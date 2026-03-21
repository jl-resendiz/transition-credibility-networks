"""Volatility mediation for EIA-860-style announcement shocks.

Computes daily abnormal return CAR and pre/post volatility around announcement dates,
then estimates mediation-style regressions:
  (1) vol_change ~ exposure
  (2) CAR ~ exposure
  (3) CAR ~ exposure + vol_change

Outputs: results/metrics/strategy2_eia860_vol_mediation.md
"""
import csv
import math
import os
import random
import hashlib
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

EVENTS_PATH = os.getenv('EVENTS_PATH', '')
EVENT_SCOPE = os.getenv('EVENT_SCOPE', 'all_matched')
CONTROL_MULT = int(os.getenv('CONTROL_MULT', '5'))
TRANSFORM_SET = os.getenv('TRANSFORM_SET', 'base')  # base | log1p | zscore
WRITE_METRICS = os.getenv('WRITE_METRICS', '1') == '1'
METRICS_SUFFIX = os.getenv('METRICS_SUFFIX', '').strip()
TIER_FILTER = os.getenv('TIER_FILTER', '')  # optional, e.g. "1" or "1,2"
BINDING_ONLY = os.getenv('BINDING_ONLY', '0') == '1'

PRE_START = -120
PRE_END = -20
POST_START = 0
POST_END = 60
CAR_START = -1
CAR_END = 20


def load_ff_factors_daily(path):
    if not os.path.exists(path):
        return None
    vwretd = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('This file') or line.startswith('The ') or line.startswith(','):
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 5:
                continue
            date = parts[0]
            if not date.isdigit() or len(date) != 8:
                continue
            try:
                mktrf_val = float(parts[1])
                rf_val = float(parts[4])
            except ValueError:
                continue
            vw = (mktrf_val + rf_val) / 100.0
            date_fmt = f"{date[:4]}-{date[4:6]}-{date[6:]}"
            vwretd[date_fmt] = vw
    return vwretd if vwretd else None


def load_daily_returns(path):
    data = defaultdict(dict)
    with open(path, 'r', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            gvkey = row['gvkey']
            date = row['datadate']
            try:
                ret = float(row['ret_daily'])
            except (ValueError, TypeError):
                continue
            data[gvkey][date] = ret
    return data


def window_stats(dates, ar, event_idx, start, end):
    if not dates or not ar:
        return None
    idx = event_idx
    if idx is None:
        return None
    vals = []
    for offset in range(start, end + 1):
        j = idx + offset
        if 0 <= j < len(dates):
            vals.append(ar[j])
    if len(vals) < (end - start + 1) * 0.4:
        return None
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / max(len(vals) - 1, 1)
    return mean, math.sqrt(var), sum(vals)


def apply_exposure_transform(obs, transform):
    if transform == 'base':
        return obs
    if transform == 'log1p':
        for row in obs:
            w = row.get('w_ij', 0.0)
            row['w_ij'] = math.log1p(w) if w > 0 else 0.0
        return obs
    if transform == 'zscore':
        vals = [row.get('w_ij', 0.0) for row in obs]
        if not vals:
            return obs
        mean = sum(vals) / len(vals)
        var = sum((v - mean) ** 2 for v in vals) / len(vals)
        std = math.sqrt(var) if var > 1e-12 else 1.0
        for row in obs:
            row['w_ij'] = (row.get('w_ij', 0.0) - mean) / std
        return obs
    return obs


def ols_cluster(y, X, cluster_ids):
    n = len(y)
    k = len(X[0])
    # OLS
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

    # invert
    def invert_matrix(A):
        m = len(A)
        aug = [[A[i][j] for j in range(m)] + [1.0 if i == j else 0.0 for j in range(m)] for i in range(m)]
        for i in range(m):
            pivot = aug[i][i]
            if abs(pivot) < 1e-12:
                return None
            inv_p = 1.0 / pivot
            for j in range(2*m):
                aug[i][j] *= inv_p
            for r in range(m):
                if r == i:
                    continue
                factor = aug[r][i]
                for c in range(2*m):
                    aug[r][c] -= factor * aug[i][c]
        return [row[m:] for row in aug]

    inv = invert_matrix(XtX)
    if inv is None:
        return None
    beta = [sum(inv[i][j] * Xty[j] for j in range(k)) for i in range(k)]
    yhat = [sum(beta[j] * X[i][j] for j in range(k)) for i in range(n)]
    resid = [y[i] - yhat[i] for i in range(n)]

    # cluster covariance
    clusters = defaultdict(list)
    for i, cid in enumerate(cluster_ids):
        clusters[cid].append(i)
    S = [[0.0 for _ in range(k)] for _ in range(k)]
    for _, idxs in clusters.items():
        xu = [0.0 for _ in range(k)]
        for i in idxs:
            for a in range(k):
                xu[a] += X[i][a] * resid[i]
        for a in range(k):
            for b in range(k):
                S[a][b] += xu[a] * xu[b]
    cov = [[sum(inv[i][a] * S[a][b] * inv[j][b] for a in range(k) for b in range(k))
            for j in range(k)] for i in range(k)]
    se = [math.sqrt(cov[i][i]) if cov[i][i] > 0 else float('nan') for i in range(k)]
    return beta, se, len(clusters)


# Load inputs
ff_daily = raw_path('factors', 'F-F_Research_Data_Factors_daily.csv')
market_ret = load_ff_factors_daily(ff_daily)
if not market_ret:
    raise RuntimeError('Missing F-F daily factors for vwretd.')

daily_ret = load_daily_returns(derived_path('returns', 'daily_returns.csv'))

# Precompute AR series per gvkey for speed
daily_ar = {}
for gk, dct in daily_ret.items():
    dates = [d for d in dct.keys() if d in market_ret]
    if not dates:
        continue
    dates.sort()
    ar = [dct[d] - market_ret[d] for d in dates]
    daily_ar[gk] = (dates, ar)

W = defaultdict(dict)
geo_path = derived_path('networks', 'weight_matrix_W_geo.csv')
with open(geo_path, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

# Events
events_path = EVENTS_PATH if EVENTS_PATH else derived_path('events', 'eia860_announcement_events.csv')
all_events = []
with open(events_path, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        ann_date = row.get('announcement_date', '').strip()
        event_date = row.get('event_date', '').strip()
        if not ann_date and not event_date:
            continue
        if TIER_FILTER:
            tiers = set([t.strip() for t in TIER_FILTER.split(',') if t.strip()])
            if row.get('exogeneity_tier', '').strip() not in tiers:
                continue
        if BINDING_ONLY and row.get('binding', '').strip().lower() != 'yes':
            continue
        effective_date = ann_date if ann_date else event_date
        if len(effective_date) == 7 and effective_date[4] == '-':
            effective_date = f"{effective_date}-15"
        all_events.append({
            'event_date': effective_date,
            'gvkeys': row['matched_gvkeys'].split(';'),
            'is_first_mover': row.get('is_first_mover') == 'True',
        })

if EVENT_SCOPE == 'first_mover':
    events = [e for e in all_events if e['is_first_mover']]
else:
    events = list(all_events)

obs = []
for event_id, event in enumerate(events):
    event_date = event['event_date']
    event_gvkeys = set([g for g in event['gvkeys'] if g])
    for fm_gk in event_gvkeys:
        if fm_gk not in W:
            continue
        neighbors = W[fm_gk]
        neighbor_gks = set(neighbors.keys()) - event_gvkeys
        non_connected = [gk for gk in fundamentals if gk not in event_gvkeys and gk not in neighbors]
        stable_seed = int(hashlib.md5(str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
        random.seed(stable_seed)
        n_ctrl = min(len(non_connected), max(CONTROL_MULT * len(neighbor_gks), 20))
        ctrl_sample = random.sample(non_connected, n_ctrl) if len(non_connected) > n_ctrl else non_connected
        candidate_firms = list(neighbor_gks) + ctrl_sample

        for gk in candidate_firms:
            w_ij = neighbors.get(gk, 0.0)
            if gk not in daily_ar:
                continue
            dates, ar = daily_ar[gk]
            # locate event index by first date >= event_date
            event_idx = None
            for i, d in enumerate(dates):
                if d >= event_date:
                    event_idx = i
                    break
            if event_idx is None:
                continue
            pre = window_stats(dates, ar, event_idx, PRE_START, PRE_END)
            post = window_stats(dates, ar, event_idx, POST_START, POST_END)
            car = window_stats(dates, ar, event_idx, CAR_START, CAR_END)
            if pre is None or post is None or car is None:
                continue
            vol_pre = pre[1]
            vol_post = post[1]
            vol_chg = vol_post - vol_pre
            car_sum = car[2]
            obs.append({
                'event_id': event_id,
                'gvkey': gk,
                'w_ij': w_ij,
                'vol_chg': vol_chg,
                'car': car_sum,
            })

apply_exposure_transform(obs, TRANSFORM_SET)

# Handle empty sample
if not obs:
    if WRITE_METRICS:
        suffix = f'_{TRANSFORM_SET}' if TRANSFORM_SET else ''
        extra = f'_{METRICS_SUFFIX}' if METRICS_SUFFIX else ''
        out_path = results_path('metrics', f'strategy2_eia860_vol_mediation{extra}{suffix}.md')
        os.makedirs(os.path.dirname(out_path), exist_ok=True)
        with open(out_path, 'w', encoding='utf-8') as f:
            f.write('# EIA860 Volatility Mediation\n\n- N: 0 (no usable observations)\n')
        print(f'Wrote: {out_path}')
    raise SystemExit('No usable observations for volatility mediation.')

# Regression 1: vol_change ~ exposure
X1 = [[1.0, o['w_ij']] for o in obs]
y1 = [o['vol_chg'] for o in obs]
res1 = ols_cluster(y1, X1, [o['event_id'] for o in obs])

# Regression 2: CAR ~ exposure
X2 = [[1.0, o['w_ij']] for o in obs]
y2 = [o['car'] for o in obs]
res2 = ols_cluster(y2, X2, [o['event_id'] for o in obs])

# Regression 3: CAR ~ exposure + vol_change
X3 = [[1.0, o['w_ij'], o['vol_chg']] for o in obs]
y3 = [o['car'] for o in obs]
res3 = ols_cluster(y3, X3, [o['event_id'] for o in obs])

if WRITE_METRICS:
    suffix = f'_{TRANSFORM_SET}' if TRANSFORM_SET else ''
    extra = f'_{METRICS_SUFFIX}' if METRICS_SUFFIX else ''
    out_path = results_path('metrics', f'strategy2_eia860_vol_mediation{extra}{suffix}.md')
    lines = [
        '# EIA860 Volatility Mediation',
        '',
        f'- event_scope: {EVENT_SCOPE}',
        f'- transform: {TRANSFORM_SET}',
        f'- N: {len(obs)}',
        '',
    ]
    if res1:
        b, se, g = res1
        lines.append('## Volatility change ~ exposure')
        lines.append(f'beta_w={b[1]:+.4f}, se={se[1]:.4f}, t={b[1]/se[1]:.2f}, clusters={g}')
        lines.append('')
    if res2:
        b, se, g = res2
        lines.append('## CAR ~ exposure')
        lines.append(f'beta_w={b[1]:+.4f}, se={se[1]:.4f}, t={b[1]/se[1]:.2f}, clusters={g}')
        lines.append('')
    if res3:
        b, se, g = res3
        lines.append('## CAR ~ exposure + vol_change')
        lines.append(f'beta_w={b[1]:+.4f}, se={se[1]:.4f}, t={b[1]/se[1]:.2f}, clusters={g}')
        lines.append(f'beta_vol={b[2]:+.4f}, se={se[2]:.4f}, t={b[2]/se[2]:.2f}, clusters={g}')
        lines.append('')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'Wrote: {out_path}')
