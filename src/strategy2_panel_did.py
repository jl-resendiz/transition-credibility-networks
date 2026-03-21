"""Strategy 2 panel DiD: firm-month panel with firm + month FE.

Estimate: AR_{j,t} = alpha_j + gamma_t + beta * (Exposure_{j,e} * Post_{e,t}) + eps_{j,t}

Key transformations:
  - Monthly abnormal returns (vwretd) from Fama-French
  - Event-time windows around first-mover retirements
  - Exposure per event using geographic W
  - Overlap handling: keep nearest event per (firm, month)
  - Two-way FE by demeaning (firm + month)
"""
import csv
import math
import os
import random
import hashlib
from collections import defaultdict

from _paths import raw_path, derived_path

# Settings (override via environment variables)
def _env_int(name, default):
    v = os.getenv(name)
    if v is None or v == '':
        return default
    try:
        return int(v)
    except ValueError:
        return default


def _env_bool(name, default):
    v = os.getenv(name)
    if v is None or v == '':
        return default
    return v.strip().lower() in ('1', 'true', 'yes', 'y')


EVENT_SCOPE = os.getenv('EVENT_SCOPE', 'first_mover')  # 'first_mover' or 'all_matched'
EVENTS_PATH = os.getenv('EVENTS_PATH', '')             # optional override for events file
EXACT_ONLY = _env_bool('EXACT_ONLY', False)            # restrict to researched announcement dates only
TAU_START = _env_int('TAU_START', -6)
TAU_END = _env_int('TAU_END', 12)
POST_START = _env_int('POST_START', 0)
POST_END = _env_int('POST_END', 12)
CONTROL_MULT = _env_int('CONTROL_MULT', 5)             # controls per neighbor
OVERLAP_RULE = os.getenv('OVERLAP_RULE', 'nearest')    # 'nearest' or 'drop'
FOREIGN_ONLY = _env_bool('FOREIGN_ONLY', False)        # if True, zero out same-country exposure
TRANSFORM_SET = os.getenv('TRANSFORM_SET', 'base')     # base | log1p | zscore
WRITE_METRICS = _env_bool('WRITE_METRICS', False)
METRICS_SUFFIX = os.getenv('METRICS_SUFFIX', '').strip()
TIER_FILTER = os.getenv('TIER_FILTER', '').strip()    # e.g., "1" or "1,2"
BINDING_ONLY = _env_bool('BINDING_ONLY', False)
EXCLUDE_YEARS = os.getenv('EXCLUDE_YEARS', '').strip()  # e.g., "2019"


def add_months(ym, delta):
    y, m = ym.split('-')
    y = int(y)
    m = int(m)
    m += delta
    while m > 12:
        y += 1
        m -= 12
    while m < 1:
        y -= 1
        m += 12
    return f"{y:04d}-{m:02d}"


def apply_exposure_transform(obs, transform):
    if transform == 'base':
        return obs
    if transform == 'log1p':
        for row in obs:
            w = row.get('w_ij', 0.0)
            w_t = math.log1p(w) if w > 0 else 0.0
            row['w_ij'] = w_t
            row['exp_post'] = w_t * row.get('post', 0.0)
        return obs
    if transform == 'zscore':
        vals = [row.get('w_ij', 0.0) for row in obs]
        if not vals:
            return obs
        mean = sum(vals) / len(vals)
        var = sum((v - mean) ** 2 for v in vals) / len(vals)
        std = math.sqrt(var) if var > 1e-12 else 1.0
        for row in obs:
            w = row.get('w_ij', 0.0)
            w_t = (w - mean) / std
            row['w_ij'] = w_t
            row['exp_post'] = w_t * row.get('post', 0.0)
        return obs
    return obs


def write_metrics(path, res_event, res_firm, res_tw, meta):
    def _fmt(res):
        if not res:
            return 'NA'
        b = res['beta'].get('exp_post', float('nan'))
        se = res['se'].get('exp_post', float('nan'))
        t = res['t'].get('exp_post', float('nan'))
        return f'{b:+.4f} (se {se:.4f}, t {t:.2f}), N={res["n"]}'
    lines = [
        '# Strategy 2 Panel DiD Metrics',
        '',
        f'- transform: {meta.get("transform")}',
        f'- event_scope: {meta.get("event_scope")}',
        f'- exact_only: {meta.get("exact_only")}',
        f'- foreign_only: {meta.get("foreign_only")}',
        f'- tau: [{meta.get("tau_start")},{meta.get("tau_end")}]',
        f'- post: [{meta.get("post_start")},{meta.get("post_end")}]',
        f'- overlap: {meta.get("overlap_rule")}',
        '',
        '## exp_post coefficient',
        f'- event-clustered: {_fmt(res_event)}',
        f'- firm-clustered: {_fmt(res_firm)}',
        f'- two-way: {_fmt(res_tw)}',
        '',
    ]
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))


def load_ff_factors_monthly(path):
    if not os.path.exists(path):
        return None
    mktrf = {}
    rf = {}
    vwretd = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('This file') or line.startswith('The '):
                continue
            if line.startswith(','):
                continue
            parts = [p.strip() for p in line.split(',')]
            if len(parts) < 5:
                continue
            date = parts[0]
            if not date.isdigit() or len(date) != 6:
                continue
            try:
                mktrf_val = float(parts[1])
                rf_val = float(parts[4])
            except ValueError:
                continue
            mktrf_dec = mktrf_val / 100.0
            rf_dec = rf_val / 100.0
            vwretd_dec = mktrf_dec + rf_dec
            date_fmt = f'{date[:4]}-{date[4:6]}'
            mktrf[date_fmt] = mktrf_dec
            rf[date_fmt] = rf_dec
            vwretd[date_fmt] = vwretd_dec
    return (mktrf, rf, vwretd) if mktrf else None


def is_exact_source(src):
    if not src:
        return False
    s = src.lower()
    if 'proxy' in s or 'approx' in s or 'mid' in s or 'month' in s:
        return False
    return True


def invert_matrix(mat):
    n = len(mat)
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(mat)]
    for col in range(n):
        max_row = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[max_row][col]) < 1e-20:
            return None
        aug[col], aug[max_row] = aug[max_row], aug[col]
        pivot = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= pivot
        for row in range(n):
            if row != col:
                factor = aug[row][col]
                for j in range(2 * n):
                    aug[row][j] -= factor * aug[col][j]
    return [row[n:] for row in aug]


def mat_mul(a, b):
    rows = len(a)
    cols = len(b[0])
    mid = len(b)
    out = [[0.0 for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for k in range(mid):
            aik = a[i][k]
            if aik == 0:
                continue
            for j in range(cols):
                out[i][j] += aik * b[k][j]
    return out


def _cluster_cov(X, resid, clusters):
    n = len(X)
    k = len(X[0])
    S = [[0.0 for _ in range(k)] for _ in range(k)]
    for _, idxs in clusters.items():
        xu = [0.0 for _ in range(k)]
        for i in idxs:
            for a in range(k):
                xu[a] += X[i][a] * resid[i]
        for a in range(k):
            for b in range(k):
                S[a][b] += xu[a] * xu[b]
    return S, len(clusters), n


def ols(data, y_var, x_vars, cluster_var=None):
    n = len(data)
    k = len(x_vars) + 1
    if n <= k + 1:
        return None
    y = [d[y_var] for d in data]
    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]
    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None
    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]
    y_hat = [sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]

    if cluster_var:
        if isinstance(cluster_var, (list, tuple)) and len(cluster_var) == 2:
            c1, c2 = cluster_var
            clusters1, clusters2, clusters12 = {}, {}, {}
            for i, d in enumerate(data):
                k1 = d.get(c1)
                k2 = d.get(c2)
                clusters1.setdefault(k1, []).append(i)
                clusters2.setdefault(k2, []).append(i)
                clusters12.setdefault((k1, k2), []).append(i)
            S1, G1, _ = _cluster_cov(X, resid, clusters1)
            S2, G2, _ = _cluster_cov(X, resid, clusters2)
            S12, G12, _ = _cluster_cov(X, resid, clusters12)
            S = [[S1[a][b] + S2[a][b] - S12[a][b] for b in range(k)] for a in range(k)]
            V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
            G = min(G1, G2)
            if G > 1:
                scale = (G / (G - 1)) * ((n - 1) / (n - k))
                for a in range(k):
                    for b in range(k):
                        V[a][b] *= scale
            se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
        else:
            clusters = {}
            for i, d in enumerate(data):
                cid = d.get(cluster_var)
                clusters.setdefault(cid, []).append(i)
            S, G, _ = _cluster_cov(X, resid, clusters)
            V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
            if G > 1:
                scale = (G / (G - 1)) * ((n - 1) / (n - k))
                for a in range(k):
                    for b in range(k):
                        V[a][b] *= scale
            se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
    else:
        se = [0.0 for _ in range(k)]

    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0.0 for a in range(k)]
    names = ['intercept'] + x_vars
    return {
        'beta': dict(zip(names, beta)),
        'se': dict(zip(names, se)),
        't': dict(zip(names, t_stats)),
        'n': n,
    }


def two_way_demean(data, y_var, x_vars, firm_key='gvkey', time_key='ym'):
    vars_all = [y_var] + x_vars
    overall = {v: 0.0 for v in vars_all}
    firm_sums = defaultdict(lambda: defaultdict(float))
    firm_counts = defaultdict(int)
    time_sums = defaultdict(lambda: defaultdict(float))
    time_counts = defaultdict(int)

    for row in data:
        f = row[firm_key]
        t = row[time_key]
        firm_counts[f] += 1
        time_counts[t] += 1
        for v in vars_all:
            val = row[v]
            overall[v] += val
            firm_sums[f][v] += val
            time_sums[t][v] += val

    n = len(data)
    if n == 0:
        return []
    for v in vars_all:
        overall[v] /= n

    firm_means = {f: {v: firm_sums[f][v] / firm_counts[f] for v in vars_all}
                  for f in firm_counts}
    time_means = {t: {v: time_sums[t][v] / time_counts[t] for v in vars_all}
                  for t in time_counts}

    out = []
    for row in data:
        f = row[firm_key]
        t = row[time_key]
        d = {firm_key: f, time_key: t}
        for v in vars_all:
            d[v] = row[v] - firm_means[f][v] - time_means[t][v] + overall[v]
        d['event_id'] = row.get('event_id')
        d['tau'] = row.get('tau')
        out.append(d)
    return out


# ----------------------------- Load data -----------------------------

monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except ValueError:
            pass

ff_monthly_path = raw_path('factors', 'F-F_Research_Data_Factors.csv')
ff_factors = load_ff_factors_monthly(ff_monthly_path)
if not ff_factors:
    raise RuntimeError('Missing F-F monthly factors for vwretd.')
mktrf_monthly, rf_monthly, market_ret = ff_factors

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


# ----------------------------- Build panel -----------------------------

all_events = []
events_path = EVENTS_PATH if EVENTS_PATH else derived_path('events', 'coal_retirement_events.csv')
with open(events_path, 'r', encoding='utf-8') as f:
    tiers = set([t.strip() for t in TIER_FILTER.split(',') if t.strip()]) if TIER_FILTER else None
    exclude_years = set([y.strip() for y in EXCLUDE_YEARS.split(',') if y.strip()]) if EXCLUDE_YEARS else set()
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        ann_date = row.get('announcement_date', '').strip()
        ret_date = row.get('event_date', '').strip()
        ann_src = row.get('announcement_source', '').strip()
        if EXACT_ONLY:
            if not ann_date:
                continue
            if ann_src and not is_exact_source(ann_src):
                continue
        effective_date = ann_date if ann_date else ret_date
        if effective_date and len(effective_date) >= 7:
            event_month = effective_date[:7]
        else:
            event_month = f'{row["ret_year"]}-07'
        event_year = int(event_month[:4]) if event_month[:4].isdigit() else None

        # optional filters for regulatory shocks
        if tiers is not None and row.get('exogeneity_tier', ''):
            if row.get('exogeneity_tier', '').strip() not in tiers:
                continue
        if BINDING_ONLY and row.get('binding', ''):
            if row.get('binding', '').strip().lower() != 'yes':
                continue
        if exclude_years and event_month[:4] in exclude_years:
            continue

        all_events.append({
            'event_month': event_month,
            'event_year': event_year,
            'gvkeys': row['matched_gvkeys'].split(';'),
            'is_first_mover': row.get('is_first_mover') == 'True',
        })

if EVENT_SCOPE == 'first_mover':
    events = [e for e in all_events if e['is_first_mover']]
else:
    events = list(all_events)

obs = []
for event_id, event in enumerate(events):
    event_month = event['event_month']
    if not event_month:
        continue
    event_gvkeys = set(event['gvkeys'])
    for fm_gk in event_gvkeys:
        if fm_gk not in W:
            continue
        fm_fic = None
        frow_fm = fundamentals.get(fm_gk)
        if frow_fm and frow_fm.get('fic'):
            fm_fic = frow_fm.get('fic')
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
            if FOREIGN_ONLY and fm_fic:
                frow_j = fundamentals.get(gk)
                fic_j = frow_j.get('fic') if frow_j else None
                if fic_j and fic_j == fm_fic:
                    w_ij = 0.0
            if gk not in monthly_ret:
                continue
            for tau in range(TAU_START, TAU_END + 1):
                ym = add_months(event_month, tau)
                if ym not in monthly_ret[gk] or ym not in market_ret:
                    continue
                ar = monthly_ret[gk][ym] - market_ret[ym]
                post = 1.0 if (tau >= POST_START and tau <= POST_END) else 0.0
                obs.append({
                    'gvkey': gk,
                    'ym': ym,
                    'event_id': event_id,
                    'tau': tau,
                    'ar': ar,
                    'post': post,
                    'exp_post': w_ij * post,
                    'w_ij': w_ij,
                })

print(f'Raw observations: {len(obs)}')

# Overlap handling (firm-month appears in multiple events)
if OVERLAP_RULE in ('nearest', 'drop'):
    by_key = defaultdict(list)
    for row in obs:
        by_key[(row['gvkey'], row['ym'])].append(row)
    filtered = []
    overlaps = 0
    for key, rows in by_key.items():
        if len(rows) == 1:
            filtered.append(rows[0])
            continue
        overlaps += 1
        if OVERLAP_RULE == 'drop':
            continue
        # keep nearest event by |tau|
        rows_sorted = sorted(rows, key=lambda r: abs(r['tau']))
        filtered.append(rows_sorted[0])
    obs = filtered
    print(f'Overlap groups: {overlaps}, kept obs: {len(obs)}')

# Apply exposure transformation after overlap handling
obs = apply_exposure_transform(obs, TRANSFORM_SET)

# Two-way demean (firm + month FE)
obs_dm = two_way_demean(obs, 'ar', ['exp_post'], firm_key='gvkey', time_key='ym')

# Regress with different clustering
res_event = ols(obs_dm, 'ar', ['exp_post'], cluster_var='event_id')
res_firm = ols(obs_dm, 'ar', ['exp_post'], cluster_var='gvkey')
res_tw = ols(obs_dm, 'ar', ['exp_post'], cluster_var=['event_id', 'gvkey'])

def _print_res(label, res):
    if not res:
        print(f'{label}: regression failed')
        return
    b = res['beta'].get('exp_post', 0.0)
    se = res['se'].get('exp_post', 0.0)
    t = res['t'].get('exp_post', 0.0)
    print(f'{label}: beta={b:+.4f}  se={se:.4f}  t={t:.2f}  N={res["n"]}')

print('\n=== Panel DiD (firm + month FE) ===')
print(f'Event scope={EVENT_SCOPE}, tau=[{TAU_START},{TAU_END}], post=[{POST_START},{POST_END}], overlap={OVERLAP_RULE}')
_print_res('Event-clustered', res_event)
_print_res('Firm-clustered', res_firm)
_print_res('Two-way clustered', res_tw)

if WRITE_METRICS:
    suffix = f'_{METRICS_SUFFIX}' if METRICS_SUFFIX else ''
    metrics_path = os.path.join('JEEM_submission_package', 'JEEM_outputs', 'metrics',
                                f'strategy2_panel_did_metrics{suffix}_{TRANSFORM_SET}.md')
    write_metrics(metrics_path, res_event, res_firm, res_tw, {
        'transform': TRANSFORM_SET,
        'event_scope': EVENT_SCOPE,
        'exact_only': EXACT_ONLY,
        'foreign_only': FOREIGN_ONLY,
        'tau_start': TAU_START,
        'tau_end': TAU_END,
        'post_start': POST_START,
        'post_end': POST_END,
        'overlap_rule': OVERLAP_RULE,
    })
