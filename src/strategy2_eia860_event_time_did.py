"""Event-time DiD for EIA-860-style announcement shocks.

For each event-time tau, estimate:
  AR_{j,tau} = a_tau + b_tau * exposure_{j,e} + eps
with event-clustered SEs. This produces an event-time profile of the exposure effect.

Outputs: results/summaries/strategy2_eia860_event_time.csv
"""
import csv
import math
import os
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

EVENTS_PATH = os.getenv('EVENTS_PATH', '')
EVENT_SCOPE = os.getenv('EVENT_SCOPE', 'all_matched')
CONTROL_MULT = int(os.getenv('CONTROL_MULT', '5'))
TAU_START = int(os.getenv('TAU_START', '-6'))
TAU_END = int(os.getenv('TAU_END', '12'))
TRANSFORM_SET = os.getenv('TRANSFORM_SET', 'base')  # base | log1p | zscore


def load_ff_factors_monthly(path):
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
            if not date.isdigit() or len(date) != 6:
                continue
            try:
                mktrf_val = float(parts[1])
                rf_val = float(parts[4])
            except ValueError:
                continue
            vw = (mktrf_val + rf_val) / 100.0
            date_fmt = f'{date[:4]}-{date[4:6]}'
            vwretd[date_fmt] = vw
    return vwretd


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


def apply_transform(vals, transform):
    if transform == 'base':
        return vals
    if transform == 'log1p':
        return [math.log1p(v) if v > 0 else 0.0 for v in vals]
    if transform == 'zscore':
        mean = sum(vals) / len(vals) if vals else 0.0
        var = sum((v - mean) ** 2 for v in vals) / len(vals) if vals else 0.0
        std = math.sqrt(var) if var > 1e-12 else 1.0
        return [(v - mean) / std for v in vals]
    return vals


def ols_cluster(y, x, clusters):
    n = len(y)
    if n == 0:
        return None
    # x is list of [1, w]
    k = 2
    XtX = [[0.0, 0.0], [0.0, 0.0]]
    Xty = [0.0, 0.0]
    for i in range(n):
        xi0, xi1 = x[i]
        XtX[0][0] += xi0 * xi0
        XtX[0][1] += xi0 * xi1
        XtX[1][0] += xi1 * xi0
        XtX[1][1] += xi1 * xi1
        Xty[0] += xi0 * y[i]
        Xty[1] += xi1 * y[i]
    det = XtX[0][0] * XtX[1][1] - XtX[0][1] * XtX[1][0]
    if abs(det) < 1e-12:
        return None
    inv = [[XtX[1][1] / det, -XtX[0][1] / det],
           [-XtX[1][0] / det, XtX[0][0] / det]]
    beta0 = inv[0][0] * Xty[0] + inv[0][1] * Xty[1]
    beta1 = inv[1][0] * Xty[0] + inv[1][1] * Xty[1]
    resid = [y[i] - beta0 * x[i][0] - beta1 * x[i][1] for i in range(n)]

    # cluster-robust SE
    clus = defaultdict(list)
    for i, cid in enumerate(clusters):
        clus[cid].append(i)
    S = [[0.0, 0.0], [0.0, 0.0]]
    for _, idxs in clus.items():
        xu0 = 0.0
        xu1 = 0.0
        for i in idxs:
            xu0 += x[i][0] * resid[i]
            xu1 += x[i][1] * resid[i]
        S[0][0] += xu0 * xu0
        S[0][1] += xu0 * xu1
        S[1][0] += xu1 * xu0
        S[1][1] += xu1 * xu1
    cov = [
        [inv[0][0] * S[0][0] * inv[0][0] + inv[0][1] * S[1][0] * inv[0][0] + inv[0][0] * S[0][1] * inv[0][1] + inv[0][1] * S[1][1] * inv[0][1],
         inv[0][0] * S[0][0] * inv[1][0] + inv[0][1] * S[1][0] * inv[1][0] + inv[0][0] * S[0][1] * inv[1][1] + inv[0][1] * S[1][1] * inv[1][1]],
        [inv[1][0] * S[0][0] * inv[0][0] + inv[1][1] * S[1][0] * inv[0][0] + inv[1][0] * S[0][1] * inv[0][1] + inv[1][1] * S[1][1] * inv[0][1],
         inv[1][0] * S[0][0] * inv[1][0] + inv[1][1] * S[1][0] * inv[1][0] + inv[1][0] * S[0][1] * inv[1][1] + inv[1][1] * S[1][1] * inv[1][1]]
    ]
    se1 = math.sqrt(cov[1][1]) if cov[1][1] > 0 else float('nan')
    t1 = beta1 / se1 if se1 and se1 > 0 else float('nan')
    return beta1, se1, t1, len(clus), n


# Load monthly returns
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except (ValueError, TypeError):
            continue

vwretd = load_ff_factors_monthly(raw_path('factors', 'F-F_Research_Data_Factors.csv'))
if not vwretd:
    raise RuntimeError('Missing Fama-French monthly factors for vwretd.')

# Load W_geo
W = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gi = row['gvkey_i']
        gj = row['gvkey_j']
        try:
            W[gi][gj] = float(row['w_ij'])
        except (ValueError, TypeError):
            continue

# Fundamentals for control pool
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
        if not ann_date:
            continue
        effective_date = ann_date
        if len(effective_date) >= 7:
            event_month = effective_date[:7]
        else:
            continue
        all_events.append({
            'event_month': event_month,
            'gvkeys': [g.strip() for g in row['matched_gvkeys'].split(';') if g.strip()],
            'is_first_mover': row.get('is_first_mover') == 'True',
        })

events = [e for e in all_events if e['is_first_mover']] if EVENT_SCOPE == 'first_mover' else list(all_events)

# Build event-time observations
rows_by_tau = {tau: [] for tau in range(TAU_START, TAU_END + 1)}
cluster_by_tau = {tau: [] for tau in range(TAU_START, TAU_END + 1)}

for event_id, e in enumerate(events):
    event_month = e['event_month']
    event_gvkeys = set(e['gvkeys'])
    for fm_gk in event_gvkeys:
        if fm_gk not in W:
            continue
        neighbors = W[fm_gk]
        neighbor_gks = set(neighbors.keys()) - event_gvkeys
        non_connected = [gk for gk in fundamentals if gk not in event_gvkeys and gk not in neighbors]
        n_ctrl = min(len(non_connected), max(CONTROL_MULT * len(neighbor_gks), 20))
        ctrl_sample = non_connected[:n_ctrl]
        candidate_firms = list(neighbor_gks) + ctrl_sample
        for gk in candidate_firms:
            w = neighbors.get(gk, 0.0)
            for tau in range(TAU_START, TAU_END + 1):
                ym = add_months(event_month, tau)
                if gk not in monthly_ret or ym not in monthly_ret[gk]:
                    continue
                ar = monthly_ret[gk][ym] - vwretd.get(ym, 0.0)
                rows_by_tau[tau].append((ar, w))
                cluster_by_tau[tau].append(event_id)

# Apply exposure transform within each tau (consistent with global transform)
out_rows = []
for tau in range(TAU_START, TAU_END + 1):
    data = rows_by_tau[tau]
    if not data:
        continue
    y = [d[0] for d in data]
    w_raw = [d[1] for d in data]
    w = apply_transform(w_raw, TRANSFORM_SET)
    X = [[1.0, w[i]] for i in range(len(w))]
    res = ols_cluster(y, X, cluster_by_tau[tau])
    if res is None:
        continue
    beta, se, tstat, n_clusters, n_obs = res
    out_rows.append({
        'tau': tau,
        'beta_w': beta,
        'se_w': se,
        't_w': tstat,
        'n_obs': n_obs,
        'n_clusters': n_clusters,
    })

suffix = f'_{TRANSFORM_SET}' if TRANSFORM_SET else ''
out_path = results_path('summaries', f'strategy2_eia860_event_time{suffix}.csv')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['tau', 'beta_w', 'se_w', 't_w', 'n_obs', 'n_clusters'])
    w.writeheader()
    for row in out_rows:
        w.writerow(row)

print(f'Wrote: {out_path} (transform={TRANSFORM_SET}, events={len(events)})')
