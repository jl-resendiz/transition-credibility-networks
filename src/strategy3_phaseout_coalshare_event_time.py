"""Event-time DiD for coal phase-out shocks using coal-share intensity.

Estimates AR_{i,tau} = a_tau + b_tau * coal_share_i * Treat_i + eps
with event-clustered SEs. coal_share is matched to event year.

Outputs: results/summaries/strategy3_phaseout_coalshare_event_time_{transform}.csv
"""
import csv
import math
import os
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

EVENTS_PATH = os.getenv('EVENTS_PATH', derived_path('events', 'coal_phaseout_shocks_events.csv'))
TIER_FILTER = os.getenv('TIER_FILTER', '').strip()
BINDING_ONLY = os.getenv('BINDING_ONLY', '0') == '1'
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
            vwretd[f'{date[:4]}-{date[4:6]}'] = vw
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


def load_coal_share(panel_path):
    years_by_gvkey = defaultdict(list)
    coal_by_year = defaultdict(dict)
    with open(panel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            yr = row['year']
            cs = row.get('coal_share', '')
            if not cs:
                continue
            try:
                cs_val = float(cs)
            except (ValueError, TypeError):
                continue
            coal_by_year[gk][int(yr)] = cs_val
            years_by_gvkey[gk].append(int(yr))
    for gk in years_by_gvkey:
        years_by_gvkey[gk] = sorted(set(years_by_gvkey[gk]))

    def get_share(gk, year):
        if gk not in coal_by_year or not years_by_gvkey[gk]:
            return None
        if year in coal_by_year[gk]:
            return coal_by_year[gk][year]
        years = years_by_gvkey[gk]
        prior = [y for y in years if y <= year]
        if prior:
            return coal_by_year[gk][max(prior)]
        return coal_by_year[gk][years[0]]

    return get_share


# Load returns
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
    raise RuntimeError('Missing F-F monthly factors for vwretd.')

# Fundamentals for control pool
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

get_coal_share = load_coal_share(derived_path('fundamentals', 'firm_alpha_panel.csv'))

# Events
tiers = set([t.strip() for t in TIER_FILTER.split(',') if t.strip()]) if TIER_FILTER else None
events = []
with open(EVENTS_PATH, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        if tiers and row.get('exogeneity_tier', '') not in tiers:
            continue
        if BINDING_ONLY and row.get('binding', '').strip().lower() != 'yes':
            continue
        event_date = row.get('event_date', '')
        if not event_date:
            continue
        event_month = event_date[:7]
        event_year = int(event_date[:4]) if event_date[:4].isdigit() else None
        events.append({
            'event_month': event_month,
            'event_year': event_year,
            'gvkeys': [g for g in row['matched_gvkeys'].split(';') if g],
        })

# Build event-time rows
rows_by_tau = {tau: [] for tau in range(TAU_START, TAU_END + 1)}
cluster_by_tau = {tau: [] for tau in range(TAU_START, TAU_END + 1)}

for event_id, e in enumerate(events):
    event_month = e['event_month']
    event_year = e['event_year']
    if event_year is None:
        continue
    event_gvkeys = set(e['gvkeys'])
    non_treated = [gk for gk in fundamentals if gk not in event_gvkeys]
    n_ctrl = min(len(non_treated), max(CONTROL_MULT * len(event_gvkeys), 20))
    ctrl_sample = non_treated[:n_ctrl]
    candidate_firms = list(event_gvkeys) + ctrl_sample
    for gk in candidate_firms:
        cs = get_coal_share(gk, event_year)
        if cs is None:
            continue
        exp = cs if gk in event_gvkeys else 0.0
        for tau in range(TAU_START, TAU_END + 1):
            ym = add_months(event_month, tau)
            if gk not in monthly_ret or ym not in monthly_ret[gk]:
                continue
            ar = monthly_ret[gk][ym] - vwretd.get(ym, 0.0)
            rows_by_tau[tau].append((ar, exp))
            cluster_by_tau[tau].append(event_id)

out_rows = []
for tau in range(TAU_START, TAU_END + 1):
    data = rows_by_tau[tau]
    if not data:
        continue
    y = [d[0] for d in data]
    exp_raw = [d[1] for d in data]
    exp = apply_transform(exp_raw, TRANSFORM_SET)
    X = [[1.0, exp[i]] for i in range(len(exp))]
    res = ols_cluster(y, X, cluster_by_tau[tau])
    if res is None:
        continue
    beta, se, tstat, n_clusters, n_obs = res
    out_rows.append({
        'tau': tau,
        'beta_exp': beta,
        'se_exp': se,
        't_exp': tstat,
        'n_obs': n_obs,
        'n_clusters': n_clusters,
    })

suffix = f'_{TRANSFORM_SET}' if TRANSFORM_SET else ''
tf = f"_tier{TIER_FILTER.replace(',', '')}" if TIER_FILTER else ''
bo = "_binding" if BINDING_ONLY else ''
out_path = results_path('summaries', f'strategy3_phaseout_coalshare_event_time{suffix}{tf}{bo}.csv')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['tau', 'beta_exp', 'se_exp', 't_exp', 'n_obs', 'n_clusters'])
    w.writeheader()
    for row in out_rows:
        w.writerow(row)

print(f'Wrote: {out_path} (events={len(events)}, transform={TRANSFORM_SET})')
