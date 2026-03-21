"""Wild cluster bootstrap for panel DiD (exp_post).

Target: Tier1-binding + log1p exposure (primary spec).
Clusters: event_id.
"""
import csv
import math
import os
import random
import hashlib
from collections import defaultdict

from _paths import raw_path, derived_path, results_path

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


EVENTS_PATH = os.getenv('EVENTS_PATH', '')
EVENT_SCOPE = os.getenv('EVENT_SCOPE', 'all_matched')
TAU_START = _env_int('TAU_START', -6)
TAU_END = _env_int('TAU_END', 12)
POST_START = _env_int('POST_START', 0)
POST_END = _env_int('POST_END', 12)
CONTROL_MULT = _env_int('CONTROL_MULT', 5)
OVERLAP_RULE = os.getenv('OVERLAP_RULE', 'nearest')
TRANSFORM_SET = os.getenv('TRANSFORM_SET', 'log1p')
TIER_FILTER = os.getenv('TIER_FILTER', '1')
BINDING_ONLY = _env_bool('BINDING_ONLY', True)
B = _env_int('B', 999)
SEED = _env_int('SEED', 42)


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


def two_way_demean(data, y_var, x_var, firm_key='gvkey', time_key='ym'):
    vars_all = [y_var, x_var]
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
        out.append(d)
    return out


def cluster_se(x, resid, clusters):
    n = len(x)
    if n == 0:
        return None
    x2 = sum(v * v for v in x)
    if x2 <= 1e-12:
        return None
    # meat
    clus = defaultdict(list)
    for i, cid in enumerate(clusters):
        clus[cid].append(i)
    S = 0.0
    for _, idxs in clus.items():
        xu = sum(x[i] * resid[i] for i in idxs)
        S += xu * xu
    se = math.sqrt(S / (x2 * x2))
    return se, len(clus)


# Load monthly returns
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except ValueError:
            pass

vwretd = load_ff_factors_monthly(raw_path('factors', 'F-F_Research_Data_Factors.csv'))
if not vwretd:
    raise RuntimeError('Missing F-F monthly factors for vwretd.')

W = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
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
events_path = EVENTS_PATH if EVENTS_PATH else derived_path('events', 'coal_phaseout_shocks_events.csv')
tiers = set([t.strip() for t in TIER_FILTER.split(',') if t.strip()]) if TIER_FILTER else None
events = []
with open(events_path, 'r', encoding='utf-8') as f:
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
        events.append({
            'event_month': event_month,
            'gvkeys': [g for g in row['matched_gvkeys'].split(';') if g],
        })

obs = []
for event_id, event in enumerate(events):
    event_month = event['event_month']
    if not event_month:
        continue
    event_gvkeys = set(event['gvkeys'])
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
            if gk not in monthly_ret:
                continue
            for tau in range(TAU_START, TAU_END + 1):
                ym = add_months(event_month, tau)
                if ym not in monthly_ret[gk]:
                    continue
                ar = monthly_ret[gk][ym] - vwretd.get(ym, 0.0)
                post = 1.0 if (tau >= POST_START and tau <= POST_END) else 0.0
                if TRANSFORM_SET == 'log1p':
                    w_t = math.log1p(w_ij) if w_ij > 0 else 0.0
                elif TRANSFORM_SET == 'zscore':
                    w_t = w_ij  # temporary, will zscore later if needed
                else:
                    w_t = w_ij
                obs.append({
                    'gvkey': gk,
                    'ym': ym,
                    'ar': ar,
                    'event_id': event_id,
                    'exp_post': w_t * post,
                })

# Overlap: keep one obs per firm-month
if OVERLAP_RULE == 'nearest':
    grouped = defaultdict(list)
    for row in obs:
        key = (row['gvkey'], row['ym'])
        grouped[key].append(row)
    obs = []
    for _, rows in grouped.items():
        rows_sorted = sorted(rows, key=lambda r: abs(r['exp_post']), reverse=True)
        obs.append(rows_sorted[0])

# zscore transform (if selected)
if TRANSFORM_SET == 'zscore' and obs:
    vals = [row['exp_post'] for row in obs]
    mean = sum(vals) / len(vals)
    var = sum((v - mean) ** 2 for v in vals) / len(vals)
    std = math.sqrt(var) if var > 1e-12 else 1.0
    for row in obs:
        row['exp_post'] = (row['exp_post'] - mean) / std

# Two-way demean
obs_dm = two_way_demean(obs, 'ar', 'exp_post', firm_key='gvkey', time_key='ym')
x = [row['exp_post'] for row in obs_dm]
y = [row['ar'] for row in obs_dm]
clusters = [row['event_id'] for row in obs_dm]

# OLS slope
x2 = sum(v * v for v in x)
beta = sum(x[i] * y[i] for i in range(len(x))) / x2 if x2 > 0 else 0.0
resid = [y[i] - beta * x[i] for i in range(len(x))]
se_res = cluster_se(x, resid, clusters)
se = se_res[0] if se_res else float('nan')
t_stat = beta / se if se and se > 0 else float('nan')

random.seed(SEED)
betas = []
cluster_ids = list(set(clusters))
cluster_map = defaultdict(list)
for i, cid in enumerate(clusters):
    cluster_map[cid].append(i)

for _ in range(B):
    # Rademacher weights per cluster
    weight = {cid: (1 if random.random() < 0.5 else -1) for cid in cluster_ids}
    y_star = [beta * x[i] + resid[i] * weight[clusters[i]] for i in range(len(x))]
    b_star = sum(x[i] * y_star[i] for i in range(len(x))) / x2 if x2 > 0 else 0.0
    betas.append(b_star)

abs_b = abs(beta)
p_val = sum(1 for b in betas if abs(b) >= abs_b) / len(betas) if betas else float('nan')

out_path = results_path('metrics', f'strategy2_panel_did_wild_bootstrap_tier1_{TRANSFORM_SET}.md')
lines = [
    '# Wild Cluster Bootstrap (event clusters)',
    '',
    f'- B: {B}',
    f'- transform: {TRANSFORM_SET}',
    f'- events: {len(events)} (tier_filter={TIER_FILTER}, binding_only={BINDING_ONLY})',
    f'- N: {len(obs_dm)}',
    '',
    f'- beta: {beta:+.4f}',
    f'- se(cluster): {se:.4f}',
    f'- t: {t_stat:.2f}',
    f'- wild p-value: {p_val:.4f}',
    '',
]
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print(f'Wrote: {out_path}')
