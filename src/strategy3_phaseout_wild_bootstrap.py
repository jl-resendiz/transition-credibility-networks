"""Wild cluster bootstrap for the coal phase-out DiD specification.

The paper "Pricing Transition Credibility on Spatial Networks" has only
14 Tier-1 binding phase-out events, making standard cluster-robust
standard errors unreliable when clustering at the event level. This
script implements a wild cluster bootstrap with Webb 6-point weights
to obtain valid inference under few clusters.

Specification:
    AR_{j,t} = firm_FE + month_FE + beta * (coal_share_j * Post_t) + eps

where coal_share_j is the firm-level coal capacity share from the alpha
panel, Post_t = 1 for months [0, +12] relative to the event, and the
sample includes treated firms (matched to the phase-out jurisdiction)
plus a control pool. Firm and month fixed effects are removed via
iterative demeaning.

Bootstrap procedure (Cameron, Gelbach & Miller 2008):
    - Cluster at the event level (G = number of events)
    - B = 999 replications using Webb 6-point weights
    - Under H0 (beta=0): y* = w_g * residual_i  (restricted residuals)
    - Bootstrap t-statistic: t*_b = beta*_b / se*_b
    - p-value = (1 + #{|t*_b| >= |t_obs|}) / (1 + B)

Outputs:
    results/metrics/strategy3_phaseout_wild_bootstrap.md
"""
import csv
import math
import os
import random
import hashlib
from collections import defaultdict

from _paths import raw_path, derived_path, results_path

# ----------------------------- Configuration -----------------------------

TIER_FILTER = os.getenv('TIER_FILTER', '1')
BINDING_ONLY = os.getenv('BINDING_ONLY', '1') == '1'
TAU_START = int(os.getenv('TAU_START', '-6'))
TAU_END = int(os.getenv('TAU_END', '12'))
POST_START = int(os.getenv('POST_START', '0'))
POST_END = int(os.getenv('POST_END', '12'))
CONTROL_MULT = int(os.getenv('CONTROL_MULT', '5'))
OVERLAP_RULE = os.getenv('OVERLAP_RULE', 'nearest')
B = int(os.getenv('B', '999'))
SEED = int(os.getenv('SEED', '42'))

# Webb 6-point weights
WEBB_WEIGHTS = [
    -math.sqrt(3.0 / 2.0),
    -math.sqrt(2.0 / 2.0),
    -math.sqrt(1.0 / 2.0),
    math.sqrt(1.0 / 2.0),
    math.sqrt(2.0 / 2.0),
    math.sqrt(3.0 / 2.0),
]


# ----------------------------- Helpers -----------------------------

def add_months(ym, delta):
    y, m = ym.split('-')
    y, m = int(y), int(m)
    m += delta
    while m > 12:
        y += 1
        m -= 12
    while m < 1:
        y -= 1
        m += 12
    return f"{y:04d}-{m:02d}"


def load_ff_factors_monthly(path):
    """Parse Fama-French monthly factors to get vwretd = Mkt-RF + RF."""
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


def load_coal_share(panel_path):
    """Load coal_share by gvkey-year from firm_alpha_panel.csv.

    Returns a function get_share(gvkey, year) that finds the nearest
    year <= the requested year, or the earliest available year.
    """
    coal_by_year = defaultdict(dict)
    years_by_gvkey = defaultdict(list)
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


def two_way_demean(data, y_var, x_var, firm_key='gvkey', time_key='ym'):
    """Remove firm and month fixed effects via additive demeaning."""
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
    """One-way clustered standard error for univariate OLS."""
    n = len(x)
    if n == 0:
        return None
    x2 = sum(v * v for v in x)
    if x2 <= 1e-12:
        return None
    clus = defaultdict(list)
    for i, cid in enumerate(clusters):
        clus[cid].append(i)
    G = len(clus)
    S = 0.0
    for _, idxs in clus.items():
        xu = sum(x[i] * resid[i] for i in idxs)
        S += xu * xu
    # Small-sample correction: G/(G-1) * (N-1)/(N-1) ~ G/(G-1)
    if G > 1:
        scale = G / (G - 1.0)
        S *= scale
    se = math.sqrt(S / (x2 * x2))
    return se, G


# ----------------------------- Load data -----------------------------

print('Loading monthly returns ...')
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except (ValueError, TypeError):
            pass
print(f'  {len(monthly_ret)} firms loaded')

print('Loading Fama-French factors ...')
vwretd = load_ff_factors_monthly(raw_path('factors', 'F-F_Research_Data_Factors.csv'))
if not vwretd:
    raise RuntimeError('Missing F-F monthly factors for vwretd.')
print(f'  {len(vwretd)} month factors loaded')

print('Loading geographic weight matrix ...')
W = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
print(f'  {len(W)} firms in weight matrix')

print('Loading firm fundamentals ...')
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row
print(f'  {len(fundamentals)} firms loaded')

print('Loading coal share from alpha panel ...')
get_coal_share = load_coal_share(derived_path('fundamentals', 'firm_alpha_panel.csv'))

# ----------------------------- Load events -----------------------------

print('Loading phase-out events ...')
events_path = derived_path('events', 'coal_phaseout_shocks_events.csv')
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
        event_date = row.get('event_date', '').strip()
        if not event_date:
            continue
        event_month = event_date[:7]
        event_year = int(event_date[:4]) if event_date[:4].isdigit() else None
        events.append({
            'event_month': event_month,
            'event_year': event_year,
            'gvkeys': [g for g in row['matched_gvkeys'].split(';') if g],
            'shock_name': row.get('shock_name', ''),
        })

print(f'  {len(events)} Tier-1 binding events')

# ----------------------------- Build panel -----------------------------

print('Building panel observations ...')
obs = []
for event_id, event in enumerate(events):
    event_month = event['event_month']
    event_year = event['event_year']
    if not event_month or event_year is None:
        continue
    treated = set(event['gvkeys'])
    # Control pool: all firms not treated in this event
    non_treated = [gk for gk in fundamentals if gk not in treated]
    stable_seed = int(hashlib.md5(str(event_id).encode('utf-8')).hexdigest()[:8], 16)
    random.seed(stable_seed)
    n_ctrl = min(len(non_treated), max(CONTROL_MULT * len(treated), 20))
    ctrl_sample = random.sample(non_treated, n_ctrl) if len(non_treated) > n_ctrl else non_treated
    candidate_firms = list(treated) + ctrl_sample

    for gk in candidate_firms:
        cs = get_coal_share(gk, event_year)
        if cs is None:
            continue
        # Exposure: coal_share for treated firms, 0 for controls
        exp = cs if gk in treated else 0.0
        if gk not in monthly_ret:
            continue
        for tau in range(TAU_START, TAU_END + 1):
            ym = add_months(event_month, tau)
            if ym not in monthly_ret[gk] or ym not in vwretd:
                continue
            ar = monthly_ret[gk][ym] - vwretd[ym]
            post = 1.0 if (tau >= POST_START and tau <= POST_END) else 0.0
            obs.append({
                'gvkey': gk,
                'ym': ym,
                'event_id': event_id,
                'ar': ar,
                'exp_post': exp * post,
            })

print(f'  Raw observations: {len(obs)}')

# ----------------------------- Overlap handling -----------------------------

if OVERLAP_RULE == 'nearest':
    grouped = defaultdict(list)
    for row in obs:
        key = (row['gvkey'], row['ym'])
        grouped[key].append(row)
    obs_clean = []
    overlaps = 0
    for _, rows in grouped.items():
        if len(rows) == 1:
            obs_clean.append(rows[0])
            continue
        overlaps += 1
        # Keep the observation with highest absolute exposure*post
        rows_sorted = sorted(rows, key=lambda r: abs(r['exp_post']), reverse=True)
        obs_clean.append(rows_sorted[0])
    obs = obs_clean
    print(f'  Overlap groups: {overlaps}, kept obs: {len(obs)}')

# ----------------------------- Two-way demean -----------------------------

print('Demeaning (firm + month FE) ...')
obs_dm = two_way_demean(obs, 'ar', 'exp_post', firm_key='gvkey', time_key='ym')
x = [row['exp_post'] for row in obs_dm]
y = [row['ar'] for row in obs_dm]
clusters = [row['event_id'] for row in obs_dm]
N = len(x)
print(f'  N after demeaning: {N}')

# ----------------------------- OLS point estimate -----------------------------

x2 = sum(v * v for v in x)
if x2 <= 1e-12:
    raise RuntimeError('No variation in treatment variable after demeaning.')

beta_hat = sum(x[i] * y[i] for i in range(N)) / x2
resid = [y[i] - beta_hat * x[i] for i in range(N)]

se_res = cluster_se(x, resid, clusters)
if se_res is None:
    raise RuntimeError('Cannot compute clustered SE.')
se_hat, G = se_res
t_obs = beta_hat / se_hat if se_hat > 1e-15 else float('nan')

print()
print('=== Observed estimates ===')
print(f'  beta(exp_post)  = {beta_hat:+.6f}')
print(f'  se(cluster)     = {se_hat:.6f}')
print(f'  t-stat          = {t_obs:.3f}')
print(f'  clusters (G)    = {G}')
print(f'  N               = {N}')

# ----------------------------- Wild cluster bootstrap -----------------------------

print()
print(f'Running wild cluster bootstrap (B={B}, Webb 6-point weights) ...')
random.seed(SEED)

# Build cluster index: event_id -> list of row indices
cluster_ids = sorted(set(clusters))
cluster_map = defaultdict(list)
for i, cid in enumerate(clusters):
    cluster_map[cid].append(i)

# Restricted residuals under H0: beta = 0 => y* uses full y as residual
# Under H0 the restricted model is: y = FE + eps, so restricted residual = y (demeaned)
resid_h0 = y[:]

n_exceed = 0
boot_t_stats = []

for b in range(B):
    if (b + 1) % 200 == 0:
        print(f'  replication {b + 1}/{B} ...')

    # Draw Webb weight per cluster
    weights = {cid: random.choice(WEBB_WEIGHTS) for cid in cluster_ids}

    # Construct bootstrap dependent variable: y*_i = w_g * resid_h0_i
    # (imposing H0: beta=0, so no X*beta_restricted term)
    y_star = [resid_h0[i] * weights[clusters[i]] for i in range(N)]

    # Bootstrap OLS: beta*
    beta_star = sum(x[i] * y_star[i] for i in range(N)) / x2

    # Bootstrap residuals and clustered SE
    resid_star = [y_star[i] - beta_star * x[i] for i in range(N)]
    se_star_res = cluster_se(x, resid_star, clusters)
    if se_star_res is None or se_star_res[0] <= 1e-15:
        # Degenerate draw: skip (conservative)
        boot_t_stats.append(0.0)
        continue
    se_star = se_star_res[0]
    t_star = beta_star / se_star

    boot_t_stats.append(t_star)
    if abs(t_star) >= abs(t_obs):
        n_exceed += 1

# p-value = (1 + #{|t*| >= |t_obs|}) / (1 + B)
p_val = (1.0 + n_exceed) / (1.0 + B)

print()
print('=== Wild Cluster Bootstrap Results ===')
print(f'  B replications   = {B}')
print(f'  Webb 6-pt weights')
print(f'  Clusters (G)     = {G}')
print(f'  N                = {N}')
print(f'  Observed beta    = {beta_hat:+.6f}')
print(f'  Observed t-stat  = {t_obs:.3f}')
print(f'  Bootstrap p-val  = {p_val:.4f}')
print(f'  |t*| >= |t_obs|  = {n_exceed} / {B}')

# Percentiles of bootstrap t distribution
boot_t_sorted = sorted(boot_t_stats)
def percentile(sorted_vals, p):
    idx = p / 100.0 * (len(sorted_vals) - 1)
    lo = int(math.floor(idx))
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = idx - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac

pcts = [2.5, 5.0, 50.0, 95.0, 97.5]
print()
print('  Bootstrap t-stat distribution:')
for p in pcts:
    print(f'    {p:5.1f}th pctile = {percentile(boot_t_sorted, p):+.3f}')

# ----------------------------- Write output -----------------------------

out_path = results_path('metrics', 'strategy3_phaseout_wild_bootstrap.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

lines = [
    '# Wild Cluster Bootstrap: Coal Phase-Out DiD',
    '',
    '## Specification',
    '',
    'AR_{j,t} = firm_FE + month_FE + beta * (coal_share_j * Post_t) + eps',
    '',
    '- Exposure: coal_share from firm_alpha_panel (treated firms only; controls = 0)',
    '- Post: months [0, +12] relative to event',
    f'- Event window: tau in [{TAU_START}, {TAU_END}]',
    f'- Overlap rule: {OVERLAP_RULE}',
    f'- Tier filter: {TIER_FILTER} (binding only: {BINDING_ONLY})',
    '',
    '## Observed estimates',
    '',
    f'- beta(exp_post): {beta_hat:+.6f}',
    f'- se(cluster): {se_hat:.6f}',
    f'- t-stat: {t_obs:.3f}',
    f'- N: {N}',
    f'- Clusters (G): {G}',
    '',
    '## Wild cluster bootstrap',
    '',
    f'- B: {B}',
    '- Weight distribution: Webb 6-point',
    f'- Seed: {SEED}',
    f'- Bootstrap p-value: {p_val:.4f}',
    f'- |t*| >= |t_obs|: {n_exceed} / {B}',
    '',
    '### Bootstrap t-stat distribution',
    '',
]
for p in pcts:
    lines.append(f'- {p:.1f}th percentile: {percentile(boot_t_sorted, p):+.4f}')
lines.append('')

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

print()
print(f'Wrote: {out_path}')
