"""Fisher randomization inference for the coal phase-out DiD specification.

The paper "Pricing Transition Credibility on Spatial Networks" has only
14 Tier-1 binding phase-out events, making standard asymptotic inference
questionable. This script implements Fisher randomization inference as a
non-parametric alternative.

Under the sharp null H0: no treatment effect for any unit at any time,
the assignment of events to dates is exchangeable. We permute the event
dates across the 14 events, recompute the full panel (post indicators,
exp_post), re-demean, re-estimate OLS, and collect the permutation
t-statistics.

Specification:
    AR_{j,t} = firm_FE + month_FE + beta * (coal_share_j * Post_t) + eps

where coal_share_j is the firm-level coal capacity share from the alpha
panel, Post_t = 1 for months [0, +12] relative to the event, and the
sample includes treated firms (matched to the phase-out jurisdiction)
plus a control pool. Firm and month fixed effects are removed via
iterative demeaning.

Fisher p-value = (1 + #{|t*_b| >= |t_obs|}) / (1 + B)

Outputs:
    results/metrics/strategy3_fisher_randomization.md
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


def two_way_demean(data, y_var, x_vars, fe1='gvkey', fe2='ym', max_iter=100, tol=1e-8):
    """Remove firm and month fixed effects via iterative demeaning."""
    vars_to_demean = [y_var] + x_vars
    # Work on copies to avoid mutating originals
    working = []
    for d in data:
        row = {fe1: d[fe1], fe2: d[fe2], 'event_id': d.get('event_id')}
        for v in vars_to_demean:
            row[v] = d[v]
        working.append(row)

    for iteration in range(max_iter):
        max_change = 0.0
        for var in vars_to_demean:
            # Demean by fe1
            means = {}
            counts = {}
            for d in working:
                k = d[fe1]
                means[k] = means.get(k, 0.0) + d[var]
                counts[k] = counts.get(k, 0) + 1
            for k in means:
                means[k] /= counts[k]
            for d in working:
                old = d[var]
                d[var] -= means[d[fe1]]
                max_change = max(max_change, abs(d[var] - old))
            # Demean by fe2
            means = {}
            counts = {}
            for d in working:
                k = d[fe2]
                means[k] = means.get(k, 0.0) + d[var]
                counts[k] = counts.get(k, 0) + 1
            for k in means:
                means[k] /= counts[k]
            for d in working:
                old = d[var]
                d[var] -= means[d[fe2]]
                max_change = max(max_change, abs(d[var] - old))
        if max_change < tol:
            break
    return working


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
    # Small-sample correction: G/(G-1)
    if G > 1:
        scale = G / (G - 1.0)
        S *= scale
    se = math.sqrt(S / (x2 * x2))
    return se, G


def ols_estimate(obs_list):
    """Demean, run OLS, return (beta, se, t, N, G) or None if degenerate."""
    obs_dm = two_way_demean(obs_list, 'ar', ['exp_post'],
                            fe1='gvkey', fe2='ym')
    x = [row['exp_post'] for row in obs_dm]
    y = [row['ar'] for row in obs_dm]
    clusters = [row['event_id'] for row in obs_dm]
    N = len(x)
    if N == 0:
        return None

    x2 = sum(v * v for v in x)
    if x2 <= 1e-12:
        return None

    beta = sum(x[i] * y[i] for i in range(N)) / x2
    resid = [y[i] - beta * x[i] for i in range(N)]

    se_res = cluster_se(x, resid, clusters)
    if se_res is None:
        return None
    se, G = se_res
    t = beta / se if se > 1e-15 else float('nan')
    return beta, se, t, N, G


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


# ----------------------------- Panel builder -----------------------------

def build_panel(events_list):
    """Build panel observations for a given list of events.

    Each event carries its own event_month. This function is called once
    for the observed data and once per permutation (with shuffled dates).
    Returns the list of observation dicts.
    """
    obs = []
    for event_id, event in enumerate(events_list):
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
    return obs


def handle_overlaps(obs):
    """Resolve overlapping observations using the nearest rule."""
    if OVERLAP_RULE != 'nearest':
        return obs
    grouped = defaultdict(list)
    for row in obs:
        key = (row['gvkey'], row['ym'])
        grouped[key].append(row)
    obs_clean = []
    for _, rows in grouped.items():
        if len(rows) == 1:
            obs_clean.append(rows[0])
            continue
        rows_sorted = sorted(rows, key=lambda r: abs(r['exp_post']), reverse=True)
        obs_clean.append(rows_sorted[0])
    return obs_clean


# ----------------------------- Build observed panel -----------------------------

print('Building observed panel ...')
obs = build_panel(events)
obs = handle_overlaps(obs)
print(f'  Observations after overlap handling: {len(obs)}')

# ----------------------------- Observed OLS -----------------------------

print('Estimating observed specification ...')
result_obs = ols_estimate(obs)
if result_obs is None:
    raise RuntimeError('Cannot estimate observed specification (degenerate).')
beta_hat, se_hat, t_obs, N, G = result_obs

print()
print('=== Observed estimates ===')
print(f'  beta(exp_post)  = {beta_hat:+.6f}')
print(f'  se(cluster)     = {se_hat:.6f}')
print(f'  t-stat          = {t_obs:.3f}')
print(f'  clusters (G)    = {G}')
print(f'  N               = {N}')

# ----------------------------- Fisher randomization -----------------------------

print()
print(f'Running Fisher randomization inference (B={B} permutations) ...')
random.seed(SEED)

# Extract the observed event months for permutation
observed_event_months = [e['event_month'] for e in events]
observed_event_years = [e['event_year'] for e in events]
n_events = len(events)

n_exceed = 0
perm_t_stats = []

for b in range(B):
    if (b + 1) % 50 == 0:
        print(f'  permutation {b + 1}/{B} ...')

    # Permute event dates across events
    perm_indices = list(range(n_events))
    random.shuffle(perm_indices)

    # Build permuted event list: same firms/gvkeys per event, but shuffled dates
    events_perm = []
    for i in range(n_events):
        events_perm.append({
            'event_month': observed_event_months[perm_indices[i]],
            'event_year': observed_event_years[perm_indices[i]],
            'gvkeys': events[i]['gvkeys'],
            'shock_name': events[i]['shock_name'],
        })

    # Build the full panel with permuted dates
    obs_perm = build_panel(events_perm)
    obs_perm = handle_overlaps(obs_perm)

    # Estimate OLS with demeaning
    result_perm = ols_estimate(obs_perm)
    if result_perm is None:
        perm_t_stats.append(0.0)
        continue

    _, _, t_perm, _, _ = result_perm
    perm_t_stats.append(t_perm)

    if abs(t_perm) >= abs(t_obs):
        n_exceed += 1

# Fisher p-value = (1 + #{|t*_b| >= |t_obs|}) / (1 + B)
fisher_p = (1.0 + n_exceed) / (1.0 + B)

print()
print('=== Fisher Randomization Results ===')
print(f'  B permutations   = {B}')
print(f'  Clusters (G)     = {G}')
print(f'  N                = {N}')
print(f'  Observed beta    = {beta_hat:+.6f}')
print(f'  Observed t-stat  = {t_obs:.3f}')
print(f'  Fisher p-value   = {fisher_p:.4f}')
print(f'  |t*| >= |t_obs|  = {n_exceed} / {B}')

# Percentiles of permutation t distribution
perm_t_sorted = sorted(perm_t_stats)
def percentile(sorted_vals, p):
    idx = p / 100.0 * (len(sorted_vals) - 1)
    lo = int(math.floor(idx))
    hi = min(lo + 1, len(sorted_vals) - 1)
    frac = idx - lo
    return sorted_vals[lo] * (1 - frac) + sorted_vals[hi] * frac

pcts = [2.5, 5.0, 50.0, 95.0, 97.5]
print()
print('  Permutation t-stat distribution:')
for p in pcts:
    print(f'    {p:5.1f}th pctile = {percentile(perm_t_sorted, p):+.3f}')

# ----------------------------- Write output -----------------------------

out_path = results_path('metrics', 'strategy3_fisher_randomization.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

lines = [
    '# Fisher Randomization Inference: Phase-Out DiD',
    '',
    f'Observed coefficient: {beta_hat:+.6f}',
    f'Observed t-statistic: {t_obs:.3f}',
    f'Fisher p-value: {fisher_p:.4f}',
    f'Permutations: {B}',
    f'Tier-1 binding events: {n_events}',
    f'N: {N}',
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
    '## Permutation t-stat distribution',
    '',
]
for p in pcts:
    lines.append(f'- {p:.1f}th percentile: {percentile(perm_t_sorted, p):+.4f}')
lines += [
    '',
    f'|t*| >= |t_obs|: {n_exceed} / {B}',
    '',
    '## Notes',
    '',
    'Under H0 of no treatment effect, event dates are permuted',
    f'across the {n_events} Tier-1 binding phase-out events. Firm + month FE',
    'absorbed via iterative demeaning. Event-clustered SEs.',
    f'Seed: {SEED}.',
    '',
]

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

print()
print(f'Wrote: {out_path}')
