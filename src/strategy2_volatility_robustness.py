"""Volatility robustness: alternative volatility measures.

Tests three volatility measures around event windows:
  (a) SD-based (baseline): vol_change = SD(AR post) - SD(AR pre)
  (b) |AR|-based (Beaver 1968): mean(|AR| post) / mean(|AR| pre) - 1
  (c) Squared-return-based: mean(AR^2 post) / mean(AR^2 pre) - 1

Runs regressions for each measure:
  - EIA-860 events: vol_measure ~ geographic exposure (W)
  - Phase-out events (Tier-1 binding): vol_measure ~ coal_share
  - CAR ~ exposure / coal_share as comparison

Outputs: results/metrics/strategy2_volatility_robustness.md
"""
import csv
import math
import os
import bisect
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

# ── Windows ──────────────────────────────────────────────────────────
PRE_START = -21
PRE_END = -1
POST_START = 1
POST_END = 20
CAR_START = -1
CAR_END = 20
MIN_COVERAGE = 0.4

CONTROL_MULT = 5


# ── Helpers ──────────────────────────────────────────────────────────

def load_ff_factors_daily(path):
    """Load Fama-French daily factors and return vwretd dict {date: vw}."""
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
            vwretd[f'{date[:4]}-{date[4:6]}-{date[6:]}'] = vw
    return vwretd if vwretd else None


def load_daily_returns(path):
    """Load daily returns into {gvkey: {date: ret}}."""
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


def load_coal_share(panel_path):
    """Load coal_share from alpha panel, return lookup function."""
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


def window_ar_values(dates, ar, event_idx, start, end):
    """Extract AR values in [event_idx + start, event_idx + end].

    Returns list of values if coverage >= MIN_COVERAGE, else None.
    """
    required = end - start + 1
    vals = []
    for offset in range(start, end + 1):
        j = event_idx + offset
        if 0 <= j < len(dates):
            vals.append(ar[j])
    if len(vals) < required * MIN_COVERAGE:
        return None
    return vals


def compute_vol_measures(pre_vals, post_vals):
    """Compute three volatility change measures from pre/post AR lists.

    Returns dict with keys: vol_sd, vol_abs, vol_sq (or None values if
    denominators are zero).
    """
    # SD-based: SD(post) - SD(pre)
    def sd(vals):
        n = len(vals)
        if n < 2:
            return None
        mean = sum(vals) / n
        var = sum((v - mean) ** 2 for v in vals) / (n - 1)
        return math.sqrt(var)

    sd_pre = sd(pre_vals)
    sd_post = sd(post_vals)
    vol_sd = (sd_post - sd_pre) if (sd_pre is not None and sd_post is not None) else None

    # |AR|-based (Beaver 1968): mean(|AR| post) / mean(|AR| pre) - 1
    abs_pre = sum(abs(v) for v in pre_vals) / len(pre_vals)
    abs_post = sum(abs(v) for v in post_vals) / len(post_vals)
    vol_abs = (abs_post / abs_pre - 1.0) if abs_pre > 1e-15 else None

    # Squared-return-based: mean(AR^2 post) / mean(AR^2 pre) - 1
    sq_pre = sum(v * v for v in pre_vals) / len(pre_vals)
    sq_post = sum(v * v for v in post_vals) / len(post_vals)
    vol_sq = (sq_post / sq_pre - 1.0) if sq_pre > 1e-15 else None

    return {'vol_sd': vol_sd, 'vol_abs': vol_abs, 'vol_sq': vol_sq}


def ols_cluster(y, x, clusters):
    """OLS with event-clustered SEs. x is list of [1, exposure].

    Returns (beta1, se1, t1, n_clusters, n) or None.
    """
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
        [inv[0][0] * S[0][0] * inv[0][0] + inv[0][1] * S[1][0] * inv[0][0]
         + inv[0][0] * S[0][1] * inv[0][1] + inv[0][1] * S[1][1] * inv[0][1],
         inv[0][0] * S[0][0] * inv[1][0] + inv[0][1] * S[1][0] * inv[1][0]
         + inv[0][0] * S[0][1] * inv[1][1] + inv[0][1] * S[1][1] * inv[1][1]],
        [inv[1][0] * S[0][0] * inv[0][0] + inv[1][1] * S[1][0] * inv[0][0]
         + inv[1][0] * S[0][1] * inv[0][1] + inv[1][1] * S[1][1] * inv[0][1],
         inv[1][0] * S[0][0] * inv[1][0] + inv[1][1] * S[1][0] * inv[1][0]
         + inv[1][0] * S[0][1] * inv[1][1] + inv[1][1] * S[1][1] * inv[1][1]]
    ]
    se1 = math.sqrt(cov[1][1]) if cov[1][1] > 0 else float('nan')
    t1 = beta1 / se1 if se1 and se1 > 0 else float('nan')
    return beta1, se1, t1, len(clus), n


# ── Load data ────────────────────────────────────────────────────────

print('Loading Fama-French daily factors...')
ff_path = raw_path('factors', 'F-F_Research_Data_Factors_daily.csv')
market_ret = load_ff_factors_daily(ff_path)
if not market_ret:
    raise RuntimeError('Missing F-F daily factors for vwretd.')
print(f'  {len(market_ret)} trading days loaded.')

print('Loading daily returns...')
daily_ret = load_daily_returns(derived_path('returns', 'daily_returns.csv'))
print(f'  {len(daily_ret)} firms loaded.')

print('Building daily AR series...')
daily_ar = {}
for gk, dct in daily_ret.items():
    dates = [d for d in dct.keys() if d in market_ret]
    if not dates:
        continue
    dates.sort()
    ar = [dct[d] - market_ret[d] for d in dates]
    daily_ar[gk] = (dates, ar)
print(f'  {len(daily_ar)} firms with AR series.')

print('Loading geographic weight matrix...')
W = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
print(f'  {len(W)} source firms in W.')

print('Loading firm fundamentals...')
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row
print(f'  {len(fundamentals)} firms loaded.')

print('Loading coal share from alpha panel...')
get_coal_share = load_coal_share(derived_path('fundamentals', 'firm_alpha_panel.csv'))
print('  Done.')


# ── EIA-860 event observations ──────────────────────────────────────

print('\n=== EIA-860 Retirement Announcements ===')

eia_events = []
with open(derived_path('events', 'eia860_announcement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        ann_date = row.get('announcement_date', '').strip()
        event_date = row.get('event_date', '').strip()
        if not ann_date and not event_date:
            continue
        effective_date = ann_date if ann_date else event_date
        if len(effective_date) == 7 and effective_date[4] == '-':
            effective_date = f"{effective_date}-15"
        eia_events.append({
            'event_date': effective_date,
            'gvkeys': [g for g in row['matched_gvkeys'].split(';') if g],
        })

print(f'  {len(eia_events)} EIA-860 events loaded.')

eia_obs = []
for event_id, event in enumerate(eia_events):
    event_date = event['event_date']
    event_gvkeys = set(event['gvkeys'])
    for fm_gk in event_gvkeys:
        if fm_gk not in W:
            continue
        neighbors = W[fm_gk]
        neighbor_gks = set(neighbors.keys()) - event_gvkeys
        non_connected = [gk for gk in fundamentals
                         if gk not in event_gvkeys and gk not in neighbors]
        n_ctrl = min(len(non_connected), max(CONTROL_MULT * len(neighbor_gks), 20))
        ctrl_sample = non_connected[:n_ctrl]
        candidate_firms = list(neighbor_gks) + ctrl_sample

        for gk in candidate_firms:
            w_ij = neighbors.get(gk, 0.0)
            if gk not in daily_ar:
                continue
            dates, ar = daily_ar[gk]
            # locate event index
            event_idx = bisect.bisect_left(dates, event_date)
            if event_idx >= len(dates):
                continue

            pre_vals = window_ar_values(dates, ar, event_idx, PRE_START, PRE_END)
            post_vals = window_ar_values(dates, ar, event_idx, POST_START, POST_END)
            car_vals = window_ar_values(dates, ar, event_idx, CAR_START, CAR_END)
            if pre_vals is None or post_vals is None or car_vals is None:
                continue

            measures = compute_vol_measures(pre_vals, post_vals)
            car_sum = sum(car_vals)

            eia_obs.append({
                'event_id': event_id,
                'gvkey': gk,
                'exposure': w_ij,
                'vol_sd': measures['vol_sd'],
                'vol_abs': measures['vol_abs'],
                'vol_sq': measures['vol_sq'],
                'car': car_sum,
            })

print(f'  {len(eia_obs)} event-firm observations for EIA-860.')


# ── Phase-out event observations (Tier-1 binding) ───────────────────

print('\n=== Binding Phase-Out Laws (Tier-1) ===')

phaseout_events = []
with open(derived_path('events', 'coal_phaseout_shocks_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        if row.get('exogeneity_tier', '').strip() != '1':
            continue
        if row.get('binding', '').strip().lower() != 'yes':
            continue
        event_date = row.get('event_date', '').strip()
        if not event_date or len(event_date) < 10:
            continue
        event_year = int(event_date[:4]) if event_date[:4].isdigit() else None
        phaseout_events.append({
            'event_date': event_date,
            'event_year': event_year,
            'gvkeys': [g for g in row['matched_gvkeys'].split(';') if g],
        })

print(f'  {len(phaseout_events)} Tier-1 binding phase-out events loaded.')

phaseout_obs = []
for event_id, event in enumerate(phaseout_events):
    event_date = event['event_date']
    event_year = event['event_year']
    if event_year is None:
        continue
    event_gvkeys = set(event['gvkeys'])
    # All matched firms are treated; add controls from fundamentals pool
    non_treated = [gk for gk in fundamentals if gk not in event_gvkeys]
    n_ctrl = min(len(non_treated), max(CONTROL_MULT * len(event_gvkeys), 20))
    ctrl_sample = non_treated[:n_ctrl]
    candidate_firms = list(event_gvkeys) + ctrl_sample

    for gk in candidate_firms:
        cs = get_coal_share(gk, event_year)
        if cs is None:
            continue
        # Exposure is coal_share for treated firms, 0 for controls
        exposure = cs if gk in event_gvkeys else 0.0
        if gk not in daily_ar:
            continue
        dates, ar = daily_ar[gk]
        event_idx = bisect.bisect_left(dates, event_date)
        if event_idx >= len(dates):
            continue

        pre_vals = window_ar_values(dates, ar, event_idx, PRE_START, PRE_END)
        post_vals = window_ar_values(dates, ar, event_idx, POST_START, POST_END)
        car_vals = window_ar_values(dates, ar, event_idx, CAR_START, CAR_END)
        if pre_vals is None or post_vals is None or car_vals is None:
            continue

        measures = compute_vol_measures(pre_vals, post_vals)
        car_sum = sum(car_vals)

        phaseout_obs.append({
            'event_id': event_id,
            'gvkey': gk,
            'exposure': exposure,
            'vol_sd': measures['vol_sd'],
            'vol_abs': measures['vol_abs'],
            'vol_sq': measures['vol_sq'],
            'car': car_sum,
        })

print(f'  {len(phaseout_obs)} event-firm observations for phase-out.')


# ── Regressions ─────────────────────────────────────────────────────

def run_regressions(obs_list, label):
    """Run regressions for all four measures, return list of result dicts."""
    measures = [
        ('SD-based (baseline)', 'vol_sd'),
        ('|AR|-based (Beaver)', 'vol_abs'),
        ('Squared-return', 'vol_sq'),
        ('CAR (level)', 'car'),
    ]
    results = []
    for measure_name, key in measures:
        # Filter to obs where the measure is not None
        valid = [o for o in obs_list if o[key] is not None]
        if not valid:
            results.append({
                'name': measure_name,
                'coeff': None, 'se': None, 't': None, 'n': 0,
            })
            continue
        y = [o[key] for o in valid]
        X = [[1.0, o['exposure']] for o in valid]
        clusters = [o['event_id'] for o in valid]
        res = ols_cluster(y, X, clusters)
        if res is None:
            results.append({
                'name': measure_name,
                'coeff': None, 'se': None, 't': None, 'n': len(valid),
            })
            continue
        beta1, se1, t1, n_clus, n_obs = res
        results.append({
            'name': measure_name,
            'coeff': beta1, 'se': se1, 't': t1, 'n': n_obs,
        })
        print(f'  {label} | {measure_name}: beta={beta1:+.6f}, se={se1:.6f}, '
              f't={t1:.2f}, N={n_obs}')
    return results


print('\nRunning EIA-860 regressions (vol_measure ~ geographic exposure)...')
eia_results = run_regressions(eia_obs, 'EIA-860')

print('\nRunning phase-out regressions (vol_measure ~ coal_share)...')
phaseout_results = run_regressions(phaseout_obs, 'Phase-out')


# ── Write output ─────────────────────────────────────────────────────

def fmt(val, decimals=6):
    if val is None:
        return '---'
    return f'{val:+.{decimals}f}' if decimals > 0 else f'{val}'


def fmt_unsigned(val, decimals=6):
    if val is None:
        return '---'
    return f'{val:.{decimals}f}'


def fmt_t(val):
    if val is None:
        return '---'
    return f'{val:.2f}'


def fmt_n(val):
    if val is None or val == 0:
        return '---'
    return str(val)


lines = [
    '# Volatility Robustness: Alternative Measures',
    '',
    f'Pre-event window: [{PRE_START}, {PRE_END}] trading days.',
    f'Post-event window: [{POST_START}, {POST_END}] trading days.',
    f'CAR window: [{CAR_START}, {CAR_END}] trading days.',
    f'Minimum data coverage: {MIN_COVERAGE:.0%} in each window.',
    '',
    '## EIA-860 Retirement Announcements',
    '',
    '| Measure | Coeff | SE | t-stat | N |',
    '|---|---|---|---|---|',
]

for r in eia_results:
    lines.append(f'| {r["name"]} | {fmt(r["coeff"])} | {fmt_unsigned(r["se"])} '
                 f'| {fmt_t(r["t"])} | {fmt_n(r["n"])} |')

lines += [
    '',
    '## Binding Phase-Out Laws (Tier-1)',
    '',
    '| Measure | Coeff | SE | t-stat | N |',
    '|---|---|---|---|---|',
]

for r in phaseout_results:
    lines.append(f'| {r["name"]} | {fmt(r["coeff"])} | {fmt_unsigned(r["se"])} '
                 f'| {fmt_t(r["t"])} | {fmt_n(r["n"])} |')

lines += [
    '',
    '## Notes',
    '',
    'All regressions use event-clustered standard errors.',
    '',
    'The |AR|-based measure follows Beaver (1968) and is robust to '
    'GARCH-type volatility clustering. The squared-return measure captures '
    'variance changes without assuming normality.',
    '',
    'SD-based: vol_change = SD(AR post) - SD(AR pre).',
    '|AR|-based: vol_change = mean(|AR| post) / mean(|AR| pre) - 1.',
    'Squared-return: vol_change = mean(AR^2 post) / mean(AR^2 pre) - 1.',
    'CAR: cumulative abnormal return over [-1, +20] window.',
    '',
    'EIA-860 regressions: outcome ~ geographic exposure (W_ij).',
    'Phase-out regressions: outcome ~ coal_share (treated) vs 0 (controls).',
    '',
]

out_path = results_path('metrics', 'strategy2_volatility_robustness.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

print(f'\nWrote: {out_path}')
