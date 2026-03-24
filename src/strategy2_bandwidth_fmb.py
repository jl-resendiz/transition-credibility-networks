"""Bandwidth sensitivity for W_geo under Fama-MacBeth + Newey-West inference.

Rebuilds W_geo at multiple bandwidths (250, 500, 750, 1000, 1500 km half-life)
and runs the full Fama-MacBeth procedure event-by-event for each bandwidth.
Reports Newey-West HAC SEs on the beta time series.

Economic motivation: the competitive benefit of coal retirement operates through
interconnected transmission grids (ENTSO-E ~1500km, US ISOs ~1000km), not at
plant-level proximity. The 500km default half-life may be too narrow.

Output: results/metrics/strategy2_bandwidth_fmb.md
"""
import csv
import os
import sys
import math
import random
import hashlib
from collections import defaultdict

from _paths import derived_path, raw_path, results_path


def _print(msg=''):
    print(msg)
    sys.stdout.flush()


# ── Matrix utilities ────────────────────────────────────────────────

def invert_matrix(mat):
    n = len(mat)
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)]
           for i, row in enumerate(mat)]
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


# ── OLS (simple, no clustering) ─────────────────────────────────────

def ols_simple(y, X_mat):
    """OLS returning betas, residuals, R2. No SEs (computed externally)."""
    n = len(y)
    k = len(X_mat[0])
    if n <= k:
        return None
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    if ss_tot < 1e-15:
        return None
    XtX = [[sum(X_mat[i][a] * X_mat[i][b] for i in range(n))
            for b in range(k)] for a in range(k)]
    Xty = [sum(X_mat[i][a] * y[i] for i in range(n)) for a in range(k)]
    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None
    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]
    y_hat = [sum(X_mat[i][a] * beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]
    ss_res = sum(r ** 2 for r in resid)
    r2 = 1 - ss_res / ss_tot
    return {'beta': beta, 'resid': resid, 'r2': r2, 'n': n}


def _normal_cdf(x):
    if x < -8:
        return 0.0
    if x > 8:
        return 1.0
    ax = abs(x)
    b0 = 0.2316419
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    t_val = 1.0 / (1.0 + b0 * ax)
    phi = (1.0 / math.sqrt(2.0 * math.pi)) * math.exp(-0.5 * ax * ax)
    cdf = 1.0 - phi * (b1 * t_val + b2 * t_val**2 + b3 * t_val**3
                        + b4 * t_val**4 + b5 * t_val**5)
    if x < 0:
        return 1.0 - cdf
    return cdf


def p_from_t(t_stat):
    return 2.0 * (1.0 - _normal_cdf(abs(t_stat)))


# ── Newey-West HAC standard errors ──────────────────────────────────

def newey_west_se(series, max_lag=None):
    """Newey-West (1987) HAC standard error for a time series of scalars."""
    T = len(series)
    if T < 3:
        return float('inf')
    mean = sum(series) / T
    demean = [x - mean for x in series]

    if max_lag is None:
        max_lag = max(1, int(4 * (T / 100) ** (2 / 9)))
    max_lag = min(max_lag, T - 1)

    gamma_0 = sum(d * d for d in demean) / T

    nw_var = gamma_0
    for lag in range(1, max_lag + 1):
        weight = 1.0 - lag / (max_lag + 1.0)
        gamma_lag = sum(demean[t] * demean[t - lag] for t in range(lag, T)) / T
        nw_var += 2.0 * weight * gamma_lag

    var_mean = nw_var / T
    if var_mean < 0:
        var_mean = gamma_0 / T
    return math.sqrt(var_mean)


# ── Haversine ───────────────────────────────────────────────────────

def haversine(lat1, lon1, lat2, lon2):
    """Distance in km between two GPS points."""
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(dlon / 2) ** 2)
    return R * 2 * math.asin(math.sqrt(a))


# ── Load centroids ──────────────────────────────────────────────────

_print('Loading firm centroids...')
centroids = {}
with open(derived_path('networks', 'firm_centroids.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        try:
            lat = float(row['centroid_lat'])
            lon = float(row['centroid_lon'])
        except (ValueError, TypeError):
            continue
        centroids[gk] = (lat, lon)
_print(f'  Firms with centroids: {len(centroids)}')


# ── Build W_geo for a given half-life ───────────────────────────────

def build_W_geo(half_life_km):
    """Build row-normalized geographic weight matrix with given half-life.

    w_ij = exp(-d_ij / DECAY_KM) / d_ij, then row-normalize.
    DECAY_KM = half_life_km / ln(2)
    """
    decay_km = half_life_km / math.log(2)
    gvkeys = sorted(centroids.keys())
    W = {}
    for gi in gvkeys:
        lat_i, lon_i = centroids[gi]
        neighbors = {}
        row_sum = 0.0
        for gj in gvkeys:
            if gi == gj:
                continue
            lat_j, lon_j = centroids[gj]
            d = haversine(lat_i, lon_i, lat_j, lon_j)
            if d > 0:
                w = math.exp(-d / decay_km) / d
                neighbors[gj] = w
                row_sum += w
        if row_sum > 0:
            W[gi] = {gj: w / row_sum for gj, w in neighbors.items()}
        else:
            W[gi] = {}
    return W


# ── Load shared data (same as strategy2_robust_inference.py) ────────

_print('Loading monthly returns...')
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except ValueError:
            pass
_print(f'  Monthly: {len(monthly_ret)} firms')


_print('Loading Fama-French factors...')


def load_ff_factors_monthly(path):
    if not os.path.exists(path):
        return None
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
            vw = (mktrf_val + rf_val) / 100.0
            vwretd[f'{date[:4]}-{date[4:6]}'] = vw
    return vwretd


market_ret_monthly = load_ff_factors_monthly(
    raw_path('factors', 'F-F_Research_Data_Factors.csv')
)
if not market_ret_monthly:
    raise RuntimeError('Missing F-F monthly factors.')
_print(f'  Market months: {len(market_ret_monthly)}')


# Load fuel weight matrix (fixed across bandwidths)
_print('Loading W_fuel...')
W_fuel = defaultdict(dict)
fuel_path = derived_path('networks', 'weight_matrix_W_fuel.csv')
if os.path.exists(fuel_path):
    with open(fuel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            W_fuel[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
_print(f'  W_fuel firms: {len(W_fuel)}')


# Load fundamentals (for SIC codes)
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row


def get_sic4(gvkey):
    f = fundamentals.get(gvkey)
    if f and f.get('sic'):
        return f['sic'][:4]
    return None


# Load events
_print('Loading events...')
all_events = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        if row.get('is_first_mover') != 'True':
            continue
        ann_date = row.get('announcement_date', '').strip()
        ret_date = row.get('event_date', '').strip()
        effective_date = ann_date if ann_date else ret_date
        event_year = None
        if effective_date and len(effective_date) >= 4 and effective_date[:4].isdigit():
            event_year = int(effective_date[:4])
        else:
            event_year = int(row['ret_year']) if row.get('ret_year') else None
        all_events.append({
            'plant': row['plant_name'],
            'year': event_year,
            'event_date': effective_date,
            'gvkeys': row['matched_gvkeys'].split(';'),
        })
_print(f'  First-mover events: {len(all_events)}')


POST_MONTHS = 3
PRE_MONTHS = 24


def compute_monthly_car(gvkey, event_month, post=3):
    if gvkey not in monthly_ret:
        return None
    months = sorted(monthly_ret[gvkey].keys())
    event_idx = None
    for i, m in enumerate(months):
        if m >= event_month:
            event_idx = i
            break
    if event_idx is None:
        return None
    pre_rets = [monthly_ret[gvkey][months[i]]
                for i in range(max(0, event_idx - PRE_MONTHS), event_idx)
                if months[i] in monthly_ret[gvkey]]
    if len(pre_rets) < 12:
        return None
    ar_list = []
    for i in range(max(0, event_idx - PRE_MONTHS), event_idx):
        m = months[i]
        if m in monthly_ret[gvkey] and m in market_ret_monthly:
            ar_list.append(monthly_ret[gvkey][m] - market_ret_monthly[m])
    pre_mean_ar = (sum(ar_list) / len(ar_list)) if ar_list else 0.0
    car = 0.0
    for offset in range(-1, post + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months) and months[idx] in monthly_ret[gvkey]:
            m = months[idx]
            r_it = monthly_ret[gvkey][m]
            if m in market_ret_monthly:
                ar = r_it - market_ret_monthly[m]
                car += ar - pre_mean_ar
    return car


# ══════════════════════════════════════════════════════════════════════
# BANDWIDTH SENSITIVITY ANALYSIS
# ══════════════════════════════════════════════════════════════════════

BANDWIDTHS = [250, 500, 750, 1000, 1500]
MIN_OBS_PER_EVENT = 20

bandwidth_results = []  # list of dicts, one per bandwidth

for bw in BANDWIDTHS:
    _print(f'\n{"=" * 70}')
    _print(f'BANDWIDTH: {bw} km half-life')
    _print(f'{"=" * 70}')

    # 1. Build W_geo for this bandwidth
    _print(f'  Building W_geo (half-life = {bw} km)...')
    W_geo_bw = build_W_geo(bw)
    n_connected = sum(1 for gk in W_geo_bw if len(W_geo_bw[gk]) > 0)
    _print(f'  Connected firms: {n_connected}')

    # 2. Build per-event datasets using this W_geo
    _print(f'  Building per-event datasets...')
    event_datasets = {}

    for event_id, event in enumerate(all_events):
        event_gvkeys = set(event['gvkeys'])
        year = event['year']
        event_date = event.get('event_date', '')
        if event_date and len(event_date) >= 7:
            event_month = event_date[:7]
        else:
            event_month = f'{year}-07' if year else None
        if not event_month:
            continue

        fm_sic4 = None
        for gk in event_gvkeys:
            fm_sic4 = get_sic4(gk)
            if fm_sic4:
                break

        obs = []
        for fm_gk in event_gvkeys:
            if fm_gk not in W_geo_bw:
                continue
            neighbors = W_geo_bw[fm_gk]
            neighbor_gks = set(neighbors.keys()) - event_gvkeys
            non_connected = [gk for gk in fundamentals
                             if gk not in event_gvkeys and gk not in neighbors]
            stable_seed = int(hashlib.md5(
                str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
            random.seed(stable_seed)
            n_ctrl = min(len(non_connected),
                         max(5 * len(neighbor_gks), 20))
            ctrl_sample = (random.sample(non_connected, n_ctrl)
                           if len(non_connected) > n_ctrl
                           else non_connected)
            candidate_firms = list(neighbor_gks) + ctrl_sample

            for gk in candidate_firms:
                w_geo = neighbors.get(gk, 0.0)
                w_fuel = W_fuel.get(fm_gk, {}).get(gk, 0.0)
                car = compute_monthly_car(gk, event_month, post=POST_MONTHS)
                if car is None:
                    continue
                obs.append({
                    'car': car,
                    'w_geo': w_geo,
                    'w_fuel': w_fuel,
                    'gvkey': gk,
                })

        if len(obs) >= MIN_OBS_PER_EVENT:
            event_datasets[event_id] = obs

    n_valid = len(event_datasets)
    total_obs = sum(len(v) for v in event_datasets.values())
    _print(f'  Valid events (>= {MIN_OBS_PER_EVENT} obs): {n_valid}')
    _print(f'  Total obs: {total_obs}')

    # 3. Fama-MacBeth: cross-sectional OLS per event
    event_betas_geo = []
    event_betas_fuel = []
    event_r2s = []
    event_ns = []
    event_ids_used = []

    SPEC_VARS = ['w_geo', 'w_fuel']

    for eid in sorted(event_datasets.keys()):
        obs = event_datasets[eid]
        y = [o['car'] for o in obs]
        X = [[1.0, o['w_geo'], o['w_fuel']] for o in obs]

        result = ols_simple(y, X)
        if result is None:
            continue

        event_betas_geo.append(result['beta'][1])
        event_betas_fuel.append(result['beta'][2])
        event_r2s.append(result['r2'])
        event_ns.append(result['n'])
        event_ids_used.append(eid)

    T_fm = len(event_ids_used)
    if T_fm < 3:
        _print(f'  WARNING: Only {T_fm} valid FM events, skipping bandwidth {bw}')
        continue

    avg_n = sum(event_ns) / T_fm
    avg_r2 = sum(event_r2s) / T_fm
    _print(f'  FM events: {T_fm}, avg N: {avg_n:.1f}, avg R2: {avg_r2:.4f}')

    # 4. Newey-West SEs on beta time series
    mean_geo = sum(event_betas_geo) / T_fm
    nw_se_geo = newey_west_se(event_betas_geo)
    t_geo = mean_geo / nw_se_geo if nw_se_geo > 1e-15 else 0.0
    p_geo = p_from_t(t_geo)

    mean_fuel = sum(event_betas_fuel) / T_fm
    nw_se_fuel = newey_west_se(event_betas_fuel)
    t_fuel = mean_fuel / nw_se_fuel if nw_se_fuel > 1e-15 else 0.0
    p_fuel = p_from_t(t_fuel)

    # 5. Difference test (geo - fuel)
    diff_series = [event_betas_geo[t] - event_betas_fuel[t] for t in range(T_fm)]
    mean_diff = sum(diff_series) / T_fm
    nw_se_diff = newey_west_se(diff_series)
    t_diff = mean_diff / nw_se_diff if nw_se_diff > 1e-15 else 0.0
    p_diff = p_from_t(t_diff)

    stars_geo = '***' if p_geo < 0.01 else '**' if p_geo < 0.05 else '*' if p_geo < 0.10 else ''
    stars_fuel = '***' if p_fuel < 0.01 else '**' if p_fuel < 0.05 else '*' if p_fuel < 0.10 else ''
    stars_diff = '***' if p_diff < 0.01 else '**' if p_diff < 0.05 else '*' if p_diff < 0.10 else ''

    _print(f'\n  beta_geo  = {mean_geo:+.6f}  SE = {nw_se_geo:.6f}  t = {t_geo:.3f}  p = {p_geo:.4f}{stars_geo}')
    _print(f'  beta_fuel = {mean_fuel:+.6f}  SE = {nw_se_fuel:.6f}  t = {t_fuel:.3f}  p = {p_fuel:.4f}{stars_fuel}')
    _print(f'  diff(g-f) = {mean_diff:+.6f}  SE = {nw_se_diff:.6f}  t = {t_diff:.3f}  p = {p_diff:.4f}{stars_diff}')

    bandwidth_results.append({
        'bw': bw,
        'T_fm': T_fm,
        'avg_n': avg_n,
        'avg_r2': avg_r2,
        'beta_geo': mean_geo,
        'se_geo': nw_se_geo,
        't_geo': t_geo,
        'p_geo': p_geo,
        'stars_geo': stars_geo,
        'beta_fuel': mean_fuel,
        'se_fuel': nw_se_fuel,
        't_fuel': t_fuel,
        'p_fuel': p_fuel,
        'stars_fuel': stars_fuel,
        'diff': mean_diff,
        'se_diff': nw_se_diff,
        't_diff': t_diff,
        'p_diff': p_diff,
        'stars_diff': stars_diff,
    })


# ══════════════════════════════════════════════════════════════════════
# WRITE OUTPUT
# ══════════════════════════════════════════════════════════════════════

out_path = results_path('metrics', 'strategy2_bandwidth_fmb.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

lines = [
    '# Bandwidth Sensitivity: W_geo under Fama-MacBeth + Newey-West Inference',
    '',
    'Tests whether the geographic proximity channel (w_geo) strengthens at',
    'wider bandwidths. Economic rationale: competitive benefits from coal',
    'retirement propagate through interconnected transmission grids',
    '(ENTSO-E ~1500km, US ISOs ~1000km), not at plant-level proximity.',
    '',
    f'Events: {len(all_events)} first-mover coal retirements',
    f'Window: [-1, +{POST_MONTHS}] months, vwretd market-adjusted returns',
    f'Minimum obs per event: {MIN_OBS_PER_EVENT}',
    '',
    'Weight formula: w_ij = exp(-d_ij / DECAY_KM) / d_ij, row-normalized',
    'DECAY_KM = half_life / ln(2)',
    '',
    '## Main Results: Bandwidth x Channel',
    '',
    '| Bandwidth (km) | Events | Avg N | beta_geo | SE(NW) | t | p | beta_fuel | SE(NW) | t | p | diff(g-f) | SE(NW) | t | p |',
    '|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
]

for r in bandwidth_results:
    lines.append(
        f'| {r["bw"]} | {r["T_fm"]} | {r["avg_n"]:.0f} '
        f'| {r["beta_geo"]:+.6f} | {r["se_geo"]:.6f} | {r["t_geo"]:.3f} | {r["p_geo"]:.4f}{r["stars_geo"]} '
        f'| {r["beta_fuel"]:+.6f} | {r["se_fuel"]:.6f} | {r["t_fuel"]:.3f} | {r["p_fuel"]:.4f}{r["stars_fuel"]} '
        f'| {r["diff"]:+.6f} | {r["se_diff"]:.6f} | {r["t_diff"]:.3f} | {r["p_diff"]:.4f}{r["stars_diff"]} |'
    )

# Add compact summary
lines += [
    '',
    '## Summary: How geo significance varies with bandwidth',
    '',
    '| Bandwidth | t(geo) | p(geo) | t(fuel) | p(fuel) | Interpretation |',
    '|---:|---:|---:|---:|---:|---|',
]

for r in bandwidth_results:
    if r['p_geo'] < 0.05:
        interp = 'geo significant at 5%'
    elif r['p_geo'] < 0.10:
        interp = 'geo significant at 10%'
    else:
        interp = 'geo not significant'
    if r['p_fuel'] < 0.05:
        interp += '; fuel significant'
    lines.append(
        f'| {r["bw"]} km | {r["t_geo"]:.3f} | {r["p_geo"]:.4f}{r["stars_geo"]} '
        f'| {r["t_fuel"]:.3f} | {r["p_fuel"]:.4f}{r["stars_fuel"]} | {interp} |'
    )

# Identify optimal bandwidth
if bandwidth_results:
    best = max(bandwidth_results, key=lambda r: abs(r['t_geo']))
    lines += [
        '',
        f'**Strongest geo channel: {best["bw"]} km half-life** '
        f'(t = {best["t_geo"]:.3f}, p = {best["p_geo"]:.4f})',
    ]
    # Check if 1000km is better than 500km
    bw500 = [r for r in bandwidth_results if r['bw'] == 500]
    bw1000 = [r for r in bandwidth_results if r['bw'] == 1000]
    if bw500 and bw1000:
        lines.append('')
        lines.append(f'Comparison: 500km t(geo) = {bw500[0]["t_geo"]:.3f} vs '
                     f'1000km t(geo) = {bw1000[0]["t_geo"]:.3f}')

lines.append('')

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

_print(f'\nWrote: {out_path}')
_print('Done.')
