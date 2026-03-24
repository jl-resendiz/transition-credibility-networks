"""Robust inference alternatives for the channel decomposition.

The two-way clustered SEs (event + firm) show that the original event-only
clustered t-stats were inflated. This script implements three alternatives
that properly account for within-firm correlation:

1. Fama-MacBeth (1973): Run cross-sectional regression event-by-event,
   average coefficients, and compute Newey-West (1987) SEs on the time
   series of betas. This is the gold standard for repeated cross-sections.

2. Event-level collapse: For each event, compute mean CAR by exposure
   quintile, then test the quintile spread over the series of events.

3. Portfolio sorts with proper inference: Shanken (1992) / Newey-West
   corrected t-stats on the time series of event-level portfolio spreads.

Output: results/metrics/strategy2_robust_inference.md
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


def mat_mul(a, b):
    rows_a = len(a)
    cols_b = len(b[0])
    mid = len(b)
    out = [[0.0 for _ in range(cols_b)] for _ in range(rows_a)]
    for i in range(rows_a):
        for k in range(mid):
            aik = a[i][k]
            if aik == 0:
                continue
            for j in range(cols_b):
                out[i][j] += aik * b[k][j]
    return out


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
    """Newey-West (1987) HAC standard error for a time series of scalars.

    Computes the variance of the sample mean accounting for autocorrelation.
    Uses Bartlett kernel with automatic lag selection (floor(4*(T/100)^{2/9})).
    """
    T = len(series)
    if T < 3:
        return float('inf')
    mean = sum(series) / T
    demean = [x - mean for x in series]

    if max_lag is None:
        max_lag = max(1, int(4 * (T / 100) ** (2 / 9)))
    max_lag = min(max_lag, T - 1)

    # Gamma_0
    gamma_0 = sum(d * d for d in demean) / T

    # Sum weighted autocovariances
    nw_var = gamma_0
    for lag in range(1, max_lag + 1):
        weight = 1.0 - lag / (max_lag + 1.0)  # Bartlett kernel
        gamma_lag = sum(demean[t] * demean[t - lag] for t in range(lag, T)) / T
        nw_var += 2.0 * weight * gamma_lag

    # Variance of the mean = nw_var / T
    var_mean = nw_var / T
    if var_mean < 0:
        # Can happen with aggressive lag; fall back to simple
        var_mean = gamma_0 / T
    return math.sqrt(var_mean)


def newey_west_cov(series_list, max_lag=None):
    """Newey-West covariance matrix for a list of K time series.

    Returns KxK matrix where [i][j] is the HAC covariance of means.
    """
    K = len(series_list)
    T = len(series_list[0])
    if T < 3:
        return None

    means = [sum(s) / T for s in series_list]
    demean = [[series_list[k][t] - means[k] for t in range(T)] for k in range(K)]

    if max_lag is None:
        max_lag = max(1, int(4 * (T / 100) ** (2 / 9)))
    max_lag = min(max_lag, T - 1)

    # Gamma_0 matrix
    S = [[0.0] * K for _ in range(K)]
    for i in range(K):
        for j in range(i, K):
            g0 = sum(demean[i][t] * demean[j][t] for t in range(T)) / T
            S[i][j] = g0
            S[j][i] = g0

    # Add weighted autocovariance lags
    for lag in range(1, max_lag + 1):
        weight = 1.0 - lag / (max_lag + 1.0)
        for i in range(K):
            for j in range(K):
                g_lag = sum(demean[i][t] * demean[j][t - lag]
                            for t in range(lag, T)) / T
                g_lag_rev = sum(demean[j][t] * demean[i][t - lag]
                                for t in range(lag, T)) / T
                S[i][j] += weight * (g_lag + g_lag_rev)

    # Variance of means = S / T
    V = [[S[i][j] / T for j in range(K)] for i in range(K)]
    return V


# ── Load data ────────────────────────────────────────────────────────

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

_print('Loading weight matrices...')
W_geo = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_geo[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

W_fuel = defaultdict(dict)
fuel_path = derived_path('networks', 'weight_matrix_W_fuel.csv')
if os.path.exists(fuel_path):
    with open(fuel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            W_fuel[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

W_reg = defaultdict(dict)
reg_path = derived_path('networks', 'weight_matrix_W_regulatory.csv')
if os.path.exists(reg_path):
    with open(reg_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            wval = row.get('w_ij')
            if wval in (None, ''):
                wval = row.get('w_reg')
            try:
                W_reg[row['gvkey_i']][row['gvkey_j']] = float(wval)
            except (ValueError, TypeError):
                continue

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


# ── Build per-event datasets ─────────────────────────────────────────

_print('\nBuilding per-event datasets...')

SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']
MIN_OBS_PER_EVENT = 20  # minimum firms per event for Fama-MacBeth

event_datasets = {}  # event_id -> list of obs dicts

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
        if fm_gk not in W_geo:
            continue
        neighbors = W_geo[fm_gk]
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
            w_reg = W_reg.get(fm_gk, {}).get(gk, 0.0)
            j_sic4 = get_sic4(gk)
            same_sector = 1.0 if (fm_sic4 and j_sic4
                                  and fm_sic4 == j_sic4) else 0.0
            car = compute_monthly_car(gk, event_month, post=POST_MONTHS)
            if car is None:
                continue
            obs.append({
                'car': car,
                'w_geo': w_geo,
                'w_fuel': w_fuel,
                'w_reg': w_reg,
                'same_sector': same_sector,
                'gvkey': gk,
            })

    if len(obs) >= MIN_OBS_PER_EVENT:
        event_datasets[event_id] = obs

n_valid_events = len(event_datasets)
_print(f'  Valid events (>= {MIN_OBS_PER_EVENT} obs): {n_valid_events}')
_print(f'  Total obs: {sum(len(v) for v in event_datasets.values())}')


# ══════════════════════════════════════════════════════════════════════
# APPROACH 1: FAMA-MACBETH
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('APPROACH 1: FAMA-MACBETH (1973)')
_print('Cross-sectional regression per event, then average betas')
_print('Newey-West (1987) HAC standard errors on the beta time series')
_print('=' * 70)

# Run cross-sectional OLS for each event
event_betas = defaultdict(list)  # var_name -> list of betas (one per event)
event_r2s = []
event_ns = []
event_ids_used = []

for event_id in sorted(event_datasets.keys()):
    obs = event_datasets[event_id]
    n_obs = len(obs)

    # Check same_sector variation
    ss_vals = set(o['same_sector'] for o in obs)
    use_vars = SPEC_VARS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']

    y = [o['car'] for o in obs]
    X = [[1.0] + [o[v] for v in use_vars] for o in obs]

    result = ols_simple(y, X)
    if result is None:
        continue

    names = ['intercept'] + use_vars
    for i, name in enumerate(names):
        event_betas[name].append(result['beta'][i])

    # Pad missing variables with NaN
    for v in SPEC_VARS:
        if v not in use_vars:
            event_betas[v].append(float('nan'))

    event_r2s.append(result['r2'])
    event_ns.append(result['n'])
    event_ids_used.append(event_id)

T_fm = len(event_ids_used)
_print(f'\n  Events with valid regressions: {T_fm}')
_print(f'  Avg N per event: {sum(event_ns) / T_fm:.1f}')
_print(f'  Avg R2 per event: {sum(event_r2s) / T_fm:.4f}')

# Export event-level betas to CSV for figure generation
_event_beta_path = results_path('summaries', 'event_level_betas.csv')
os.makedirs(os.path.dirname(_event_beta_path), exist_ok=True)
with open(_event_beta_path, 'w', newline='', encoding='utf-8') as _f:
    _w = csv.writer(_f)
    _w.writerow(['event_id', 'event_n', 'event_r2', 'beta_fuel', 'beta_geo',
                 'beta_reg', 'beta_same_sector'])
    for _idx in range(T_fm):
        _w.writerow([
            event_ids_used[_idx],
            event_ns[_idx],
            f'{event_r2s[_idx]:.6f}',
            f'{event_betas["w_fuel"][_idx]:.6f}' if not math.isnan(event_betas['w_fuel'][_idx]) else '',
            f'{event_betas["w_geo"][_idx]:.6f}' if not math.isnan(event_betas['w_geo'][_idx]) else '',
            f'{event_betas["w_reg"][_idx]:.6f}' if not math.isnan(event_betas['w_reg'][_idx]) else '',
            f'{event_betas["same_sector"][_idx]:.6f}' if not math.isnan(event_betas['same_sector'][_idx]) else '',
        ])
_print(f'  Exported event-level betas to: {_event_beta_path}')

# Compute Fama-MacBeth averages and Newey-West SEs
_print('\n  Fama-MacBeth coefficients with Newey-West SEs:')
_print(f'  {"Variable":<15} {"Mean beta":>12} {"NW SE":>10} {"t":>8} {"p":>8}')
_print('  ' + '-' * 55)

fm_results = {}
for v in ['intercept'] + SPEC_VARS:
    betas = event_betas[v]
    # Remove NaN
    clean = [b for b in betas if not math.isnan(b)]
    if len(clean) < 3:
        continue
    mean_b = sum(clean) / len(clean)
    nw_se = newey_west_se(clean)
    t_stat = mean_b / nw_se if nw_se > 1e-15 else 0.0
    p_val = p_from_t(t_stat)
    stars = '***' if p_val < 0.01 else '**' if p_val < 0.05 else '*' if p_val < 0.10 else ''
    _print(f'  {v:<15} {mean_b:+12.6f} {nw_se:10.6f} {t_stat:8.3f} {p_val:8.4f}{stars}')
    fm_results[v] = {
        'mean': mean_b, 'se': nw_se, 't': t_stat, 'p': p_val,
        'n_events': len(clean),
    }

# FM difference test (geo - fuel)
if 'w_geo' in fm_results and 'w_fuel' in fm_results:
    geo_betas = [b for b in event_betas['w_geo'] if not math.isnan(b)]
    fuel_betas = [b for b in event_betas['w_fuel'] if not math.isnan(b)]
    # Align series (same events)
    T_diff = min(len(geo_betas), len(fuel_betas))
    diff_series = [geo_betas[t] - fuel_betas[t] for t in range(T_diff)]
    mean_diff = sum(diff_series) / T_diff
    nw_se_diff = newey_west_se(diff_series)
    t_diff = mean_diff / nw_se_diff if nw_se_diff > 1e-15 else 0.0
    p_diff = p_from_t(t_diff)
    stars_d = '***' if p_diff < 0.01 else '**' if p_diff < 0.05 else '*' if p_diff < 0.10 else ''
    _print(f'\n  Difference test (beta_geo - beta_fuel):')
    _print(f'    Mean diff = {mean_diff:+.6f}, NW SE = {nw_se_diff:.6f}, '
           f't = {t_diff:.3f}, p = {p_diff:.4f}{stars_d}')

    # Also: joint Wald test on [beta_geo, beta_fuel, beta_reg]
    # using Newey-West covariance matrix of the beta series
    channels_present = [v for v in ['w_geo', 'w_fuel', 'w_reg']
                        if v in fm_results and fm_results[v]['n_events'] == T_fm]
    if len(channels_present) >= 2:
        series_list = [[b for b in event_betas[v] if not math.isnan(b)]
                       for v in channels_present]
        # Ensure all same length
        min_len = min(len(s) for s in series_list)
        series_list = [s[:min_len] for s in series_list]
        means = [sum(s) / len(s) for s in series_list]
        V_nw = newey_west_cov(series_list)
        if V_nw:
            inv_V = invert_matrix(V_nw)
            if inv_V:
                q = len(channels_present)
                quad = sum(means[a] * sum(inv_V[a][b] * means[b]
                                          for b in range(q))
                           for a in range(q))
                f_fm = quad / q
                _print(f'\n  Joint Wald F-test (FM + NW): F = {f_fm:.4f}, q = {q}')


# ══════════════════════════════════════════════════════════════════════
# APPROACH 2: EVENT-LEVEL PORTFOLIO SPREADS
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('APPROACH 2: EVENT-LEVEL PORTFOLIO SORTS')
_print('For each event, form quintiles on exposure, compute Q5-Q1 spread.')
_print('Then test the time series of spreads with Newey-West SEs.')
_print('=' * 70)

fuel_spreads = []
geo_spreads = []
channel_split_spreads = []

for event_id in sorted(event_datasets.keys()):
    obs = event_datasets[event_id]
    n_obs = len(obs)
    if n_obs < 25:  # need at least 5 per quintile
        continue

    # Fuel quintiles
    sorted_fuel = sorted(obs, key=lambda o: o['w_fuel'])
    q_size = n_obs // 5
    if q_size < 1:
        continue
    q1_fuel = sorted_fuel[:q_size]
    q5_fuel = sorted_fuel[-q_size:]
    mean_q1_fuel = sum(o['car'] for o in q1_fuel) / len(q1_fuel)
    mean_q5_fuel = sum(o['car'] for o in q5_fuel) / len(q5_fuel)
    fuel_spread = mean_q5_fuel - mean_q1_fuel
    fuel_spreads.append(fuel_spread)

    # Geo quintiles
    sorted_geo = sorted(obs, key=lambda o: o['w_geo'])
    q1_geo = sorted_geo[:q_size]
    q5_geo = sorted_geo[-q_size:]
    mean_q1_geo = sum(o['car'] for o in q1_geo) / len(q1_geo)
    mean_q5_geo = sum(o['car'] for o in q5_geo) / len(q5_geo)
    geo_spread = mean_q5_geo - mean_q1_geo
    geo_spreads.append(geo_spread)

    channel_split_spreads.append(geo_spread - fuel_spread)

T_ps = len(fuel_spreads)
_print(f'\n  Events with valid quintile sorts: {T_ps}')

# Fuel spread
mean_fuel = sum(fuel_spreads) / T_ps
nw_fuel = newey_west_se(fuel_spreads)
t_fuel = mean_fuel / nw_fuel if nw_fuel > 1e-15 else 0.0
p_fuel = p_from_t(t_fuel)

# Geo spread
mean_geo = sum(geo_spreads) / T_ps
nw_geo = newey_west_se(geo_spreads)
t_geo = mean_geo / nw_geo if nw_geo > 1e-15 else 0.0
p_geo = p_from_t(t_geo)

# Channel split
mean_split = sum(channel_split_spreads) / T_ps
nw_split = newey_west_se(channel_split_spreads)
t_split = mean_split / nw_split if nw_split > 1e-15 else 0.0
p_split = p_from_t(t_split)

_print(f'\n  {"Spread":<25} {"Mean":>10} {"NW SE":>10} {"t":>8} {"p":>8}')
_print('  ' + '-' * 60)
for label, mn, se, t, p in [
    ('Fuel Q5-Q1', mean_fuel, nw_fuel, t_fuel, p_fuel),
    ('Geo Q5-Q1', mean_geo, nw_geo, t_geo, p_geo),
    ('Channel split (G-F)', mean_split, nw_split, t_split, p_split),
]:
    stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
    _print(f'  {label:<25} {mn:+10.4f} {se:10.4f} {t:8.3f} {p:8.4f}{stars}')

# Also: simple t-test (no NW correction) for comparison
simple_se_fuel = math.sqrt(sum((x - mean_fuel)**2 for x in fuel_spreads)
                           / (T_ps * (T_ps - 1)))
simple_t_fuel = mean_fuel / simple_se_fuel if simple_se_fuel > 0 else 0
simple_se_geo = math.sqrt(sum((x - mean_geo)**2 for x in geo_spreads)
                          / (T_ps * (T_ps - 1)))
simple_t_geo = mean_geo / simple_se_geo if simple_se_geo > 0 else 0
simple_se_split = math.sqrt(sum((x - mean_split)**2 for x in channel_split_spreads)
                            / (T_ps * (T_ps - 1)))
simple_t_split = mean_split / simple_se_split if simple_se_split > 0 else 0

_print(f'\n  For comparison (simple SEs, no HAC correction):')
_print(f'  {"Spread":<25} {"t(simple)":>10} {"t(NW)":>10}')
_print('  ' + '-' * 45)
_print(f'  {"Fuel Q5-Q1":<25} {simple_t_fuel:10.3f} {t_fuel:10.3f}')
_print(f'  {"Geo Q5-Q1":<25} {simple_t_geo:10.3f} {t_geo:10.3f}')
_print(f'  {"Channel split":<25} {simple_t_split:10.3f} {t_split:10.3f}')


# ══════════════════════════════════════════════════════════════════════
# APPROACH 3: LONG-SHORT PORTFOLIO WITH PROPER INFERENCE
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('APPROACH 3: LONG-SHORT PORTFOLIO (NW-corrected)')
_print('Long high-geo/low-fuel, short low-geo/high-fuel')
_print('=' * 70)

ls_returns = []

for event_id in sorted(event_datasets.keys()):
    obs = event_datasets[event_id]
    n_obs = len(obs)
    if n_obs < 25:
        continue
    q_size = n_obs // 5

    # Long: high geo + low fuel
    sorted_geo = sorted(obs, key=lambda o: o['w_geo'])
    sorted_fuel = sorted(obs, key=lambda o: o['w_fuel'])
    high_geo = set(o['gvkey'] for o in sorted_geo[-q_size:])
    low_fuel = set(o['gvkey'] for o in sorted_fuel[:q_size])
    long_firms = high_geo | low_fuel
    long_cars = [o['car'] for o in obs if o['gvkey'] in long_firms]

    # Short: low geo + high fuel
    low_geo = set(o['gvkey'] for o in sorted_geo[:q_size])
    high_fuel = set(o['gvkey'] for o in sorted_fuel[-q_size:])
    short_firms = low_geo | high_fuel
    short_cars = [o['car'] for o in obs if o['gvkey'] in short_firms]

    if long_cars and short_cars:
        ls_ret = sum(long_cars) / len(long_cars) - sum(short_cars) / len(short_cars)
        ls_returns.append(ls_ret)

T_ls = len(ls_returns)
mean_ls = sum(ls_returns) / T_ls
nw_ls = newey_west_se(ls_returns)
t_ls = mean_ls / nw_ls if nw_ls > 1e-15 else 0.0
p_ls = p_from_t(t_ls)
simple_se_ls = math.sqrt(sum((x - mean_ls)**2 for x in ls_returns)
                         / (T_ls * (T_ls - 1)))
simple_t_ls = mean_ls / simple_se_ls if simple_se_ls > 0 else 0

stars_ls = '***' if p_ls < 0.01 else '**' if p_ls < 0.05 else '*' if p_ls < 0.10 else ''
_print(f'\n  N events: {T_ls}')
_print(f'  Mean L/S return: {mean_ls:+.4f} ({mean_ls * 100:+.2f}%)')
_print(f'  NW SE: {nw_ls:.4f}')
_print(f'  t (NW): {t_ls:.3f}, p = {p_ls:.4f}{stars_ls}')
_print(f'  t (simple): {simple_t_ls:.3f} (for comparison)')


# ══════════════════════════════════════════════════════════════════════
# APPROACH 4: EVENT WINDOW SENSITIVITY
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('APPROACH 4: EVENT WINDOW SENSITIVITY')
_print('Pooled OLS with event-clustered SEs at multiple windows')
_print('=' * 70)

window_sensitivity = []  # list of dicts with results per window

for post_m in [1, 2, 3]:
    window_label = f'[-1, +{post_m}]'
    _print(f'\n  Window {window_label} ...')

    # Rebuild CARs for this window
    w_obs_all = []
    w_event_ids = []
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
            if fm_gk not in W_geo:
                continue
            neighbors = W_geo[fm_gk]
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
                w_geo_val = neighbors.get(gk, 0.0)
                w_fuel_val = W_fuel.get(fm_gk, {}).get(gk, 0.0)
                car = compute_monthly_car(gk, event_month, post=post_m)
                if car is None:
                    continue
                obs.append({
                    'car': car,
                    'w_geo': w_geo_val,
                    'w_fuel': w_fuel_val,
                    'event_id': event_id,
                })

        if len(obs) >= MIN_OBS_PER_EVENT:
            w_obs_all.extend(obs)
            w_event_ids.append(event_id)

    if len(w_obs_all) < 10 or len(w_event_ids) < 2:
        _print(f'    Skipping: too few obs ({len(w_obs_all)}) or events ({len(w_event_ids)})')
        continue

    # Pooled OLS: CAR = a + b1*w_fuel + b2*w_geo + e
    y_w = [o['car'] for o in w_obs_all]
    X_w = [[1.0, o['w_fuel'], o['w_geo']] for o in w_obs_all]
    result_w = ols_simple(y_w, X_w)
    if result_w is None:
        _print('    OLS failed')
        continue

    beta_w = result_w['beta']
    resid_w = result_w['resid']
    n_w = result_w['n']
    k_w = 3  # intercept + 2 channels

    # Event-clustered standard errors
    # V = (X'X)^-1 * ( sum_g X_g' e_g e_g' X_g ) * (X'X)^-1
    XtX = [[sum(X_w[i][a] * X_w[i][b] for i in range(n_w))
            for b in range(k_w)] for a in range(k_w)]
    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        _print('    Singular X\'X')
        continue

    # Group residuals by event
    event_groups = defaultdict(list)
    for idx, o in enumerate(w_obs_all):
        event_groups[o['event_id']].append(idx)

    G = len(event_groups)
    meat = [[0.0] * k_w for _ in range(k_w)]
    for eid, indices in event_groups.items():
        # X_g' * e_g  (k x 1 vector)
        score = [sum(X_w[i][a] * resid_w[i] for i in indices) for a in range(k_w)]
        for a in range(k_w):
            for b in range(k_w):
                meat[a][b] += score[a] * score[b]

    # Small-sample correction: G/(G-1) * (n-1)/(n-k)
    correction = (G / (G - 1.0)) * ((n_w - 1.0) / (n_w - k_w))
    for a in range(k_w):
        for b in range(k_w):
            meat[a][b] *= correction

    V_cl = mat_mul(mat_mul(inv_XtX, meat), inv_XtX)

    se_fuel = math.sqrt(max(V_cl[1][1], 0.0))
    se_geo = math.sqrt(max(V_cl[2][2], 0.0))
    t_fuel_w = beta_w[1] / se_fuel if se_fuel > 1e-15 else 0.0
    t_geo_w = beta_w[2] / se_geo if se_geo > 1e-15 else 0.0

    row = {
        'window': window_label,
        'post_months': post_m,
        'n_obs': n_w,
        'n_events': G,
        'fuel_beta': beta_w[1],
        'fuel_se': se_fuel,
        'fuel_t': t_fuel_w,
        'geo_beta': beta_w[2],
        'geo_se': se_geo,
        'geo_t': t_geo_w,
        'r2': result_w['r2'],
    }
    window_sensitivity.append(row)
    _print(f'    N={n_w}, events={G}, fuel beta={beta_w[1]:+.6f} (t={t_fuel_w:.3f}), '
           f'geo beta={beta_w[2]:+.6f} (t={t_geo_w:.3f})')

# Also add [0, +1] window (no pre-event month)
_print(f'\n  Window [0, +1] (custom) ...')
w_obs_01 = []
w_event_ids_01 = []
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

    for fm_gk in event_gvkeys:
        if fm_gk not in W_geo:
            continue
        neighbors = W_geo[fm_gk]
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
            w_geo_val = neighbors.get(gk, 0.0)
            w_fuel_val = W_fuel.get(fm_gk, {}).get(gk, 0.0)
            # Custom CAR for [0, +1]: no pre-event month
            if gk not in monthly_ret:
                continue
            months_gk = sorted(monthly_ret[gk].keys())
            ev_idx = None
            for mi, m in enumerate(months_gk):
                if m >= event_month:
                    ev_idx = mi
                    break
            if ev_idx is None:
                continue
            # Pre-period AR for adjustment
            pre_rets_01 = [monthly_ret[gk][months_gk[pi]]
                           for pi in range(max(0, ev_idx - PRE_MONTHS), ev_idx)
                           if months_gk[pi] in monthly_ret[gk]]
            if len(pre_rets_01) < 12:
                continue
            ar_list_01 = []
            for pi in range(max(0, ev_idx - PRE_MONTHS), ev_idx):
                m = months_gk[pi]
                if m in monthly_ret[gk] and m in market_ret_monthly:
                    ar_list_01.append(monthly_ret[gk][m] - market_ret_monthly[m])
            pre_mean_ar_01 = (sum(ar_list_01) / len(ar_list_01)) if ar_list_01 else 0.0
            car_01 = 0.0
            for offset in range(0, 2):  # month 0 and month +1
                idx = ev_idx + offset
                if 0 <= idx < len(months_gk) and months_gk[idx] in monthly_ret[gk]:
                    m = months_gk[idx]
                    r_it = monthly_ret[gk][m]
                    if m in market_ret_monthly:
                        ar = r_it - market_ret_monthly[m]
                        car_01 += ar - pre_mean_ar_01
            w_obs_01.append({
                'car': car_01,
                'w_geo': w_geo_val,
                'w_fuel': w_fuel_val,
                'event_id': event_id,
            })

    if event_id not in [o['event_id'] for o in w_obs_01]:
        continue
    # Count obs for this event
    n_this = sum(1 for o in w_obs_01 if o['event_id'] == event_id)
    if n_this >= MIN_OBS_PER_EVENT and event_id not in w_event_ids_01:
        w_event_ids_01.append(event_id)

# Filter to events with enough obs
w_obs_01_filtered = [o for o in w_obs_01 if o['event_id'] in w_event_ids_01]

if len(w_obs_01_filtered) >= 10 and len(w_event_ids_01) >= 2:
    y_01 = [o['car'] for o in w_obs_01_filtered]
    X_01 = [[1.0, o['w_fuel'], o['w_geo']] for o in w_obs_01_filtered]
    result_01 = ols_simple(y_01, X_01)
    if result_01 is not None:
        beta_01 = result_01['beta']
        resid_01 = result_01['resid']
        n_01 = result_01['n']
        k_01 = 3

        XtX_01 = [[sum(X_01[i][a] * X_01[i][b] for i in range(n_01))
                    for b in range(k_01)] for a in range(k_01)]
        inv_XtX_01 = invert_matrix(XtX_01)
        if inv_XtX_01 is not None:
            event_groups_01 = defaultdict(list)
            for idx, o in enumerate(w_obs_01_filtered):
                event_groups_01[o['event_id']].append(idx)
            G_01 = len(event_groups_01)
            meat_01 = [[0.0] * k_01 for _ in range(k_01)]
            for eid, indices in event_groups_01.items():
                score = [sum(X_01[i][a] * resid_01[i] for i in indices)
                         for a in range(k_01)]
                for a in range(k_01):
                    for b in range(k_01):
                        meat_01[a][b] += score[a] * score[b]
            correction_01 = (G_01 / (G_01 - 1.0)) * ((n_01 - 1.0) / (n_01 - k_01))
            for a in range(k_01):
                for b in range(k_01):
                    meat_01[a][b] *= correction_01
            V_01 = mat_mul(mat_mul(inv_XtX_01, meat_01), inv_XtX_01)
            se_fuel_01 = math.sqrt(max(V_01[1][1], 0.0))
            se_geo_01 = math.sqrt(max(V_01[2][2], 0.0))
            t_fuel_01 = beta_01[1] / se_fuel_01 if se_fuel_01 > 1e-15 else 0.0
            t_geo_01 = beta_01[2] / se_geo_01 if se_geo_01 > 1e-15 else 0.0

            window_sensitivity.append({
                'window': '[0, +1]',
                'post_months': 1,
                'n_obs': n_01,
                'n_events': G_01,
                'fuel_beta': beta_01[1],
                'fuel_se': se_fuel_01,
                'fuel_t': t_fuel_01,
                'geo_beta': beta_01[2],
                'geo_se': se_geo_01,
                'geo_t': t_geo_01,
                'r2': result_01['r2'],
            })
            _print(f'    N={n_01}, events={G_01}, fuel beta={beta_01[1]:+.6f} (t={t_fuel_01:.3f}), '
                   f'geo beta={beta_01[2]:+.6f} (t={t_geo_01:.3f})')

# Sort sensitivity results by window for display
window_order = {'[-1, +1]': 0, '[-1, +2]': 1, '[-1, +3]': 2, '[0, +1]': 3}
window_sensitivity.sort(key=lambda r: window_order.get(r['window'], 99))


# ══════════════════════════════════════════════════════════════════════
# APPROACH 5: OUTLIER DIAGNOSTICS (Cook's Distance)
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('APPROACH 5: OUTLIER DIAGNOSTICS (Cook\'s Distance)')
_print('Pooled OLS on baseline window, flag obs with Cook\'s D > 4/N')
_print('=' * 70)

# Build pooled dataset from the baseline event_datasets
pooled_obs = []
for event_id in sorted(event_datasets.keys()):
    for o in event_datasets[event_id]:
        pooled_obs.append({
            'car': o['car'],
            'w_fuel': o['w_fuel'],
            'w_geo': o['w_geo'],
            'event_id': event_id,
        })

y_pool = [o['car'] for o in pooled_obs]
X_pool = [[1.0, o['w_fuel'], o['w_geo']] for o in pooled_obs]
n_pool = len(y_pool)
k_pool = 3

result_pool = ols_simple(y_pool, X_pool)
cook_results = None

if result_pool is not None:
    beta_pool = result_pool['beta']
    resid_pool = result_pool['resid']

    # Compute (X'X)^{-1}
    XtX_pool = [[sum(X_pool[i][a] * X_pool[i][b] for i in range(n_pool))
                  for b in range(k_pool)] for a in range(k_pool)]
    inv_XtX_pool = invert_matrix(XtX_pool)

    if inv_XtX_pool is not None:
        # MSE = sum(e^2) / (N - k)
        ss_res_pool = sum(r ** 2 for r in resid_pool)
        mse_pool = ss_res_pool / (n_pool - k_pool)

        # Leverage h_ii = X_i' (X'X)^{-1} X_i
        # Cook's D_i = (e_i^2 / (k * MSE)) * (h_ii / (1 - h_ii)^2)
        cook_d = []
        for i in range(n_pool):
            h_ii = 0.0
            for a in range(k_pool):
                for b in range(k_pool):
                    h_ii += X_pool[i][a] * inv_XtX_pool[a][b] * X_pool[i][b]
            denom = (1.0 - h_ii) ** 2
            if denom < 1e-15:
                cook_d.append(float('inf'))
            else:
                d_i = (resid_pool[i] ** 2 / (k_pool * mse_pool)) * (h_ii / denom)
                cook_d.append(d_i)

        threshold = 4.0 / n_pool
        max_cook = max(cook_d)
        n_high = sum(1 for d in cook_d if d > threshold)

        _print(f'\n  Pooled N: {n_pool}, k: {k_pool}')
        _print(f'  MSE: {mse_pool:.8f}')
        _print(f'  Max Cook\'s D: {max_cook:.6f}')
        _print(f'  Threshold (4/N): {threshold:.6f}')
        _print(f'  Obs with Cook\'s D > 4/N: {n_high} ({100.0 * n_high / n_pool:.1f}%)')
        _print(f'  Fuel beta (full sample): {beta_pool[1]:+.6f}')

        # Re-run dropping high-Cook's-D observations
        y_trim = [y_pool[i] for i in range(n_pool) if cook_d[i] <= threshold]
        X_trim = [X_pool[i] for i in range(n_pool) if cook_d[i] <= threshold]
        n_trim = len(y_trim)

        result_trim = ols_simple(y_trim, X_trim) if n_trim > k_pool else None
        if result_trim is not None:
            _print(f'\n  After dropping {n_high} high-influence obs (N={n_trim}):')
            _print(f'  Fuel beta (trimmed): {result_trim["beta"][1]:+.6f}')
            _print(f'  Geo beta (trimmed):  {result_trim["beta"][2]:+.6f}')
            _print(f'  R2 (trimmed): {result_trim["r2"]:.4f}')

            cook_results = {
                'n_pool': n_pool,
                'max_cook': max_cook,
                'threshold': threshold,
                'n_high': n_high,
                'fuel_beta_full': beta_pool[1],
                'geo_beta_full': beta_pool[2],
                'r2_full': result_pool['r2'],
                'fuel_beta_trim': result_trim['beta'][1],
                'geo_beta_trim': result_trim['beta'][2],
                'r2_trim': result_trim['r2'],
                'n_trim': n_trim,
            }


# ══════════════════════════════════════════════════════════════════════
# WRITE OUTPUT
# ══════════════════════════════════════════════════════════════════════

out_path = results_path('metrics', 'strategy2_robust_inference.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

lines = [
    '# Robust Inference: Three Alternative Approaches',
    '',
    'Two-way clustering (event + firm) showed that event-only clustered',
    't-stats were inflated by within-firm correlation across events.',
    'These three approaches provide valid inference under this structure.',
    '',
    f'Events: {len(all_events)} first-mover coal retirements',
    f'Window: [-1, +{POST_MONTHS}] months, vwretd market-adjusted returns',
    '',
    '## Approach 1: Fama-MacBeth (1973) with Newey-West SEs',
    '',
    'Cross-sectional regression per event, then average betas.',
    f'Valid events: {T_fm} (min {MIN_OBS_PER_EVENT} firms per event)',
    f'Avg firms per event: {sum(event_ns) / T_fm:.1f}',
    f'Avg within-event R2: {sum(event_r2s) / T_fm:.4f}',
    '',
    '| Variable | Mean beta | NW SE | t | p |',
    '|---|---:|---:|---:|---:|',
]

for v in ['intercept'] + SPEC_VARS:
    if v in fm_results:
        r = fm_results[v]
        stars = '***' if r['p'] < 0.01 else '**' if r['p'] < 0.05 else '*' if r['p'] < 0.10 else ''
        lines.append(f'| {v} | {r["mean"]:+.6f} | {r["se"]:.6f} | {r["t"]:.3f} | {r["p"]:.4f}{stars} |')

lines.append('')
if 'w_geo' in fm_results and 'w_fuel' in fm_results:
    lines.append(f'Difference test (FM): beta_geo - beta_fuel = {mean_diff:+.6f} '
                 f'(NW SE = {nw_se_diff:.6f}, t = {t_diff:.3f}, p = {p_diff:.4f})')
    if 'f_fm' in dir():
        lines.append(f'Joint Wald F-test (FM + NW): F = {f_fm:.4f}')

lines += [
    '',
    '## Approach 2: Event-Level Portfolio Sorts (Newey-West SEs)',
    '',
    f'Valid events: {T_ps} (min 25 firms per event for quintile formation)',
    '',
    '| Spread | Mean | NW SE | t(NW) | p | t(simple) |',
    '|---|---:|---:|---:|---:|---:|',
    f'| Fuel Q5-Q1 | {mean_fuel:+.4f} | {nw_fuel:.4f} | {t_fuel:.3f} | {p_fuel:.4f} | {simple_t_fuel:.3f} |',
    f'| Geo Q5-Q1 | {mean_geo:+.4f} | {nw_geo:.4f} | {t_geo:.3f} | {p_geo:.4f} | {simple_t_geo:.3f} |',
    f'| Channel split (G-F) | {mean_split:+.4f} | {nw_split:.4f} | {t_split:.3f} | {p_split:.4f} | {simple_t_split:.3f} |',
    '',
    '## Approach 3: Long-Short Portfolio (Newey-West)',
    '',
    f'Events: {T_ls}',
    f'Mean L/S return: {mean_ls:+.4f} ({mean_ls * 100:+.2f}%)',
    f'NW SE: {nw_ls:.4f}, t(NW) = {t_ls:.3f}, p = {p_ls:.4f}',
    f't(simple) = {simple_t_ls:.3f} (for comparison)',
    '',
    '## Summary: Inference Comparison',
    '',
    '| Method | geo t | fuel t | diff t | Note |',
    '|---|---:|---:|---:|---|',
]

# Pooled event-clustered (from Approach 4, baseline [-1, +3] window)
_baseline_ws = [r for r in window_sensitivity if r['window'] == '[-1, +3]']
geo_t_ec = _baseline_ws[0]['geo_t'] if _baseline_ws else 0.0
fuel_t_ec = _baseline_ws[0]['fuel_t'] if _baseline_ws else 0.0
lines.append(f'| Pooled, event-clustered | {geo_t_ec:.3f} | {fuel_t_ec:.3f} | -- | Primary |')
# Fama-MacBeth
geo_t_fm = fm_results.get('w_geo', {}).get('t', 0)
fuel_t_fm = fm_results.get('w_fuel', {}).get('t', 0)
lines.append(f'| Fama-MacBeth + NW | {geo_t_fm:.3f} | {fuel_t_fm:.3f} | {t_diff:.3f} | Gold standard |')
# Portfolio sorts
lines.append(f'| Portfolio sorts + NW | {t_geo:.3f} | {t_fuel:.3f} | {t_split:.3f} | Non-parametric |')

lines += [
    '',
    '## Interpretation',
    '',
    'The comparison reveals the extent to which the original event-only',
    'clustered inference was inflated by within-firm serial correlation.',
    'The Fama-MacBeth approach is the appropriate gold standard for this',
    'repeated cross-section design (Petersen 2009, Table 5).',
]

# Window sensitivity section
if window_sensitivity:
    lines += [
        '',
        '## Approach 4: Event Window Sensitivity',
        '',
        'Pooled OLS with event-clustered SEs at alternative windows.',
        'Referee-requested robustness check for the baseline [-1, +3] window.',
        '',
        '| Window | N | Events | Fuel beta | Fuel SE | Fuel t | Geo beta | Geo SE | Geo t | R2 |',
        '|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|',
    ]
    for r in window_sensitivity:
        lines.append(
            f'| {r["window"]} | {r["n_obs"]} | {r["n_events"]} '
            f'| {r["fuel_beta"]:+.6f} | {r["fuel_se"]:.6f} | {r["fuel_t"]:.3f} '
            f'| {r["geo_beta"]:+.6f} | {r["geo_se"]:.6f} | {r["geo_t"]:.3f} '
            f'| {r["r2"]:.4f} |'
        )

# Cook's distance section
if cook_results is not None:
    cr = cook_results
    lines += [
        '',
        '## Approach 5: Outlier Diagnostics',
        '',
        'Cook\'s distance on the pooled OLS baseline specification.',
        f'Threshold: 4/N = {cr["threshold"]:.6f}',
        '',
        '| Metric | Value |',
        '|---|---:|',
        f'| N (full sample) | {cr["n_pool"]} |',
        f'| Max Cook\'s D | {cr["max_cook"]:.6f} |',
        f'| Obs with D > 4/N | {cr["n_high"]} ({100.0 * cr["n_high"] / cr["n_pool"]:.1f}%) |',
        f'| Fuel beta (full) | {cr["fuel_beta_full"]:+.6f} |',
        f'| Geo beta (full) | {cr["geo_beta_full"]:+.6f} |',
        f'| R2 (full) | {cr["r2_full"]:.4f} |',
        f'| N (trimmed) | {cr["n_trim"]} |',
        f'| Fuel beta (trimmed) | {cr["fuel_beta_trim"]:+.6f} |',
        f'| Geo beta (trimmed) | {cr["geo_beta_trim"]:+.6f} |',
        f'| R2 (trimmed) | {cr["r2_trim"]:.4f} |',
    ]

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

_print(f'\nWrote: {out_path}')
_print('Done.')
