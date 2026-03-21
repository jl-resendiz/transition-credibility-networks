"""Spatial Transition Score (STS): a computable, publicly-available
alternative to ESG ratings for measuring transition risk exposure.

For each firm j relative to a retirement event i:

    STS_j = w_fuel_ij * has_ets_j - w_geo_ij

A more negative STS means more transition risk (stranding dominates benefit).
All inputs are publicly available (GEM plant data, ETS membership, GPS coordinates).

Part 1 -- In-sample validation: STS vs ESG in predicting CARs around
    coal retirement events. Portfolio sorts on STS quintiles.

Part 2 -- Out-of-sample validation: temporal split (pre/post 2020).
    Estimate on training set, predict on test set.
    Compare MSPE for STS, ESG, and naive (mean CAR) models.

Part 3 -- Firm-level STS ranking: average exposure across all events.
"""
import csv
import os
import sys
import math
import random
import hashlib
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

# ── Configuration ────────────────────────────────────────────────────

POST_MONTHS = 3   # [-1, +3] window
PRE_MONTHS = 24   # pre-event months for AR demeaning


def _print(msg=''):
    print(msg)
    sys.stdout.flush()


# ── Matrix utilities ────────────────────────────────────────────────

def invert_matrix(mat):
    """Gauss-Jordan inversion of a square matrix."""
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
    """Multiply two matrices."""
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


# ── OLS with event-clustered SEs ─────────────────────────────────────

def _cluster_cov(X, resid, cluster_map, k):
    """Compute clustered meat matrix. cluster_map: {cid: [indices]}."""
    S = [[0.0] * k for _ in range(k)]
    for idxs in cluster_map.values():
        xu = [0.0] * k
        for i in idxs:
            ri = resid[i]
            Xi = X[i]
            for a in range(k):
                xu[a] += Xi[a] * ri
        for a in range(k):
            xua = xu[a]
            for b_idx in range(a, k):
                v = xua * xu[b_idx]
                S[a][b_idx] += v
                if a != b_idx:
                    S[b_idx][a] += v
    return S


def ols_full(data, y_var, x_vars, cluster_var=None):
    """OLS regression returning betas, cluster-robust SEs, t-stats, R2,
    and the FULL variance-covariance matrix V for inference.

    Returns dict with keys: beta, se, t, r2, n, V (full vcov matrix),
    ss_res, X, resid, y, y_hat, inv_XtX, clusters.
    """
    n = len(data)
    k = len(x_vars) + 1
    if n <= k + 1:
        return None

    y = [d[y_var] for d in data]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    if ss_tot < 1e-15:
        return None

    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]

    XtX = [[sum(X[i][a] * X[i][b_idx] for i in range(n))
            for b_idx in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    beta = [sum(inv_XtX[a][b_idx] * Xty[b_idx] for b_idx in range(k))
            for a in range(k)]

    y_hat = [sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]
    ss_res = sum(r ** 2 for r in resid)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0

    V = None
    G = None
    if cluster_var:
        cluster_map = {}
        for i, d in enumerate(data):
            cid = d.get(cluster_var, None)
            cluster_map.setdefault(cid, []).append(i)
        G = len(cluster_map)
        S = _cluster_cov(X, resid, cluster_map, k)
        V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
        if G > 1:
            scale = (G / (G - 1)) * ((n - 1) / (n - k))
            for a in range(k):
                for b_idx in range(k):
                    V[a][b_idx] *= scale
        se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
    else:
        s2 = ss_res / (n - k) if n > k else 0
        V = [[s2 * inv_XtX[a][b_idx] for b_idx in range(k)] for a in range(k)]
        se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]

    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0 for a in range(k)]

    names = ['intercept'] + x_vars
    return {
        'beta': dict(zip(names, beta)),
        'se': dict(zip(names, se)),
        't': dict(zip(names, t_stats)),
        'r2': r2,
        'n': n,
        'V': V,
        'ss_res': ss_res,
        'X': X,
        'resid': resid,
        'y': y,
        'y_hat': y_hat,
        'inv_XtX': inv_XtX,
        'clusters': G,
    }


# ── p-value from t-statistic (two-sided, normal approximation) ──────

def _normal_cdf(x):
    """Standard normal CDF approximation (Abramowitz & Stegun 26.2.17)."""
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
    """Two-sided p-value from t-stat using normal CDF approximation."""
    return 2.0 * (1.0 - _normal_cdf(abs(t_stat)))


def sig_stars(p):
    """Return significance stars for a p-value."""
    if p < 0.01:
        return '***'
    if p < 0.05:
        return '**'
    if p < 0.10:
        return '*'
    return ''


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

# Fama-French monthly factors (vwretd = Mkt-RF + RF)
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

# Weight matrices
_print('Loading weight matrices...')
W_geo = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_geo[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
_print(f'  W_geo firms: {len(W_geo)}')

W_fuel = defaultdict(dict)
fuel_path = derived_path('networks', 'weight_matrix_W_fuel.csv')
if os.path.exists(fuel_path):
    with open(fuel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            W_fuel[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
    _print(f'  W_fuel edges: {sum(len(v) for v in W_fuel.values())}')
else:
    _print('  W_fuel: NOT FOUND')

# Fundamentals (latest record per firm for sector classification + names)
_print('Loading fundamentals...')
fundamentals = {}
fundamentals_by_year = defaultdict(dict)
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        fundamentals_by_year[gk][fy] = row
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row
_print(f'  Fundamentals: {len(fundamentals)} firms')


def get_sic4(gvkey):
    f = fundamentals.get(gvkey)
    if f and f.get('sic'):
        return f['sic'][:4]
    return None


def get_firm_name(gvkey):
    f = fundamentals.get(gvkey)
    if f and f.get('conm'):
        return f['conm']
    return gvkey


def get_country(gvkey):
    f = fundamentals.get(gvkey)
    if f and f.get('fic'):
        return f['fic']
    return '—'


# ESG scores (Refinitiv Environmental Score)
_print('Loading ESG scores...')
esg_scores = {}
esg_path = raw_path('refinitiv', 'refinitiv_esg.csv')
if os.path.exists(esg_path):
    with open(esg_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row.get('gvkey', '').strip()
            env = row.get('env_score', '').strip()
            if gk and env:
                try:
                    esg_scores[gk] = float(env) / 100.0  # normalize to [0,1]
                except ValueError:
                    pass
    _print(f'  ESG scores: {len(esg_scores)} firms (normalized to [0,1])')
else:
    raise RuntimeError(f'Missing ESG file: {esg_path}')

# ETS membership
_print('Loading ETS membership...')
ets_membership = {}
ets_path = derived_path('networks', 'firm_ets_membership.csv')
if os.path.exists(ets_path):
    with open(ets_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row.get('gvkey', '').strip()
            has = row.get('has_ets', '0').strip()
            ets_membership[gk] = 1.0 if has == '1' else 0.0
    _print(f'  ETS membership: {len(ets_membership)} firms '
           f'({sum(1 for v in ets_membership.values() if v == 1.0)} with ETS)')
else:
    _print('  ETS membership: NOT FOUND; defaulting all to 0')

# Events (first-mover coal retirements only)
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


# ── CAR computation ──────────────────────────────────────────────────

def compute_monthly_car(gvkey, event_month, post=3):
    """Monthly CAR[-1, +post] using vwretd model with pre-event demeaning."""
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

    # Require enough pre-event data
    pre_rets = [monthly_ret[gvkey][months[i]]
                for i in range(max(0, event_idx - PRE_MONTHS), event_idx)
                if months[i] in monthly_ret[gvkey]]
    if len(pre_rets) < 12:
        return None

    # Pre-demean ARs by pre-window mean
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


# ── Build regression dataset ──────────────────────────────────────────

def build_obs():
    """Build cross-sectional regression dataset for the [-1,+3] window.
    Computes the STS for each event-firm pair and includes ESG where available.
    Returns (all_obs, esg_obs) where esg_obs is the subset with ESG scores."""
    all_obs = []
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

        # Get first-mover SIC4
        fm_sic4 = None
        for gk in event_gvkeys:
            fm_sic4 = get_sic4(gk)
            if fm_sic4:
                break

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
                has_ets = ets_membership.get(gk, 0.0)

                # STS = w_fuel * has_ets - w_geo
                sts = w_fuel * has_ets - w_geo

                car = compute_monthly_car(gk, event_month, post=POST_MONTHS)
                if car is None:
                    continue

                rec = {
                    'car': car,
                    'sts': sts,
                    'w_geo': w_geo,
                    'w_fuel': w_fuel,
                    'has_ets': has_ets,
                    'event_id': event_id,
                    'event_date': event_date,
                    'gvkey': gk,
                }

                # Add ESG if available
                if gk in esg_scores:
                    rec['esg_score'] = esg_scores[gk]

                all_obs.append(rec)

    esg_obs = [o for o in all_obs if 'esg_score' in o]
    return all_obs, esg_obs


# ── Correlation helper ────────────────────────────────────────────────

def correlation(xs, ys):
    """Pearson correlation between two lists."""
    n = len(xs)
    if n < 3:
        return 0.0
    mx = sum(xs) / n
    my = sum(ys) / n
    cov = sum((xs[i] - mx) * (ys[i] - my) for i in range(n))
    sx = math.sqrt(sum((x - mx) ** 2 for x in xs))
    sy = math.sqrt(sum((y - my) ** 2 for y in ys))
    if sx < 1e-15 or sy < 1e-15:
        return 0.0
    return cov / (sx * sy)


# ── Main analysis ────────────────────────────────────────────────────

_print()
_print('=' * 70)
_print('SPATIAL TRANSITION SCORE (STS)')
_print('STS_j = w_fuel_ij * has_ets_j - w_geo_ij')
_print(f'Window: [-1, +{POST_MONTHS}] months')
_print('=' * 70)

# Build dataset
_print('\nBuilding dataset...')
all_obs, esg_obs = build_obs()
n_all = len(all_obs)
n_esg = len(esg_obs)
n_events_all = len(set(o['event_id'] for o in all_obs))
n_events_esg = len(set(o['event_id'] for o in esg_obs))
n_firms_all = len(set(o['gvkey'] for o in all_obs))
n_firms_esg = len(set(o['gvkey'] for o in esg_obs))
_print(f'  All obs: N = {n_all}, events = {n_events_all}, '
       f'firms = {n_firms_all}')
_print(f'  With ESG: N = {n_esg}, events = {n_events_esg}, '
       f'firms = {n_firms_esg}')

# STS summary stats
sts_vals = [o['sts'] for o in all_obs]
sts_mean = sum(sts_vals) / len(sts_vals) if sts_vals else 0
sts_min = min(sts_vals) if sts_vals else 0
sts_max = max(sts_vals) if sts_vals else 0
sts_sorted = sorted(sts_vals)
sts_median = sts_sorted[len(sts_sorted) // 2] if sts_sorted else 0
_print(f'  STS: mean={sts_mean:.4f}, median={sts_median:.4f}, '
       f'min={sts_min:.4f}, max={sts_max:.4f}')

# ══════════════════════════════════════════════════════════════════════
# PART 1: In-Sample Validation
# ══════════════════════════════════════════════════════════════════════

_print()
_print('=' * 70)
_print('PART 1: IN-SAMPLE VALIDATION')
_print('=' * 70)

# --- STS regression (all obs) ---
_print('\nRegression: CAR = a + b*STS (all observations)...')
res_sts_all = ols_full(all_obs, 'car', ['sts'], cluster_var='event_id')
if res_sts_all:
    _print(f'  N = {res_sts_all["n"]}, R2 = {res_sts_all["r2"]:.6f}')
    b_sts = res_sts_all['beta']['sts']
    se_sts = res_sts_all['se']['sts']
    t_sts = res_sts_all['t']['sts']
    p_sts = p_from_t(t_sts)
    _print(f'  STS: beta={b_sts:+.6f}, se={se_sts:.6f}, '
           f't={t_sts:.3f}, p={p_sts:.4f}{sig_stars(p_sts)}')
else:
    _print('  ERROR: STS regression failed.')

# --- ESG regression (ESG subsample) ---
_print('\nRegression: CAR = a + b*ESG (ESG subsample)...')
res_esg_sub = ols_full(esg_obs, 'car', ['esg_score'], cluster_var='event_id')
if res_esg_sub:
    _print(f'  N = {res_esg_sub["n"]}, R2 = {res_esg_sub["r2"]:.6f}')
    b_esg = res_esg_sub['beta']['esg_score']
    se_esg = res_esg_sub['se']['esg_score']
    t_esg = res_esg_sub['t']['esg_score']
    p_esg = p_from_t(t_esg)
    _print(f'  ESG: beta={b_esg:+.6f}, se={se_esg:.6f}, '
           f't={t_esg:.3f}, p={p_esg:.4f}{sig_stars(p_esg)}')
else:
    _print('  ERROR: ESG regression failed.')

# --- STS regression on ESG subsample (for apples-to-apples comparison) ---
_print('\nRegression: CAR = a + b*STS (ESG subsample, for comparison)...')
res_sts_sub = ols_full(esg_obs, 'car', ['sts'], cluster_var='event_id')
if res_sts_sub:
    _print(f'  N = {res_sts_sub["n"]}, R2 = {res_sts_sub["r2"]:.6f}')
    b_sts_sub = res_sts_sub['beta']['sts']
    se_sts_sub = res_sts_sub['se']['sts']
    t_sts_sub = res_sts_sub['t']['sts']
    p_sts_sub = p_from_t(t_sts_sub)
    _print(f'  STS: beta={b_sts_sub:+.6f}, se={se_sts_sub:.6f}, '
           f't={t_sts_sub:.3f}, p={p_sts_sub:.4f}{sig_stars(p_sts_sub)}')
else:
    _print('  ERROR: STS regression on ESG subsample failed.')

# --- Portfolio sorts on STS quintiles ---
_print('\nPortfolio sorts on STS quintiles...')
obs_sorted_sts = sorted(all_obs, key=lambda o: o['sts'])
n_per_q = len(obs_sorted_sts) // 5
quintile_results = []
quintile_cars = {}

for q in range(5):
    start = q * n_per_q
    end = start + n_per_q if q < 4 else len(obs_sorted_sts)
    q_obs = obs_sorted_sts[start:end]
    cars = [o['car'] for o in q_obs]
    mean_car = sum(cars) / len(cars)
    mean_sts = sum(o['sts'] for o in q_obs) / len(q_obs)
    quintile_results.append({
        'quintile': q + 1,
        'n': len(q_obs),
        'mean_car': mean_car,
        'mean_sts': mean_sts,
    })
    quintile_cars[q + 1] = cars
    _print(f'  Q{q+1}: N={len(q_obs)}, mean_STS={mean_sts:.4f}, '
           f'mean_CAR={mean_car:.4f}')

# Q5-Q1 spread and t-stat
q1_cars = quintile_cars[1]
q5_cars = quintile_cars[5]
mean_q1 = sum(q1_cars) / len(q1_cars)
mean_q5 = sum(q5_cars) / len(q5_cars)
spread = mean_q5 - mean_q1

# t-stat for Q5-Q1 spread (unequal variance t-test)
var_q1 = sum((c - mean_q1) ** 2 for c in q1_cars) / (len(q1_cars) - 1) if len(q1_cars) > 1 else 0
var_q5 = sum((c - mean_q5) ** 2 for c in q5_cars) / (len(q5_cars) - 1) if len(q5_cars) > 1 else 0
se_spread = math.sqrt(var_q1 / len(q1_cars) + var_q5 / len(q5_cars)) if (var_q1 + var_q5) > 0 else 0
t_spread = spread / se_spread if se_spread > 1e-15 else 0
p_spread = p_from_t(t_spread)
_print(f'  Q5-Q1 spread: {spread:+.4f} (t = {t_spread:.3f}, p = {p_spread:.4f})')


# ══════════════════════════════════════════════════════════════════════
# PART 2: Out-of-Sample Validation (temporal split at 2020-01-01)
# ══════════════════════════════════════════════════════════════════════

_print()
_print('=' * 70)
_print('PART 2: OUT-OF-SAMPLE VALIDATION')
_print('Split: training = events before 2020, test = events 2020+')
_print('=' * 70)

SPLIT_DATE = '2020-01'

# Split all_obs by event date
train_obs = [o for o in all_obs if o['event_date'][:7] < SPLIT_DATE]
test_obs = [o for o in all_obs if o['event_date'][:7] >= SPLIT_DATE]

# Also split ESG subsample
train_esg = [o for o in esg_obs if o['event_date'][:7] < SPLIT_DATE]
test_esg = [o for o in esg_obs if o['event_date'][:7] >= SPLIT_DATE]

n_train = len(train_obs)
n_test = len(test_obs)
n_train_esg = len(train_esg)
n_test_esg = len(test_esg)
n_train_events = len(set(o['event_id'] for o in train_obs))
n_test_events = len(set(o['event_id'] for o in test_obs))

_print(f'  Training: N = {n_train}, events = {n_train_events}')
_print(f'  Test:     N = {n_test}, events = {n_test_events}')
_print(f'  Training (ESG sub): N = {n_train_esg}')
_print(f'  Test (ESG sub):     N = {n_test_esg}')

# --- STS: estimate on training, predict on test ---
_print('\nSTS model: estimate on training...')
res_sts_train = ols_full(train_obs, 'car', ['sts'], cluster_var='event_id')

sts_mspe = None
sts_corr = None
sts_dir_acc = None
sts_predicted = []
sts_actual = []

if res_sts_train and n_test > 0:
    alpha_sts = res_sts_train['beta']['intercept']
    beta_sts_train = res_sts_train['beta']['sts']
    _print(f'  Training: alpha={alpha_sts:+.6f}, beta={beta_sts_train:+.6f}, '
           f'R2={res_sts_train["r2"]:.6f}')

    # Predict on test set
    for o in test_obs:
        pred = alpha_sts + beta_sts_train * o['sts']
        actual = o['car']
        sts_predicted.append(pred)
        sts_actual.append(actual)

    sts_mspe = sum((sts_predicted[i] - sts_actual[i]) ** 2
                   for i in range(len(sts_actual))) / len(sts_actual)
    sts_corr = correlation(sts_predicted, sts_actual)
    sts_dir_acc = sum(1 for i in range(len(sts_actual))
                      if (sts_predicted[i] >= 0) == (sts_actual[i] >= 0)) / len(sts_actual)
    _print(f'  Test MSPE:  {sts_mspe:.8f}')
    _print(f'  Test corr:  {sts_corr:.4f}')
    _print(f'  Dir. acc.:  {sts_dir_acc:.4f}')
else:
    _print('  ERROR: STS training regression failed or no test observations.')

# --- ESG: estimate on training, predict on test ---
_print('\nESG model: estimate on training (ESG subsample)...')
res_esg_train = ols_full(train_esg, 'car', ['esg_score'], cluster_var='event_id')

esg_mspe = None
esg_corr = None
esg_dir_acc = None
esg_predicted = []
esg_actual = []

if res_esg_train and n_test_esg > 0:
    alpha_esg_train = res_esg_train['beta']['intercept']
    beta_esg_train = res_esg_train['beta']['esg_score']
    _print(f'  Training: alpha={alpha_esg_train:+.6f}, '
           f'beta={beta_esg_train:+.6f}, R2={res_esg_train["r2"]:.6f}')

    for o in test_esg:
        pred = alpha_esg_train + beta_esg_train * o['esg_score']
        actual = o['car']
        esg_predicted.append(pred)
        esg_actual.append(actual)

    esg_mspe = sum((esg_predicted[i] - esg_actual[i]) ** 2
                   for i in range(len(esg_actual))) / len(esg_actual)
    esg_corr = correlation(esg_predicted, esg_actual)
    esg_dir_acc = sum(1 for i in range(len(esg_actual))
                      if (esg_predicted[i] >= 0) == (esg_actual[i] >= 0)) / len(esg_actual)
    _print(f'  Test MSPE:  {esg_mspe:.8f}')
    _print(f'  Test corr:  {esg_corr:.4f}')
    _print(f'  Dir. acc.:  {esg_dir_acc:.4f}')
else:
    _print('  ERROR: ESG training regression failed or no test observations.')

# --- Naive model: predict mean CAR from training set ---
_print('\nNaive model (mean CAR from training)...')
train_cars = [o['car'] for o in train_obs]
naive_mean = sum(train_cars) / len(train_cars) if train_cars else 0

naive_mspe = None
naive_dir_acc = None
if n_test > 0:
    naive_mspe = sum((naive_mean - o['car']) ** 2 for o in test_obs) / n_test
    naive_dir_acc = sum(1 for o in test_obs
                        if (naive_mean >= 0) == (o['car'] >= 0)) / n_test
    _print(f'  Training mean CAR: {naive_mean:+.6f}')
    _print(f'  Test MSPE:  {naive_mspe:.8f}')
    _print(f'  Dir. acc.:  {naive_dir_acc:.4f}')

# --- Also compute naive on ESG subsample for fair comparison ---
train_esg_cars = [o['car'] for o in train_esg]
naive_esg_mean = sum(train_esg_cars) / len(train_esg_cars) if train_esg_cars else 0
naive_esg_mspe = None
naive_esg_dir_acc = None
if n_test_esg > 0:
    naive_esg_mspe = sum((naive_esg_mean - o['car']) ** 2
                         for o in test_esg) / n_test_esg
    naive_esg_dir_acc = sum(1 for o in test_esg
                            if (naive_esg_mean >= 0) == (o['car'] >= 0)) / n_test_esg


# ══════════════════════════════════════════════════════════════════════
# PART 3: Firm-Level STS Ranking
# ══════════════════════════════════════════════════════════════════════

_print()
_print('=' * 70)
_print('PART 3: FIRM-LEVEL STS RANKING')
_print('=' * 70)

# Average STS across all events for each firm
firm_sts = defaultdict(list)
for o in all_obs:
    firm_sts[o['gvkey']].append(o['sts'])

firm_avg_sts = []
for gk, sts_list in firm_sts.items():
    avg = sum(sts_list) / len(sts_list)
    firm_avg_sts.append({
        'gvkey': gk,
        'avg_sts': avg,
        'n_events': len(sts_list),
        'name': get_firm_name(gk),
        'country': get_country(gk),
    })

# Sort: most negative first (most exposed)
firm_avg_sts.sort(key=lambda x: x['avg_sts'])

_print(f'\n  Total firms with STS: {len(firm_avg_sts)}')

_print('\n  Most exposed (most negative average STS):')
for i, f in enumerate(firm_avg_sts[:10]):
    _print(f'    {i+1:2d}. {f["name"]:<30s} ({f["country"]})  '
           f'avg_STS={f["avg_sts"]:+.4f}  N={f["n_events"]}')

_print('\n  Least exposed (most positive average STS):')
least_exposed = list(reversed(firm_avg_sts[-10:]))
for i, f in enumerate(least_exposed):
    _print(f'    {i+1:2d}. {f["name"]:<30s} ({f["country"]})  '
           f'avg_STS={f["avg_sts"]:+.4f}  N={f["n_events"]}')


# ══════════════════════════════════════════════════════════════════════
# Write output
# ══════════════════════════════════════════════════════════════════════

_print('\nWriting results...')
out_path = results_path('metrics', 'strategy2_spatial_score.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

lines = []
lines.append('# Spatial Transition Score (STS): A Computable Alternative to ESG')
lines.append('')
lines.append('## Construction')
lines.append('')
lines.append('STS_j = w_fuel_ij * has_ets_j - w_geo_ij')
lines.append('')
lines.append('Inputs: GEM plant-level data (public), ETS membership (public), '
             'GPS coordinates (public).')
lines.append('No corporate disclosure required.')
lines.append('')

# ── Part 1 ──

lines.append('## Part 1: In-Sample Validation')
lines.append('')
lines.append(f'Full sample: N = {n_all}, events = {n_events_all}, '
             f'firms = {n_firms_all}.')
lines.append(f'ESG subsample: N = {n_esg}, events = {n_events_esg}, '
             f'firms = {n_firms_esg}.')
lines.append(f'Standard errors: event-clustered.')
lines.append(f'Window: [-1, +{POST_MONTHS}] months.')
lines.append('')

lines.append('### Univariate regressions: CAR = a + b * Predictor')
lines.append('')
lines.append('| Predictor | Beta | SE | t-stat | R-squared | N |')
lines.append('|---|---:|---:|---:|---:|---:|')

if res_sts_all:
    lines.append(f'| STS (full sample) | '
                 f'{res_sts_all["beta"]["sts"]:+.6f} | '
                 f'{res_sts_all["se"]["sts"]:.6f} | '
                 f'{res_sts_all["t"]["sts"]:.3f} | '
                 f'{res_sts_all["r2"]:.6f} | '
                 f'{res_sts_all["n"]} |')
if res_sts_sub:
    lines.append(f'| STS (ESG subsample) | '
                 f'{res_sts_sub["beta"]["sts"]:+.6f} | '
                 f'{res_sts_sub["se"]["sts"]:.6f} | '
                 f'{res_sts_sub["t"]["sts"]:.3f} | '
                 f'{res_sts_sub["r2"]:.6f} | '
                 f'{res_sts_sub["n"]} |')
if res_esg_sub:
    lines.append(f'| ESG score | '
                 f'{res_esg_sub["beta"]["esg_score"]:+.6f} | '
                 f'{res_esg_sub["se"]["esg_score"]:.6f} | '
                 f'{res_esg_sub["t"]["esg_score"]:.3f} | '
                 f'{res_esg_sub["r2"]:.6f} | '
                 f'{res_esg_sub["n"]} |')
lines.append('')

# Portfolio sorts
lines.append('### STS Portfolio Sorts')
lines.append('')
lines.append('| Quintile | Mean STS | Mean CAR | N_firms |')
lines.append('|---|---:|---:|---:|')
for qr in quintile_results:
    lines.append(f'| Q{qr["quintile"]} | {qr["mean_sts"]:+.4f} | '
                 f'{qr["mean_car"]:+.4f} | {qr["n"]} |')
lines.append(f'| Q5-Q1 | | {spread:+.4f} (t = {t_spread:.3f}{sig_stars(p_spread)}) | |')
lines.append('')

# ── Part 2 ──

lines.append('## Part 2: Out-of-Sample Validation')
lines.append('')
lines.append(f'Training: events before 2020 (N = {n_train}, '
             f'events = {n_train_events})')
lines.append(f'Test: events 2020+ (N = {n_test}, events = {n_test_events})')
lines.append('')

lines.append('### Full sample comparison (STS vs Naive)')
lines.append('')
lines.append('| Metric | STS | Naive (mean) |')
lines.append('|---|---:|---:|')
lines.append(f'| MSPE | {sts_mspe:.8f} | {naive_mspe:.8f} |'
             if sts_mspe is not None and naive_mspe is not None
             else '| MSPE | -- | -- |')
lines.append(f'| Correlation (pred, actual) | {sts_corr:.4f} | -- |'
             if sts_corr is not None else '| Correlation (pred, actual) | -- | -- |')
lines.append(f'| Directional accuracy | {sts_dir_acc:.4f} | {naive_dir_acc:.4f} |'
             if sts_dir_acc is not None and naive_dir_acc is not None
             else '| Directional accuracy | -- | -- |')
lines.append('')

lines.append(f'### ESG subsample comparison (STS vs ESG vs Naive)')
lines.append('')
lines.append('| Metric | STS | ESG | Naive (mean) |')
lines.append('|---|---:|---:|---:|')

# For fair comparison, run STS on ESG subsample split
train_esg_sts = ols_full(train_esg, 'car', ['sts'], cluster_var='event_id')
sts_esg_mspe = None
sts_esg_corr = None
sts_esg_dir_acc = None
if train_esg_sts and n_test_esg > 0:
    a_sts_e = train_esg_sts['beta']['intercept']
    b_sts_e = train_esg_sts['beta']['sts']
    pred_sts_e = [a_sts_e + b_sts_e * o['sts'] for o in test_esg]
    act_sts_e = [o['car'] for o in test_esg]
    sts_esg_mspe = sum((pred_sts_e[i] - act_sts_e[i]) ** 2
                       for i in range(len(act_sts_e))) / len(act_sts_e)
    sts_esg_corr = correlation(pred_sts_e, act_sts_e)
    sts_esg_dir_acc = sum(1 for i in range(len(act_sts_e))
                          if (pred_sts_e[i] >= 0) == (act_sts_e[i] >= 0)) / len(act_sts_e)


def _fmt_val(v, fmt='.8f'):
    return f'{v:{fmt}}' if v is not None else '--'


lines.append(f'| MSPE | {_fmt_val(sts_esg_mspe)} | '
             f'{_fmt_val(esg_mspe)} | {_fmt_val(naive_esg_mspe)} |')
lines.append(f'| Correlation (pred, actual) | {_fmt_val(sts_esg_corr, ".4f")} | '
             f'{_fmt_val(esg_corr, ".4f")} | -- |')
lines.append(f'| Directional accuracy | {_fmt_val(sts_esg_dir_acc, ".4f")} | '
             f'{_fmt_val(esg_dir_acc, ".4f")} | {_fmt_val(naive_esg_dir_acc, ".4f")} |')
lines.append('')

# ── Part 3 ──

lines.append('## Part 3: Firm-Level Exposure Ranking')
lines.append('')

lines.append('### Most Exposed (most negative average STS)')
lines.append('')
lines.append('| Rank | Firm | Country | Avg STS | N events |')
lines.append('|---|---|---|---:|---:|')
for i, f in enumerate(firm_avg_sts[:10]):
    lines.append(f'| {i+1} | {f["name"]} | {f["country"]} | '
                 f'{f["avg_sts"]:+.4f} | {f["n_events"]} |')
lines.append('')

lines.append('### Least Exposed (most positive average STS)')
lines.append('')
lines.append('| Rank | Firm | Country | Avg STS | N events |')
lines.append('|---|---|---|---:|---:|')
for i, f in enumerate(least_exposed):
    lines.append(f'| {i+1} | {f["name"]} | {f["country"]} | '
                 f'{f["avg_sts"]:+.4f} | {f["n_events"]} |')
lines.append('')

# ── Key Finding ──

lines.append('## Key Finding')
lines.append('')

# Build the key finding paragraph based on results
sts_sig_insample = (res_sts_all is not None
                    and abs(res_sts_all['t']['sts']) > 1.96)
esg_sig_insample = (res_esg_sub is not None
                    and abs(res_esg_sub['t']['esg_score']) > 1.96)
sts_beats_esg_oos = (sts_esg_mspe is not None and esg_mspe is not None
                     and sts_esg_mspe < esg_mspe)

if sts_sig_insample and not esg_sig_insample and sts_beats_esg_oos:
    finding = (
        'The Spatial Transition Score is a freely computable measure, constructed '
        'entirely from publicly available plant-level data (GEM trackers, ETS membership '
        'records, GPS coordinates), that outperforms Refinitiv ESG environmental scores '
        'in predicting transition-related repricing around coal retirement events. '
        f'In-sample, STS predicts CARs with a t-statistic of '
        f'{res_sts_all["t"]["sts"]:.2f} '
        f'(R-squared = {res_sts_all["r2"]:.4f}), while the ESG score is not significant '
        f'(t = {res_esg_sub["t"]["esg_score"]:.2f}, '
        f'R-squared = {res_esg_sub["r2"]:.4f}). '
        f'Out-of-sample, STS achieves lower MSPE ({sts_esg_mspe:.6f}) '
        f'than ESG ({esg_mspe:.6f}). Investors can construct this score without '
        'purchasing proprietary ESG ratings, using only publicly observable information '
        'about power plant locations, fuel mix, and carbon pricing jurisdiction.'
    )
elif sts_beats_esg_oos:
    sts_t_str = (f'{res_sts_all["t"]["sts"]:.2f}'
                 if res_sts_all else 'N/A')
    esg_t_str = (f'{res_esg_sub["t"]["esg_score"]:.2f}'
                 if res_esg_sub else 'N/A')
    sts_r2_str = (f'{res_sts_all["r2"]:.4f}'
                  if res_sts_all else 'N/A')
    esg_r2_str = (f'{res_esg_sub["r2"]:.4f}'
                  if res_esg_sub else 'N/A')
    finding = (
        'The Spatial Transition Score, constructed entirely from publicly available '
        'plant-level data, achieves lower out-of-sample prediction error than Refinitiv '
        f'ESG environmental scores (MSPE: {_fmt_val(sts_esg_mspe)} vs '
        f'{_fmt_val(esg_mspe)}). '
        f'In-sample, STS has t = {sts_t_str} (R-squared = {sts_r2_str}), '
        f'while ESG has t = {esg_t_str} (R-squared = {esg_r2_str}). '
        'The score requires no corporate disclosure and can be computed from GEM '
        'plant trackers, ETS membership records, and GPS coordinates alone.'
    )
else:
    sts_t_str = (f'{res_sts_all["t"]["sts"]:.2f}'
                 if res_sts_all else 'N/A')
    esg_t_str = (f'{res_esg_sub["t"]["esg_score"]:.2f}'
                 if res_esg_sub else 'N/A')
    sts_r2_str = (f'{res_sts_all["r2"]:.4f}'
                  if res_sts_all else 'N/A')
    esg_r2_str = (f'{res_esg_sub["r2"]:.4f}'
                  if res_esg_sub else 'N/A')
    finding = (
        'The Spatial Transition Score is a freely computable spatial measure '
        'constructed from publicly available plant-level data. '
        f'In-sample, STS has t = {sts_t_str} (R-squared = {sts_r2_str}), '
        f'while ESG has t = {esg_t_str} (R-squared = {esg_r2_str}). '
        f'Out-of-sample MSPE: STS = {_fmt_val(sts_esg_mspe)}, '
        f'ESG = {_fmt_val(esg_mspe)}, Naive = {_fmt_val(naive_esg_mspe)}. '
        'The score can be constructed using only GEM plant trackers, ETS membership '
        'records, and GPS coordinates, without purchasing proprietary ESG ratings.'
    )

lines.append(finding)
lines.append('')

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

_print(f'\nWrote: {out_path}')
_print('Done.')
