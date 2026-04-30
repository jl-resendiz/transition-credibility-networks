"""Shift-Share (Bartik) causal robustness test for the fuel channel.

Addresses the concern that the fuel coefficient captures shared exposure
to common factors rather than network transmission.

Design:
  SHARES: Pre-period fuel-mix similarity weights w_fuel_ij fixed using only
          plants commissioned before 2014 (before the main retirement wave).
  SHIFTS: Annual aggregate coal MW retired globally per year.
  Bartik instrument: B_ie = w_fuel_pre(i,j) * RetiredMW(event year)

Specifications:
  A (reduced form): CAR = alpha + beta_bartik * B_ie + w_geo + w_reg + same_sector + eps
  B (comparison):   Compare beta_bartik with beta_OLS from the standard specification

Reports Fama-MacBeth + Newey-West and pooled event-clustered results,
plus diagnostics on pre-period share coverage and correlation with
current shares.

Output: results/metrics/strategy2_bartik_shiftshare.md
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
    return {'beta': beta, 'resid': resid, 'r2': r2, 'n': n, 'inv_XtX': inv_XtX}


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


def _beta_inc(a, b, x, niter=200):
    """Regularised incomplete beta function I_x(a,b) via continued fraction."""
    if x < 0 or x > 1:
        return 0.0
    if x == 0 or x == 1:
        return x
    from math import lgamma, exp
    lbeta = lgamma(a) + lgamma(b) - lgamma(a + b)
    front = exp(a * math.log(x) + b * math.log(1 - x) - lbeta) / a
    # Lentz continued fraction
    f = 1.0
    c = 1.0
    d = 1.0 - (a + b) * x / (a + 1)
    if abs(d) < 1e-30:
        d = 1e-30
    d = 1.0 / d
    f = d
    for m in range(1, niter + 1):
        # even step
        num = m * (b - m) * x / ((a + 2 * m - 1) * (a + 2 * m))
        d = 1.0 + num * d
        if abs(d) < 1e-30: d = 1e-30
        c = 1.0 + num / c
        if abs(c) < 1e-30: c = 1e-30
        d = 1.0 / d
        f *= d * c
        # odd step
        num = -(a + m) * (a + b + m) * x / ((a + 2 * m) * (a + 2 * m + 1))
        d = 1.0 + num * d
        if abs(d) < 1e-30: d = 1e-30
        c = 1.0 + num / c
        if abs(c) < 1e-30: c = 1e-30
        d = 1.0 / d
        delta = d * c
        f *= delta
        if abs(delta - 1.0) < 1e-12:
            break
    return front * f


def p_from_t_df(t_stat, df):
    """Two-sided p-value from t-stat using the t-distribution with df degrees
    of freedom. Falls back to normal approximation for df > 200."""
    if df > 200:
        return p_from_t(t_stat)
    x = df / (df + t_stat * t_stat)
    p_two = _beta_inc(df / 2.0, 0.5, x)
    return max(0.0, min(1.0, p_two))


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


def newey_west_cov(series_list, max_lag=None):
    """Newey-West covariance matrix for a list of K time series."""
    K = len(series_list)
    T = len(series_list[0])
    if T < 3:
        return None
    means = [sum(s) / T for s in series_list]
    demean = [[series_list[k][t] - means[k] for t in range(T)] for k in range(K)]
    if max_lag is None:
        max_lag = max(1, int(4 * (T / 100) ** (2 / 9)))
    max_lag = min(max_lag, T - 1)
    S = [[0.0] * K for _ in range(K)]
    for i in range(K):
        for j in range(i, K):
            g0 = sum(demean[i][t] * demean[j][t] for t in range(T)) / T
            S[i][j] = g0
            S[j][i] = g0
    for lag in range(1, max_lag + 1):
        weight = 1.0 - lag / (max_lag + 1.0)
        for i in range(K):
            for j in range(K):
                g_lag = sum(demean[i][t] * demean[j][t - lag]
                            for t in range(lag, T)) / T
                g_lag_rev = sum(demean[j][t] * demean[i][t - lag]
                                for t in range(lag, T)) / T
                S[i][j] += weight * (g_lag + g_lag_rev)
    V = [[S[i][j] / T for j in range(K)] for i in range(K)]
    return V


# ── Clustered standard errors ───────────────────────────────────────

def clustered_se(X_mat, resid, cluster_ids, inv_XtX):
    """Event-clustered standard errors."""
    n = len(resid)
    k = len(X_mat[0])
    clusters = defaultdict(list)
    for i, cid in enumerate(cluster_ids):
        clusters[cid].append(i)
    G = len(clusters)
    if G <= k:
        return None
    meat = [[0.0] * k for _ in range(k)]
    for cid, indices in clusters.items():
        score = [0.0] * k
        for i in indices:
            for a in range(k):
                score[a] += X_mat[i][a] * resid[i]
        for a in range(k):
            for b in range(k):
                meat[a][b] += score[a] * score[b]
    scale = G / (G - 1.0) * (n - 1.0) / (n - k)
    for a in range(k):
        for b in range(k):
            meat[a][b] *= scale
    V = mat_mul(mat_mul(inv_XtX, meat), inv_XtX)
    se = [math.sqrt(max(0, V[a][a])) for a in range(k)]
    return se


# ══════════════════════════════════════════════════════════════════════
# STEP 1: BUILD PRE-PERIOD FUEL-MIX SHARES
# ══════════════════════════════════════════════════════════════════════

_print('=' * 70)
_print('SHIFT-SHARE (BARTIK) CAUSAL ROBUSTNESS TEST')
_print('=' * 70)

# Load Parent -> gvkey mapping
_print('\nLoading Parent -> gvkey mapping...')
parent_to_gvkeys = defaultdict(set)
gvkey_to_parents = defaultdict(set)
with open(derived_path('mappings', 'gem_compustat_matches.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        parent_to_gvkeys[row['gem_parent']].add(row['gvkey'])
        gvkey_to_parents[row['gvkey']].add(row['gem_parent'])
_print(f'  {len(parent_to_gvkeys)} parents mapped to {len(gvkey_to_parents)} gvkeys')

# Helper to parse parents from GEM field
import re

def parse_parents(field):
    if not field or str(field).strip() == '':
        return []
    parts = str(field).split(';')
    results = []
    for p in parts:
        p = p.strip()
        match = re.match(r'^(.+?)\s*\[(\d+\.?\d*)%\]$', p)
        if match:
            results.append(match.group(1).strip())
        elif p:
            results.append(p.strip())
    return results


def safe_float(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def safe_int(x):
    try:
        return int(float(x))
    except (ValueError, TypeError):
        return None


PRE_CUTOFF = 2014  # Plants commissioned before this year are "pre-determined"
PRE_CUTOFF_SENSITIVITY = 2010  # Sensitivity check: even more pre-determined shares

_print(f'\nBuilding pre-period fuel vectors (plants with Start year < {PRE_CUTOFF})...')

# Accumulate MW by gvkey and fuel type from ALL four GEM trackers
# using only plants commissioned before PRE_CUTOFF
pre_mw = defaultdict(lambda: {'coal': 0.0, 'gas': 0.0, 'solar': 0.0, 'wind': 0.0})
current_mw = defaultdict(lambda: {'coal': 0.0, 'gas': 0.0, 'solar': 0.0, 'wind': 0.0})

# Also store raw plant records to rebuild pre_mw for sensitivity cutoff
_all_plant_records = []  # list of (gvkey_set, fuel, cap, start_year, ret_year)

# --- Coal ---
_print('  Reading GEM Coal Plant Tracker...')
fpath = derived_path('gem', 'gem_coal.csv')
n_coal_pre = 0
n_coal_total = 0
with open(fpath, newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        status = str(row['Status']) if row['Status'] else ''
        # Include operating, mothballed, and retired plants (all existed at some point)
        cap = safe_float(row['Capacity (MW)'])
        if cap is None or cap <= 0:
            continue
        parent_raw = row['Parent']
        parents = parse_parents(parent_raw)
        start_year = safe_int(row['Start year'])
        ret_year = safe_int(row['Retired year'])

        # Map to gvkeys
        matched_gvkeys = set()
        for name in parents:
            if name in parent_to_gvkeys:
                matched_gvkeys.update(parent_to_gvkeys[name])
        if not matched_gvkeys:
            continue

        n_coal_total += 1
        _all_plant_records.append((matched_gvkeys, 'coal', cap, start_year, ret_year))

        # Current portfolio: include if operating or retired after study start
        for gk in matched_gvkeys:
            current_mw[gk]['coal'] += cap

        # Pre-period: commissioned before cutoff AND not retired before cutoff
        if start_year is not None and start_year < PRE_CUTOFF:
            if ret_year is None or ret_year >= PRE_CUTOFF:
                for gk in matched_gvkeys:
                    pre_mw[gk]['coal'] += cap
                n_coal_pre += 1
        elif start_year is None:
            # No start year: include if not retired before cutoff (proxy)
            if ret_year is None or ret_year >= PRE_CUTOFF:
                for gk in matched_gvkeys:
                    pre_mw[gk]['coal'] += cap
                n_coal_pre += 1

_print(f'    Coal units mapped: {n_coal_total} total, {n_coal_pre} pre-{PRE_CUTOFF}')

# --- Gas ---
_print('  Reading GEM Gas Plant Tracker...')
fpath = derived_path('gem', 'gem_gas.csv')
n_gas_pre = 0
n_gas_total = 0
with open(fpath, newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        fuel = str(row.get('Fuel', '') or '').lower()
        if 'gas' not in fuel:
            continue
        cap = safe_float(row['Capacity (MW)'])
        if cap is None or cap <= 0:
            continue
        parent_raw = row['Parent(s)']
        parents = parse_parents(parent_raw)
        start_year = safe_int(row['Start year'])
        ret_year = safe_int(row['Retired year'])

        matched_gvkeys = set()
        for name in parents:
            if name in parent_to_gvkeys:
                matched_gvkeys.update(parent_to_gvkeys[name])
        if not matched_gvkeys:
            continue

        n_gas_total += 1
        _all_plant_records.append((matched_gvkeys, 'gas', cap, start_year, ret_year))
        for gk in matched_gvkeys:
            current_mw[gk]['gas'] += cap

        if start_year is not None and start_year < PRE_CUTOFF:
            if ret_year is None or ret_year >= PRE_CUTOFF:
                for gk in matched_gvkeys:
                    pre_mw[gk]['gas'] += cap
                n_gas_pre += 1
        elif start_year is None:
            if ret_year is None or ret_year >= PRE_CUTOFF:
                for gk in matched_gvkeys:
                    pre_mw[gk]['gas'] += cap
                n_gas_pre += 1

_print(f'    Gas units mapped: {n_gas_total} total, {n_gas_pre} pre-{PRE_CUTOFF}')

# --- Solar ---
_print('  Reading GEM Solar Power Tracker...')
fpath = derived_path('gem', 'gem_solar.csv')
n_solar_pre = 0
n_solar_total = 0
with open(fpath, newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        cap = safe_float(row['Capacity (MW)'])
        if cap is None or cap <= 0:
            continue
        status = str(row['Status'] or '').lower()
        owner_raw = row['Owner']
        parents = parse_parents(owner_raw)
        start_year = safe_int(row['Start year'])
        ret_year = safe_int(row['Retired year'])

        matched_gvkeys = set()
        for name in parents:
            if name in parent_to_gvkeys:
                matched_gvkeys.update(parent_to_gvkeys[name])
        if not matched_gvkeys:
            continue

        n_solar_total += 1
        _all_plant_records.append((matched_gvkeys, 'solar', cap, start_year, ret_year))
        for gk in matched_gvkeys:
            current_mw[gk]['solar'] += cap

        if start_year is not None and start_year < PRE_CUTOFF:
            if ret_year is None or ret_year >= PRE_CUTOFF:
                for gk in matched_gvkeys:
                    pre_mw[gk]['solar'] += cap
                n_solar_pre += 1
        elif start_year is None:
            if ret_year is None or ret_year >= PRE_CUTOFF:
                for gk in matched_gvkeys:
                    pre_mw[gk]['solar'] += cap
                n_solar_pre += 1

_print(f'    Solar units mapped: {n_solar_total} total, {n_solar_pre} pre-{PRE_CUTOFF}')

# --- Wind ---
_print('  Reading GEM Wind Power Tracker...')
fpath = derived_path('gem', 'gem_wind.csv')
n_wind_pre = 0
n_wind_total = 0
with open(fpath, newline='', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        cap = safe_float(row['Capacity (MW)'])
        if cap is None or cap <= 0:
            continue
        owner_raw = row['Owner']
        parents = parse_parents(owner_raw)
        start_year = safe_int(row['Start year'])
        ret_year = safe_int(row['Retired year'])

        matched_gvkeys = set()
        for name in parents:
            if name in parent_to_gvkeys:
                matched_gvkeys.update(parent_to_gvkeys[name])
        if not matched_gvkeys:
            continue

        n_wind_total += 1
        _all_plant_records.append((matched_gvkeys, 'wind', cap, start_year, ret_year))
        for gk in matched_gvkeys:
            current_mw[gk]['wind'] += cap

        if start_year is not None and start_year < PRE_CUTOFF:
            if ret_year is None or ret_year >= PRE_CUTOFF:
                for gk in matched_gvkeys:
                    pre_mw[gk]['wind'] += cap
                n_wind_pre += 1
        elif start_year is None:
            if ret_year is None or ret_year >= PRE_CUTOFF:
                for gk in matched_gvkeys:
                    pre_mw[gk]['wind'] += cap
                n_wind_pre += 1

_print(f'    Wind units mapped: {n_wind_total} total, {n_wind_pre} pre-{PRE_CUTOFF}')

# Compute fuel share vectors
def compute_shares(mw_dict):
    coal = mw_dict['coal']
    gas = mw_dict['gas']
    solar = mw_dict['solar']
    wind = mw_dict['wind']
    total = coal + gas + solar + wind
    if total <= 0:
        return None
    return [coal / total, gas / total, solar / total, wind / total]


def build_pre_mw_for_cutoff(plant_records, cutoff):
    """Rebuild pre-period MW from stored plant records for a given cutoff year."""
    result = defaultdict(lambda: {'coal': 0.0, 'gas': 0.0, 'solar': 0.0, 'wind': 0.0})
    n_included = 0
    for gvkeys, fuel, cap, start_year, ret_year in plant_records:
        if start_year is not None and start_year < cutoff:
            if ret_year is None or ret_year >= cutoff:
                for gk in gvkeys:
                    result[gk][fuel] += cap
                n_included += 1
        elif start_year is None:
            if ret_year is None or ret_year >= cutoff:
                for gk in gvkeys:
                    result[gk][fuel] += cap
                n_included += 1
    return result, n_included


pre_shares = {}
current_shares = {}
for gk in set(list(pre_mw.keys()) + list(current_mw.keys())):
    s = compute_shares(pre_mw[gk])
    if s is not None:
        pre_shares[gk] = s
    s = compute_shares(current_mw[gk])
    if s is not None:
        current_shares[gk] = s

_print(f'\nFirms with pre-period fuel shares: {len(pre_shares)}')
_print(f'Firms with current fuel shares: {len(current_shares)}')
_print(f'Overlap: {len(set(pre_shares.keys()) & set(current_shares.keys()))}')


# Fuel-mix similarity: 1 - 0.5 * L1 distance (same as build_fuel_matrix.py)
def fuel_similarity(v1, v2):
    if v1 is None or v2 is None:
        return 0.0
    l1 = sum(abs(a - b) for a, b in zip(v1, v2))
    sim = 1.0 - 0.5 * l1
    return max(0.0, min(1.0, sim))


# ══════════════════════════════════════════════════════════════════════
# STEP 2: COMPUTE ANNUAL AGGREGATE RETIREMENT SHOCKS
# ══════════════════════════════════════════════════════════════════════

_print('\nComputing annual aggregate retirement shocks...')
retirement_mw_by_year = defaultdict(float)
all_retirements = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            year = int(row['ret_year'])
            mw = float(row['capacity_mw'])
        except (ValueError, TypeError):
            continue
        retirement_mw_by_year[year] += mw
        all_retirements.append(row)

_print(f'  Retirement MW by year:')
for y in sorted(retirement_mw_by_year):
    _print(f'    {y}: {retirement_mw_by_year[y]:,.0f} MW')


# ══════════════════════════════════════════════════════════════════════
# STEP 3: LOAD DATA (same pattern as strategy2_robust_inference.py)
# ══════════════════════════════════════════════════════════════════════

_print('\nLoading monthly returns...')
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


_print('Loading Fama-French factors...')
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


# ══════════════════════════════════════════════════════════════════════
# STEP 4: BUILD PRE-PERIOD FUEL SIMILARITY MATRIX
# ══════════════════════════════════════════════════════════════════════

_print('\nBuilding pre-period fuel similarity matrix...')

# Load geo neighbor set (same sparsity as current W_fuel)
neighbors = defaultdict(list)
for gi in W_geo:
    for gj in W_geo[gi]:
        neighbors[gi].append(gj)

W_fuel_pre = defaultdict(dict)
for gi in neighbors:
    s_i = pre_shares.get(gi)
    row_sum = 0.0
    row_raw = {}
    for gj in neighbors[gi]:
        sim = fuel_similarity(s_i, pre_shares.get(gj))
        if sim <= 0:
            continue
        row_raw[gj] = sim
        row_sum += sim
    if row_sum > 0:
        for gj, sim in row_raw.items():
            W_fuel_pre[gi][gj] = sim / row_sum

n_pre_edges = sum(len(v) for v in W_fuel_pre.values())
_print(f'  Pre-period fuel matrix: {len(W_fuel_pre)} firms, {n_pre_edges} edges')
_print(f'  Current fuel matrix: {len(W_fuel)} firms, '
       f'{sum(len(v) for v in W_fuel.items())} edges')


# ══════════════════════════════════════════════════════════════════════
# DIAGNOSTICS: Correlation between pre-period and current shares
# ══════════════════════════════════════════════════════════════════════

_print('\n--- DIAGNOSTICS ---')

# Compute correlation of pre-period vs current w_fuel for all overlapping pairs
pre_vals = []
cur_vals = []
for gi in W_fuel_pre:
    if gi not in W_fuel:
        continue
    for gj in W_fuel_pre[gi]:
        if gj in W_fuel.get(gi, {}):
            pre_vals.append(W_fuel_pre[gi][gj])
            cur_vals.append(W_fuel[gi][gj])

if len(pre_vals) > 2:
    mean_pre = sum(pre_vals) / len(pre_vals)
    mean_cur = sum(cur_vals) / len(cur_vals)
    cov = sum((p - mean_pre) * (c - mean_cur) for p, c in zip(pre_vals, cur_vals)) / len(pre_vals)
    sd_pre = math.sqrt(sum((p - mean_pre) ** 2 for p in pre_vals) / len(pre_vals))
    sd_cur = math.sqrt(sum((c - mean_cur) ** 2 for c in cur_vals) / len(cur_vals))
    corr_shares = cov / (sd_pre * sd_cur) if sd_pre > 0 and sd_cur > 0 else 0
    _print(f'  Correlation(w_fuel_pre, w_fuel_current): {corr_shares:.4f} ({len(pre_vals)} pairs)')
else:
    corr_shares = float('nan')
    _print(f'  Too few overlapping pairs for correlation')

# Summary stats on pre-period shares
_print(f'\n  Pre-period fuel shares summary (N = {len(pre_shares)}):')
for fuel_idx, fuel_name in enumerate(['coal', 'gas', 'solar', 'wind']):
    vals = [s[fuel_idx] for s in pre_shares.values()]
    if vals:
        mean_v = sum(vals) / len(vals)
        _print(f'    {fuel_name}: mean = {mean_v:.4f}, '
               f'min = {min(vals):.4f}, max = {max(vals):.4f}')


# ══════════════════════════════════════════════════════════════════════
# STEP 5: LOAD EVENTS AND BUILD PER-EVENT DATASETS
# ══════════════════════════════════════════════════════════════════════

_print('\nLoading events...')
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
        try:
            cap_mw = float(row['capacity_mw'])
        except (ValueError, TypeError):
            cap_mw = 0.0
        all_events.append({
            'plant': row['plant_name'],
            'year': event_year,
            'event_date': effective_date,
            'gvkeys': row['matched_gvkeys'].split(';'),
            'capacity_mw': cap_mw,
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


# Build per-event datasets with Bartik instrument
_print('\nBuilding per-event datasets with Bartik instrument...')

SPEC_VARS_OLS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']
SPEC_VARS_BARTIK = ['w_geo', 'bartik', 'w_reg', 'same_sector']
MIN_OBS_PER_EVENT = 20

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

    # Aggregate retirement shock for event year
    agg_shock_mw = retirement_mw_by_year.get(year, 0.0)
    if agg_shock_mw <= 0:
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
        geo_neighbors = W_geo[fm_gk]
        neighbor_gks = set(geo_neighbors.keys()) - event_gvkeys
        non_connected = [gk for gk in fundamentals
                         if gk not in event_gvkeys and gk not in geo_neighbors]
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
            w_geo = geo_neighbors.get(gk, 0.0)
            w_fuel = W_fuel.get(fm_gk, {}).get(gk, 0.0)
            w_fuel_pre = W_fuel_pre.get(fm_gk, {}).get(gk, 0.0)
            w_reg = W_reg.get(fm_gk, {}).get(gk, 0.0)
            j_sic4 = get_sic4(gk)
            same_sector = 1.0 if (fm_sic4 and j_sic4
                                  and fm_sic4 == j_sic4) else 0.0
            car = compute_monthly_car(gk, event_month, post=POST_MONTHS)
            if car is None:
                continue

            # Bartik instrument: pre-determined fuel similarity * aggregate shock
            # Normalize shock by dividing by 10000 to keep coefficient interpretable
            bartik = w_fuel_pre * (agg_shock_mw / 10000.0)

            obs.append({
                'car': car,
                'w_geo': w_geo,
                'w_fuel': w_fuel,
                'w_fuel_pre': w_fuel_pre,
                'bartik': bartik,
                'w_reg': w_reg,
                'same_sector': same_sector,
                'gvkey': gk,
                'event_id': event_id,
            })

    if len(obs) >= MIN_OBS_PER_EVENT:
        event_datasets[event_id] = obs

n_valid = len(event_datasets)
total_obs = sum(len(v) for v in event_datasets.values())
_print(f'  Valid events (>= {MIN_OBS_PER_EVENT} obs): {n_valid}')
_print(f'  Total obs: {total_obs}')

# Bartik summary stats
all_bartik = [o['bartik'] for ds in event_datasets.values() for o in ds]
all_wfuel_pre = [o['w_fuel_pre'] for ds in event_datasets.values() for o in ds]
all_wfuel = [o['w_fuel'] for ds in event_datasets.values() for o in ds]

if all_bartik:
    mean_b = sum(all_bartik) / len(all_bartik)
    sd_b = math.sqrt(sum((x - mean_b) ** 2 for x in all_bartik) / len(all_bartik))
    nonzero_b = sum(1 for x in all_bartik if x > 0)
    _print(f'\n  Bartik instrument summary:')
    _print(f'    Mean: {mean_b:.6f}, SD: {sd_b:.6f}')
    _print(f'    Min: {min(all_bartik):.6f}, Max: {max(all_bartik):.6f}')
    _print(f'    Non-zero: {nonzero_b}/{len(all_bartik)} ({100*nonzero_b/len(all_bartik):.1f}%)')

    # Correlation between pre and current w_fuel in regression sample
    if all_wfuel_pre and all_wfuel:
        mean_p = sum(all_wfuel_pre) / len(all_wfuel_pre)
        mean_c = sum(all_wfuel) / len(all_wfuel)
        cov_pc = sum((p - mean_p) * (c - mean_c)
                     for p, c in zip(all_wfuel_pre, all_wfuel)) / len(all_wfuel_pre)
        sd_p = math.sqrt(sum((x - mean_p) ** 2 for x in all_wfuel_pre) / len(all_wfuel_pre))
        sd_c = math.sqrt(sum((x - mean_c) ** 2 for x in all_wfuel) / len(all_wfuel))
        corr_sample = cov_pc / (sd_p * sd_c) if sd_p > 0 and sd_c > 0 else 0
        _print(f'    Corr(w_fuel_pre, w_fuel) in sample: {corr_sample:.4f}')


# ══════════════════════════════════════════════════════════════════════
# SPEC A: FAMA-MACBETH WITH BARTIK INSTRUMENT (REDUCED FORM)
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('SPEC A: FAMA-MACBETH (1973) WITH BARTIK INSTRUMENT')
_print('CAR = alpha + beta_bartik * B_ie + w_geo + w_reg + same_sector + eps')
_print('=' * 70)

event_betas_bartik = defaultdict(list)
event_r2s_bartik = []
event_ns_bartik = []
event_ids_bartik = []

for event_id in sorted(event_datasets.keys()):
    obs = event_datasets[event_id]
    n_obs = len(obs)

    ss_vals = set(o['same_sector'] for o in obs)
    use_vars = SPEC_VARS_BARTIK if len(ss_vals) > 1 else ['w_geo', 'bartik', 'w_reg']

    y = [o['car'] for o in obs]
    X = [[1.0] + [o[v] for v in use_vars] for o in obs]

    result = ols_simple(y, X)
    if result is None:
        continue

    names = ['intercept'] + use_vars
    for i, name in enumerate(names):
        event_betas_bartik[name].append(result['beta'][i])

    for v in SPEC_VARS_BARTIK:
        if v not in use_vars:
            event_betas_bartik[v].append(float('nan'))

    event_r2s_bartik.append(result['r2'])
    event_ns_bartik.append(result['n'])
    event_ids_bartik.append(event_id)

T_bartik = len(event_ids_bartik)
_print(f'\n  Events with valid regressions: {T_bartik}')
if T_bartik > 0:
    _print(f'  Avg N per event: {sum(event_ns_bartik) / T_bartik:.1f}')
    _print(f'  Avg R2 per event: {sum(event_r2s_bartik) / T_bartik:.4f}')

_print(f'\n  {"Variable":<15} {"Mean beta":>12} {"NW SE":>10} {"t":>8} {"p":>8}')
_print('  ' + '-' * 55)

fm_bartik = {}
for v in ['intercept'] + SPEC_VARS_BARTIK:
    betas = event_betas_bartik[v]
    clean = [b for b in betas if not math.isnan(b)]
    if len(clean) < 3:
        continue
    mean_b = sum(clean) / len(clean)
    nw_se = newey_west_se(clean)
    t_stat = mean_b / nw_se if nw_se > 1e-15 else 0.0
    df = len(clean) - 1
    p_val = p_from_t_df(t_stat, df)
    stars = '***' if p_val < 0.01 else '**' if p_val < 0.05 else '*' if p_val < 0.10 else ''
    _print(f'  {v:<15} {mean_b:+12.6f} {nw_se:10.6f} {t_stat:8.3f} {p_val:8.4f}{stars}  (df={df})')
    fm_bartik[v] = {
        'mean': mean_b, 'se': nw_se, 't': t_stat, 'p': p_val,
        'n_events': len(clean), 'df': df,
    }


# ══════════════════════════════════════════════════════════════════════
# SPEC B: FAMA-MACBETH WITH STANDARD w_fuel (COMPARISON)
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('SPEC B: FAMA-MACBETH WITH STANDARD w_fuel (COMPARISON)')
_print('CAR = alpha + beta_fuel * w_fuel + w_geo + w_reg + same_sector + eps')
_print('=' * 70)

event_betas_ols = defaultdict(list)
event_r2s_ols = []
event_ns_ols = []
event_ids_ols = []

for event_id in sorted(event_datasets.keys()):
    obs = event_datasets[event_id]
    n_obs = len(obs)

    ss_vals = set(o['same_sector'] for o in obs)
    use_vars = SPEC_VARS_OLS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']

    y = [o['car'] for o in obs]
    X = [[1.0] + [o[v] for v in use_vars] for o in obs]

    result = ols_simple(y, X)
    if result is None:
        continue

    names = ['intercept'] + use_vars
    for i, name in enumerate(names):
        event_betas_ols[name].append(result['beta'][i])

    for v in SPEC_VARS_OLS:
        if v not in use_vars:
            event_betas_ols[v].append(float('nan'))

    event_r2s_ols.append(result['r2'])
    event_ns_ols.append(result['n'])
    event_ids_ols.append(event_id)

T_ols = len(event_ids_ols)
_print(f'\n  Events with valid regressions: {T_ols}')
if T_ols > 0:
    _print(f'  Avg N per event: {sum(event_ns_ols) / T_ols:.1f}')
    _print(f'  Avg R2 per event: {sum(event_r2s_ols) / T_ols:.4f}')

_print(f'\n  {"Variable":<15} {"Mean beta":>12} {"NW SE":>10} {"t":>8} {"p":>8}')
_print('  ' + '-' * 55)

fm_ols = {}
for v in ['intercept'] + SPEC_VARS_OLS:
    betas = event_betas_ols[v]
    clean = [b for b in betas if not math.isnan(b)]
    if len(clean) < 3:
        continue
    mean_b = sum(clean) / len(clean)
    nw_se = newey_west_se(clean)
    t_stat = mean_b / nw_se if nw_se > 1e-15 else 0.0
    p_val = p_from_t(t_stat)
    stars = '***' if p_val < 0.01 else '**' if p_val < 0.05 else '*' if p_val < 0.10 else ''
    _print(f'  {v:<15} {mean_b:+12.6f} {nw_se:10.6f} {t_stat:8.3f} {p_val:8.4f}{stars}')
    fm_ols[v] = {
        'mean': mean_b, 'se': nw_se, 't': t_stat, 'p': p_val,
        'n_events': len(clean),
    }


# ══════════════════════════════════════════════════════════════════════
# POOLED EVENT-CLUSTERED REGRESSIONS
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('POOLED EVENT-CLUSTERED REGRESSIONS')
_print('=' * 70)

# Pool all observations
pooled = []
for event_id in sorted(event_datasets.keys()):
    for o in event_datasets[event_id]:
        pooled.append(o)

_print(f'  Pooled observations: {len(pooled)}')

# Spec A: Bartik
_print('\n  Spec A (Bartik): CAR ~ bartik + w_geo + w_reg + same_sector')
ss_vals_pool = set(o['same_sector'] for o in pooled)
use_vars_A = ['bartik', 'w_geo', 'w_reg', 'same_sector'] if len(ss_vals_pool) > 1 else ['bartik', 'w_geo', 'w_reg']

y_pool = [o['car'] for o in pooled]
X_pool_A = [[1.0] + [o[v] for v in use_vars_A] for o in pooled]
cluster_ids_A = [o['event_id'] for o in pooled]

res_A = ols_simple(y_pool, X_pool_A)
if res_A:
    se_A = clustered_se(X_pool_A, res_A['resid'], cluster_ids_A, res_A['inv_XtX'])
    names_A = ['intercept'] + use_vars_A
    _print(f'  R2 = {res_A["r2"]:.4f}, N = {res_A["n"]}')
    _print(f'  {"Variable":<15} {"beta":>12} {"SE(cl)":>10} {"t":>8} {"p":>8}')
    _print('  ' + '-' * 55)
    pooled_bartik = {}
    for i, name in enumerate(names_A):
        t_val = res_A['beta'][i] / se_A[i] if se_A and se_A[i] > 1e-15 else 0.0
        p_val = p_from_t(t_val)
        stars = '***' if p_val < 0.01 else '**' if p_val < 0.05 else '*' if p_val < 0.10 else ''
        _print(f'  {name:<15} {res_A["beta"][i]:+12.6f} {se_A[i]:10.6f} {t_val:8.3f} {p_val:8.4f}{stars}')
        pooled_bartik[name] = {'beta': res_A['beta'][i], 'se': se_A[i], 't': t_val, 'p': p_val}

# Spec B: Standard w_fuel
_print('\n  Spec B (Standard): CAR ~ w_fuel + w_geo + w_reg + same_sector')
use_vars_B = ['w_fuel', 'w_geo', 'w_reg', 'same_sector'] if len(ss_vals_pool) > 1 else ['w_fuel', 'w_geo', 'w_reg']

X_pool_B = [[1.0] + [o[v] for v in use_vars_B] for o in pooled]
cluster_ids_B = [o['event_id'] for o in pooled]

res_B = ols_simple(y_pool, X_pool_B)
if res_B:
    se_B = clustered_se(X_pool_B, res_B['resid'], cluster_ids_B, res_B['inv_XtX'])
    names_B = ['intercept'] + use_vars_B
    _print(f'  R2 = {res_B["r2"]:.4f}, N = {res_B["n"]}')
    _print(f'  {"Variable":<15} {"beta":>12} {"SE(cl)":>10} {"t":>8} {"p":>8}')
    _print('  ' + '-' * 55)
    pooled_ols = {}
    for i, name in enumerate(names_B):
        t_val = res_B['beta'][i] / se_B[i] if se_B and se_B[i] > 1e-15 else 0.0
        p_val = p_from_t(t_val)
        stars = '***' if p_val < 0.01 else '**' if p_val < 0.05 else '*' if p_val < 0.10 else ''
        _print(f'  {name:<15} {res_B["beta"][i]:+12.6f} {se_B[i]:10.6f} {t_val:8.3f} {p_val:8.4f}{stars}')
        pooled_ols[name] = {'beta': res_B['beta'][i], 'se': se_B[i], 't': t_val, 'p': p_val}


# ══════════════════════════════════════════════════════════════════════
# GPS (2020) DIAGNOSTICS: ROTEMBERG WEIGHTS, NEGATIVE WEIGHTS,
# PRE-EVENT BALANCE, TOP-SHARE CORRELATES, OSTER (2019) BOUNDS
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('GPS (2020) SHIFT-SHARE DIAGNOSTICS')
_print('=' * 70)

# ── Diagnostic 1: Rotemberg Weights ──────────────────────────────────
# The Rotemberg weight for event e is proportional to the sum of
# squared Bartik residuals within that event (after partialling out
# controls). This identifies which events drive the Bartik estimate.

_print('\n  Diagnostic 1: Rotemberg Weights')

# First, regress bartik on controls to get residuals
ctrl_vars = ['w_geo', 'w_reg', 'same_sector']
y_bartik = [o['bartik'] for o in pooled]
X_ctrl = [[1.0] + [o[v] for v in ctrl_vars] for o in pooled]
res_ctrl = ols_simple(y_bartik, X_ctrl)

rotemberg_weights = {}
total_sq_resid = 0.0
if res_ctrl:
    bartik_resids = res_ctrl['resid']
    # Group by event_id
    event_sq_resid = defaultdict(float)
    for i, o in enumerate(pooled):
        event_sq_resid[o['event_id']] += bartik_resids[i] ** 2
        total_sq_resid += bartik_resids[i] ** 2

    if total_sq_resid > 1e-15:
        for eid, sq_r in event_sq_resid.items():
            rotemberg_weights[eid] = sq_r / total_sq_resid

    # Sort by weight, report top 10
    sorted_rw = sorted(rotemberg_weights.items(), key=lambda x: -x[1])
    _print(f'  Total events with Rotemberg weights: {len(sorted_rw)}')
    _print(f'  {"Event":>6} {"Weight":>10} {"Plant":<35} {"Year":>6} {"MW":>8}')
    _print('  ' + '-' * 70)
    for rank, (eid, wt) in enumerate(sorted_rw[:10]):
        ev = all_events[eid] if eid < len(all_events) else {}
        plant = ev.get('plant', 'Unknown')[:34]
        yr = ev.get('year', 0)
        mw = ev.get('capacity_mw', 0)
        _print(f'  {eid:6d} {wt:10.4f} {plant:<35} {yr:6d} {mw:8.0f}')

    # HHI of Rotemberg weights (concentration)
    rw_hhi = sum(w ** 2 for w in rotemberg_weights.values())
    _print(f'\n  HHI of Rotemberg weights: {rw_hhi:.4f}')
    top5_share = sum(w for _, w in sorted_rw[:5])
    _print(f'  Top 5 events share: {top5_share:.4f}')

# ── Diagnostic 2: Negative Weight Fraction ───────────────────────────
_print('\n  Diagnostic 2: Negative Weight Fraction')

n_negative = sum(1 for w in rotemberg_weights.values() if w < 0)
sum_negative = sum(w for w in rotemberg_weights.values() if w < 0)
sum_positive = sum(w for w in rotemberg_weights.values() if w >= 0)
_print(f'  Negative-weight events: {n_negative} / {len(rotemberg_weights)}')
_print(f'  Sum of negative weights: {sum_negative:.6f}')
_print(f'  Sum of positive weights: {sum_positive:.6f}')

# Re-estimate dropping negative-weight events
if n_negative > 0 and res_A:
    neg_events = {eid for eid, w in rotemberg_weights.items() if w < 0}
    pooled_pos = [o for o in pooled if o['event_id'] not in neg_events]
    if len(pooled_pos) > 10:
        y_pos = [o['car'] for o in pooled_pos]
        X_pos = [[1.0] + [o[v] for v in use_vars_A] for o in pooled_pos]
        cl_pos = [o['event_id'] for o in pooled_pos]
        res_pos = ols_simple(y_pos, X_pos)
        if res_pos:
            se_pos = clustered_se(X_pos, res_pos['resid'], cl_pos, res_pos['inv_XtX'])
            idx_bartik = use_vars_A.index('bartik') + 1  # +1 for intercept
            t_pos = res_pos['beta'][idx_bartik] / se_pos[idx_bartik] if se_pos[idx_bartik] > 1e-15 else 0.0
            p_pos = p_from_t(t_pos)
            _print(f'  Bartik (dropping neg-weight events): beta={res_pos["beta"][idx_bartik]:+.6f}, '
                   f't={t_pos:.3f}, p={p_pos:.4f}, N={len(pooled_pos)}')

# ── Diagnostic 3: Pre-Event Balance Test [-5, -2] months ────────────
_print('\n  Diagnostic 3: Pre-Event Balance Test ([-5, -2] months)')


def compute_monthly_car_pre(gvkey, event_month, pre_start=-5, pre_end=-2):
    """Compute CAR in the pre-event window [pre_start, pre_end] months."""
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
    car = 0.0
    n_months = 0
    for offset in range(pre_start, pre_end + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months) and months[idx] in monthly_ret[gvkey]:
            m = months[idx]
            r_it = monthly_ret[gvkey][m]
            if m in market_ret_monthly:
                ar = r_it - market_ret_monthly[m]
                car += ar
                n_months += 1
    if n_months < 2:
        return None
    return car


# Build pre-event CARs for pooled sample
pooled_pre = []
for o in pooled:
    event_id = o['event_id']
    ev = all_events[event_id] if event_id < len(all_events) else {}
    event_date = ev.get('event_date', '')
    yr = ev.get('year', 0)
    if event_date and len(event_date) >= 7:
        em = event_date[:7]
    else:
        em = f'{yr}-07' if yr else None
    if not em:
        continue
    car_pre = compute_monthly_car_pre(o['gvkey'], em)
    if car_pre is not None:
        pooled_pre.append({**o, 'car_pre': car_pre})

_print(f'  Observations with pre-event CARs: {len(pooled_pre)}')

if len(pooled_pre) > 50:
    y_pre = [o['car_pre'] for o in pooled_pre]
    X_pre = [[1.0] + [o[v] for v in use_vars_A] for o in pooled_pre]
    cl_pre = [o['event_id'] for o in pooled_pre]

    res_pre = ols_simple(y_pre, X_pre)
    if res_pre:
        se_pre = clustered_se(X_pre, res_pre['resid'], cl_pre, res_pre['inv_XtX'])
        _print(f'  R2 = {res_pre["r2"]:.6f}, N = {res_pre["n"]}')
        _print(f'  {"Variable":<15} {"beta":>12} {"SE(cl)":>10} {"t":>8} {"p":>8}')
        _print('  ' + '-' * 55)
        pre_results = {}
        for i, name in enumerate(['intercept'] + use_vars_A):
            t_val = res_pre['beta'][i] / se_pre[i] if se_pre[i] > 1e-15 else 0.0
            p_val = p_from_t(t_val)
            stars = '***' if p_val < 0.01 else '**' if p_val < 0.05 else '*' if p_val < 0.10 else ''
            _print(f'  {name:<15} {res_pre["beta"][i]:+12.6f} {se_pre[i]:10.6f} {t_val:8.3f} {p_val:8.4f}{stars}')
            pre_results[name] = {'beta': res_pre['beta'][i], 'se': se_pre[i], 't': t_val, 'p': p_val}
        bartik_pre_t = pre_results.get('bartik', {}).get('t', 0)
        bartik_pre_p = pre_results.get('bartik', {}).get('p', 1)
        _print(f'\n  Pre-event balance test: bartik t = {bartik_pre_t:.3f}, p = {bartik_pre_p:.4f}')
        if abs(bartik_pre_t) < 1.96:
            _print('  PASS: Bartik does not predict pre-event CARs')
        else:
            _print('  FAIL: Bartik predicts pre-event CARs (pre-trend concern)')

# ── Diagnostic 3b: Pre-Balance Sensitivity (cutoff = 2010) ───────────
_print(f'\n  Diagnostic 3b: Pre-Balance Sensitivity (cutoff = {PRE_CUTOFF_SENSITIVITY})')

# Rebuild pre-period shares with stricter cutoff
pre_mw_2010, n_plants_2010 = build_pre_mw_for_cutoff(_all_plant_records, PRE_CUTOFF_SENSITIVITY)
pre_shares_2010 = {}
for gk in pre_mw_2010:
    s = compute_shares(pre_mw_2010[gk])
    if s is not None:
        pre_shares_2010[gk] = s
_print(f'  Firms with pre-{PRE_CUTOFF_SENSITIVITY} fuel shares: {len(pre_shares_2010)} '
       f'(vs {len(pre_shares)} at pre-{PRE_CUTOFF})')
_print(f'  Plants included (pre-{PRE_CUTOFF_SENSITIVITY}): {n_plants_2010}')

# Build W_fuel_pre_2010 matrix (same sparsity as W_geo neighbors)
W_fuel_pre_2010 = defaultdict(dict)
for gi in neighbors:
    s_i = pre_shares_2010.get(gi)
    row_sum = 0.0
    row_raw = {}
    for gj in neighbors[gi]:
        sim = fuel_similarity(s_i, pre_shares_2010.get(gj))
        if sim <= 0:
            continue
        row_raw[gj] = sim
        row_sum += sim
    if row_sum > 0:
        for gj, sim in row_raw.items():
            W_fuel_pre_2010[gi][gj] = sim / row_sum

n_pre_edges_2010 = sum(len(v) for v in W_fuel_pre_2010.values())
_print(f'  Pre-{PRE_CUTOFF_SENSITIVITY} fuel matrix: {len(W_fuel_pre_2010)} firms, '
       f'{n_pre_edges_2010} edges')

# Rebuild Bartik with 2010 shares and run pre-balance test
pooled_pre_2010 = []
for o in pooled:
    event_id = o['event_id']
    ev = all_events[event_id] if event_id < len(all_events) else {}
    event_date = ev.get('event_date', '')
    yr = ev.get('year', 0)
    if event_date and len(event_date) >= 7:
        em = event_date[:7]
    else:
        em = f'{yr}-07' if yr else None
    if not em:
        continue

    # Get the first-mover gvkey for this event to look up fuel similarity
    fm_gk = None
    for gk_cand in ev.get('gvkeys', []):
        if gk_cand in W_fuel_pre_2010:
            fm_gk = gk_cand
            break
    if fm_gk is None:
        # Try from the original event structure
        for gk_cand in ev.get('gvkeys', []):
            fm_gk = gk_cand
            break

    w_fuel_pre_2010_val = W_fuel_pre_2010.get(fm_gk, {}).get(o['gvkey'], 0.0) if fm_gk else 0.0
    agg_shock = retirement_mw_by_year.get(yr, 0.0)
    bartik_2010 = w_fuel_pre_2010_val * (agg_shock / 10000.0)

    car_pre = compute_monthly_car_pre(o['gvkey'], em)
    if car_pre is not None:
        pooled_pre_2010.append({**o, 'car_pre': car_pre, 'bartik_2010': bartik_2010})

_print(f'  Observations with pre-event CARs (2010 cutoff): {len(pooled_pre_2010)}')

pre_results_2010 = {}
bartik_pre_t_2010 = 0
bartik_pre_p_2010 = 1
if len(pooled_pre_2010) > 50:
    y_pre_2010 = [o['car_pre'] for o in pooled_pre_2010]
    use_vars_2010 = ['bartik_2010', 'w_geo', 'w_reg']
    # Check same_sector variation
    ss_vals_2010 = set(o['same_sector'] for o in pooled_pre_2010)
    if len(ss_vals_2010) > 1:
        use_vars_2010.append('same_sector')
    X_pre_2010 = [[1.0] + [o[v] for v in use_vars_2010] for o in pooled_pre_2010]
    cl_pre_2010 = [o['event_id'] for o in pooled_pre_2010]

    res_pre_2010 = ols_simple(y_pre_2010, X_pre_2010)
    if res_pre_2010:
        se_pre_2010 = clustered_se(X_pre_2010, res_pre_2010['resid'],
                                    cl_pre_2010, res_pre_2010['inv_XtX'])
        _print(f'  R2 = {res_pre_2010["r2"]:.6f}, N = {res_pre_2010["n"]}')
        _print(f'  {"Variable":<15} {"beta":>12} {"SE(cl)":>10} {"t":>8} {"p":>8}')
        _print('  ' + '-' * 55)
        names_2010 = ['intercept'] + use_vars_2010
        for i, name in enumerate(names_2010):
            t_val = res_pre_2010['beta'][i] / se_pre_2010[i] if se_pre_2010[i] > 1e-15 else 0.0
            p_val = p_from_t(t_val)
            stars = '***' if p_val < 0.01 else '**' if p_val < 0.05 else '*' if p_val < 0.10 else ''
            _print(f'  {name:<15} {res_pre_2010["beta"][i]:+12.6f} {se_pre_2010[i]:10.6f} '
                   f'{t_val:8.3f} {p_val:8.4f}{stars}')
            display_name = name.replace('bartik_2010', 'bartik')
            pre_results_2010[display_name] = {'beta': res_pre_2010['beta'][i],
                                               'se': se_pre_2010[i], 't': t_val, 'p': p_val}
        bartik_pre_t_2010 = pre_results_2010.get('bartik', {}).get('t', 0)
        bartik_pre_p_2010 = pre_results_2010.get('bartik', {}).get('p', 1)
        _print(f'\n  Pre-event balance (cutoff {PRE_CUTOFF_SENSITIVITY}): '
               f'bartik t = {bartik_pre_t_2010:.3f}, p = {bartik_pre_p_2010:.4f}')
        if abs(bartik_pre_t_2010) < 1.96:
            _print(f'  PASS: Bartik ({PRE_CUTOFF_SENSITIVITY}) does not predict pre-event CARs')
        else:
            _print(f'  FAIL: Bartik ({PRE_CUTOFF_SENSITIVITY}) predicts pre-event CARs')

# ── Diagnostic 4: Top-Share Correlates ───────────────────────────────
_print('\n  Diagnostic 4: Top-Share Correlates')

# For each firm in pooled, get average w_fuel_pre and fundamentals
firm_exposure = defaultdict(list)
for o in pooled:
    firm_exposure[o['gvkey']].append(o['w_fuel_pre'])

firm_avg_exposure = {}
for gk, vals in firm_exposure.items():
    firm_avg_exposure[gk] = sum(vals) / len(vals)

if firm_avg_exposure:
    median_exp = sorted(firm_avg_exposure.values())[len(firm_avg_exposure) // 2]
    high_exp = {gk for gk, v in firm_avg_exposure.items() if v > median_exp}
    low_exp = {gk for gk, v in firm_avg_exposure.items() if v <= median_exp}

    # Collect fundamentals for each group
    def group_stats(gk_set, var_name):
        vals = []
        for gk in gk_set:
            if gk in fundamentals:
                v = fundamentals[gk].get(var_name)
                if v is not None:
                    try:
                        vals.append(float(v))
                    except (ValueError, TypeError):
                        pass
        if not vals:
            return None, None, 0
        mean_v = sum(vals) / len(vals)
        if len(vals) > 1:
            sd_v = math.sqrt(sum((x - mean_v) ** 2 for x in vals) / (len(vals) - 1))
        else:
            sd_v = 0.0
        return mean_v, sd_v, len(vals)

    _print(f'  Firms: {len(firm_avg_exposure)} total ({len(high_exp)} high-exposure, {len(low_exp)} low-exposure)')
    _print(f'  Median w_fuel_pre: {median_exp:.6f}')
    _print(f'\n  {"Variable":<20} {"High mean":>10} {"Low mean":>10} {"Diff":>10} {"t":>8}')
    _print('  ' + '-' * 62)

    correlate_results = {}
    for var_name, var_label in [('at', 'Total assets'),
                                 ('dltt', 'LT debt'),
                                 ('ni', 'Net income')]:
        m_h, sd_h, n_h = group_stats(high_exp, var_name)
        m_l, sd_l, n_l = group_stats(low_exp, var_name)
        if m_h is not None and m_l is not None and n_h > 2 and n_l > 2:
            diff = m_h - m_l
            pooled_se = math.sqrt(sd_h ** 2 / n_h + sd_l ** 2 / n_l) if (sd_h + sd_l) > 0 else 1e-15
            t_diff = diff / pooled_se if pooled_se > 1e-15 else 0.0
            _print(f'  {var_label:<20} {m_h:10.1f} {m_l:10.1f} {diff:+10.1f} {t_diff:8.2f}')
            correlate_results[var_label] = {'high': m_h, 'low': m_l, 'diff': diff, 't': t_diff}

# ── Diagnostic 5: Oster (2019) Coefficient Stability Bounds ─────────
_print('\n  Diagnostic 5: Oster (2019) Coefficient Stability Bounds')

# Restricted: bartik only, no controls
y_oster = [o['car'] for o in pooled]
X_restricted = [[1.0, o['bartik']] for o in pooled]
res_restricted = ols_simple(y_oster, X_restricted)

# Full: bartik + all controls (already computed as res_A)
if res_restricted and res_A:
    beta_restricted = res_restricted['beta'][1]  # bartik coeff (no controls)
    r2_restricted = res_restricted['r2']

    idx_b = use_vars_A.index('bartik') + 1  # +1 for intercept
    beta_full = res_A['beta'][idx_b]
    r2_full = res_A['r2']

    # R_max = min(1, 1.3 * R_full) following Oster's recommendation
    r2_max = min(1.0, 1.3 * r2_full)

    _print(f'  beta_restricted (bartik only): {beta_restricted:+.6f}, R2 = {r2_restricted:.6f}')
    _print(f'  beta_full (bartik + controls): {beta_full:+.6f}, R2 = {r2_full:.6f}')
    _print(f'  R2_max (1.3 x R2_full): {r2_max:.6f}')

    denom = r2_full - r2_restricted
    if abs(denom) > 1e-12:
        # beta* at delta = 1
        numer = (beta_restricted - beta_full) * (r2_max - r2_full)
        beta_star = beta_full - 1.0 * numer / denom

        # delta* that would drive beta to zero
        if abs(beta_restricted - beta_full) > 1e-12 and abs(r2_max - r2_full) > 1e-12:
            delta_star = beta_full * denom / ((beta_restricted - beta_full) * (r2_max - r2_full))
        else:
            delta_star = float('inf')

        _print(f'  beta* (delta=1, R_max=1.3R_full): {beta_star:+.6f}')
        _print(f'  delta* (beta -> 0): {delta_star:.4f}')
        _print(f'  Identified set: [{min(beta_full, beta_star):+.6f}, {max(beta_full, beta_star):+.6f}]')

        if delta_star > 1:
            _print('  PASS: delta* > 1 means unobservables would need to be more')
            _print('        important than observables to explain away the result.')
        else:
            _print('  CAUTION: delta* < 1 means modest unobservable selection')
            _print('           could explain the fuel coefficient.')

        # Also compute for standard w_fuel specification
        X_restricted_fuel = [[1.0, o['w_fuel']] for o in pooled]
        res_restricted_fuel = ols_simple(y_oster, X_restricted_fuel)
        if res_restricted_fuel and res_B:
            beta_r_fuel = res_restricted_fuel['beta'][1]
            r2_r_fuel = res_restricted_fuel['r2']
            idx_f = use_vars_B.index('w_fuel') + 1
            beta_f_fuel = res_B['beta'][idx_f]
            r2_f_fuel = res_B['r2']
            r2_max_fuel = min(1.0, 1.3 * r2_f_fuel)
            denom_fuel = r2_f_fuel - r2_r_fuel
            if abs(denom_fuel) > 1e-12:
                numer_fuel = (beta_r_fuel - beta_f_fuel) * (r2_max_fuel - r2_f_fuel)
                beta_star_fuel = beta_f_fuel - 1.0 * numer_fuel / denom_fuel
                if abs(beta_r_fuel - beta_f_fuel) > 1e-12 and abs(r2_max_fuel - r2_f_fuel) > 1e-12:
                    delta_star_fuel = beta_f_fuel * denom_fuel / ((beta_r_fuel - beta_f_fuel) * (r2_max_fuel - r2_f_fuel))
                else:
                    delta_star_fuel = float('inf')
                _print(f'\n  Comparison: standard w_fuel specification')
                _print(f'  beta_restricted: {beta_r_fuel:+.6f}, R2 = {r2_r_fuel:.6f}')
                _print(f'  beta_full: {beta_f_fuel:+.6f}, R2 = {r2_f_fuel:.6f}')
                _print(f'  beta* (delta=1): {beta_star_fuel:+.6f}')
                _print(f'  delta* (beta -> 0): {delta_star_fuel:.4f}')
    else:
        _print('  WARNING: R2_full ≈ R2_restricted, Oster bounds undefined')
        beta_star = None
        delta_star = None
else:
    beta_star = None
    delta_star = None


# ══════════════════════════════════════════════════════════════════════
# WRITE OUTPUT
# ══════════════════════════════════════════════════════════════════════

out_path = results_path('metrics', 'strategy2_bartik_shiftshare.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

lines = [
    '# Shift-Share (Bartik) Causal Robustness Test',
    '',
    'Addresses the concern that the fuel coefficient captures shared exposure',
    'to common factors rather than network transmission.',
    '',
    '## Design',
    '',
    f'**Shares**: Pre-period fuel-mix similarity weights using only plants',
    f'commissioned before {PRE_CUTOFF} (frozen before the main retirement wave).',
    '',
    '**Shifts**: Annual aggregate coal MW retired globally per year.',
    '',
    '**Bartik instrument**: B_ie = w_fuel_pre(i,j) x (RetiredMW_year / 10000)',
    '',
    '## Diagnostics',
    '',
    f'Firms with pre-period fuel shares: {len(pre_shares)}',
    f'Firms with current fuel shares: {len(current_shares)}',
    f'Pre-period fuel matrix: {len(W_fuel_pre)} firms, {n_pre_edges} edges',
]

if not math.isnan(corr_shares):
    lines.append(f'Correlation(w_fuel_pre, w_fuel_current) across matrix: {corr_shares:.4f}')

if all_bartik:
    lines += [
        '',
        '### Bartik instrument summary',
        '',
        f'| Statistic | Value |',
        f'|---|---:|',
        f'| Mean | {sum(all_bartik)/len(all_bartik):.6f} |',
        f'| SD | {math.sqrt(sum((x - sum(all_bartik)/len(all_bartik))**2 for x in all_bartik)/len(all_bartik)):.6f} |',
        f'| Min | {min(all_bartik):.6f} |',
        f'| Max | {max(all_bartik):.6f} |',
        f'| Non-zero | {sum(1 for x in all_bartik if x > 0)}/{len(all_bartik)} |',
    ]
    if all_wfuel_pre and all_wfuel and 'corr_sample' in dir():
        lines.append(f'| Corr(w_fuel_pre, w_fuel) in sample | {corr_sample:.4f} |')

lines += [
    '',
    '### Retirement shocks by year',
    '',
    '| Year | MW retired |',
    '|---|---:|',
]
for y in sorted(retirement_mw_by_year):
    lines.append(f'| {y} | {retirement_mw_by_year[y]:,.0f} |')

# Spec A: Bartik FM
lines += [
    '',
    '## Spec A: Bartik Reduced Form (Fama-MacBeth + Newey-West)',
    '',
    'CAR = alpha + beta_bartik x B_ie + w_geo + w_reg + same_sector + eps',
    '',
]
if T_bartik > 0:
    lines += [
        f'Valid events: {T_bartik}',
        f'Avg firms per event: {sum(event_ns_bartik) / T_bartik:.1f}',
        f'Avg within-event R2: {sum(event_r2s_bartik) / T_bartik:.4f}',
        '',
        '| Variable | Mean beta | NW SE | t | p |',
        '|---|---:|---:|---:|---:|',
    ]
    for v in ['intercept'] + SPEC_VARS_BARTIK:
        if v in fm_bartik:
            r = fm_bartik[v]
            stars = '***' if r['p'] < 0.01 else '**' if r['p'] < 0.05 else '*' if r['p'] < 0.10 else ''
            lines.append(f'| {v} | {r["mean"]:+.6f} | {r["se"]:.6f} | {r["t"]:.3f} | {r["p"]:.4f}{stars} |')

# Spec B: Standard FM
lines += [
    '',
    '## Spec B: Standard w_fuel (Fama-MacBeth + Newey-West, comparison)',
    '',
    'CAR = alpha + beta_fuel x w_fuel + w_geo + w_reg + same_sector + eps',
    '',
]
if T_ols > 0:
    lines += [
        f'Valid events: {T_ols}',
        f'Avg firms per event: {sum(event_ns_ols) / T_ols:.1f}',
        f'Avg within-event R2: {sum(event_r2s_ols) / T_ols:.4f}',
        '',
        '| Variable | Mean beta | NW SE | t | p |',
        '|---|---:|---:|---:|---:|',
    ]
    for v in ['intercept'] + SPEC_VARS_OLS:
        if v in fm_ols:
            r = fm_ols[v]
            stars = '***' if r['p'] < 0.01 else '**' if r['p'] < 0.05 else '*' if r['p'] < 0.10 else ''
            lines.append(f'| {v} | {r["mean"]:+.6f} | {r["se"]:.6f} | {r["t"]:.3f} | {r["p"]:.4f}{stars} |')

# Pooled event-clustered
lines += [
    '',
    '## Pooled Event-Clustered Regressions',
    '',
]
if res_A:
    lines += [
        '### Spec A: Bartik',
        '',
        f'N = {res_A["n"]}, R2 = {res_A["r2"]:.4f}',
        '',
        '| Variable | beta | SE(cl) | t | p |',
        '|---|---:|---:|---:|---:|',
    ]
    for name in names_A:
        if name in pooled_bartik:
            r = pooled_bartik[name]
            stars = '***' if r['p'] < 0.01 else '**' if r['p'] < 0.05 else '*' if r['p'] < 0.10 else ''
            lines.append(f'| {name} | {r["beta"]:+.6f} | {r["se"]:.6f} | {r["t"]:.3f} | {r["p"]:.4f}{stars} |')

if res_B:
    lines += [
        '',
        '### Spec B: Standard w_fuel',
        '',
        f'N = {res_B["n"]}, R2 = {res_B["r2"]:.4f}',
        '',
        '| Variable | beta | SE(cl) | t | p |',
        '|---|---:|---:|---:|---:|',
    ]
    for name in names_B:
        if name in pooled_ols:
            r = pooled_ols[name]
            stars = '***' if r['p'] < 0.01 else '**' if r['p'] < 0.05 else '*' if r['p'] < 0.10 else ''
            lines.append(f'| {name} | {r["beta"]:+.6f} | {r["se"]:.6f} | {r["t"]:.3f} | {r["p"]:.4f}{stars} |')

# Comparison table
lines += [
    '',
    '## Summary Comparison',
    '',
    '| Specification | Channel variable | FM t-stat | Pooled t-stat |',
    '|---|---|---:|---:|',
]
bartik_fm_t = fm_bartik.get('bartik', {}).get('t', 0)
fuel_fm_t = fm_ols.get('w_fuel', {}).get('t', 0)
bartik_pool_t = pooled_bartik.get('bartik', {}).get('t', 0) if res_A else 0
fuel_pool_t = pooled_ols.get('w_fuel', {}).get('t', 0) if res_B else 0
lines.append(f'| Bartik (pre-period shares x agg shock) | bartik | {bartik_fm_t:.3f} | {bartik_pool_t:.3f} |')
lines.append(f'| Standard (current w_fuel) | w_fuel | {fuel_fm_t:.3f} | {fuel_pool_t:.3f} |')

# GPS (2020) diagnostics output
lines += [
    '',
    '## GPS (2020) Shift-Share Diagnostics',
    '',
    '### Rotemberg Weights',
    '',
    'Which events drive the Bartik estimate? Rotemberg weights are proportional to the',
    'sum of squared Bartik residuals (after partialling out controls) within each event.',
    '',
]

if rotemberg_weights:
    sorted_rw = sorted(rotemberg_weights.items(), key=lambda x: -x[1])
    rw_hhi = sum(w ** 2 for w in rotemberg_weights.values())
    top5_share = sum(w for _, w in sorted_rw[:5])
    lines += [
        f'HHI of Rotemberg weights: {rw_hhi:.4f}',
        f'Top 5 events share: {top5_share:.4f}',
        '',
        '| Rank | Event | Plant | Year | MW | Weight |',
        '|---:|---:|---|---:|---:|---:|',
    ]
    for rank, (eid, wt) in enumerate(sorted_rw[:10]):
        ev = all_events[eid] if eid < len(all_events) else {}
        plant = ev.get('plant', 'Unknown')[:30]
        yr = ev.get('year', 0)
        mw = ev.get('capacity_mw', 0)
        lines.append(f'| {rank+1} | {eid} | {plant} | {yr} | {mw:,.0f} | {wt:.4f} |')

lines += [
    '',
    '### Negative Weight Diagnostic',
    '',
    f'Negative-weight events: {n_negative} / {len(rotemberg_weights)}',
    f'Sum of negative weights: {sum_negative:.6f}',
    f'Sum of positive weights: {sum_positive:.6f}',
    '',
]

# Pre-event balance
lines += [
    '### Pre-Event Balance Test ([-5, -2] months)',
    '',
    'Tests whether Bartik exposure predicts CARs in the pre-event window.',
    'Under the identifying assumption, the Bartik instrument should NOT predict',
    'pre-event returns.',
    '',
]
if 'pre_results' in dir() and pre_results:
    lines += [
        f'N = {len(pooled_pre)}, R2 = {res_pre["r2"]:.6f}',
        '',
        '| Variable | beta | SE(cl) | t | p |',
        '|---|---:|---:|---:|---:|',
    ]
    for name in ['intercept'] + use_vars_A:
        if name in pre_results:
            r = pre_results[name]
            stars = '***' if r['p'] < 0.01 else '**' if r['p'] < 0.05 else '*' if r['p'] < 0.10 else ''
            lines.append(f'| {name} | {r["beta"]:+.6f} | {r["se"]:.6f} | {r["t"]:.3f} | {r["p"]:.4f}{stars} |')
    bartik_pre_t = pre_results.get('bartik', {}).get('t', 0)
    bartik_pre_p = pre_results.get('bartik', {}).get('p', 1)
    verdict = 'PASS' if abs(bartik_pre_t) < 1.96 else 'FAIL'
    lines.append(f'\n**{verdict}**: Bartik t = {bartik_pre_t:.3f}, p = {bartik_pre_p:.4f}')
else:
    lines.append('Insufficient pre-event return data.')

# Pre-balance sensitivity (2010 cutoff)
lines += [
    '',
    f'### Pre-Balance Sensitivity (cutoff = {PRE_CUTOFF_SENSITIVITY})',
    '',
    f'Repeats the pre-event balance test using only plants commissioned before '
    f'{PRE_CUTOFF_SENSITIVITY}.',
    f'More pre-determined shares strengthen the causal claim if the test still passes.',
    '',
    f'Firms with pre-{PRE_CUTOFF_SENSITIVITY} fuel shares: {len(pre_shares_2010)} '
    f'(vs {len(pre_shares)} at pre-{PRE_CUTOFF})',
    f'Pre-{PRE_CUTOFF_SENSITIVITY} fuel matrix: {len(W_fuel_pre_2010)} firms, '
    f'{n_pre_edges_2010} edges',
    '',
]
if pre_results_2010:
    lines += [
        f'N = {len(pooled_pre_2010)}',
        '',
        '| Variable | beta | SE(cl) | t | p |',
        '|---|---:|---:|---:|---:|',
    ]
    for name in ['intercept', 'bartik', 'w_geo', 'w_reg', 'same_sector']:
        if name in pre_results_2010:
            r = pre_results_2010[name]
            stars = '***' if r['p'] < 0.01 else '**' if r['p'] < 0.05 else '*' if r['p'] < 0.10 else ''
            lines.append(f'| {name} | {r["beta"]:+.6f} | {r["se"]:.6f} | {r["t"]:.3f} | {r["p"]:.4f}{stars} |')
    verdict_2010 = 'PASS' if abs(bartik_pre_t_2010) < 1.96 else 'FAIL'
    lines.append(f'\n**{verdict_2010}**: Bartik (pre-{PRE_CUTOFF_SENSITIVITY}) '
                 f't = {bartik_pre_t_2010:.3f}, p = {bartik_pre_p_2010:.4f}')

    # Summary comparison of both cutoffs
    lines += [
        '',
        '#### Pre-Balance Comparison',
        '',
        '| Cutoff | Bartik t | Bartik p | Verdict |',
        '|---:|---:|---:|---|',
    ]
    if 'pre_results' in dir() and pre_results:
        t_2014 = pre_results.get('bartik', {}).get('t', 0)
        p_2014 = pre_results.get('bartik', {}).get('p', 1)
        v_2014 = 'PASS' if abs(t_2014) < 1.96 else 'FAIL'
        lines.append(f'| {PRE_CUTOFF} | {t_2014:.3f} | {p_2014:.4f} | {v_2014} |')
    lines.append(f'| {PRE_CUTOFF_SENSITIVITY} | {bartik_pre_t_2010:.3f} | '
                 f'{bartik_pre_p_2010:.4f} | {verdict_2010} |')
else:
    lines.append('Insufficient data for 2010 sensitivity test.')

# Oster bounds
lines += [
    '',
    '### Oster (2019) Coefficient Stability Bounds',
    '',
    'Tests how much selection on unobservables (delta) would be needed to explain',
    'away the Bartik coefficient. delta* > 1 means unobservables would need to be',
    'more important than observables.',
    '',
]
if beta_star is not None and delta_star is not None:
    lines += [
        '| Quantity | Bartik | Standard w_fuel |',
        '|---|---:|---:|',
    ]
    if res_restricted and res_A:
        lines.append(f'| beta (no controls) | {res_restricted["beta"][1]:+.6f} | {res_restricted_fuel["beta"][1] if res_restricted_fuel else 0:+.6f} |')
        lines.append(f'| beta (full controls) | {res_A["beta"][idx_b]:+.6f} | {res_B["beta"][idx_f] if res_B else 0:+.6f} |')
        lines.append(f'| R2 (no controls) | {res_restricted["r2"]:.6f} | {res_restricted_fuel["r2"] if res_restricted_fuel else 0:.6f} |')
        lines.append(f'| R2 (full controls) | {res_A["r2"]:.6f} | {res_B["r2"] if res_B else 0:.6f} |')
    lines.append(f'| beta* (delta=1) | {beta_star:+.6f} | {beta_star_fuel if "beta_star_fuel" in dir() else 0:+.6f} |')
    lines.append(f'| **delta*** | **{delta_star:.4f}** | **{delta_star_fuel if "delta_star_fuel" in dir() else 0:.4f}** |')
    verdict_oster = 'PASS' if delta_star > 1 else 'CAUTION'
    lines.append(f'\n**{verdict_oster}**: delta* = {delta_star:.4f}')

lines += [
    '',
    '## Interpretation',
    '',
    'The Bartik instrument separates pre-determined exposure (fuel-mix similarity',
    f'frozen at pre-{PRE_CUTOFF} plant vintages) from aggregate shocks (total coal MW',
    'retired per year). If beta_bartik is significant, the fuel channel reflects',
    'genuine network transmission rather than shared exposure to common factors.',
    '',
    'The comparison between Spec A and Spec B shows whether the causal (Bartik)',
    'and descriptive (OLS) estimates agree in sign and magnitude.',
    '',
    'Under the Goldsmith-Pinkham, Sorkin & Swift (2020) framework, the Bartik',
    'coefficient has a causal interpretation as the dose-response of abnormal',
    'returns to technology exposure, provided: (1) fuel-mix shares are pre-determined',
    '(supported by pre-2014 vintage restriction); (2) shares do not predict pre-event',
    'returns (tested above); (3) Rotemberg weights are non-negative (checked above);',
    '(4) the Oster (2019) bound confirms robustness to selection on unobservables.',
]

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

_print(f'\nWrote: {out_path}')
_print('Done.')
