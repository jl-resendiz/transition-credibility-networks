"""Geographic diversification test: does multinational scope attenuate w_geo?

The geographic channel (w_geo) is insignificant in the full sample. This script
tests the hypothesis that geographic competitive effects are real at the plant
level but get diversified away for multinationals operating in many countries.
Single-country utilities should show a stronger geographic effect.

Three specifications:
  1. Full sample baseline (fuel + geo + reg + same_sector)
  2. Single-country subsample (firms with all plants in one country)
  3. Diversification interaction (w_geo x n_countries, w_geo x log(n_countries))

Each specification reports:
  - Pooled OLS with event-clustered SEs
  - Fama-MacBeth + Newey-West

Output: results/metrics/strategy2_geo_diversification.md
"""
import csv
import os
import sys
import math
import random
import hashlib
import re
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
    return {'beta': beta, 'resid': resid, 'r2': r2, 'n': n,
            'inv_XtX': inv_XtX, 'X_mat': X_mat, 'y': y}


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


def sig_stars(p):
    if p < 0.01:
        return '***'
    if p < 0.05:
        return '**'
    if p < 0.10:
        return '*'
    return ''


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


# ── Event-clustered SEs ─────────────────────────────────────────────

def event_clustered_se(X_mat, resid, inv_XtX, cluster_ids):
    """Compute cluster-robust SEs (clustered by event)."""
    n = len(resid)
    k = len(X_mat[0])
    clusters = defaultdict(list)
    for i in range(n):
        clusters[cluster_ids[i]].append(i)
    G = len(clusters)
    if G < 2:
        return [float('inf')] * k

    # Meat: sum of (X_g' e_g)(X_g' e_g)' over clusters
    meat = [[0.0] * k for _ in range(k)]
    for cid, indices in clusters.items():
        # X_g' e_g
        xe = [sum(X_mat[i][a] * resid[i] for i in indices) for a in range(k)]
        for a in range(k):
            for b in range(k):
                meat[a][b] += xe[a] * xe[b]

    # Scale factor: G/(G-1) * (n-1)/(n-k)
    scale = (G / (G - 1.0)) * ((n - 1.0) / (n - k))

    # Sandwich: inv_XtX * meat * inv_XtX
    # First: meat * inv_XtX
    temp = [[sum(meat[a][c] * inv_XtX[c][b] for c in range(k))
             for b in range(k)] for a in range(k)]
    # Then: inv_XtX * temp
    V = [[sum(inv_XtX[a][c] * temp[c][b] for c in range(k))
          for b in range(k)] for a in range(k)]

    ses = []
    for a in range(k):
        v = V[a][a] * scale
        ses.append(math.sqrt(v) if v > 0 else float('inf'))
    return ses


# ══════════════════════════════════════════════════════════════════════
# STEP 1: Count plant-countries per gvkey from GEM trackers
# ══════════════════════════════════════════════════════════════════════

_print('Counting plant countries per firm from GEM trackers...')


def parse_parents(field):
    if not field or str(field).strip() == '':
        return []
    parts = str(field).split(';')
    results = []
    for p in parts:
        p = p.strip()
        match = re.match(r'^(.+?)\s*\[(\d+\.?\d*)%\]$', p)
        if match:
            results.append((match.group(1).strip(), float(match.group(2))))
        elif p:
            results.append((p.strip(), None))
    return results


# Load GEM->Compustat matches
parent_to_gvkeys = defaultdict(set)
with open(derived_path('mappings', 'gem_compustat_matches.csv'), 'r',
          encoding='utf-8') as f:
    for row in csv.DictReader(f):
        parent_to_gvkeys[row['gem_parent']].add(row['gvkey'])

# gvkey -> {country: total_mw}
gvkey_country_mw = defaultdict(lambda: defaultdict(float))

trackers = [
    ('gem_coal.csv', 'Parent'),
    ('gem_gas.csv', 'Parent(s)'),
    ('gem_solar.csv', 'Owner'),
    ('gem_wind.csv', 'Owner'),
]

for fname, parent_col in trackers:
    fpath = derived_path('gem', fname)
    _print(f'  Reading {fname}...')
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            status = row.get('Status', '')
            if status != 'operating':
                continue
            try:
                cap = float(row['Capacity (MW)'])
            except (ValueError, TypeError):
                continue
            country = row.get('Country/Area', '').strip()
            if not country:
                continue

            parsed = parse_parents(row.get(parent_col, ''))
            for name, pct in parsed:
                if name in parent_to_gvkeys:
                    share = ((pct / 100.0) if pct
                             else 1.0 / len(parsed) if len(parsed) > 1
                             else 1.0)
                    for gvkey in parent_to_gvkeys[name]:
                        gvkey_country_mw[gvkey][country] += cap * share

# Compute n_countries and single_country flag per gvkey
firm_n_countries = {}   # gvkey -> int
firm_is_single = {}     # gvkey -> bool (True if >90% MW in one country)
firm_total_mw = {}      # gvkey -> float

for gvkey, country_mw in gvkey_country_mw.items():
    total = sum(country_mw.values())
    firm_total_mw[gvkey] = total
    n_countries = len(country_mw)
    firm_n_countries[gvkey] = n_countries
    # Single-country: all plants in one country OR >90% MW in one country
    max_country_mw = max(country_mw.values())
    firm_is_single[gvkey] = (max_country_mw / total >= 0.90) if total > 0 else True

n_single = sum(1 for v in firm_is_single.values() if v)
n_multi = sum(1 for v in firm_is_single.values() if not v)
_print(f'\n  Firms with plant data: {len(gvkey_country_mw)}')
_print(f'  Single-country (>=90% MW in one country): {n_single}')
_print(f'  Multi-country: {n_multi}')

# Distribution of n_countries
from collections import Counter
nc_dist = Counter(firm_n_countries.values())
_print(f'\n  Distribution of n_countries:')
for nc in sorted(nc_dist.keys()):
    _print(f'    n_countries = {nc}: {nc_dist[nc]} firms')


# ══════════════════════════════════════════════════════════════════════
# STEP 2: Load shared data
# ══════════════════════════════════════════════════════════════════════

_print('\nLoading monthly returns...')
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r',
          encoding='utf-8') as f:
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
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r',
          encoding='utf-8') as f:
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
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r',
          encoding='utf-8') as f:
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
with open(derived_path('events', 'coal_retirement_events.csv'), 'r',
          encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        if row.get('is_first_mover') != 'True':
            continue
        ann_date = row.get('announcement_date', '').strip()
        ret_date = row.get('event_date', '').strip()
        effective_date = ann_date if ann_date else ret_date
        event_year = None
        if (effective_date and len(effective_date) >= 4
                and effective_date[:4].isdigit()):
            event_year = int(effective_date[:4])
        else:
            event_year = (int(row['ret_year'])
                          if row.get('ret_year') else None)
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
# STEP 3: Build per-event datasets (with n_countries attached)
# ══════════════════════════════════════════════════════════════════════

_print('\nBuilding per-event datasets...')

SPEC_VARS = ['w_fuel', 'w_geo', 'w_reg', 'same_sector']
MIN_OBS_PER_EVENT = 20

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
            nc = firm_n_countries.get(gk, 1)
            obs.append({
                'car': car,
                'w_fuel': w_fuel,
                'w_geo': w_geo,
                'w_reg': w_reg,
                'same_sector': same_sector,
                'gvkey': gk,
                'n_countries': nc,
                'is_single': firm_is_single.get(gk, True),
                'event_id': event_id,
            })

    if len(obs) >= MIN_OBS_PER_EVENT:
        event_datasets[event_id] = obs

n_valid = len(event_datasets)
total_obs = sum(len(v) for v in event_datasets.values())
_print(f'  Valid events (>= {MIN_OBS_PER_EVENT} obs): {n_valid}')
_print(f'  Total obs: {total_obs}')


# ══════════════════════════════════════════════════════════════════════
# HELPER: Run pooled OLS with event-clustered SEs
# ══════════════════════════════════════════════════════════════════════

def run_pooled_ols(obs_list, var_names, label):
    """Run pooled OLS on obs_list with given variable names.

    Returns dict with beta, se, t, p per variable, plus r2 and n.
    """
    _print(f'\n  --- Pooled OLS ({label}) ---')
    # Check same_sector variation
    if 'same_sector' in var_names:
        ss_vals = set(o['same_sector'] for o in obs_list)
        if len(ss_vals) <= 1:
            var_names = [v for v in var_names if v != 'same_sector']

    y = [o['car'] for o in obs_list]
    X = [[1.0] + [o[v] for v in var_names] for o in obs_list]
    cluster_ids = [o['event_id'] for o in obs_list]

    result = ols_simple(y, X)
    if result is None:
        _print('    OLS failed (singular)')
        return None

    ses = event_clustered_se(X, result['resid'], result['inv_XtX'],
                             cluster_ids)

    names = ['intercept'] + var_names
    out = {'r2': result['r2'], 'n': result['n'],
           'n_clusters': len(set(cluster_ids)), 'vars': {}}

    _print(f'    N = {result["n"]}, R2 = {result["r2"]:.4f}, '
           f'clusters = {len(set(cluster_ids))}')
    _print(f'    {"Variable":<25} {"beta":>12} {"SE(cl)":>10} '
           f'{"t":>8} {"p":>8}')
    _print('    ' + '-' * 65)

    for i, name in enumerate(names):
        b = result['beta'][i]
        se = ses[i]
        t = b / se if se > 1e-15 else 0.0
        p = p_from_t(t)
        stars = sig_stars(p)
        _print(f'    {name:<25} {b:+12.6f} {se:10.6f} {t:8.3f} '
               f'{p:8.4f}{stars}')
        out['vars'][name] = {'beta': b, 'se': se, 't': t, 'p': p}

    return out


# ══════════════════════════════════════════════════════════════════════
# HELPER: Run Fama-MacBeth with Newey-West SEs
# ══════════════════════════════════════════════════════════════════════

def run_fama_macbeth(datasets, var_names, label):
    """Run Fama-MacBeth: cross-sectional OLS per event, average betas.

    Returns dict with mean_beta, se, t, p per variable, plus T and avg_n.
    """
    _print(f'\n  --- Fama-MacBeth ({label}) ---')

    event_betas = defaultdict(list)
    event_r2s = []
    event_ns = []
    event_ids_used = []

    for eid in sorted(datasets.keys()):
        obs = datasets[eid]
        use_vars = list(var_names)
        if 'same_sector' in use_vars:
            ss_vals = set(o['same_sector'] for o in obs)
            if len(ss_vals) <= 1:
                use_vars = [v for v in use_vars if v != 'same_sector']

        y = [o['car'] for o in obs]
        X = [[1.0] + [o[v] for v in use_vars] for o in obs]

        result = ols_simple(y, X)
        if result is None:
            continue

        names = ['intercept'] + use_vars
        for i, name in enumerate(names):
            event_betas[name].append(result['beta'][i])
        # Pad missing variables with NaN
        for v in var_names:
            if v not in use_vars:
                event_betas[v].append(float('nan'))

        event_r2s.append(result['r2'])
        event_ns.append(result['n'])
        event_ids_used.append(eid)

    T_fm = len(event_ids_used)
    if T_fm < 3:
        _print(f'    Too few events ({T_fm})')
        return None

    avg_n = sum(event_ns) / T_fm
    avg_r2 = sum(event_r2s) / T_fm
    _print(f'    T = {T_fm} events, avg N = {avg_n:.1f}, '
           f'avg R2 = {avg_r2:.4f}')
    _print(f'    {"Variable":<25} {"Mean beta":>12} {"NW SE":>10} '
           f'{"t":>8} {"p":>8}')
    _print('    ' + '-' * 65)

    out = {'T': T_fm, 'avg_n': avg_n, 'avg_r2': avg_r2, 'vars': {}}

    for v in ['intercept'] + list(var_names):
        betas = event_betas.get(v, [])
        clean = [b for b in betas if not math.isnan(b)]
        if len(clean) < 3:
            continue
        mean_b = sum(clean) / len(clean)
        nw_se = newey_west_se(clean)
        t_stat = mean_b / nw_se if nw_se > 1e-15 else 0.0
        p_val = p_from_t(t_stat)
        stars = sig_stars(p_val)
        _print(f'    {v:<25} {mean_b:+12.6f} {nw_se:10.6f} {t_stat:8.3f} '
               f'{p_val:8.4f}{stars}')
        out['vars'][v] = {'beta': mean_b, 'se': nw_se, 't': t_stat,
                          'p': p_val, 'n_events': len(clean)}

    return out


# ══════════════════════════════════════════════════════════════════════
# SPEC 1: Full sample baseline
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('SPECIFICATION 1: FULL SAMPLE BASELINE')
_print('CAR = b1*w_fuel + b2*w_geo + b3*w_reg + b4*same_sector + e')
_print('=' * 70)

all_obs = []
for eid in sorted(event_datasets.keys()):
    all_obs.extend(event_datasets[eid])

pooled_full = run_pooled_ols(all_obs, SPEC_VARS, 'Full sample')
fm_full = run_fama_macbeth(event_datasets, SPEC_VARS, 'Full sample')


# ══════════════════════════════════════════════════════════════════════
# SPEC 2: Single-country subsample
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('SPECIFICATION 2: SINGLE-COUNTRY SUBSAMPLE')
_print('Same spec, but only firms with >=90% MW in one country')
_print('=' * 70)

# Filter event datasets to single-country firms only
single_datasets = {}
for eid, obs in event_datasets.items():
    filtered = [o for o in obs if o['is_single']]
    if len(filtered) >= MIN_OBS_PER_EVENT:
        single_datasets[eid] = filtered

n_single_events = len(single_datasets)
n_single_obs = sum(len(v) for v in single_datasets.values())
_print(f'  Events with >= {MIN_OBS_PER_EVENT} single-country obs: '
       f'{n_single_events}')
_print(f'  Total single-country obs: {n_single_obs}')

# Count single-country firms in each event
single_counts = []
for eid, obs in event_datasets.items():
    n_s = sum(1 for o in obs if o['is_single'])
    single_counts.append(n_s)
if single_counts:
    _print(f'  Across all events: mean single-country firms = '
           f'{sum(single_counts)/len(single_counts):.1f}, '
           f'min = {min(single_counts)}, max = {max(single_counts)}')

all_single_obs = []
for eid in sorted(single_datasets.keys()):
    all_single_obs.extend(single_datasets[eid])

pooled_single = run_pooled_ols(all_single_obs, SPEC_VARS,
                                'Single-country subsample')
fm_single = run_fama_macbeth(single_datasets, SPEC_VARS,
                              'Single-country subsample')


# Also run multi-country subsample for comparison
_print('\n  --- Multi-country subsample (for comparison) ---')
multi_datasets = {}
for eid, obs in event_datasets.items():
    filtered = [o for o in obs if not o['is_single']]
    if len(filtered) >= MIN_OBS_PER_EVENT:
        multi_datasets[eid] = filtered

n_multi_events = len(multi_datasets)
n_multi_obs = sum(len(v) for v in multi_datasets.values())
_print(f'  Events with >= {MIN_OBS_PER_EVENT} multi-country obs: '
       f'{n_multi_events}')
_print(f'  Total multi-country obs: {n_multi_obs}')

all_multi_obs = []
for eid in sorted(multi_datasets.keys()):
    all_multi_obs.extend(multi_datasets[eid])

pooled_multi = None
fm_multi = None
if all_multi_obs:
    pooled_multi = run_pooled_ols(all_multi_obs, SPEC_VARS,
                                   'Multi-country subsample')
if multi_datasets:
    fm_multi = run_fama_macbeth(multi_datasets, SPEC_VARS,
                                 'Multi-country subsample')


# ══════════════════════════════════════════════════════════════════════
# SPEC 3: Diversification interaction
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('SPECIFICATION 3: DIVERSIFICATION INTERACTION')
_print('CAR = b1*w_fuel + b2*w_geo + b3*w_reg + b4*same_sector '
       '+ b5*(w_geo x n_countries) + e')
_print('=' * 70)

# Add interaction terms to observations
for eid, obs in event_datasets.items():
    for o in obs:
        nc = o['n_countries']
        o['w_geo_x_nc'] = o['w_geo'] * nc
        o['w_geo_x_log_nc'] = o['w_geo'] * math.log(max(nc, 1))

# Spec 3a: w_geo x n_countries
SPEC_3A_VARS = ['w_fuel', 'w_geo', 'w_reg', 'same_sector', 'w_geo_x_nc']

_print('\n--- Spec 3a: w_geo x n_countries ---')
pooled_3a = run_pooled_ols(all_obs, SPEC_3A_VARS, 'Interaction: n_countries')
fm_3a = run_fama_macbeth(event_datasets, SPEC_3A_VARS,
                          'Interaction: n_countries')

# Spec 3b: w_geo x log(n_countries)
SPEC_3B_VARS = ['w_fuel', 'w_geo', 'w_reg', 'same_sector', 'w_geo_x_log_nc']

_print('\n--- Spec 3b: w_geo x log(n_countries) ---')
pooled_3b = run_pooled_ols(all_obs, SPEC_3B_VARS,
                            'Interaction: log(n_countries)')
fm_3b = run_fama_macbeth(event_datasets, SPEC_3B_VARS,
                          'Interaction: log(n_countries)')


# ══════════════════════════════════════════════════════════════════════
# ADDITIONAL DIAGNOSTICS
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('DIAGNOSTICS')
_print('=' * 70)

# Distribution of n_countries across event obs
all_nc = [o['n_countries'] for o in all_obs]
nc_counter = Counter(all_nc)
_print(f'\n  n_countries distribution across event observations:')
for nc in sorted(nc_counter.keys())[:15]:
    pct = 100.0 * nc_counter[nc] / len(all_nc)
    _print(f'    n_countries = {nc}: {nc_counter[nc]} obs ({pct:.1f}%)')

# Mean geo exposure by single/multi
single_geo = [o['w_geo'] for o in all_obs if o['is_single']]
multi_geo = [o['w_geo'] for o in all_obs if not o['is_single']]
if single_geo and multi_geo:
    _print(f'\n  Mean w_geo (single-country): {sum(single_geo)/len(single_geo):.6f}')
    _print(f'  Mean w_geo (multi-country):  {sum(multi_geo)/len(multi_geo):.6f}')

# Mean CAR by single/multi
single_car = [o['car'] for o in all_obs if o['is_single']]
multi_car = [o['car'] for o in all_obs if not o['is_single']]
if single_car and multi_car:
    _print(f'\n  Mean CAR (single-country): {sum(single_car)/len(single_car):.6f}')
    _print(f'  Mean CAR (multi-country):  {sum(multi_car)/len(multi_car):.6f}')

# Correlation between n_countries and w_geo across obs
if all_nc:
    mean_nc = sum(all_nc) / len(all_nc)
    mean_wg = sum(o['w_geo'] for o in all_obs) / len(all_obs)
    cov_nc_wg = sum((o['n_countries'] - mean_nc) * (o['w_geo'] - mean_wg)
                     for o in all_obs) / len(all_obs)
    var_nc = sum((nc - mean_nc) ** 2 for nc in all_nc) / len(all_nc)
    var_wg = sum((o['w_geo'] - mean_wg) ** 2 for o in all_obs) / len(all_obs)
    if var_nc > 0 and var_wg > 0:
        corr = cov_nc_wg / math.sqrt(var_nc * var_wg)
        _print(f'\n  Correlation(n_countries, w_geo): {corr:.4f}')


# ══════════════════════════════════════════════════════════════════════
# WRITE OUTPUT
# ══════════════════════════════════════════════════════════════════════

out_path = results_path('metrics', 'strategy2_geo_diversification.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)


def fmt_row(var_name, pooled_res, fm_res):
    """Format a table row for one variable across pooled and FM results."""
    pv = pooled_res['vars'].get(var_name, {}) if pooled_res else {}
    fv = fm_res['vars'].get(var_name, {}) if fm_res else {}
    pb = pv.get('beta', float('nan'))
    pse = pv.get('se', float('nan'))
    pt = pv.get('t', float('nan'))
    pp = pv.get('p', float('nan'))
    fb = fv.get('beta', float('nan'))
    fse = fv.get('se', float('nan'))
    ft = fv.get('t', float('nan'))
    fp = fv.get('p', float('nan'))

    def f(x, fmt='+.4f'):
        if math.isnan(x) or math.isinf(x):
            return '--'
        return format(x, fmt)

    def stars(p_val):
        if math.isnan(p_val) or math.isinf(p_val):
            return ''
        return sig_stars(p_val)

    return (f'| {var_name:<22} '
            f'| {f(pb)} | ({f(pse, ".4f")}) | {f(pt, ".2f")} | {f(pp, ".3f")}{stars(pp)} '
            f'| {f(fb)} | ({f(fse, ".4f")}) | {f(ft, ".2f")} | {f(fp, ".3f")}{stars(fp)} |')


lines = [
    '# Geographic Diversification Test',
    '',
    'Tests whether multinational diversification attenuates the geographic',
    'competitive channel (w_geo). If geographic effects are real at the plant',
    'level, single-country utilities should show a stronger geo coefficient.',
    '',
    f'Events: {len(all_events)} first-mover coal retirements',
    f'Window: [-1, +{POST_MONTHS}] months, vwretd market-adjusted returns',
    f'Minimum obs per event: {MIN_OBS_PER_EVENT}',
    f'Single-country definition: >=90% of MW in one country',
    '',
    '## Firm-Level Diversification Summary',
    '',
    f'- Firms with plant data: {len(gvkey_country_mw)}',
    f'- Single-country firms: {n_single}',
    f'- Multi-country firms: {n_multi}',
    '',
    '| n_countries | Firms |',
    '|---:|---:|',
]

for nc in sorted(nc_dist.keys()):
    lines.append(f'| {nc} | {nc_dist[nc]} |')

lines += [
    '',
    '## Specification 1: Full Sample Baseline',
    '',
    'CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector + e',
    '',
]

if pooled_full and fm_full:
    lines.append(f'N = {pooled_full["n"]}, Events = {fm_full["T"]}, '
                 f'R2(pooled) = {pooled_full["r2"]:.4f}, '
                 f'R2(FM avg) = {fm_full["avg_r2"]:.4f}')
    lines.append('')
    lines.append('| Variable | beta(OLS) | SE(cl) | t | p '
                 '| beta(FM) | SE(NW) | t | p |')
    lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for v in SPEC_VARS:
        lines.append(fmt_row(v, pooled_full, fm_full))

lines += [
    '',
    '## Specification 2: Single-Country Subsample',
    '',
    'Same specification, restricted to firms with >=90% MW in one country.',
    '',
    f'Events with >= {MIN_OBS_PER_EVENT} single-country obs: '
    f'{n_single_events}',
    f'Total single-country obs: {n_single_obs}',
    '',
]

if pooled_single and fm_single:
    lines.append(f'N = {pooled_single["n"]}, Events = {fm_single["T"]}, '
                 f'R2(pooled) = {pooled_single["r2"]:.4f}, '
                 f'R2(FM avg) = {fm_single["avg_r2"]:.4f}')
    lines.append('')
    lines.append('| Variable | beta(OLS) | SE(cl) | t | p '
                 '| beta(FM) | SE(NW) | t | p |')
    lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for v in SPEC_VARS:
        lines.append(fmt_row(v, pooled_single, fm_single))

lines += [
    '',
    '### Multi-Country Subsample (for comparison)',
    '',
]

if pooled_multi and fm_multi:
    lines.append(f'N = {pooled_multi["n"]}, Events = {fm_multi["T"]}, '
                 f'R2(pooled) = {pooled_multi["r2"]:.4f}, '
                 f'R2(FM avg) = {fm_multi["avg_r2"]:.4f}')
    lines.append('')
    lines.append('| Variable | beta(OLS) | SE(cl) | t | p '
                 '| beta(FM) | SE(NW) | t | p |')
    lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for v in SPEC_VARS:
        lines.append(fmt_row(v, pooled_multi, fm_multi))
elif n_multi_events == 0:
    lines.append(f'Too few events with >= {MIN_OBS_PER_EVENT} '
                 'multi-country obs.')

lines += [
    '',
    '## Specification 3a: Diversification Interaction (n_countries)',
    '',
    'CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector '
    '+ b5 (w_geo x n_countries) + e',
    '',
    'If b2 < 0 and b5 > 0: geo benefit exists but weakens with diversification',
    '',
]

if pooled_3a and fm_3a:
    lines.append(f'N = {pooled_3a["n"]}, Events = {fm_3a["T"]}, '
                 f'R2(pooled) = {pooled_3a["r2"]:.4f}, '
                 f'R2(FM avg) = {fm_3a["avg_r2"]:.4f}')
    lines.append('')
    lines.append('| Variable | beta(OLS) | SE(cl) | t | p '
                 '| beta(FM) | SE(NW) | t | p |')
    lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for v in SPEC_3A_VARS:
        lines.append(fmt_row(v, pooled_3a, fm_3a))

lines += [
    '',
    '## Specification 3b: Diversification Interaction (log n_countries)',
    '',
    'CAR = b1 w_fuel + b2 w_geo + b3 w_reg + b4 same_sector '
    '+ b5 (w_geo x log(n_countries)) + e',
    '',
    'Log version: captures concave attenuation '
    '(first additional country matters most)',
    '',
]

if pooled_3b and fm_3b:
    lines.append(f'N = {pooled_3b["n"]}, Events = {fm_3b["T"]}, '
                 f'R2(pooled) = {pooled_3b["r2"]:.4f}, '
                 f'R2(FM avg) = {fm_3b["avg_r2"]:.4f}')
    lines.append('')
    lines.append('| Variable | beta(OLS) | SE(cl) | t | p '
                 '| beta(FM) | SE(NW) | t | p |')
    lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for v in SPEC_3B_VARS:
        lines.append(fmt_row(v, pooled_3b, fm_3b))

# Interpretation
lines += [
    '',
    '## Interpretation',
    '',
]

# Summarize key findings
if fm_full and fm_single:
    fg = fm_full['vars'].get('w_geo', {})
    sg = fm_single['vars'].get('w_geo', {})
    if fg and sg:
        lines.append(f'- Full sample w_geo: beta = {fg.get("beta", 0):+.6f}, '
                     f't = {fg.get("t", 0):.3f}, p = {fg.get("p", 1):.4f}')
        lines.append(f'- Single-country w_geo: beta = {sg.get("beta", 0):+.6f}, '
                     f't = {sg.get("t", 0):.3f}, p = {sg.get("p", 1):.4f}')
        if abs(sg.get('t', 0)) > abs(fg.get('t', 0)):
            lines.append('- **w_geo is stronger in the single-country subsample**, '
                         'consistent with aggregation attenuation.')
        else:
            lines.append('- w_geo is NOT stronger in the single-country subsample.')

if fm_3a:
    ia = fm_3a['vars'].get('w_geo_x_nc', {})
    g3 = fm_3a['vars'].get('w_geo', {})
    if ia and g3:
        lines.append(f'- Interaction (w_geo x n_countries): '
                     f'beta = {ia.get("beta", 0):+.6f}, '
                     f't = {ia.get("t", 0):.3f}, p = {ia.get("p", 1):.4f}')
        if g3.get('beta', 0) < 0 and ia.get('beta', 0) > 0:
            lines.append('  Pattern: negative base geo + positive interaction '
                         '-> geo benefit weakens with more countries.')

if fm_3b:
    ib = fm_3b['vars'].get('w_geo_x_log_nc', {})
    if ib:
        lines.append(f'- Interaction (w_geo x log(n_countries)): '
                     f'beta = {ib.get("beta", 0):+.6f}, '
                     f't = {ib.get("t", 0):.3f}, p = {ib.get("p", 1):.4f}')

# Diagnostics section
lines += [
    '',
    '## Diagnostics',
    '',
]

if single_geo and multi_geo:
    lines.append(f'- Mean w_geo (single-country): '
                 f'{sum(single_geo)/len(single_geo):.6f}')
    lines.append(f'- Mean w_geo (multi-country): '
                 f'{sum(multi_geo)/len(multi_geo):.6f}')

if single_car and multi_car:
    lines.append(f'- Mean CAR (single-country): '
                 f'{sum(single_car)/len(single_car):.6f}')
    lines.append(f'- Mean CAR (multi-country): '
                 f'{sum(multi_car)/len(multi_car):.6f}')

if all_nc:
    mean_nc_val = sum(all_nc) / len(all_nc)
    lines.append(f'- Mean n_countries across obs: {mean_nc_val:.2f}')

lines.append('')

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

_print(f'\nWrote: {out_path}')
_print('Done.')
