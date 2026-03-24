"""ESG horse race and ETS interaction with Fama-MacBeth + Newey-West inference.

Re-tests two key results from the channel decomposition paper:
1. ESG horse race: spatial exposure (t=-2.86) vs ESG score (t=-0.59)
2. ETS interaction: w_fuel x has_ets = -3.242 (t=-3.41)

Both were originally tested with event-clustered SEs. This script applies
Fama-MacBeth (1973) cross-sectional regressions with Newey-West (1987)
HAC standard errors on the time series of betas, which properly accounts
for cross-event correlation.

Output: results/metrics/strategy2_esg_ets_fmb.md
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


# ── Load ESG scores ──────────────────────────────────────────────────

_print('Loading ESG scores...')
esg_scores = {}  # gvkey -> env_score (normalized to 0-1)
esg_path = raw_path('refinitiv', 'refinitiv_esg.csv')
if os.path.exists(esg_path):
    with open(esg_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row.get('gvkey', '').strip()
            score_str = row.get('env_score', '').strip()
            if not gk or not score_str:
                continue
            try:
                score = float(score_str)
            except ValueError:
                continue
            # Normalize to [0,1] if on 0-100 scale
            if score > 1.0:
                score = score / 100.0
            esg_scores[gk] = score
    _print(f'  ESG scores: {len(esg_scores)} firms')
else:
    _print('  WARNING: ESG file not found')

# ── Load ETS membership ──────────────────────────────────────────────

_print('Loading ETS membership...')
ets_member = {}  # gvkey -> 0 or 1
ets_path = derived_path('networks', 'firm_ets_membership.csv')
if os.path.exists(ets_path):
    with open(ets_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row.get('gvkey', '').strip()
            has = row.get('has_ets', '0').strip()
            if gk:
                ets_member[gk] = 1 if has == '1' else 0
    _print(f'  ETS members: {sum(v for v in ets_member.values())} of {len(ets_member)} firms')
else:
    _print('  WARNING: ETS file not found')


# ── Load events ──────────────────────────────────────────────────────

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

MIN_OBS_PER_EVENT = 20

event_datasets = {}  # event_id -> list of obs dicts (ALL firms)
event_datasets_esg = {}  # event_id -> list of obs dicts (ESG subset)

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

    obs_all = []
    obs_esg = []
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
            has_ets = float(ets_member.get(gk, 0))
            car = compute_monthly_car(gk, event_month, post=POST_MONTHS)
            if car is None:
                continue

            ob = {
                'car': car,
                'w_geo': w_geo,
                'w_fuel': w_fuel,
                'w_reg': w_reg,
                'same_sector': same_sector,
                'has_ets': has_ets,
                'w_fuel_ets': w_fuel * has_ets,
                'w_geo_ets': w_geo * has_ets,
                'gvkey': gk,
            }
            obs_all.append(ob)

            # ESG subset
            if gk in esg_scores:
                ob_esg = dict(ob)
                ob_esg['esg_score'] = esg_scores[gk]
                obs_esg.append(ob_esg)

    if len(obs_all) >= MIN_OBS_PER_EVENT:
        event_datasets[event_id] = obs_all
    if len(obs_esg) >= MIN_OBS_PER_EVENT:
        event_datasets_esg[event_id] = obs_esg

_print(f'  Valid events (all firms, >= {MIN_OBS_PER_EVENT} obs): {len(event_datasets)}')
_print(f'  Valid events (ESG subset, >= {MIN_OBS_PER_EVENT} obs): {len(event_datasets_esg)}')
_print(f'  Total obs (all): {sum(len(v) for v in event_datasets.values())}')
_print(f'  Total obs (ESG): {sum(len(v) for v in event_datasets_esg.values())}')


# ── Helper: run FM across events for a given spec ────────────────────

def run_fama_macbeth(datasets, var_names, label):
    """Run FM regressions across events. Returns dict of results per variable."""
    event_betas = defaultdict(list)
    event_r2s = []
    event_ns = []
    events_used = []

    for event_id in sorted(datasets.keys()):
        obs = datasets[event_id]
        n_obs = len(obs)

        # Check which vars have variation
        use_vars = []
        for v in var_names:
            vals = set(o[v] for o in obs)
            if len(vals) > 1:
                use_vars.append(v)
            # If no variation, skip this variable for this event

        if len(use_vars) == 0:
            continue

        y = [o['car'] for o in obs]
        X = [[1.0] + [o[v] for v in use_vars] for o in obs]

        result = ols_simple(y, X)
        if result is None:
            continue

        names = ['intercept'] + use_vars
        for i, name in enumerate(names):
            event_betas[name].append(result['beta'][i])

        # Pad missing vars with NaN
        for v in var_names:
            if v not in use_vars:
                event_betas[v].append(float('nan'))

        event_r2s.append(result['r2'])
        event_ns.append(result['n'])
        events_used.append(event_id)

    T_fm = len(events_used)
    if T_fm < 3:
        _print(f'  [{label}] Too few events ({T_fm}), skipping.')
        return None

    _print(f'\n  [{label}] Events: {T_fm}, Avg N: {sum(event_ns)/T_fm:.1f}, '
           f'Avg R2: {sum(event_r2s)/T_fm:.4f}')

    fm_results = {}
    for v in ['intercept'] + var_names:
        betas = event_betas[v]
        clean = [b for b in betas if not math.isnan(b)]
        if len(clean) < 3:
            continue
        mean_b = sum(clean) / len(clean)
        nw_se = newey_west_se(clean)
        t_stat = mean_b / nw_se if nw_se > 1e-15 else 0.0
        p_val = p_from_t(t_stat)
        stars = '***' if p_val < 0.01 else '**' if p_val < 0.05 else '*' if p_val < 0.10 else ''
        _print(f'    {v:<18} {mean_b:+12.6f} {nw_se:10.6f} {t_stat:8.3f} {p_val:8.4f}{stars}')
        fm_results[v] = {
            'mean': mean_b, 'se': nw_se, 't': t_stat, 'p': p_val,
            'n_events': len(clean), 'betas': clean,
        }

    return {'results': fm_results, 'T': T_fm, 'avg_n': sum(event_ns)/T_fm,
            'avg_r2': sum(event_r2s)/T_fm, 'event_betas': event_betas}


# ══════════════════════════════════════════════════════════════════════
# PANEL A: ESG HORSE RACE (restricted to ESG-covered firms)
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('PANEL A: ESG HORSE RACE (Fama-MacBeth + Newey-West)')
_print('Sample restricted to firms with Refinitiv ESG scores')
_print('=' * 70)

_print(f'\n  {"Variable":<18} {"Mean beta":>12} {"NW SE":>10} {"t":>8} {"p":>8}')
_print('  ' + '-' * 60)

# Spec (1): CAR = a + b1*esg_score
esg_spec1 = run_fama_macbeth(event_datasets_esg,
                              ['esg_score'],
                              'ESG-1: ESG only')

# Spec (2): CAR = a + b1*w_fuel + b2*w_geo + b3*same_sector
esg_spec2 = run_fama_macbeth(event_datasets_esg,
                              ['w_fuel', 'w_geo', 'same_sector'],
                              'ESG-2: Spatial only')

# Spec (3): CAR = a + b1*w_fuel + b2*w_geo + b3*esg_score + b4*same_sector
esg_spec3 = run_fama_macbeth(event_datasets_esg,
                              ['w_fuel', 'w_geo', 'esg_score', 'same_sector'],
                              'ESG-3: Horse race')

# Spec (4): Full model with ETS interaction
esg_spec4 = run_fama_macbeth(event_datasets_esg,
                              ['w_fuel', 'w_geo', 'w_reg', 'esg_score',
                               'w_fuel_ets', 'same_sector'],
                              'ESG-4: Full + ETS')


# ══════════════════════════════════════════════════════════════════════
# PANEL B: ETS INTERACTION (all firms)
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('PANEL B: ETS INTERACTION (Fama-MacBeth + Newey-West)')
_print('Full sample (all firms)')
_print('=' * 70)

_print(f'\n  {"Variable":<18} {"Mean beta":>12} {"NW SE":>10} {"t":>8} {"p":>8}')
_print('  ' + '-' * 60)

# Spec (1): CAR = a + b1*w_fuel + b2*(w_fuel x has_ets) + b3*w_geo + b4*w_reg + b5*same_sector
ets_spec1 = run_fama_macbeth(event_datasets,
                              ['w_fuel', 'w_fuel_ets', 'w_geo', 'w_reg', 'same_sector'],
                              'ETS-1: Fuel x ETS')

# Spec (2): Placebo: also include w_geo x has_ets
ets_spec2 = run_fama_macbeth(event_datasets,
                              ['w_fuel', 'w_fuel_ets', 'w_geo', 'w_geo_ets',
                               'w_reg', 'same_sector'],
                              'ETS-2: Placebo (Geo x ETS)')


# ══════════════════════════════════════════════════════════════════════
# DIFFERENCE TEST: (fuel x ETS) - (geo x ETS)
# ══════════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('DIFFERENCE TEST: beta(w_fuel x ETS) - beta(w_geo x ETS)')
_print('=' * 70)

if ets_spec2 and 'w_fuel_ets' in ets_spec2['results'] and 'w_geo_ets' in ets_spec2['results']:
    fuel_ets_betas = ets_spec2['results']['w_fuel_ets']['betas']
    geo_ets_betas = ets_spec2['results']['w_geo_ets']['betas']
    T_diff = min(len(fuel_ets_betas), len(geo_ets_betas))
    diff_series = [fuel_ets_betas[t] - geo_ets_betas[t] for t in range(T_diff)]
    mean_diff = sum(diff_series) / T_diff
    nw_se_diff = newey_west_se(diff_series)
    t_diff = mean_diff / nw_se_diff if nw_se_diff > 1e-15 else 0.0
    p_diff = p_from_t(t_diff)
    stars_d = '***' if p_diff < 0.01 else '**' if p_diff < 0.05 else '*' if p_diff < 0.10 else ''
    _print(f'  Mean diff = {mean_diff:+.6f}, NW SE = {nw_se_diff:.6f}, '
           f't = {t_diff:.3f}, p = {p_diff:.4f}{stars_d}')
else:
    _print('  Could not compute difference test (missing betas).')
    mean_diff = t_diff = p_diff = None
    stars_d = ''


# ══════════════════════════════════════════════════════════════════════
# WRITE OUTPUT
# ══════════════════════════════════════════════════════════════════════

_print('\nWriting output...')

out_path = results_path('metrics', 'strategy2_esg_ets_fmb.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)


def fmt_row(var, res):
    """Format a result row for markdown table."""
    if res is None:
        return f'| {var} | | | | |'
    stars = '***' if res['p'] < 0.01 else '**' if res['p'] < 0.05 else '*' if res['p'] < 0.10 else ''
    return (f'| {var} | {res["mean"]:+.4f} | {res["se"]:.4f} | '
            f'{res["t"]:.2f}{stars} | {res["p"]:.4f} |')


def write_spec_table(f, spec, var_names, spec_label, spec_info):
    """Write a spec's results as a markdown table."""
    if spec is None:
        f.write(f'\n**{spec_label}**: Too few events to estimate.\n\n')
        return
    res = spec['results']
    f.write(f'\n**{spec_label}** (T={spec["T"]}, '
            f'avg N={spec["avg_n"]:.0f}, avg R2={spec["avg_r2"]:.4f})\n\n')
    f.write('| Variable | Mean beta | NW SE | t-stat | p-value |\n')
    f.write('|:---------|----------:|------:|-------:|--------:|\n')
    for v in ['intercept'] + var_names:
        if v in res:
            f.write(fmt_row(v, res[v]) + '\n')
    f.write('\n')


with open(out_path, 'w', encoding='utf-8') as f:
    f.write('# ESG Horse Race and ETS Interaction: Fama-MacBeth + Newey-West\n\n')
    f.write('Re-tests two key results with robust inference that properly\n')
    f.write('accounts for cross-event correlation in the time series of betas.\n\n')
    f.write('**Method**: For each event, run cross-sectional OLS. Average betas\n')
    f.write('across events. Compute Newey-West (1987) HAC SEs with Bartlett kernel\n')
    f.write('and automatic lag selection (floor(4*(T/100)^{2/9})).\n\n')
    f.write(f'Significance: \\*p<0.10, \\*\\*p<0.05, \\*\\*\\*p<0.01\n\n')

    # Panel A
    f.write('---\n\n')
    f.write('## Panel A: ESG Horse Race\n\n')
    f.write('Sample restricted to firms with Refinitiv environmental scores.\n')
    f.write(f'ESG coverage: {len(esg_scores)} firms.\n\n')

    f.write('**Original results** (event-clustered SEs):\n')
    f.write('- Spatial exposure w_fuel: t = -2.86\n')
    f.write('- ESG score: t = -0.59 (insignificant)\n\n')

    write_spec_table(f, esg_spec1, ['esg_score'],
                     'Spec (1): ESG score only', None)
    write_spec_table(f, esg_spec2, ['w_fuel', 'w_geo', 'same_sector'],
                     'Spec (2): Spatial exposure only', None)
    write_spec_table(f, esg_spec3, ['w_fuel', 'w_geo', 'esg_score', 'same_sector'],
                     'Spec (3): Horse race (spatial + ESG)', None)
    write_spec_table(f, esg_spec4, ['w_fuel', 'w_geo', 'w_reg', 'esg_score',
                                     'w_fuel_ets', 'same_sector'],
                     'Spec (4): Full model + ETS interaction', None)

    # Summary for Panel A
    f.write('### Panel A Summary\n\n')
    if esg_spec3 and 'w_fuel' in esg_spec3['results'] and 'esg_score' in esg_spec3['results']:
        fuel_r = esg_spec3['results']['w_fuel']
        esg_r = esg_spec3['results']['esg_score']
        f.write(f'In the horse race (Spec 3), spatial fuel exposure has '
                f't={fuel_r["t"]:.2f} (p={fuel_r["p"]:.4f}) while ESG score '
                f'has t={esg_r["t"]:.2f} (p={esg_r["p"]:.4f}).\n\n')
        if abs(fuel_r['t']) > 1.96 and abs(esg_r['t']) < 1.96:
            f.write('**Result**: Spatial exposure survives FM+NW; ESG does not. '
                    'The original finding is confirmed.\n\n')
        elif abs(fuel_r['t']) < 1.96 and abs(esg_r['t']) < 1.96:
            f.write('**Result**: Neither survives FM+NW inference. The original '
                    'spatial t-stat was inflated by event clustering.\n\n')
        elif abs(fuel_r['t']) < 1.96:
            f.write('**Result**: Spatial exposure no longer significant under FM+NW, '
                    'but ESG remains insignificant too.\n\n')
        else:
            f.write('**Result**: Both survive FM+NW inference.\n\n')
    else:
        f.write('Insufficient data for summary.\n\n')

    # Panel B
    f.write('---\n\n')
    f.write('## Panel B: ETS Interaction\n\n')
    f.write('Full sample (all firms, not restricted to ESG coverage).\n')
    f.write(f'ETS coverage: {sum(v for v in ets_member.values())} firms under ETS '
            f'of {len(ets_member)} total.\n\n')

    f.write('**Original result** (event-clustered SEs):\n')
    f.write('- w_fuel x has_ets: coeff = -3.242, t = -3.41\n\n')

    write_spec_table(f, ets_spec1, ['w_fuel', 'w_fuel_ets', 'w_geo', 'w_reg', 'same_sector'],
                     'Spec (1): Fuel x ETS interaction', None)
    write_spec_table(f, ets_spec2, ['w_fuel', 'w_fuel_ets', 'w_geo', 'w_geo_ets',
                                     'w_reg', 'same_sector'],
                     'Spec (2): Placebo (Geo x ETS added)', None)

    # Difference test
    f.write('### Difference Test: beta(w_fuel x ETS) - beta(w_geo x ETS)\n\n')
    if t_diff is not None:
        f.write(f'Mean difference = {mean_diff:+.4f}, NW SE = {nw_se_diff:.4f}, '
                f't = {t_diff:.2f}, p = {p_diff:.4f}{stars_d}\n\n')
        if abs(t_diff) > 1.96:
            f.write('The fuel-mix channel through ETS is significantly stronger '
                    'than the geographic channel through ETS.\n\n')
        else:
            f.write('The difference is not statistically significant.\n\n')
    else:
        f.write('Could not compute difference test.\n\n')

    # Panel B summary
    f.write('### Panel B Summary\n\n')
    if ets_spec1 and 'w_fuel_ets' in ets_spec1['results']:
        ets_r = ets_spec1['results']['w_fuel_ets']
        f.write(f'The fuel x ETS interaction has FM+NW t={ets_r["t"]:.2f} '
                f'(p={ets_r["p"]:.4f})')
        orig_t = -3.41
        reduction = (1 - abs(ets_r['t']) / abs(orig_t)) * 100
        f.write(f', a {reduction:.0f}% reduction from the original t={orig_t:.2f}.\n\n')
        if abs(ets_r['t']) > 1.96:
            f.write('**Result**: The ETS interaction survives FM+NW inference.\n\n')
        elif abs(ets_r['t']) > 1.645:
            f.write('**Result**: The ETS interaction is marginally significant (p<0.10) '
                    'under FM+NW.\n\n')
        else:
            f.write('**Result**: The ETS interaction does NOT survive FM+NW inference.\n\n')
    else:
        f.write('Insufficient data for summary.\n\n')

    # Placebo
    if ets_spec2 and 'w_geo_ets' in ets_spec2['results']:
        geo_ets_r = ets_spec2['results']['w_geo_ets']
        f.write(f'**Placebo**: w_geo x has_ets has t={geo_ets_r["t"]:.2f} '
                f'(p={geo_ets_r["p"]:.4f}). ')
        if abs(geo_ets_r['t']) < 1.96:
            f.write('The placebo is insignificant, as expected: ETS membership '
                    'amplifies the fuel-mix channel, not geographic proximity.\n\n')
        else:
            f.write('The placebo is significant, which complicates the story.\n\n')

    # Overall
    f.write('---\n\n')
    f.write('## Overall Assessment\n\n')
    f.write('| Test | Original t | FM+NW t | Survives? |\n')
    f.write('|:-----|----------:|---------:|:---------:|\n')

    if esg_spec3 and 'w_fuel' in esg_spec3['results']:
        fuel_t = esg_spec3['results']['w_fuel']['t']
        surv = 'Yes' if abs(fuel_t) > 1.96 else 'Marginal' if abs(fuel_t) > 1.645 else 'No'
        f.write(f'| Spatial exposure (ESG race) | -2.86 | {fuel_t:.2f} | {surv} |\n')
    if esg_spec3 and 'esg_score' in esg_spec3['results']:
        esg_t = esg_spec3['results']['esg_score']['t']
        surv = 'Yes' if abs(esg_t) > 1.96 else 'Marginal' if abs(esg_t) > 1.645 else 'No'
        f.write(f'| ESG score (horse race) | -0.59 | {esg_t:.2f} | {surv} |\n')
    if ets_spec1 and 'w_fuel_ets' in ets_spec1['results']:
        ets_t = ets_spec1['results']['w_fuel_ets']['t']
        surv = 'Yes' if abs(ets_t) > 1.96 else 'Marginal' if abs(ets_t) > 1.645 else 'No'
        f.write(f'| w_fuel x has_ets | -3.41 | {ets_t:.2f} | {surv} |\n')
    if ets_spec2 and 'w_geo_ets' in ets_spec2['results']:
        geoets_t = ets_spec2['results']['w_geo_ets']['t']
        surv = 'Yes' if abs(geoets_t) > 1.96 else 'Marginal' if abs(geoets_t) > 1.645 else 'No'
        f.write(f'| w_geo x has_ets (placebo) | n/a | {geoets_t:.2f} | {surv} |\n')

    f.write('\n')

_print(f'\nOutput written to {out_path}')
_print('Done.')
