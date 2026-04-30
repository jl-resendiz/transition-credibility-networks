"""Alternative specifications for the cascading revelation finding.

The baseline test (learning_order_test.py) finds that fuel
contagion STRENGTHENS with successive retirements — opposite to the
Bayesian learning prediction — but the result is not statistically
significant (FM+NW t=-1.57, p=0.117).

This script implements four alternative specifications to diagnose
whether the weak significance reflects:
  (A) US dominance diluting a strong non-US effect
  (C) The cascade operating specifically in ETS jurisdictions
  (D) Calendar time rather than within-country order as the learning dim
  (E) Non-parametric tests that avoid linearity assumptions

All use Fama-MacBeth event-by-event cross-sectional regressions with
Newey-West HAC standard errors, matching the baseline methodology.

Output: results/metrics/learning_alternatives.md
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


def sig_stars(p):
    if p < 0.01:
        return '***'
    if p < 0.05:
        return '**'
    if p < 0.10:
        return '*'
    return ''


# ── Mann-Whitney U test (stdlib) ────────────────────────────────────

def mann_whitney_u(x, y):
    """Two-sided Mann-Whitney U test.  Returns (U, z, p).

    Normal approximation (valid for n1, n2 >= 8).
    """
    n1, n2 = len(x), len(y)
    if n1 < 2 or n2 < 2:
        return None, None, None
    combined = [(v, 0) for v in x] + [(v, 1) for v in y]
    combined.sort(key=lambda t: t[0])
    # Assign ranks (handle ties with average rank)
    N = n1 + n2
    ranks = [0.0] * N
    i = 0
    while i < N:
        j = i
        while j < N - 1 and combined[j + 1][0] == combined[i][0]:
            j += 1
        avg_rank = (i + j) / 2.0 + 1.0  # 1-based
        for k in range(i, j + 1):
            ranks[k] = avg_rank
        i = j + 1
    R1 = sum(ranks[k] for k in range(N) if combined[k][1] == 0)
    U1 = R1 - n1 * (n1 + 1) / 2.0
    mu = n1 * n2 / 2.0
    # Tie correction
    tie_counts = defaultdict(int)
    for k in range(N):
        tie_counts[ranks[k]] += 1
    tie_correction = sum(t ** 3 - t for t in tie_counts.values()) / (12.0 * (N * (N - 1)))
    sigma = math.sqrt(n1 * n2 * ((N + 1) / 12.0 - tie_correction))
    if sigma < 1e-15:
        return U1, 0.0, 1.0
    z = (U1 - mu) / sigma
    p = p_from_t(z)  # two-sided
    return U1, z, p


# ── Kolmogorov-Smirnov test (stdlib) ────────────────────────────────

def ks_two_sample(x, y):
    """Two-sample KS test.  Returns (D, p_approx).

    Uses the asymptotic formula p ~ 2 * exp(-2 * lambda^2) where
    lambda = sqrt(n_eff) * D and n_eff = n1*n2/(n1+n2).
    """
    n1, n2 = len(x), len(y)
    if n1 < 2 or n2 < 2:
        return None, None
    xs = sorted(x)
    ys = sorted(y)
    # Merge and compute ECDF difference
    all_vals = sorted(set(xs + ys))
    D = 0.0
    for v in all_vals:
        # ECDF of x at v
        f1 = sum(1 for xi in xs if xi <= v) / n1
        f2 = sum(1 for yi in ys if yi <= v) / n2
        D = max(D, abs(f1 - f2))
    n_eff = n1 * n2 / (n1 + n2)
    lam = math.sqrt(n_eff) * D
    # Asymptotic p-value (Smirnov 1948)
    p = 2.0 * math.exp(-2.0 * lam * lam) if lam > 0 else 1.0
    p = min(max(p, 0.0), 1.0)
    return D, p


# ═══════════════════════════════════════════════════════════════════
# LOAD DATA (same as learning_order_test.py)
# ═══════════════════════════════════════════════════════════════════

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


# Load ETS membership
_print('Loading ETS membership...')
firm_ets = {}
ets_path = derived_path('networks', 'firm_ets_membership.csv')
if os.path.exists(ets_path):
    with open(ets_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            firm_ets[row['gvkey']] = {
                'has_ets': int(row.get('has_ets', '0')),
                'ets_names': row.get('ets_names', ''),
                'fic': row.get('fic', ''),
            }
_print(f'  Firms with ETS info: {len(firm_ets)}')
n_ets = sum(1 for v in firm_ets.values() if v['has_ets'] == 1)
_print(f'  Firms in ETS: {n_ets}')


def get_sic4(gvkey):
    f = fundamentals.get(gvkey)
    if f and f.get('sic'):
        return f['sic'][:4]
    return None


def get_country(gvkey):
    f = fundamentals.get(gvkey)
    if f and f.get('fic'):
        return f['fic']
    return None


def get_has_ets(gvkey):
    """Return 1 if firm is in an ETS jurisdiction, 0 otherwise."""
    info = firm_ets.get(gvkey)
    if info:
        return info['has_ets']
    return 0


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
            'country': row.get('country', ''),
        })
_print(f'  First-mover events: {len(all_events)}')

POST_MONTHS = 3
PRE_MONTHS = 24
MIN_OBS_PER_EVENT = 20


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


# ═══════════════════════════════════════════════════════════════════
# STEP 1: Compute event order & metadata
# ═══════════════════════════════════════════════════════════════════

_print('\nComputing event order within each country...')

for ev in all_events:
    country = ev.get('country', '').strip()
    if not country:
        for gk in ev['gvkeys']:
            c = get_country(gk)
            if c:
                country = c
                break
    ev['event_country'] = country

    # Determine if the first-mover firm is in an ETS jurisdiction
    ev_has_ets = 0
    for gk in ev['gvkeys']:
        if get_has_ets(gk):
            ev_has_ets = 1
            break
    ev['has_ets'] = ev_has_ets


def sort_key(ev):
    d = ev.get('event_date', '')
    if d and len(d) >= 7:
        return d[:7]
    y = ev.get('year')
    if y:
        return f'{y}-07'
    return '9999-99'


all_events_sorted = sorted(all_events, key=sort_key)

country_event_count = defaultdict(int)
for ev in all_events_sorted:
    c = ev['event_country']
    ev['event_order'] = country_event_count[c]
    country_event_count[c] += 1
    ev['is_early'] = 1.0 if ev['event_order'] <= 2 else 0.0
    ev['log_order'] = math.log(1.0 + ev['event_order'])

# Determine if event is US
for ev in all_events:
    ev['is_us'] = 1.0 if ev['event_country'] in ('United States', 'USA') else 0.0

# Post-Paris dummy (Paris Agreement entered force 4 Nov 2016)
for ev in all_events:
    y = ev.get('year')
    ev['post_paris'] = 1.0 if (y and y >= 2016) else 0.0

n_us = sum(1 for ev in all_events if ev['is_us'] == 1.0)
n_nonus = sum(1 for ev in all_events if ev['is_us'] == 0.0)
n_ets = sum(1 for ev in all_events if ev['has_ets'] == 1)
n_noets = sum(1 for ev in all_events if ev['has_ets'] == 0)
n_post = sum(1 for ev in all_events if ev['post_paris'] == 1.0)
n_pre = sum(1 for ev in all_events if ev['post_paris'] == 0.0)

_print(f'  US events: {n_us}, Non-US: {n_nonus}')
_print(f'  ETS events: {n_ets}, Non-ETS: {n_noets}')
_print(f'  Post-Paris: {n_post}, Pre-Paris: {n_pre}')

# ═══════════════════════════════════════════════════════════════════
# STEP 2: Build per-event datasets
# ═══════════════════════════════════════════════════════════════════

_print('\nBuilding per-event datasets...')

event_datasets = {}
event_meta = {}

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
        event_meta[event_id] = {
            'event_order': event['event_order'],
            'is_early': event['is_early'],
            'log_order': event['log_order'],
            'country': event['event_country'],
            'plant': event['plant'],
            'year': event['year'],
            'is_us': event['is_us'],
            'has_ets': event['has_ets'],
            'post_paris': event['post_paris'],
        }

n_valid = len(event_datasets)
_print(f'  Valid events (>= {MIN_OBS_PER_EVENT} obs): {n_valid}')
_print(f'  Total obs: {sum(len(v) for v in event_datasets.values())}')


# ═══════════════════════════════════════════════════════════════════
# HELPER: Run FM on a subset and report results
# ═══════════════════════════════════════════════════════════════════

def run_fm_simple(event_ids, label='', vars_list=None):
    """Run Fama-MacBeth on a subset of events.

    Returns dict of {varname: {mean, se, t, p, n}} plus meta keys
    'n_events', 'avg_n', 'avg_r2'.
    """
    if vars_list is None:
        vars_list = ['w_fuel', 'w_geo', 'same_sector']
    betas = defaultdict(list)
    r2s = []
    ns = []
    ids_used = []

    for eid in sorted(event_ids):
        obs = event_datasets[eid]
        ss_vals = set(o['same_sector'] for o in obs)
        use_vars = vars_list if ('same_sector' not in vars_list or len(ss_vals) > 1) else [v for v in vars_list if v != 'same_sector']

        y = [o['car'] for o in obs]
        X = [[1.0] + [o[v] for v in use_vars] for o in obs]
        result = ols_simple(y, X)
        if result is None:
            continue
        names = ['intercept'] + use_vars
        for i, name in enumerate(names):
            betas[name].append(result['beta'][i])
        for v in vars_list:
            if v not in use_vars:
                betas[v].append(float('nan'))
        r2s.append(result['r2'])
        ns.append(result['n'])
        ids_used.append(eid)

    T = len(ids_used)
    out = {
        'n_events': T,
        'avg_n': sum(ns) / T if T > 0 else 0,
        'avg_r2': sum(r2s) / T if T > 0 else 0,
        'event_ids': ids_used,
        'betas': betas,
    }
    for v in ['intercept'] + vars_list:
        clean = [b for b in betas.get(v, []) if not math.isnan(b)]
        if len(clean) < 3:
            out[v] = {'mean': float('nan'), 'se': float('nan'),
                      't': float('nan'), 'p': float('nan'), 'n': 0}
            continue
        mean_b = sum(clean) / len(clean)
        nw_se = newey_west_se(clean)
        t_stat = mean_b / nw_se if nw_se > 1e-15 else 0.0
        p_val = p_from_t(t_stat)
        out[v] = {'mean': mean_b, 'se': nw_se, 't': t_stat, 'p': p_val,
                  'n': len(clean)}
    return out


def run_fm_interaction(event_ids, interaction_var_name, get_interaction_val,
                       base_vars=None):
    """Run FM with an interaction term: w_fuel * interaction_value.

    get_interaction_val(event_id) returns the scalar interaction value
    for each event (e.g., log_order, post_paris, has_ets).
    The interaction is constant within each event, so w_fuel * val
    is perfectly collinear with w_fuel unless val varies across events.

    Returns dict with variable-level FM results.
    """
    if base_vars is None:
        base_vars = ['w_fuel', 'w_geo']
    full_vars = base_vars + [interaction_var_name]
    betas = defaultdict(list)
    r2s = []
    ns = []
    ids_used = []

    for eid in sorted(event_ids):
        obs = event_datasets[eid]
        int_val = get_interaction_val(eid)

        # Build interaction variable for each observation
        for o in obs:
            o[interaction_var_name] = o['w_fuel'] * int_val

        y = [o['car'] for o in obs]
        X = [[1.0] + [o[v] for v in full_vars] for o in obs]
        result = ols_simple(y, X)
        if result is None:
            continue

        names = ['intercept'] + full_vars
        for i, name in enumerate(names):
            betas[name].append(result['beta'][i])
        r2s.append(result['r2'])
        ns.append(result['n'])
        ids_used.append(eid)

    T = len(ids_used)
    out = {
        'n_events': T,
        'avg_n': sum(ns) / T if T > 0 else 0,
        'avg_r2': sum(r2s) / T if T > 0 else 0,
        'event_ids': ids_used,
        'betas': betas,
    }
    for v in ['intercept'] + full_vars:
        clean = [b for b in betas.get(v, []) if not math.isnan(b)]
        if len(clean) < 3:
            out[v] = {'mean': float('nan'), 'se': float('nan'),
                      't': float('nan'), 'p': float('nan'), 'n': 0}
            continue
        mean_b = sum(clean) / len(clean)
        nw_se = newey_west_se(clean)
        t_stat = mean_b / nw_se if nw_se > 1e-15 else 0.0
        p_val = p_from_t(t_stat)
        out[v] = {'mean': mean_b, 'se': nw_se, 't': t_stat, 'p': p_val,
                  'n': len(clean)}
    return out


def fmt_row(v, r):
    """Format a variable result as a markdown table row."""
    s = sig_stars(r['p'])
    return f'| {v:<25} | {r["mean"]:+12.6f} | {r["se"]:10.6f} | {r["t"]:8.3f} | {r["p"]:8.4f}{s:>3s} |'


def table_header():
    return (f'| {"Variable":<25} | {"Mean beta":>12} | {"NW SE":>10} | {"t":>8} | {"p":>8}    |\n'
            f'|{"-"*27}|{"-"*14}|{"-"*12}|{"-"*10}|{"-"*12}|')


# ═══════════════════════════════════════════════════════════════════
# ALT A: US vs non-US split
# ═══════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('ALT A: US vs Non-US Split')
_print('=' * 70)

us_events = [eid for eid in event_datasets if event_meta[eid]['is_us'] == 1.0]
nonus_events = [eid for eid in event_datasets if event_meta[eid]['is_us'] == 0.0]

_print(f'\n  US events: {len(us_events)}, Non-US events: {len(nonus_events)}')

# A1: Simple FM on each subsample
_print('\n  --- A1: Simple fuel beta by subsample ---')
for label, eids in [('US', us_events), ('Non-US', nonus_events)]:
    res = run_fm_simple(eids)
    _print(f'\n  {label}: {res["n_events"]} events, avg N={res["avg_n"]:.0f}')
    for v in ['w_fuel', 'w_geo', 'same_sector']:
        r = res[v]
        if not math.isnan(r['mean']):
            _print(f'    {v:<15} {r["mean"]:+10.4f}  NW SE={r["se"]:.4f}  t={r["t"]:.3f}  p={r["p"]:.4f}{sig_stars(r["p"])}')

# A2: Early vs late WITHIN each subsample
_print('\n  --- A2: Early vs Late within US and Non-US ---')
results_altA = {}
for label, eids in [('US', us_events), ('Non-US', nonus_events)]:
    early_ids = [e for e in eids if event_meta[e]['is_early'] == 1.0]
    late_ids = [e for e in eids if event_meta[e]['is_early'] == 0.0]
    _print(f'\n  {label}: early={len(early_ids)}, late={len(late_ids)}')

    res_early = run_fm_simple(early_ids)
    res_late = run_fm_simple(late_ids)

    fb_early = res_early['w_fuel']['mean'] if not math.isnan(res_early['w_fuel']['mean']) else None
    fb_late = res_late['w_fuel']['mean'] if not math.isnan(res_late['w_fuel']['mean']) else None

    _print(f'    Early fuel beta: {fb_early:+.4f}' if fb_early is not None else '    Early fuel beta: N/A')
    _print(f'    Late  fuel beta: {fb_late:+.4f}' if fb_late is not None else '    Late  fuel beta: N/A')
    if fb_early is not None and fb_late is not None:
        diff = fb_early - fb_late
        _print(f'    Difference (early - late): {diff:+.4f}')

    results_altA[label] = {
        'early': res_early, 'late': res_late,
        'n_early': len(early_ids), 'n_late': len(late_ids),
    }

# A3: Continuous log_order interaction within each subsample
_print('\n  --- A3: w_fuel x log_order within US and Non-US ---')
results_altA_cont = {}
for label, eids in [('US', us_events), ('Non-US', nonus_events)]:
    res = run_fm_interaction(
        eids, 'w_fuel_x_logorder',
        lambda eid: event_meta[eid]['log_order']
    )
    _print(f'\n  {label}: {res["n_events"]} events')
    for v in ['w_fuel', 'w_geo', 'w_fuel_x_logorder']:
        r = res[v]
        if not math.isnan(r['mean']):
            _print(f'    {v:<25} {r["mean"]:+10.4f}  t={r["t"]:.3f}  p={r["p"]:.4f}{sig_stars(r["p"])}')
    results_altA_cont[label] = res


# ═══════════════════════════════════════════════════════════════════
# ALT C: Triple interaction w_fuel x early x has_ets
# ═══════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('ALT C: ETS x Learning Order Interaction')
_print('=' * 70)

ets_events = [eid for eid in event_datasets if event_meta[eid]['has_ets'] == 1]
noets_events = [eid for eid in event_datasets if event_meta[eid]['has_ets'] == 0]
_print(f'\n  ETS events: {len(ets_events)}, Non-ETS events: {len(noets_events)}')

# C1: Early vs late WITHIN ETS and non-ETS
_print('\n  --- C1: Early vs Late within ETS and Non-ETS ---')
results_altC = {}
for label, eids in [('ETS', ets_events), ('Non-ETS', noets_events)]:
    early_ids = [e for e in eids if event_meta[e]['is_early'] == 1.0]
    late_ids = [e for e in eids if event_meta[e]['is_early'] == 0.0]
    _print(f'\n  {label}: early={len(early_ids)}, late={len(late_ids)}')

    res_early = run_fm_simple(early_ids)
    res_late = run_fm_simple(late_ids)

    fb_early = res_early['w_fuel']['mean'] if not math.isnan(res_early['w_fuel']['mean']) else None
    fb_late = res_late['w_fuel']['mean'] if not math.isnan(res_late['w_fuel']['mean']) else None

    _print(f'    Early fuel beta: {fb_early:+.4f}' if fb_early is not None else '    Early fuel beta: N/A')
    _print(f'    Late  fuel beta: {fb_late:+.4f}' if fb_late is not None else '    Late  fuel beta: N/A')
    if fb_early is not None and fb_late is not None:
        diff = fb_early - fb_late
        _print(f'    Difference (early - late): {diff:+.4f}')

    results_altC[label] = {
        'early': res_early, 'late': res_late,
        'n_early': len(early_ids), 'n_late': len(late_ids),
    }

# C2: Continuous log_order interaction within ETS and non-ETS
_print('\n  --- C2: w_fuel x log_order within ETS and Non-ETS ---')
results_altC_cont = {}
for label, eids in [('ETS', ets_events), ('Non-ETS', noets_events)]:
    res = run_fm_interaction(
        eids, 'w_fuel_x_logorder',
        lambda eid: event_meta[eid]['log_order']
    )
    _print(f'\n  {label}: {res["n_events"]} events')
    for v in ['w_fuel', 'w_geo', 'w_fuel_x_logorder']:
        r = res[v]
        if not math.isnan(r['mean']):
            _print(f'    {v:<25} {r["mean"]:+10.4f}  t={r["t"]:.3f}  p={r["p"]:.4f}{sig_stars(r["p"])}')
    results_altC_cont[label] = res

# ── Extract per-event fuel betas (needed by C3, D2, and Alt E) ──
_print('\n  Extracting per-event fuel betas (simple FM)...')
all_fuel_betas = {}  # eid -> fuel_beta
for eid in sorted(event_datasets.keys()):
    obs = event_datasets[eid]
    ss_vals = set(o['same_sector'] for o in obs)
    use_vars = ['w_fuel', 'w_geo', 'same_sector'] if len(ss_vals) > 1 else ['w_fuel', 'w_geo']
    y = [o['car'] for o in obs]
    X = [[1.0] + [o[v] for v in use_vars] for o in obs]
    result = ols_simple(y, X)
    if result is None:
        continue
    all_fuel_betas[eid] = result['beta'][1]
_print(f'  Events with fuel betas: {len(all_fuel_betas)}')

# C3: Second-stage test: regress per-event fuel betas on has_ets
# (Since has_ets is constant within event, we cannot include it as a
# first-stage interaction — it would be collinear with w_fuel.
# Instead we use the two-stage approach: first extract fuel betas,
# then regress them on has_ets in a second-stage cross-section.)
_print('\n  --- C3: Second-stage: fuel_beta on has_ets ---')
all_valid_ids = list(event_datasets.keys())
fb_ets_pairs = []
for eid in sorted(all_fuel_betas.keys()):
    fb_ets_pairs.append((all_fuel_betas[eid], float(event_meta[eid]['has_ets'])))

y_c3 = [p[0] for p in fb_ets_pairs]
X_c3 = [[1.0, p[1]] for p in fb_ets_pairs]
res_c3_ols = ols_simple(y_c3, X_c3)
res_ets_pooled = {}
if res_c3_ols and len(fb_ets_pairs) > 2:
    b_const, b_ets = res_c3_ols['beta']
    n_c3 = len(fb_ets_pairs)
    ss_res = sum(r**2 for r in res_c3_ols['resid'])
    mse = ss_res / (n_c3 - 2)
    XtX = [[sum(X_c3[i][a] * X_c3[i][b] for i in range(n_c3))
            for b in range(2)] for a in range(2)]
    inv_c3 = invert_matrix(XtX)
    if inv_c3:
        se_ets = math.sqrt(mse * inv_c3[1][1])
        t_ets = b_ets / se_ets if se_ets > 1e-15 else 0.0
        p_ets = p_from_t(t_ets)
        _print(f'  N={n_c3}, gamma_1(has_ets) = {b_ets:+.4f}, SE={se_ets:.4f}, t={t_ets:.3f}, p={p_ets:.4f}{sig_stars(p_ets)}')
        _print(f'  R2={res_c3_ols["r2"]:.4f}')
        res_ets_pooled['w_fuel_x_ets'] = {'mean': b_ets, 'se': se_ets,
                                           't': t_ets, 'p': p_ets, 'n': n_c3}
    # Also split by has_ets and compare means with NW
    ets_fb = [fb for fb, e in fb_ets_pairs if e == 1.0]
    noets_fb = [fb for fb, e in fb_ets_pairs if e == 0.0]
    if len(ets_fb) >= 3 and len(noets_fb) >= 3:
        m_ets = sum(ets_fb) / len(ets_fb)
        m_noets = sum(noets_fb) / len(noets_fb)
        nw_ets = newey_west_se(ets_fb)
        nw_noets = newey_west_se(noets_fb)
        diff = m_ets - m_noets
        se_diff = math.sqrt(nw_ets**2 + nw_noets**2)
        t_diff = diff / se_diff if se_diff > 1e-15 else 0.0
        p_diff = p_from_t(t_diff)
        _print(f'  ETS mean fuel beta: {m_ets:+.4f} (N={len(ets_fb)})')
        _print(f'  Non-ETS mean fuel beta: {m_noets:+.4f} (N={len(noets_fb)})')
        _print(f'  Diff (ETS - non-ETS): {diff:+.4f}, NW t={t_diff:.3f}, p={p_diff:.4f}{sig_stars(p_diff)}')
        res_ets_pooled['ets_diff'] = {'mean': diff, 'se': se_diff,
                                       't': t_diff, 'p': p_diff}


# ═══════════════════════════════════════════════════════════════════
# ALT D: Calendar time instead of order
# ═══════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('ALT D: Calendar Time as Learning Dimension')
_print('=' * 70)

# D1: Pre-Paris vs post-Paris subsample
_print('\n  --- D1: Fuel beta pre-Paris vs post-Paris ---')
pre_paris_ids = [eid for eid in event_datasets if event_meta[eid]['post_paris'] == 0.0]
post_paris_ids = [eid for eid in event_datasets if event_meta[eid]['post_paris'] == 1.0]
_print(f'  Pre-Paris: {len(pre_paris_ids)}, Post-Paris: {len(post_paris_ids)}')

results_altD = {}
for label, eids in [('Pre-Paris (before 2016)', pre_paris_ids),
                     ('Post-Paris (2016+)', post_paris_ids)]:
    res = run_fm_simple(eids)
    _print(f'\n  {label}: {res["n_events"]} events, avg N={res["avg_n"]:.0f}')
    for v in ['w_fuel', 'w_geo', 'same_sector']:
        r = res[v]
        if not math.isnan(r['mean']):
            _print(f'    {v:<15} {r["mean"]:+10.4f}  t={r["t"]:.3f}  p={r["p"]:.4f}{sig_stars(r["p"])}')
    results_altD[label] = res

# D2: Second-stage test: regress per-event fuel betas on post_paris
# (Same logic as C3: post_paris is event-level, so second-stage.)
_print('\n  --- D2: Second-stage: fuel_beta on post_paris ---')
fb_paris_pairs = []
for eid in sorted(all_fuel_betas.keys()):
    fb_paris_pairs.append((all_fuel_betas[eid], event_meta[eid]['post_paris']))

y_d2 = [p[0] for p in fb_paris_pairs]
X_d2 = [[1.0, p[1]] for p in fb_paris_pairs]
res_d2_ols = ols_simple(y_d2, X_d2)
res_paris_pooled = {}
if res_d2_ols and len(fb_paris_pairs) > 2:
    b_const_d2, b_paris = res_d2_ols['beta']
    n_d2 = len(fb_paris_pairs)
    ss_res_d2 = sum(r**2 for r in res_d2_ols['resid'])
    mse_d2 = ss_res_d2 / (n_d2 - 2)
    XtX_d2 = [[sum(X_d2[i][a] * X_d2[i][b] for i in range(n_d2))
                for b in range(2)] for a in range(2)]
    inv_d2 = invert_matrix(XtX_d2)
    if inv_d2:
        se_paris = math.sqrt(mse_d2 * inv_d2[1][1])
        t_paris = b_paris / se_paris if se_paris > 1e-15 else 0.0
        p_paris = p_from_t(t_paris)
        _print(f'  N={n_d2}, gamma_1(post_paris) = {b_paris:+.4f}, SE={se_paris:.4f}, t={t_paris:.3f}, p={p_paris:.4f}{sig_stars(p_paris)}')
        _print(f'  R2={res_d2_ols["r2"]:.4f}')
        res_paris_pooled['w_fuel_x_post_paris'] = {'mean': b_paris, 'se': se_paris,
                                                    't': t_paris, 'p': p_paris, 'n': n_d2}
    # Also compare means
    pre_fb_list = [fb for fb, pp in fb_paris_pairs if pp == 0.0]
    post_fb_list = [fb for fb, pp in fb_paris_pairs if pp == 1.0]
    if len(pre_fb_list) >= 3 and len(post_fb_list) >= 3:
        m_pre = sum(pre_fb_list) / len(pre_fb_list)
        m_post = sum(post_fb_list) / len(post_fb_list)
        nw_pre = newey_west_se(pre_fb_list)
        nw_post = newey_west_se(post_fb_list)
        diff_pp = m_post - m_pre
        se_diff_pp = math.sqrt(nw_pre**2 + nw_post**2)
        t_diff_pp = diff_pp / se_diff_pp if se_diff_pp > 1e-15 else 0.0
        p_diff_pp = p_from_t(t_diff_pp)
        _print(f'  Pre-Paris mean fuel beta: {m_pre:+.4f} (N={len(pre_fb_list)})')
        _print(f'  Post-Paris mean fuel beta: {m_post:+.4f} (N={len(post_fb_list)})')
        _print(f'  Diff (Post - Pre): {diff_pp:+.4f}, NW t={t_diff_pp:.3f}, p={p_diff_pp:.4f}{sig_stars(p_diff_pp)}')
        res_paris_pooled['paris_diff'] = {'mean': diff_pp, 'se': se_diff_pp,
                                           't': t_diff_pp, 'p': p_diff_pp}

# D3: w_fuel x event_year (continuous calendar time)
_print('\n  --- D3: w_fuel x event_year (continuous) ---')
# Center year around median to reduce collinearity
years = [event_meta[eid]['year'] for eid in event_datasets
         if event_meta[eid]['year'] is not None]
years.sort()
median_year = years[len(years) // 2] if years else 2018
_print(f'  Median event year: {median_year}')

res_year_cont = run_fm_interaction(
    [eid for eid in all_valid_ids if event_meta[eid]['year'] is not None],
    'w_fuel_x_year',
    lambda eid: float(event_meta[eid]['year'] - median_year)
)
_print(f'  Events: {res_year_cont["n_events"]}')
for v in ['w_fuel', 'w_geo', 'w_fuel_x_year']:
    r = res_year_cont[v]
    if not math.isnan(r['mean']):
        _print(f'    {v:<25} {r["mean"]:+10.4f}  t={r["t"]:.3f}  p={r["p"]:.4f}{sig_stars(r["p"])}')

# D4: Tercile analysis by event year
_print('\n  --- D4: Fuel beta by year tercile ---')
year_ids = [(event_meta[eid]['year'], eid) for eid in event_datasets
            if event_meta[eid]['year'] is not None]
year_ids.sort()
n_total = len(year_ids)
t1_end = n_total // 3
t2_end = 2 * n_total // 3
terciles = [
    ('T1 (earliest)', [eid for _, eid in year_ids[:t1_end]]),
    ('T2 (middle)', [eid for _, eid in year_ids[t1_end:t2_end]]),
    ('T3 (latest)', [eid for _, eid in year_ids[t2_end:]]),
]
results_altD_tercile = {}
for label, eids in terciles:
    if not eids:
        continue
    yrs = [event_meta[e]['year'] for e in eids]
    res = run_fm_simple(eids)
    _print(f'\n  {label} (years {min(yrs)}-{max(yrs)}): {res["n_events"]} events')
    r = res['w_fuel']
    if not math.isnan(r['mean']):
        _print(f'    w_fuel: {r["mean"]:+.4f}  t={r["t"]:.3f}  p={r["p"]:.4f}{sig_stars(r["p"])}')
    results_altD_tercile[label] = res


# ═══════════════════════════════════════════════════════════════════
# ALT E: Non-parametric comparison of fuel beta distributions
# ═══════════════════════════════════════════════════════════════════

_print('\n' + '=' * 70)
_print('ALT E: Non-Parametric Tests on Fuel Beta Distributions')
_print('=' * 70)

# Per-event fuel betas already computed above (before C3)
_print(f'\n  Using pre-computed per-event fuel betas: {len(all_fuel_betas)}')

# E1: Early vs Late (Mann-Whitney + KS)
early_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
            if event_meta[eid]['is_early'] == 1.0]
late_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
           if event_meta[eid]['is_early'] == 0.0]

_print(f'\n  --- E1: Early ({len(early_fb)}) vs Late ({len(late_fb)}) fuel betas ---')
_print(f'  Early: mean={sum(early_fb)/len(early_fb):+.4f}, median={sorted(early_fb)[len(early_fb)//2]:+.4f}')
_print(f'  Late:  mean={sum(late_fb)/len(late_fb):+.4f}, median={sorted(late_fb)[len(late_fb)//2]:+.4f}')

U, z_mw, p_mw = mann_whitney_u(early_fb, late_fb)
_print(f'  Mann-Whitney U={U:.0f}, z={z_mw:.3f}, p={p_mw:.4f}{sig_stars(p_mw)}')

D_ks, p_ks = ks_two_sample(early_fb, late_fb)
_print(f'  Kolmogorov-Smirnov D={D_ks:.4f}, p={p_ks:.4f}{sig_stars(p_ks)}')

# Two-sample t-test (Welch) for early vs late fuel betas
_mean_e = sum(early_fb) / len(early_fb)
_mean_l = sum(late_fb) / len(late_fb)
_var_e = sum((x - _mean_e)**2 for x in early_fb) / (len(early_fb) - 1) if len(early_fb) > 1 else 0
_var_l = sum((x - _mean_l)**2 for x in late_fb) / (len(late_fb) - 1) if len(late_fb) > 1 else 0
_se_diff = (_var_e / len(early_fb) + _var_l / len(late_fb)) ** 0.5
_t_welch = (_mean_e - _mean_l) / _se_diff if _se_diff > 1e-15 else 0.0
_p_welch = p_from_t(_t_welch)
_print(f'  Welch t-test: diff={_mean_e - _mean_l:+.4f}, t={_t_welch:.3f}, p={_p_welch:.4f}{sig_stars(_p_welch)}')

results_altE = {
    'early_vs_late': {
        'n_early': len(early_fb), 'n_late': len(late_fb),
        'mean_early': sum(early_fb) / len(early_fb),
        'mean_late': sum(late_fb) / len(late_fb),
        'median_early': sorted(early_fb)[len(early_fb) // 2],
        'median_late': sorted(late_fb)[len(late_fb) // 2],
        'mw_U': U, 'mw_z': z_mw, 'mw_p': p_mw,
        'ks_D': D_ks, 'ks_p': p_ks,
        'welch_t': _t_welch, 'welch_p': _p_welch,
        'welch_diff': _mean_e - _mean_l,
    }
}

# E2: US early vs late (non-parametric)
us_early_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
               if event_meta[eid]['is_us'] == 1.0 and event_meta[eid]['is_early'] == 1.0]
us_late_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
              if event_meta[eid]['is_us'] == 1.0 and event_meta[eid]['is_early'] == 0.0]
nonus_early_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
                  if event_meta[eid]['is_us'] == 0.0 and event_meta[eid]['is_early'] == 1.0]
nonus_late_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
                 if event_meta[eid]['is_us'] == 0.0 and event_meta[eid]['is_early'] == 0.0]

_print(f'\n  --- E2: US early ({len(us_early_fb)}) vs late ({len(us_late_fb)}) ---')
if len(us_early_fb) >= 2 and len(us_late_fb) >= 2:
    U_us, z_us, p_us = mann_whitney_u(us_early_fb, us_late_fb)
    _print(f'  Mann-Whitney z={z_us:.3f}, p={p_us:.4f}{sig_stars(p_us)}')
    results_altE['us_early_vs_late'] = {'mw_z': z_us, 'mw_p': p_us,
                                         'n_early': len(us_early_fb), 'n_late': len(us_late_fb)}
else:
    _print('  Insufficient observations')

_print(f'\n  --- E3: Non-US early ({len(nonus_early_fb)}) vs late ({len(nonus_late_fb)}) ---')
if len(nonus_early_fb) >= 2 and len(nonus_late_fb) >= 2:
    U_nonus, z_nonus, p_nonus = mann_whitney_u(nonus_early_fb, nonus_late_fb)
    _print(f'  Mann-Whitney z={z_nonus:.3f}, p={p_nonus:.4f}{sig_stars(p_nonus)}')
    results_altE['nonus_early_vs_late'] = {'mw_z': z_nonus, 'mw_p': p_nonus,
                                            'n_early': len(nonus_early_fb), 'n_late': len(nonus_late_fb)}
else:
    _print('  Insufficient observations')

# E4: ETS early vs late (non-parametric)
ets_early_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
                if event_meta[eid]['has_ets'] == 1 and event_meta[eid]['is_early'] == 1.0]
ets_late_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
               if event_meta[eid]['has_ets'] == 1 and event_meta[eid]['is_early'] == 0.0]

_print(f'\n  --- E4: ETS early ({len(ets_early_fb)}) vs late ({len(ets_late_fb)}) ---')
if len(ets_early_fb) >= 2 and len(ets_late_fb) >= 2:
    U_ets, z_ets, p_ets = mann_whitney_u(ets_early_fb, ets_late_fb)
    _print(f'  Mann-Whitney z={z_ets:.3f}, p={p_ets:.4f}{sig_stars(p_ets)}')
    results_altE['ets_early_vs_late'] = {'mw_z': z_ets, 'mw_p': p_ets,
                                          'n_early': len(ets_early_fb), 'n_late': len(ets_late_fb)}
else:
    _print('  Insufficient observations')

# E5: Pre-Paris vs Post-Paris fuel betas (non-parametric)
pre_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
          if event_meta[eid]['post_paris'] == 0.0]
post_fb = [all_fuel_betas[eid] for eid in all_fuel_betas
           if event_meta[eid]['post_paris'] == 1.0]

_print(f'\n  --- E5: Pre-Paris ({len(pre_fb)}) vs Post-Paris ({len(post_fb)}) fuel betas ---')
if len(pre_fb) >= 2 and len(post_fb) >= 2:
    _print(f'  Pre:  mean={sum(pre_fb)/len(pre_fb):+.4f}, median={sorted(pre_fb)[len(pre_fb)//2]:+.4f}')
    _print(f'  Post: mean={sum(post_fb)/len(post_fb):+.4f}, median={sorted(post_fb)[len(post_fb)//2]:+.4f}')
    U_pp, z_pp, p_pp = mann_whitney_u(pre_fb, post_fb)
    D_pp, pk_pp = ks_two_sample(pre_fb, post_fb)
    _print(f'  Mann-Whitney z={z_pp:.3f}, p={p_pp:.4f}{sig_stars(p_pp)}')
    _print(f'  KS D={D_pp:.4f}, p={pk_pp:.4f}{sig_stars(pk_pp)}')
    results_altE['pre_vs_post_paris'] = {
        'mw_z': z_pp, 'mw_p': p_pp, 'ks_D': D_pp, 'ks_p': pk_pp,
        'n_pre': len(pre_fb), 'n_post': len(post_fb),
        'mean_pre': sum(pre_fb) / len(pre_fb),
        'mean_post': sum(post_fb) / len(post_fb),
    }
else:
    _print('  Insufficient observations')


# ═══════════════════════════════════════════════════════════════════
# WRITE RESULTS
# ═══════════════════════════════════════════════════════════════════

_print('\nWriting results...')

out_path = results_path('metrics', 'learning_alternatives.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

L = []  # output lines


def add(line=''):
    L.append(line)


add('# Learning Order / Cascading Revelation: Alternative Specifications')
add()
add('**Motivation**: The baseline test finds fuel contagion strengthens with')
add('successive retirements (opposite to Bayesian learning), but the result')
add('is not statistically significant (FM+NW t=-1.57, p=0.117). These')
add('alternatives diagnose whether significance improves under different')
add('sample splits, functional forms, and non-parametric tests.')
add()

# ── Alt A ──

add('## Alt A: US vs Non-US Split')
add()
add('**Rationale**: 93 of 179 events are US. US retirements span order 0-92,')
add('so log_order is dominated by within-US variation. The cascade might')
add('be stronger outside the US, where each retirement is genuinely novel.')
add()

add('### A1: Simple fuel beta by subsample')
add()
add(table_header())
for label, eids in [('US', us_events), ('Non-US', nonus_events)]:
    res = run_fm_simple(eids)
    r = res['w_fuel']
    if not math.isnan(r['mean']):
        add(fmt_row(f'w_fuel ({label}, N={res["n_events"]})', r))
add()

add('### A2: Early vs Late within each subsample')
add()
for label in ['US', 'Non-US']:
    ra = results_altA[label]
    add(f'**{label}**: early={ra["n_early"]}, late={ra["n_late"]}')
    add()
    add(table_header())
    re = ra['early']['w_fuel']
    rl = ra['late']['w_fuel']
    if not math.isnan(re['mean']):
        add(fmt_row(f'w_fuel (early)', re))
    if not math.isnan(rl['mean']):
        add(fmt_row(f'w_fuel (late)', rl))
    add()
    if not math.isnan(re['mean']) and not math.isnan(rl['mean']):
        diff = re['mean'] - rl['mean']
        add(f'Difference (early - late): {diff:+.4f}')
        add()

add('### A3: Continuous w_fuel x log_order within each subsample')
add()
for label in ['US', 'Non-US']:
    res = results_altA_cont[label]
    add(f'**{label}**: {res["n_events"]} events')
    add()
    add(table_header())
    for v in ['w_fuel', 'w_fuel_x_logorder']:
        r = res[v]
        if not math.isnan(r['mean']):
            add(fmt_row(v, r))
    add()

# ── Alt C ──

add('## Alt C: ETS x Learning Order Interaction')
add()
add('**Rationale**: In ETS jurisdictions, retirements reinforce regime')
add('credibility. In non-ETS jurisdictions, retirements may be idiosyncratic.')
add('The cascade should operate specifically in ETS jurisdictions.')
add()

add(f'ETS events: {len(ets_events)}, Non-ETS events: {len(noets_events)}')
add()

add('### C1: Early vs Late within ETS and Non-ETS')
add()
for label in ['ETS', 'Non-ETS']:
    ra = results_altC[label]
    add(f'**{label}**: early={ra["n_early"]}, late={ra["n_late"]}')
    add()
    add(table_header())
    re = ra['early']['w_fuel']
    rl = ra['late']['w_fuel']
    if not math.isnan(re['mean']):
        add(fmt_row(f'w_fuel (early)', re))
    if not math.isnan(rl['mean']):
        add(fmt_row(f'w_fuel (late)', rl))
    add()
    if not math.isnan(re['mean']) and not math.isnan(rl['mean']):
        diff = re['mean'] - rl['mean']
        add(f'Difference (early - late): {diff:+.4f}')
        add()

add('### C2: Continuous w_fuel x log_order within ETS and Non-ETS')
add()
for label in ['ETS', 'Non-ETS']:
    res = results_altC_cont[label]
    add(f'**{label}**: {res["n_events"]} events')
    add()
    add(table_header())
    for v in ['w_fuel', 'w_fuel_x_logorder']:
        r = res[v]
        if not math.isnan(r['mean']):
            add(fmt_row(v, r))
    add()

add('### C3: Second-stage fuel_beta on has_ets')
add()
add('Since has_ets is constant within each event, the interaction w_fuel x has_ets')
add('is collinear with w_fuel in the first-stage cross-section. Instead, we extract')
add('per-event fuel betas and regress them on has_ets in a second stage.')
add()
if 'w_fuel_x_ets' in res_ets_pooled:
    r = res_ets_pooled['w_fuel_x_ets']
    add(f'gamma_1(has_ets) = {r["mean"]:+.4f}, SE = {r["se"]:.4f}, t = {r["t"]:.3f}, p = {r["p"]:.4f}{sig_stars(r["p"])}')
if 'ets_diff' in res_ets_pooled:
    r = res_ets_pooled['ets_diff']
    add(f'ETS - Non-ETS mean fuel beta: {r["mean"]:+.4f}, NW t = {r["t"]:.3f}, p = {r["p"]:.4f}{sig_stars(r["p"])}')
add()

# ── Alt D ──

add('## Alt D: Calendar Time as Learning Dimension')
add()
add('**Rationale**: Instead of within-country order, the learning dimension')
add('might be calendar time. The Paris Agreement (2015) and subsequent COP')
add('commitments may have shifted the baseline belief about transition')
add('probability, making each post-Paris retirement more informative.')
add()

add('### D1: Pre-Paris vs Post-Paris subsample')
add()
for label in ['Pre-Paris (before 2016)', 'Post-Paris (2016+)']:
    res = results_altD[label]
    add(f'**{label}**: {res["n_events"]} events')
    add()
    add(table_header())
    for v in ['w_fuel', 'w_geo', 'same_sector']:
        r = res[v]
        if not math.isnan(r['mean']):
            add(fmt_row(v, r))
    add()

add('### D2: Second-stage fuel_beta on post_paris')
add()
add('Same second-stage approach as C3: post_paris is event-level.')
add()
if 'w_fuel_x_post_paris' in res_paris_pooled:
    r = res_paris_pooled['w_fuel_x_post_paris']
    add(f'gamma_1(post_paris) = {r["mean"]:+.4f}, SE = {r["se"]:.4f}, t = {r["t"]:.3f}, p = {r["p"]:.4f}{sig_stars(r["p"])}')
if 'paris_diff' in res_paris_pooled:
    r = res_paris_pooled['paris_diff']
    add(f'Post - Pre Paris mean fuel beta: {r["mean"]:+.4f}, NW t = {r["t"]:.3f}, p = {r["p"]:.4f}{sig_stars(r["p"])}')
add()

add('### D3: Continuous w_fuel x event_year')
add()
add(f'Year centered at median = {median_year}')
add()
add(table_header())
for v in ['w_fuel', 'w_fuel_x_year']:
    r = res_year_cont[v]
    if not math.isnan(r['mean']):
        add(fmt_row(v, r))
add()

add('### D4: Fuel beta by year tercile')
add()
add(table_header())
for label, eids in terciles:
    if not eids:
        continue
    res = results_altD_tercile.get(label)
    if res:
        r = res['w_fuel']
        yrs = [event_meta[e]['year'] for e in eids]
        if not math.isnan(r['mean']):
            add(fmt_row(f'w_fuel {label} ({min(yrs)}-{max(yrs)}, N={res["n_events"]})', r))
add()

# ── Alt E ──

add('## Alt E: Non-Parametric Tests')
add()
add('**Rationale**: Regression-based tests assume linearity. Non-parametric')
add('tests compare the DISTRIBUTION of event-level fuel betas without')
add('functional form assumptions.')
add()

r = results_altE['early_vs_late']
add('### E1: Early vs Late (full sample)')
add()
add(f'| Statistic | Early (N={r["n_early"]}) | Late (N={r["n_late"]}) |')
add('|---|---|---|')
add(f'| Mean fuel beta | {r["mean_early"]:+.4f} | {r["mean_late"]:+.4f} |')
add(f'| Median fuel beta | {r["median_early"]:+.4f} | {r["median_late"]:+.4f} |')
add()
add(f'Welch t-test: diff={r["welch_diff"]:+.4f}, t={r["welch_t"]:.3f}, p={r["welch_p"]:.4f}{sig_stars(r["welch_p"])}')
add(f'Mann-Whitney U: z={r["mw_z"]:.3f}, p={r["mw_p"]:.4f}{sig_stars(r["mw_p"])}')
add(f'Kolmogorov-Smirnov: D={r["ks_D"]:.4f}, p={r["ks_p"]:.4f}{sig_stars(r["ks_p"])}')
add()

if 'us_early_vs_late' in results_altE:
    r = results_altE['us_early_vs_late']
    add(f'### E2: US early ({r["n_early"]}) vs late ({r["n_late"]})')
    add(f'Mann-Whitney z={r["mw_z"]:.3f}, p={r["mw_p"]:.4f}{sig_stars(r["mw_p"])}')
    add()

if 'nonus_early_vs_late' in results_altE:
    r = results_altE['nonus_early_vs_late']
    add(f'### E3: Non-US early ({r["n_early"]}) vs late ({r["n_late"]})')
    add(f'Mann-Whitney z={r["mw_z"]:.3f}, p={r["mw_p"]:.4f}{sig_stars(r["mw_p"])}')
    add()

if 'ets_early_vs_late' in results_altE:
    r = results_altE['ets_early_vs_late']
    add(f'### E4: ETS early ({r["n_early"]}) vs late ({r["n_late"]})')
    add(f'Mann-Whitney z={r["mw_z"]:.3f}, p={r["mw_p"]:.4f}{sig_stars(r["mw_p"])}')
    add()

if 'pre_vs_post_paris' in results_altE:
    r = results_altE['pre_vs_post_paris']
    add(f'### E5: Pre-Paris ({r["n_pre"]}) vs Post-Paris ({r["n_post"]}) fuel betas')
    add(f'Pre-Paris mean: {r["mean_pre"]:+.4f}, Post-Paris mean: {r["mean_post"]:+.4f}')
    add(f'Mann-Whitney z={r["mw_z"]:.3f}, p={r["mw_p"]:.4f}{sig_stars(r["mw_p"])}')
    add(f'KS D={r["ks_D"]:.4f}, p={r["ks_p"]:.4f}{sig_stars(r["ks_p"])}')
    add()

# ── Summary ──

add('## Summary: Which Alternatives Strengthen the Finding?')
add()
add('| Alternative | Key test | t-stat | p-value | Verdict |')
add('|---|---|---|---|---|')

# Baseline
add(f'| Baseline (log_order) | w_fuel x log_order | -1.569 | 0.1167 | Marginal |')

# Alt A
for label in ['US', 'Non-US']:
    res = results_altA_cont[label]
    r = res.get('w_fuel_x_logorder', {})
    if not math.isnan(r.get('mean', float('nan'))):
        verdict = 'Significant' if r['p'] < 0.10 else 'Not significant'
        add(f'| Alt A: {label} log_order | w_fuel x log_order | {r["t"]:.3f} | {r["p"]:.4f} | {verdict} |')

# Alt C
for label in ['ETS', 'Non-ETS']:
    res = results_altC_cont[label]
    r = res.get('w_fuel_x_logorder', {})
    if not math.isnan(r.get('mean', float('nan'))):
        verdict = 'Significant' if r['p'] < 0.10 else 'Not significant'
        add(f'| Alt C: {label} log_order | w_fuel x log_order | {r["t"]:.3f} | {r["p"]:.4f} | {verdict} |')

# Alt C3
r = res_ets_pooled.get('ets_diff', res_ets_pooled.get('w_fuel_x_ets', {}))
if r and not math.isnan(r.get('mean', float('nan'))):
    verdict = 'Significant' if r['p'] < 0.10 else 'Not significant'
    add(f'| Alt C: ETS vs non-ETS | fuel_beta ~ has_ets | {r["t"]:.3f} | {r["p"]:.4f} | {verdict} |')

# Alt D
r = res_paris_pooled.get('paris_diff', res_paris_pooled.get('w_fuel_x_post_paris', {}))
if r and not math.isnan(r.get('mean', float('nan'))):
    verdict = 'Significant' if r['p'] < 0.10 else 'Not significant'
    add(f'| Alt D: post-Paris | fuel_beta ~ post_paris | {r["t"]:.3f} | {r["p"]:.4f} | {verdict} |')

r = res_year_cont.get('w_fuel_x_year', {})
if not math.isnan(r.get('mean', float('nan'))):
    verdict = 'Significant' if r['p'] < 0.10 else 'Not significant'
    add(f'| Alt D: calendar year | w_fuel x year | {r["t"]:.3f} | {r["p"]:.4f} | {verdict} |')

# Alt E
r = results_altE.get('early_vs_late', {})
if r.get('mw_p') is not None:
    verdict = 'Significant' if r['mw_p'] < 0.10 else 'Not significant'
    add(f'| Alt E: Mann-Whitney | early vs late | {r["mw_z"]:.3f} | {r["mw_p"]:.4f} | {verdict} |')

add()

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(L))

_print(f'\nResults written to {out_path}')
_print('Done.')
