"""Firm-level collapsed specification and corrected Wald F-test.

Addresses two issues identified in the econometric review:

1. The pooled R^2 = 0.0012 is mechanically depressed because the same firm
   appears ~179 times with nearly identical spatial weights but different CARs.
   This script collapses to firm-level by averaging CARs and weights across
   events, then runs the cross-sectional regression on N ~ 389 firms.

2. The original F-test in strategy2_joint_tests.py uses the SSR-based
   (homoskedastic) F-statistic, which is invalid under clustered errors.
   This script computes the correct Wald F-statistic using the cluster-robust
   variance-covariance matrix: F = (R*beta)' [R V R']^{-1} (R*beta) / q.

Output: results/metrics/strategy2_firm_level_test.md
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

POST_MONTHS = 3
PRE_MONTHS = 24
CHANNEL_VARS = ['w_geo', 'w_fuel', 'w_reg']
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']


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


# ── OLS with optional clustering and full vcov ────────────────────

def _cluster_cov(X, resid, cluster_map, k):
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
            for b in range(a, k):
                v = xua * xu[b]
                S[a][b] += v
                if a != b:
                    S[b][a] += v
    return S


def ols_full(data, y_var, x_vars, cluster_var=None):
    """OLS with HC1 or cluster-robust SEs. Returns full vcov matrix."""
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

    XtX = [[sum(X[i][a] * X[i][b] for i in range(n))
            for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]

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
                for b in range(k):
                    V[a][b] *= scale
        se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
    else:
        # HC1 (White) robust standard errors for cross-sectional case
        S = [[0.0] * k for _ in range(k)]
        for i in range(n):
            ri = resid[i]
            for a in range(k):
                for b in range(a, k):
                    v = X[i][a] * X[i][b] * ri * ri
                    S[a][b] += v
                    if a != b:
                        S[b][a] += v
        V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
        hc1_scale = n / (n - k)
        for a in range(k):
            for b in range(k):
                V[a][b] *= hc1_scale
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
        'clusters': G,
    }


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


def wald_f_test(res, restriction_vars):
    """Wald F-test using cluster-robust vcov: F = (Rb)' [RVR']^{-1} (Rb) / q.

    This is the correct F-test under heteroskedasticity/clustering,
    unlike the SSR-based F-test which assumes homoskedastic errors.
    """
    names = ['intercept'] + [v for v in res['beta'] if v != 'intercept']
    q = len(restriction_vars)
    V = res['V']

    # Build R matrix (q x k): R selects the restricted coefficients
    # Under H0: R*beta = 0
    idx_map = {name: i for i, name in enumerate(names)}
    Rb = [res['beta'][v] for v in restriction_vars]

    # Extract RVR' (q x q submatrix of V)
    RVR = [[V[idx_map[restriction_vars[a]]][idx_map[restriction_vars[b]]]
            for b in range(q)] for a in range(q)]

    inv_RVR = invert_matrix(RVR)
    if inv_RVR is None:
        return None, None

    # F = (Rb)' inv(RVR') (Rb) / q
    quad = sum(Rb[a] * sum(inv_RVR[a][b] * Rb[b] for b in range(q))
               for a in range(q))
    f_stat = quad / q
    return f_stat, q


# ── Load data (same as strategy2_joint_tests.py) ────────────────────

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
    _print(f'  W_reg edges: {sum(len(v) for v in W_reg.values())}')

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


# ── CAR computation ──────────────────────────────────────────────────

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


# ── Build pooled dataset ─────────────────────────────────────────────

def build_obs():
    obs = []
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
                    'event_id': event_id,
                    'gvkey': gk,
                })
    return obs


# ── Collapse to firm level ───────────────────────────────────────────

def collapse_to_firm(obs):
    """Average CARs and weights across events for each firm."""
    firm_data = defaultdict(lambda: {
        'car_sum': 0.0, 'w_geo_sum': 0.0, 'w_fuel_sum': 0.0,
        'w_reg_sum': 0.0, 'same_sector_sum': 0.0, 'count': 0
    })
    for o in obs:
        gk = o['gvkey']
        firm_data[gk]['car_sum'] += o['car']
        firm_data[gk]['w_geo_sum'] += o['w_geo']
        firm_data[gk]['w_fuel_sum'] += o['w_fuel']
        firm_data[gk]['w_reg_sum'] += o['w_reg']
        firm_data[gk]['same_sector_sum'] += o['same_sector']
        firm_data[gk]['count'] += 1

    collapsed = []
    for gk, d in firm_data.items():
        n = d['count']
        collapsed.append({
            'car': d['car_sum'] / n,
            'w_geo': d['w_geo_sum'] / n,
            'w_fuel': d['w_fuel_sum'] / n,
            'w_reg': d['w_reg_sum'] / n,
            'same_sector': d['same_sector_sum'] / n,
            'gvkey': gk,
            'n_events': n,
        })
    return collapsed


# ══════════════════════════════════════════════════════════════════════
# MAIN
# ══════════════════════════════════════════════════════════════════════

_print()
_print('=' * 70)
_print('FIRM-LEVEL COLLAPSED SPECIFICATION + CORRECTED WALD F-TEST')
_print('=' * 70)

# Build pooled dataset
_print('\nBuilding pooled dataset...')
obs = build_obs()
n_pooled = len(obs)
n_events = len(set(o['event_id'] for o in obs))
n_firms_pooled = len(set(o['gvkey'] for o in obs))
_print(f'  N_pooled = {n_pooled} (event-firm pairs)')
_print(f'  Events: {n_events}')
_print(f'  Unique firms: {n_firms_pooled}')
_print(f'  Avg appearances per firm: {n_pooled / n_firms_pooled:.1f}')

# Check same_sector variation
ss_vals = set(o['same_sector'] for o in obs)
spec_vars = SPEC_VARS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']

# ── Panel A: Pooled specification (original, event-clustered) ────────

_print('\n' + '-' * 60)
_print('PANEL A: POOLED SPECIFICATION (event-clustered)')
_print('-' * 60)

res_pooled = ols_full(obs, 'car', spec_vars, cluster_var='event_id')
if res_pooled:
    _print(f'  N = {res_pooled["n"]}, R2 = {res_pooled["r2"]:.6f}, '
           f'Clusters = {res_pooled["clusters"]}')
    for v in ['intercept'] + spec_vars:
        b = res_pooled['beta'][v]
        s = res_pooled['se'][v]
        t = res_pooled['t'][v]
        p = p_from_t(t)
        stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
        _print(f'  {v:<15} {b:+.6f}  ({s:.6f})  t={t:.3f}  p={p:.4f}{stars}')

    # Wald F-test (cluster-robust)
    f_wald, q = wald_f_test(res_pooled, CHANNEL_VARS)
    if f_wald is not None:
        # Approximate p-value using F(q, G-1) via normal for large G
        # For G=175, normal is fine
        chi2 = f_wald * q
        # Chi-squared p-value approximation via normal
        z_chi2 = (chi2 / q - 1) * math.sqrt(q / 2.0)  # Wilson-Hilferty
        p_wald = 2 * (1 - _normal_cdf(abs(z_chi2)))
        # Better: use permutation for exact p-value
        _print(f'\n  Wald F-test (cluster-robust): F = {f_wald:.4f}, q = {q}')
        _print(f'  Chi2 = {chi2:.4f}')

# ── Panel B: Pooled specification (two-way clustered) ────────────────

_print('\n' + '-' * 60)
_print('PANEL B: POOLED SPECIFICATION (two-way clustered: event + firm)')
_print('-' * 60)

# Implement two-way clustering
n = len(obs)
k = len(spec_vars) + 1
y = [o['car'] for o in obs]
y_mean = sum(y) / n
ss_tot = sum((yi - y_mean) ** 2 for yi in y)
X = [[1.0] + [o[xv] for xv in spec_vars] for o in obs]
XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]
inv_XtX = invert_matrix(XtX)
beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]
y_hat = [sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
resid = [y[i] - y_hat[i] for i in range(n)]
ss_res = sum(r ** 2 for r in resid)
r2_tw = 1 - ss_res / ss_tot

# Cluster maps
clust_event = defaultdict(list)
clust_firm = defaultdict(list)
clust_both = defaultdict(list)
for i, o in enumerate(obs):
    clust_event[o['event_id']].append(i)
    clust_firm[o['gvkey']].append(i)
    clust_both[(o['event_id'], o['gvkey'])].append(i)

S_event = _cluster_cov(X, resid, clust_event, k)
S_firm = _cluster_cov(X, resid, clust_firm, k)
S_both = _cluster_cov(X, resid, clust_both, k)
S_tw = [[S_event[a][b] + S_firm[a][b] - S_both[a][b]
         for b in range(k)] for a in range(k)]
V_tw = mat_mul(mat_mul(inv_XtX, S_tw), inv_XtX)
G_min = min(len(clust_event), len(clust_firm))
if G_min > 1:
    scale = (G_min / (G_min - 1)) * ((n - 1) / (n - k))
    for a in range(k):
        for b in range(k):
            V_tw[a][b] *= scale
se_tw = [math.sqrt(V_tw[a][a]) if V_tw[a][a] > 0 else 0.0 for a in range(k)]
t_tw = [beta[a] / se_tw[a] if se_tw[a] > 1e-15 else 0.0 for a in range(k)]

names_tw = ['intercept'] + spec_vars
_print(f'  N = {n}, R2 = {r2_tw:.6f}')
_print(f'  Event clusters: {len(clust_event)}, Firm clusters: {len(clust_firm)}')
for idx, v in enumerate(names_tw):
    p = p_from_t(t_tw[idx])
    stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
    _print(f'  {v:<15} {beta[idx]:+.6f}  ({se_tw[idx]:.6f})  t={t_tw[idx]:.3f}  p={p:.4f}{stars}')

# Wald F-test (two-way clustered)
idx_map_tw = {name: i for i, name in enumerate(names_tw)}
Rb_tw = [beta[idx_map_tw[v]] for v in CHANNEL_VARS]
q_tw = len(CHANNEL_VARS)
RVR_tw = [[V_tw[idx_map_tw[CHANNEL_VARS[a]]][idx_map_tw[CHANNEL_VARS[b]]]
           for b in range(q_tw)] for a in range(q_tw)]
inv_RVR_tw = invert_matrix(RVR_tw)
if inv_RVR_tw:
    quad_tw = sum(Rb_tw[a] * sum(inv_RVR_tw[a][b] * Rb_tw[b] for b in range(q_tw))
                  for a in range(q_tw))
    f_wald_tw = quad_tw / q_tw
    _print(f'\n  Wald F-test (two-way clustered): F = {f_wald_tw:.4f}, q = {q_tw}')

# Difference test with two-way SEs
idx_geo = idx_map_tw['w_geo']
idx_fuel = idx_map_tw['w_fuel']
diff_tw = beta[idx_geo] - beta[idx_fuel]
se_diff_tw = math.sqrt(V_tw[idx_geo][idx_geo] + V_tw[idx_fuel][idx_fuel]
                        - 2 * V_tw[idx_geo][idx_fuel])
t_diff_tw = diff_tw / se_diff_tw if se_diff_tw > 1e-15 else 0.0
p_diff_tw = p_from_t(t_diff_tw)
_print(f'\n  Difference test (beta_geo - beta_fuel): {diff_tw:+.6f}')
_print(f'  SE = {se_diff_tw:.6f}, t = {t_diff_tw:.3f}, p = {p_diff_tw:.4f}')

# ── Panel C: Firm-level collapsed specification ──────────────────────

_print('\n' + '-' * 60)
_print('PANEL C: FIRM-LEVEL COLLAPSED SPECIFICATION (HC1 robust SEs)')
_print('-' * 60)

collapsed = collapse_to_firm(obs)
_print(f'  N_firms = {len(collapsed)}')
_print(f'  Avg events per firm: {sum(c["n_events"] for c in collapsed) / len(collapsed):.1f}')
_print(f'  Min events: {min(c["n_events"] for c in collapsed)}')
_print(f'  Max events: {max(c["n_events"] for c in collapsed)}')

# Descriptive stats for collapsed data
for v in ['car', 'w_geo', 'w_fuel', 'w_reg']:
    vals = [c[v] for c in collapsed]
    mean_v = sum(vals) / len(vals)
    sd_v = math.sqrt(sum((x - mean_v) ** 2 for x in vals) / (len(vals) - 1))
    _print(f'  {v}: mean={mean_v:.4f}, sd={sd_v:.4f}')

# Check same_sector variation in collapsed
ss_collapsed = set(round(c['same_sector'], 4) for c in collapsed)
spec_vars_c = spec_vars if len(ss_collapsed) > 1 else ['w_geo', 'w_fuel', 'w_reg']

res_firm = ols_full(collapsed, 'car', spec_vars_c)
if res_firm:
    _print(f'\n  N = {res_firm["n"]}, R2 = {res_firm["r2"]:.6f}')
    for v in ['intercept'] + spec_vars_c:
        b = res_firm['beta'][v]
        s = res_firm['se'][v]
        t = res_firm['t'][v]
        p = p_from_t(t)
        stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
        _print(f'  {v:<15} {b:+.6f}  ({s:.6f})  t={t:.3f}  p={p:.4f}{stars}')

    # Wald F-test (HC1)
    channel_vars_c = [v for v in CHANNEL_VARS if v in spec_vars_c]
    f_wald_firm, q_firm = wald_f_test(res_firm, channel_vars_c)
    if f_wald_firm is not None:
        _print(f'\n  Wald F-test (HC1): F = {f_wald_firm:.4f}, q = {q_firm}')

    # Difference test
    if 'w_geo' in res_firm['beta'] and 'w_fuel' in res_firm['beta']:
        names_c = ['intercept'] + spec_vars_c
        idx_geo_c = names_c.index('w_geo')
        idx_fuel_c = names_c.index('w_fuel')
        V_c = res_firm['V']
        diff_c = res_firm['beta']['w_geo'] - res_firm['beta']['w_fuel']
        se_diff_c = math.sqrt(V_c[idx_geo_c][idx_geo_c] + V_c[idx_fuel_c][idx_fuel_c]
                              - 2 * V_c[idx_geo_c][idx_fuel_c])
        t_diff_c = diff_c / se_diff_c if se_diff_c > 1e-15 else 0.0
        p_diff_c = p_from_t(t_diff_c)
        _print(f'\n  Difference test (beta_geo - beta_fuel): {diff_c:+.6f}')
        _print(f'  SE = {se_diff_c:.6f}, t = {t_diff_c:.3f}, p = {p_diff_c:.4f}')

# ── Panel D: Weighted firm-level (by sqrt of event count) ────────────

_print('\n' + '-' * 60)
_print('PANEL D: WEIGHTED FIRM-LEVEL (WLS, weight = sqrt(n_events))')
_print('-' * 60)

# WLS: multiply each obs by sqrt(weight)
for c in collapsed:
    c['wt'] = math.sqrt(c['n_events'])

weighted = []
for c in collapsed:
    w = c['wt']
    weighted.append({
        'car': c['car'] * w,
        'w_geo': c['w_geo'] * w,
        'w_fuel': c['w_fuel'] * w,
        'w_reg': c['w_reg'] * w,
        'same_sector': c['same_sector'] * w,
        'gvkey': c['gvkey'],
        '_weight': w,  # for intercept adjustment
    })

# For WLS, the intercept absorbs the weight, so we use the weighted data directly
# but we need to transform the intercept column too (it becomes w_i * 1 = w_i)
# This is handled by the OLS function since X[i][0] = 1.0 in the unweighted case
# For proper WLS: y* = w*y, X* = w*X, then OLS on (y*, X*) with no intercept offset
# Since our OLS adds intercept=1.0, we need a different approach:
# Just run OLS on the weighted data but replace the intercept column
n_w = len(weighted)
k_w = len(spec_vars_c) + 1
y_w = [c['car'] * c['wt'] for c in collapsed]
X_w = [[c['wt']] + [c[xv] * c['wt'] for xv in spec_vars_c] for c in collapsed]
y_w_mean = sum(y_w) / n_w
ss_tot_w = sum((yi - y_w_mean) ** 2 for yi in y_w)

XtX_w = [[sum(X_w[i][a] * X_w[i][b] for i in range(n_w)) for b in range(k_w)] for a in range(k_w)]
Xty_w = [sum(X_w[i][a] * y_w[i] for i in range(n_w)) for a in range(k_w)]
inv_XtX_w = invert_matrix(XtX_w)
if inv_XtX_w:
    beta_w = [sum(inv_XtX_w[a][b] * Xty_w[b] for b in range(k_w)) for a in range(k_w)]
    y_hat_w = [sum(X_w[i][a] * beta_w[a] for a in range(k_w)) for i in range(n_w)]
    resid_w = [y_w[i] - y_hat_w[i] for i in range(n_w)]
    ss_res_w = sum(r ** 2 for r in resid_w)

    # Compute R2 on the ORIGINAL (unweighted) scale
    y_hat_orig = [sum([1.0] + [collapsed[i][xv] for xv in spec_vars_c])[0]
                  if False else  # placeholder
                  beta_w[0] + sum(beta_w[j + 1] * collapsed[i][spec_vars_c[j]]
                                  for j in range(len(spec_vars_c)))
                  for i in range(n_w)]
    y_orig = [c['car'] for c in collapsed]
    y_orig_mean = sum(y_orig) / n_w
    ss_tot_orig = sum((yi - y_orig_mean) ** 2 for yi in y_orig)
    ss_res_orig = sum((y_orig[i] - y_hat_orig[i]) ** 2 for i in range(n_w))
    r2_wls = 1 - ss_res_orig / ss_tot_orig if ss_tot_orig > 0 else 0

    # HC1 SEs on weighted regression
    S_w = [[0.0] * k_w for _ in range(k_w)]
    for i in range(n_w):
        ri = resid_w[i]
        for a in range(k_w):
            for b in range(a, k_w):
                v = X_w[i][a] * X_w[i][b] * ri * ri
                S_w[a][b] += v
                if a != b:
                    S_w[b][a] += v
    V_w = mat_mul(mat_mul(inv_XtX_w, S_w), inv_XtX_w)
    hc1_w = n_w / (n_w - k_w)
    for a in range(k_w):
        for b in range(k_w):
            V_w[a][b] *= hc1_w
    se_w = [math.sqrt(V_w[a][a]) if V_w[a][a] > 0 else 0.0 for a in range(k_w)]
    t_w = [beta_w[a] / se_w[a] if se_w[a] > 1e-15 else 0.0 for a in range(k_w)]

    _print(f'  N = {n_w}, R2 (unweighted scale) = {r2_wls:.6f}')
    names_w = ['intercept'] + spec_vars_c
    for idx, v in enumerate(names_w):
        p = p_from_t(t_w[idx])
        stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
        _print(f'  {v:<15} {beta_w[idx]:+.6f}  ({se_w[idx]:.6f})  t={t_w[idx]:.3f}  p={p:.4f}{stars}')

# ── Summary comparison ───────────────────────────────────────────────

_print('\n' + '=' * 70)
_print('SUMMARY: R-SQUARED COMPARISON')
_print('=' * 70)
_print(f'  Pooled (N={n_pooled}):         R2 = {res_pooled["r2"]:.6f}')
_print(f'  Firm-level (N={len(collapsed)}): R2 = {res_firm["r2"]:.6f}')
if inv_XtX_w:
    _print(f'  WLS firm-level (N={n_w}):  R2 = {r2_wls:.6f}')
ratio = res_firm['r2'] / res_pooled['r2'] if res_pooled['r2'] > 0 else float('inf')
_print(f'  Ratio (firm / pooled): {ratio:.1f}x')

# ── Write output ─────────────────────────────────────────────────────

out_path = results_path('metrics', 'strategy2_firm_level_test.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

lines = [
    '# Firm-Level Collapsed Specification and Corrected F-Tests',
    '',
    'This analysis addresses two issues:',
    '1. The pooled R^2 (0.12%) is mechanically depressed because the same firm appears',
    '   ~179 times with nearly identical spatial weights but event-varying CARs.',
    '2. The original F-test used SSR-based (homoskedastic) statistic; this uses the',
    '   correct Wald statistic with cluster-robust variance-covariance matrix.',
    '',
    f'Pooled: N = {n_pooled} event-firm pairs, {n_events} events, {n_firms_pooled} firms',
    f'Collapsed: N = {len(collapsed)} firms',
    '',
    '## Panel A: Pooled Specification (Event-Clustered)',
    '',
    '| Variable | Beta | SE | t | p |',
    '|---|---:|---:|---:|---:|',
]

if res_pooled:
    for v in ['intercept'] + spec_vars:
        b = res_pooled['beta'][v]
        s = res_pooled['se'][v]
        t = res_pooled['t'][v]
        p = p_from_t(t)
        stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
        lines.append(f'| {v} | {b:+.6f} | {s:.6f} | {t:.3f} | {p:.4f}{stars} |')
    lines.append('')
    lines.append(f'R-squared = {res_pooled["r2"]:.6f}, N = {res_pooled["n"]}, '
                 f'Clusters = {res_pooled["clusters"]}')
    if f_wald is not None:
        lines.append(f'Wald F-test (event-clustered): F = {f_wald:.4f}, q = {q}')

lines += [
    '',
    '## Panel B: Pooled Specification (Two-Way Clustered: Event + Firm)',
    '',
    '| Variable | Beta | SE(two-way) | t | p |',
    '|---|---:|---:|---:|---:|',
]
for idx, v in enumerate(names_tw):
    p = p_from_t(t_tw[idx])
    stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
    lines.append(f'| {v} | {beta[idx]:+.6f} | {se_tw[idx]:.6f} | {t_tw[idx]:.3f} | {p:.4f}{stars} |')

lines.append('')
lines.append(f'Event clusters: {len(clust_event)}, Firm clusters: {len(clust_firm)}')
if inv_RVR_tw:
    lines.append(f'Wald F-test (two-way): F = {f_wald_tw:.4f}, q = {q_tw}')
lines.append(f'Difference test: beta_geo - beta_fuel = {diff_tw:+.6f} '
             f'(t = {t_diff_tw:.3f}, p = {p_diff_tw:.4f})')

lines += [
    '',
    '## Panel C: Firm-Level Collapsed Specification (HC1 Robust SEs)',
    '',
    '| Variable | Beta | SE(HC1) | t | p |',
    '|---|---:|---:|---:|---:|',
]
if res_firm:
    for v in ['intercept'] + spec_vars_c:
        b = res_firm['beta'][v]
        s = res_firm['se'][v]
        t = res_firm['t'][v]
        p = p_from_t(t)
        stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
        lines.append(f'| {v} | {b:+.6f} | {s:.6f} | {t:.3f} | {p:.4f}{stars} |')
    lines.append('')
    lines.append(f'R-squared = {res_firm["r2"]:.6f}, N = {res_firm["n"]}')
    if f_wald_firm is not None:
        lines.append(f'Wald F-test (HC1): F = {f_wald_firm:.4f}, q = {q_firm}')
    lines.append(f'Difference test: beta_geo - beta_fuel = {diff_c:+.6f} '
                 f'(t = {t_diff_c:.3f}, p = {p_diff_c:.4f})')

lines += [
    '',
    '## R-Squared Comparison',
    '',
    '| Specification | N | R-squared |',
    '|---|---:|---:|',
    f'| Pooled (event-firm pairs) | {n_pooled} | {res_pooled["r2"]:.6f} |',
    f'| Firm-level collapsed (OLS) | {len(collapsed)} | {res_firm["r2"]:.6f} |',
]
if inv_XtX_w:
    lines.append(f'| Firm-level collapsed (WLS) | {n_w} | {r2_wls:.6f} |')
lines.append('')
lines.append(f'Ratio (firm-level / pooled): {ratio:.1f}x')
lines.append('')
lines.append('## Interpretation')
lines.append('')
lines.append('The pooled R^2 is mechanically depressed by the panel structure: ')
lines.append('spatial weights vary across firms but not across events, while CARs ')
lines.append('vary substantially event-to-event. Collapsing to firm-level removes ')
lines.append('this noise inflation and provides a more interpretable R^2.')

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

_print(f'\nWrote: {out_path}')
_print('Done.')
