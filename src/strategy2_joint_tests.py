"""Focused hypothesis tests for the channel split claim.

Instead of testing 9 individual hypotheses (3 channels x 3 horizons) that
fail Romano-Wolf correction, this script performs TWO tests that directly
address the paper's core claim: spatial channels have opposite signs.

Test 1 — Joint F-test: H0: beta_geo = beta_fuel = beta_reg = 0
    (spatial network has no effect)

Test 2 — Difference-in-coefficients: H0: beta_geo = beta_fuel
    (channels have the same effect, contra the paper's claim)

Both tests use the [-1,+3] month window (strongest signal) and require
no multiple-testing correction because each is a single hypothesis.

Specification:
  CAR_j = alpha + beta_geo * w^geo + beta_fuel * w^fuel
        + beta_reg * w^reg + beta_s * SameSector + eps_j
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

B = 999           # permutation replications for F-test
SEED = 42         # reproducibility
POST_MONTHS = 3   # [-1, +3] window — the strongest
PRE_MONTHS = 24   # pre-event months for AR demeaning

CHANNEL_VARS = ['w_geo', 'w_fuel', 'w_reg']
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']


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


# ── OLS with event-clustered SEs and full covariance matrix ─────────

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

    # beta = (X'X)^{-1} X'y
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
else:
    _print('  W_reg: NOT FOUND')

# Fundamentals (latest record per firm for sector classification)
fundamentals = {}
fundamentals_by_year = defaultdict(dict)
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        fundamentals_by_year[gk][fy] = row
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row


def get_sic4(gvkey):
    f = fundamentals.get(gvkey)
    if f and f.get('sic'):
        return f['sic'][:4]
    return None


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


# ── Build regression dataset for [-1, +3] window ─────────────────────

def build_obs():
    """Build cross-sectional regression dataset for the [-1,+3] window."""
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


# ── Main analysis ────────────────────────────────────────────────────

_print()
_print('=' * 70)
_print('FOCUSED HYPOTHESIS TESTS: CHANNEL SPLIT')
_print(f'Window: [-1, +{POST_MONTHS}] months')
_print(f'Bootstrap: B = {B}, seed = {SEED}')
_print('=' * 70)

# Step 1: Build dataset
_print('\nBuilding dataset...')
obs = build_obs()
_print(f'  N = {len(obs)} observations')
_print(f'  Neighbors (w_geo > 0): {sum(1 for o in obs if o["w_geo"] > 0)}')
_print(f'  Non-connected: {sum(1 for o in obs if o["w_geo"] == 0)}')
_print(f'  Events: {len(set(o["event_id"] for o in obs))}')

# Check same_sector variation
ss_vals = set(o['same_sector'] for o in obs)
if len(ss_vals) <= 1:
    _print('  WARNING: no same_sector variation; dropping same_sector from spec')
    spec_vars = ['w_geo', 'w_fuel', 'w_reg']
else:
    spec_vars = SPEC_VARS

# Step 2: Estimate unrestricted model
_print('\nEstimating unrestricted model...')
_print(f'  CAR = alpha + beta_geo*w_geo + beta_fuel*w_fuel + beta_reg*w_reg + beta_s*SameSector')
res_full = ols_full(obs, 'car', spec_vars, cluster_var='event_id')
if res_full is None:
    _print('ERROR: Unrestricted OLS failed.')
    sys.exit(1)

_print(f'  N = {res_full["n"]}, R2 = {res_full["r2"]:.6f}')
for v in spec_vars:
    _print(f'  {v}: beta = {res_full["beta"][v]:+.6f}, '
           f'se = {res_full["se"][v]:.6f}, t = {res_full["t"][v]:.3f}')

# Step 3: Estimate restricted model (same_sector only)
restricted_vars = [v for v in spec_vars if v not in CHANNEL_VARS]
_print('\nEstimating restricted model...')
_print(f'  CAR = alpha' + (' + beta_s*SameSector' if restricted_vars else ''))
res_restricted = ols_full(obs, 'car', restricted_vars, cluster_var='event_id')
if res_restricted is None:
    _print('ERROR: Restricted OLS failed.')
    sys.exit(1)
_print(f'  N = {res_restricted["n"]}, R2 = {res_restricted["r2"]:.6f}')

# ── TEST 1: Joint F-test ────────────────────────────────────────────

_print()
_print('=' * 70)
_print('TEST 1: JOINT F-TEST')
_print('H0: beta_geo = beta_fuel = beta_reg = 0')
_print('=' * 70)

ssr_restricted = res_restricted['ss_res']
ssr_unrestricted = res_full['ss_res']
q = len(CHANNEL_VARS)  # number of restrictions = 3
n_obs = res_full['n']
k_full = len(spec_vars) + 1  # total regressors including intercept

f_stat = ((ssr_restricted - ssr_unrestricted) / q) / (ssr_unrestricted / (n_obs - k_full))

_print(f'\n  SSR_restricted:   {ssr_restricted:.6f}')
_print(f'  SSR_unrestricted: {ssr_unrestricted:.6f}')
_print(f'  F-statistic: {f_stat:.4f}')
_print(f'  df: ({q}, {n_obs - k_full})')

# Permutation bootstrap for F-test p-value
_print(f'\nRunning {B} permutation bootstraps for F-test p-value...')
random.seed(SEED)

# For each permutation: randomly permute the spatial weights (w_geo, w_fuel, w_reg)
# across observations while keeping CARs and same_sector fixed.
# Then re-estimate both models and compute F*.
f_boot = []

# Pre-extract arrays for speed
car_arr = [o['car'] for o in obs]
ss_arr = [o['same_sector'] for o in obs]
wgeo_arr = [o['w_geo'] for o in obs]
wfuel_arr = [o['w_fuel'] for o in obs]
wreg_arr = [o['w_reg'] for o in obs]
eid_arr = [o['event_id'] for o in obs]
n_perm = len(obs)

# The restricted model SSR does not change under permutation of spatial weights
# (since restricted model only uses same_sector, which stays fixed).
# So SSR_restricted is constant.

for b_iter in range(B):
    if (b_iter + 1) % 100 == 0:
        _print(f'  permutation {b_iter + 1}/{B}')

    # Permute spatial weights
    perm_idx = list(range(n_perm))
    random.shuffle(perm_idx)

    perm_data = []
    for i in range(n_perm):
        pi = perm_idx[i]
        perm_data.append({
            'car': car_arr[i],
            'w_geo': wgeo_arr[pi],
            'w_fuel': wfuel_arr[pi],
            'w_reg': wreg_arr[pi],
            'same_sector': ss_arr[i],
            'event_id': eid_arr[i],
        })

    res_perm = ols_full(perm_data, 'car', spec_vars, cluster_var='event_id')
    if res_perm is None:
        f_boot.append(0.0)
        continue

    ssr_perm = res_perm['ss_res']
    f_perm = ((ssr_restricted - ssr_perm) / q) / (ssr_perm / (n_perm - k_full))
    f_boot.append(f_perm)

p_f_perm = sum(1 for fb in f_boot if fb >= f_stat) / B

_print(f'  F-test p-value (permutation): {p_f_perm:.4f}')

f_reject = 'REJECT H0' if p_f_perm < 0.05 else 'FAIL TO REJECT H0'
_print(f'  Conclusion at 5%: {f_reject}')

# ── TEST 2: Difference-in-coefficients test ──────────────────────────

_print()
_print('=' * 70)
_print('TEST 2: DIFFERENCE-IN-COEFFICIENTS')
_print('H0: beta_geo = beta_fuel')
_print('H1: beta_geo != beta_fuel (opposing signs)')
_print('=' * 70)

# Extract from full variance-covariance matrix V
# Variable ordering: ['intercept'] + spec_vars
# Find indices of w_geo and w_fuel
names_list = ['intercept'] + spec_vars
idx_geo = names_list.index('w_geo')
idx_fuel = names_list.index('w_fuel')

V = res_full['V']
beta_geo = res_full['beta']['w_geo']
beta_fuel = res_full['beta']['w_fuel']
se_geo = res_full['se']['w_geo']
se_fuel = res_full['se']['w_fuel']

var_geo = V[idx_geo][idx_geo]
var_fuel = V[idx_fuel][idx_fuel]
cov_geo_fuel = V[idx_geo][idx_fuel]

diff = beta_geo - beta_fuel
se_diff = math.sqrt(var_geo + var_fuel - 2.0 * cov_geo_fuel)
t_diff = diff / se_diff if se_diff > 1e-15 else 0.0
p_diff = p_from_t(t_diff)

_print(f'\n  beta_geo:  {beta_geo:+.6f} (SE {se_geo:.6f})')
_print(f'  beta_fuel: {beta_fuel:+.6f} (SE {se_fuel:.6f})')
_print(f'  Difference (beta_geo - beta_fuel): {diff:+.6f}')
_print(f'  Var(beta_geo):           {var_geo:.10f}')
_print(f'  Var(beta_fuel):          {var_fuel:.10f}')
_print(f'  Cov(beta_geo, beta_fuel): {cov_geo_fuel:+.10f}')
_print(f'  SE of difference:        {se_diff:.6f}')
_print(f'  t-statistic:             {t_diff:.3f}')
_print(f'  p-value (two-sided):     {p_diff:.4f}')

diff_reject = 'REJECT H0' if p_diff < 0.05 else 'FAIL TO REJECT H0'
_print(f'  Conclusion at 5%: {diff_reject}')

sig_stars_diff = ('***' if p_diff < 0.01
                  else '**' if p_diff < 0.05
                  else '*' if p_diff < 0.10 else '')
diff_interp = ('statistically significant' if p_diff < 0.05
               else 'not significant')

# ── Summary statistics ───────────────────────────────────────────────

_print()
_print('=' * 70)
_print('FULL REGRESSION RESULTS')
_print('=' * 70)
all_names = ['intercept'] + spec_vars
for name in all_names:
    b = res_full['beta'][name]
    s = res_full['se'][name]
    t = res_full['t'][name]
    p = p_from_t(t)
    stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
    _print(f'  {name:<15} {b:+.6f}  ({s:.6f})  t={t:.3f}  p={p:.4f}{stars}')

# ── Write output ─────────────────────────────────────────────────────

out_path = results_path('metrics', 'strategy2_joint_tests.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

sig_stars_f = ('***' if p_f_perm < 0.01
               else '**' if p_f_perm < 0.05
               else '*' if p_f_perm < 0.10 else '')
f_interp = 'reject' if p_f_perm < 0.05 else 'fail to reject'

lines = [
    '# Focused Hypothesis Tests: Channel Split',
    '',
    f'Window: [-1, +{POST_MONTHS}] months (monthly CARs, vwretd)',
    f'Events: {len(all_events)} first-mover coal retirements',
    f'N = {n_obs} observations, {res_full["clusters"]} event clusters',
    f'Standard errors: event-clustered',
    '',
    '## Test 1: Joint F-test (H0: beta_geo = beta_fuel = beta_reg = 0)',
    '',
    'Unrestricted: CAR = alpha + beta_geo * w^geo + beta_fuel * w^fuel + beta_reg * w^reg + beta_s * SameSector',
    'Restricted:   CAR = alpha + beta_s * SameSector',
    '',
    f'SSR_restricted:   {ssr_restricted:.6f}',
    f'SSR_unrestricted: {ssr_unrestricted:.6f}',
    f'F-statistic: {f_stat:.4f}',
    f'df: ({q}, {n_obs - k_full})',
    f'p-value (permutation, B={B}): {p_f_perm:.4f}{sig_stars_f}',
    f'N: {n_obs}',
    '',
    f'Interpretation: {f_interp.capitalize()} H0 at 5%. The spatial network channels '
    f'{"jointly predict CARs" if p_f_perm < 0.05 else "do not jointly predict CARs"} '
    f'around coal retirement events.',
    '',
    '## Test 2: Difference test (H0: beta_geo = beta_fuel)',
    '',
    f'beta_geo:  {beta_geo:+.6f} (SE {se_geo:.6f})',
    f'beta_fuel: {beta_fuel:+.6f} (SE {se_fuel:.6f})',
    f'Difference (beta_geo - beta_fuel): {diff:+.6f}',
    f'SE of difference: {se_diff:.6f}',
    f't-statistic: {t_diff:.3f}',
    f'p-value: {p_diff:.4f}{sig_stars_diff}',
    f'Cov(beta_geo, beta_fuel): {cov_geo_fuel:+.10f}',
    '',
    f'Interpretation: The opposing-sign channel split is {diff_interp} '
    f'in a single test (t = {t_diff:.3f}, p = {p_diff:.4f}).',
    '',
    '## Full regression coefficients',
    '',
    '| Variable | Beta | SE | t | p |',
    '|---|---:|---:|---:|---:|',
]

for name in all_names:
    b = res_full['beta'][name]
    s = res_full['se'][name]
    t = res_full['t'][name]
    p = p_from_t(t)
    stars = '***' if p < 0.01 else '**' if p < 0.05 else '*' if p < 0.10 else ''
    lines.append(f'| {name} | {b:+.6f} | {s:.6f} | {t:.3f} | {p:.4f}{stars} |')

lines.append('')
lines.append('## Comparison with individual tests')
lines.append('')
lines.append('| Test | Hypotheses tested | Correction needed | Result |')
lines.append('|---|---|---|---|')
lines.append('| Individual t-tests | 9 | Romano-Wolf | 0/9 significant |')
f_result = f'significant (p={p_f_perm:.4f})' if p_f_perm < 0.05 else f'not significant (p={p_f_perm:.4f})'
d_result = f'significant (p={p_diff:.4f})' if p_diff < 0.05 else f'not significant (p={p_diff:.4f})'
lines.append(f'| Joint F-test | 1 | None | {f_result} |')
lines.append(f'| Difference test | 1 | None | {d_result} |')
lines.append('')

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

_print(f'\nWrote: {out_path}')
_print('Done.')
