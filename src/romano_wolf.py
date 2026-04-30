"""Romano-Wolf (2005, 2016) stepdown correction for multiple hypothesis testing.

Tests 3 hypotheses simultaneously (3 spatial channels at the primary [-1,+3]
window) from Table 2:
  CAR_j = alpha + beta_geo * w^geo_ij + beta_fuel * w^fuel_ij
        + beta_reg * w^reg_ij + beta_s * SameSector_j + eps_j

for window [-1,+3] months.

Implements:
  1. Bonferroni correction: p_adj = min(p_raw * 3, 1)
  2. Bootstrap max-t (Westfall-Young): adjusted p = P(max|t*| >= |t_j|)
  3. Romano-Wolf stepdown: after rejecting most significant, recompute max
     over remaining hypotheses and enforce monotonicity

Reference: Romano, J.P. and Wolf, M. (2005). 'Stepwise Multiple Testing as
Formalized Data Snooping', Econometrica, 73(4), 1237-1282.
"""
import csv
import os
import sys
import math
import random
import hashlib
import shutil
import subprocess
import tempfile
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

# ── Configuration ────────────────────────────────────────────────────

B = 999           # bootstrap replications
SEED = 42         # reproducibility
MONTH_POSTS = [3]
PRE_MONTHS = 24   # pre-event months for AR demeaning

# Regression variables: 3 spatial channels + same_sector control
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


# ── OLS with event-clustered SEs ────────────────────────────────────

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


def ols(data, y_var, x_vars, cluster_var=None):
    """OLS regression returning betas, cluster-robust SEs, t-stats, R2."""
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
        se = [math.sqrt(s2 * inv_XtX[a][a]) if inv_XtX[a][a] > 0 else 0.0
              for a in range(k)]

    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0 for a in range(k)]

    names = ['intercept'] + x_vars
    return {
        'beta': dict(zip(names, beta)),
        'se': dict(zip(names, se)),
        't': dict(zip(names, t_stats)),
        'r2': r2,
        'n': n,
        'X': X,
        'resid': resid,
        'y': y,
        'y_hat': y_hat,
        'inv_XtX': inv_XtX,
    }


def bootstrap_t_stats(y_star, X, inv_XtX, cluster_map, x_vars, k, n,
                      beta_obs=None, se_obs=None):
    """Compute centred bootstrap t-stats: t*_j = (beta*_j - beta_j) / se_j.

    Uses the OBSERVED betas and SEs as centering and scaling, following
    Cameron, Gelbach & Miller (2008). This ensures the bootstrap distribution
    is centred at zero under H0 and uses the correct scale.

    If beta_obs/se_obs are None, falls back to beta*/se* (uncorrected).
    """
    # X'y*
    Xty = [0.0] * k
    for i in range(n):
        yi = y_star[i]
        Xi = X[i]
        for a in range(k):
            Xty[a] += Xi[a] * yi

    # beta* = inv_XtX @ X'y*
    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]

    names = ['intercept'] + x_vars

    if beta_obs is not None and se_obs is not None:
        # Centred t*: (beta* - beta_obs) / se_obs
        t_stats = []
        for a in range(k):
            se_a = se_obs.get(names[a], 0.0) if isinstance(se_obs, dict) else se_obs[a]
            b_a = beta_obs.get(names[a], 0.0) if isinstance(beta_obs, dict) else beta_obs[a]
            if se_a > 1e-15:
                t_stats.append((beta[a] - b_a) / se_a)
            else:
                t_stats.append(0.0)
        return dict(zip(names, t_stats))

    # Fallback: uncorrected beta*/se*
    resid = [0.0] * n
    for i in range(n):
        Xi = X[i]
        fitted = sum(Xi[a] * beta[a] for a in range(k))
        resid[i] = y_star[i] - fitted

    S = _cluster_cov(X, resid, cluster_map, k)
    V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
    G = len(cluster_map)
    if G > 1:
        scale = (G / (G - 1)) * ((n - 1) / (n - k))
        for a in range(k):
            for b in range(k):
                V[a][b] *= scale

    se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0.0 for a in range(k)]

    return dict(zip(names, t_stats))


# ── p-value from t-statistic (two-sided, normal approximation) ──────

def p_from_t(t_stat):
    """Two-sided p-value from t-stat using normal CDF approximation."""
    x = abs(t_stat)
    if x > 8:
        return 0.0
    b0 = 0.2316419
    b1 = 0.319381530
    b2 = -0.356563782
    b3 = 1.781477937
    b4 = -1.821255978
    b5 = 1.330274429
    t_val = 1.0 / (1.0 + b0 * x)
    phi = (1.0 / math.sqrt(2.0 * math.pi)) * math.exp(-0.5 * x * x)
    cdf = 1.0 - phi * (b1 * t_val + b2 * t_val**2 + b3 * t_val**3
                        + b4 * t_val**4 + b5 * t_val**5)
    return 2.0 * (1.0 - cdf)


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

# Alpha panel (time-varying fossil intensity)
alpha_panel = defaultdict(dict)
alpha_panel_path = derived_path('fundamentals', 'firm_alpha_panel.csv')
if os.path.exists(alpha_panel_path):
    with open(alpha_panel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            fy = row.get('fyear') or row.get('year')
            alpha = row.get('alpha', '')
            if gk and fy and alpha not in ('', None):
                alpha_panel[gk][str(fy)] = alpha

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

def compute_monthly_car(gvkey, event_month, post=12):
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


# ── Build regression datasets for all 3 windows ─────────────────────

def build_obs_for_window(post):
    """Build cross-sectional regression dataset for one monthly window."""
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

                car = compute_monthly_car(gk, event_month, post=post)
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
_print('ROMANO-WOLF MULTIPLE HYPOTHESIS TESTING CORRECTION')
_print(f'B = {B} bootstrap replications, seed = {SEED}')
_print(f'3 hypotheses: 3 channels x 1 window [-1,+3]')
_print('=' * 70)

# Step 1: Build datasets and estimate all 9 regressions
datasets = {}
results = {}
hypothesis_labels = []

for post in MONTH_POSTS:
    _print(f'\nBuilding dataset for window [-1, +{post}]...')
    obs = build_obs_for_window(post)
    datasets[post] = obs

    # Check same_sector variation
    ss_vals = set(o['same_sector'] for o in obs)
    if len(ss_vals) <= 1:
        spec_vars = ['w_geo', 'w_fuel', 'w_reg']
    else:
        spec_vars = SPEC_VARS

    res = ols(obs, 'car', spec_vars, cluster_var='event_id')
    if res is None:
        _print(f'  WARNING: OLS failed for post={post}')
        continue
    results[post] = res
    _print(f'  N={res["n"]}, R2={res["r2"]:.4f}')

    for ch in CHANNEL_VARS:
        hypothesis_labels.append((ch, post))
        t_val = res['t'].get(ch, 0.0)
        b_val = res['beta'].get(ch, 0.0)
        se_val = res['se'].get(ch, 0.0)
        _print(f'  {ch}: beta={b_val:+.6f}, se={se_val:.6f}, t={t_val:.2f}')

# Collect the 9 observed t-statistics
obs_t = []
for ch, post in hypothesis_labels:
    obs_t.append(results[post]['t'].get(ch, 0.0))

n_hyp = len(obs_t)
_print(f'\n{n_hyp} hypotheses collected')
for i, (ch, post) in enumerate(hypothesis_labels):
    _print(f'  H{i+1}: {ch} at [-1,+{post}], |t| = {abs(obs_t[i]):.3f}')

# Step 2: Raw p-values and Bonferroni correction
raw_p = [p_from_t(t) for t in obs_t]
bonf_p = [min(p * n_hyp, 1.0) for p in raw_p]

# Step 3: Bootstrap max-t (Westfall-Young) and Romano-Wolf stepdown

# Pre-compute per-window structures (X, inv_XtX, cluster_map, etc.)
window_cache = {}
for post in MONTH_POSTS:
    if post not in results:
        continue
    obs = datasets[post]
    res = results[post]

    ss_vals = set(o['same_sector'] for o in obs)
    spec_vars = SPEC_VARS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']

    cmap = defaultdict(list)
    cids = []
    for i, o in enumerate(obs):
        eid = o['event_id']
        cmap[eid].append(i)
        cids.append(eid)

    window_cache[post] = {
        'X': res['X'],
        'inv_XtX': res['inv_XtX'],
        'y_hat': res['y_hat'],
        'resid': res['resid'],
        'cluster_map': dict(cmap),
        'cluster_ids': cids,
        'spec_vars': spec_vars,
        'n': res['n'],
        'k': len(spec_vars) + 1,
        'beta_obs': res['beta'],   # observed betas for centred bootstrap
        'se_obs': res['se'],       # observed SEs for centred bootstrap
    }


# ── Julia hybrid bootstrap ──────────────────────────────────────────

def _find_julia():
    """Locate Julia executable. Returns path or None."""
    jl = shutil.which('julia')
    if jl:
        return jl
    fallback = os.path.expanduser(r'~\AppData\Local\Microsoft\WindowsApps\julia.exe')
    if os.path.isfile(fallback):
        return fallback
    return None


def _export_for_julia(tmpdir):
    """Export OLS matrices to CSV files that romano_wolf_bootstrap.jl expects."""
    data_dir = os.path.join(tmpdir, 'data')
    os.makedirs(data_dir, exist_ok=True)

    for post in MONTH_POSTS:
        if post not in window_cache:
            continue
        wc = window_cache[post]

        # X matrix with header
        k = wc['k']
        sv = wc['spec_vars']
        x_header = ['intercept'] + sv
        with open(os.path.join(data_dir, f'X_{post}.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(x_header)
            for row in wc['X']:
                w.writerow(row)

        # vectors: obs_idx, y_hat, resid, cluster_id
        with open(os.path.join(data_dir, f'vectors_{post}.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['obs_idx', 'y_hat', 'resid', 'cluster_id'])
            for i in range(wc['n']):
                w.writerow([i + 1, wc['y_hat'][i], wc['resid'][i],
                            wc['cluster_ids'][i]])

        # inv_XtX with header
        inv_header = [f'c{j}' for j in range(k)]
        with open(os.path.join(data_dir, f'inv_XtX_{post}.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(inv_header)
            for row in wc['inv_XtX']:
                w.writerow(row)

        # spec_vars (one per line, no header)
        with open(os.path.join(data_dir, f'spec_vars_{post}.txt'), 'w') as f:
            for v in sv:
                f.write(v + '\n')

        # Observed betas and SEs for centred bootstrap
        beta_obs = wc['beta_obs']
        se_obs = wc['se_obs']
        names = ['intercept'] + sv
        with open(os.path.join(data_dir, f'obs_beta_se_{post}.csv'), 'w', newline='') as f:
            w = csv.writer(f)
            w.writerow(['beta_obs', 'se_obs'])
            for name in names:
                b = beta_obs.get(name, 0.0) if isinstance(beta_obs, dict) else 0.0
                s = se_obs.get(name, 0.0) if isinstance(se_obs, dict) else 0.0
                w.writerow([b, s])

    return data_dir


def _run_julia_bootstrap():
    """Try to run bootstrap via Julia. Returns (maxt_p, rw_p) or None."""
    julia_exe = _find_julia()
    if julia_exe is None:
        _print('  Julia not found, falling back to Python bootstrap.')
        return None

    jl_script = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'romano_wolf_bootstrap.jl')
    if not os.path.isfile(jl_script):
        _print(f'  Julia script not found at {jl_script}, falling back.')
        return None

    tmpdir = tempfile.mkdtemp(prefix='rw_bootstrap_')
    try:
        _print(f'  Exporting matrices to {tmpdir}...')
        data_dir = _export_for_julia(tmpdir)
        out_csv = os.path.join(tmpdir, 'julia_rw_results.csv')

        _print(f'  Calling Julia: {julia_exe}')
        proc = subprocess.run(
            [julia_exe, jl_script, data_dir, out_csv],
            capture_output=True, text=True, timeout=600,
        )
        if proc.returncode != 0:
            _print(f'  Julia failed (rc={proc.returncode}):')
            for line in (proc.stderr or '').strip().split('\n')[:10]:
                _print(f'    {line}')
            return None

        # Print Julia stdout for diagnostics
        for line in (proc.stdout or '').strip().split('\n'):
            _print(f'  [julia] {line}')

        # Read results CSV: channel,window,obs_t,maxt_p,rw_p
        if not os.path.isfile(out_csv):
            _print('  Julia output CSV not found.')
            return None

        julia_maxt = {}
        julia_rw = {}
        with open(out_csv, 'r', encoding='utf-8') as f:
            for row in csv.DictReader(f):
                key = (row['channel'], int(row['window']))
                julia_maxt[key] = float(row['maxt_p'])
                julia_rw[key] = float(row['rw_p'])

        # Build p-value lists in hypothesis_labels order
        maxt_p = []
        rw_p = []
        for ch, post in hypothesis_labels:
            key = (ch, post)
            if key not in julia_maxt:
                _print(f'  WARNING: Julia missing result for {key}')
                return None
            maxt_p.append(julia_maxt[key])
            rw_p.append(julia_rw[key])

        return maxt_p, rw_p

    except subprocess.TimeoutExpired:
        _print('  Julia timed out after 600s, falling back.')
        return None
    except Exception as e:
        _print(f'  Julia error: {e}')
        return None
    finally:
        # Clean up temp directory
        try:
            shutil.rmtree(tmpdir, ignore_errors=True)
        except Exception:
            pass


def _python_bootstrap():
    """Pure-Python fallback bootstrap. Returns (maxt_p, rw_p)."""
    _print(f'\nRunning {B} Python bootstrap replications...')
    random.seed(SEED)

    # Get all unique cluster IDs across windows for consistent Rademacher draws
    all_cluster_ids = set()
    for post in MONTH_POSTS:
        if post in window_cache:
            all_cluster_ids.update(window_cache[post]['cluster_map'].keys())
    all_cluster_ids = sorted(all_cluster_ids)

    import gc
    gc.collect()

    # Pre-allocate y_star buffers per window
    y_star_buf = {}
    for post in MONTH_POSTS:
        if post in window_cache:
            y_star_buf[post] = [0.0] * window_cache[post]['n']

    boot_max_t = [0.0] * B
    boot_t_all = [[0.0] * n_hyp for _ in range(B)]

    for b_iter in range(B):
        if (b_iter + 1) % 100 == 0:
            _print(f'  bootstrap {b_iter + 1}/{B}')

        rademacher = {cid: (1 if random.random() < 0.5 else -1)
                      for cid in all_cluster_ids}

        col = 0
        for post in MONTH_POSTS:
            if post not in window_cache:
                col += len(CHANNEL_VARS)
                continue

            wc = window_cache[post]
            y_hat = wc['y_hat']
            resid = wc['resid']
            cids = wc['cluster_ids']
            n = wc['n']

            buf = y_star_buf[post]
            for i in range(n):
                buf[i] = y_hat[i] + rademacher.get(cids[i], 1) * resid[i]

            t_dict = bootstrap_t_stats(
                buf, wc['X'], wc['inv_XtX'], wc['cluster_map'],
                wc['spec_vars'], wc['k'], n,
                beta_obs=wc['beta_obs'], se_obs=wc['se_obs']
            )
            if t_dict is None:
                col += len(CHANNEL_VARS)
                continue

            for ch in CHANNEL_VARS:
                boot_t_all[b_iter][col] = abs(t_dict.get(ch, 0.0))
                col += 1

        boot_max_t[b_iter] = max(boot_t_all[b_iter])

    _print('  Bootstrap complete.')

    # Westfall-Young max-t adjusted p-values
    maxt_p = []
    for j in range(n_hyp):
        abs_t_j = abs(obs_t[j])
        count = sum(1 for mt in boot_max_t if mt >= abs_t_j)
        maxt_p.append(count / B)

    # Romano-Wolf stepdown p-values
    order = sorted(range(n_hyp), key=lambda j: abs(obs_t[j]), reverse=True)
    rw_p = [0.0] * n_hyp
    remaining = list(range(n_hyp))

    for step_idx, j in enumerate(order):
        abs_t_j = abs(obs_t[j])
        count = 0
        for b_iter in range(B):
            step_max = max(boot_t_all[b_iter][h] for h in remaining)
            if step_max >= abs_t_j:
                count += 1
        p_step = count / B

        if step_idx > 0:
            p_step = max(p_step, rw_p[order[step_idx - 1]])
        rw_p[j] = p_step

        remaining.remove(j)
        if not remaining:
            break

    return maxt_p, rw_p


# Try Julia first, fall back to Python
_print(f'\nRunning {B} bootstrap replications...')
julia_result = _run_julia_bootstrap()

if julia_result is not None:
    maxt_p, rw_p = julia_result
    _print('  Using Julia bootstrap p-values (max-t and Romano-Wolf).')
else:
    maxt_p, rw_p = _python_bootstrap()

# Step 6: Output results
_print()
_print('Results:')
_print(f'  {"Variable":<10} {"Window":<12} {"t":>8} {"raw p":>8} '
       f'{"Bonf p":>8} {"Max-t p":>8} {"RW p":>8}')
_print(f'  {"-" * 66}')
for j, (ch, post) in enumerate(hypothesis_labels):
    window_str = f'[-1,+{post}]'
    _print(f'  {ch:<10} {window_str:<12} {obs_t[j]:>8.3f} {raw_p[j]:>8.4f} '
           f'{bonf_p[j]:>8.4f} {maxt_p[j]:>8.4f} {rw_p[j]:>8.4f}')

# Write markdown output
out_path = results_path('metrics', 'romano_wolf.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)

lines = [
    '# Multiple Hypothesis Testing Correction',
    '',
    'Romano-Wolf (2005, 2016) stepdown correction for the 3 coefficients',
    'in Table 2: 3 spatial channels (w_geo, w_fuel, w_reg) at the primary',
    '[-1,+3] month window.',
    '',
    'Specification for each window:',
    '  CAR_j = alpha + beta_geo * w^geo_ij + beta_fuel * w^fuel_ij',
    '        + beta_reg * w^reg_ij + beta_s * SameSector_j + eps_j',
    '',
    f'Bootstrap replications: B = {B}',
    f'Seed: {SEED}',
    f'Events: {len(all_events)} first-mover coal retirements',
    f'Return model: vwretd (market-adjusted, Fama-French)',
    '',
    'Methods:',
    '- Raw p: two-sided p-value from asymptotic normal',
    '- Bonferroni: p_adj = min(p_raw x 3, 1)',
    '- Max-t (Westfall-Young): P(max|t*| >= |t_j|) using Rademacher',
    '  cluster bootstrap under H0',
    '- Romano-Wolf: stepdown refinement of max-t; after rejecting the',
    '  most significant hypothesis, recompute max over remaining',
    '  hypotheses and enforce monotonicity',
    '',
    '## Results',
    '',
    '| Variable | Window | Beta | SE | Raw t | Raw p | Bonferroni p | Max-t p | RW p |',
    '|---|---|---:|---:|---:|---:|---:|---:|---:|',
]

for j, (ch, post) in enumerate(hypothesis_labels):
    window_str = f'[-1,+{post}]'
    b_val = results[post]['beta'].get(ch, 0.0)
    se_val = results[post]['se'].get(ch, 0.0)
    sig_raw = ('***' if raw_p[j] < 0.01
               else '**' if raw_p[j] < 0.05
               else '*' if raw_p[j] < 0.10 else '')
    sig_rw = ('***' if rw_p[j] < 0.01
              else '**' if rw_p[j] < 0.05
              else '*' if rw_p[j] < 0.10 else '')
    lines.append(
        f'| {ch} | {window_str} | {b_val:+.6f} | {se_val:.6f} '
        f'| {obs_t[j]:.3f} | {raw_p[j]:.4f}{sig_raw} '
        f'| {bonf_p[j]:.4f} | {maxt_p[j]:.4f} | {rw_p[j]:.4f}{sig_rw} |'
    )

lines.append('')
lines.append('## Interpretation')
lines.append('')

# Count rejections at 5%
rej_raw = sum(1 for p in raw_p if p < 0.05)
rej_bonf = sum(1 for p in bonf_p if p < 0.05)
rej_maxt = sum(1 for p in maxt_p if p < 0.05)
rej_rw = sum(1 for p in rw_p if p < 0.05)

lines.append(f'Rejections at 5% level:')
lines.append(f'- Raw: {rej_raw}/{n_hyp}')
lines.append(f'- Bonferroni: {rej_bonf}/{n_hyp}')
lines.append(f'- Max-t (Westfall-Young): {rej_maxt}/{n_hyp}')
lines.append(f'- Romano-Wolf stepdown: {rej_rw}/{n_hyp}')
lines.append('')

# Note on FWER
lines.append(f'Without correction, testing {n_hyp} hypotheses at 5% yields a')
lines.append(f'family-wise error rate of 1 - (1-0.05)^{n_hyp} = {1 - 0.95**n_hyp:.1%}.')
lines.append('The corrections above control the FWER at the nominal level.')
lines.append('')

# Sample sizes per window
lines.append('## Sample sizes')
lines.append('')
for post in MONTH_POSTS:
    if post in results:
        lines.append(f'- Window [-1,+{post}]: N = {results[post]["n"]}, '
                      f'R2 = {results[post]["r2"]:.4f}')
lines.append('')

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

_print(f'\nWrote: {out_path}')
