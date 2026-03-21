"""Strategy 2 regression: Spatial event study specification from paper.

For each first-mover event i, pool all firms j and estimate:
  CAR_j = alpha + beta_1 * w_ij + beta_2 * (w_ij * ConditionalTransformer_j) + beta_3 * SameSector_j + eps_j

A positive beta_1 with positive beta_2 increment is the coordination signature.
Uses geographic-only W and both binary and continuous conditionality measures.
Supports first-movers only or all matched retirements, and event-year fundamentals.
Runs on both daily and monthly CARs under multiple return models for robustness:
  - vwretd: market-adjusted using CRSP value-weighted return (Fama-French)
  - capm: CAPM-adjusted using pre-event beta and Fama-French factors
  - constant_mean: firm's own pre-event mean return subtracted
"""
import csv, os, math, hashlib, datetime
from collections import defaultdict

from _paths import raw_path, derived_path

# Classification modes:
# - "strict": alpha in [0.3,0.7] AND lambda<=median AND kappa>=median
# - "alpha_only": alpha in [0.3,0.7] only
CLASS_MODES = ['alpha_only', 'strict']

# Event scope:
# - "first_mover": only first-mover retirements (clean identification)
# - "all_matched": all matched retirements (power)
EVENT_SCOPES = ['first_mover', 'all_matched']

# Restrict to events with researched announcement dates only
EXACT_ONLY = os.getenv('EXACT_ONLY', '0') == '1'

# Transformations for exposure variables (diagnostics / robustness)
TRANSFORM_SET = os.getenv('TRANSFORM_SET', 'base')  # base | log1p | zscore
WRITE_METRICS = os.getenv('WRITE_METRICS', '0') == '1'

# Monthly windows for robustness
MONTH_POSTS = [12, 6, 3]

# Return models for CAR computation:
# - "vwretd": AR_t = R_it - R_mt (CRSP value-weighted, from Fama-French)
# - "capm": AR_t = R_it - [RF + beta*(Mkt-RF)] using pre-event beta
# - "constant_mean": AR_t = R_it - mean(R_i, pre-event window)
RETURN_MODELS = ['vwretd', 'capm', 'constant_mean']

# Light-run mode to speed up iteration (set RUN_LIGHT=1)
if os.getenv('RUN_LIGHT', '0') == '1':
    CLASS_MODES = ['alpha_only']
    EVENT_SCOPES = ['first_mover']
    RETURN_MODELS = ['vwretd']

# If RUN_FE_ONLY=1, skip all but the primary spec (alpha_only, first_mover, vwretd, post=12)
RUN_FE_ONLY = os.getenv('RUN_FE_ONLY', '0') == '1'

# Channel decomposition (geo/reg/fuel) for mechanism mapping
# Default: run only for primary spec (alpha_only + vwretd + first_mover)
RUN_CHANNEL_DECOMP = os.getenv('RUN_CHANNEL_DECOMP', '1') != '0'

# Lead-time interaction (announcement -> retirement) for credibility timing
RUN_LEAD_TIME = os.getenv('RUN_LEAD_TIME', '0') == '1'

# Pre-trend adjustment (demean CAR by pre-event AR mean)
PRE_DEMEAN_DAILY = True
PRE_DEMEAN_MONTHLY = True
PRE_DAYS = 250
PRE_MONTHS = 24

# Clustering modes for inference
CLUSTER_MODES = ['event', 'twoway']


# -- OLS with cluster-robust standard errors --

def invert_matrix(mat):
    n = len(mat)
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(mat)]
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


def summarize(vals):
    vals = [v for v in vals if v is not None and not math.isnan(v)]
    if not vals:
        return None
    vals.sort()
    n = len(vals)
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / (n - 1) if n > 1 else 0
    sd = math.sqrt(var)
    p = lambda q: vals[int(q * (n - 1))]
    if sd > 0 and n > 2:
        m3 = sum((v - mean) ** 3 for v in vals) / n
        skew = m3 / (sd ** 3)
    else:
        skew = 0.0
    if sd > 0 and n > 3:
        m4 = sum((v - mean) ** 4 for v in vals) / n
        kurt = m4 / (sd ** 4) - 3
    else:
        kurt = 0.0
    return {
        'N': n, 'Mean': mean, 'SD': sd, 'Min': vals[0], 'P1': p(0.01), 'P5': p(0.05),
        'P25': p(0.25), 'Median': p(0.50), 'P75': p(0.75), 'P95': p(0.95), 'P99': p(0.99),
        'Max': vals[-1], 'Skew': skew, 'Kurt': kurt,
    }


def write_metrics_md(path, metrics, extra_lines=None):
    lines = []
    lines.append('| Variable | N | Mean | SD | Min | P1 | P5 | P25 | Median | P75 | P95 | P99 | Max | Skew | Kurt |')
    lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|')
    for var, s in metrics.items():
        if not s:
            continue
        lines.append(
            f"| {var} | {s['N']} | {s['Mean']:.3f} | {s['SD']:.3f} | {s['Min']:.3f} | {s['P1']:.3f} | {s['P5']:.3f} | "
            f"{s['P25']:.3f} | {s['Median']:.3f} | {s['P75']:.3f} | {s['P95']:.3f} | {s['P99']:.3f} | "
            f"{s['Max']:.3f} | {s['Skew']:.2f} | {s['Kurt']:.2f} |"
        )
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        if extra_lines:
            lines.append('')
            lines.extend(extra_lines)
        f.write('\n'.join(lines))


def apply_exposure_transform(obs_list, transform_set):
    if not obs_list:
        return
    if transform_set == 'base':
        return
    keys = ['w_ij', 'w_geo', 'w_reg', 'w_fuel']
    if transform_set == 'log1p':
        for o in obs_list:
            for k in keys:
                if k in o and o[k] is not None:
                    o[k] = math.log1p(max(o[k], 0.0))
    elif transform_set == 'zscore':
        stats = {}
        for k in keys:
            vals = [o[k] for o in obs_list if k in o and o[k] is not None]
            if len(vals) < 2:
                continue
            mean = sum(vals) / len(vals)
            var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
            sd = math.sqrt(var) if var > 0 else None
            stats[k] = (mean, sd)
        for o in obs_list:
            for k, (mean, sd) in stats.items():
                if sd and k in o and o[k] is not None:
                    o[k] = (o[k] - mean) / sd
    # refresh interactions using transformed exposures
    for o in obs_list:
        w = o.get('w_ij', 0.0)
        w_geo = o.get('w_geo', 0.0)
        w_reg = o.get('w_reg', 0.0)
        w_fuel = o.get('w_fuel', 0.0)
        cond = o.get('conditional', 0.0)
        cs = o.get('cond_score', 0.0)
        csr = o.get('cond_score_rank', 0.0)
        csk = o.get('cond_score_kernel', 0.0)
        csf = o.get('cond_score_fin', 0.0)
        lt = o.get('lead_time_z', None)
        o['w_x_cond'] = w * cond
        o['w_x_condscore'] = w * cs
        o['w_x_condscore_rank'] = w * csr
        o['w_x_condscore_kernel'] = w * csk
        o['w_x_condscore_fin'] = w * csf
        o['w_geo_x_cond'] = w_geo * cond
        o['w_reg_x_cond'] = w_reg * cond
        o['w_fuel_x_cond'] = w_fuel * cond
        o['w_geo_x_condscore'] = w_geo * cs
        o['w_reg_x_condscore'] = w_reg * cs
        o['w_fuel_x_condscore'] = w_fuel * cs
        if lt is not None:
            o['w_x_lead_time'] = w * lt


def mat_mul(a, b):
    rows = len(a)
    cols = len(b[0])
    mid = len(b)
    out = [[0.0 for _ in range(cols)] for _ in range(rows)]
    for i in range(rows):
        for k in range(mid):
            aik = a[i][k]
            if aik == 0:
                continue
            for j in range(cols):
                out[i][j] += aik * b[k][j]
    return out


def _cluster_cov(X, resid, clusters):
    n = len(X)
    k = len(X[0])
    S = [[0.0 for _ in range(k)] for _ in range(k)]
    for _, idxs in clusters.items():
        xu = [0.0 for _ in range(k)]
        for i in idxs:
            for a in range(k):
                xu[a] += X[i][a] * resid[i]
        for a in range(k):
            for b in range(k):
                S[a][b] += xu[a] * xu[b]
    return S, len(clusters), n


def ols(data, y_var, x_vars, cluster_var=None, debug=False):
    """OLS regression returning betas, cluster-robust standard errors, t-stats, R2.

    cluster_var can be a string (single cluster) or a list/tuple of two cluster keys
    for two-way clustering (Cameron-Gelbach-Miller).
    """
    n = len(data)
    k = len(x_vars) + 1
    if n <= k + 1:
        if debug:
            print(f'  [OLS] n={n} <= k+1={k+1}, skipping')
        return None

    y = [d[y_var] for d in data]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    if ss_tot < 1e-15:
        if debug:
            print(f'  [OLS] ss_tot={ss_tot:.2e} too small')
        return None

    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]

    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

    if debug:
        print(f'  [OLS] n={n}, k={k}, ss_tot={ss_tot:.4f}')
        print(f"  [OLS] X'X diagonal: {[XtX[a][a] for a in range(k)]}")

    aug = [row[:] + [Xty[a]] for a, row in enumerate(XtX)]
    for col in range(k):
        max_row = max(range(col, k), key=lambda r: abs(aug[r][col]))
        aug[col], aug[max_row] = aug[max_row], aug[col]
        pivot = aug[col][col]
        if abs(pivot) < 1e-20:
            if debug:
                print(f'  [OLS] Singular at col {col}, pivot={pivot:.2e}')
            return None
        for row in range(k):
            if row != col:
                factor = aug[row][col] / pivot
                for j in range(k + 1):
                    aug[row][j] -= factor * aug[col][j]

    beta = [aug[a][k] / aug[a][a] for a in range(k)]

    y_hat = [sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]
    ss_res = sum(r ** 2 for r in resid)
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    adj_r2 = 1 - (ss_res / (n - k)) / (ss_tot / (n - 1)) if n > k else 0

    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    if cluster_var:
        if isinstance(cluster_var, (list, tuple)) and len(cluster_var) == 2:
            c1, c2 = cluster_var
            clusters1 = {}
            clusters2 = {}
            clusters12 = {}
            for i, d in enumerate(data):
                k1 = d.get(c1, None)
                k2 = d.get(c2, None)
                clusters1.setdefault(k1, []).append(i)
                clusters2.setdefault(k2, []).append(i)
                clusters12.setdefault((k1, k2), []).append(i)
            S1, G1, _ = _cluster_cov(X, resid, clusters1)
            S2, G2, _ = _cluster_cov(X, resid, clusters2)
            S12, G12, _ = _cluster_cov(X, resid, clusters12)
            # V = invX * (S1 + S2 - S12) * invX
            S = [[S1[a][b] + S2[a][b] - S12[a][b] for b in range(k)] for a in range(k)]
            V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
            # Small-sample correction using min cluster count
            G = min(G1, G2)
            if G > 1:
                scale = (G / (G - 1)) * ((n - 1) / (n - k))
                for a in range(k):
                    for b in range(k):
                        V[a][b] *= scale
            se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
            clusters = (G1, G2)
        else:
            clusters = {}
            for i, d in enumerate(data):
                cid = d.get(cluster_var, None)
                clusters.setdefault(cid, []).append(i)
            S, G, _ = _cluster_cov(X, resid, clusters)
            V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
            if G > 1:
                scale = (G / (G - 1)) * ((n - 1) / (n - k))
                for a in range(k):
                    for b in range(k):
                        V[a][b] *= scale
            se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
    else:
        s2 = ss_res / (n - k) if n > k else 0
        se = [math.sqrt(s2 * inv_XtX[a][a]) if inv_XtX[a][a] > 0 else 0.0 for a in range(k)]

    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0 for a in range(k)]

    names = ['intercept'] + x_vars
    return {
        'beta': dict(zip(names, beta)),
        'se': dict(zip(names, se)),
        't': dict(zip(names, t_stats)),
        'r2': r2,
        'adj_r2': adj_r2,
        'n': n,
        'clusters': clusters if cluster_var else None,
    }


# ── Load data ────────────────────────────────────────────────────────

# Daily returns
print('Loading daily returns...')
daily_ret = defaultdict(dict)
with open(derived_path('returns', 'daily_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        daily_ret[row['gvkey']][row['datadate']] = float(row['ret_daily'])
print(f'  Daily: {len(daily_ret)} firms')

# Market return (daily): prefer Fama-French vwretd; fallback to equal-weighted
def load_ff_factors_daily(path):
    if not os.path.exists(path):
        return None
    mktrf = {}
    rf = {}
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
            if not date.isdigit() or len(date) != 8:
                continue
            try:
                mktrf_val = float(parts[1])
                rf_val = float(parts[4])
            except ValueError:
                continue
            mktrf_dec = mktrf_val / 100.0
            rf_dec = rf_val / 100.0
            vwretd_dec = mktrf_dec + rf_dec
            date_fmt = f'{date[:4]}-{date[4:6]}-{date[6:]}'
            mktrf[date_fmt] = mktrf_dec
            rf[date_fmt] = rf_dec
            vwretd[date_fmt] = vwretd_dec
    return (mktrf, rf, vwretd) if mktrf else None


ff_daily_path = raw_path('factors', 'F-F_Research_Data_Factors_daily.csv')
ff_daily = load_ff_factors_daily(ff_daily_path)
if ff_daily:
    mktrf_daily, rf_daily, market_ret_daily = ff_daily
    print(f'  Market daily dates (F-F vwretd): {len(market_ret_daily)}')
else:
    _mkt_sum = defaultdict(float)
    _mkt_cnt = defaultdict(int)
    for gk, series in daily_ret.items():
        for d, r in series.items():
            _mkt_sum[d] += r
            _mkt_cnt[d] += 1
    market_ret_daily = {d: (_mkt_sum[d] / _mkt_cnt[d]) for d in _mkt_sum if _mkt_cnt[d] > 0}
    mktrf_daily = {}
    rf_daily = {}
    print(f'  Market daily dates (equal-weighted): {len(market_ret_daily)}')

# Monthly returns
print('Loading monthly returns...')
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        monthly_ret[row['gvkey']][row['datadate'][:7]] = float(row['ret_monthly'])
print(f'  Monthly: {len(monthly_ret)} firms')

# Market return (monthly): prefer Fama-French vwretd; fallback to equal-weighted
def load_ff_factors_monthly(path):
    if not os.path.exists(path):
        return None
    mktrf = {}
    rf = {}
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
            mktrf_dec = mktrf_val / 100.0
            rf_dec = rf_val / 100.0
            vwretd_dec = mktrf_dec + rf_dec
            date_fmt = f'{date[:4]}-{date[4:6]}'
            mktrf[date_fmt] = mktrf_dec
            rf[date_fmt] = rf_dec
            vwretd[date_fmt] = vwretd_dec
    return (mktrf, rf, vwretd) if mktrf else None


ff_monthly_path = raw_path('factors', 'F-F_Research_Data_Factors.csv')
ff_monthly = load_ff_factors_monthly(ff_monthly_path)
if ff_monthly:
    mktrf_monthly, rf_monthly, market_ret_monthly = ff_monthly
    print(f'  Market monthly months (F-F vwretd): {len(market_ret_monthly)}')
else:
    _mkt_sum = defaultdict(float)
    _mkt_cnt = defaultdict(int)
    for gk, series in monthly_ret.items():
        for m, r in series.items():
            _mkt_sum[m] += r
            _mkt_cnt[m] += 1
    market_ret_monthly = {m: (_mkt_sum[m] / _mkt_cnt[m]) for m in _mkt_sum if _mkt_cnt[m] > 0}
    mktrf_monthly = {}
    rf_monthly = {}
    print(f'  Market monthly months (equal-weighted): {len(market_ret_monthly)}')

# If FF factors missing, drop CAPM
if not mktrf_daily or not rf_daily or not mktrf_monthly or not rf_monthly:
    RETURN_MODELS = [m for m in RETURN_MODELS if m != 'capm']

# Weight matrices (geo primary; optional reg/fuel components)
print('Loading weight matrix...')
W = defaultdict(dict)
geo_path = derived_path('networks', 'weight_matrix_W_geo.csv')
with open(geo_path, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
print(f'  Firms in W: {len(W)}')

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
    print(f'  Regulatory W: {sum(len(v) for v in W_reg.values())} edges')
else:
    print('  Regulatory W: NOT FOUND (weight_matrix_W_regulatory.csv)')

W_fuel = defaultdict(dict)
fuel_path = derived_path('networks', 'weight_matrix_W_fuel.csv')
if os.path.exists(fuel_path):
    with open(fuel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            W_fuel[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])
    print(f'  Fuel W: {sum(len(v) for v in W_fuel.values())} edges')
else:
    print('  Fuel W: NOT FOUND (weight_matrix_W_fuel.csv). Run build_fuel_matrix.py')

# Fundamentals (panel)
fundamentals = {}
fundamentals_by_year = defaultdict(dict)
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        fundamentals_by_year[gk][fy] = row
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

# Winsorize kappa (interest coverage) at 1st/99th percentiles
kappa_vals = []
for rows in fundamentals_by_year.values():
    for row in rows.values():
        try:
            kap = float(row.get('kappa')) if row.get('kappa') not in (None, '') else None
        except (ValueError, TypeError):
            kap = None
        if kap is not None:
            kappa_vals.append(kap)
if kappa_vals:
    kappa_vals.sort()
    n_k = len(kappa_vals)
    kappa_p1 = kappa_vals[int(0.01 * (n_k - 1))]
    kappa_p99 = kappa_vals[int(0.99 * (n_k - 1))]
    for rows in fundamentals_by_year.values():
        for row in rows.values():
            if row.get('kappa') in (None, ''):
                continue
            try:
                kap = float(row['kappa'])
            except (ValueError, TypeError):
                continue
            kap = min(max(kap, kappa_p1), kappa_p99)
            row['kappa'] = f'{kap:.6f}'
    print(f'  Winsorized kappa at 1/99 pct: {kappa_p1:.3f}, {kappa_p99:.3f}')

# Time-varying alpha (panel)
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

# Alpha trajectory (pre-event slope of log CO2/Revenue)
alpha_traj = defaultdict(dict)  # gvkey -> year -> slope_log
traj_path = derived_path('trajectories', 'alpha_trajectory_panel.csv')
if os.path.exists(traj_path):
    with open(traj_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row.get('gvkey')
            if not gk:
                continue
            try:
                year = int(row.get('year'))
            except (ValueError, TypeError):
                continue
            try:
                slope = float(row.get('slope_log')) if row.get('slope_log') not in ('', None) else None
            except (ValueError, TypeError):
                slope = None
            if slope is None:
                continue
            alpha_traj[gk][year] = slope
    print(f'  Alpha trajectories: {sum(len(v) for v in alpha_traj.values())} firm-years '
          f'({len(alpha_traj)} firms)')
else:
    print('  Alpha trajectories: NOT FOUND (alpha_trajectory_panel.csv)')

# Classification thresholds (use all firm-years)
_lambdas_all = []
_kappas_all = []
for gk, rows in fundamentals_by_year.items():
    for fy, row in rows.items():
        try:
            lam = float(row['lambda']) if row.get('lambda') else None
            kap = float(row['kappa']) if row.get('kappa') else None
            if lam is not None:
                _lambdas_all.append(lam)
            if kap is not None:
                _kappas_all.append(kap)
        except (ValueError, TypeError):
            continue

_lambdas_all.sort()
_kappas_all.sort()
MEDIAN_LAMBDA = _lambdas_all[len(_lambdas_all) // 2] if _lambdas_all else 0.5
MEDIAN_KAPPA = _kappas_all[len(_kappas_all) // 2] if _kappas_all else 1.0
MAX_LAMBDA = max(_lambdas_all) if _lambdas_all else 1.0
MAX_KAPPA = max(_kappas_all) if _kappas_all else 1.0

# Alpha distribution for rank/kernel conditionality
alpha_vals = []
for gk, rows in fundamentals_by_year.items():
    for fy, row in rows.items():
        aval = row.get('alpha', None)
        if aval not in (None, ''):
            try:
                alpha_vals.append(float(aval))
            except ValueError:
                continue
alpha_vals.sort()
if alpha_vals:
    p25 = alpha_vals[len(alpha_vals) // 4]
    p75 = alpha_vals[(3 * len(alpha_vals)) // 4]
    iqr = p75 - p25
else:
    iqr = 0.2
KERNEL_SIGMA = max(0.05, iqr / 2.0)


def get_fundamentals_for_year(gvkey, year):
    rows = fundamentals_by_year.get(gvkey, {})
    if not rows:
        return fundamentals.get(gvkey)
    # pick latest year <= event year, else latest available
    years = [int(y) for y in rows.keys() if str(y).isdigit()]
    if not years:
        return fundamentals.get(gvkey)
    years_le = [y for y in years if y <= year]
    chosen = max(years_le) if years_le else max(years)
    return rows.get(str(chosen))


def get_alpha_for_year(gvkey, year):
    # prefer time-varying alpha panel if available
    if gvkey in alpha_panel:
        years = [int(y) for y in alpha_panel[gvkey].keys() if str(y).isdigit()]
        if years:
            years_le = [y for y in years if y <= year]
            chosen = max(years_le) if years_le else max(years)
            aval = alpha_panel[gvkey].get(str(chosen))
            if aval not in (None, ''):
                try:
                    return float(aval)
                except ValueError:
                    pass
    # fallback to fundamentals alpha
    f = get_fundamentals_for_year(gvkey, year)
    if f and f.get('alpha', '') != '':
        try:
            return float(f['alpha'])
        except ValueError:
            return None
    return None


def classify_firm(gvkey, mode='strict', event_year=None):
    f = fundamentals.get(gvkey)
    if event_year is None and f:
        event_year = int(f['fyear'])
    if event_year is None:
        return 'unknown'

    alpha = get_alpha_for_year(gvkey, event_year)
    if alpha is None:
        return 'unknown'
    if alpha < 0.3:
        return 'always_transform'
    elif alpha > 0.7:
        return 'never_transform'
    else:
        if mode == 'alpha_only':
            return 'conditional'
        f = get_fundamentals_for_year(gvkey, event_year)
        try:
            lam = float(f['lambda']) if f and f.get('lambda') else None
            kap = float(f['kappa']) if f and f.get('kappa') else None
        except (ValueError, TypeError):
            return 'conditional'
        if lam is not None and kap is not None:
            if lam <= MEDIAN_LAMBDA and kap >= MEDIAN_KAPPA:
                return 'conditional'
            else:
                return 'constrained'
        return 'conditional'


def conditional_score(alpha):
    """Continuous conditionality score in [0,1], peaks at alpha=0.5."""
    if alpha is None:
        return None
    score = 1.0 - abs(2.0 * alpha - 1.0)
    if score < 0:
        return 0.0
    if score > 1:
        return 1.0
    return score


def conditional_score_fin(alpha, lam, kap):
    """Financially-adjusted conditionality score in [0,1]."""
    base = conditional_score(alpha)
    if base is None:
        return None
    if lam is None or kap is None or MAX_LAMBDA <= 0 or MAX_KAPPA <= 0:
        return base
    lam_term = 1.0 - (lam / MAX_LAMBDA)
    kap_term = kap / MAX_KAPPA
    lam_term = max(0.0, min(1.0, lam_term))
    kap_term = max(0.0, min(1.0, kap_term))
    return base * lam_term * kap_term


def conditional_score_kernel(alpha, sigma):
    """Kernel-smoothed conditionality centered at 0.5."""
    if alpha is None:
        return None
    if sigma <= 0:
        return conditional_score(alpha)
    x = alpha - 0.5
    return math.exp(-(x * x) / (2.0 * sigma * sigma))


def conditional_score_rank(alpha, alpha_sorted):
    """Rank-based conditionality using percentile rank of alpha."""
    if alpha is None or not alpha_sorted:
        return None
    import bisect
    n = len(alpha_sorted)
    if n <= 1:
        return conditional_score(alpha)
    r = bisect.bisect_right(alpha_sorted, alpha)
    p = (r - 1) / (n - 1)
    return 1.0 - abs(2.0 * p - 1.0)


def get_sic4(gvkey):
    """Return 4-digit SIC for finer sector classification."""
    f = fundamentals.get(gvkey)
    if f and f.get('sic'):
        return f['sic'][:4]
    return None


def get_alpha_slope_for_year(gvkey, year):
    """Return pre-event slope of log(CO2/Revenue) for the event year."""
    if gvkey not in alpha_traj or year is None:
        return None
    years = [int(y) for y in alpha_traj[gvkey].keys()]
    if not years:
        return None
    years_le = [y for y in years if y <= year]
    chosen = max(years_le) if years_le else max(years)
    try:
        return float(alpha_traj[gvkey].get(chosen))
    except (ValueError, TypeError):
        return None


def is_exact_source(src):
    if not src:
        return False
    s = src.lower()
    if 'proxy' in s or 'approx' in s or 'mid' in s or 'month' in s:
        return False
    return True


# Events (load once, filter by scope later)
all_events = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        # Prefer announcement_date (when market learned) over event_date (physical closure)
        ann_date = row.get('announcement_date', '').strip()
        ret_date = row.get('event_date', '').strip()
        ann_src = row.get('announcement_source', '').strip()
        if EXACT_ONLY:
            if not ann_date:
                continue
            if ann_src and not is_exact_source(ann_src):
                continue
        effective_date = ann_date if ann_date else ret_date
        event_year = None
        if effective_date and len(effective_date) >= 4 and effective_date[:4].isdigit():
            event_year = int(effective_date[:4])
        else:
            event_year = int(row['ret_year']) if row.get('ret_year') else None
        # Lead time in days (announcement -> retirement), if both are day-exact
        lead_days = None
        if ann_date and ret_date and len(ann_date) >= 10 and len(ret_date) >= 10:
            try:
                ann_dt = datetime.date.fromisoformat(ann_date[:10])
                ret_dt = datetime.date.fromisoformat(ret_date[:10])
                lead_days = (ret_dt - ann_dt).days
            except ValueError:
                lead_days = None
        all_events.append({
            'plant': row['plant_name'],
            'year': event_year,
            'event_date': effective_date,
            'date_type': 'announcement' if ann_date else 'retirement',
            'announcement_source': ann_src,
            'gvkeys': row['matched_gvkeys'].split(';'),
            'is_first_mover': row.get('is_first_mover') == 'True',
            'lead_days': lead_days,
        })
print(f'All matched events: {len(all_events)}')


# ── CAR computation ──────────────────────────────────────────────────

def estimate_beta_daily(gvkey, dates, event_idx, pre=250):
    if not mktrf_daily or not rf_daily:
        return None
    pre_start = max(0, event_idx - pre)
    pre_end = event_idx
    xs = []
    ys = []
    for i in range(pre_start, pre_end):
        d = dates[i]
        if d in daily_ret[gvkey] and d in mktrf_daily and d in rf_daily:
            xs.append(mktrf_daily[d])
            ys.append(daily_ret[gvkey][d] - rf_daily[d])
    if len(xs) < 60:
        return None
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    denom = sum((x - mean_x) ** 2 for x in xs)
    if denom < 1e-12:
        return None
    beta = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / denom
    return beta


def compute_daily_car(gvkey, event_year, event_date='', pre=250, post=20, model='vwretd', pre_demean=False):
    """Daily CAR[-1, +20]. model in {'vwretd','capm','constant_mean'}."""
    if gvkey not in daily_ret:
        return None
    dates = sorted(daily_ret[gvkey].keys())
    event_idx = None
    if event_date and len(event_date) >= 10:
        ed = event_date[:10]
        for i, d in enumerate(dates):
            if d >= ed:
                event_idx = i
                break
    if event_idx is None:
        yr_s = str(event_year)
        for i, d in enumerate(dates):
            if d.startswith(yr_s):
                event_idx = i
                break
    if event_idx is None:
        return None

    pre_rets = [daily_ret[gvkey][dates[i]] for i in range(max(0, event_idx - pre), event_idx)
                if dates[i] in daily_ret[gvkey]]
    if len(pre_rets) < 60:
        return None

    mu = None
    beta = None
    if model == 'constant_mean':
        mu = sum(pre_rets) / len(pre_rets)
    elif model == 'capm':
        beta = estimate_beta_daily(gvkey, dates, event_idx, pre=pre)
        if beta is None:
            return None
    # Pre-demean ARs by pre-window mean for vwretd/capm
    pre_mean_ar = 0.0
    if pre_demean and model in ('vwretd', 'capm'):
        ar_list = []
        for i in range(max(0, event_idx - pre), event_idx):
            d = dates[i]
            if d not in daily_ret[gvkey]:
                continue
            r_it = daily_ret[gvkey][d]
            if model == 'capm':
                if d in mktrf_daily and d in rf_daily:
                    ar_list.append(r_it - (rf_daily[d] + beta * mktrf_daily[d]))
            else:
                if d in market_ret_daily:
                    ar_list.append(r_it - market_ret_daily[d])
        if ar_list:
            pre_mean_ar = sum(ar_list) / len(ar_list)

    car = 0.0
    for offset in range(-1, post + 1):
        idx = event_idx + offset
        if 0 <= idx < len(dates) and dates[idx] in daily_ret[gvkey]:
            d = dates[idx]
            r_it = daily_ret[gvkey][d]
            if model == 'constant_mean':
                car += r_it - mu
            elif model == 'capm':
                if d not in mktrf_daily or d not in rf_daily:
                    continue
                ar = r_it - (rf_daily[d] + beta * mktrf_daily[d])
                car += ar - pre_mean_ar if pre_demean else ar
            else:  # vwretd
                if d in market_ret_daily:
                    ar = r_it - market_ret_daily[d]
                    car += ar - pre_mean_ar if pre_demean else ar
    return car


def estimate_beta_monthly(gvkey, months, event_idx, pre=24):
    if not mktrf_monthly or not rf_monthly:
        return None
    pre_start = max(0, event_idx - pre)
    pre_end = event_idx
    xs = []
    ys = []
    for i in range(pre_start, pre_end):
        m = months[i]
        if m in monthly_ret[gvkey] and m in mktrf_monthly and m in rf_monthly:
            xs.append(mktrf_monthly[m])
            ys.append(monthly_ret[gvkey][m] - rf_monthly[m])
    if len(xs) < 12:
        return None
    mean_x = sum(xs) / len(xs)
    mean_y = sum(ys) / len(ys)
    denom = sum((x - mean_x) ** 2 for x in xs)
    if denom < 1e-12:
        return None
    beta = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys)) / denom
    return beta


def compute_monthly_car(gvkey, event_month, pre=24, post=12, model='vwretd', pre_demean=False):
    """Monthly CAR[-1, +k]. model in {'vwretd','capm','constant_mean'}."""
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

    pre_rets = [monthly_ret[gvkey][months[i]] for i in range(max(0, event_idx - pre), event_idx)
                if months[i] in monthly_ret[gvkey]]
    if len(pre_rets) < 12:
        return None

    mu = None
    beta = None
    if model == 'constant_mean':
        mu = sum(pre_rets) / len(pre_rets)
    elif model == 'capm':
        beta = estimate_beta_monthly(gvkey, months, event_idx, pre=pre)
        if beta is None:
            return None
    # Pre-demean ARs by pre-window mean for vwretd/capm
    pre_mean_ar = 0.0
    if pre_demean and model in ('vwretd', 'capm'):
        ar_list = []
        for i in range(max(0, event_idx - pre), event_idx):
            m = months[i]
            if m not in monthly_ret[gvkey]:
                continue
            r_it = monthly_ret[gvkey][m]
            if model == 'capm':
                if m in mktrf_monthly and m in rf_monthly:
                    ar_list.append(r_it - (rf_monthly[m] + beta * mktrf_monthly[m]))
            else:
                if m in market_ret_monthly:
                    ar_list.append(r_it - market_ret_monthly[m])
        if ar_list:
            pre_mean_ar = sum(ar_list) / len(ar_list)

    car = 0.0
    for offset in range(-1, post + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months) and months[idx] in monthly_ret[gvkey]:
            m = months[idx]
            r_it = monthly_ret[gvkey][m]
            if model == 'constant_mean':
                car += r_it - mu
            elif model == 'capm':
                if m not in mktrf_monthly or m not in rf_monthly:
                    continue
                ar = r_it - (rf_monthly[m] + beta * mktrf_monthly[m])
                car += ar - pre_mean_ar if pre_demean else ar
            else:  # vwretd
                if m in market_ret_monthly:
                    ar = r_it - market_ret_monthly[m]
                    car += ar - pre_mean_ar if pre_demean else ar
    return car


# -- Build regression dataset --


def build_regression_dataset(mode='strict', ret_model='vwretd', month_posts=None, events=None):
    print()
    print(f'=== BUILDING REGRESSION DATASET ({mode}, {ret_model}) ===')

    # For each event i, for each firm j (neighbor + non-connected controls),
    # one observation: (CAR_j, w_ij, ConditionalTransformer_j, SameSector_j)
    daily_obs = []
    if month_posts is None:
        month_posts = [12]
    monthly_obs_by_post = {p: [] for p in month_posts}

    if events is None:
        events = []

    for event_id, event in enumerate(events):
        event_gvkeys = set(event['gvkeys'])
        year = event['year']
        event_date = event.get('event_date', '')
        lead_time_z = event.get('lead_time_z')
        if event_date and len(event_date) >= 7:
            event_month = event_date[:7]
        else:
            event_month = f'{year}-07' if year else None

        # Get SIC4 of first-mover
        fm_sic4 = None
        for gk in event_gvkeys:
            fm_sic4 = get_sic4(gk)
            if fm_sic4:
                break

        for fm_gk in event_gvkeys:
            if fm_gk not in W:
                continue
            neighbors = W[fm_gk]

            # Neighbors + matched non-connected controls (5x neighbors to keep balanced)
            neighbor_gks = set(neighbors.keys()) - event_gvkeys
            non_connected = [gk for gk in fundamentals if gk not in event_gvkeys and gk not in neighbors]
            import random as _rng
            stable_seed = int(hashlib.md5(str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
            _rng.seed(stable_seed)
            n_ctrl = min(len(non_connected), max(5 * len(neighbor_gks), 20))
            ctrl_sample = _rng.sample(non_connected, n_ctrl) if len(non_connected) > n_ctrl else non_connected
            candidate_firms = list(neighbor_gks) + ctrl_sample

            for gk in candidate_firms:
                w_ij = neighbors.get(gk, 0.0)
                w_geo = w_ij
                w_reg = W_reg.get(fm_gk, {}).get(gk, 0.0)
                w_fuel = W_fuel.get(fm_gk, {}).get(gk, 0.0)
                ftype = classify_firm(gk, mode=mode, event_year=year)
                is_conditional = 1.0 if ftype == 'conditional' else 0.0
                alpha_j = get_alpha_for_year(gk, year) if year else None
                frow = get_fundamentals_for_year(gk, year) if year else fundamentals.get(gk)
                lam_j = None
                kap_j = None
                if frow:
                    try:
                        lam_j = float(frow['lambda']) if frow.get('lambda') else None
                        kap_j = float(frow['kappa']) if frow.get('kappa') else None
                    except (ValueError, TypeError):
                        lam_j = None
                        kap_j = None
                cond_score = conditional_score(alpha_j)
                cond_score_rank = conditional_score_rank(alpha_j, alpha_vals)
                cond_score_kernel = conditional_score_kernel(alpha_j, KERNEL_SIGMA)
                cond_score_fin = conditional_score_fin(alpha_j, lam_j, kap_j)
                alpha_slope = get_alpha_slope_for_year(gk, year) if year else None
                j_sic4 = get_sic4(gk)
                fic = frow.get('fic') if frow else None
                fe_cy = f'{fic}_{year}' if (fic and year) else None
                fe_ty = f'{j_sic4}_{year}' if (j_sic4 and year) else None
                same_sector = 1.0 if (fm_sic4 and j_sic4 and fm_sic4 == j_sic4) else 0.0

                # Daily CAR
                car_d = compute_daily_car(
                    gk, year, event_date, model=ret_model, pre_demean=PRE_DEMEAN_DAILY
                ) if year else None
                if car_d is not None:
                    daily_obs.append({
                        'car': car_d,
                        'w_ij': w_ij,
                        'w_geo': w_geo,
                        'w_reg': w_reg,
                        'w_fuel': w_fuel,
                        'w_x_cond': w_ij * is_conditional,
                        'w_x_lead_time': (w_ij * lead_time_z) if lead_time_z is not None else None,
                        'w_x_condscore': w_ij * (cond_score if cond_score is not None else 0.0),
                        'w_x_condscore_rank': w_ij * (cond_score_rank if cond_score_rank is not None else 0.0),
                        'w_x_condscore_kernel': w_ij * (cond_score_kernel if cond_score_kernel is not None else 0.0),
                        'w_x_condscore_fin': w_ij * (cond_score_fin if cond_score_fin is not None else 0.0),
                        'w_geo_x_cond': w_geo * is_conditional,
                        'w_reg_x_cond': w_reg * is_conditional,
                        'w_fuel_x_cond': w_fuel * is_conditional,
                        'w_geo_x_condscore': w_geo * (cond_score if cond_score is not None else 0.0),
                        'w_reg_x_condscore': w_reg * (cond_score if cond_score is not None else 0.0),
                        'w_fuel_x_condscore': w_fuel * (cond_score if cond_score is not None else 0.0),
                        'same_sector': same_sector,
                        'conditional': is_conditional,
                        'cond_score': cond_score if cond_score is not None else 0.0,
                        'cond_score_rank': cond_score_rank if cond_score_rank is not None else 0.0,
                        'cond_score_kernel': cond_score_kernel if cond_score_kernel is not None else 0.0,
                        'cond_score_fin': cond_score_fin if cond_score_fin is not None else 0.0,
                        'alpha_slope': alpha_slope,
                        'lead_time_z': lead_time_z,
                        'is_neighbor': 1.0 if w_ij > 0 else 0.0,
                        'event_id': event_id,
                        'gvkey': gk,
                        'fe_cy': fe_cy,
                        'fe_ty': fe_ty,
                    })

                # Monthly CARs for multiple windows
                for post in month_posts:
                    car_m = compute_monthly_car(
                        gk, event_month, post=post, model=ret_model, pre_demean=PRE_DEMEAN_MONTHLY
                    ) if event_month else None
                    if car_m is not None:
                        monthly_obs_by_post[post].append({
                            'car': car_m,
                            'w_ij': w_ij,
                            'w_geo': w_geo,
                            'w_reg': w_reg,
                            'w_fuel': w_fuel,
                            'w_x_cond': w_ij * is_conditional,
                            'w_x_lead_time': (w_ij * lead_time_z) if lead_time_z is not None else None,
                            'w_x_condscore': w_ij * (cond_score if cond_score is not None else 0.0),
                            'w_x_condscore_rank': w_ij * (cond_score_rank if cond_score_rank is not None else 0.0),
                            'w_x_condscore_kernel': w_ij * (cond_score_kernel if cond_score_kernel is not None else 0.0),
                            'w_x_condscore_fin': w_ij * (cond_score_fin if cond_score_fin is not None else 0.0),
                            'w_geo_x_cond': w_geo * is_conditional,
                            'w_reg_x_cond': w_reg * is_conditional,
                            'w_fuel_x_cond': w_fuel * is_conditional,
                            'w_geo_x_condscore': w_geo * (cond_score if cond_score is not None else 0.0),
                            'w_reg_x_condscore': w_reg * (cond_score if cond_score is not None else 0.0),
                            'w_fuel_x_condscore': w_fuel * (cond_score if cond_score is not None else 0.0),
                            'same_sector': same_sector,
                            'conditional': is_conditional,
                            'cond_score': cond_score if cond_score is not None else 0.0,
                            'cond_score_rank': cond_score_rank if cond_score_rank is not None else 0.0,
                            'cond_score_kernel': cond_score_kernel if cond_score_kernel is not None else 0.0,
                            'cond_score_fin': cond_score_fin if cond_score_fin is not None else 0.0,
                            'alpha_slope': alpha_slope,
                            'lead_time_z': lead_time_z,
                            'is_neighbor': 1.0 if w_ij > 0 else 0.0,
                            'event_id': event_id,
                            'gvkey': gk,
                            'fe_cy': fe_cy,
                            'fe_ty': fe_ty,
                        })

    print(f'Daily regression obs: {len(daily_obs)}')
    print(f'  Neighbors (w>0): {sum(1 for o in daily_obs if o["w_ij"] > 0)}')
    print(f'  Non-connected: {sum(1 for o in daily_obs if o["w_ij"] == 0)}')

    daily_obs_slope = [o for o in daily_obs if o.get('alpha_slope') is not None]
    monthly_obs_slope_by_post = {
        post: [o for o in obs if o.get('alpha_slope') is not None]
        for post, obs in monthly_obs_by_post.items()
    }

    print(f'  Alpha slope obs (daily): {len(daily_obs_slope)}')

    daily_obs_lead = [o for o in daily_obs if o.get('lead_time_z') is not None]
    monthly_obs_lead_by_post = {
        post: [o for o in obs if o.get('lead_time_z') is not None]
        for post, obs in monthly_obs_by_post.items()
    }
    if RUN_LEAD_TIME:
        print(f'  Lead-time obs (daily): {len(daily_obs_lead)}')

    # Apply exposure transformations (if requested)
    apply_exposure_transform(daily_obs, TRANSFORM_SET)
    for post, obs in monthly_obs_by_post.items():
        apply_exposure_transform(obs, TRANSFORM_SET)
    if RUN_LEAD_TIME:
        apply_exposure_transform(daily_obs_lead, TRANSFORM_SET)
        for post, obs in monthly_obs_lead_by_post.items():
            apply_exposure_transform(obs, TRANSFORM_SET)

    if os.getenv('PRINT_EXPOSURE_STATS', '0') == '1' and 12 in monthly_obs_by_post:
        obs = monthly_obs_by_post[12]
        if obs:
            vals = [o.get('w_ij', 0.0) for o in obs]
            mean_all = sum(vals) / len(vals)
            var_all = sum((v - mean_all) ** 2 for v in vals) / (len(vals) - 1) if len(vals) > 1 else 0.0
            sd_all = math.sqrt(var_all) if var_all > 0 else 0.0

            exposed = [v for v in vals if v > 0]
            if exposed:
                exposed.sort()
                p25 = exposed[int(0.25 * (len(exposed) - 1))]
                p75 = exposed[int(0.75 * (len(exposed) - 1))]
                iqr_exposed = p75 - p25
            else:
                iqr_exposed = 0.0

            print('EXPOSURE STATS (monthly +12)')
            print(f'  N={len(vals)}, mean={mean_all:.6f}, sd={sd_all:.6f}')
            print(f'  exposed N={len(exposed)}, exposed IQR={iqr_exposed:.6f}')

    if WRITE_METRICS and 12 in monthly_obs_by_post:
        obs = monthly_obs_by_post[12]
        metrics = {
            'car': summarize([o.get('car') for o in obs]),
            'w_ij': summarize([o.get('w_ij') for o in obs]),
            'w_geo': summarize([o.get('w_geo') for o in obs]),
            'w_reg': summarize([o.get('w_reg') for o in obs]),
            'w_fuel': summarize([o.get('w_fuel') for o in obs]),
        }
        out_path = os.path.join('JEEM_submission_package', 'JEEM_outputs', 'metrics',
                                f'strategy2_spatial_metrics_{TRANSFORM_SET}_monthly12.md')
        write_metrics_md(out_path, metrics)

    def add_zscore(obs_list, key):
        vals = [o[key] for o in obs_list if key in o and o[key] is not None]
        if len(vals) < 2:
            return
        mean = sum(vals) / len(vals)
        var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1)
        sd = math.sqrt(var) if var > 1e-12 else 1.0
        zkey = f'{key}_z'
        wkey = f'w_x_{key}_z'
        for o in obs_list:
            v = o.get(key, None)
            if v is None:
                o[zkey] = 0.0
                o[wkey] = 0.0
            else:
                z = (v - mean) / sd
                o[zkey] = z
                o[wkey] = o.get('w_ij', 0.0) * z

    def add_maxnorm(obs_list, key):
        vals = [o[key] for o in obs_list if key in o and o[key] is not None]
        if not vals:
            return
        vmax = max(vals)
        if vmax <= 0:
            return
        mkey = f'{key}_max'
        wkey = f'w_x_{key}_max'
        for o in obs_list:
            v = o.get(key, None)
            if v is None:
                o[mkey] = 0.0
                o[wkey] = 0.0
            else:
                m = v / vmax
                o[mkey] = m
                o[wkey] = o.get('w_ij', 0.0) * m

    def zscore_field(obs_list, src_key, dst_key):
        """Z-score a single field in-place, storing result under dst_key."""
        vals = [o[src_key] for o in obs_list if src_key in o and o[src_key] is not None]
        if len(vals) < 2:
            return
        mean_v = sum(vals) / len(vals)
        var_v = sum((v - mean_v) ** 2 for v in vals) / (len(vals) - 1)
        sd_v = math.sqrt(var_v) if var_v > 1e-12 else 1.0
        for o in obs_list:
            v = o.get(src_key, None)
            o[dst_key] = ((v - mean_v) / sd_v) if v is not None else 0.0

    # Add z-scored CondScoreFin to stabilize scale
    add_zscore(daily_obs, 'cond_score_fin')
    for post in month_posts:
        add_zscore(monthly_obs_by_post[post], 'cond_score_fin')

    # Add z-scored alpha trajectory slope (pre-event)
    add_zscore(daily_obs_slope, 'alpha_slope')
    for post in month_posts:
        add_zscore(monthly_obs_slope_by_post[post], 'alpha_slope')

    # Add max-normalized CondScoreFin (divide by sample max)
    add_maxnorm(daily_obs, 'cond_score_fin')
    for post in month_posts:
        add_maxnorm(monthly_obs_by_post[post], 'cond_score_fin')

    # Z-score the interaction itself: z(w_ij x CondScoreFin)
    # Puts b2 in "CAR per 1 SD of the interaction" units
    zscore_field(daily_obs, 'w_x_condscore_fin', 'w_x_condscore_fin_iz')
    for post in month_posts:
        zscore_field(monthly_obs_by_post[post], 'w_x_condscore_fin', 'w_x_condscore_fin_iz')

    for post in month_posts:
        obs = monthly_obs_by_post[post]
        print(f'Monthly regression obs (post={post}): {len(obs)}')
        print(f'  Neighbors (w>0): {sum(1 for o in obs if o["w_ij"] > 0)}')
        print(f'  Non-connected: {sum(1 for o in obs if o["w_ij"] == 0)}')
        if monthly_obs_slope_by_post.get(post) is not None:
            print(f'  Alpha slope obs (post={post}): {len(monthly_obs_slope_by_post[post])}')

    return (
        daily_obs,
        monthly_obs_by_post,
        daily_obs_slope,
        monthly_obs_slope_by_post,
        daily_obs_lead,
        monthly_obs_lead_by_post,
    )


# -- Run regressions --

def residualize_by_fe(obs_list, y_key, x_keys, fe_keys, iters=6):
    """Residualize y and x by multiple fixed effects via alternating demeaning."""
    keys_all = [y_key] + x_keys
    # Initialize working copy
    work = []
    for o in obs_list:
        row = {k: o.get(k, 0.0) for k in keys_all}
        for fe in fe_keys:
            row[fe] = o.get(fe)
        row['event_id'] = o.get('event_id')
        row['gvkey'] = o.get('gvkey')
        work.append(row)

    for _ in range(iters):
        for fe in fe_keys:
            sums = defaultdict(lambda: defaultdict(float))
            counts = defaultdict(int)
            for r in work:
                key = r.get(fe)
                if key is None:
                    continue
                counts[key] += 1
                for k in keys_all:
                    sums[key][k] += r[k]
            for r in work:
                key = r.get(fe)
                if key is None or counts.get(key, 0) == 0:
                    continue
                for k in keys_all:
                    r[k] -= sums[key][k] / counts[key]

    # Build output with residualized variables
    out = []
    for r in work:
        d = {
            'event_id': r.get('event_id'),
            'gvkey': r.get('gvkey'),
        }
        d[f'{y_key}_resid'] = r[y_key]
        for k in x_keys:
            d[f'{k}_resid'] = r[k]
        out.append(d)
    return out

def print_regression(label, result, x_vars):
    print()
    print(label)
    if result is None:
        print('  Regression failed (insufficient data or singular)')
        return
    clusters = result.get('clusters')
    if clusters:
        if isinstance(clusters, tuple):
            print(f'  N={result["n"]}, clusters={clusters}, R2={result["r2"]:.4f}, Adj R2={result["adj_r2"]:.4f}')
            print('  SEs two-way clustered (event x firm)')
        else:
            cluster_n = len(clusters) if hasattr(clusters, '__len__') else clusters
            print(f'  N={result["n"]}, clusters={cluster_n}, R2={result["r2"]:.4f}, Adj R2={result["adj_r2"]:.4f}')
            print('  SEs clustered by event')
    else:
        print(f'  N={result["n"]}, R2={result["r2"]:.4f}, Adj R2={result["adj_r2"]:.4f}')
    print(f'  {"Variable":<20} {"Beta":>10} {"SE":>10} {"t":>8} {"Sig":>5}')
    print(f'  {"-"*55}')
    for v in ['intercept'] + x_vars:
        b = result['beta'].get(v, 0)
        se = result['se'].get(v, 0)
        t = result['t'].get(v, 0)
        sig = '***' if abs(t) > 2.576 else '**' if abs(t) > 1.96 else '*' if abs(t) > 1.645 else ''
        print(f'  {v:<20} {b:>10.6f} {se:>10.6f} {t:>8.2f} {sig:>5}')


for scope in EVENT_SCOPES:
    if scope == 'first_mover':
        events = [e for e in all_events if e['is_first_mover']]
    else:
        events = list(all_events)
    if RUN_LEAD_TIME:
        lead_vals = [e.get('lead_days') for e in events if e.get('lead_days') is not None and e.get('lead_days') > 0]
        if len(lead_vals) >= 5:
            mean_lead = sum(lead_vals) / len(lead_vals)
            var_lead = sum((v - mean_lead) ** 2 for v in lead_vals) / (len(lead_vals) - 1)
            sd_lead = math.sqrt(var_lead) if var_lead > 0 else None
        else:
            mean_lead = None
            sd_lead = None
        for e in events:
            ld = e.get('lead_days')
            if sd_lead and ld is not None and ld > 0:
                e['lead_time_z'] = (ld - mean_lead) / sd_lead
            else:
                e['lead_time_z'] = None
    n_ann = sum(1 for e in events if e['date_type'] == 'announcement')
    print()
    print('=' * 70)
    print(f'EVENT SCOPE: {scope} (n={len(events)}, announcements={n_ann})')
    print('=' * 70)

    for ret_model in RETURN_MODELS:
        if RUN_FE_ONLY and ret_model != 'vwretd':
            continue
        print()
        print('#' * 70)
        print(f'#  RETURN MODEL: {ret_model.upper()}')
        if ret_model == 'vwretd':
            print('#  AR_t = R_it - R_mt  (CRSP value-weighted market return)')
        elif ret_model == 'capm':
            print('#  AR_t = R_it - [RF + beta*(Mkt-RF)]  (beta from pre-event window)')
        else:
            print('#  AR_t = R_it - mean(R_i, pre-event window)')
        print('#' * 70)

        for mode in CLASS_MODES:
            if RUN_FE_ONLY and not (mode == 'alpha_only' and scope == 'first_mover' and ret_model == 'vwretd'):
                continue
            daily_obs, monthly_obs_by_post, daily_obs_slope, monthly_obs_slope_by_post, daily_obs_lead, monthly_obs_lead_by_post = build_regression_dataset(
                mode, ret_model=ret_model, month_posts=MONTH_POSTS, events=events
            )

            # Paper specification: CAR_j = a + b1*w_ij + b2*(w_ij x Conditional_j) + b3*SameSector_j
            # Check if same_sector has variation; if not, drop it
            ss_source = daily_obs if daily_obs else monthly_obs_by_post.get(MONTH_POSTS[0], [])
            ss_vals = set(o['same_sector'] for o in ss_source) if ss_source else set()
            if len(ss_vals) <= 1:
                print()
                print(f'WARNING: same_sector has no variation (all={ss_vals}). Dropping from spec.')
                print('  (All firms are in the same 4-digit SIC; using reduced specification)')
                spec_vars = ['w_ij', 'w_x_cond']
                spec_vars_score = ['w_ij', 'w_x_condscore']
                spec_vars_score_rank = ['w_ij', 'w_x_condscore_rank']
                spec_vars_score_kernel = ['w_ij', 'w_x_condscore_kernel']
                spec_vars_score_fin = ['w_ij', 'w_x_cond_score_fin_z']
                spec_vars_score_fin_max = ['w_ij', 'w_x_cond_score_fin_max']
                spec_vars_score_fin_iz = ['w_ij', 'w_x_condscore_fin_iz']
                spec_vars_slope = ['w_ij', 'w_x_alpha_slope_z']
                spec_vars_chan = ['w_geo', 'w_reg', 'w_fuel']
                spec_vars_chan_cond = [
                    'w_geo', 'w_reg', 'w_fuel',
                    'w_geo_x_condscore', 'w_reg_x_condscore', 'w_fuel_x_condscore',
                ]
            else:
                spec_vars = ['w_ij', 'w_x_cond', 'same_sector']
                spec_vars_score = ['w_ij', 'w_x_condscore', 'same_sector']
                spec_vars_score_rank = ['w_ij', 'w_x_condscore_rank', 'same_sector']
                spec_vars_score_kernel = ['w_ij', 'w_x_condscore_kernel', 'same_sector']
                spec_vars_score_fin = ['w_ij', 'w_x_cond_score_fin_z', 'same_sector']
                spec_vars_score_fin_max = ['w_ij', 'w_x_cond_score_fin_max', 'same_sector']
                spec_vars_score_fin_iz = ['w_ij', 'w_x_condscore_fin_iz', 'same_sector']
                spec_vars_slope = ['w_ij', 'w_x_alpha_slope_z', 'same_sector']
                spec_vars_bivar = ['w_ij']
                spec_vars_chan = ['w_geo', 'w_reg', 'w_fuel', 'same_sector']
                spec_vars_chan2 = ['w_geo', 'w_fuel', 'same_sector']
                spec_vars_chan_cond = [
                    'w_geo', 'w_reg', 'w_fuel',
                    'w_geo_x_condscore', 'w_reg_x_condscore', 'w_fuel_x_condscore',
                    'same_sector',
                ]
                print()
                print(f'Same-sector variation: {ss_vals}')

            print()
            print('=' * 60)
            print(f'STRATEGY 2 SPATIAL REGRESSION ({mode}, {ret_model})')
            print('CAR_j = a + b1*w_ij + b2*(w_ij x Cond_j) + b3*SameSector_j')
            print('=' * 60)

            # Daily
            res_daily = ols(daily_obs, 'car', spec_vars, cluster_var='event_id')
            print_regression('DAILY CARs [-1, +20]', res_daily, spec_vars)
            res_daily_tw = ols(daily_obs, 'car', spec_vars, cluster_var=['event_id', 'gvkey'])
            print_regression('DAILY CARs [-1, +20] (two-way)', res_daily_tw, spec_vars)

            # Channel decomposition (geo/reg/fuel) for primary spec
            if RUN_CHANNEL_DECOMP and mode == 'alpha_only' and ret_model == 'vwretd' and scope == 'first_mover':
                print('\nCHANNEL DECOMPOSITION (geo/reg/fuel)')
                res_chan_d = ols(daily_obs, 'car', spec_vars_chan, cluster_var='event_id')
                print_regression('DAILY (channels)', res_chan_d, spec_vars_chan)
                res_chan_d_tw = ols(daily_obs, 'car', spec_vars_chan, cluster_var=['event_id', 'gvkey'])
                print_regression('DAILY (channels, two-way)', res_chan_d_tw, spec_vars_chan)
                print('CHANNEL x CondScore: add interactions with CondScore')
                res_chan_cd = ols(daily_obs, 'car', spec_vars_chan_cond, cluster_var='event_id')
                print_regression('DAILY (channels x CondScore)', res_chan_cd, spec_vars_chan_cond)
                res_chan_cd_tw = ols(daily_obs, 'car', spec_vars_chan_cond, cluster_var=['event_id', 'gvkey'])
                print_regression('DAILY (channels x CondScore, two-way)', res_chan_cd_tw, spec_vars_chan_cond)

            # Monthly (multiple windows)
            res_monthly_by_post = {}
            res_monthly_tw_by_post = {}
            for post in MONTH_POSTS:
                monthly_obs = monthly_obs_by_post.get(post, [])
                res_monthly = ols(monthly_obs, 'car', spec_vars, cluster_var='event_id')
                res_monthly_by_post[post] = res_monthly
                print_regression(f'MONTHLY CARs [-1, +{post}]', res_monthly, spec_vars)
                if post == 12:
                    # Bivariate (exposure only) for transparency
                    res_monthly_bi = ols(monthly_obs, 'car', spec_vars_bivar, cluster_var='event_id')
                    print_regression(f'MONTHLY CARs [-1, +{post}] (bivariate)', res_monthly_bi, spec_vars_bivar)
                    res_monthly_tw = ols(monthly_obs, 'car', spec_vars, cluster_var=['event_id', 'gvkey'])
                    res_monthly_tw_by_post[post] = res_monthly_tw
                    print_regression(f'MONTHLY CARs [-1, +{post}] (two-way)', res_monthly_tw, spec_vars)
                    # Common-shock controls: country-year and technology-year FE (primary spec only)
                    if mode == 'alpha_only' and ret_model == 'vwretd' and scope == 'first_mover':
                        fe_x = ['w_ij']
                        if 'same_sector' in spec_vars:
                            fe_x.append('same_sector')
                        fe_obs = residualize_by_fe(monthly_obs, 'car', fe_x, ['fe_cy', 'fe_ty'])
                        fe_x_resid = [f'{k}_resid' for k in fe_x]
                        res_fe = ols(fe_obs, 'car_resid', fe_x_resid, cluster_var=['event_id', 'gvkey'])
                        print_regression('MONTHLY (country-year + tech-year FE)', res_fe, fe_x_resid)

                if RUN_CHANNEL_DECOMP and mode == 'alpha_only' and ret_model == 'vwretd' and scope == 'first_mover':
                    print('\nCHANNEL DECOMPOSITION (geo/reg/fuel): monthly')
                    res_chan_m = ols(monthly_obs, 'car', spec_vars_chan, cluster_var='event_id')
                    print_regression(f'MONTHLY (channels, +{post})', res_chan_m, spec_vars_chan)
                    # Two-channel (geo+fuel) robustness for short-horizon split
                    if post == 3:
                        res_chan_m2 = ols(monthly_obs, 'car', spec_vars_chan2, cluster_var='event_id')
                        print_regression('MONTHLY (channels: geo+fuel, +3)', res_chan_m2, spec_vars_chan2)
                    if post == 12:
                        res_chan_m_tw = ols(monthly_obs, 'car', spec_vars_chan, cluster_var=['event_id', 'gvkey'])
                        print_regression(f'MONTHLY (channels, +{post}, two-way)', res_chan_m_tw, spec_vars_chan)
                    res_chan_cm = ols(monthly_obs, 'car', spec_vars_chan_cond, cluster_var='event_id')
                    print_regression(f'MONTHLY (channels x CondScore, +{post})', res_chan_cm, spec_vars_chan_cond)
                    if post == 12:
                        res_chan_cm_tw = ols(monthly_obs, 'car', spec_vars_chan_cond, cluster_var=['event_id', 'gvkey'])
                        print_regression(f'MONTHLY (channels x CondScore, +{post}, two-way)', res_chan_cm_tw, spec_vars_chan_cond)

            # Continuous conditionality score (alpha-only)
            print()
            print('=' * 60)
            print(f'CONTINUOUS CONDSCORE (alpha-only) ({mode}, {ret_model})')
            print('CAR_j = a + b1*w_ij + b2*(w_ij x CondScore_j) + b3*SameSector_j')
            print('=' * 60)
            res_daily_score = ols(daily_obs, 'car', spec_vars_score, cluster_var='event_id')
            print_regression('DAILY (condscore)', res_daily_score, spec_vars_score)
            res_daily_score_tw = ols(daily_obs, 'car', spec_vars_score, cluster_var=['event_id', 'gvkey'])
            print_regression('DAILY (condscore, two-way)', res_daily_score_tw, spec_vars_score)
            res_monthly_score_by_post = {}
            for post in MONTH_POSTS:
                monthly_obs = monthly_obs_by_post.get(post, [])
                res_monthly_score = ols(monthly_obs, 'car', spec_vars_score, cluster_var='event_id')
                res_monthly_score_by_post[post] = res_monthly_score
                print_regression(f'MONTHLY (condscore, +{post})', res_monthly_score, spec_vars_score)
                if post == 12:
                    res_monthly_score_tw = ols(monthly_obs, 'car', spec_vars_score, cluster_var=['event_id', 'gvkey'])
                    print_regression(f'MONTHLY (condscore, +{post}, two-way)', res_monthly_score_tw, spec_vars_score)

            # Rank-based conditionality score
            print()
            print('=' * 60)
            print(f'RANK-BASED CONDSCORE ({mode}, {ret_model})')
            print('CAR_j = a + b1*w_ij + b2*(w_ij x CondScore_rank_j) + b3*SameSector_j')
            print('=' * 60)
            res_daily_score_rank = ols(daily_obs, 'car', spec_vars_score_rank, cluster_var='event_id')
            print_regression('DAILY (condscore_rank)', res_daily_score_rank, spec_vars_score_rank)
            res_daily_score_rank_tw = ols(daily_obs, 'car', spec_vars_score_rank, cluster_var=['event_id', 'gvkey'])
            print_regression('DAILY (condscore_rank, two-way)', res_daily_score_rank_tw, spec_vars_score_rank)
            for post in MONTH_POSTS:
                monthly_obs = monthly_obs_by_post.get(post, [])
                res_monthly_score_rank = ols(monthly_obs, 'car', spec_vars_score_rank, cluster_var='event_id')
                print_regression(f'MONTHLY (condscore_rank, +{post})', res_monthly_score_rank, spec_vars_score_rank)
                if post == 12:
                    res_monthly_score_rank_tw = ols(monthly_obs, 'car', spec_vars_score_rank, cluster_var=['event_id', 'gvkey'])
                    print_regression(f'MONTHLY (condscore_rank, +{post}, two-way)', res_monthly_score_rank_tw, spec_vars_score_rank)

            # Kernel-smoothed conditionality score
            print()
            print('=' * 60)
            print(f'KERNEL CONDSCORE (sigma={KERNEL_SIGMA:.3f}) ({mode}, {ret_model})')
            print('CAR_j = a + b1*w_ij + b2*(w_ij x CondScore_kernel_j) + b3*SameSector_j')
            print('=' * 60)
            res_daily_score_kernel = ols(daily_obs, 'car', spec_vars_score_kernel, cluster_var='event_id')
            print_regression('DAILY (condscore_kernel)', res_daily_score_kernel, spec_vars_score_kernel)
            res_daily_score_kernel_tw = ols(daily_obs, 'car', spec_vars_score_kernel, cluster_var=['event_id', 'gvkey'])
            print_regression('DAILY (condscore_kernel, two-way)', res_daily_score_kernel_tw, spec_vars_score_kernel)
            for post in MONTH_POSTS:
                monthly_obs = monthly_obs_by_post.get(post, [])
                res_monthly_score_kernel = ols(monthly_obs, 'car', spec_vars_score_kernel, cluster_var='event_id')
                print_regression(f'MONTHLY (condscore_kernel, +{post})', res_monthly_score_kernel, spec_vars_score_kernel)
                if post == 12:
                    res_monthly_score_kernel_tw = ols(monthly_obs, 'car', spec_vars_score_kernel, cluster_var=['event_id', 'gvkey'])
                    print_regression(f'MONTHLY (condscore_kernel, +{post}, two-way)', res_monthly_score_kernel_tw, spec_vars_score_kernel)

            # Continuous conditionality score (financially adjusted)
            print()
            print('=' * 60)
            print(f'CONTINUOUS CONDSCORE (alpha+financials, z-scored) ({mode}, {ret_model})')
            print('CAR_j = a + b1*w_ij + b2*(w_ij x z(CondScoreFin_j)) + b3*SameSector_j')
            print('=' * 60)
            res_daily_score_fin = ols(daily_obs, 'car', spec_vars_score_fin, cluster_var='event_id')
            print_regression('DAILY (condscore_fin)', res_daily_score_fin, spec_vars_score_fin)
            res_daily_score_fin_tw = ols(daily_obs, 'car', spec_vars_score_fin, cluster_var=['event_id', 'gvkey'])
            print_regression('DAILY (condscore_fin, two-way)', res_daily_score_fin_tw, spec_vars_score_fin)
            for post in MONTH_POSTS:
                monthly_obs = monthly_obs_by_post.get(post, [])
                res_monthly_score_fin = ols(monthly_obs, 'car', spec_vars_score_fin, cluster_var='event_id')
                print_regression(f'MONTHLY (condscore_fin, +{post})', res_monthly_score_fin, spec_vars_score_fin)
                if post == 12:
                    res_monthly_score_fin_tw = ols(monthly_obs, 'car', spec_vars_score_fin, cluster_var=['event_id', 'gvkey'])
                    print_regression(f'MONTHLY (condscore_fin, +{post}, two-way)', res_monthly_score_fin_tw, spec_vars_score_fin)

            # Continuous conditionality score (financially adjusted, max-normalized)
            print()
            print('=' * 60)
            print(f'CONTINUOUS CONDSCORE (alpha+financials, max-normalized) ({mode}, {ret_model})')
            print('CAR_j = a + b1*w_ij + b2*(w_ij x max(CondScoreFin_j)) + b3*SameSector_j')
            print('=' * 60)
            res_daily_score_fin_max = ols(daily_obs, 'car', spec_vars_score_fin_max, cluster_var='event_id')
            print_regression('DAILY (condscore_fin_max)', res_daily_score_fin_max, spec_vars_score_fin_max)
            res_daily_score_fin_max_tw = ols(daily_obs, 'car', spec_vars_score_fin_max, cluster_var=['event_id', 'gvkey'])
            print_regression('DAILY (condscore_fin_max, two-way)', res_daily_score_fin_max_tw, spec_vars_score_fin_max)
            for post in MONTH_POSTS:
                monthly_obs = monthly_obs_by_post.get(post, [])
                res_monthly_score_fin_max = ols(monthly_obs, 'car', spec_vars_score_fin_max, cluster_var='event_id')
                print_regression(f'MONTHLY (condscore_fin_max, +{post})', res_monthly_score_fin_max, spec_vars_score_fin_max)
                if post == 12:
                    res_monthly_score_fin_max_tw = ols(monthly_obs, 'car', spec_vars_score_fin_max, cluster_var=['event_id', 'gvkey'])
                    print_regression(f'MONTHLY (condscore_fin_max, +{post}, two-way)', res_monthly_score_fin_max_tw, spec_vars_score_fin_max)

            # Interaction z-scored: z(w_ij x CondScoreFin) — b2 in CAR per 1 SD units
            print()
            print('=' * 60)
            print(f'INTERACTION Z-SCORED: z(w_ij x CondScoreFin) ({mode}, {ret_model})')
            print('CAR_j = a + b1*w_ij + b2*z(w_ij x CondScoreFin_j) + b3*SameSector_j')
            print('=' * 60)
            res_daily_iz = ols(daily_obs, 'car', spec_vars_score_fin_iz, cluster_var='event_id')
            print_regression('DAILY (interaction z-scored)', res_daily_iz, spec_vars_score_fin_iz)
            res_daily_iz_tw = ols(daily_obs, 'car', spec_vars_score_fin_iz, cluster_var=['event_id', 'gvkey'])
            print_regression('DAILY (interaction z-scored, two-way)', res_daily_iz_tw, spec_vars_score_fin_iz)
            for post in MONTH_POSTS:
                monthly_obs = monthly_obs_by_post.get(post, [])
                res_monthly_iz = ols(monthly_obs, 'car', spec_vars_score_fin_iz, cluster_var='event_id')
                print_regression(f'MONTHLY (interaction z-scored, +{post})', res_monthly_iz, spec_vars_score_fin_iz)
                if post == 12:
                    res_monthly_iz_tw = ols(monthly_obs, 'car', spec_vars_score_fin_iz, cluster_var=['event_id', 'gvkey'])
                    print_regression(f'MONTHLY (interaction z-scored, +{post}, two-way)', res_monthly_iz_tw, spec_vars_score_fin_iz)

            # Alpha trajectory slope (pre-event CO2/Revenue trend)
            print()
            print('=' * 60)
            print(f'ALPHA TRAJECTORY (pre-event slope) ({mode}, {ret_model})')
            print('CAR_j = a + b1*w_ij + b2*(w_ij x z(alpha_slope_j)) + b3*SameSector_j')
            print('=' * 60)
            res_daily_slope = ols(daily_obs_slope, 'car', spec_vars_slope, cluster_var='event_id')
            print_regression('DAILY (alpha_slope)', res_daily_slope, spec_vars_slope)
            res_daily_slope_tw = ols(daily_obs_slope, 'car', spec_vars_slope, cluster_var=['event_id', 'gvkey'])
            print_regression('DAILY (alpha_slope, two-way)', res_daily_slope_tw, spec_vars_slope)
            for post in MONTH_POSTS:
                monthly_obs_slope = monthly_obs_slope_by_post.get(post, [])
                res_monthly_slope = ols(monthly_obs_slope, 'car', spec_vars_slope, cluster_var='event_id')
                print_regression(f'MONTHLY (alpha_slope, +{post})', res_monthly_slope, spec_vars_slope)
                if post == 12:
                    res_monthly_slope_tw = ols(monthly_obs_slope, 'car', spec_vars_slope, cluster_var=['event_id', 'gvkey'])
                    print_regression(f'MONTHLY (alpha_slope, +{post}, two-way)', res_monthly_slope_tw, spec_vars_slope)

            # Lead-time interaction (announcement -> retirement)
            if RUN_LEAD_TIME:
                spec_vars_lead = ['w_ij', 'w_x_lead_time']
                if 'same_sector' in spec_vars:
                    spec_vars_lead.append('same_sector')
                print()
                print('=' * 60)
                print(f'LEAD-TIME INTERACTION (announcement -> retirement) ({mode}, {ret_model})')
                print('CAR_j = a + b1*w_ij + b2*(w_ij x z(lead_time)) + b3*SameSector_j')
                print('=' * 60)
                res_daily_lead = ols(daily_obs_lead, 'car', spec_vars_lead, cluster_var='event_id')
                print_regression('DAILY (lead_time)', res_daily_lead, spec_vars_lead)
                res_daily_lead_tw = ols(daily_obs_lead, 'car', spec_vars_lead, cluster_var=['event_id', 'gvkey'])
                print_regression('DAILY (lead_time, two-way)', res_daily_lead_tw, spec_vars_lead)
                for post in MONTH_POSTS:
                    monthly_obs_lead = monthly_obs_lead_by_post.get(post, [])
                    res_monthly_lead = ols(monthly_obs_lead, 'car', spec_vars_lead, cluster_var='event_id')
                    print_regression(f'MONTHLY (lead_time, +{post})', res_monthly_lead, spec_vars_lead)
                    if post == 12:
                        res_monthly_lead_tw = ols(monthly_obs_lead, 'car', spec_vars_lead, cluster_var=['event_id', 'gvkey'])
                        print_regression(f'MONTHLY (lead_time, +{post}, two-way)', res_monthly_lead_tw, spec_vars_lead)

            # Also run extended spec with Conditional main effect (binary)
            spec_ext = ['w_ij', 'conditional', 'w_x_cond']
            if 'same_sector' in spec_vars:
                spec_ext.append('same_sector')
            print()
            print('=' * 60)
            print(f'EXTENDED SPEC: + Conditional main effect ({ret_model})')
            print('=' * 60)
            res_daily_ext = ols(daily_obs, 'car', spec_ext, cluster_var='event_id')
            print_regression('DAILY (extended)', res_daily_ext, spec_ext)
            res_daily_ext_tw = ols(daily_obs, 'car', spec_ext, cluster_var=['event_id', 'gvkey'])
            print_regression('DAILY (extended, two-way)', res_daily_ext_tw, spec_ext)
            for post in MONTH_POSTS:
                monthly_obs = monthly_obs_by_post.get(post, [])
                res_monthly_ext = ols(monthly_obs, 'car', spec_ext, cluster_var='event_id')
                print_regression(f'MONTHLY (extended, +{post})', res_monthly_ext, spec_ext)
                if post == 12:
                    res_monthly_ext_tw = ols(monthly_obs, 'car', spec_ext, cluster_var=['event_id', 'gvkey'])
                    print_regression(f'MONTHLY (extended, +{post}, two-way)', res_monthly_ext_tw, spec_ext)

            # Coordination signature interpretation (binary)
            print()
            print(f'=== COORDINATION SIGNATURE ({ret_model}) ===')
            if res_daily and res_daily['beta'].get('w_ij', 0) > 0:
                print('Daily: b1 (w_ij) > 0  -> spatial spillover EXISTS')
            else:
                print('Daily: b1 (w_ij) <= 0  -> no spatial spillover detected')
            if res_daily and res_daily['beta'].get('w_x_cond', 0) > 0:
                print('Daily: b2 (w x Cond) > 0 -> conditional transformers respond MORE (coordination)')
            else:
                print('Daily: b2 (w x Cond) <= 0 -> no extra response for conditional transformers')
            for post in MONTH_POSTS:
                res_monthly = res_monthly_by_post.get(post)
                if res_monthly and res_monthly['beta'].get('w_ij', 0) > 0:
                    print(f'Monthly(+{post}): b1 (w_ij) > 0 -> spatial spillover EXISTS')
                else:
                    print(f'Monthly(+{post}): b1 (w_ij) <= 0 -> no spatial spillover detected')
                if res_monthly and res_monthly['beta'].get('w_x_cond', 0) > 0:
                    print(f'Monthly(+{post}): b2 (w x Cond) > 0 -> conditional transformers respond MORE')
                else:
                    print(f'Monthly(+{post}): b2 (w x Cond) <= 0 -> no extra response for conditional transformers')

