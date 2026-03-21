"""Strategy 1: Kernel Sensitivity for Regime-Dependent R^2.

Compute the R^2 split using continuous geographic kernels with different half-lives
without re-building the full weight matrix. This avoids overwriting the baseline
network and keeps the analysis reproducible.

Half-life values (km): 250, 500, 1000
"""
import csv
import math
import random
from collections import defaultdict

from _paths import derived_path

HALF_LIVES = [250, 500, 1000]
N_BOOT = 1000
SEED = 42

def haversine(lat1, lon1, lat2, lon2):
    """Distance in km between two GPS points."""
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(a))

def ols_r2(data, y_var, x_vars):
    """Simple OLS R^2 using normal equations (no numpy)."""
    n = len(data)
    k = len(x_vars) + 1  # intercept
    if n <= k:
        return None, None
    y = [d[y_var] for d in data]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]
    # X'X
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    # X'y
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]
    # Solve via Gaussian elimination
    aug = [XtX[a] + [Xty[a]] for a in range(k)]
    for col in range(k):
        max_row = max(range(col, k), key=lambda r: abs(aug[r][col]))
        aug[col], aug[max_row] = aug[max_row], aug[col]
        if abs(aug[col][col]) < 1e-12:
            return None, None
        for row in range(k):
            if row != col:
                factor = aug[row][col] / aug[col][col]
                for j in range(k + 1):
                    aug[row][j] -= factor * aug[col][j]
    beta = [aug[a][k] / aug[a][a] for a in range(k)]
    y_hat = [sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
    ss_res = sum((y[i] - y_hat[i]) ** 2 for i in range(n))
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    adj_r2 = 1 - (ss_res / (n - k)) / (ss_tot / (n - 1)) if n > k and ss_tot > 0 else 0
    return r2, adj_r2

def bootstrap_p(high_d, low_d, x_vars, y_var='delivery', n_boot=1000):
    diffs = []
    for _ in range(n_boot):
        boot_high = random.choices(high_d, k=len(high_d))
        boot_low = random.choices(low_d, k=len(low_d))
        r2_h, _ = ols_r2(boot_high, y_var, x_vars)
        r2_l, _ = ols_r2(boot_low, y_var, x_vars)
        if r2_h is not None and r2_l is not None:
            diffs.append(r2_l - r2_h)
    diffs.sort()
    if not diffs:
        return None
    p_val = sum(1 for d in diffs if d <= 0) / len(diffs)
    return p_val

# 1) Load fundamentals (latest year per gvkey)
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

# Winsorize kappa at 1/99 pct
kappa_vals = []
for f in fundamentals.values():
    try:
        kap = float(f.get('kappa')) if f.get('kappa') not in (None, '') else None
    except (ValueError, TypeError):
        kap = None
    if kap is not None:
        kappa_vals.append(kap)
if kappa_vals:
    kappa_vals.sort()
    n_k = len(kappa_vals)
    kappa_p1 = kappa_vals[int(0.01 * (n_k - 1))]
    kappa_p99 = kappa_vals[int(0.99 * (n_k - 1))]
    for f in fundamentals.values():
        if f.get('kappa') in (None, ''):
            continue
        try:
            kap = float(f['kappa'])
        except (ValueError, TypeError):
            continue
        kap = min(max(kap, kappa_p1), kappa_p99)
        f['kappa'] = f'{kap:.6f}'

# 2) Load centroids (lat/lon)
centroids = {}
with open(derived_path('networks', 'firm_centroids.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        try:
            lat = float(row['centroid_lat'])
            lon = float(row['centroid_lon'])
        except (ValueError, TypeError):
            continue
        centroids[gk] = (lat, lon)

gvkeys = sorted(centroids.keys())
n = len(gvkeys)

# Pre-compute distance matrix
dist = [[0.0] * n for _ in range(n)]
for i in range(n):
    lat_i, lon_i = centroids[gvkeys[i]]
    for j in range(i + 1, n):
        lat_j, lon_j = centroids[gvkeys[j]]
        d = haversine(lat_i, lon_i, lat_j, lon_j)
        dist[i][j] = d
        dist[j][i] = d

# 3) Region mapping for controls (same as Strategy 1)
REGION_MAP = {}
for c in ['AT','BE','BG','CH','CZ','DE','DK','ES','FI','FR','GB','GR','HR','HU',
          'IE','IT','LT','LU','NL','NO','PL','PT','RO','SE','SI','SK','UA','RS','BA','MK','EE','LV','CY','MT']:
    REGION_MAP[c] = 'europe'
for c in ['US','CA']:
    REGION_MAP[c] = 'north_america'
for c in ['BR','CL','MX','CO','AR','PE','VE','EC','UY','PY','BO','PA','CR','DO','JM','TT','GT','HN','SV','NI']:
    REGION_MAP[c] = 'latin_america'
for c in ['CN','JP','KR','TW','HK','MO']:
    REGION_MAP[c] = 'east_asia'
for c in ['IN','PH','TH','MY','ID','SG','VN','BD','PK','LK','MM','KH','LA','NP']:
    REGION_MAP[c] = 'south_se_asia'
for c in ['AU','NZ','FJ','PG']:
    REGION_MAP[c] = 'oceania'

# 4) Build base sample (no density yet)
base_sample = {}
for gk, f in fundamentals.items():
    if gk not in centroids:
        continue
    theta = {}
    skip = False
    for var in ['alpha', 'lambda', 'rho', 'kappa', 'delta']:
        val = f.get(var, '')
        if val == '':
            skip = True
            break
        theta[var] = float(val)
    if skip:
        continue
    try:
        log_at = math.log(max(float(f['at']), 1.0)) if f.get('at') else None
    except (ValueError, TypeError):
        log_at = None
    if log_at is None:
        continue
    fic = f.get('fic', '')
    region = REGION_MAP.get(fic, 'other')
    base_sample[gk] = {
        'gvkey': gk,
        'density': None,
        'delivery': 1.0 - theta['alpha'],
        'log_at': log_at,
        'reg_europe': 1.0 if region == 'europe' else 0.0,
        'reg_north_am': 1.0 if region == 'north_america' else 0.0,
        'reg_latin_am': 1.0 if region == 'latin_america' else 0.0,
        'reg_east_asia': 1.0 if region == 'east_asia' else 0.0,
        'reg_south_se_asia': 1.0 if region == 'south_se_asia' else 0.0,
        'reg_oceania': 1.0 if region == 'oceania' else 0.0,
        **theta,
    }

print(f'Base sample: {len(base_sample)} firms with fundamentals + centroids')

random.seed(SEED)

results = []
for half_life in HALF_LIVES:
    decay = half_life / math.log(2)
    density = {}
    for i, gi in enumerate(gvkeys):
        if gi not in base_sample:
            continue
        row_sum = 0.0
        for j in range(n):
            if i == j:
                continue
            d = dist[i][j]
            if d > 0:
                row_sum += math.exp(-d / decay) / d
        density[gi] = row_sum

    # Build sample with density
    sample = []
    for gk, row in base_sample.items():
        if gk not in density:
            continue
        rec = dict(row)
        rec['density'] = density[gk]
        sample.append(rec)

    densities = sorted(s['density'] for s in sample)
    median_d = densities[len(densities) // 2]
    high_d = [s for s in sample if s['density'] > median_d]
    low_d = [s for s in sample if s['density'] <= median_d]

    # Baseline: lambda, rho, kappa, delta
    x_vars_base = ['lambda', 'rho', 'kappa', 'delta']
    r2_high, _ = ols_r2(high_d, 'delivery', x_vars_base)
    r2_low, _ = ols_r2(low_d, 'delivery', x_vars_base)
    p_val = bootstrap_p(high_d, low_d, x_vars_base, n_boot=N_BOOT)

    # Controls: size + region dummies (drop all-zero)
    active_controls = ['log_at']
    for reg in ['reg_europe', 'reg_north_am', 'reg_latin_am', 'reg_east_asia', 'reg_south_se_asia', 'reg_oceania']:
        has_h = any(s[reg] > 0 for s in high_d)
        has_l = any(s[reg] > 0 for s in low_d)
        if has_h and has_l:
            active_controls.append(reg)
    x_vars_ctrl = x_vars_base + active_controls
    r2_high_c, _ = ols_r2(high_d, 'delivery', x_vars_ctrl)
    r2_low_c, _ = ols_r2(low_d, 'delivery', x_vars_ctrl)
    p_val_c = bootstrap_p(high_d, low_d, x_vars_ctrl, n_boot=N_BOOT)

    results.append({
        'half_life': half_life,
        'median_d': median_d,
        'n_low': len(low_d),
        'n_high': len(high_d),
        'r2_low': r2_low,
        'r2_high': r2_high,
        'p': p_val,
        'r2_low_c': r2_low_c,
        'r2_high_c': r2_high_c,
        'p_c': p_val_c,
    })

print('\n=== KERNEL SENSITIVITY RESULTS ===')
for r in results:
    print(f'Half-life {r["half_life"]} km (median d={r["median_d"]:.3f}): '
          f'R2 low/high = {r["r2_low"]:.4f}/{r["r2_high"]:.4f} (p={r["p"]:.3f}); '
          f'with controls = {r["r2_low_c"]:.4f}/{r["r2_high_c"]:.4f} (p={r["p_c"]:.3f})')

# Write summary table
out_path = derived_path('results', 'summaries', 'strategy1_kernel_sensitivity_summary.md')
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('# Strategy 1: Kernel Half-Life Sensitivity\n\n')
    f.write('| Half-life (km) | Median density | R2 low | R2 high | p (baseline) | R2 low (ctrl) | R2 high (ctrl) | p (ctrl) | N low | N high |\n')
    f.write('|---|---|---|---|---|---|---|---|---|---|\n')
    for r in results:
        f.write(f'| {r["half_life"]} | {r["median_d"]:.3f} | {r["r2_low"]:.4f} | {r["r2_high"]:.4f} | {r["p"]:.3f} | '
                f'{r["r2_low_c"]:.4f} | {r["r2_high_c"]:.4f} | {r["p_c"]:.3f} | {r["n_low"]} | {r["n_high"]} |\n')

print(f'\nWrote summary: {out_path}')
