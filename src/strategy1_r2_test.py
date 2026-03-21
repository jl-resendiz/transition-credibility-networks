"""Strategy 1: Regime-Dependent R^2 Test.

Test whether fundamentals predict delivery less in high-density networks.
  - Delivery proxy = clean share = 1 - alpha (cross-sectional)
  - theta = (alpha, lambda, rho, kappa, delta)
  - Split by network density (above/below median neighbors)
  - H_coord: R^2_high < R^2_low
  - Refinement: interact theta with density, expect delta < 0
"""
import csv, os, math
from collections import defaultdict

from _paths import derived_path

# 1. Load firm fundamentals (latest available year per gvkey)
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

print(f'Firms with fundamentals: {len(fundamentals)}')

# Winsorize kappa (interest coverage) at 1st/99th percentiles
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
    print(f'Winsorized kappa at 1/99 pct: {kappa_p1:.3f}, {kappa_p99:.3f}')

# 2. Load network density from firm centroids (raw weight sum preferred)
density = {}
centroids_path = derived_path('networks', 'firm_centroids.csv')
if os.path.exists(centroids_path):
    with open(centroids_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            wsum = row.get('w_sum', '')
            if wsum not in ('', None):
                try:
                    density[gk] = float(wsum)
                except (ValueError, TypeError):
                    density[gk] = float(row.get('n_neighbors', 0))
            else:
                density[gk] = float(row.get('n_neighbors', 0))
    print(f'Firms with density (w_sum): {len(density)}')
else:
    raise SystemExit('Missing firm_centroids.csv for density.')

# 3. Load GEM data to compute delivery outcome (clean MW share change)
# We need MW by fuel at two points in time. For now, use alpha as a cross-sectional proxy:
# delivery = 1 - alpha (higher clean share = more transformation)
# In the full implementation, we'd compute alpha at t and t+k for time-series variation.

# 3b. Region mapping for country controls
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
# Everything else -> 'other' (Middle East, Africa, Central Asia)

# 4. Build analysis sample: firms with all five theta components + density + alpha
sample = []
for gk, f in fundamentals.items():
    if gk not in density:
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

    # Controls
    fic = f.get('fic', '')
    region = REGION_MAP.get(fic, 'other')
    try:
        log_at = math.log(max(float(f['at']), 1.0)) if f.get('at') else None
    except (ValueError, TypeError):
        log_at = None
    if log_at is None:
        continue  # need size control

    sample.append({
        'gvkey': gk,
        'conm': f['conm'],
        'fic': fic,
        'density': density[gk],
        'delivery': 1.0 - theta['alpha'],
        'log_at': log_at,
        # Region dummies (omit 'other' as reference)
        'reg_europe': 1.0 if region == 'europe' else 0.0,
        'reg_north_am': 1.0 if region == 'north_america' else 0.0,
        'reg_latin_am': 1.0 if region == 'latin_america' else 0.0,
        'reg_east_asia': 1.0 if region == 'east_asia' else 0.0,
        'reg_south_se_asia': 1.0 if region == 'south_se_asia' else 0.0,
        'reg_oceania': 1.0 if region == 'oceania' else 0.0,
        **theta,
    })

print(f'Analysis sample: {len(sample)} firms')

# 5. Split by network density median
densities = sorted(s['density'] for s in sample)
median_d = densities[len(densities) // 2]
print(f'Median network density: {median_d}')

high_d = [s for s in sample if s['density'] > median_d]
low_d = [s for s in sample if s['density'] <= median_d]
print(f'High-density firms: {len(high_d)}, Low-density firms: {len(low_d)}')

# 6. OLS regression: Delivery = beta' * theta + epsilon
def ols_r2(data, y_var, x_vars):
    """Simple OLS R^2 using normal equations."""
    n = len(data)
    k = len(x_vars) + 1  # +1 for intercept
    if n <= k:
        return None, None, None

    # Build matrices manually (no numpy)
    y = [d[y_var] for d in data]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)

    # X matrix: intercept + x_vars
    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]

    # X'X
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    # X'y
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

    # Solve via Gaussian elimination
    aug = [XtX[a] + [Xty[a]] for a in range(k)]
    for col in range(k):
        # Pivot
        max_row = max(range(col, k), key=lambda r: abs(aug[r][col]))
        aug[col], aug[max_row] = aug[max_row], aug[col]
        if abs(aug[col][col]) < 1e-12:
            return None, None, None
        for row in range(k):
            if row != col:
                factor = aug[row][col] / aug[col][col]
                for j in range(k + 1):
                    aug[row][j] -= factor * aug[col][j]

    beta = [aug[a][k] / aug[a][a] for a in range(k)]

    # Residuals and R^2
    y_hat = [sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
    ss_res = sum((y[i] - y_hat[i]) ** 2 for i in range(n))

    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0
    adj_r2 = 1 - (ss_res / (n - k)) / (ss_tot / (n - 1)) if n > k and ss_tot > 0 else 0

    return r2, adj_r2, dict(zip(['intercept'] + x_vars, beta))

x_vars_base = ['lambda', 'rho', 'kappa', 'delta']
controls = ['log_at', 'reg_europe', 'reg_north_am', 'reg_latin_am', 'reg_east_asia', 'reg_south_se_asia', 'reg_oceania']
x_vars = x_vars_base  # no-controls spec first

print(f'\n=== STRATEGY 1: REGIME-DEPENDENT R^2 (no controls) ===')
print(f'Y = Delivery (clean MW share)')
print(f'X = lambda, rho, kappa, delta')

r2_high, adj_r2_high, betas_high = ols_r2(high_d, 'delivery', x_vars)
r2_low, adj_r2_low, betas_low = ols_r2(low_d, 'delivery', x_vars)

print(f'\nHigh-density subsample (n={len(high_d)}):')
print(f'  R^2 = {r2_high:.4f}, Adj R^2 = {adj_r2_high:.4f}')
if betas_high:
    for var, b in betas_high.items():
        print(f'  beta_{var} = {b:.4f}')

print(f'\nLow-density subsample (n={len(low_d)}):')
print(f'  R^2 = {r2_low:.4f}, Adj R^2 = {adj_r2_low:.4f}')
if betas_low:
    for var, b in betas_low.items():
        print(f'  beta_{var} = {b:.4f}')

if r2_high is not None and r2_low is not None:
    print(f'\nR^2 difference: {r2_low - r2_high:.4f} (low - high)')
    if r2_low > r2_high:
        print('  --> CONSISTENT with coordination model (fundamentals predict less in dense networks)')
    else:
        print('  --> NOT consistent with coordination model')

# 7. Refinement: Pooled regression with density interaction
print(f'\n=== REFINEMENT: POOLED WITH DENSITY INTERACTION ===')

# Add interaction terms
for s in sample:
    for var in x_vars:
        s[f'{var}_x_d'] = s[var] * s['density']

x_vars_pooled = x_vars + ['density'] + [f'{v}_x_d' for v in x_vars]
r2_pooled, adj_r2_pooled, betas_pooled = ols_r2(sample, 'delivery', x_vars_pooled)

print(f'Pooled sample (n={len(sample)}):')
print(f'  R^2 = {r2_pooled:.4f}, Adj R^2 = {adj_r2_pooled:.4f}')
if betas_pooled:
    print(f'\n  Main effects:')
    for var in x_vars:
        print(f'    beta_{var} = {betas_pooled.get(var, 0):.4f}')
    print(f'    beta_density = {betas_pooled.get("density", 0):.4f}')
    print(f'\n  Interaction effects (theta x density):')
    for var in x_vars:
        b = betas_pooled.get(f'{var}_x_d', 0)
        sign = '-' if b < 0 else '+'
        print(f'    beta_{var}_x_d = {b:.4f}  ({sign} {"CONSISTENT" if b < 0 else "not consistent"} with coordination)')

# 8. Bootstrap R^2 difference
print(f'\n=== BOOTSTRAP TEST (R^2 difference) ===')
import random
random.seed(42)
N_BOOT = 1000
r2_diffs = []

for b in range(N_BOOT):
    boot_high = random.choices(high_d, k=len(high_d))
    boot_low = random.choices(low_d, k=len(low_d))
    r2_h, _, _ = ols_r2(boot_high, 'delivery', x_vars)
    r2_l, _, _ = ols_r2(boot_low, 'delivery', x_vars)
    if r2_h is not None and r2_l is not None:
        r2_diffs.append(r2_l - r2_h)

r2_diffs.sort()
n_boot = len(r2_diffs)
if n_boot > 0:
    mean_diff = sum(r2_diffs) / n_boot
    p025 = r2_diffs[int(0.025 * n_boot)]
    p975 = r2_diffs[int(0.975 * n_boot)]
    p_value = sum(1 for d in r2_diffs if d <= 0) / n_boot
    print(f'Bootstrap iterations: {n_boot}')
    print(f'Mean R^2 difference (low - high): {mean_diff:.4f}')
    print(f'95% CI: [{p025:.4f}, {p975:.4f}]')
    print(f'P-value (one-sided, H0: diff <= 0): {p_value:.4f}')
    if p_value < 0.05:
        print(f'  --> Significant at 5% level (p={p_value:.3f})')
    elif p_value < 0.10:
        print(f'  --> Significant at 10% level (p={p_value:.3f})')
    else:
        print(f'  --> NOT significant (p={p_value:.3f})')

# 9. WITH CONTROLS: size + region dummies
# Drop region dummies that are all-zero in either subsample to avoid collinearity
print(f'\n=== STRATEGY 1 WITH CONTROLS (size, region) ===')
active_controls = ['log_at']
for reg in ['reg_europe', 'reg_north_am', 'reg_latin_am', 'reg_east_asia', 'reg_south_se_asia', 'reg_oceania']:
    has_h = any(s[reg] > 0 for s in high_d)
    has_l = any(s[reg] > 0 for s in low_d)
    if has_h and has_l:
        active_controls.append(reg)
    else:
        print(f'  Dropping {reg} (all-zero in {"high" if not has_h else "low"} subsample)')

x_vars_ctrl = x_vars_base + active_controls

r2_high_c, adj_r2_high_c, _ = ols_r2(high_d, 'delivery', x_vars_ctrl)
r2_low_c, adj_r2_low_c, _ = ols_r2(low_d, 'delivery', x_vars_ctrl)

if r2_high_c is not None and r2_low_c is not None:
    print(f'High-density (n={len(high_d)}): R^2={r2_high_c:.4f}, Adj R^2={adj_r2_high_c:.4f}')
    print(f'Low-density  (n={len(low_d)}):  R^2={r2_low_c:.4f}, Adj R^2={adj_r2_low_c:.4f}')
else:
    print(f'High-density R^2: {"N/A" if r2_high_c is None else f"{r2_high_c:.4f}"}')
    print(f'Low-density  R^2: {"N/A" if r2_low_c is None else f"{r2_low_c:.4f}"}')
if r2_high_c is not None and r2_low_c is not None:
    diff_c = r2_low_c - r2_high_c
    print(f'R^2 difference: {diff_c:.4f} (low - high)')
    print(f'  {"CONSISTENT" if diff_c > 0 else "NOT consistent"} with coordination model')

# Bootstrap with controls
r2_diffs_c = []
for b in range(N_BOOT):
    boot_high = random.choices(high_d, k=len(high_d))
    boot_low = random.choices(low_d, k=len(low_d))
    r2_h, _, _ = ols_r2(boot_high, 'delivery', x_vars_ctrl)
    r2_l, _, _ = ols_r2(boot_low, 'delivery', x_vars_ctrl)
    if r2_h is not None and r2_l is not None:
        r2_diffs_c.append(r2_l - r2_h)

r2_diffs_c.sort()
n_bc = len(r2_diffs_c)
if n_bc > 0:
    p_val_c = sum(1 for d in r2_diffs_c if d <= 0) / n_bc
    mean_dc = sum(r2_diffs_c) / n_bc
    print(f'Bootstrap: mean diff={mean_dc:.4f}, p={p_val_c:.3f}')
    if p_val_c < 0.05:
        print(f'  --> Significant at 5% level')
    elif p_val_c < 0.10:
        print(f'  --> Significant at 10% level')
    else:
        print(f'  --> NOT significant')

# Interaction with controls
for s in sample:
    for var in x_vars_base:
        s[f'{var}_x_d_c'] = s[var] * s['density']

x_vars_int_ctrl = x_vars_ctrl + ['density'] + [f'{v}_x_d_c' for v in x_vars_base]
r2_ic, adj_ic, betas_ic = ols_r2(sample, 'delivery', x_vars_int_ctrl)
if betas_ic:
    print(f'\nPooled with controls + interaction (n={len(sample)}, R^2={r2_ic:.4f}):')
    print(f'  Interaction effects (theta x density):')
    for var in x_vars_base:
        b = betas_ic.get(f'{var}_x_d_c', 0)
        sign = '-' if b < 0 else '+'
        print(f'    beta_{var}_x_d = {b:.4f}  ({sign} {"CONSISTENT" if b < 0 else "not consistent"})')

# 10. Summary stats
print(f'\n=== SAMPLE SUMMARY STATISTICS ===')
for var in ['delivery', 'alpha', 'lambda', 'rho', 'kappa', 'delta', 'density', 'log_at']:
    vals = sorted(s[var] for s in sample)
    n = len(vals)
    mean = sum(vals) / n
    std = (sum((v - mean)**2 for v in vals) / (n - 1)) ** 0.5
    print(f'  {var:>10}: mean={mean:>8.3f}  std={std:>8.3f}  min={vals[0]:>8.3f}  median={vals[n//2]:>8.3f}  max={vals[-1]:>8.3f}')
