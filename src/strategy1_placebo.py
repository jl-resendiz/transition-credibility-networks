"""Placebo weight matrix test for Strategy 1.

Randomly permute firm identities in the geographic weight matrix, preserving
the sparsity structure (same degree distribution, same weights, but random
assignment of which firm gets which position). Re-run the regime-dependent
R² test. If the observed R² difference is also present with random networks,
the spatial structure is not doing the work.

Also tests alternative distance thresholds (250km, 1000km).
"""
import csv, os, math, random
from collections import defaultdict, Counter

from _paths import derived_path
random.seed(42)


# ── OLS ──────────────────────────────────────────────────────────────
def ols_r2(data, y_var, x_vars):
    n = len(data)
    k = len(x_vars) + 1
    if n <= k:
        return None, None

    y = [d[y_var] for d in data]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    if ss_tot < 1e-15:
        return 0.0, 0.0

    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

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
    adj_r2 = 1 - (ss_res / (n - k)) / (ss_tot / (n - 1)) if n > k else 0
    return r2, adj_r2


# ── Load fundamentals (latest year per firm) ─────────────────────────
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

x_vars = ['lambda', 'rho', 'kappa', 'delta']


def build_sample(density_map):
    """Build analysis sample given a density mapping."""
    sample = []
    for gk, f in fundamentals.items():
        if gk not in density_map:
            continue
        skip = False
        theta = {}
        for var in ['alpha'] + x_vars:
            val = f.get(var, '')
            if val == '':
                skip = True
                break
            theta[var] = float(val)
        if skip:
            continue
        sample.append({
            'gvkey': gk,
            'density': density_map[gk],
            'delivery': 1.0 - theta['alpha'],
            **{v: theta[v] for v in x_vars},
        })
    return sample


def run_r2_test(sample):
    """Split by median density and return (r2_low - r2_high)."""
    if len(sample) < 20:
        return None
    densities = sorted(s['density'] for s in sample)
    median_d = densities[len(densities) // 2]
    high_d = [s for s in sample if s['density'] > median_d]
    low_d = [s for s in sample if s['density'] <= median_d]
    if len(high_d) < 10 or len(low_d) < 10:
        return None
    r2_h, _ = ols_r2(high_d, 'delivery', x_vars)
    r2_l, _ = ols_r2(low_d, 'delivery', x_vars)
    if r2_h is None or r2_l is None:
        return None
    return r2_l - r2_h


# ── Load actual geographic W and compute true density ────────────────
print('Loading firm centroids for density...')
geo_firms = set()
true_density = {}
with open(derived_path('networks', 'firm_centroids.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        geo_firms.add(gk)
        wsum = row.get('w_sum', '')
        if wsum not in ('', None):
            try:
                true_density[gk] = float(wsum)
            except (ValueError, TypeError):
                true_density[gk] = float(row.get('n_neighbors', 0))
        else:
            true_density[gk] = float(row.get('n_neighbors', 0))

# True R² difference
true_sample = build_sample(true_density)
true_diff = run_r2_test(true_sample)
print(f'True R² difference (low - high): {true_diff:.4f}')
print(f'Sample size: {len(true_sample)}')


# ══════════════════════════════════════════════════════════════════════
# PART 1: PLACEBO WEIGHT MATRICES
# ══════════════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('PLACEBO WEIGHT MATRIX TEST')
print('Permute firm identities, preserve edge structure')
print('=' * 60)

N_PLACEBO = 500
placebo_diffs = []
firm_list = sorted(geo_firms)
density_vals = [true_density[gk] for gk in firm_list]

for p in range(N_PLACEBO):
    # Random permutation of firm identities
    perm = firm_list[:]
    random.shuffle(perm)
    mapping = dict(zip(firm_list, perm))

    # Build permuted density (shuffle density values across firms)
    perm_density = {gk: density_vals[i] for i, gk in enumerate(perm)}
    perm_sample = build_sample(perm_density)
    diff = run_r2_test(perm_sample)
    if diff is not None:
        placebo_diffs.append(diff)

placebo_diffs.sort()
n_p = len(placebo_diffs)
if n_p > 0:
    mean_placebo = sum(placebo_diffs) / n_p
    p_value = sum(1 for d in placebo_diffs if d >= true_diff) / n_p
    p025 = placebo_diffs[int(0.025 * n_p)]
    p975 = placebo_diffs[int(0.975 * n_p)]

    print(f'\nPlacebo iterations: {n_p}')
    print(f'True R² diff:    {true_diff:.4f}')
    print(f'Mean placebo:    {mean_placebo:.4f}')
    print(f'Placebo 95% CI:  [{p025:.4f}, {p975:.4f}]')
    print(f'P-value (frac placebo >= true): {p_value:.4f}')

    if p_value < 0.05:
        print(f'  --> TRUE spatial structure significant at 5% (p={p_value:.3f})')
        print(f'  --> Placebo networks do NOT replicate the pattern')
    elif p_value < 0.10:
        print(f'  --> Marginally significant at 10% (p={p_value:.3f})')
    else:
        print(f'  --> NOT significant (p={p_value:.3f})')
        print(f'  --> Random networks produce similar R² differences')


# ══════════════════════════════════════════════════════════════════════
# PART 2: ALTERNATIVE DISTANCE THRESHOLDS
# ══════════════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('ALTERNATIVE DISTANCE THRESHOLDS')
print('=' * 60)

# Load centroids
centroids = {}
with open(derived_path('networks', 'firm_centroids.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            centroids[row['gvkey']] = (float(row['centroid_lat']), float(row['centroid_lon']))
        except (ValueError, TypeError):
            pass

print(f'Firms with centroids: {len(centroids)}')


def haversine(lat1, lon1, lat2, lon2):
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2) ** 2
    return R * 2 * math.asin(math.sqrt(a))


def build_density_at_threshold(threshold_km, decay_km=None):
    """Build density map at given distance threshold."""
    if decay_km is None:
        decay_km = threshold_km
    gvkeys = sorted(centroids.keys())
    density_map = {}
    for gi in gvkeys:
        lat_i, lon_i = centroids[gi]
        count = 0
        for gj in gvkeys:
            if gi == gj:
                continue
            lat_j, lon_j = centroids[gj]
            d = haversine(lat_i, lon_i, lat_j, lon_j)
            if 0 < d <= threshold_km:
                count += 1
        density_map[gi] = count
    return density_map


for threshold in [250, 500, 1000]:
    print(f'\n--- Threshold: {threshold} km ---')
    density_map = build_density_at_threshold(threshold)
    n_connected = sum(1 for d in density_map.values() if d > 0)
    avg_d = sum(density_map.values()) / len(density_map) if density_map else 0
    densities = sorted(density_map.values())
    median_d = densities[len(densities) // 2] if densities else 0
    print(f'  Connected firms: {n_connected}/{len(density_map)}')
    print(f'  Avg density: {avg_d:.1f}, Median: {median_d}')

    sample = build_sample(density_map)
    diff = run_r2_test(sample)
    if diff is not None:
        print(f'  R² diff (low - high): {diff:.4f}')
        print(f'  {"CONSISTENT" if diff > 0 else "NOT consistent"} with coordination model')

        # Quick bootstrap for this threshold
        if len(sample) >= 20:
            ds = sorted(s['density'] for s in sample)
            med = ds[len(ds) // 2]
            high = [s for s in sample if s['density'] > med]
            low = [s for s in sample if s['density'] <= med]
            boot_diffs = []
            for _ in range(500):
                bh = random.choices(high, k=len(high))
                bl = random.choices(low, k=len(low))
                r2h, _ = ols_r2(bh, 'delivery', x_vars)
                r2l, _ = ols_r2(bl, 'delivery', x_vars)
                if r2h is not None and r2l is not None:
                    boot_diffs.append(r2l - r2h)
            if boot_diffs:
                p_val = sum(1 for d in boot_diffs if d <= 0) / len(boot_diffs)
                print(f'  Bootstrap p-value: {p_val:.3f}')
    else:
        print(f'  Could not compute R² difference')


# ══════════════════════════════════════════════════════════════════════
# PART 3: RESTRICT TO CONDITIONAL TRANSFORMERS
# ══════════════════════════════════════════════════════════════════════
print('\n' + '=' * 60)
print('RESTRICT TO CONDITIONAL TRANSFORMERS')
print('(Effects should strengthen per paper)')
print('=' * 60)

conditional_sample = []
for s in true_sample:
    f = fundamentals.get(s['gvkey'])
    if not f or f.get('alpha', '') == '':
        continue
    alpha = float(f['alpha'])
    if 0.3 <= alpha <= 0.7:
        conditional_sample.append(s)

print(f'Conditional transformers: {len(conditional_sample)}')
if len(conditional_sample) >= 20:
    diff_cond = run_r2_test(conditional_sample)
    if diff_cond is not None:
        print(f'R² diff (conditional only): {diff_cond:.4f}')
        print(f'R² diff (full sample):      {true_diff:.4f}')
        if diff_cond > true_diff:
            print('  --> STRENGTHENED (consistent with coordination model)')
        else:
            print('  --> WEAKENED (not consistent)')
    else:
        print('  Could not compute (insufficient variation)')
else:
    print('  Too few conditional transformers for split regression')
