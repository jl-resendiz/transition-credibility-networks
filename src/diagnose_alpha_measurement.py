"""Diagnose measurement differences between GEM physical delivery and Refinitiv CO2/Revenue.

Outputs:
  - Summary correlations between within-firm changes in GEM clean share and Refinitiv CO2/Revenue.
  - OLS regression of Refinitiv change on GEM change with density interaction (firm-clustered SEs).

Interpretation:
  If correlation is low, within-firm Refinitiv variation is capturing reporting/operational noise
  rather than physical transformation. Density interaction checks whether measurement quality varies
  with network density.
"""
import csv
import math
from collections import defaultdict

from _paths import raw_path, derived_path


def safe_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def ols_clustered(data, y_var, x_vars, cluster_var):
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
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

    # Solve XtX * beta = Xty via Gaussian elimination
    aug = [row[:] + [Xty[a]] for a, row in enumerate(XtX)]
    for col in range(k):
        max_row = max(range(col, k), key=lambda r: abs(aug[r][col]))
        aug[col], aug[max_row] = aug[max_row], aug[col]
        pivot = aug[col][col]
        if abs(pivot) < 1e-20:
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

    # Cluster-robust SEs
    # Compute (X'X)^-1
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

    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    clusters = defaultdict(list)
    for i, d in enumerate(data):
        clusters[d[cluster_var]].append(i)

    S = [[0.0 for _ in range(k)] for _ in range(k)]
    for idxs in clusters.values():
        xu = [0.0 for _ in range(k)]
        for i in idxs:
            for a in range(k):
                xu[a] += X[i][a] * resid[i]
        for a in range(k):
            for b in range(k):
                S[a][b] += xu[a] * xu[b]

    # V = (X'X)^-1 S (X'X)^-1
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

    V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
    G = len(clusters)
    if G > 1:
        scale = (G / (G - 1)) * ((n - 1) / (n - k))
        for a in range(k):
            for b in range(k):
                V[a][b] *= scale
    se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]

    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0.0 for a in range(k)]
    names = ['intercept'] + x_vars
    return {
        'beta': dict(zip(names, beta)),
        'se': dict(zip(names, se)),
        't': dict(zip(names, t_stats)),
        'r2': r2,
        'n': n,
        'clusters': len(clusters),
    }


def corr(x, y):
    n = len(x)
    if n < 2:
        return None
    xbar = sum(x) / n
    ybar = sum(y) / n
    num = sum((xi - xbar) * (yi - ybar) for xi, yi in zip(x, y))
    denx = math.sqrt(sum((xi - xbar) ** 2 for xi in x))
    deny = math.sqrt(sum((yi - ybar) ** 2 for yi in y))
    if denx == 0 or deny == 0:
        return None
    return num / (denx * deny)


# Load GEM alpha panel (physical clean share)
gem_path = derived_path('fundamentals', 'firm_alpha_panel.csv')
gem = defaultdict(dict)  # gvkey -> year -> clean_share
with open(gem_path, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        try:
            year = int(row['year'])
        except (ValueError, TypeError):
            continue
        alpha = safe_float(row.get('alpha'))
        if alpha is None:
            continue
        gem[gk][year] = 1.0 - alpha

print(f'GEM firms: {len(gem)}')

# Load Refinitiv panel (CO2/Revenue)
ref_path = raw_path('refinitiv', 'refinitiv_panel.csv')
ref = defaultdict(dict)  # gvkey -> year -> co2_to_revenue
with open(ref_path, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row.get('gvkey')
        if not gk:
            continue
        try:
            year = int(row.get('year'))
        except (ValueError, TypeError):
            continue
        val = safe_float(row.get('co2_to_revenue'))
        if val is None or val <= 0:
            continue
        ref[gk][year] = val

print(f'Refinitiv firms: {len(ref)}')

# Load density (raw weight sum preferred)
density = {}
with open(derived_path('networks', 'firm_centroids.csv'), 'r', encoding='utf-8') as f:
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

dens_vals = sorted(density.values())
median_density = dens_vals[len(dens_vals) // 2] if dens_vals else 0
print(f'Median density: {median_density}')

# Build year-over-year deltas for matched firm-years
rows = []
for gk in set(gem.keys()) & set(ref.keys()):
    years = sorted(set(gem[gk].keys()) & set(ref[gk].keys()))
    for y in years:
        if (y - 1) not in gem[gk] or (y - 1) not in ref[gk]:
            continue
        clean_t = gem[gk][y]
        clean_lag = gem[gk][y - 1]
        d_clean = clean_t - clean_lag

        ref_t = ref[gk][y]
        ref_lag = ref[gk][y - 1]
        if ref_t <= 0 or ref_lag <= 0:
            continue
        # Use negative log-change so positive means cleaner (lower CO2/Revenue)
        d_ref = -(math.log(ref_t) - math.log(ref_lag))

        rows.append({
            'gvkey': gk,
            'year': y,
            'd_clean': d_clean,
            'd_ref': d_ref,
            'density': density.get(gk, 0),
        })

print(f'Matched firm-year deltas: {len(rows)}')

if rows:
    xs = [r['d_clean'] for r in rows]
    ys = [r['d_ref'] for r in rows]
    rho = corr(xs, ys)
    print(f'Overall corr(d_clean, -Delta log CO2/Rev): {rho:.3f}' if rho is not None else 'Overall corr: n/a')

    low = [r for r in rows if r['density'] <= median_density]
    high = [r for r in rows if r['density'] > median_density]
    if low:
        rho_low = corr([r['d_clean'] for r in low], [r['d_ref'] for r in low])
        print(f'Low-density corr: {rho_low:.3f}' if rho_low is not None else 'Low-density corr: n/a')
    if high:
        rho_high = corr([r['d_clean'] for r in high], [r['d_ref'] for r in high])
        print(f'High-density corr: {rho_high:.3f}' if rho_high is not None else 'High-density corr: n/a')

    # Regression: d_ref on d_clean and d_clean*density (standardized density)
    dens = [r['density'] for r in rows]
    dens_mean = sum(dens) / len(dens)
    dens_sd = math.sqrt(sum((d - dens_mean) ** 2 for d in dens) / (len(dens) - 1)) if len(dens) > 1 else 1.0
    for r in rows:
        z_d = (r['density'] - dens_mean) / dens_sd if dens_sd > 0 else 0.0
        r['z_density'] = z_d
        r['d_clean_x_zd'] = r['d_clean'] * z_d

    reg = ols_clustered(rows, 'd_ref', ['d_clean', 'd_clean_x_zd'], 'gvkey')
    if reg:
        print('\nRegression: d_ref = a + b1*d_clean + b2*(d_clean x z_density)')
        for name in ['intercept', 'd_clean', 'd_clean_x_zd']:
            print(f'  {name:15s} beta={reg["beta"][name]: .4f}  se={reg["se"][name]: .4f}  t={reg["t"][name]: .2f}')
        print(f'  R^2={reg["r2"]:.4f}, N={reg["n"]}, clusters={reg["clusters"]}')
    else:
        print('Regression could not be estimated (insufficient variation).')
