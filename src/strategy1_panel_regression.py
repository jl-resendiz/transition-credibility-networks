"""Strategy 1b: Panel regression of delivery on theta with firm and year fixed effects.

Uses time-varying alpha from firm_alpha_panel.csv (built from GEM start/retirement years).
Delivery_it = 1 - alpha_it (clean MW share in year t).
Theta_it = (lambda_it, rho_it, kappa_it, delta_it) from Compustat fundamentals.
Two-way demeaning (firm + year FE) with firm-clustered standard errors.

Also runs density-split panel regression (high vs low geographic density).
"""
import csv, os, math
from collections import defaultdict

from _paths import derived_path

WRITE_METRICS = os.getenv('WRITE_METRICS', '').strip().lower() in ('1', 'true', 'yes', 'y')

# 1. Load time-varying alpha
alpha_panel = {}  # (gvkey, year) -> alpha
alpha_path = derived_path('fundamentals', 'firm_alpha_panel.csv')
if os.path.exists(alpha_path):
    with open(alpha_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            alpha_panel[(row['gvkey'], int(row['year']))] = float(row['alpha'])
    print(f'Time-varying alpha: {len(alpha_panel)} firm-year observations')
else:
    print('WARNING: firm_alpha_panel.csv not found. Run build_time_varying_alpha.py first.')

# 2. Load fundamentals (all years)
print('Loading fundamentals...')
fund_rows = []
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = int(row['fyear'])
        try:
            lam = float(row['lambda']) if row.get('lambda') else None
            rho = float(row['rho']) if row.get('rho') else None
            kappa = float(row['kappa']) if row.get('kappa') else None
            delta = float(row['delta']) if row.get('delta') else None
        except (ValueError, TypeError):
            continue
        if None in (lam, rho, kappa, delta):
            continue

        # Get time-varying alpha for this firm-year
        alpha = alpha_panel.get((gk, fy))
        if alpha is None:
            continue

        fund_rows.append({
            'gvkey': gk,
            'fyear': fy,
            'delivery': 1.0 - alpha,
            'alpha': alpha,
            'lambda': lam,
            'rho': rho,
            'kappa': kappa,
            'delta': delta,
        })

# Winsorize kappa (interest coverage) at 1st/99th percentiles
kappa_vals = [r['kappa'] for r in fund_rows if r.get('kappa') is not None]
if kappa_vals:
    kappa_vals.sort()
    n_k = len(kappa_vals)
    kappa_p1 = kappa_vals[int(0.01 * (n_k - 1))]
    kappa_p99 = kappa_vals[int(0.99 * (n_k - 1))]
    for r in fund_rows:
        kap = r.get('kappa')
        if kap is None:
            continue
        r['kappa'] = min(max(kap, kappa_p1), kappa_p99)
    print(f'Winsorized kappa at 1/99 pct: {kappa_p1:.3f}, {kappa_p99:.3f}')

print(f'Panel rows (with time-varying alpha + complete theta): {len(fund_rows)}')
print(f'Firms: {len(set(r["gvkey"] for r in fund_rows))}')
print(f'Years: {sorted(set(r["fyear"] for r in fund_rows))}')

# Check within-firm variation in delivery
firm_deliveries = defaultdict(set)
for r in fund_rows:
    firm_deliveries[r['gvkey']].add(round(r['delivery'], 4))
n_varying = sum(1 for gk, vals in firm_deliveries.items() if len(vals) > 1)
print(f'Firms with within-firm delivery variation: {n_varying} / {len(firm_deliveries)} ({100*n_varying/len(firm_deliveries):.1f}%)')

# 3. Load geographic density (raw weight sum preferred)
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
    print(f'Firms with geographic density (w_sum): {len(density)}')

# 4. OLS with cluster-robust SEs (firm-level)
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


def ols(data, y_var, x_vars, cluster_var=None):
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
    adj_r2 = 1 - (ss_res / (n - k)) / (ss_tot / (n - 1)) if n > k else 0

    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    if cluster_var:
        clusters = {}
        for i, d in enumerate(data):
            cid = d.get(cluster_var)
            if cid not in clusters:
                clusters[cid] = []
            clusters[cid].append(i)
        S = [[0.0 for _ in range(k)] for _ in range(k)]
        for _, idxs in clusters.items():
            xu = [0.0 for _ in range(k)]
            for i in idxs:
                for a in range(k):
                    xu[a] += X[i][a] * resid[i]
            for a in range(k):
                for b in range(k):
                    S[a][b] += xu[a] * xu[b]
        V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
        G = len(clusters)
        if G > 1:
            scale = (G / (G - 1)) * ((n - 1) / (n - k))
            for a in range(k):
                for b in range(k):
                    V[a][b] *= scale
        se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]
    else:
        s2 = ss_res / (n - k) if n > k else 0
        se = [math.sqrt(s2 * inv_XtX[a][a]) if inv_XtX[a][a] > 0 else 0.0 for a in range(k)]

    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0.0 for a in range(k)]

    names = ['intercept'] + x_vars
    return {
        'beta': dict(zip(names, beta)),
        'se': dict(zip(names, se)),
        't': dict(zip(names, t_stats)),
        'r2': r2,
        'adj_r2': adj_r2,
        'n': n,
        'clusters': len(set(d.get(cluster_var) for d in data)) if cluster_var else None,
    }


x_vars = ['lambda', 'rho', 'kappa', 'delta']

# 5. Two-way demeaning (firm + year FE)
def two_way_demean(data, y_var, x_vars, firm_key='gvkey', time_key='fyear'):
    vars_all = [y_var] + x_vars
    overall = {v: 0.0 for v in vars_all}
    firm_sums = defaultdict(lambda: defaultdict(float))
    firm_counts = defaultdict(int)
    time_sums = defaultdict(lambda: defaultdict(float))
    time_counts = defaultdict(int)

    for row in data:
        f = row[firm_key]
        t = row[time_key]
        firm_counts[f] += 1
        time_counts[t] += 1
        for v in vars_all:
            val = row[v]
            overall[v] += val
            firm_sums[f][v] += val
            time_sums[t][v] += val

    n = len(data)
    if n == 0:
        return []
    for v in vars_all:
        overall[v] /= n

    firm_means = {f: {v: firm_sums[f][v] / firm_counts[f] for v in vars_all}
                  for f in firm_counts}
    time_means = {t: {v: time_sums[t][v] / time_counts[t] for v in vars_all}
                  for t in time_counts}

    out = []
    for row in data:
        f = row[firm_key]
        t = row[time_key]
        d = {
            firm_key: f,
            time_key: t,
            'density': row.get('density', 0),
        }
        for v in vars_all:
            d[v] = row[v] - firm_means[f][v] - time_means[t][v] + overall[v]
        out.append(d)
    return out


sample = []
for r in fund_rows:
    sample.append({
        'gvkey': r['gvkey'],
        'fyear': r['fyear'],
        'delivery': r['delivery'],
        'lambda': r['lambda'],
        'rho': r['rho'],
        'kappa': r['kappa'],
        'delta': r['delta'],
        'density': density.get(r['gvkey'], 0),
    })

sample_tw = two_way_demean(sample, 'delivery', x_vars)

# 6. Pooled panel regression
print('\n=== PANEL REGRESSION (TWO-WAY FE VIA DEMEANING, TIME-VARYING ALPHA) ===')
res = ols(sample_tw, 'delivery', x_vars, cluster_var='gvkey')
res_main = res
if res is not None:
    print(f'Within (two-way FE) R^2: {res["r2"]:.4f}, Adj R^2: {res["adj_r2"]:.4f}, N={res["n"]}, clusters={res.get("clusters")}')
    for var in ['intercept'] + x_vars:
        b = res['beta'].get(var, 0)
        se = res['se'].get(var, 0)
        t = res['t'].get(var, 0)
        sig = '***' if abs(t) > 2.576 else '**' if abs(t) > 1.96 else '*' if abs(t) > 1.645 else ''
        print(f'  {var:<10} beta={b:.4f}  se={se:.4f}  t={t:.2f} {sig}')

# 7. Density-split panel regression (geographic density)
res_h = None
res_l = None
diff = None
if density:
    densities_in_sample = [s['density'] for s in sample if s['density'] > 0]
    if densities_in_sample:
        densities_in_sample.sort()
        median_d = densities_in_sample[len(densities_in_sample) // 2]
        print(f'\nMedian geographic density: {median_d}')

        sample_with_d = [s for s in sample_tw if s['gvkey'] in density]
        high_d = [s for s in sample_with_d if s['density'] > median_d]
        low_d = [s for s in sample_with_d if s['density'] <= median_d]

        print(f'\n--- High-density panel (n={len(high_d)}) ---')
        res_h = ols(high_d, 'delivery', x_vars, cluster_var='gvkey')
        if res_h is not None:
            print(f'  Within (two-way FE) R^2: {res_h["r2"]:.4f}, Adj R^2: {res_h["adj_r2"]:.4f}')
            for var in ['intercept'] + x_vars:
                b = res_h['beta'].get(var, 0)
                print(f'  beta_{var} = {b:.4f}')

        print(f'\n--- Low-density panel (n={len(low_d)}) ---')
        res_l = ols(low_d, 'delivery', x_vars, cluster_var='gvkey')
        if res_l is not None:
            print(f'  Within (two-way FE) R^2: {res_l["r2"]:.4f}, Adj R^2: {res_l["adj_r2"]:.4f}')
            for var in ['intercept'] + x_vars:
                b = res_l['beta'].get(var, 0)
                print(f'  beta_{var} = {b:.4f}')

        if res_h is not None and res_l is not None:
            diff = res_l['r2'] - res_h['r2']
            print(f'\n  R^2 difference (low - high): {diff:.4f}')
            if diff > 0:
                print('  --> CONSISTENT with coordination model')
            else:
                print('  --> NOT consistent with coordination model')

# 8. Pooled panel with density interaction
print('\n=== PANEL WITH DENSITY INTERACTION ===')
for s in sample_tw:
    for v in x_vars:
        s[f'{v}_x_d'] = s[v] * s['density']

# Density main effect is time-invariant and absorbed by firm FE
x_vars_int = x_vars + [f'{v}_x_d' for v in x_vars]
sample_with_d = [s for s in sample_tw if s['gvkey'] in density]
res_int = ols(sample_with_d, 'delivery', x_vars_int, cluster_var='gvkey')
if res_int is not None:
    print(f'Within (two-way FE) R^2: {res_int["r2"]:.4f}, Adj R^2: {res_int["adj_r2"]:.4f}, N={res_int["n"]}, clusters={res_int.get("clusters")}')
    print(f'\n  Main effects:')
    for v in x_vars:
        print(f'    beta_{v} = {res_int["beta"].get(v, 0):.4f}')
    print(f'\n  Interaction effects (theta x density):')
    for v in x_vars:
        b = res_int['beta'].get(f'{v}_x_d', 0)
        sign = '-' if b < 0 else '+'
        print(f'    beta_{v}_x_d = {b:.6f}  ({sign} {"CONSISTENT" if b < 0 else "not consistent"})')

if WRITE_METRICS:
    metrics_path = os.path.join('JEEM_submission_package', 'JEEM_outputs', 'metrics',
                                'strategy1_panel_metrics.md')
    lines = [
        '# Strategy 1 Panel Metrics',
        '',
        f'- N: {res_main["n"] if res_main else "NA"}',
        '',
        '## Main (two-way FE)',
        f'- R2: {res_main["r2"]:.4f}  Adj R2: {res_main["adj_r2"]:.4f}' if res_main else '- R2: NA',
        '',
        '## Density split',
    ]
    if res_h and res_l:
        lines += [
            f'- High density R2: {res_h["r2"]:.4f}',
            f'- Low density R2: {res_l["r2"]:.4f}',
            f'- Diff (low - high): {diff:.4f}' if diff is not None else '- Diff (low - high): NA',
        ]
    else:
        lines.append('- Density split not available')
    if res_int:
        lines += [
            '',
            '## Density interactions',
            f'- R2: {res_int["r2"]:.4f}  Adj R2: {res_int["adj_r2"]:.4f}',
        ]
    os.makedirs(os.path.dirname(metrics_path), exist_ok=True)
    with open(metrics_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
