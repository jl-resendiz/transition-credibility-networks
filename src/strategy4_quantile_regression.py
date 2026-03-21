"""Strategy 4: Quantile regression of delivery on theta.

Estimate:
  Q_q(Delivery_i | theta_i) = beta_q' * theta_i
at q in {0.10, 0.25, 0.50, 0.75, 0.90}.

Under the coordination model, fundamentals set a ceiling (feasibility frontier)
but coordination determines whether the ceiling is reached. This predicts
|beta_q| increasing in q: fundamentals bind at upper quantiles (where the
feasibility constraint is active) but not at lower quantiles (where coordination
failure or low effort dominates).

Implementation: iteratively reweighted least squares (IRLS) for quantile regression.
"""
import csv, os, math
from collections import defaultdict

from _paths import derived_path


# ── Quantile regression via IRLS ─────────────────────────────────────

def quantile_regression(data, y_var, x_vars, tau, max_iter=100, tol=1e-6):
    """Quantile regression using iteratively reweighted least squares.

    Uses the Koenker-Bassett (1978) check function rho_tau(u) = u*(tau - I(u<0)).
    IRLS: weight w_i = tau / |u_i| if u_i > 0, (1-tau) / |u_i| if u_i < 0.
    """
    n = len(data)
    k = len(x_vars) + 1
    if n <= k:
        return None

    y = [d[y_var] for d in data]
    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]

    # Initialize with OLS
    # X'X
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]
    aug = [XtX[a] + [Xty[a]] for a in range(k)]
    for col in range(k):
        max_row = max(range(col, k), key=lambda r: abs(aug[r][col]))
        aug[col], aug[max_row] = aug[max_row], aug[col]
        if abs(aug[col][col]) < 1e-12:
            return None
        for row in range(k):
            if row != col:
                factor = aug[row][col] / aug[col][col]
                for j in range(k + 1):
                    aug[row][j] -= factor * aug[col][j]
    beta = [aug[a][k] / aug[a][a] for a in range(k)]

    # IRLS iterations
    for iteration in range(max_iter):
        # Compute residuals
        resid = [y[i] - sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]

        # Compute weights
        eps = 1e-6  # small constant to avoid division by zero
        weights = []
        for i in range(n):
            u = resid[i]
            abs_u = max(abs(u), eps)
            if u >= 0:
                weights.append(tau / abs_u)
            else:
                weights.append((1 - tau) / abs_u)

        # Weighted least squares: (X'WX)^{-1} X'Wy
        WX = [[weights[i] * X[i][a] for a in range(k)] for i in range(n)]
        XtWX = [[sum(X[i][a] * WX[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
        XtWy = [sum(X[i][a] * weights[i] * y[i] for i in range(n)) for a in range(k)]

        aug2 = [XtWX[a] + [XtWy[a]] for a in range(k)]
        for col in range(k):
            max_row = max(range(col, k), key=lambda r: abs(aug2[r][col]))
            aug2[col], aug2[max_row] = aug2[max_row], aug2[col]
            if abs(aug2[col][col]) < 1e-12:
                break
            for row in range(k):
                if row != col:
                    factor = aug2[row][col] / aug2[col][col]
                    for j in range(k + 1):
                        aug2[row][j] -= factor * aug2[col][j]
        else:
            new_beta = [aug2[a][k] / aug2[a][a] for a in range(k)]
            # Check convergence
            change = sum((new_beta[a] - beta[a]) ** 2 for a in range(k))
            beta = new_beta
            if change < tol:
                break
            continue
        break  # singular matrix

    # Final residuals and quantile loss
    resid = [y[i] - sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
    q_loss = sum(r * (tau - (1.0 if r < 0 else 0.0)) for r in resid)

    # Null model (intercept only): quantile of y
    y_sorted = sorted(y)
    q_idx = max(0, min(n - 1, int(tau * n)))
    y_q = y_sorted[q_idx]
    null_resid = [yi - y_q for yi in y]
    null_loss = sum(r * (tau - (1.0 if r < 0 else 0.0)) for r in null_resid)

    # Pseudo-R² (Koenker-Machado)
    pseudo_r2 = 1 - q_loss / null_loss if null_loss > 0 else 0

    # Bootstrap standard errors (simple residual bootstrap)
    import random
    random.seed(42 + int(tau * 1000))
    N_BOOT = 200
    boot_betas = [[] for _ in range(k)]

    for _ in range(N_BOOT):
        indices = random.choices(range(n), k=n)
        by = [y[i] for i in indices]
        bX = [X[i] for i in indices]

        # Quick WLS solve
        b_beta = beta[:]
        for _ in range(30):
            b_resid = [by[i] - sum(bX[i][a] * b_beta[a] for a in range(k)) for i in range(n)]
            b_weights = []
            for i in range(n):
                u = b_resid[i]
                abs_u = max(abs(u), eps)
                b_weights.append(tau / abs_u if u >= 0 else (1 - tau) / abs_u)

            bWX = [[b_weights[i] * bX[i][a] for a in range(k)] for i in range(n)]
            bXtWX = [[sum(bX[i][a] * bWX[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
            bXtWy = [sum(bX[i][a] * b_weights[i] * by[i] for i in range(n)) for a in range(k)]

            baug = [bXtWX[a] + [bXtWy[a]] for a in range(k)]
            ok = True
            for col in range(k):
                max_row = max(range(col, k), key=lambda r: abs(baug[r][col]))
                baug[col], baug[max_row] = baug[max_row], baug[col]
                if abs(baug[col][col]) < 1e-12:
                    ok = False
                    break
                for row in range(k):
                    if row != col:
                        factor = baug[row][col] / baug[col][col]
                        for j in range(k + 1):
                            baug[row][j] -= factor * baug[col][j]
            if not ok:
                break
            new_b = [baug[a][k] / baug[a][a] for a in range(k)]
            change = sum((new_b[a] - b_beta[a]) ** 2 for a in range(k))
            b_beta = new_b
            if change < tol:
                break

        for a in range(k):
            boot_betas[a].append(b_beta[a])

    # SE from bootstrap
    se = []
    for a in range(k):
        vals = boot_betas[a]
        mean_b = sum(vals) / len(vals)
        var_b = sum((v - mean_b) ** 2 for v in vals) / (len(vals) - 1) if len(vals) > 1 else 0
        se.append(math.sqrt(var_b))

    t_stats = [beta[a] / se[a] if se[a] > 1e-15 else 0 for a in range(k)]

    names = ['intercept'] + x_vars
    return {
        'beta': dict(zip(names, beta)),
        'se': dict(zip(names, se)),
        't': dict(zip(names, t_stats)),
        'pseudo_r2': pseudo_r2,
        'q_loss': q_loss,
        'n': n,
    }


# ── Load data ────────────────────────────────────────────────────────
print('Loading fundamentals...')
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

x_vars = ['lambda', 'rho', 'kappa', 'delta']

sample = []
for gk, f in fundamentals.items():
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
        'delivery': 1.0 - theta['alpha'],
        **{v: theta[v] for v in x_vars},
    })

print(f'Sample: {len(sample)} firms')

# ── Run quantile regressions ─────────────────────────────────────────
quantiles = [0.10, 0.25, 0.50, 0.75, 0.90]

print('\n' + '=' * 60)
print('STRATEGY 4: QUANTILE REGRESSION')
print('Q_q(Delivery | theta) = beta_q\' * theta')
print('Prediction: |beta_q| increasing in q')
print('=' * 60)

results = {}
for q in quantiles:
    print(f'\n--- Quantile tau = {q:.2f} ---')
    res = quantile_regression(sample, 'delivery', x_vars, q)
    if res is None:
        print('  Failed to converge')
        continue
    results[q] = res
    print(f'  N={res["n"]}, Pseudo-R²={res["pseudo_r2"]:.4f}')
    print(f'  {"Variable":<12} {"Beta":>10} {"SE":>10} {"t":>8} {"Sig":>5}')
    print(f'  {"-"*47}')
    for v in ['intercept'] + x_vars:
        b = res['beta'].get(v, 0)
        se = res['se'].get(v, 0)
        t = res['t'].get(v, 0)
        sig = '***' if abs(t) > 2.576 else '**' if abs(t) > 1.96 else '*' if abs(t) > 1.645 else ''
        print(f'  {v:<12} {b:>10.4f} {se:>10.4f} {t:>8.2f} {sig:>5}')

# ── Summary: coefficient magnitudes across quantiles ─────────────────
print('\n' + '=' * 60)
print('COEFFICIENT MAGNITUDES ACROSS QUANTILES')
print('=' * 60)

if len(results) >= 3:
    print(f'\n  {"Variable":<12}', end='')
    for q in quantiles:
        print(f'  q={q:.2f}', end='')
    print(f'  {"Trend":>10}')
    print(f'  {"-"*70}')

    for v in x_vars:
        print(f'  {v:<12}', end='')
        betas_across = []
        for q in quantiles:
            if q in results:
                b = results[q]['beta'].get(v, 0)
                betas_across.append(abs(b))
                print(f'  {b:>6.4f}', end='')
            else:
                betas_across.append(None)
                print(f'  {"N/A":>6}', end='')

        # Check if magnitudes increase across quantiles
        valid = [b for b in betas_across if b is not None]
        if len(valid) >= 3:
            increasing = all(valid[i] <= valid[i + 1] + 0.001 for i in range(len(valid) - 1))
            ratio = valid[-1] / valid[0] if valid[0] > 0.001 else float('inf')
            if increasing:
                print(f'  {"INCREASING":>10} (ratio={ratio:.1f}x)')
            elif valid[-1] > valid[0]:
                print(f'  {"MOSTLY UP":>10} (ratio={ratio:.1f}x)')
            else:
                print(f'  {"FLAT/DOWN":>10}')
        else:
            print(f'  {"N/A":>10}')

    # Pseudo-R² across quantiles
    print(f'\n  {"Pseudo-R²":<12}', end='')
    for q in quantiles:
        if q in results:
            print(f'  {results[q]["pseudo_r2"]:>6.4f}', end='')
        else:
            print(f'  {"N/A":>6}', end='')
    print()

    # Overall verdict
    print('\n=== COORDINATION MODEL PREDICTION ===')
    print('Prediction: |beta_q| increasing in q (fundamentals bind at upper quantiles)')
    n_increasing = 0
    for v in x_vars:
        valid = []
        for q in quantiles:
            if q in results:
                valid.append(abs(results[q]['beta'].get(v, 0)))
        if len(valid) >= 3 and valid[-1] > valid[0]:
            n_increasing += 1
    print(f'Variables with |beta| higher at q=0.9 than q=0.1: {n_increasing}/{len(x_vars)}')
    if n_increasing >= 3:
        print('  --> CONSISTENT with coordination model (fundamentals ceiling)')
    elif n_increasing >= 2:
        print('  --> PARTIALLY consistent')
    else:
        print('  --> NOT consistent')

    # Write metrics report
    metrics_path = os.path.join('JEEM_submission_package', 'JEEM_outputs', 'metrics',
                                'strategy4_quantile_metrics.md')
    lines = [
        '# Strategy 4 Quantile Regression Metrics',
        '',
        f'- N: {results[quantiles[0]]["n"] if quantiles and quantiles[0] in results else "NA"}',
        '',
        '| Quantile | Pseudo-R2 | ' + ' | '.join(x_vars) + ' |',
        '|---|---:|' + '|'.join(['---:'] * len(x_vars)) + '|',
    ]
    for q in quantiles:
        if q not in results:
            continue
        row = results[q]
        coeffs = [f'{row["beta"].get(v, 0.0):+.4f}' for v in x_vars]
        lines.append(f"| {q:.2f} | {row['pseudo_r2']:.4f} | " + ' | '.join(coeffs) + ' |')
    os.makedirs(os.path.dirname(metrics_path), exist_ok=True)
    with open(metrics_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
