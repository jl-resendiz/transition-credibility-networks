"""Strategy 1 forward-delivery R^2 split (GEM physical delivery change).

Tests whether fundamentals at time t predict future delivery change (clean share)
over horizons h in {3,5} years, and whether predictive power differs by network density.

Also reports coverage by density and a common-shock control via country x year demeaning.
"""
import csv
import math
from collections import defaultdict

from _paths import derived_path

HORIZONS = [3, 5]


def safe_float(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return None


def ols_r2(data, y_var, x_vars):
    n = len(data)
    k = len(x_vars) + 1
    if n <= k:
        return None, None
    y = [d[y_var] for d in data]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
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
    r2 = 1 - ss_res / ss_tot if ss_tot > 0 else 0.0
    adj_r2 = 1 - (ss_res / (n - k)) / (ss_tot / (n - 1)) if n > k and ss_tot > 0 else 0.0
    return r2, adj_r2


def demean_by_group(data, group_key, vars_list):
    sums = defaultdict(lambda: defaultdict(float))
    counts = defaultdict(int)
    for row in data:
        g = row[group_key]
        counts[g] += 1
        for v in vars_list:
            sums[g][v] += row[v]
    out = []
    for row in data:
        g = row[group_key]
        if counts[g] < 2:
            # drop singletons to avoid zero-variance artifacts
            continue
        new_row = row.copy()
        for v in vars_list:
            new_row[v] = row[v] - sums[g][v] / counts[g]
        out.append(new_row)
    return out, len(out)


# Load GEM alpha panel (physical clean share)
alpha_panel = defaultdict(dict)
with open(derived_path('fundamentals', 'firm_alpha_panel.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        yr = row.get('year') or row.get('fyear')
        a = safe_float(row.get('alpha'))
        if gk and yr and a is not None:
            alpha_panel[gk][int(yr)] = a

# Load fundamentals by year
fund_by_year = defaultdict(dict)
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        fund_by_year[gk][int(fy)] = row

# Density from centroids (raw weight sum preferred)
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


def build_panel(h):
    panel = []
    for gk, yr_map in alpha_panel.items():
        for y, a0 in yr_map.items():
            y1 = y + h
            if y1 not in yr_map:
                continue
            f = fund_by_year.get(gk, {}).get(y)
            if not f:
                continue
            lamb = safe_float(f.get('lambda'))
            rho = safe_float(f.get('rho'))
            kappa = safe_float(f.get('kappa'))
            delta = safe_float(f.get('delta'))
            if None in (lamb, rho, kappa, delta):
                continue
            # Forward delivery change (clean share)
            delivery0 = 1.0 - a0
            delivery1 = 1.0 - yr_map[y1]
            d_deliv = delivery1 - delivery0

            fic = f.get('fic', '')
            try:
                log_at = math.log(max(float(f['at']), 1.0)) if f.get('at') else None
            except (ValueError, TypeError):
                log_at = None
            if log_at is None:
                continue

            panel.append({
                'gvkey': gk,
                'year': y,
                'fic': fic,
                'density': density.get(gk, 0),
                'd_delivery': d_deliv,
                'alpha': a0,
                'lambda': lamb,
                'rho': rho,
                'kappa': kappa,
                'delta': delta,
                'log_at': log_at,
            })
    return panel


def run_split(panel, label, x_vars):
    if not panel:
        print(f'{label}: no data')
        return
    dens_vals = sorted(p['density'] for p in panel)
    median_d = dens_vals[len(dens_vals) // 2]
    low = [p for p in panel if p['density'] <= median_d]
    high = [p for p in panel if p['density'] > median_d]
    r2_h, adj_h = ols_r2(high, 'd_delivery', x_vars)
    r2_l, adj_l = ols_r2(low, 'd_delivery', x_vars)
    print(f'\n{label}')
    print(f'  Median density: {median_d}')
    print(f'  N low={len(low)}, N high={len(high)}, firms={len(set(p["gvkey"] for p in panel))}')
    print(f'  R2 low={r2_l:.4f}, high={r2_h:.4f}, diff={r2_l - r2_h:+.4f}')
    return r2_l, r2_h, median_d, len(low), len(high)


print('=== Forward Delivery R2 Split (GEM clean-share change) ===')

for h in HORIZONS:
    panel = build_panel(h)
    print(f'\n--- Horizon: {h} years ---')
    print(f'Panel rows: {len(panel)}')
    # Coverage by density
    dens_vals = sorted(p['density'] for p in panel)
    if dens_vals:
        median_d = dens_vals[len(dens_vals) // 2]
        low = [p for p in panel if p['density'] <= median_d]
        high = [p for p in panel if p['density'] > median_d]
        print(f'Coverage: low={len(low)} rows, high={len(high)} rows; firms_low={len(set(p["gvkey"] for p in low))}, firms_high={len(set(p["gvkey"] for p in high))}')

    # Base spec
    x_base = ['alpha', 'lambda', 'rho', 'kappa', 'delta']
    run_split(panel, f'Base spec (theta), h={h}', x_base)

    # With size control
    x_ctrl = x_base + ['log_at']
    run_split(panel, f'With size control, h={h}', x_ctrl)

    # Country x year demeaning (common-shock control)
    vars_all = ['d_delivery'] + x_ctrl
    panel_cy, n_cy = demean_by_group(panel, 'fic', vars_all)
    # Within country groups per year
    panel_cy2 = []
    for row in panel_cy:
        row['cy'] = f'{row["fic"]}_{row["year"]}'
        panel_cy2.append(row)
    panel_cy2, n_cy2 = demean_by_group(panel_cy2, 'cy', vars_all)
    run_split(panel_cy2, f'Country x year FE (demeaned), h={h}', x_ctrl)
