"""Compute regression results for reviewer-requested robustness tables.

This is the heavy-computation step. It loads raw data, runs all OLS
regressions, the placebo permutation test, VIF diagnostics, and
bandwidth sensitivity. Results are exported to results/json/referee_tables.json
for fast table generation by referee_tables.py.

Tasks:
1) Correlation matrix among w_geo, w_fuel, w_reg (3x3).
2) Channel decomposition with firm controls (Size, lambda, rho) for 3-month CARs.
3) Strong placebo: shuffle exposure networks (permute gvkeys) for 3-month CARs.
4) Specification progression (w_geo only, w_fuel only, both, full + controls).
5) Bandwidth sensitivity: w_geo at half-lives {250, 500, 750, 1000, 1500} km.
6) Controls sensitivity: size-only, leverage-only, firm FE.
7) VIF diagnostics.

Output: results/json/referee_tables.json
"""
import csv
import os
import json
import math
import random
from collections import defaultdict

from _paths import raw_path, derived_path, results_path

POST = 3
PRE_MONTHS = 24
RANDOM_SEED = 20260222
N_PLACEBO = 200


def load_monthly_returns(path):
    monthly = defaultdict(dict)
    with open(path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            d = row['datadate']
            if not d or len(d) < 7:
                continue
            ym = d[:7]
            try:
                r = float(row['ret_monthly'])
            except (ValueError, TypeError):
                continue
            monthly[gk][ym] = r
    return monthly


def load_ff_factors_monthly(path):
    mktrf = {}
    rf = {}
    vwretd = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('This file') or line.startswith('The '):
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
            date_fmt = f"{date[:4]}-{date[4:6]}"
            mktrf[date_fmt] = mktrf_dec
            rf[date_fmt] = rf_dec
            vwretd[date_fmt] = vwretd_dec
    return vwretd


def compute_monthly_car(gvkey, event_month, post, monthly_ret, market_ret):
    if gvkey not in monthly_ret:
        return None
    months = sorted(monthly_ret[gvkey].keys())
    if event_month not in months:
        return None
    event_idx = months.index(event_month)
    # pre-demean
    pre_mean_ar = 0.0
    ar_list = []
    for i in range(max(0, event_idx - PRE_MONTHS), event_idx):
        m = months[i]
        if m in monthly_ret[gvkey] and m in market_ret:
            ar_list.append(monthly_ret[gvkey][m] - market_ret[m])
    if ar_list:
        pre_mean_ar = sum(ar_list) / len(ar_list)
    car = 0.0
    for offset in range(-1, post + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months):
            m = months[idx]
            if m in monthly_ret[gvkey] and m in market_ret:
                ar = monthly_ret[gvkey][m] - market_ret[m]
                car += ar - pre_mean_ar
    return car


def load_weight_matrix(path, value_key='w_ij'):
    W = defaultdict(dict)
    if not os.path.exists(path):
        return W
    with open(path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gi = row['gvkey_i']
            gj = row['gvkey_j']
            try:
                w = float(row.get(value_key, row.get('w_ij', 0.0)))
            except (ValueError, TypeError):
                continue
            W[gi][gj] = w
    return W


def haversine(lat1, lon1, lat2, lon2):
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2) ** 2
    return R * 2 * math.asin(math.sqrt(a))


def build_geo_weights_from_centroids(centroid_path, half_life_km=500):
    centroids = {}
    with open(centroid_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                centroids[row['gvkey']] = (float(row['centroid_lat']), float(row['centroid_lon']))
            except (ValueError, TypeError):
                continue
    gvkeys = sorted(centroids.keys())
    decay_km = half_life_km / math.log(2)
    W_geo = defaultdict(dict)
    for i, gi in enumerate(gvkeys):
        lat_i, lon_i = centroids[gi]
        neighbors = {}
        row_sum = 0.0
        for j, gj in enumerate(gvkeys):
            if gi == gj:
                continue
            lat_j, lon_j = centroids[gj]
            d = haversine(lat_i, lon_i, lat_j, lon_j)
            if d <= 0:
                continue
            w = math.exp(-d / decay_km) / d
            neighbors[gj] = w
            row_sum += w
        if row_sum > 0:
            for gj, w in neighbors.items():
                W_geo[gi][gj] = w / row_sum
    return W_geo


def get_fundamentals_by_year(path):
    fundamentals_by_year = defaultdict(dict)
    fundamentals_latest = {}
    with open(path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            fy = row['fyear']
            fundamentals_by_year[gk][fy] = row
            if gk not in fundamentals_latest or fy > fundamentals_latest[gk]['fyear']:
                fundamentals_latest[gk] = row
    return fundamentals_by_year, fundamentals_latest


def get_fundamentals_for_year(gvkey, year, fundamentals_by_year, fundamentals_latest):
    rows = fundamentals_by_year.get(gvkey, {})
    if not rows:
        return fundamentals_latest.get(gvkey)
    years = [int(y) for y in rows.keys() if str(y).isdigit()]
    if not years:
        return fundamentals_latest.get(gvkey)
    years_le = [y for y in years if y <= year]
    chosen = max(years_le) if years_le else max(years)
    return rows.get(str(chosen)) or rows.get(chosen) or fundamentals_latest.get(gvkey)


def build_events(path):
    events = []
    with open(path, 'r', encoding='utf-8') as f:
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
            events.append({
                'event_date': effective_date,
                'year': event_year,
                'gvkeys': row['matched_gvkeys'].split(';')
            })
    return events


def ols(data, y_var, x_vars, cluster_var=None):
    if not data:
        return None
    y = [d[y_var] for d in data]
    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]
    n = len(y)
    k = len(x_vars) + 1
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]
    inv = invert_matrix(XtX)
    if inv is None:
        return None
    beta = [sum(inv[a][b] * Xty[b] for b in range(k)) for a in range(k)]
    yhat = [sum(beta[j] * X[i][j] for j in range(k)) for i in range(n)]
    resid = [y[i] - yhat[i] for i in range(n)]

    # cluster-robust variance
    if cluster_var is None:
        s2 = sum(r * r for r in resid) / (n - k)
        vcov = [[s2 * inv[i][j] for j in range(k)] for i in range(k)]
    else:
        if isinstance(cluster_var, list):
            clusters = list(zip(*[[d[c] for d in data] for c in cluster_var]))
        else:
            clusters = [d[cluster_var] for d in data]
        meat = [[0.0 for _ in range(k)] for _ in range(k)]
        scores = defaultdict(lambda: [0.0] * k)
        for i in range(n):
            key = clusters[i]
            xi = X[i]
            for a in range(k):
                scores[key][a] += xi[a] * resid[i]
        for s in scores.values():
            for a in range(k):
                for b in range(k):
                    meat[a][b] += s[a] * s[b]
        temp = mat_mul(inv, meat)
        vcov = mat_mul(temp, inv)

    se = [math.sqrt(max(vcov[i][i], 0.0)) for i in range(k)]
    out = {
        'beta': {'intercept': beta[0]},
        'se': {'intercept': se[0]},
        't': {'intercept': beta[0] / se[0] if se[0] > 0 else 0.0},
        'n': n
    }
    for i, xv in enumerate(x_vars, start=1):
        out['beta'][xv] = beta[i]
        out['se'][xv] = se[i]
        out['t'][xv] = beta[i] / se[i] if se[i] > 0 else 0.0
    return out


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


def corr(x, y):
    n = len(x)
    if n < 2:
        return None
    mean_x = sum(x) / n
    mean_y = sum(y) / n
    num = sum((x[i] - mean_x) * (y[i] - mean_y) for i in range(n))
    denx = sum((x[i] - mean_x) ** 2 for i in range(n))
    deny = sum((y[i] - mean_y) ** 2 for i in range(n))
    if denx <= 0 or deny <= 0:
        return None
    return num / math.sqrt(denx * deny)


def compute_vif(data, x_vars):
    """Return dict var -> VIF from the diagonal of (X'X)^{-1} * diag(X'X).

    VIF_j = [(X'X)^{-1}]_jj * (X'X)_jj, which avoids running k separate
    auxiliary regressions and is O(nk^2 + k^3) instead of O(nk^3).
    """
    vifs = {}
    if not data:
        return vifs
    n = len(data)
    k = len(x_vars)
    if k == 0:
        return vifs
    # Build X'X and invert once (k x k, no intercept needed for VIF)
    # Center variables for numerical stability
    means = {v: sum(d[v] for d in data) / n for v in x_vars}
    XtX = [[0.0] * k for _ in range(k)]
    for d in data:
        row = [d[v] - means[v] for v in x_vars]
        for a in range(k):
            for b in range(a, k):
                XtX[a][b] += row[a] * row[b]
    for a in range(k):
        for b in range(a + 1, k):
            XtX[b][a] = XtX[a][b]
    inv = invert_matrix(XtX)
    if inv is None:
        return {v: None for v in x_vars}
    for j, v in enumerate(x_vars):
        vifs[v] = inv[j][j] * XtX[j][j] if XtX[j][j] > 0 else 1.0
    return vifs


def serialize_ols_result(res, vars_list):
    """Convert an OLS result dict to a JSON-serializable dict."""
    if res is None:
        return None
    coefficients = {}
    for v in vars_list:
        coefficients[v] = {
            'beta': res['beta'][v],
            'se': res['se'][v],
            't': res['t'][v],
        }
    return {
        'n': res['n'],
        'vars': vars_list,
        'coefficients': coefficients,
    }


def main():
    random.seed(RANDOM_SEED)

    print('Loading data...', flush=True)
    monthly_ret = load_monthly_returns(derived_path('returns', 'monthly_returns.csv'))
    market_ret = load_ff_factors_monthly(raw_path('factors', 'F-F_Research_Data_Factors.csv'))
    if not market_ret:
        raise RuntimeError('Missing F-F factors for monthly vwretd')

    fundamentals_by_year, fundamentals_latest = get_fundamentals_by_year(
        derived_path('fundamentals', 'firm_fundamentals.csv')
    )

    W_geo = load_weight_matrix(derived_path('networks', 'weight_matrix_W_geo.csv'))
    W_fuel = load_weight_matrix(derived_path('networks', 'weight_matrix_W_fuel.csv'))
    W_reg = load_weight_matrix(derived_path('networks', 'weight_matrix_W_regulatory.csv'), value_key='w_reg')

    events = build_events(derived_path('events', 'coal_retirement_events.csv'))

    print('Building observation dataset...', flush=True)
    obs = []
    all_event_months = []
    car_cache = {}
    for event_id, event in enumerate(events):
        event_year = event['year']
        event_date = event['event_date']
        if event_date and len(event_date) >= 7:
            event_month = event_date[:7]
        else:
            event_month = f"{event_year}-07" if event_year else None
        if event_month:
            all_event_months.append(event_month)
        fm_sic4 = None
        for gk in event['gvkeys']:
            frow = get_fundamentals_for_year(gk, event_year, fundamentals_by_year, fundamentals_latest)
            if frow and frow.get('sic'):
                fm_sic4 = frow['sic'][:4]
                break
        for fm_gk in event['gvkeys']:
            neighbors = W_geo.get(fm_gk, {})
            neighbor_gks = set(neighbors.keys()) - set(event['gvkeys'])
            non_connected = [gk for gk in fundamentals_latest if gk not in event['gvkeys'] and gk not in neighbors]
            stable_seed = int.from_bytes(fm_gk.encode('utf-8'), 'little', signed=False) % (2**32)
            rng = random.Random(stable_seed)
            n_ctrl = min(len(non_connected), max(5 * len(neighbor_gks), 20))
            ctrl_sample = rng.sample(non_connected, n_ctrl) if len(non_connected) > n_ctrl else non_connected
            candidate_firms = list(neighbor_gks) + ctrl_sample
            for gk in candidate_firms:
                w_geo = neighbors.get(gk, 0.0)
                w_fuel = W_fuel.get(fm_gk, {}).get(gk, 0.0)
                w_reg = W_reg.get(fm_gk, {}).get(gk, 0.0)
                frow = get_fundamentals_for_year(gk, event_year, fundamentals_by_year, fundamentals_latest)
                if not frow:
                    continue
                sic = frow.get('sic')
                same_sector = 1.0 if (fm_sic4 and sic and sic[:4] == fm_sic4) else 0.0
                try:
                    at = float(frow['at']) if frow.get('at') else None
                except (ValueError, TypeError):
                    at = None
                log_assets = math.log(at) if at and at > 0 else None
                try:
                    lam = float(frow['lambda']) if frow.get('lambda') else None
                except (ValueError, TypeError):
                    lam = None
                try:
                    rho = float(frow['rho']) if frow.get('rho') else None
                except (ValueError, TypeError):
                    rho = None
                if event_month:
                    key = (gk, event_month)
                    if key in car_cache:
                        car = car_cache[key]
                    else:
                        car = compute_monthly_car(gk, event_month, POST, monthly_ret, market_ret)
                        car_cache[key] = car
                else:
                    car = None
                if car is None:
                    continue
                obs.append({
                    'car': car,
                    'w_geo': w_geo,
                    'w_fuel': w_fuel,
                    'w_reg': w_reg,
                    'same_sector': same_sector,
                    'log_assets': log_assets,
                    'lambda': lam,
                    'rho': rho,
                    'event_id': event_id,
                    'gvkey': gk,
                    'fm_gk': fm_gk,
                })

    print(f'  Total observations: {len(obs)}', flush=True)

    # ── Task 1: correlation matrix ──
    print('Computing correlation matrix...', flush=True)
    x_geo = [o['w_geo'] for o in obs]
    x_fuel = [o['w_fuel'] for o in obs]
    x_reg = [o['w_reg'] for o in obs]
    corrs = [
        [1.0, corr(x_geo, x_fuel), corr(x_geo, x_reg)],
        [corr(x_fuel, x_geo), 1.0, corr(x_fuel, x_reg)],
        [corr(x_reg, x_geo), corr(x_reg, x_fuel), 1.0],
    ]

    # ── Task 2: channel decomposition with controls ──
    print('Running channel decomposition with controls...', flush=True)
    obs_controls = [o for o in obs if o['log_assets'] is not None and o['lambda'] is not None and o['rho'] is not None]
    spec_controls = ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'log_assets', 'lambda', 'rho']
    res_controls = ols(obs_controls, 'car', spec_controls, cluster_var='event_id')

    # ── VIF diagnostics ──
    print('Computing VIF diagnostics...', flush=True)
    vif_vars = ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'log_assets', 'lambda', 'rho']
    vifs = compute_vif(obs_controls, vif_vars)

    # ── Task 3: strong placebo (shuffle exposure networks) ──
    print(f'Running placebo test ({N_PLACEBO} permutations)...', flush=True)
    placebo_betas = []
    gvkeys_all = sorted({o['gvkey'] for o in obs} | {o['fm_gk'] for o in obs})
    cars = [o['car'] for o in obs]
    w_reg_list = [o['w_reg'] for o in obs]
    same_list = [o['same_sector'] for o in obs]
    fm_list = [o['fm_gk'] for o in obs]
    gv_list = [o['gvkey'] for o in obs]

    for _ in range(N_PLACEBO):
        perm = gvkeys_all[:]
        random.shuffle(perm)
        perm_map = {g: p for g, p in zip(gvkeys_all, perm)}

        # accumulate XtX and Xty for regressors: [1, w_geo, w_fuel, w_reg, same]
        n = len(obs)
        s00 = n
        s01 = s02 = s03 = s04 = 0.0
        s11 = s12 = s13 = s14 = 0.0
        s22 = s23 = s24 = 0.0
        s33 = s34 = 0.0
        s44 = 0.0
        t0 = t1 = t2 = t3 = t4 = 0.0

        for i in range(n):
            fm = fm_list[i]
            gk = gv_list[i]
            fm_p = perm_map.get(fm, fm)
            gk_p = perm_map.get(gk, gk)
            w_geo_p = W_geo.get(fm_p, {}).get(gk_p, 0.0)
            w_fuel_p = W_fuel.get(fm_p, {}).get(gk_p, 0.0)
            w_reg = W_reg.get(fm_p, {}).get(gk_p, 0.0)
            same = same_list[i]
            y = cars[i]

            s01 += w_geo_p
            s02 += w_fuel_p
            s03 += w_reg
            s04 += same

            s11 += w_geo_p * w_geo_p
            s12 += w_geo_p * w_fuel_p
            s13 += w_geo_p * w_reg
            s14 += w_geo_p * same

            s22 += w_fuel_p * w_fuel_p
            s23 += w_fuel_p * w_reg
            s24 += w_fuel_p * same

            s33 += w_reg * w_reg
            s34 += w_reg * same

            s44 += same * same

            t0 += y
            t1 += y * w_geo_p
            t2 += y * w_fuel_p
            t3 += y * w_reg
            t4 += y * same

        XtX = [
            [s00, s01, s02, s03, s04],
            [s01, s11, s12, s13, s14],
            [s02, s12, s22, s23, s24],
            [s03, s13, s23, s33, s34],
            [s04, s14, s24, s34, s44],
        ]
        Xty = [t0, t1, t2, t3, t4]
        inv = invert_matrix(XtX)
        if inv is None:
            continue
        beta = [sum(inv[a][b] * Xty[b] for b in range(5)) for a in range(5)]
        placebo_betas.append((beta[1], beta[2], beta[3]))

    def mean_sd(vals):
        if not vals:
            return (0.0, 0.0)
        mean = sum(vals) / len(vals)
        var = sum((v - mean) ** 2 for v in vals) / (len(vals) - 1) if len(vals) > 1 else 0.0
        return mean, math.sqrt(var)

    geo_mean, geo_sd = mean_sd([b[0] for b in placebo_betas])
    fuel_mean, fuel_sd = mean_sd([b[1] for b in placebo_betas])
    reg_mean, reg_sd = mean_sd([b[2] for b in placebo_betas])

    # ── Task 4: spec progression ──
    print('Running specification progression...', flush=True)
    res_geo_only = ols(obs, 'car', ['w_geo'], cluster_var='event_id')
    res_fuel_only = ols(obs, 'car', ['w_fuel'], cluster_var='event_id')
    res_geo_fuel = ols(obs, 'car', ['w_geo', 'w_fuel'], cluster_var='event_id')
    res_full = ols(obs_controls, 'car', ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'log_assets', 'lambda', 'rho'], cluster_var='event_id')

    # Controls sensitivity
    print('Running controls sensitivity...', flush=True)
    res_size_only = ols([o for o in obs if o['log_assets'] is not None],
                        'car', ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'log_assets'], cluster_var='event_id')
    res_lev_only = ols([o for o in obs if o['lambda'] is not None],
                       'car', ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'lambda'], cluster_var='event_id')

    # Firm FE via within transformation by gvkey
    def demean_by_gvkey(data, y_var, x_vars):
        sums = defaultdict(lambda: defaultdict(float))
        counts = defaultdict(int)
        for d in data:
            g = d['gvkey']
            counts[g] += 1
            sums[g][y_var] += d[y_var]
            for v in x_vars:
                sums[g][v] += d[v]
        out = []
        for d in data:
            g = d['gvkey']
            if counts[g] <= 1:
                continue
            dd = {'gvkey': g, 'event_id': d['event_id']}
            dd[y_var] = d[y_var] - sums[g][y_var] / counts[g]
            for v in x_vars:
                dd[v] = d[v] - sums[g][v] / counts[g]
            out.append(dd)
        return out

    fe_vars = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']
    fe_data = demean_by_gvkey(obs, 'car', fe_vars)
    res_fe = ols(fe_data, 'car', fe_vars, cluster_var='event_id')

    # ── Task 5: bandwidth sensitivity (250, 500, 750, 1000, 1500 km) ──
    print('Building bandwidth sensitivity (250, 500, 750, 1000, 1500 km)...',
          flush=True)
    centroids_path = derived_path('networks', 'firm_centroids.csv')
    bandwidths = [250, 500, 750, 1000, 1500]
    W_geo_alt = {h: build_geo_weights_from_centroids(centroids_path, half_life_km=h)
                 for h in bandwidths}

    def rebuild_obs_with_geo(W_geo_alt):
        obs_alt = []
        for event_id, event in enumerate(events):
            event_year = event['year']
            event_date = event['event_date']
            if event_date and len(event_date) >= 7:
                event_month = event_date[:7]
            else:
                event_month = f"{event_year}-07" if event_year else None
            if not event_month:
                continue
            fm_sic4 = None
            for gk in event['gvkeys']:
                frow = get_fundamentals_for_year(gk, event_year, fundamentals_by_year, fundamentals_latest)
                if frow and frow.get('sic'):
                    fm_sic4 = frow['sic'][:4]
                    break
            for fm_gk in event['gvkeys']:
                neighbors = W_geo_alt.get(fm_gk, {})
                neighbor_gks = set(neighbors.keys()) - set(event['gvkeys'])
                non_connected = [gk for gk in fundamentals_latest if gk not in event['gvkeys'] and gk not in neighbors]
                stable_seed = int.from_bytes(fm_gk.encode('utf-8'), 'little', signed=False) % (2**32)
                rng = random.Random(stable_seed)
                n_ctrl = min(len(non_connected), max(5 * len(neighbor_gks), 20))
                ctrl_sample = rng.sample(non_connected, n_ctrl) if len(non_connected) > n_ctrl else non_connected
                candidate_firms = list(neighbor_gks) + ctrl_sample
                for gk in candidate_firms:
                    w_geo = neighbors.get(gk, 0.0)
                    w_fuel = W_fuel.get(fm_gk, {}).get(gk, 0.0)
                    w_reg = W_reg.get(fm_gk, {}).get(gk, 0.0)
                    frow = get_fundamentals_for_year(gk, event_year, fundamentals_by_year, fundamentals_latest)
                    if not frow:
                        continue
                    sic = frow.get('sic')
                    same_sector = 1.0 if (fm_sic4 and sic and sic[:4] == fm_sic4) else 0.0
                    key = (gk, event_month)
                    if key in car_cache:
                        car = car_cache[key]
                    else:
                        car = compute_monthly_car(gk, event_month, POST, monthly_ret, market_ret)
                        car_cache[key] = car
                    if car is None:
                        continue
                    obs_alt.append({
                        'car': car,
                        'w_geo': w_geo,
                        'w_fuel': w_fuel,
                        'w_reg': w_reg,
                        'same_sector': same_sector,
                        'event_id': event_id,
                        'gvkey': gk,
                    })
        return obs_alt

    bw_results = {}
    for h in bandwidths:
        obs_h = rebuild_obs_with_geo(W_geo_alt[h])
        bw_results[h] = ols(
            obs_h, 'car',
            ['w_geo', 'w_fuel', 'w_reg', 'same_sector'],
            cluster_var='event_id',
        )

    # ── Assemble JSON output ──
    print('Writing JSON...', flush=True)
    output = {
        'correlations': corrs,
        'channel_controls': serialize_ols_result(res_controls, spec_controls),
        'vifs': {v: vifs[v] for v in vif_vars},
        'placebo': {
            'n_permutations': N_PLACEBO,
            'w_geo': {'mean': geo_mean, 'sd': geo_sd},
            'w_fuel': {'mean': fuel_mean, 'sd': fuel_sd},
            'w_reg': {'mean': reg_mean, 'sd': reg_sd},
        },
        'spec_progression': {
            'geo_only': serialize_ols_result(res_geo_only, ['w_geo']),
            'fuel_only': serialize_ols_result(res_fuel_only, ['w_fuel']),
            'geo_fuel': serialize_ols_result(res_geo_fuel, ['w_geo', 'w_fuel']),
            'full_controls': serialize_ols_result(res_full, ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'log_assets', 'lambda', 'rho']),
        },
        'controls_sensitivity': {
            'size_only': serialize_ols_result(res_size_only, ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'log_assets']),
            'leverage_only': serialize_ols_result(res_lev_only, ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'lambda']),
            'firm_fe': serialize_ols_result(res_fe, ['w_geo', 'w_fuel', 'w_reg', 'same_sector']),
        },
        'bandwidth_sensitivity': {
            f'half_life_{h}km': serialize_ols_result(
                bw_results[h], ['w_geo', 'w_fuel', 'w_reg', 'same_sector'])
            for h in bandwidths
        },
    }

    json_dir = results_path('json')
    os.makedirs(json_dir, exist_ok=True)
    json_path = results_path('json', 'referee_tables.json')
    with open(json_path, 'w', encoding='utf-8') as f:
        json.dump(output, f, indent=2)

    print(f'Wrote: {json_path}', flush=True)
    print('Done (compute step).', flush=True)


if __name__ == '__main__':
    main()
