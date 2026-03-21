"""Strategy 2: Conley (1999) spatial standard errors for channel decomposition.

Replicates the baseline channel decomposition regression (Table 2):
  CAR_j = alpha + beta_geo * w^geo_ij + beta_fuel * w^fuel_ij
        + beta_reg * w^reg_ij + beta_s * SameSector_j + eps_j

Then computes Conley (1999) standard errors using the Bartlett kernel
and great-circle (Haversine) distances between firm centroids. Tests
three spatial cutoffs (250 km, 500 km, 1000 km) alongside event-clustered
SEs for comparison.

Reference: Conley, T.G. (1999). 'GMM Estimation with Cross Sectional
Dependence', Journal of Econometrics, 92(1), 1-45.
"""
import csv
import os
import math
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

# ── Configuration ────────────────────────────────────────────────────

POST = 3            # months after event for CAR window
PRE_MONTHS = 24     # pre-event months for demeaning
CONLEY_CUTOFFS = [250, 500, 1000]  # km


# ── Haversine distance ──────────────────────────────────────────────

def haversine(lat1, lon1, lat2, lon2):
    """Great-circle distance in km between two points."""
    R = 6371  # km
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = (math.sin(dlat / 2) ** 2
         + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2))
         * math.sin(dlon / 2) ** 2)
    return 2 * R * math.asin(math.sqrt(a))


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


# ── Data loaders ────────────────────────────────────────────────────

def load_monthly_returns(path):
    """Load monthly returns into {gvkey: {YYYY-MM: ret}}."""
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
    """Load Fama-French monthly factors; return vwretd = Mkt-RF + RF."""
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
            vwretd_dec = (mktrf_val + rf_val) / 100.0
            date_fmt = f"{date[:4]}-{date[4:6]}"
            vwretd[date_fmt] = vwretd_dec
    return vwretd


def load_weight_matrix(path, value_key='w_ij'):
    """Load sparse weight matrix from CSV into {gvkey_i: {gvkey_j: w}}."""
    W = defaultdict(dict)
    if not os.path.exists(path):
        print(f'  WARNING: weight matrix not found: {path}')
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


def load_centroids(path):
    """Load firm centroids into {gvkey: (lat, lon)}."""
    centroids = {}
    with open(path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            try:
                centroids[row['gvkey']] = (
                    float(row['centroid_lat']),
                    float(row['centroid_lon'])
                )
            except (ValueError, TypeError, KeyError):
                continue
    return centroids


def load_fundamentals(path):
    """Load firm fundamentals; return by-year dict and latest-year dict."""
    by_year = defaultdict(dict)
    latest = {}
    with open(path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            fy = row['fyear']
            by_year[gk][fy] = row
            if gk not in latest or fy > latest[gk]['fyear']:
                latest[gk] = row
    return by_year, latest


def load_alpha_panel(path):
    """Load alpha panel into {gvkey: {year_str: alpha_float}}."""
    panel = defaultdict(dict)
    if not os.path.exists(path):
        return panel
    with open(path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            yr = row.get('year') or row.get('fyear')
            alpha = row.get('alpha', '')
            if gk and yr and alpha not in ('', None):
                try:
                    panel[gk][str(yr)] = float(alpha)
                except (ValueError, TypeError):
                    continue
    return panel


def load_events(path):
    """Load coal retirement events (first-movers only)."""
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
            if (effective_date and len(effective_date) >= 4
                    and effective_date[:4].isdigit()):
                event_year = int(effective_date[:4])
            else:
                event_year = (int(row['ret_year'])
                              if row.get('ret_year') else None)
            events.append({
                'event_date': effective_date,
                'year': event_year,
                'gvkeys': row['matched_gvkeys'].split(';'),
            })
    return events


def get_fundamentals_for_year(gvkey, year, by_year, latest):
    """Get fundamentals row for gvkey at or before year."""
    rows = by_year.get(gvkey, {})
    if not rows:
        return latest.get(gvkey)
    years = [int(y) for y in rows.keys() if str(y).isdigit()]
    if not years:
        return latest.get(gvkey)
    years_le = [y for y in years if y <= year]
    chosen = max(years_le) if years_le else max(years)
    return rows.get(str(chosen)) or latest.get(gvkey)


# ── CAR computation ─────────────────────────────────────────────────

def compute_monthly_car(gvkey, event_month, post, monthly_ret, market_ret):
    """Cumulative abnormal return over [-1, +post] months, pre-demeaned."""
    if gvkey not in monthly_ret:
        return None
    months = sorted(monthly_ret[gvkey].keys())
    if event_month not in months:
        return None
    event_idx = months.index(event_month)
    # pre-event mean AR for demeaning
    ar_list = []
    for i in range(max(0, event_idx - PRE_MONTHS), event_idx):
        m = months[i]
        if m in monthly_ret[gvkey] and m in market_ret:
            ar_list.append(monthly_ret[gvkey][m] - market_ret[m])
    pre_mean_ar = sum(ar_list) / len(ar_list) if ar_list else 0.0
    car = 0.0
    for offset in range(-1, post + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months):
            m = months[idx]
            if m in monthly_ret[gvkey] and m in market_ret:
                ar = monthly_ret[gvkey][m] - market_ret[m]
                car += ar - pre_mean_ar
    return car


# ── OLS regression ──────────────────────────────────────────────────

def ols_fit(data, y_var, x_vars):
    """OLS fit returning beta, residuals, X matrix, XtX_inv, and diagnostics."""
    n = len(data)
    k = len(x_vars) + 1  # +1 for intercept
    if n <= k + 1:
        return None

    y = [d[y_var] for d in data]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    if ss_tot < 1e-15:
        return None

    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]

    XtX = [[sum(X[i][a] * X[i][b] for i in range(n))
            for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]
    y_hat = [sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]
    ss_res = sum(r ** 2 for r in resid)
    r2 = 1 - ss_res / ss_tot

    return {
        'beta': beta,
        'resid': resid,
        'X': X,
        'inv_XtX': inv_XtX,
        'n': n,
        'k': k,
        'r2': r2,
        'ss_res': ss_res,
    }


def homoskedastic_se(fit):
    """Classical (non-robust) standard errors."""
    n, k = fit['n'], fit['k']
    s2 = fit['ss_res'] / (n - k)
    inv_XtX = fit['inv_XtX']
    return [math.sqrt(max(s2 * inv_XtX[a][a], 0.0)) for a in range(k)]


def event_clustered_se(fit, data, cluster_var='event_id'):
    """Cluster-robust SEs (clustered by event), with small-sample correction."""
    n, k = fit['n'], fit['k']
    X = fit['X']
    resid = fit['resid']
    inv_XtX = fit['inv_XtX']

    # Accumulate score by cluster
    clusters = defaultdict(list)
    for i, d in enumerate(data):
        clusters[d[cluster_var]].append(i)
    G = len(clusters)

    meat = [[0.0 for _ in range(k)] for _ in range(k)]
    for idxs in clusters.values():
        xu = [0.0] * k
        for i in idxs:
            for a in range(k):
                xu[a] += X[i][a] * resid[i]
        for a in range(k):
            for b in range(k):
                meat[a][b] += xu[a] * xu[b]

    V = mat_mul(mat_mul(inv_XtX, meat), inv_XtX)

    # Small-sample correction
    if G > 1:
        scale = (G / (G - 1)) * ((n - 1) / (n - k))
        for a in range(k):
            for b in range(k):
                V[a][b] *= scale

    return [math.sqrt(max(V[a][a], 0.0)) for a in range(k)]


def conley_se(fit, data, centroids, cutoff_km):
    """Conley (1999) spatial SEs with Bartlett kernel.

    Omega = sum_i sum_j K(d_ij/cutoff) * (X_i' e_i)(X_j' e_j)'
    V = (X'X)^{-1} Omega (X'X)^{-1}

    where K(z) = max(1 - z, 0) is the Bartlett kernel.
    """
    n, k = fit['n'], fit['k']
    X = fit['X']
    resid = fit['resid']
    inv_XtX = fit['inv_XtX']

    # Precompute score vectors: s_i = X_i * e_i  (k-vector for each obs)
    scores = []
    for i in range(n):
        s = [X[i][a] * resid[i] for a in range(k)]
        scores.append(s)

    # Look up centroid for each observation's gvkey
    obs_coords = []
    for d in data:
        gk = d['gvkey']
        if gk in centroids:
            obs_coords.append(centroids[gk])
        else:
            obs_coords.append(None)

    # Build Omega via pairwise kernel weighting
    omega = [[0.0 for _ in range(k)] for _ in range(k)]
    pairs_used = 0

    for i in range(n):
        if obs_coords[i] is None:
            continue
        lat_i, lon_i = obs_coords[i]
        # Diagonal term (distance = 0, kernel = 1)
        for a in range(k):
            for b in range(k):
                omega[a][b] += scores[i][a] * scores[i][b]
        for j in range(i + 1, n):
            if obs_coords[j] is None:
                continue
            lat_j, lon_j = obs_coords[j]
            d = haversine(lat_i, lon_i, lat_j, lon_j)
            kernel = max(1.0 - d / cutoff_km, 0.0)
            if kernel <= 0:
                continue
            pairs_used += 1
            for a in range(k):
                for b in range(k):
                    cross = kernel * (scores[i][a] * scores[j][b]
                                      + scores[j][a] * scores[i][b])
                    omega[a][b] += cross

    V = mat_mul(mat_mul(inv_XtX, omega), inv_XtX)
    se = [math.sqrt(max(V[a][a], 0.0)) for a in range(k)]
    return se, pairs_used


# ── Significance stars ──────────────────────────────────────────────

def stars(t):
    at = abs(t)
    if at >= 2.58:
        return '***'
    if at >= 1.96:
        return '**'
    if at >= 1.65:
        return '*'
    return ''


# ── Main ────────────────────────────────────────────────────────────

def main():
    print('=' * 70)
    print('Strategy 2: Conley (1999) Spatial Standard Errors')
    print('  Channel decomposition with spatial residual correlation')
    print('=' * 70)

    # Load data
    print('\nLoading monthly returns...')
    monthly_ret = load_monthly_returns(
        derived_path('returns', 'monthly_returns.csv'))
    print(f'  {len(monthly_ret)} firms')

    print('Loading Fama-French monthly factors...')
    market_ret = load_ff_factors_monthly(
        raw_path('factors', 'F-F_Research_Data_Factors.csv'))
    if not market_ret:
        raise RuntimeError('Missing Fama-French factors')
    print(f'  {len(market_ret)} months')

    print('Loading retirement events...')
    events = load_events(
        derived_path('events', 'coal_retirement_events.csv'))
    print(f'  {len(events)} first-mover events')

    print('Loading weight matrices...')
    W_geo = load_weight_matrix(
        derived_path('networks', 'weight_matrix_W_geo.csv'))
    W_fuel = load_weight_matrix(
        derived_path('networks', 'weight_matrix_W_fuel.csv'))
    W_reg = load_weight_matrix(
        derived_path('networks', 'weight_matrix_W_regulatory.csv'),
        value_key='w_reg')
    print(f'  W_geo:  {len(W_geo)} firms')
    print(f'  W_fuel: {len(W_fuel)} firms')
    print(f'  W_reg:  {len(W_reg)} firms')

    print('Loading firm centroids...')
    centroids = load_centroids(
        derived_path('networks', 'firm_centroids.csv'))
    print(f'  {len(centroids)} firms with centroids')

    print('Loading fundamentals...')
    fundamentals_by_year, fundamentals_latest = load_fundamentals(
        derived_path('fundamentals', 'firm_fundamentals.csv'))
    print(f'  {len(fundamentals_latest)} firms')

    print('Loading alpha panel...')
    alpha_panel = load_alpha_panel(
        derived_path('fundamentals', 'firm_alpha_panel.csv'))
    print(f'  {len(alpha_panel)} firms with alpha data')

    # ── Build event-firm observation pool ────────────────────────────

    print('\nBuilding event-firm observation pool...')
    obs = []
    car_cache = {}

    for event_id, event in enumerate(events):
        event_year = event['year']
        event_date = event['event_date']
        if event_date and len(event_date) >= 7:
            event_month = event_date[:7]
        else:
            event_month = f"{event_year}-07" if event_year else None
        if not event_month:
            continue

        # Identify first-mover SIC4 for same-sector indicator
        fm_sic4 = None
        for gk in event['gvkeys']:
            frow = get_fundamentals_for_year(
                gk, event_year, fundamentals_by_year, fundamentals_latest)
            if frow and frow.get('sic'):
                fm_sic4 = frow['sic'][:4]
                break

        # For each first-mover gvkey, find all neighbor/control firms
        for fm_gk in event['gvkeys']:
            neighbors_geo = W_geo.get(fm_gk, {})
            # Pool all firms that appear in any weight matrix
            candidate_gks = set()
            candidate_gks.update(neighbors_geo.keys())
            candidate_gks.update(W_fuel.get(fm_gk, {}).keys())
            candidate_gks.update(W_reg.get(fm_gk, {}).keys())
            # Also include firms with centroids but zero weight (controls)
            for gk in fundamentals_latest:
                if gk not in event['gvkeys']:
                    candidate_gks.add(gk)
            # Remove event firms themselves
            candidate_gks -= set(event['gvkeys'])

            for gk in candidate_gks:
                w_geo = W_geo.get(fm_gk, {}).get(gk, 0.0)
                w_fuel = W_fuel.get(fm_gk, {}).get(gk, 0.0)
                w_reg = W_reg.get(fm_gk, {}).get(gk, 0.0)

                frow = get_fundamentals_for_year(
                    gk, event_year, fundamentals_by_year,
                    fundamentals_latest)
                if not frow:
                    continue
                sic = frow.get('sic')
                same_sector = (1.0 if (fm_sic4 and sic
                                       and sic[:4] == fm_sic4)
                               else 0.0)

                # Compute CAR (with caching)
                key = (gk, event_month)
                if key in car_cache:
                    car = car_cache[key]
                else:
                    car = compute_monthly_car(
                        gk, event_month, POST, monthly_ret, market_ret)
                    car_cache[key] = car
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
                    'fm_gk': fm_gk,
                })

    print(f'  Total event-firm observations: {len(obs)}')
    n_events = len(set(o['event_id'] for o in obs))
    n_firms = len(set(o['gvkey'] for o in obs))
    print(f'  Unique events: {n_events}')
    print(f'  Unique firms:  {n_firms}')

    if not obs:
        print('ERROR: No observations built. Check data files.')
        return

    # ── Run OLS ─────────────────────────────────────────────────────

    x_vars = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']
    var_labels = ['intercept'] + x_vars

    print('\nFitting OLS...')
    fit = ols_fit(obs, 'car', x_vars)
    if fit is None:
        print('ERROR: OLS failed (singular matrix or insufficient data).')
        return

    beta = fit['beta']
    r2 = fit['r2']
    n = fit['n']
    k = fit['k']
    print(f'  N = {n}, k = {k}, R2 = {r2:.4f}')
    for i, name in enumerate(var_labels):
        print(f'  {name:15s}: beta = {beta[i]:+.6f}')

    # ── Event-clustered SEs ─────────────────────────────────────────

    print('\nComputing event-clustered standard errors...')
    se_event = event_clustered_se(fit, obs, cluster_var='event_id')
    print('  Done.')
    for i, name in enumerate(var_labels):
        t = beta[i] / se_event[i] if se_event[i] > 1e-15 else 0.0
        print(f'  {name:15s}: SE = {se_event[i]:.6f}, t = {t:+.3f}')

    # ── Conley SEs at each cutoff ───────────────────────────────────

    conley_results = {}
    for cutoff in CONLEY_CUTOFFS:
        print(f'\nComputing Conley SEs (cutoff = {cutoff} km)...')
        se_conley, pairs = conley_se(fit, obs, centroids, cutoff)
        conley_results[cutoff] = se_conley
        print(f'  Kernel-weighted pairs: {pairs}')
        for i, name in enumerate(var_labels):
            t = beta[i] / se_conley[i] if se_conley[i] > 1e-15 else 0.0
            print(f'  {name:15s}: SE = {se_conley[i]:.6f}, t = {t:+.3f}')

    # ── Write results markdown ──────────────────────────────────────

    out_path = results_path('metrics', 'strategy2_conley_se.md')
    os.makedirs(os.path.dirname(out_path), exist_ok=True)

    lines = []
    lines.append('# Conley (1999) Spatial Standard Errors')
    lines.append('')
    lines.append('Channel decomposition regression (monthly CAR [-1, +3]):')
    lines.append('```')
    lines.append('CAR_j = alpha + beta_geo * w^geo_ij + beta_fuel * w^fuel_ij'
                 ' + beta_reg * w^reg_ij + beta_s * SameSector_j + eps_j')
    lines.append('```')
    lines.append('')
    lines.append(f'N = {n}, R2 = {r2:.4f}, Events = {n_events}, '
                 f'Firms = {n_firms}')
    lines.append('')

    # Header row
    header = ('| Variable | Coef '
              '| Event-clust SE | t '
              '| Conley 250km SE | t '
              '| Conley 500km SE | t '
              '| Conley 1000km SE | t |')
    sep = ('|---|---:'
           '|---:|---:'
           '|---:|---:'
           '|---:|---:'
           '|---:|---:|')
    lines.append(header)
    lines.append(sep)

    display_names = {
        'intercept': 'Intercept',
        'w_geo': 'w_geo',
        'w_fuel': 'w_fuel',
        'w_reg': 'w_reg',
        'same_sector': 'SameSector',
    }

    for i, name in enumerate(var_labels):
        b = beta[i]
        # Event-clustered
        se_e = se_event[i]
        t_e = b / se_e if se_e > 1e-15 else 0.0
        # Conley at each cutoff
        se_250 = conley_results[250][i]
        t_250 = b / se_250 if se_250 > 1e-15 else 0.0
        se_500 = conley_results[500][i]
        t_500 = b / se_500 if se_500 > 1e-15 else 0.0
        se_1000 = conley_results[1000][i]
        t_1000 = b / se_1000 if se_1000 > 1e-15 else 0.0

        label = display_names.get(name, name)
        row = (f'| {label} | {b:+.4f}{stars(t_e)} '
               f'| {se_e:.4f} | {t_e:+.2f} '
               f'| {se_250:.4f} | {t_250:+.2f}{stars(t_250)} '
               f'| {se_500:.4f} | {t_500:+.2f}{stars(t_500)} '
               f'| {se_1000:.4f} | {t_1000:+.2f}{stars(t_1000)} |')
        lines.append(row)

    lines.append('')
    lines.append('Significance: \\*\\*\\* p<0.01, \\*\\* p<0.05, \\* p<0.10')
    lines.append('')
    lines.append('Notes:')
    lines.append('- Conley (1999) SEs use the Bartlett kernel: '
                 'K(d) = max(1 - d/cutoff, 0)')
    lines.append('- Distance is great-circle (Haversine) between firm '
                 'centroids')
    lines.append('- Event-clustered SEs include small-sample '
                 'correction: G/(G-1) * (N-1)/(N-k)')
    lines.append('- CARs are market-adjusted (ret - vwretd) and '
                 'pre-demeaned using 24-month pre-event window')

    # Ratio table: how much do Conley SEs differ from event-clustered?
    lines.append('')
    lines.append('## SE Ratios (Conley / Event-clustered)')
    lines.append('')
    ratio_header = '| Variable | 250km | 500km | 1000km |'
    ratio_sep = '|---|---:|---:|---:|'
    lines.append(ratio_header)
    lines.append(ratio_sep)
    for i, name in enumerate(var_labels):
        if name == 'intercept':
            continue
        label = display_names.get(name, name)
        r250 = (conley_results[250][i] / se_event[i]
                if se_event[i] > 1e-15 else 0.0)
        r500 = (conley_results[500][i] / se_event[i]
                if se_event[i] > 1e-15 else 0.0)
        r1000 = (conley_results[1000][i] / se_event[i]
                 if se_event[i] > 1e-15 else 0.0)
        lines.append(f'| {label} | {r250:.3f} | {r500:.3f} | {r1000:.3f} |')

    with open(out_path, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))

    print(f'\nResults written to: {out_path}')
    print('Done.')


if __name__ == '__main__':
    main()
