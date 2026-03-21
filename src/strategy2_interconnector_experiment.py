import csv
import math
import os
from collections import defaultdict
from datetime import datetime

from _paths import derived_path, results_path, raw_path

DECAY_KM = 500 / math.log(2)  # consistent with build_weight_matrix.py
NEAR_THRESHOLDS = [500, 1000]
PRE_MONTHS = (-6, -1)
POST_MONTHS = [3, 6, 12]
DAILY_WINDOW = (-1, 20)

COUNTRY_MAP = {
    'Australia (Tasmania)': 'AUS',
    'Australia (Victoria)': 'AUS',
    'Bangladesh': 'BGD',
    'Belgium': 'BEL',
    'Bhutan': 'BTN',
    'Denmark': 'DNK',
    'Estonia': 'EST',
    'Finland': 'FIN',
    'France': 'FRA',
    'Germany': 'DEU',
    'Greece (Crete)': 'GRC',
    'Greece (mainland)': 'GRC',
    'India': 'IND',
    'Indonesia (W.Kalimantan)': 'IDN',
    'Ireland': 'IRL',
    'Italy': 'ITA',
    'Italy (Sardinia)': 'ITA',
    'Italy (mainland)': 'ITA',
    'Japan (Hokkaido)': 'JPN',
    'Japan (Honshu)': 'JPN',
    'Laos': 'LAO',
    'Lithuania': 'LTU',
    'Malaysia (Sarawak)': 'MYS',
    'Montenegro': 'MNE',
    'Myanmar': 'MMR',
    'Nepal': 'NPL',
    'Netherlands': 'NLD',
    'Norway': 'NOR',
    'Poland': 'POL',
    'Spain': 'ESP',
    'Sweden': 'SWE',
    'UK': 'GBR',
    'UK (Wales)': 'GBR',
    # ambiguous multi-country endpoint; skip
    'Thailand-Malaysia-Singapore': None,
}


def haversine(lat1, lon1, lat2, lon2):
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(a))


def weight_from_distance(d):
    if d <= 0:
        return 0.0
    return math.exp(-d / DECAY_KM) / d


def load_ff_factors_daily(path):
    if not os.path.exists(path):
        return None
    vwretd = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('This file') or line.startswith('The ') or line.startswith(','):
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
            vw = (mktrf_val + rf_val) / 100.0
            date_fmt = f"{date[:4]}-{date[4:6]}-{date[6:]}"
            vwretd[date_fmt] = vw
    return vwretd if vwretd else None


def load_ff_factors_monthly(path):
    if not os.path.exists(path):
        return None
    vwretd = {}
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith('This file') or line.startswith('The ') or line.startswith(','):
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
            vw = (mktrf_val + rf_val) / 100.0
            date_fmt = f"{date[:4]}-{date[4:6]}"
            vwretd[date_fmt] = vw
    return vwretd if vwretd else None


def load_daily_returns(path):
    data = defaultdict(dict)
    with open(path, 'r', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            gvkey = row['gvkey']
            date = row['datadate']
            try:
                ret = float(row['ret_daily'])
            except (ValueError, TypeError):
                continue
            data[gvkey][date] = ret
    return data


def load_monthly_returns(path):
    data = defaultdict(dict)
    with open(path, 'r', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            gvkey = row['gvkey']
            date = row['datadate']
            try:
                ret = float(row['ret_monthly'])
            except (ValueError, TypeError):
                continue
            if len(date) >= 7:
                ym = date[:7]
            else:
                continue
            data[gvkey][ym] = ret
    return data


def compute_window_car_daily(gvkey, event_date, start, end, daily_ret, market_ret):
    if gvkey not in daily_ret or not market_ret:
        return None
    dates = sorted(daily_ret[gvkey].keys())
    if not dates:
        return None
    event_idx = None
    if event_date:
        for i, d in enumerate(dates):
            if d >= event_date:
                event_idx = i
                break
    if event_idx is None:
        return None
    ar_list = []
    for offset in range(start, end + 1):
        idx = event_idx + offset
        if 0 <= idx < len(dates):
            d = dates[idx]
            if d in market_ret:
                ar_list.append(daily_ret[gvkey][d] - market_ret[d])
    if len(ar_list) < (end - start + 1) * 0.4:
        return None
    return sum(ar_list)


def compute_window_car_monthly(gvkey, event_month, start, end, monthly_ret, market_ret):
    if gvkey not in monthly_ret or not market_ret:
        return None
    months = sorted(monthly_ret[gvkey].keys())
    if event_month not in months:
        return None
    event_idx = months.index(event_month)
    ar_list = []
    for offset in range(start, end + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months):
            m = months[idx]
            if m in market_ret:
                ar_list.append(monthly_ret[gvkey][m] - market_ret[m])
    if len(ar_list) < (end - start + 1) * 0.4:
        return None
    return sum(ar_list)


def ols_pooled(data, y_key, x_keys, cluster_keys):
    n = len(data)
    k = len(x_keys) + 1
    if n <= k + 1:
        return None
    y = [d[y_key] for d in data]
    X = [[1.0] + [d[x] for x in x_keys] for d in data]

    # XtX and Xty
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

    # invert XtX
    def invert_matrix(A):
        m = len(A)
        aug = [[A[i][j] for j in range(m)] + [1.0 if i == j else 0.0 for j in range(m)] for i in range(m)]
        for i in range(m):
            pivot = aug[i][i]
            if abs(pivot) < 1e-12:
                return None
            inv_p = 1.0 / pivot
            for j in range(2*m):
                aug[i][j] *= inv_p
            for r in range(m):
                if r == i:
                    continue
                factor = aug[r][i]
                for c in range(2*m):
                    aug[r][c] -= factor * aug[i][c]
        return [row[m:] for row in aug]

    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    beta = [sum(inv_XtX[i][j] * Xty[j] for j in range(k)) for i in range(k)]
    y_hat = [sum(beta[j] * X[i][j] for j in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]

    # cluster robust
    def cluster_cov(cluster_name):
        clusters = defaultdict(list)
        for idx, d in enumerate(data):
            clusters[d[cluster_name]].append(idx)
        S = [[0.0 for _ in range(k)] for _ in range(k)]
        for _, idxs in clusters.items():
            xu = [0.0 for _ in range(k)]
            for i in idxs:
                for a in range(k):
                    xu[a] += X[i][a] * resid[i]
            for a in range(k):
                for b in range(k):
                    S[a][b] += xu[a] * xu[b]
        return S, len(clusters)

    # single cluster (event)
    covs = {}
    for cname in cluster_keys:
        S, g = cluster_cov(cname)
        cov = [[sum(inv_XtX[i][a] * S[a][b] * inv_XtX[j][b] for a in range(k) for b in range(k))
                for j in range(k)] for i in range(k)]
        covs[cname] = (cov, g)

    return beta, covs


def main():
    # Load firm centroids
    centroids = {}
    w_sum = {}
    with open(derived_path('networks', 'firm_centroids.csv'), 'r', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            gvkey = row['gvkey']
            centroids[gvkey] = (float(row['centroid_lat']), float(row['centroid_lon']), float(row['total_mw']))
            try:
                w_sum[gvkey] = float(row['w_sum'])
            except (ValueError, TypeError):
                w_sum[gvkey] = None

    # Load gvkey->country (latest year)
    country_by_gvkey = {}
    latest_year = {}
    with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            gvkey = row['gvkey']
            try:
                year = int(row['fyear'])
            except (ValueError, TypeError):
                continue
            if gvkey not in latest_year or year > latest_year[gvkey]:
                latest_year[gvkey] = year
                country_by_gvkey[gvkey] = row['fic']

    # Country centroids (weighted by MW)
    country_acc = defaultdict(lambda: [0.0, 0.0, 0.0])
    for gvkey, (lat, lon, mw) in centroids.items():
        c = country_by_gvkey.get(gvkey)
        if not c:
            continue
        country_acc[c][0] += lat * mw
        country_acc[c][1] += lon * mw
        country_acc[c][2] += mw
    country_centroids = {}
    for c, (lat_sum, lon_sum, mw_sum) in country_acc.items():
        if mw_sum > 0:
            country_centroids[c] = (lat_sum / mw_sum, lon_sum / mw_sum)

    # Load returns + market
    daily_ret = load_daily_returns(derived_path('returns', 'daily_returns.csv'))
    monthly_ret = load_monthly_returns(derived_path('returns', 'monthly_returns.csv'))
    # Fama-French files may live in raw/fama_french or directly under finance_data
    ff_daily_path = raw_path('factors', 'F-F_Research_Data_Factors_daily.csv')
    ff_monthly_path = raw_path('factors', 'F-F_Research_Data_Factors.csv')
    if not os.path.exists(ff_daily_path):
        ff_daily_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), '..', 'F-F_Research_Data_Factors_daily.csv')
    if not os.path.exists(ff_monthly_path):
        ff_monthly_path = os.path.join(os.path.dirname(os.path.dirname(__file__)), '..', 'F-F_Research_Data_Factors.csv')
    market_daily = load_ff_factors_daily(ff_daily_path)
    market_monthly = load_ff_factors_monthly(ff_monthly_path)

    # Return coverage window
    all_months = sorted({m for g in monthly_ret.values() for m in g.keys()})
    min_month = all_months[0] if all_months else None
    max_month = all_months[-1] if all_months else None

    # Build interconnector events
    inter_path = os.path.join(os.path.dirname(__file__), '..', '..', 'interconnectors', 'cross_border_interconnectors.csv')
    inter_events = []
    with open(inter_path, 'r', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            if row.get('status', '').strip().lower() != 'operational':
                continue
            if not row.get('commissioning_year'):
                continue
            iso_a = COUNTRY_MAP.get(row['country_a'])
            iso_b = COUNTRY_MAP.get(row['country_b'])
            if not iso_a or not iso_b:
                continue
            if iso_a == iso_b:
                continue
            if iso_a not in country_centroids or iso_b not in country_centroids:
                continue
            year = int(row['commissioning_year'])
            event_date = row.get('commissioning_date') or f"{year}-07-01"
            event_month = event_date[:7]
            if min_month and (event_month < min_month or event_month > max_month):
                continue
            inter_events.append({
                'event_id': row['project_name'],
                'country_a': iso_a,
                'country_b': iso_b,
                'event_date': event_date,
                'event_month': event_month,
                'lat_a': country_centroids[iso_a][0],
                'lon_a': country_centroids[iso_a][1],
                'lat_b': country_centroids[iso_b][0],
                'lon_b': country_centroids[iso_b][1],
            })

    # Build retirement events (first movers) with plant lat/lon
    retire_events = []
    with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
        r = csv.DictReader(f)
        for row in r:
            if row.get('is_first_mover', '').lower() != 'true':
                continue
            try:
                lat = float(row['lat'])
                lon = float(row['lon'])
            except (ValueError, TypeError):
                continue
            date = row.get('announcement_date') or row.get('event_date') or ''
            if not date:
                try:
                    year = int(row.get('ret_year'))
                    date = f"{year}-07-01"
                except (ValueError, TypeError):
                    continue
            event_month = date[:7]
            if min_month and (event_month < min_month or event_month > max_month):
                continue
            retire_events.append({
                'event_id': row['gem_id'],
                'event_date': date,
                'event_month': event_month,
                'lat': lat,
                'lon': lon,
            })

    # Exposure mapping + CARs for interconnectors
    out_dir = results_path('interconnectors')
    os.makedirs(out_dir, exist_ok=True)

    mapping_rows = []
    inter_obs = []  # firm-event obs for interconnectors

    density_vals = [v for v in w_sum.values() if v is not None]
    density_median = sorted(density_vals)[len(density_vals)//2] if density_vals else None

    for ev in inter_events:
        exposures = {}
        near_counts = {th: 0 for th in NEAR_THRESHOLDS}
        for gvkey, (lat, lon, mw) in centroids.items():
            d_a = haversine(lat, lon, ev['lat_a'], ev['lon_a'])
            d_b = haversine(lat, lon, ev['lat_b'], ev['lon_b'])
            d_min = min(d_a, d_b)
            for th in NEAR_THRESHOLDS:
                if d_min <= th:
                    near_counts[th] += 1
            w = weight_from_distance(d_a) + weight_from_distance(d_b)
            exposures[gvkey] = w

        vals = list(exposures.values())
        vals_sorted = sorted(vals)
        p50 = vals_sorted[len(vals_sorted)//2] if vals_sorted else 0.0

        mapping_rows.append({
            'event_id': ev['event_id'],
            'n_firms': len(exposures),
            'n_near_500km': near_counts[500],
            'n_near_1000km': near_counts[1000],
            'exposure_min': min(vals) if vals else 0.0,
            'exposure_median': p50,
            'exposure_max': max(vals) if vals else 0.0,
        })

        for gvkey, w in exposures.items():
            # CARs
            car_d = compute_window_car_daily(gvkey, ev['event_date'], DAILY_WINDOW[0], DAILY_WINDOW[1], daily_ret, market_daily)
            car_pre = compute_window_car_monthly(gvkey, ev['event_month'], PRE_MONTHS[0], PRE_MONTHS[1], monthly_ret, market_monthly)
            cars_post = {}
            for post in POST_MONTHS:
                cars_post[post] = compute_window_car_monthly(gvkey, ev['event_month'], -1, post, monthly_ret, market_monthly)
            inter_obs.append({
                'event_id': ev['event_id'],
                'gvkey': gvkey,
                'exposure': w,
                'car_d': car_d,
                'car_pre': car_pre,
                'car_m3': cars_post[3],
                'car_m6': cars_post[6],
                'car_m12': cars_post[12],
                'density': w_sum.get(gvkey),
            })

    # Write mapping
    mapping_path = os.path.join(out_dir, 'interconnector_exposure_mapping.csv')
    with open(mapping_path, 'w', newline='', encoding='utf-8') as f:
        w = csv.DictWriter(f, fieldnames=list(mapping_rows[0].keys())) if mapping_rows else None
        if w:
            w.writeheader(); w.writerows(mapping_rows)

    # Pre-trend + raw CAR diff summaries (event-level high vs low exposure)
    def event_level_diff(obs, car_key):
        by_event = defaultdict(list)
        for o in obs:
            if o[car_key] is None:
                continue
            by_event[o['event_id']].append(o)
        diffs = []
        for eid, rows in by_event.items():
            exps = sorted(r['exposure'] for r in rows)
            if not exps:
                continue
            med = exps[len(exps)//2]
            high = [r[car_key] for r in rows if r['exposure'] >= med]
            low = [r[car_key] for r in rows if r['exposure'] < med]
            if not high or not low:
                continue
            diffs.append(sum(high)/len(high) - sum(low)/len(low))
        if not diffs:
            return None
        mean = sum(diffs)/len(diffs)
        var = sum((d-mean)**2 for d in diffs)/max(len(diffs)-1,1)
        se = math.sqrt(var/len(diffs))
        t = mean/se if se>0 else 0
        return mean, se, t, len(diffs)

    summary = []
    pre = event_level_diff(inter_obs, 'car_pre')
    summary.append(('Pre-trend CAR[-6,-1] (monthly)', pre))
    for post_key, label in [('car_m3', 'CAR[-1,+3] (monthly)'), ('car_m6','CAR[-1,+6] (monthly)'), ('car_m12','CAR[-1,+12] (monthly)')]:
        summary.append((label, event_level_diff(inter_obs, post_key)))
    summary.append(('CAR[-1,+20] (daily)', event_level_diff(inter_obs, 'car_d')))

    # Coverage diagnostics
    def count_nonnull(obs, key):
        return sum(1 for o in obs if o.get(key) is not None)
    inter_counts = {
        'car_pre': count_nonnull(inter_obs, 'car_pre'),
        'car_m3': count_nonnull(inter_obs, 'car_m3'),
        'car_m6': count_nonnull(inter_obs, 'car_m6'),
        'car_m12': count_nonnull(inter_obs, 'car_m12'),
        'car_d': count_nonnull(inter_obs, 'car_d'),
    }

    # Density-group comparison (interconnectors)
    if density_median is not None:
        by_event = defaultdict(list)
        for o in inter_obs:
            if o['car_m12'] is None or o['density'] is None:
                continue
            by_event[o['event_id']].append(o)
        diffs = []
        for eid, rows in by_event.items():
            high = [r['car_m12'] for r in rows if r['density'] >= density_median]
            low = [r['car_m12'] for r in rows if r['density'] < density_median]
            if high and low:
                diffs.append(sum(high)/len(high) - sum(low)/len(low))
        if diffs:
            mean = sum(diffs)/len(diffs)
            var = sum((d-mean)**2 for d in diffs)/max(len(diffs)-1,1)
            se = math.sqrt(var/len(diffs))
            t = mean/se if se>0 else 0
            density_summary = (mean, se, t, len(diffs))
        else:
            density_summary = None
    else:
        density_summary = None

    # Overlap diagnostic: firms exposed to any interconnector vs any retirement
    inter_firms = {o['gvkey'] for o in inter_obs if o['exposure'] is not None}

    # Build retirement obs (exposure to plant location)
    retire_obs = []
    for ev in retire_events:
        exposures = {}
        for gvkey, (lat, lon, mw) in centroids.items():
            d = haversine(lat, lon, ev['lat'], ev['lon'])
            w = weight_from_distance(d)
            exposures[gvkey] = w
        for gvkey, w in exposures.items():
            car_m12 = compute_window_car_monthly(gvkey, ev['event_month'], -1, 12, monthly_ret, market_monthly)
            if car_m12 is None:
                continue
            retire_obs.append({
                'event_id': ev['event_id'],
                'gvkey': gvkey,
                'exposure': w,
                'car_m12': car_m12,
            })

    retire_firms = {o['gvkey'] for o in retire_obs}
    overlap = len(inter_firms & retire_firms)

    # Pooled regression: CAR_m12 ~ exposure_z + shock + exposure_z*shock
    pooled = []
    def add_pooled(obs, shock_type):
        by_event = defaultdict(list)
        for o in obs:
            by_event[o['event_id']].append(o)
        for eid, rows in by_event.items():
            exps = [r['exposure'] for r in rows]
            mean = sum(exps)/len(exps)
            var = sum((x-mean)**2 for x in exps)/max(len(exps)-1,1)
            sd = math.sqrt(var) if var>0 else 1.0
            for r in rows:
                ez = (r['exposure'] - mean) / sd if sd>0 else 0.0
                pooled.append({
                    'event_id': eid,
                    'gvkey': r['gvkey'],
                    'car': r['car_m12'],
                    'exposure_z': ez,
                    'shock': shock_type,
                    'exposure_z_x_shock': ez * shock_type,
                })
    # use interconnector obs with car_m12
    add_pooled([o for o in inter_obs if o['car_m12'] is not None], 1)
    add_pooled(retire_obs, 0)

    pooled_result = ols_pooled(pooled, 'car', ['exposure_z', 'shock', 'exposure_z_x_shock'], ['event_id'])

    # Write summary
    summary_path = os.path.join(out_dir, 'interconnector_eventstudy_summary.md')
    with open(summary_path, 'w', encoding='utf-8') as f:
        f.write('# Interconnector Event Study Summary\n\n')
        f.write(f'Interconnector events (operational, cross-border) used: {len(inter_events)}\n')
        f.write(f'Retirement first-mover events used: {len(retire_events)}\n')
        f.write(f'Overlap firms (interconnector-exposed vs retirement-exposed): {overlap}\n\n')

        f.write('## Exposure mapping (per event)\n')
        f.write(f'Output: {mapping_path}\n\n')

        f.write('## Pre-trends and raw CAR differences (high vs low exposure, event-level)\n')
        for label, res in summary:
            if res is None:
                f.write(f'- {label}: insufficient data\n')
            else:
                mean, se, t, n = res
                f.write(f'- {label}: mean diff={mean:.4f}, se={se:.4f}, t={t:.2f}, N_events={n}\n')

        f.write('\n## Coverage diagnostics (interconnector firm-event obs)\n')
        f.write(f'- car_pre: {inter_counts["car_pre"]}\n')
        f.write(f'- car_m3: {inter_counts["car_m3"]}\n')
        f.write(f'- car_m6: {inter_counts["car_m6"]}\n')
        f.write(f'- car_m12: {inter_counts["car_m12"]}\n')
        f.write(f'- car_d: {inter_counts["car_d"]}\n')
        f.write(f'- retirement obs (car_m12): {len(retire_obs)}\n')

        if density_summary:
            mean, se, t, n = density_summary
            f.write('\n## Density-group (high vs low density, monthly +12)\n')
            f.write(f'- mean diff={mean:.4f}, se={se:.4f}, t={t:.2f}, N_events={n}\n')

        f.write('\n## Pooled regression (CAR12) interconnector vs retirement\n')
        if pooled_result is None:
            f.write('Regression failed (insufficient data).\n')
        else:
            beta, covs = pooled_result
            se_event = [math.sqrt(covs['event_id'][0][i][i]) for i in range(len(beta))]
            f.write('Spec: CAR ~ exposure_z + shock + exposure_z*shock\n')
            f.write(f'  beta_const={beta[0]:.4f}\n')
            f.write(f'  beta_exposure_z={beta[1]:.4f} (se={se_event[1]:.4f}, t={beta[1]/se_event[1]:.2f})\n')
            f.write(f'  beta_shock={beta[2]:.4f} (se={se_event[2]:.4f}, t={beta[2]/se_event[2]:.2f})\n')
            f.write(f'  beta_exposure_z_x_shock={beta[3]:.4f} (se={se_event[3]:.4f}, t={beta[3]/se_event[3]:.2f})\n')

    print(f'Wrote: {summary_path}')
    print(f'Wrote: {mapping_path}')


if __name__ == '__main__':
    main()
