"""Interconnector diagnostics: extensive pre-tests and robustness.

This script runs a full suite of diagnostics for the interconnector experiment:
  - Exposure mapping (kernel to endpoints)
  - Pre-trends (event-level high vs low exposure)
  - Raw CAR differences for multiple windows
  - Return-model sensitivity (vwretd, constant-mean, CAPM)
  - Exposure transforms (raw, zscore, log1p)
  - Foreign-only exposure variant
  - Placebo time shifts (+/- 6, 12 months)
  - Leave-one-event-out stability
  - Event-time profiles (tau=-6..+12)
  - Pooled regression with retirement events (interconnector vs retirement)

Outputs (results/interconnectors/):
  - interconnector_diagnostics_summary.md
  - interconnector_event_time_<model>_<variant>.csv
  - interconnector_leave_one_out.csv
"""
import csv
import math
import os
from collections import defaultdict

from _paths import derived_path, results_path, raw_path


# --- Configuration ---
DECAY_KM = 500 / math.log(2)
PRE_DAYS = 250
PRE_MONTHS = 24
MONTH_WINDOWS = [3, 6, 12]
DAILY_WINDOW = (-1, 20)
PLACEBO_SHIFTS = [-12, -6, 6, 12]
EVENT_TAU = list(range(-6, 13))

RETURN_MODELS = ["vwretd", "constant_mean", "capm"]
EXPOSURE_TRANSFORMS = ["raw", "zscore", "log1p"]
EXPOSURE_VARIANTS = ["base", "foreign"]


COUNTRY_MAP = {
    "Australia (Tasmania)": "AUS",
    "Australia (Victoria)": "AUS",
    "Bangladesh": "BGD",
    "Belgium": "BEL",
    "Bhutan": "BTN",
    "Denmark": "DNK",
    "Estonia": "EST",
    "Finland": "FIN",
    "France": "FRA",
    "Germany": "DEU",
    "Greece (Crete)": "GRC",
    "Greece (mainland)": "GRC",
    "India": "IND",
    "Indonesia (W.Kalimantan)": "IDN",
    "Ireland": "IRL",
    "Italy": "ITA",
    "Italy (Sardinia)": "ITA",
    "Italy (mainland)": "ITA",
    "Japan (Hokkaido)": "JPN",
    "Japan (Honshu)": "JPN",
    "Laos": "LAO",
    "Lithuania": "LTU",
    "Malaysia (Sarawak)": "MYS",
    "Montenegro": "MNE",
    "Myanmar": "MMR",
    "Nepal": "NPL",
    "Netherlands": "NLD",
    "Norway": "NOR",
    "Poland": "POL",
    "Spain": "ESP",
    "Sweden": "SWE",
    "UK": "GBR",
    "UK (Wales)": "GBR",
    # ambiguous multi-country endpoint; skip
    "Thailand-Malaysia-Singapore": None,
}


def haversine(lat1, lon1, lat2, lon2):
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat / 2) ** 2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon / 2) ** 2
    return R * 2 * math.asin(math.sqrt(a))


def weight_from_distance(d):
    if d <= 0:
        return 0.0
    return math.exp(-d / DECAY_KM) / d


def load_ff_factors_daily(path):
    if not os.path.exists(path):
        return None
    vwretd = {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("This file") or line.startswith("The ") or line.startswith(","):
                continue
            parts = [p.strip() for p in line.split(",")]
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
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("This file") or line.startswith("The ") or line.startswith(","):
                continue
            parts = [p.strip() for p in line.split(",")]
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
    with open(path, "r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            gvkey = row["gvkey"]
            date = row["datadate"]
            try:
                ret = float(row["ret_daily"])
            except (ValueError, TypeError):
                continue
            data[gvkey][date] = ret
    return data


def load_monthly_returns(path):
    data = defaultdict(dict)
    with open(path, "r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            gvkey = row["gvkey"]
            date = row["datadate"]
            try:
                ret = float(row["ret_monthly"])
            except (ValueError, TypeError):
                continue
            if len(date) >= 7:
                ym = date[:7]
            else:
                continue
            data[gvkey][ym] = ret
    return data


def compute_ar_series_monthly(gvkey, event_month, monthly_ret, market_ret, model):
    if gvkey not in monthly_ret or not market_ret:
        return None, None
    months = sorted(monthly_ret[gvkey].keys())
    if event_month not in months:
        return None, None
    idx = months.index(event_month)

    pre_start = max(0, idx - PRE_MONTHS)
    pre_months = months[pre_start:idx]
    if len(pre_months) < max(6, int(PRE_MONTHS * 0.4)):
        return None, None

    ret_pre = [monthly_ret[gvkey][m] for m in pre_months if m in market_ret]
    mkt_pre = [market_ret[m] for m in pre_months if m in market_ret]
    if len(ret_pre) < max(6, int(PRE_MONTHS * 0.4)):
        return None, None

    if model == "constant_mean":
        mean_ret = sum(ret_pre) / len(ret_pre)
        ar = {m: monthly_ret[gvkey][m] - mean_ret for m in months if m in market_ret}
        return ar, None
    if model == "vwretd":
        ar = {m: monthly_ret[gvkey][m] - market_ret[m] for m in months if m in market_ret}
        return ar, None
    if model == "capm":
        mean_ret = sum(ret_pre) / len(ret_pre)
        mean_mkt = sum(mkt_pre) / len(mkt_pre)
        cov = sum((ri - mean_ret) * (mi - mean_mkt) for ri, mi in zip(ret_pre, mkt_pre)) / max(len(ret_pre) - 1, 1)
        var = sum((mi - mean_mkt) ** 2 for mi in mkt_pre) / max(len(mkt_pre) - 1, 1)
        beta = cov / var if var > 0 else 0.0
        alpha = mean_ret - beta * mean_mkt
        ar = {m: monthly_ret[gvkey][m] - (alpha + beta * market_ret[m]) for m in months if m in market_ret}
        return ar, (alpha, beta)
    return None, None


def compute_ar_series_daily(gvkey, event_date, daily_ret, market_ret, model):
    if gvkey not in daily_ret or not market_ret:
        return None, None
    dates = sorted(daily_ret[gvkey].keys())
    if not dates:
        return None, None
    event_idx = None
    for i, d in enumerate(dates):
        if d >= event_date:
            event_idx = i
            break
    if event_idx is None:
        return None, None

    pre_start = max(0, event_idx - PRE_DAYS)
    pre_dates = dates[pre_start:event_idx]
    ret_pre = [daily_ret[gvkey][d] for d in pre_dates if d in market_ret]
    mkt_pre = [market_ret[d] for d in pre_dates if d in market_ret]
    if len(ret_pre) < max(50, int(PRE_DAYS * 0.4)):
        return None, None

    if model == "constant_mean":
        mean_ret = sum(ret_pre) / len(ret_pre)
        ar = {d: daily_ret[gvkey][d] - mean_ret for d in dates if d in market_ret}
        return ar, None
    if model == "vwretd":
        ar = {d: daily_ret[gvkey][d] - market_ret[d] for d in dates if d in market_ret}
        return ar, None
    if model == "capm":
        mean_ret = sum(ret_pre) / len(ret_pre)
        mean_mkt = sum(mkt_pre) / len(mkt_pre)
        cov = sum((ri - mean_ret) * (mi - mean_mkt) for ri, mi in zip(ret_pre, mkt_pre)) / max(len(ret_pre) - 1, 1)
        var = sum((mi - mean_mkt) ** 2 for mi in mkt_pre) / max(len(mkt_pre) - 1, 1)
        beta = cov / var if var > 0 else 0.0
        alpha = mean_ret - beta * mean_mkt
        ar = {d: daily_ret[gvkey][d] - (alpha + beta * market_ret[d]) for d in dates if d in market_ret}
        return ar, (alpha, beta)
    return None, None


def compute_car_from_ar(ar_map, event_key, window, is_monthly=True):
    if ar_map is None:
        return None
    keys = sorted(ar_map.keys())
    if event_key not in keys:
        return None
    idx = keys.index(event_key)
    start, end = window
    vals = []
    for offset in range(start, end + 1):
        j = idx + offset
        if 0 <= j < len(keys):
            k = keys[j]
            vals.append(ar_map[k])
    if len(vals) < (end - start + 1) * 0.4:
        return None
    return sum(vals)


def transform_exposure(values, transform):
    if transform == "raw":
        return values
    if transform == "log1p":
        return [math.log1p(v) for v in values]
    if transform == "zscore":
        mean = sum(values) / len(values)
        var = sum((v - mean) ** 2 for v in values) / max(len(values) - 1, 1)
        sd = math.sqrt(var) if var > 0 else 1.0
        return [(v - mean) / sd for v in values]
    return values


def event_level_diff(rows, car_key, exposure_key):
    by_event = defaultdict(list)
    for r in rows:
        if r.get(car_key) is None:
            continue
        by_event[r["event_id"]].append(r)
    diffs = []
    for eid, obs in by_event.items():
        vals = [o[exposure_key] for o in obs]
        if not vals:
            continue
        vals_sorted = sorted(vals)
        med = vals_sorted[len(vals_sorted) // 2]
        high = [o[car_key] for o in obs if o[exposure_key] >= med]
        low = [o[car_key] for o in obs if o[exposure_key] < med]
        if not high or not low:
            continue
        diffs.append(sum(high) / len(high) - sum(low) / len(low))
    if not diffs:
        return None
    mean = sum(diffs) / len(diffs)
    var = sum((d - mean) ** 2 for d in diffs) / max(len(diffs) - 1, 1)
    se = math.sqrt(var / len(diffs)) if len(diffs) > 0 else 0.0
    t = mean / se if se > 0 else 0.0
    return mean, se, t, len(diffs)


def ols_pooled(data, y_key, x_keys, cluster_key="event_id"):
    n = len(data)
    k = len(x_keys) + 1
    if n <= k + 1:
        return None
    y = [d[y_key] for d in data]
    X = [[1.0] + [d[x] for x in x_keys] for d in data]

    XtX = [[sum(X[i][a] * X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]

    def invert_matrix(A):
        m = len(A)
        aug = [[A[i][j] for j in range(m)] + [1.0 if i == j else 0.0 for j in range(m)] for i in range(m)]
        for i in range(m):
            pivot = aug[i][i]
            if abs(pivot) < 1e-12:
                return None
            inv_p = 1.0 / pivot
            for j in range(2 * m):
                aug[i][j] *= inv_p
            for r in range(m):
                if r == i:
                    continue
                factor = aug[r][i]
                for c in range(2 * m):
                    aug[r][c] -= factor * aug[i][c]
        return [row[m:] for row in aug]

    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    beta = [sum(inv_XtX[i][j] * Xty[j] for j in range(k)) for i in range(k)]
    y_hat = [sum(beta[j] * X[i][j] for j in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]

    clusters = defaultdict(list)
    for idx, d in enumerate(data):
        clusters[d[cluster_key]].append(idx)

    S = [[0.0 for _ in range(k)] for _ in range(k)]
    for _, idxs in clusters.items():
        xu = [0.0 for _ in range(k)]
        for i in idxs:
            for a in range(k):
                xu[a] += X[i][a] * resid[i]
        for a in range(k):
            for b in range(k):
                S[a][b] += xu[a] * xu[b]

    cov = [[sum(inv_XtX[i][a] * S[a][b] * inv_XtX[j][b] for a in range(k) for b in range(k))
            for j in range(k)] for i in range(k)]
    se = [math.sqrt(cov[i][i]) if cov[i][i] > 0 else float("nan") for i in range(k)]
    return beta, se, len(clusters)


def main():
    out_dir = results_path("interconnectors")
    os.makedirs(out_dir, exist_ok=True)

    # Load centroids
    centroids = {}
    w_sum = {}
    with open(derived_path("networks", "firm_centroids.csv"), "r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            gvkey = row["gvkey"]
            centroids[gvkey] = (float(row["centroid_lat"]), float(row["centroid_lon"]), float(row["total_mw"]))
            try:
                w_sum[gvkey] = float(row["w_sum"])
            except (ValueError, TypeError):
                w_sum[gvkey] = None

    # gvkey -> country (latest year)
    country_by_gvkey = {}
    latest_year = {}
    with open(derived_path("fundamentals", "firm_fundamentals.csv"), "r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            gvkey = row["gvkey"]
            try:
                year = int(row["fyear"])
            except (ValueError, TypeError):
                continue
            if gvkey not in latest_year or year > latest_year[gvkey]:
                latest_year[gvkey] = year
                country_by_gvkey[gvkey] = row["fic"]

    # Country centroids (MW-weighted)
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

    # Returns
    daily_ret = load_daily_returns(derived_path("returns", "daily_returns.csv"))
    monthly_ret = load_monthly_returns(derived_path("returns", "monthly_returns.csv"))

    ff_daily = raw_path("factors", "F-F_Research_Data_Factors_daily.csv")
    ff_monthly = raw_path("factors", "F-F_Research_Data_Factors.csv")
    market_daily = load_ff_factors_daily(ff_daily)
    market_monthly = load_ff_factors_monthly(ff_monthly)

    # Return coverage window
    all_months = sorted({m for g in monthly_ret.values() for m in g.keys()})
    min_month = all_months[0] if all_months else None
    max_month = all_months[-1] if all_months else None

    # Build interconnector events
    inter_path = os.path.join(os.path.dirname(__file__), "..", "..", "interconnectors", "cross_border_interconnectors.csv")
    inter_events = []
    with open(inter_path, "r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            if row.get("status", "").strip().lower() != "operational":
                continue
            if not row.get("commissioning_year"):
                continue
            iso_a = COUNTRY_MAP.get(row["country_a"])
            iso_b = COUNTRY_MAP.get(row["country_b"])
            if not iso_a or not iso_b:
                continue
            if iso_a == iso_b:
                continue
            if iso_a not in country_centroids or iso_b not in country_centroids:
                continue
            year = int(row["commissioning_year"])
            # Use announcement/FID date when available; fall back to commissioning date or mid-year.
            ann = (row.get("announcement_or_fid_date") or "").strip()
            event_date = ann or (row.get("commissioning_date") or f"{year}-07-01")
            # If month-only, use mid-month convention.
            if len(event_date) == 7 and event_date[4] == "-":
                event_date = f"{event_date}-15"
            event_month = event_date[:7]
            if min_month and (event_month < min_month or event_month > max_month):
                continue
            inter_events.append({
                "event_id": row["project_name"],
                "country_a": iso_a,
                "country_b": iso_b,
                "event_date": event_date,
                "event_month": event_month,
                "lat_a": country_centroids[iso_a][0],
                "lon_a": country_centroids[iso_a][1],
                "lat_b": country_centroids[iso_b][0],
                "lon_b": country_centroids[iso_b][1],
            })

    # Retirements (first mover)
    retire_events = []
    with open(derived_path("events", "coal_retirement_events.csv"), "r", encoding="utf-8") as f:
        r = csv.DictReader(f)
        for row in r:
            if row.get("is_first_mover", "").lower() != "true":
                continue
            try:
                lat = float(row["lat"])
                lon = float(row["lon"])
            except (ValueError, TypeError):
                continue
            date = row.get("announcement_date") or row.get("event_date") or ""
            if not date:
                try:
                    year = int(row.get("ret_year"))
                    date = f"{year}-07-01"
                except (ValueError, TypeError):
                    continue
            event_month = date[:7]
            if min_month and (event_month < min_month or event_month > max_month):
                continue
            retire_events.append({
                "event_id": row["gem_id"],
                "event_date": date,
                "event_month": event_month,
                "lat": lat,
                "lon": lon,
            })

    # Precompute exposures for interconnectors
    inter_obs = []
    mapping_rows = []
    density_vals = [v for v in w_sum.values() if v is not None]
    density_median = sorted(density_vals)[len(density_vals) // 2] if density_vals else None

    for ev in inter_events:
        exposures = {}
        for gvkey, (lat, lon, mw) in centroids.items():
            d_a = haversine(lat, lon, ev["lat_a"], ev["lon_a"])
            d_b = haversine(lat, lon, ev["lat_b"], ev["lon_b"])
            w_a = weight_from_distance(d_a)
            w_b = weight_from_distance(d_b)
            exposures[gvkey] = (w_a + w_b, w_a, w_b)

        vals = [v[0] for v in exposures.values()]
        vals_sorted = sorted(vals)
        p50 = vals_sorted[len(vals_sorted) // 2] if vals_sorted else 0.0

        mapping_rows.append({
            "event_id": ev["event_id"],
            "n_firms": len(exposures),
            "exposure_min": min(vals) if vals else 0.0,
            "exposure_median": p50,
            "exposure_max": max(vals) if vals else 0.0,
        })

        for gvkey, (w_tot, w_a, w_b) in exposures.items():
            firm_country = country_by_gvkey.get(gvkey)
            if firm_country == ev["country_a"]:
                w_foreign = w_b
            elif firm_country == ev["country_b"]:
                w_foreign = w_a
            else:
                w_foreign = w_tot
            inter_obs.append({
                "event_id": ev["event_id"],
                "gvkey": gvkey,
                "country": firm_country,
                "event_month": ev["event_month"],
                "event_date": ev["event_date"],
                "exposure_raw": w_tot,
                "exposure_foreign_raw": w_foreign,
                "density": w_sum.get(gvkey),
            })

    # Write mapping
    mapping_path = os.path.join(out_dir, "interconnector_exposure_mapping.csv")
    with open(mapping_path, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=list(mapping_rows[0].keys())) if mapping_rows else None
        if w:
            w.writeheader()
            w.writerows(mapping_rows)

    # Compute AR series and CARs
    obs_by_model = {m: [] for m in RETURN_MODELS}

    for ev in inter_events:
        for gvkey in centroids.keys():
            row = next((o for o in inter_obs if o["event_id"] == ev["event_id"] and o["gvkey"] == gvkey), None)
            if row is None:
                continue
            for model in RETURN_MODELS:
                ar_m, _ = compute_ar_series_monthly(gvkey, ev["event_month"], monthly_ret, market_monthly, model)
                ar_d, _ = compute_ar_series_daily(gvkey, ev["event_date"], daily_ret, market_daily, model)
                if ar_m is None and ar_d is None:
                    continue
                out = dict(row)
                out["model"] = model
                out["car_pre"] = compute_car_from_ar(ar_m, ev["event_month"], (-6, -1), True) if ar_m else None
                for post in MONTH_WINDOWS:
                    out[f"car_m{post}"] = compute_car_from_ar(ar_m, ev["event_month"], (-1, post), True) if ar_m else None
                out["car_d"] = compute_car_from_ar(ar_d, ev["event_date"], DAILY_WINDOW, False) if ar_d else None
                obs_by_model[model].append(out)

    # Build transformed exposure variants (per event)
    for model, rows in obs_by_model.items():
        by_event = defaultdict(list)
        for r in rows:
            by_event[r["event_id"]].append(r)
        for _, obs in by_event.items():
            base_vals = [o["exposure_raw"] for o in obs]
            foreign_vals = [o["exposure_foreign_raw"] for o in obs]
            for t in EXPOSURE_TRANSFORMS:
                base_t = transform_exposure(base_vals, t)
                foreign_t = transform_exposure(foreign_vals, t)
                for i, o in enumerate(obs):
                    o[f"exposure_{t}"] = base_t[i]
                    o[f"exposure_foreign_{t}"] = foreign_t[i]

    # Diagnostics summary
    summary_path = os.path.join(out_dir, "interconnector_diagnostics_summary.md")
    with open(summary_path, "w", encoding="utf-8") as f:
        f.write("# Interconnector Diagnostics Summary\n\n")
        f.write(f"Interconnector events used: {len(inter_events)}\n")
        f.write(f"Retirement first-mover events used: {len(retire_events)}\n\n")
        f.write(f"Exposure mapping: {mapping_path}\n\n")

        for model in RETURN_MODELS:
            f.write(f"## Return model: {model}\n")
            rows = obs_by_model[model]
            if not rows:
                f.write("No observations.\n\n")
                continue
            for variant in EXPOSURE_VARIANTS:
                f.write(f"### Exposure variant: {variant}\n")
                for t in EXPOSURE_TRANSFORMS:
                    key = f"exposure_{t}" if variant == "base" else f"exposure_foreign_{t}"
                    f.write(f"- Transform: {t}\n")
                    pre = event_level_diff(rows, "car_pre", key)
                    if pre:
                        mean, se, tt, n = pre
                        f.write(f"  - Pre-trend CAR[-6,-1]: mean={mean:.4f}, se={se:.4f}, t={tt:.2f}, N_events={n}\n")
                    else:
                        f.write("  - Pre-trend CAR[-6,-1]: insufficient\n")
                    for post in MONTH_WINDOWS:
                        res = event_level_diff(rows, f"car_m{post}", key)
                        if res:
                            mean, se, tt, n = res
                            f.write(f"  - CAR[-1,+{post}] (monthly): mean={mean:.4f}, se={se:.4f}, t={tt:.2f}, N_events={n}\n")
                        else:
                            f.write(f"  - CAR[-1,+{post}] (monthly): insufficient\n")
                    res_d = event_level_diff(rows, "car_d", key)
                    if res_d:
                        mean, se, tt, n = res_d
                        f.write(f"  - CAR[-1,+20] (daily): mean={mean:.4f}, se={se:.4f}, t={tt:.2f}, N_events={n}\n")
                    else:
                        f.write("  - CAR[-1,+20] (daily): insufficient\n")
                f.write("\n")
            f.write("\n")

        # Placebo shifts (base exposure, vwretd)
        f.write("## Placebo time shifts (vwretd, base exposure)\n")
        base_rows = obs_by_model.get("vwretd", [])
        for shift in PLACEBO_SHIFTS:
            placebo = []
            for r in base_rows:
                em = r["event_month"]
                try:
                    y, m = em.split("-")
                    y, m = int(y), int(m)
                except Exception:
                    continue
                m_new = m + shift
                y_new = y + (m_new - 1) // 12
                m_new = ((m_new - 1) % 12) + 1
                em_new = f"{y_new:04d}-{m_new:02d}"
                ar_m, _ = compute_ar_series_monthly(r["gvkey"], em_new, monthly_ret, market_monthly, "vwretd")
                if ar_m is None:
                    continue
                out = dict(r)
                out["car_m12"] = compute_car_from_ar(ar_m, em_new, (-1, 12), True)
                if out["car_m12"] is not None:
                    placebo.append(out)
            res = event_level_diff(placebo, "car_m12", "exposure_raw")
            if res:
                mean, se, tt, n = res
                f.write(f"- Shift {shift:+} months: mean={mean:.4f}, se={se:.4f}, t={tt:.2f}, N_events={n}\n")
            else:
                f.write(f"- Shift {shift:+} months: insufficient\n")
        f.write("\n")

        # Leave-one-out (vwretd, base, CAR12)
        f.write("## Leave-one-event-out (vwretd, base, CAR12)\n")
        loo_rows = []
        by_event = defaultdict(list)
        for r in base_rows:
            if r.get("car_m12") is not None:
                by_event[r["event_id"]].append(r)
        for eid in by_event.keys():
            rows_ = [r for r in base_rows if r["event_id"] != eid]
            res = event_level_diff(rows_, "car_m12", "exposure_raw")
            if res:
                mean, se, tt, n = res
                loo_rows.append({"event_id": eid, "mean": mean, "se": se, "t": tt, "N_events": n})
        loo_path = os.path.join(out_dir, "interconnector_leave_one_out.csv")
        if loo_rows:
            with open(loo_path, "w", newline="", encoding="utf-8") as fp:
                w = csv.DictWriter(fp, fieldnames=loo_rows[0].keys())
                w.writeheader()
                w.writerows(loo_rows)
            f.write(f"Leave-one-out written: {loo_path}\n\n")

        # Event-time profile (vwretd, base + foreign)
        for variant in ["base", "foreign"]:
            rows = base_rows
            key = "exposure_raw" if variant == "base" else "exposure_foreign_raw"
            by_event = defaultdict(list)
            for r in rows:
                by_event[r["event_id"]].append(r)
            out_rows = []
            for tau in EVENT_TAU:
                diffs = []
                for _, obs in by_event.items():
                    vals = [o[key] for o in obs]
                    if not vals:
                        continue
                    med = sorted(vals)[len(vals)//2]
                    for o in obs:
                        ar_m, _ = compute_ar_series_monthly(o["gvkey"], o["event_month"], monthly_ret, market_monthly, "vwretd")
                        if ar_m is None:
                            continue
                        car_tau = compute_car_from_ar(ar_m, o["event_month"], (tau, tau), True)
                        if car_tau is None:
                            continue
                        if o[key] >= med:
                            diffs.append(("high", car_tau))
                        else:
                            diffs.append(("low", car_tau))
                if not diffs:
                    continue
                high = [v for g, v in diffs if g == "high"]
                low = [v for g, v in diffs if g == "low"]
                if not high or not low:
                    continue
                mean = sum(high) / len(high) - sum(low) / len(low)
                out_rows.append({"tau": tau, "diff_high_low": mean})
            out_path = os.path.join(out_dir, f"interconnector_event_time_vwretd_{variant}.csv")
            with open(out_path, "w", newline="", encoding="utf-8") as fp:
                w = csv.DictWriter(fp, fieldnames=["tau", "diff_high_low"])
                w.writeheader()
                w.writerows(out_rows)
            f.write(f"Event-time profile written: {out_path}\n\n")

        # Pooled regression with retirements (vwretd, base exposure z)
        pooled = []
        def add_pooled(obs, shock_type):
            by_event = defaultdict(list)
            for o in obs:
                if o.get("car_m12") is None:
                    continue
                by_event[o["event_id"]].append(o)
            for eid, rows in by_event.items():
                exps = [r["exposure_raw"] for r in rows]
                mean = sum(exps) / len(exps)
                var = sum((x - mean) ** 2 for x in exps) / max(len(exps) - 1, 1)
                sd = math.sqrt(var) if var > 0 else 1.0
                for r in rows:
                    ez = (r["exposure_raw"] - mean) / sd if sd > 0 else 0.0
                    pooled.append({
                        "event_id": eid,
                        "car": r["car_m12"],
                        "exposure_z": ez,
                        "shock": shock_type,
                        "exposure_z_x_shock": ez * shock_type,
                    })
        add_pooled(base_rows, 1)
        retire_obs = []
        for ev in retire_events:
            for gvkey, (lat, lon, mw) in centroids.items():
                d = haversine(lat, lon, ev["lat"], ev["lon"])
                w = weight_from_distance(d)
                ar_m, _ = compute_ar_series_monthly(gvkey, ev["event_month"], monthly_ret, market_monthly, "vwretd")
                car12 = compute_car_from_ar(ar_m, ev["event_month"], (-1, 12), True) if ar_m else None
                if car12 is None:
                    continue
                retire_obs.append({"event_id": ev["event_id"], "exposure_raw": w, "car_m12": car12})
        add_pooled(retire_obs, 0)

        f.write("## Pooled regression (interconnector vs retirement, vwretd, CAR12)\n")
        pooled_res = ols_pooled(pooled, "car", ["exposure_z", "shock", "exposure_z_x_shock"], "event_id")
        if pooled_res is None:
            f.write("Regression failed (insufficient data).\n")
        else:
            beta, se, g = pooled_res
            f.write(f"Clusters (events): {g}\n")
            f.write(f"beta_exposure_z={beta[1]:.4f} (se={se[1]:.4f}, t={beta[1]/se[1]:.2f})\n")
            f.write(f"beta_shock={beta[2]:.4f} (se={se[2]:.4f}, t={beta[2]/se[2]:.2f})\n")
            f.write(f"beta_exposure_z_x_shock={beta[3]:.4f} (se={se[3]:.4f}, t={beta[3]/se[3]:.2f})\n")

    print(f"Wrote: {summary_path}")


if __name__ == "__main__":
    main()
