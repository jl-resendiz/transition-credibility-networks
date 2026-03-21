"""Strategy 2 extensions: (1) Event-level Fama-MacBeth slopes, (2) Foreign exposure test.

Primary spec mirrors the paper's Strategy 2:
  - First-mover events only
  - Announcement dates when available
  - Market-adjusted returns (vwretd)
  - Daily CAR[-1,+20] and Monthly CAR[-1,+12]
  - Pre-demean ARs by pre-event mean (consistent with strategy2_spatial_regression.py)

Outputs:
  - Event-level slope summaries (mean, SE, t, N_events)
  - Pooled foreign vs domestic exposure regressions (event-clustered, two-way)
"""
import csv
import math
import os
import hashlib
from collections import defaultdict

from _paths import raw_path, derived_path

EXACT_ONLY = False
PRE_DEMEAN_DAILY = True
PRE_DEMEAN_MONTHLY = True
PRE_DAYS = 250
PRE_MONTHS = 24


# ---- Helpers ----

def load_ff_factors_daily(path):
    if not os.path.exists(path):
        return None
    mktrf, rf, vwretd = {}, {}, {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("This file") or line.startswith("The "):
                continue
            if line.startswith(","):
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
            mktrf_dec = mktrf_val / 100.0
            rf_dec = rf_val / 100.0
            vwretd_dec = mktrf_dec + rf_dec
            date_fmt = f"{date[:4]}-{date[4:6]}-{date[6:]}"
            mktrf[date_fmt] = mktrf_dec
            rf[date_fmt] = rf_dec
            vwretd[date_fmt] = vwretd_dec
    return (mktrf, rf, vwretd) if mktrf else None


def load_ff_factors_monthly(path):
    if not os.path.exists(path):
        return None
    mktrf, rf, vwretd = {}, {}, {}
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("This file") or line.startswith("The "):
                continue
            if line.startswith(","):
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
            mktrf_dec = mktrf_val / 100.0
            rf_dec = rf_val / 100.0
            vwretd_dec = mktrf_dec + rf_dec
            date_fmt = f"{date[:4]}-{date[4:6]}"
            mktrf[date_fmt] = mktrf_dec
            rf[date_fmt] = rf_dec
            vwretd[date_fmt] = vwretd_dec
    return (mktrf, rf, vwretd) if mktrf else None


def compute_daily_car(gvkey, event_year, event_date, daily_ret, market_ret_daily):
    if gvkey not in daily_ret:
        return None
    dates = sorted(daily_ret[gvkey].keys())
    event_idx = None
    if event_date and len(event_date) >= 10:
        ed = event_date[:10]
        for i, d in enumerate(dates):
            if d >= ed:
                event_idx = i
                break
    if event_idx is None:
        yr_s = str(event_year)
        for i, d in enumerate(dates):
            if d.startswith(yr_s):
                event_idx = i
                break
    if event_idx is None:
        return None

    # pre-window AR mean for pre-demean
    pre_mean_ar = 0.0
    if PRE_DEMEAN_DAILY:
        ar_list = []
        for i in range(max(0, event_idx - PRE_DAYS), event_idx):
            d = dates[i]
            if d in daily_ret[gvkey] and d in market_ret_daily:
                ar_list.append(daily_ret[gvkey][d] - market_ret_daily[d])
        if ar_list:
            pre_mean_ar = sum(ar_list) / len(ar_list)

    car = 0.0
    for offset in range(-1, 20 + 1):
        idx = event_idx + offset
        if 0 <= idx < len(dates):
            d = dates[idx]
            if d in daily_ret[gvkey] and d in market_ret_daily:
                ar = daily_ret[gvkey][d] - market_ret_daily[d]
                car += ar - pre_mean_ar if PRE_DEMEAN_DAILY else ar
    return car


def compute_monthly_car(gvkey, event_month, post, monthly_ret, market_ret_monthly):
    if gvkey not in monthly_ret:
        return None
    months = sorted(monthly_ret[gvkey].keys())
    if not event_month or event_month not in months:
        return None
    event_idx = months.index(event_month)

    # pre-window AR mean for pre-demean
    pre_mean_ar = 0.0
    if PRE_DEMEAN_MONTHLY:
        ar_list = []
        for i in range(max(0, event_idx - PRE_MONTHS), event_idx):
            m = months[i]
            if m in monthly_ret[gvkey] and m in market_ret_monthly:
                ar_list.append(monthly_ret[gvkey][m] - market_ret_monthly[m])
        if ar_list:
            pre_mean_ar = sum(ar_list) / len(ar_list)

    car = 0.0
    for offset in range(-1, post + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months):
            m = months[idx]
            if m in monthly_ret[gvkey] and m in market_ret_monthly:
                ar = monthly_ret[gvkey][m] - market_ret_monthly[m]
                car += ar - pre_mean_ar if PRE_DEMEAN_MONTHLY else ar
    return car


def ols_simple(y, x, add_const=True):
    """Simple OLS slope (and intercept) for small event-level regressions."""
    n = len(y)
    if n < 5:
        return None
    if add_const:
        x_mean = sum(x) / n
        y_mean = sum(y) / n
        denom = sum((xi - x_mean) ** 2 for xi in x)
        if denom < 1e-12:
            return None
        b = sum((xi - x_mean) * (yi - y_mean) for xi, yi in zip(x, y)) / denom
        a = y_mean - b * x_mean
        return a, b
    else:
        denom = sum(xi * xi for xi in x)
        if denom < 1e-12:
            return None
        b = sum(xi * yi for xi, yi in zip(x, y)) / denom
        return 0.0, b


def ols_pooled(data, y_key, x_keys):
    """Pooled OLS with event-cluster and two-way cluster SEs (re-using strategy2 formula)."""
    n = len(data)
    k = len(x_keys) + 1
    if n <= k + 1:
        return None
    y = [d[y_key] for d in data]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    X = [[1.0] + [d[x] for x in x_keys] for d in data]
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

    def invert_matrix(mat):
        m = len(mat)
        aug = [row[:] + [1.0 if i == j else 0.0 for j in range(m)] for i, row in enumerate(mat)]
        for col in range(m):
            max_row = max(range(col, m), key=lambda r: abs(aug[r][col]))
            if abs(aug[max_row][col]) < 1e-20:
                return None
            aug[col], aug[max_row] = aug[max_row], aug[col]
            pivot = aug[col][col]
            for j in range(2 * m):
                aug[col][j] /= pivot
            for row in range(m):
                if row != col:
                    factor = aug[row][col]
                    for j in range(2 * m):
                        aug[row][j] -= factor * aug[col][j]
        return [row[m:] for row in aug]

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

    def cluster_cov(clusters):
        k_ = len(X[0])
        S = [[0.0 for _ in range(k_)] for _ in range(k_)]
        for _, idxs in clusters.items():
            xu = [0.0 for _ in range(k_)]
            for i in idxs:
                for a in range(k_):
                    xu[a] += X[i][a] * resid[i]
            for a in range(k_):
                for b in range(k_):
                    S[a][b] += xu[a] * xu[b]
        return S, len(clusters)

    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    # event-cluster
    clusters_event = {}
    for i, d in enumerate(data):
        clusters_event.setdefault(d['event_id'], []).append(i)
    S1, G1 = cluster_cov(clusters_event)
    V1 = mat_mul(mat_mul(inv_XtX, S1), inv_XtX)
    scale = (G1 / (G1 - 1)) * ((n - 1) / (n - k)) if G1 > 1 else 1.0
    for a in range(k):
        for b in range(k):
            V1[a][b] *= scale
    se_event = [math.sqrt(V1[a][a]) if V1[a][a] > 0 else 0.0 for a in range(k)]

    # two-way cluster (event x firm)
    clusters_firm = {}
    clusters_both = {}
    for i, d in enumerate(data):
        clusters_firm.setdefault(d['gvkey'], []).append(i)
        clusters_both.setdefault((d['event_id'], d['gvkey']), []).append(i)
    S2, G2 = cluster_cov(clusters_firm)
    S12, G12 = cluster_cov(clusters_both)
    S = [[S1[a][b] + S2[a][b] - S12[a][b] for b in range(k)] for a in range(k)]
    V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
    G = min(G1, G2)
    if G > 1:
        scale2 = (G / (G - 1)) * ((n - 1) / (n - k))
        for a in range(k):
            for b in range(k):
                V[a][b] *= scale2
    se_tw = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]

    names = ["intercept"] + x_keys
    out = {
        "beta": dict(zip(names, beta)),
        "se_event": dict(zip(names, se_event)),
        "se_tw": dict(zip(names, se_tw)),
        "n": n,
        "clusters": (G1, G2),
        "r2": 1 - sum(r * r for r in resid) / ss_tot if ss_tot > 0 else 0.0,
    }
    return out


# ---- Load data ----

print("Loading returns and factors...")
daily_ret = defaultdict(dict)
with open(derived_path("returns", "daily_returns.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        daily_ret[row["gvkey"]][row["datadate"]] = float(row["ret_daily"])

monthly_ret = defaultdict(dict)
with open(derived_path("returns", "monthly_returns.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        monthly_ret[row["gvkey"]][row["datadate"][:7]] = float(row["ret_monthly"])

ff_daily = load_ff_factors_daily(raw_path("factors", "F-F_Research_Data_Factors_daily.csv"))
ff_monthly = load_ff_factors_monthly(raw_path("factors", "F-F_Research_Data_Factors.csv"))
if not ff_daily or not ff_monthly:
    raise RuntimeError("Missing Fama-French factors for vwretd.")
mktrf_daily, rf_daily, market_ret_daily = ff_daily
mktrf_monthly, rf_monthly, market_ret_monthly = ff_monthly

print("Loading weights and fundamentals...")
W = defaultdict(dict)
with open(derived_path("networks", "weight_matrix_W_geo.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        W[row["gvkey_i"]][row["gvkey_j"]] = float(row["w_ij"])

fundamentals = {}
fundamentals_by_year = defaultdict(dict)
with open(derived_path("fundamentals", "firm_fundamentals.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        gk = row["gvkey"]
        fy = row["fyear"]
        fundamentals_by_year[gk][fy] = row
        if gk not in fundamentals or fy > fundamentals[gk]["fyear"]:
            fundamentals[gk] = row

def get_fundamentals_for_year(gvkey, year):
    rows = fundamentals_by_year.get(gvkey, {})
    if not rows:
        return fundamentals.get(gvkey)
    years = [int(y) for y in rows.keys() if str(y).isdigit()]
    if not years:
        return fundamentals.get(gvkey)
    eligible = [y for y in years if y <= year]
    if eligible:
        return rows[str(max(eligible))]
    return rows[str(max(years))]

def get_sic4(gvkey):
    row = fundamentals.get(gvkey)
    if not row:
        return None
    sic = row.get("sic", "")
    if sic and len(sic) >= 4:
        return sic[:4]
    return None


def is_exact_source(src):
    if not src:
        return False
    s = src.lower()
    if "proxy" in s or "approx" in s or "mid" in s or "month" in s:
        return False
    return True


# Events
all_events = []
with open(derived_path("events", "coal_retirement_events.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        if not row.get("matched_gvkeys"):
            continue
        ann_date = row.get("announcement_date", "").strip()
        ret_date = row.get("event_date", "").strip()
        ann_src = row.get("announcement_source", "").strip()
        if EXACT_ONLY:
            if not ann_date:
                continue
            if ann_src and not is_exact_source(ann_src):
                continue
        effective_date = ann_date if ann_date else ret_date
        event_year = None
        if effective_date and len(effective_date) >= 4 and effective_date[:4].isdigit():
            event_year = int(effective_date[:4])
        else:
            event_year = int(row["ret_year"]) if row.get("ret_year") else None
        all_events.append(
            {
                "plant": row["plant_name"],
                "year": event_year,
                "event_date": effective_date,
                "date_type": "announcement" if ann_date else "retirement",
                "announcement_source": ann_src,
                "gvkeys": row["matched_gvkeys"].split(";"),
                "is_first_mover": row.get("is_first_mover") == "True",
            }
        )
events = [e for e in all_events if e["is_first_mover"]]
print(f"Events (first mover): {len(events)}")


# ---- Build datasets ----

def build_obs(events, monthly_post=12):
    daily_obs = []
    monthly_obs = []
    for event_id, event in enumerate(events):
        event_gvkeys = set(event["gvkeys"])
        year = event["year"]
        event_date = event.get("event_date", "")
        event_month = event_date[:7] if event_date and len(event_date) >= 7 else (f"{year}-07" if year else None)

        # first mover SIC and country
        fm_sic4 = None
        fm_fic = None
        for gk in event_gvkeys:
            fm_sic4 = get_sic4(gk) or fm_sic4
            frow = get_fundamentals_for_year(gk, year) if year else fundamentals.get(gk)
            if frow and frow.get("fic"):
                fm_fic = frow.get("fic")
            if fm_sic4 and fm_fic:
                break

        for fm_gk in event_gvkeys:
            if fm_gk not in W:
                continue
            neighbors = W[fm_gk]

            neighbor_gks = set(neighbors.keys()) - event_gvkeys
            non_connected = [gk for gk in fundamentals if gk not in event_gvkeys and gk not in neighbors]
            import random as _rng
            stable_seed = int(hashlib.md5(str(fm_gk).encode("utf-8")).hexdigest()[:8], 16)
            _rng.seed(stable_seed)
            n_ctrl = min(len(non_connected), max(5 * len(neighbor_gks), 20))
            ctrl_sample = _rng.sample(non_connected, n_ctrl) if len(non_connected) > n_ctrl else non_connected
            candidate_firms = list(neighbor_gks) + ctrl_sample

            for gk in candidate_firms:
                w_ij = neighbors.get(gk, 0.0)
                frow = get_fundamentals_for_year(gk, year) if year else fundamentals.get(gk)
                fic = frow.get("fic") if frow else None
                same_sector = 1.0 if (fm_sic4 and get_sic4(gk) and fm_sic4 == get_sic4(gk)) else 0.0
                foreign = 1.0 if (fm_fic and fic and fm_fic != fic) else 0.0
                w_foreign = w_ij if foreign else 0.0
                w_domestic = w_ij if not foreign else 0.0

                # Daily
                car_d = compute_daily_car(gk, year, event_date, daily_ret, market_ret_daily) if year else None
                if car_d is not None:
                    daily_obs.append({
                        "car": car_d,
                        "w_ij": w_ij,
                        "w_foreign": w_foreign,
                        "w_domestic": w_domestic,
                        "same_sector": same_sector,
                        "event_id": event_id,
                        "gvkey": gk,
                    })

                # Monthly
                if event_month:
                    car_m = compute_monthly_car(gk, event_month, monthly_post, monthly_ret, market_ret_monthly)
                    if car_m is not None:
                        monthly_obs.append({
                            "car": car_m,
                            "w_ij": w_ij,
                            "w_foreign": w_foreign,
                            "w_domestic": w_domestic,
                            "same_sector": same_sector,
                            "event_id": event_id,
                            "gvkey": gk,
                        })
    return daily_obs, monthly_obs


daily_obs, monthly_obs = build_obs(events, monthly_post=12)
print(f"Daily obs: {len(daily_obs)}; Monthly obs (+12): {len(monthly_obs)}")


# ---- (1) Fama-MacBeth by event ----

def fmb_by_event(obs, label):
    by_event = defaultdict(list)
    for o in obs:
        by_event[o["event_id"]].append(o)

    slopes = []
    for eid, rows in by_event.items():
        # need variation in w_ij
        wvals = [r["w_ij"] for r in rows]
        if len(set(wvals)) < 3:
            continue
        y = [r["car"] for r in rows]
        x = wvals
        res = ols_simple(y, x, add_const=True)
        if res is None:
            continue
        _, b = res
        slopes.append(b)

    if not slopes:
        print(f"{label}: no usable events")
        return None
    mean_b = sum(slopes) / len(slopes)
    var_b = sum((b - mean_b) ** 2 for b in slopes) / (len(slopes) - 1) if len(slopes) > 1 else 0.0
    se = math.sqrt(var_b / len(slopes)) if len(slopes) > 1 else 0.0
    t = mean_b / se if se > 1e-12 else 0.0
    print(f"{label}: N_events={len(slopes)}, mean b={mean_b:.4f}, SE={se:.4f}, t={t:.2f}")
    return {"n_events": len(slopes), "mean_b": mean_b, "se": se, "t": t}


print("\nFAMA-MACBETH (event-level slopes)")
fmb_daily = fmb_by_event(daily_obs, "Daily [-1,+20]")
fmb_monthly = fmb_by_event(monthly_obs, "Monthly [-1,+12]")


# ---- (2) Foreign exposure test ----

def run_foreign(obs, label):
    # Regress CAR on foreign + domestic exposure + same_sector
    x_keys = ["w_foreign", "w_domestic", "same_sector"]
    res = ols_pooled(obs, "car", x_keys)
    if res is None:
        print(f"{label}: regression failed")
        return None
    b = res["beta"]
    se_e = res["se_event"]
    se_tw = res["se_tw"]
    print(f"\n{label} (pooled)")
    for k in ["w_foreign", "w_domestic", "same_sector"]:
        print(f"  {k}: {b[k]:.4f} (SE_event={se_e[k]:.4f}, SE_tw={se_tw[k]:.4f})")
    print(f"  N={res['n']}, clusters(event,firm)={res['clusters']}, R2={res['r2']:.4f}")
    return res


print("\nFOREIGN EXPOSURE TEST (pooled)")
foreign_daily = run_foreign(daily_obs, "Daily [-1,+20]")
foreign_monthly = run_foreign(monthly_obs, "Monthly [-1,+12]")


# ---- Save summary ----

base_dir = os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(derived_path("summaries", "x")))))
out_path = os.path.join(base_dir, "results", "summaries", "strategy2_event_fmb_foreign_summary.md")
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, "w", encoding="utf-8") as f:
    f.write("# Strategy 2: Event-level Fama-MacBeth and Foreign Exposure\n\n")
    f.write(f"Events (first-mover): {len(events)}\n")
    f.write(f"Daily obs: {len(daily_obs)}; Monthly obs (+12): {len(monthly_obs)}\n\n")
    if fmb_daily:
        f.write(f"## Fama-MacBeth (Daily [-1,+20])\n")
        f.write(f"- N_events: {fmb_daily['n_events']}\n")
        f.write(f"- mean b: {fmb_daily['mean_b']:.4f}\n")
        f.write(f"- SE: {fmb_daily['se']:.4f}\n")
        f.write(f"- t: {fmb_daily['t']:.2f}\n\n")
    if fmb_monthly:
        f.write(f"## Fama-MacBeth (Monthly [-1,+12])\n")
        f.write(f"- N_events: {fmb_monthly['n_events']}\n")
        f.write(f"- mean b: {fmb_monthly['mean_b']:.4f}\n")
        f.write(f"- SE: {fmb_monthly['se']:.4f}\n")
        f.write(f"- t: {fmb_monthly['t']:.2f}\n\n")
    if foreign_daily:
        f.write("## Foreign exposure (Daily)\n")
        for k in ["w_foreign", "w_domestic", "same_sector"]:
            f.write(f"- {k}: {foreign_daily['beta'][k]:.4f} (SE_event={foreign_daily['se_event'][k]:.4f}, SE_tw={foreign_daily['se_tw'][k]:.4f})\n")
        f.write(f"- N={foreign_daily['n']}, R2={foreign_daily['r2']:.4f}\n\n")
    if foreign_monthly:
        f.write("## Foreign exposure (Monthly)\n")
        for k in ["w_foreign", "w_domestic", "same_sector"]:
            f.write(f"- {k}: {foreign_monthly['beta'][k]:.4f} (SE_event={foreign_monthly['se_event'][k]:.4f}, SE_tw={foreign_monthly['se_tw'][k]:.4f})\n")
        f.write(f"- N={foreign_monthly['n']}, R2={foreign_monthly['r2']:.4f}\n\n")

print(f"\nSaved summary -> {out_path}")
