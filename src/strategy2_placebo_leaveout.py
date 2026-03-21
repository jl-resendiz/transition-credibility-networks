"""Strategy 2 robustness: placebo timing and leave-one-out checks.

1) Placebo timing: shift event dates by +/- K months and re-run the
   monthly exposure regression ([-1,+12], vwretd, first-movers).
2) Leave-one-out: drop each event (and each country) and re-estimate
   exposure coefficient distribution.
"""
import csv
import os
import hashlib
from collections import defaultdict

from _paths import raw_path, derived_path
from datetime import datetime

import math

# Settings (match primary spec)
PRE_DEMEAN_MONTHLY = True
PRE_MONTHS = 24
POST = 12
SHIFT_MONTHS = [-12, -6, 6, 12]
EXACT_ONLY = False


# ---------- Utilities ----------

def parse_month(ym):
    return datetime.strptime(ym, "%Y-%m")


def add_months(dt, k):
    y = dt.year + (dt.month - 1 + k) // 12
    m = (dt.month - 1 + k) % 12 + 1
    return datetime(y, m, 1)


def load_ff_factors_monthly(path):
    if not os.path.exists(path):
        return None
    mktrf = {}
    rf = {}
    vwretd = {}
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


def _cluster_cov(X, resid, clusters):
    k = len(X[0])
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


def ols_clustered(data, y_key, x_keys):
    n = len(data)
    k = len(x_keys) + 1
    if n <= k + 1:
        return None
    y = [d[y_key] for d in data]
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    if ss_tot < 1e-12:
        return None
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
    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None

    clusters_event = {}
    clusters_firm = {}
    clusters_both = {}
    for i, d in enumerate(data):
        clusters_event.setdefault(d["event_id"], []).append(i)
        clusters_firm.setdefault(d["gvkey"], []).append(i)
        clusters_both.setdefault((d["event_id"], d["gvkey"]), []).append(i)
    S1, G1 = _cluster_cov(X, resid, clusters_event)
    S2, G2 = _cluster_cov(X, resid, clusters_firm)
    S12, G12 = _cluster_cov(X, resid, clusters_both)
    S = [[S1[a][b] + S2[a][b] - S12[a][b] for b in range(k)] for a in range(k)]
    V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
    G = min(G1, G2)
    if G > 1:
        scale = (G / (G - 1)) * ((n - 1) / (n - k))
        for a in range(k):
            for b in range(k):
                V[a][b] *= scale
    se_tw = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0.0 for a in range(k)]

    # event-cluster only
    V1 = mat_mul(mat_mul(inv_XtX, S1), inv_XtX)
    if G1 > 1:
        scale1 = (G1 / (G1 - 1)) * ((n - 1) / (n - k))
        for a in range(k):
            for b in range(k):
                V1[a][b] *= scale1
    se_event = [math.sqrt(V1[a][a]) if V1[a][a] > 0 else 0.0 for a in range(k)]

    names = ["intercept"] + x_keys
    return {
        "beta": dict(zip(names, beta)),
        "se_event": dict(zip(names, se_event)),
        "se_tw": dict(zip(names, se_tw)),
        "n": n,
        "clusters": (G1, G2),
        "r2": 1 - sum(r ** 2 for r in resid) / ss_tot,
    }


# ---------- Load data ----------

monthly_ret = defaultdict(dict)
with open(derived_path("returns", "monthly_returns.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        monthly_ret[row["gvkey"]][row["datadate"][:7]] = float(row["ret_monthly"])

ff_monthly = load_ff_factors_monthly(raw_path("factors", "F-F_Research_Data_Factors.csv"))
if not ff_monthly:
    raise RuntimeError("Missing Fama-French monthly factors.")
mktrf_monthly, rf_monthly, market_ret_monthly = ff_monthly

W = defaultdict(dict)
with open(derived_path("networks", "weight_matrix_W_geo.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        W[row["gvkey_i"]][row["gvkey_j"]] = float(row["w_ij"])

fundamentals = {}
with open(derived_path("fundamentals", "firm_fundamentals.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        gk = row["gvkey"]
        if gk not in fundamentals or row["fyear"] > fundamentals[gk]["fyear"]:
            fundamentals[gk] = row

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


# ---------- Build dataset ----------

def build_monthly_obs(events, shift_months=0):
    obs = []
    for event_id, event in enumerate(events):
        event_gvkeys = set(event["gvkeys"])
        year = event["year"]
        event_date = event.get("event_date", "")
        if event_date and len(event_date) >= 7:
            event_month = event_date[:7]
        else:
            event_month = f"{year}-07" if year else None
        if event_month and shift_months != 0:
            try:
                dt = parse_month(event_month)
                event_month = add_months(dt, shift_months).strftime("%Y-%m")
            except Exception:
                event_month = None

        fm_sic4 = None
        fm_fic = None
        for gk in event_gvkeys:
            fm_sic4 = get_sic4(gk) or fm_sic4
            frow = fundamentals.get(gk)
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
                frow = fundamentals.get(gk)
                fic = frow.get("fic") if frow else None
                same_sector = 1.0 if (fm_sic4 and get_sic4(gk) and fm_sic4 == get_sic4(gk)) else 0.0

                car_m = compute_monthly_car(gk, event_month, POST, monthly_ret, market_ret_monthly) if event_month else None
                if car_m is not None:
                    obs.append({
                        "car": car_m,
                        "w_ij": w_ij,
                        "same_sector": same_sector,
                        "event_id": event_id,
                        "gvkey": gk,
                    })
    return obs


# ---------- Placebo timing ----------

def run_placebo():
    base_obs = build_monthly_obs(events, shift_months=0)
    res_base = ols_clustered(base_obs, "car", ["w_ij", "same_sector"])
    print("\nBASE (true dates):")
    if res_base:
        print(f"  w_ij = {res_base['beta']['w_ij']:.4f} (SE_event={res_base['se_event']['w_ij']:.4f}, SE_tw={res_base['se_tw']['w_ij']:.4f})")
    results = []
    for k in SHIFT_MONTHS:
        obs = build_monthly_obs(events, shift_months=k)
        res = ols_clustered(obs, "car", ["w_ij", "same_sector"])
        if res:
            results.append((k, res))
            print(f"SHIFT {k:+} months: w_ij={res['beta']['w_ij']:.4f} (SE_event={res['se_event']['w_ij']:.4f}, SE_tw={res['se_tw']['w_ij']:.4f}) N={res['n']}")
    return res_base, results


# ---------- Leave-one-out ----------

def run_leave_one_out():
    # Leave-one-event-out
    betas_evt = []
    for i in range(len(events)):
        ev = [e for j, e in enumerate(events) if j != i]
        obs = build_monthly_obs(ev, shift_months=0)
        res = ols_clustered(obs, "car", ["w_ij", "same_sector"])
        if res:
            betas_evt.append(res["beta"]["w_ij"])
    # Leave-one-country-out (by first-mover country)
    fm_country = {}
    for idx, e in enumerate(events):
        c = None
        for gk in e["gvkeys"]:
            frow = fundamentals.get(gk)
            if frow and frow.get("fic"):
                c = frow["fic"]
                break
        fm_country[idx] = c

    unique_c = sorted({c for c in fm_country.values() if c})
    betas_cty = []
    for c in unique_c:
        ev = [e for idx, e in enumerate(events) if fm_country.get(idx) != c]
        obs = build_monthly_obs(ev, shift_months=0)
        res = ols_clustered(obs, "car", ["w_ij", "same_sector"])
        if res:
            betas_cty.append(res["beta"]["w_ij"])

    def summarize(arr):
        if not arr:
            return None
        mean = sum(arr) / len(arr)
        p5 = sorted(arr)[int(0.05 * (len(arr)-1))]
        p95 = sorted(arr)[int(0.95 * (len(arr)-1))]
        return mean, p5, p95, len(arr)

    s_evt = summarize(betas_evt)
    s_cty = summarize(betas_cty)
    return s_evt, s_cty


if __name__ == "__main__":
    print("Running placebo timing...")
    res_base, res_shift = run_placebo()
    print("\nRunning leave-one-out...")
    s_evt, s_cty = run_leave_one_out()
    if s_evt:
        print(f"Leave-one-event-out: mean={s_evt[0]:.4f}, p5={s_evt[1]:.4f}, p95={s_evt[2]:.4f}, N={s_evt[3]}")
    if s_cty:
        print(f"Leave-one-country-out: mean={s_cty[0]:.4f}, p5={s_cty[1]:.4f}, p95={s_cty[2]:.4f}, N={s_cty[3]}")

    # Save summary
    out_path = os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(derived_path("summaries", "x"))))),
                             "results", "summaries", "strategy2_placebo_leaveout_summary.md")
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w", encoding="utf-8") as f:
        f.write("# Strategy 2: Placebo Timing and Leave-One-Out\n\n")
        if res_base:
            f.write(f"Base w_ij: {res_base['beta']['w_ij']:.4f} (SE_event={res_base['se_event']['w_ij']:.4f}, SE_tw={res_base['se_tw']['w_ij']:.4f})\n\n")
        f.write("## Placebo timing shifts\n")
        for k, res in res_shift:
            f.write(f"- Shift {k:+} months: w_ij={res['beta']['w_ij']:.4f} (SE_event={res['se_event']['w_ij']:.4f}, SE_tw={res['se_tw']['w_ij']:.4f}), N={res['n']}\n")
        f.write("\n## Leave-one-out\n")
        if s_evt:
            f.write(f"- Leave-one-event-out: mean={s_evt[0]:.4f}, p5={s_evt[1]:.4f}, p95={s_evt[2]:.4f}, N={s_evt[3]}\n")
        if s_cty:
            f.write(f"- Leave-one-country-out: mean={s_cty[0]:.4f}, p5={s_cty[1]:.4f}, p95={s_cty[2]:.4f}, N={s_cty[3]}\n")

    print(f"\nSaved summary -> {out_path}")
