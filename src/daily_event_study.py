"""Cross-sectional regression of daily CARs on the network weights.

For each daily window in {[-1,+1], [0,+5], [0,+10], [-1,+10], [0,+20]}:
  CAR_ie = alpha + gamma_geo w^geo_i + gamma_fuel w^fuel_i + gamma_reg w^reg_i
           + gamma_s SameSector_i + eps_ie
estimated by:
  (a) Pooled OLS with event-clustered SEs
  (b) Fama-MacBeth with Newey-West HAC (lag 4)
  (c) Two-way clustering (event × firm), CGM 2011

Reads `daily_car_panel.csv` from compute_daily_ar_panel.py.
Output: results/metrics/daily_event_study.md
"""
import csv
import math
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import derived_path, results_path

CAR_FF3 = os.path.join(derived_path('returns'), 'daily_car_panel.csv')
CAR_MKT = os.path.join(derived_path('returns'), 'daily_car_marketadj_panel.csv')
# Default: FF3 (legacy default). Override via CLI: `python daily_event_study.py --mkt`
import sys as _sys
if '--mkt' in _sys.argv:
    CAR_PATH = CAR_MKT
    OUT_PATH = os.path.join(results_path('metrics'), 'daily_event_study_marketadj.md')
else:
    CAR_PATH = CAR_FF3
    OUT_PATH = os.path.join(results_path('metrics'), 'daily_event_study.md')

WINDOWS = [
    ('car_m1_p1',  '[-1, +1]'),
    ('car_0_p5',   '[0, +5]'),
    ('car_0_p10',  '[0, +10]'),
    ('car_m1_p10', '[-1, +10]'),
    ('car_0_p20',  '[0, +20]'),
]


# ── Linear-algebra helpers ────────────────────────────────────────────

def invert_matrix(mat):
    n = len(mat)
    aug = [row[:] + [1.0 if i == j else 0.0 for j in range(n)] for i, row in enumerate(mat)]
    for col in range(n):
        max_r = max(range(col, n), key=lambda r: abs(aug[r][col]))
        if abs(aug[max_r][col]) < 1e-20:
            return None
        aug[col], aug[max_r] = aug[max_r], aug[col]
        piv = aug[col][col]
        for j in range(2 * n):
            aug[col][j] /= piv
        for r in range(n):
            if r != col:
                f = aug[r][col]
                for j in range(2 * n):
                    aug[r][j] -= f * aug[col][j]
    return [row[n:] for row in aug]


def mat_mul(a, b):
    rows_a, mid, cols_b = len(a), len(b), len(b[0])
    out = [[0.0] * cols_b for _ in range(rows_a)]
    for i in range(rows_a):
        for k in range(mid):
            v = a[i][k]
            if v == 0: continue
            for j in range(cols_b):
                out[i][j] += v * b[k][j]
    return out


def cluster_meat(X, resid, cluster_keys, k):
    cmap = defaultdict(list)
    for i, c in enumerate(cluster_keys):
        cmap[c].append(i)
    S = [[0.0] * k for _ in range(k)]
    for idxs in cmap.values():
        xu = [0.0] * k
        for i in idxs:
            ri = resid[i]
            for a in range(k):
                xu[a] += X[i][a] * ri
        for a in range(k):
            for b in range(a, k):
                v = xu[a] * xu[b]
                S[a][b] += v
                if a != b:
                    S[b][a] += v
    return S, len(cmap)


def ols_pooled(rows, y_key, x_keys, cluster_key='event_id'):
    n = len(rows)
    k = len(x_keys) + 1
    if n <= k + 1: return None
    y = [r[y_key] for r in rows]
    X = [[1.0] + [r[xk] for xk in x_keys] for r in rows]
    XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
    inv = invert_matrix(XtX)
    if inv is None: return None
    beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
    yhat = [sum(X[i][a]*beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - yhat[i] for i in range(n)]
    ss_tot = sum((yi - sum(y)/n)**2 for yi in y)
    ss_res = sum(r**2 for r in resid)
    r2 = 1 - ss_res/ss_tot if ss_tot > 0 else 0
    # event-clustered V
    cluster_ids = [r[cluster_key] for r in rows]
    S, G = cluster_meat(X, resid, cluster_ids, k)
    V = mat_mul(mat_mul(inv, S), inv)
    if G > 1:
        scale = (G/(G-1)) * ((n-1)/(n-k))
        for a in range(k):
            for b in range(k):
                V[a][b] *= scale
    se = [math.sqrt(V[a][a]) if V[a][a] > 0 else 0 for a in range(k)]
    t = [beta[a]/se[a] if se[a] > 1e-15 else 0 for a in range(k)]
    # Two-way (event × firm) clustering
    firm_ids = [r['gvkey'] for r in rows]
    S_e, G_e = cluster_meat(X, resid, cluster_ids, k)
    S_f, G_f = cluster_meat(X, resid, firm_ids, k)
    S_ef, G_ef = cluster_meat(X, resid, list(zip(cluster_ids, firm_ids)), k)
    def vmat(S, G):
        V_ = mat_mul(mat_mul(inv, S), inv)
        if G > 1:
            sc = (G/(G-1)) * ((n-1)/(n-k))
            for a in range(k):
                for b in range(k):
                    V_[a][b] *= sc
        return V_
    Ve = vmat(S_e, G_e); Vf = vmat(S_f, G_f); Vef = vmat(S_ef, G_ef)
    Vt = [[Ve[a][b] + Vf[a][b] - Vef[a][b] for b in range(k)] for a in range(k)]
    for a in range(k):
        if Vt[a][a] <= 0:
            Vt[a][a] = max(Ve[a][a], Vf[a][a])
    se_tw = [math.sqrt(Vt[a][a]) if Vt[a][a] > 0 else 0 for a in range(k)]
    t_tw = [beta[a]/se_tw[a] if se_tw[a] > 1e-15 else 0 for a in range(k)]
    return {
        'beta': beta, 'se': se, 't': t, 'se_tw': se_tw, 't_tw': t_tw,
        'n': n, 'r2': r2, 'G_event': G_e, 'G_firm': G_f,
        'V_event': Ve, 'V_tw': Vt, 'inv_XtX': inv, 'X': X, 'resid': resid,
    }


def fm_nw(rows, y_key, x_keys, lag=4, min_firms=20):
    """Fama-MacBeth with Newey-West HAC."""
    by_event = defaultdict(list)
    for r in rows:
        by_event[r['event_id']].append(r)
    event_betas = []  # list of beta vectors per event
    event_ids = []
    for eid in sorted(by_event.keys()):
        evrows = by_event[eid]
        if len(evrows) < min_firms:
            continue
        n = len(evrows)
        k = len(x_keys) + 1
        y = [r[y_key] for r in evrows]
        X = [[1.0] + [r[xk] for xk in x_keys] for r in evrows]
        XtX = [[sum(X[i][a]*X[i][b] for i in range(n)) for b in range(k)] for a in range(k)]
        Xty = [sum(X[i][a]*y[i] for i in range(n)) for a in range(k)]
        inv = invert_matrix(XtX)
        if inv is None:
            continue
        beta = [sum(inv[a][b]*Xty[b] for b in range(k)) for a in range(k)]
        event_betas.append(beta)
        event_ids.append(eid)

    T = len(event_betas)
    if T < 2:
        return None
    k = len(event_betas[0])
    means = [sum(b[a] for b in event_betas) / T for a in range(k)]
    # NW HAC SE
    nw_se = []
    nw_t = []
    for a in range(k):
        x = [event_betas[t][a] - means[a] for t in range(T)]
        # gamma_0
        S = sum(xt*xt for xt in x) / T
        for L in range(1, lag + 1):
            w = 1 - L/(lag + 1)
            cov = sum(x[t]*x[t-L] for t in range(L, T)) / T
            S += 2 * w * cov
        var = S / T
        se = math.sqrt(max(var, 0))
        nw_se.append(se)
        nw_t.append(means[a] / se if se > 1e-15 else 0)
    # Difference test (geo - fuel) — assumes specific x_key ordering
    if 'w_geo' in x_keys and 'w_fuel' in x_keys:
        idx_g = x_keys.index('w_geo') + 1
        idx_f = x_keys.index('w_fuel') + 1
        diffs = [b[idx_g] - b[idx_f] for b in event_betas]
        diff_mean = sum(diffs) / T
        diff_dem = [d - diff_mean for d in diffs]
        S = sum(d*d for d in diff_dem) / T
        for L in range(1, lag + 1):
            w = 1 - L/(lag + 1)
            cov = sum(diff_dem[t]*diff_dem[t-L] for t in range(L, T)) / T
            S += 2 * w * cov
        var = S / T
        diff_se = math.sqrt(max(var, 0))
        diff_t = diff_mean / diff_se if diff_se > 1e-15 else 0
    else:
        diff_mean = diff_se = diff_t = None
    return {
        'mean': means, 'nw_se': nw_se, 'nw_t': nw_t, 'T': T,
        'diff_mean': diff_mean, 'diff_se': diff_se, 'diff_t': diff_t,
    }


def stars(t):
    p = 2 * (1 - 0.5 * (1 + math.erf(abs(t) / math.sqrt(2))))
    if p < 0.01: return '***'
    if p < 0.05: return '**'
    if p < 0.10: return '*'
    return ''


# ─── Main ──────────────────────────────────────────────────────────────

def main():
    print(f'Loading {CAR_PATH}...')
    rows_all = []
    with open(CAR_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                r = {
                    'gvkey': row['gvkey'],
                    'event_id': int(row['event_id']),
                    'w_geo': float(row['w_geo']),
                    'w_fuel': float(row['w_fuel']),
                    'w_reg': float(row['w_reg']),
                    'same_sector': float(row['same_sector']),
                }
                # parse CAR fields (skip None)
                for car_key, _ in WINDOWS:
                    v = row[car_key].strip()
                    r[car_key] = float(v) if v not in ('', 'None') else None
                rows_all.append(r)
            except (ValueError, KeyError):
                continue
    print(f'  {len(rows_all):,} firm-event rows.')

    spec_vars = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']
    names = ['intercept'] + spec_vars

    # Run per-window cross-sectional regressions
    results_pool = {}
    results_fm = {}

    for car_key, label in WINDOWS:
        rows = [r for r in rows_all if r[car_key] is not None]
        print(f'\nWindow {label} (key={car_key}): N = {len(rows):,}')

        res_p = ols_pooled(rows, car_key, spec_vars)
        results_pool[car_key] = res_p
        if res_p:
            for i, nm in enumerate(names):
                print(f'  {nm:13s}: beta = {res_p["beta"][i]:+.6f}, '
                      f'SE_event = {res_p["se"][i]:.6f}, t = {res_p["t"][i]:+.3f}, '
                      f'SE_tw = {res_p["se_tw"][i]:.6f}, t_tw = {res_p["t_tw"][i]:+.3f}')
            print(f'  R^2 = {res_p["r2"]:.4f}, event clusters = {res_p["G_event"]}, '
                  f'firm clusters = {res_p["G_firm"]}')

        res_fm = fm_nw(rows, car_key, spec_vars, lag=4, min_firms=20)
        results_fm[car_key] = res_fm
        if res_fm:
            print(f'  FM (T = {res_fm["T"]}, NW lag 4):')
            for i, nm in enumerate(names):
                print(f'    {nm:13s}: mean = {res_fm["mean"][i]:+.6f}, '
                      f'SE_NW = {res_fm["nw_se"][i]:.6f}, t = {res_fm["nw_t"][i]:+.3f}')
            if res_fm['diff_t'] is not None:
                print(f'    Diff (geo-fuel): {res_fm["diff_mean"]:+.6f}, '
                      f'SE = {res_fm["diff_se"]:.6f}, t = {res_fm["diff_t"]:+.3f}')

    # ── Write output markdown ─────────────────────────────────────────
    lines = [
        '# Daily Event Study (FF3 abnormal returns) around Announcement Dates',
        '',
        'Cross-sectional regression of daily-aggregated CARs on the same channel',
        'weights as the headline monthly regression. Estimation window for',
        'firm-by-firm FF3 betas: [-252, -22] trading days before announcement_date.',
        '',
        'Spec: CAR_ie = a + gamma_geo w^geo_i + gamma_fuel w^fuel_i',
        '              + gamma_reg w^reg_i + gamma_s SameSector_i + eps_ie.',
        '',
        '## Pooled OLS (event-clustered SEs) and two-way (event × firm) clustering',
        '',
    ]

    # Build comparison table across windows for w_fuel
    fuel_idx = names.index('w_fuel')
    geo_idx = names.index('w_geo')
    reg_idx = names.index('w_reg')

    lines.append('| Window | N | gamma_fuel | t (event-cl) | t (two-way) | gamma_geo | t (event-cl) | t (two-way) | R^2 |')
    lines.append('|---|---:|---:|---:|---:|---:|---:|---:|---:|')
    for car_key, label in WINDOWS:
        r = results_pool.get(car_key)
        if r is None:
            continue
        lines.append(
            f'| {label} | {r["n"]:,} '
            f'| {r["beta"][fuel_idx]:+.4f}{stars(r["t"][fuel_idx])} '
            f'| {r["t"][fuel_idx]:+.3f} | {r["t_tw"][fuel_idx]:+.3f} '
            f'| {r["beta"][geo_idx]:+.4f}{stars(r["t"][geo_idx])} '
            f'| {r["t"][geo_idx]:+.3f} | {r["t_tw"][geo_idx]:+.3f} '
            f'| {r["r2"]:.4f} |'
        )

    lines += ['', '## Fama-MacBeth (Newey-West, lag=4)', '']
    lines.append('| Window | T (events) | gamma_fuel | NW t | gamma_geo | NW t | (geo - fuel) | NW t |')
    lines.append('|---|---:|---:|---:|---:|---:|---:|---:|')
    for car_key, label in WINDOWS:
        r = results_fm.get(car_key)
        if r is None:
            continue
        diff_str = f'{r["diff_mean"]:+.4f}' if r["diff_mean"] is not None else '—'
        diff_t = f'{r["diff_t"]:+.3f}' if r["diff_t"] is not None else '—'
        lines.append(
            f'| {label} | {r["T"]} '
            f'| {r["mean"][fuel_idx]:+.4f}{stars(r["nw_t"][fuel_idx])} '
            f'| {r["nw_t"][fuel_idx]:+.3f} '
            f'| {r["mean"][geo_idx]:+.4f}{stars(r["nw_t"][geo_idx])} '
            f'| {r["nw_t"][geo_idx]:+.3f} '
            f'| {diff_str} | {diff_t} |'
        )

    lines += [
        '',
        '## Interpretation',
        '',
        'The daily event-study uses precise YYYY-MM-DD announcement dates for',
        'all 179 first-mover events. Firm-by-firm FF3 abnormal returns are',
        'computed on a [-252, -22] daily estimation window, then aggregated to',
        'the cumulative windows above. The cross-sectional regression mirrors',
        'the headline monthly specification.',
        '',
        'A negative and statistically significant gamma_fuel in short windows',
        '([-1,+1] or [0,+5]) confirms that the fuel-mix channel transmits at',
        'announcement-day frequency, addressing the referee concern that the',
        '4-month monthly window is non-standard. A null daily effect with',
        'large monthly effect would, instead, support a gradual-diffusion',
        'mechanism (Hong-Stein 1999; Cohen-Frazzini 2008) — which is also',
        'a defensible interpretation consistent with the paper.',
        '',
    ]

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, 'w', encoding='utf-8') as f:
        f.write('\n'.join(lines))
    print(f'\nWrote {OUT_PATH}')


if __name__ == '__main__':
    main()
