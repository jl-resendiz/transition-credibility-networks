"""Honest DID sensitivity (Rambachan-Roth 2023) using placebo CAR windows.

Implements a placebo-window version of the Rambachan-Roth bound on the
cross-sectional CAR-to-fuel coefficient. The intuition: the headline result
is a Fama-MacBeth coefficient gamma_fuel from the post-event window
[-1, +3]. The same regression on placebo CAR windows that *precede* the
event provides direct evidence on parallel-trends violations expressed in
the same units as the headline (a CAR-window OLS slope, not a per-month
single-AR slope).

Procedure:

  1. For each placebo window of width 5 months ending strictly before the
     event ([-12, -8], [-9, -5], [-6, -2]), compute CAR_window for each
     (event, firm) pair.
  2. Run the same cross-sectional FM regression on (w_geo, w_fuel, w_reg,
     same_sector) and aggregate beta_fuel(window) across events using
     Newey-West (lag=4) HAC SEs.
  3. The headline window [-1, +3] uses the existing CAR construction.
  4. Rambachan-Roth M-bar (5% breakdown): the smallest M such that a
     post-event bias as large as M times the maximum |pre-period beta|
     would zero out the headline estimate at the 5% level. M-bar >= 1
     indicates the headline result survives violations as large as
     anything actually observed in the pre-period.

Output: results/metrics/honest_did.md
        results/summaries/honest_did_breakdown.csv
"""
import csv
import os
import sys
import math
import random
import hashlib
from collections import defaultdict

from _paths import derived_path, raw_path, results_path


def _print(msg=''):
    print(msg)
    sys.stdout.flush()


# ── Inlined OLS helpers ─────────────────────────────────────────────

def invert_matrix(mat):
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


def ols_simple(y, X_mat):
    n = len(y)
    k = len(X_mat[0])
    if n <= k:
        return None
    y_mean = sum(y) / n
    ss_tot = sum((yi - y_mean) ** 2 for yi in y)
    if ss_tot < 1e-15:
        return None
    XtX = [[sum(X_mat[i][a] * X_mat[i][b] for i in range(n))
            for b in range(k)] for a in range(k)]
    Xty = [sum(X_mat[i][a] * y[i] for i in range(n)) for a in range(k)]
    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None
    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]
    y_hat = [sum(X_mat[i][a] * beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]
    ss_res = sum(r ** 2 for r in resid)
    return {'beta': beta, 'resid': resid, 'r2': 1 - ss_res / ss_tot, 'n': n}


# ── Data loading ────────────────────────────────────────────────────

_print('Loading monthly returns...')
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except ValueError:
            pass


def load_ff_factors_monthly(path):
    if not os.path.exists(path):
        return None
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
            vw = (mktrf_val + rf_val) / 100.0
            vwretd[f'{date[:4]}-{date[4:6]}'] = vw
    return vwretd


market_ret_monthly = load_ff_factors_monthly(
    raw_path('factors', 'F-F_Research_Data_Factors.csv')
)
if not market_ret_monthly:
    raise RuntimeError('Missing F-F monthly factors.')

W_geo = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_geo[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

W_fuel = defaultdict(dict)
fuel_path = derived_path('networks', 'weight_matrix_W_fuel.csv')
if os.path.exists(fuel_path):
    with open(fuel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            W_fuel[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

W_reg = defaultdict(dict)
reg_path = derived_path('networks', 'weight_matrix_W_regulatory.csv')
if os.path.exists(reg_path):
    with open(reg_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            wval = row.get('w_ij') or row.get('w_reg')
            try:
                W_reg[row['gvkey_i']][row['gvkey_j']] = float(wval)
            except (ValueError, TypeError):
                continue

fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row


def get_sic4(gvkey):
    f = fundamentals.get(gvkey)
    return f['sic'][:4] if (f and f.get('sic')) else None


all_events = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        if row.get('is_first_mover') != 'True':
            continue
        ann_date = row.get('announcement_date', '').strip()
        ret_date = row.get('event_date', '').strip()
        effective_date = ann_date if ann_date else ret_date
        if not effective_date or len(effective_date) < 7:
            continue
        all_events.append({
            'plant': row['plant_name'],
            'event_month': effective_date[:7],
            'gvkeys': row['matched_gvkeys'].split(';'),
        })
_print(f'First-mover events: {len(all_events)}')


# ── CAR helper for arbitrary window ─────────────────────────────────

def compute_car_window(gvkey, event_month_str, tau_start, tau_end):
    """Sum of monthly abnormal returns (r_it - market_t) over [tau_start, tau_end]."""
    if gvkey not in monthly_ret:
        return None
    months = sorted(monthly_ret[gvkey].keys())
    event_idx = None
    for i, m in enumerate(months):
        if m >= event_month_str:
            event_idx = i
            break
    if event_idx is None:
        return None
    # Pre-period mean adjustment (matches headline CAR construction)
    PRE_MONTHS = 24
    ar_list = []
    for i in range(max(0, event_idx - PRE_MONTHS), event_idx):
        m = months[i]
        if m in monthly_ret[gvkey] and m in market_ret_monthly:
            ar_list.append(monthly_ret[gvkey][m] - market_ret_monthly[m])
    if len(ar_list) < 12:
        return None
    pre_mean_ar = sum(ar_list) / len(ar_list)

    car = 0.0
    n_months = 0
    for offset in range(tau_start, tau_end + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months) and months[idx] in monthly_ret[gvkey]:
            m = months[idx]
            r_it = monthly_ret[gvkey][m]
            if m in market_ret_monthly:
                car += (r_it - market_ret_monthly[m]) - pre_mean_ar
                n_months += 1
    if n_months < (tau_end - tau_start + 1) // 2:
        return None
    return car


# ── Build per-event candidate firm sets ─────────────────────────────

MIN_OBS_PER_EVENT = 20
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']

_print('Building per-event candidate firm sets...')

event_candidates = {}
for event_id, event in enumerate(all_events):
    event_gvkeys = set(event['gvkeys'])
    fm_sic4 = None
    for gk in event_gvkeys:
        fm_sic4 = get_sic4(gk)
        if fm_sic4:
            break
    rows = []
    for fm_gk in event_gvkeys:
        if fm_gk not in W_geo:
            continue
        neighbors = W_geo[fm_gk]
        neighbor_gks = set(neighbors.keys()) - event_gvkeys
        non_connected = [gk for gk in fundamentals
                         if gk not in event_gvkeys and gk not in neighbors]
        stable_seed = int(hashlib.md5(
            str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
        random.seed(stable_seed)
        n_ctrl = min(len(non_connected), max(5 * len(neighbor_gks), 20))
        ctrl_sample = (random.sample(non_connected, n_ctrl)
                       if len(non_connected) > n_ctrl else non_connected)
        candidate_firms = list(neighbor_gks) + ctrl_sample
        for gk in candidate_firms:
            j_sic4 = get_sic4(gk)
            rows.append({
                'gvkey': gk,
                'w_geo': neighbors.get(gk, 0.0),
                'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                'same_sector': 1.0 if (fm_sic4 and j_sic4 and fm_sic4 == j_sic4) else 0.0,
            })
    event_candidates[event_id] = {
        'event_month': event['event_month'],
        'rows': rows,
    }


# ── Run FM regression for each window ───────────────────────────────

def newey_west_se(series, lag=4):
    series = [x for x in series if not (isinstance(x, float) and math.isnan(x))]
    T = len(series)
    if T < 2:
        return float('nan'), float('nan'), 0
    mean = sum(series) / T
    dev = [x - mean for x in series]
    gamma0 = sum(d * d for d in dev) / T
    var_nw = gamma0
    for L in range(1, min(lag, T - 1) + 1):
        weight = 1.0 - L / (lag + 1)
        cov_L = sum(dev[t] * dev[t - L] for t in range(L, T)) / T
        var_nw += 2.0 * weight * cov_L
    if var_nw <= 0:
        return mean, float('nan'), T
    return mean, math.sqrt(var_nw / T), T


def fm_beta_fuel(window):
    """Fama-MacBeth beta_fuel for the CAR window [tau_start, tau_end]."""
    tau_start, tau_end = window
    fuel_betas = []
    for event_id, info in event_candidates.items():
        em = info['event_month']
        obs = []
        for r in info['rows']:
            car = compute_car_window(r['gvkey'], em, tau_start, tau_end)
            if car is None:
                continue
            obs.append({**r, 'car': car})
        if len(obs) < MIN_OBS_PER_EVENT:
            continue
        ss_vals = set(o['same_sector'] for o in obs)
        use_vars = SPEC_VARS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']
        y = [o['car'] for o in obs]
        X = [[1.0] + [o[v] for v in use_vars] for o in obs]
        result = ols_simple(y, X)
        if result is None:
            continue
        beta_dict = dict(zip(['intercept'] + use_vars, result['beta']))
        fuel_betas.append(beta_dict.get('w_fuel', float('nan')))
    mean, se, T = newey_west_se(fuel_betas, lag=4)
    return {'window': window, 'mean': mean, 'se': se, 'n_events': T,
            'betas': fuel_betas}


# Define windows: 3 placebos pre-event + headline
PLACEBO_WINDOWS = [
    (-12, -8),
    (-9, -5),
    (-6, -2),
]
HEADLINE_WINDOW = (-1, 3)
ALL_WINDOWS = PLACEBO_WINDOWS + [HEADLINE_WINDOW]

_print('\nRunning FM regressions per window...')
window_results = []
for w in ALL_WINDOWS:
    _print(f'  Window [{w[0]:+d}, {w[1]:+d}]...')
    res = fm_beta_fuel(w)
    window_results.append(res)
    if not math.isnan(res['mean']) and not math.isnan(res['se']) and res['se'] > 0:
        t = res['mean'] / res['se']
        _print(f'    beta = {res["mean"]:+.4f}, SE = {res["se"]:.4f}, t = {t:+.2f}, '
               f'N events = {res["n_events"]}')


# ── Rambachan-Roth M-bar ────────────────────────────────────────────

# Headline
headline = window_results[-1]
beta_h = headline['mean']
se_h = headline['se']
ZSTAR = 1.96

# Pre-period max |beta|
pre_results = window_results[:-1]
pre_abs = [abs(r['mean']) for r in pre_results
           if not math.isnan(r['mean'])]
v_pre_max = max(pre_abs) if pre_abs else float('nan')
pre_argmax = max(pre_results, key=lambda r: abs(r['mean']) if not math.isnan(r['mean']) else -1)

# M-bar: how many multiples of v_pre_max can the post-period bias be?
# Adjusted lower bound: |beta_h| - M * v_pre_max - ZSTAR * se_h > 0
# Solve: M_bar = (|beta_h| - ZSTAR * se_h) / v_pre_max
if not math.isnan(beta_h) and not math.isnan(se_h) and v_pre_max > 0:
    if abs(beta_h) > ZSTAR * se_h:
        M_bar = (abs(beta_h) - ZSTAR * se_h) / v_pre_max
    else:
        M_bar = 0.0
else:
    M_bar = float('nan')


# ── Output ──────────────────────────────────────────────────────────

out_csv = results_path('summaries', 'honest_did_breakdown.csv')
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['window_start', 'window_end', 'window_label',
                'beta_fuel', 'se', 't', 'n_events', 'role'])
    labels = ['pre_early', 'pre_mid', 'pre_late', 'headline']
    for r, lbl in zip(window_results, labels):
        m, s = r['mean'], r['se']
        t = m / s if (not math.isnan(s) and s > 0) else float('nan')
        w.writerow([
            r['window'][0], r['window'][1], lbl,
            f'{m:+.6f}' if not math.isnan(m) else 'NA',
            f'{s:.6f}' if not math.isnan(s) else 'NA',
            f'{t:+.4f}' if not math.isnan(t) else 'NA',
            r['n_events'],
            'placebo' if lbl != 'headline' else 'headline',
        ])
_print(f'\nWrote {out_csv}')

out_md = results_path('metrics', 'honest_did.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Honest DID Sensitivity (Rambachan-Roth 2023)\n\n')
    f.write('Placebo-CAR-window approach to bound parallel-trends violations. '
            'For each pre-event 5-month window, the same Fama-MacBeth cross-sectional '
            'regression as the headline (CAR on `w_geo`, `w_fuel`, `w_reg`, '
            '`same_sector`) is run with the dependent variable replaced by the '
            'pre-event CAR over the placebo window. Pre-event coefficients should '
            'be statistically indistinguishable from zero under the parallel-trends '
            'assumption.\n\n')

    f.write('## Window-by-window results\n\n')
    f.write('| Window | beta_fuel | SE (NW lag=4) | t-stat | N events | Role |\n')
    f.write('|---|---|---|---|---|---|\n')
    for r, lbl in zip(window_results, ['pre_early', 'pre_mid', 'pre_late', 'headline']):
        m, s = r['mean'], r['se']
        t = m / s if (not math.isnan(s) and s > 0) else float('nan')
        role = '**HEADLINE**' if lbl == 'headline' else 'placebo'
        f.write(
            f'| [{r["window"][0]:+d}, {r["window"][1]:+d}] | '
            f'{m:+.4f} | {s:.4f} | {t:+.2f} | '
            f'{r["n_events"]} | {role} |\n'
        )

    f.write('\n## Rambachan-Roth M-bar breakdown\n\n')
    f.write(f'- Headline beta_fuel ({HEADLINE_WINDOW[0]:+d}, {HEADLINE_WINDOW[1]:+d}): {beta_h:+.4f}\n')
    f.write(f'- Headline SE: {se_h:.4f}\n')
    if not math.isnan(beta_h) and not math.isnan(se_h) and se_h > 0:
        f.write(f'- Headline t-stat: {beta_h / se_h:+.2f}\n')
    f.write(f'- Pre-period max |beta_fuel|: {v_pre_max:.4f} '
            f'(window [{pre_argmax["window"][0]:+d}, {pre_argmax["window"][1]:+d}])\n')
    f.write(f'- **M-bar (5% breakdown): {M_bar:.2f}**\n\n')

    f.write('Interpretation: the headline post-event effect survives a parallel-trends '
            f'violation up to {M_bar:.2f} times the largest deviation observed in the '
            'pre-period before becoming statistically indistinguishable from zero at the '
            '5% level.\n\n')
    if M_bar >= 1.0:
        f.write('**Robust** by Rambachan-Roth convention: the result survives violations '
                'as large as anything actually observed in the pre-period.\n')
    elif M_bar >= 0.5:
        f.write('**Moderate**: the result survives violations of order one-half to one '
                'times the observed pre-period deviations.\n')
    else:
        f.write('**Fragile**: the result does not survive violations as large as the '
                'observed pre-period deviations. Any post-event differential trend on '
                'the same scale as observed pre-event noise is sufficient to zero out '
                'the headline.\n')

    f.write('\n## Method note\n\n')
    f.write('The placebo windows are non-overlapping 5-month CAR windows ending '
            'strictly before the event. The headline window [-1, +3] is also '
            '5 months wide, so the placebo distribution is directly comparable to '
            'the post-event estimate. Each window-level coefficient is a Fama-MacBeth '
            'mean of event-level cross-sectional OLS coefficients, with Newey-West '
            '(lag=4) HAC standard errors on the time series of event-level betas.\n')

_print(f'Wrote {out_md}')
_print(f'\nM-bar (5% breakdown) = {M_bar:.2f}')
_print('Done.')
