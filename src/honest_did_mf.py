"""Honest DID sensitivity (Rambachan-Roth 2023) on MULTI-FACTOR abnormal returns.

Recalibrates the placebo-CAR-window analysis from `honest_did.py`
using the 4-factor (FF3 + sample-constructed utility industry) abnormal
returns, replacing the single-factor market-adjusted CARs.

The original Honest DID gave M-bar = 1.26 on single-factor CARs. This script
re-runs the placebo windows ([-12,-8], [-9,-5], [-6,-2], [-1,+3]) with
multi-factor-adjusted CARs and reports the recalibrated M-bar.

Inputs:  same as multifactor_inference.py
Outputs: results/metrics/honest_did_mf.md
         results/summaries/honest_did_breakdown_mf.csv
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


# ── OLS helpers (inlined) ────────────────────────────────────────────

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
    XtX = [[sum(X_mat[i][a] * X_mat[i][b] for i in range(n))
            for b in range(k)] for a in range(k)]
    Xty = [sum(X_mat[i][a] * y[i] for i in range(n)) for a in range(k)]
    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None
    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]
    y_hat = [sum(X_mat[i][a] * beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]
    return {'beta': beta, 'resid': resid, 'n': n}


# ── Data loading (mirrors multifactor_inference.py) ────────

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


def load_ff3_monthly(path):
    factors = {}
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
                mkt_rf = float(parts[1]) / 100.0
                smb = float(parts[2]) / 100.0
                hml = float(parts[3]) / 100.0
                rf = float(parts[4]) / 100.0
            except ValueError:
                continue
            factors[f'{date[:4]}-{date[4:6]}'] = {
                'mkt_rf': mkt_rf, 'smb': smb, 'hml': hml, 'rf': rf,
            }
    return factors


ff3 = load_ff3_monthly(raw_path('factors', 'F-F_Research_Data_Factors.csv'))

month_returns = defaultdict(list)
for gk, dates in monthly_ret.items():
    for m, r in dates.items():
        month_returns[m].append(r)
util_factor = {m: sum(rs) / len(rs)
               for m, rs in month_returns.items() if len(rs) >= 30}
factors_panel = {
    m: {**ff, 'utl_excess': util_factor[m] - ff['rf']}
    for m, ff in ff3.items() if m in util_factor
}

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
            'event_month': effective_date[:7],
            'gvkeys': row['matched_gvkeys'].split(';'),
        })


# ── Multi-factor CAR over arbitrary window ──────────────────────────

PRE_MONTHS = 24
PRE_MIN = 12


def compute_car_multifactor_window(gvkey, event_month_str, tau_start, tau_end):
    """4-factor abnormal CAR over [tau_start, tau_end]."""
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

    # Pre-event window 4-factor regression
    y_pre, X_pre = [], []
    for i in range(max(0, event_idx - PRE_MONTHS), event_idx):
        m = months[i]
        if m in monthly_ret[gvkey] and m in factors_panel:
            f = factors_panel[m]
            y_pre.append(monthly_ret[gvkey][m] - f['rf'])
            X_pre.append([1.0, f['mkt_rf'], f['smb'], f['hml'], f['utl_excess']])
    if len(y_pre) < PRE_MIN:
        return None
    result = ols_simple(y_pre, X_pre)
    if result is None:
        return None
    alpha, b_mkt, b_smb, b_hml, b_utl = result['beta']

    width = tau_end - tau_start + 1
    car = 0.0
    n_months = 0
    for offset in range(tau_start, tau_end + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months):
            m = months[idx]
            if m in monthly_ret[gvkey] and m in factors_panel:
                f = factors_panel[m]
                expected = (alpha + b_mkt * f['mkt_rf']
                            + b_smb * f['smb'] + b_hml * f['hml']
                            + b_utl * f['utl_excess'])
                actual = monthly_ret[gvkey][m] - f['rf']
                car += actual - expected
                n_months += 1
    if n_months < width // 2 + 1:
        return None
    return car


# ── Build per-event candidate firm sets ─────────────────────────────

MIN_OBS_PER_EVENT = 20
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']

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


# ── Newey-West aggregation across events ────────────────────────────

def newey_west_se(series, lag=4):
    series = [x for x in series if not (isinstance(x, float) and math.isnan(x))]
    n = len(series)
    if n < 2:
        return float('nan'), 0
    mean = sum(series) / n
    dev = [x - mean for x in series]
    gamma0 = sum(d * d for d in dev) / n
    var_nw = gamma0
    for L in range(1, min(lag, n - 1) + 1):
        weight = 1.0 - L / (lag + 1)
        cov_L = sum(dev[t] * dev[t - L] for t in range(L, n)) / n
        var_nw += 2.0 * weight * cov_L
    if var_nw <= 0:
        return mean, float('nan')
    return mean, math.sqrt(var_nw / n)


def fm_beta_fuel(window):
    """Fama-MacBeth beta_fuel for the multi-factor CAR window."""
    tau_start, tau_end = window
    fuel_betas = []
    for event_id, info in event_candidates.items():
        em = info['event_month']
        obs = []
        for r in info['rows']:
            car = compute_car_multifactor_window(r['gvkey'], em, tau_start, tau_end)
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
    mean, se = newey_west_se(fuel_betas, lag=4)
    return {'window': window, 'mean': mean, 'se': se, 'n_events': len(fuel_betas)}


PLACEBO_WINDOWS = [(-12, -8), (-9, -5), (-6, -2)]
HEADLINE_WINDOW = (-1, 3)
ALL_WINDOWS = PLACEBO_WINDOWS + [HEADLINE_WINDOW]

_print('\nRunning multi-factor FM regressions per window...')
window_results = []
for w in ALL_WINDOWS:
    _print(f'  Window [{w[0]:+d}, {w[1]:+d}]...')
    res = fm_beta_fuel(w)
    window_results.append(res)
    if not math.isnan(res['mean']) and not math.isnan(res['se']) and res['se'] > 0:
        t = res['mean'] / res['se']
        _print(f'    beta = {res["mean"]:+.4f}, SE = {res["se"]:.4f}, '
               f't = {t:+.2f}, N events = {res["n_events"]}')


# ── M-bar breakdown ─────────────────────────────────────────────────

headline = window_results[-1]
beta_h = headline['mean']
se_h = headline['se']
ZSTAR = 1.96

pre_results = window_results[:-1]
pre_abs = [abs(r['mean']) for r in pre_results if not math.isnan(r['mean'])]
v_pre_max = max(pre_abs) if pre_abs else float('nan')
pre_argmax = max(pre_results,
                 key=lambda r: abs(r['mean']) if not math.isnan(r['mean']) else -1)

if not math.isnan(beta_h) and not math.isnan(se_h) and v_pre_max > 0:
    if abs(beta_h) > ZSTAR * se_h:
        M_bar = (abs(beta_h) - ZSTAR * se_h) / v_pre_max
    else:
        M_bar = 0.0
else:
    M_bar = float('nan')


# ── Output ──────────────────────────────────────────────────────────

out_csv = results_path('summaries', 'honest_did_breakdown_mf.csv')
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

out_md = results_path('metrics', 'honest_did_mf.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Honest DID Sensitivity (Multi-Factor CARs)\n\n')
    f.write('Recalibration of the placebo-CAR-window Rambachan-Roth analysis '
            'using 4-factor (FF3 + sample-constructed utility industry) '
            'abnormal returns, replacing the single-factor market-adjusted CARs '
            'used in `honest_did.md`.\n\n')

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

    f.write('\n## Recalibrated M-bar\n\n')
    f.write(f'- Headline beta_fuel ({HEADLINE_WINDOW[0]:+d}, {HEADLINE_WINDOW[1]:+d}): {beta_h:+.4f}\n')
    f.write(f'- Headline SE: {se_h:.4f}\n')
    if not math.isnan(beta_h) and not math.isnan(se_h) and se_h > 0:
        f.write(f'- Headline t-stat: {beta_h / se_h:+.2f}\n')
    f.write(f'- Pre-period max |beta_fuel|: {v_pre_max:.4f} '
            f'(window [{pre_argmax["window"][0]:+d}, {pre_argmax["window"][1]:+d}])\n')
    f.write(f'- **M-bar (5% breakdown, multi-factor): {M_bar:.2f}**\n\n')
    f.write(f'For comparison, M-bar on single-factor CARs was 1.26 '
            '(reported in `honest_did.md`).\n\n')
    if M_bar >= 1.0:
        f.write('**Robust** by Rambachan-Roth convention: the multi-factor result '
                'survives violations as large as anything observed in the multi-factor '
                'pre-period.\n')
    elif M_bar >= 0.5:
        f.write('**Moderate** robustness under the multi-factor CARs.\n')
    else:
        f.write('**Fragile** under multi-factor adjustment: the result does not survive '
                'violations as large as the observed pre-period deviations.\n')

_print(f'Wrote {out_md}')
_print(f'\nM-bar (multi-factor) = {M_bar:.2f}')
_print('Done.')
