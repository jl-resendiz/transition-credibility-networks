r"""Anomaly-vs-risk demarcation: post-formation decay test.

Senior referee critique: with corrected magnitudes ($-2.2$ pp/SD over 4
months, $-6.6\%$ annualized), the headline is large enough to invite
"risk premium vs mispricing" pushback. The cleanest demarcation is
post-formation alpha decay:

  - Risk premium hypothesis: CAR accumulates roughly linearly with horizon
    after event. The fuel-similarity coefficient $\beta(H)$ should INCREASE
    in absolute value with H (per-month risk compensation accumulates).

  - Mispricing hypothesis: CAR levels off as institutional arbitrageurs
    correct the mispricing. The fuel-similarity coefficient $\beta(H)$
    should plateau or decay back toward zero at longer horizons.

This script computes $\hat\beta_{\mathrm{fuel}}(H)$ for windows
$[-1, +H]$ at $H \in \{1, 3, 6, 12, 24\}$ months. The pattern of
$\beta(H)$ across $H$ discriminates the two hypotheses.

Inputs:  monthly_returns + weight matrices + events
Outputs: results/metrics/anomaly_vs_risk.md
         results/summaries/anomaly_vs_risk_horizons.csv
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
    return {'beta': beta, 'n': n}


def newey_west_se(series, lag=4):
    series = [x for x in series if not (isinstance(x, float) and math.isnan(x))]
    n = len(series)
    if n < 2:
        return float('nan'), 0
    mean = sum(series) / n
    dev = [x - mean for x in series]
    gamma0 = sum(d * d for d in dev) / n
    var = gamma0
    for L in range(1, min(lag, n - 1) + 1):
        weight = 1.0 - L / (lag + 1)
        cov_L = sum(dev[t] * dev[t - L] for t in range(L, n)) / n
        var += 2.0 * weight * cov_L
    return mean, math.sqrt(var / n) if var > 0 else float('nan'), n


# ── Data loading ────────────────────────────────────────────────────

monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except ValueError:
            pass

market_ret = {}
with open(raw_path('factors', 'F-F_Research_Data_Factors.csv'), 'r', encoding='utf-8') as f:
    for line in f:
        line = line.strip()
        if not line or line.startswith('This file') or line.startswith('The '):
            continue
        if line.startswith(','):
            continue
        parts = [p.strip() for p in line.split(',')]
        if len(parts) < 5 or not parts[0].isdigit() or len(parts[0]) != 6:
            continue
        try:
            mktrf = float(parts[1])
            rf = float(parts[4])
        except ValueError:
            continue
        market_ret[f'{parts[0][:4]}-{parts[0][4:6]}'] = (mktrf + rf) / 100.0

W_geo = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_geo[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

W_fuel = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_fuel.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_fuel[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

W_reg = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_regulatory.csv'), 'r', encoding='utf-8') as f:
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
        ann = row.get('announcement_date', '').strip()
        ret = row.get('event_date', '').strip()
        ed = ann if (ann and len(ann) >= 7) else ret
        if not ed or len(ed) < 7:
            continue
        all_events.append({'event_month': ed[:7],
                          'gvkeys': row['matched_gvkeys'].split(';')})


# ── CAR over [-1, +H] ───────────────────────────────────────────────

PRE_MONTHS = 24


def compute_car_horizon(gvkey, em, post_h):
    if gvkey not in monthly_ret:
        return None
    months = sorted(monthly_ret[gvkey].keys())
    event_idx = None
    for i, m in enumerate(months):
        if m >= em:
            event_idx = i
            break
    if event_idx is None:
        return None
    ar_list = []
    for i in range(max(0, event_idx - PRE_MONTHS), event_idx):
        m = months[i]
        if m in monthly_ret[gvkey] and m in market_ret:
            ar_list.append(monthly_ret[gvkey][m] - market_ret[m])
    if len(ar_list) < 12:
        return None
    pre_mean = sum(ar_list) / len(ar_list)
    car = 0.0
    n = 0
    for offset in range(-1, post_h + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months) and months[idx] in monthly_ret[gvkey]:
            m = months[idx]
            if m in market_ret:
                car += (monthly_ret[gvkey][m] - market_ret[m]) - pre_mean
                n += 1
    if n < (post_h + 1) // 2:
        return None
    return car


# ── Build per-event candidate sets and compute CARs at each horizon ─

MIN_OBS = 20
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']
HORIZONS = [1, 3, 6, 12, 24]


_print('Building per-event candidate sets...')
event_candidates = {}
for event_id, event in enumerate(all_events):
    em = event['event_month']
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
        candidates = list(neighbor_gks) + ctrl_sample
        for gk in candidates:
            j_sic4 = get_sic4(gk)
            rows.append({
                'gvkey': gk,
                'w_geo': neighbors.get(gk, 0.0),
                'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                'same_sector': 1.0 if (fm_sic4 and j_sic4 and fm_sic4 == j_sic4) else 0.0,
            })
    event_candidates[event_id] = {'event_month': em, 'rows': rows}


# ── Run FM regression at each horizon ───────────────────────────────

_print('\nRunning FM regressions at each post-event horizon...')
horizon_results = []
for H in HORIZONS:
    _print(f'  Horizon [-1, +{H}]...')
    fuel_betas = []
    for event_id, info in event_candidates.items():
        em = info['event_month']
        obs = []
        for r in info['rows']:
            car = compute_car_horizon(r['gvkey'], em, H)
            if car is None:
                continue
            obs.append({**r, 'car': car})
        if len(obs) < MIN_OBS:
            continue
        ss_vals = set(o['same_sector'] for o in obs)
        use_vars = SPEC_VARS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']
        y = [o['car'] for o in obs]
        X = [[1.0] + [o[v] for v in use_vars] for o in obs]
        result = ols_simple(y, X)
        if result is None:
            continue
        bd = dict(zip(['intercept'] + use_vars, result['beta']))
        fuel_betas.append(bd.get('w_fuel', float('nan')))
    mean, se, T = newey_west_se(fuel_betas, lag=4)
    horizon_results.append({
        'H': H, 'window': f'[-1, +{H}]',
        'mean': mean, 'se': se, 'n_events': T,
        't': mean / se if (not math.isnan(se) and se > 0) else float('nan'),
    })
    _print(f'    beta = {mean:+.4f}, SE = {se:.4f}, t = {mean/se:+.2f}, '
           f'N = {T}')


# ── Output ──────────────────────────────────────────────────────────

out_csv = results_path('summaries', 'anomaly_vs_risk_horizons.csv')
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['horizon_months', 'window',
                'fuel_mean', 'fuel_se', 'fuel_t',
                'fuel_per_month', 'n_events_FM'])
    for r in horizon_results:
        bp = r['mean']
        per_m = bp / (r['H'] + 1) if r['H'] >= 0 else float('nan')
        w.writerow([
            r['H'], r['window'],
            f'{r["mean"]:+.6f}',
            f'{r["se"]:.6f}' if not math.isnan(r['se']) else 'NA',
            f'{r["t"]:+.4f}' if not math.isnan(r['t']) else 'NA',
            f'{per_m:+.6f}',
            r['n_events'],
        ])
_print(f'\nWrote {out_csv}')

out_md = results_path('metrics', 'anomaly_vs_risk.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Anomaly vs Risk: Post-Formation Decay Test\n\n')
    f.write('Tests whether the fuel-similarity premium accumulates with '
            'horizon (risk-premium hypothesis) or plateaus / decays '
            '(mispricing hypothesis).\n\n')

    f.write('## Cumulative effect by post-event horizon\n\n')
    f.write('| Window | $\\hat\\beta_{\\mathrm{fuel}}$ | SE (NW lag=4) | t-stat | per-month | N events |\n')
    f.write('|---|---|---|---|---|---|\n')
    for r in horizon_results:
        per_m = r['mean'] / (r['H'] + 1) if r['H'] >= 0 else float('nan')
        f.write(f'| {r["window"]} | {r["mean"]:+.4f} | {r["se"]:.4f} | '
                f'{r["t"]:+.2f} | {per_m:+.4f} | {r["n_events"]} |\n')

    f.write('\n## Interpretation\n\n')
    if len(horizon_results) >= 3:
        b1 = horizon_results[0]['mean']  # H=1
        b3 = horizon_results[1]['mean']  # H=3
        b6 = horizon_results[2]['mean'] if len(horizon_results) > 2 else float('nan')
        b12 = horizon_results[3]['mean'] if len(horizon_results) > 3 else float('nan')
        b24 = horizon_results[4]['mean'] if len(horizon_results) > 4 else float('nan')

        f.write(f'- $\\hat\\beta(1) = {b1:+.4f}$, $\\hat\\beta(3) = {b3:+.4f}$, '
                f'$\\hat\\beta(6) = {b6:+.4f}$, $\\hat\\beta(12) = {b12:+.4f}$, '
                f'$\\hat\\beta(24) = {b24:+.4f}$.\n\n')

        if not (math.isnan(b3) or math.isnan(b24)):
            ratio = b24 / b3 if b3 != 0 else float('nan')
            f.write(f'- Ratio $\\hat\\beta(24)/\\hat\\beta(3) = {ratio:.2f}$. ')

            # Decision logic
            # Risk: ratio >> 1 (effect grows with horizon)
            # Persistent: ratio ~ 1 (effect plateaus, no decay)
            # Mispricing: ratio < 1 (effect decays)
            if ratio > 1.5:
                f.write('The effect grows with horizon, consistent with a '
                        'risk-premium interpretation: longer horizons accumulate '
                        'more of the priced compensation for transition exposure.\n')
            elif ratio > 0.7:
                f.write('The effect roughly persists with horizon (within '
                        'sampling noise), neither growing nor decaying. This is '
                        'consistent with either a slow-moving risk channel or '
                        'a non-correctable persistent mispricing.\n')
            else:
                f.write('The effect decays at longer horizons, consistent with '
                        'a mispricing interpretation: institutional arbitrage '
                        'corrects the initial under-reaction over time.\n')

    f.write('\n## Caveats\n\n')
    f.write('- This is one of three demarcation tests. The institutional '
            'ownership split and the characteristic-matched (DGTW) benchmark, '
            'flagged by referees as natural complements, require additional '
            'data not included in the replication package; these are reserved '
            'for future versions.\n')
    f.write('- Longer horizons are subject to greater confounding from '
            'subsequent retirement events, so the H=24 result should be '
            'read as suggestive rather than definitive.\n')

_print(f'Wrote {out_md}')
_print('\nDone.')
