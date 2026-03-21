"""Daily CAR robustness for phase-out shocks:
  - Main regression (CAR ~ exposure) with event-clustered SEs
  - Leave-one-out by event
  - Placebo timing shifts (+/- months)

Outputs:
  - results/metrics/strategy3_phaseout_car_robustness_{suffix}.md
"""
import csv
import math
import os
import bisect
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

EVENTS_PATH = os.getenv('EVENTS_PATH', derived_path('events', 'coal_phaseout_shocks_events.csv'))
TIER_FILTER = os.getenv('TIER_FILTER', '')
BINDING_ONLY = os.getenv('BINDING_ONLY', '0') == '1'
TRANSFORM_SET = os.getenv('TRANSFORM_SET', 'log1p')  # base | log1p | zscore
CONTROL_MULT = int(os.getenv('CONTROL_MULT', '5'))

CAR_START = int(os.getenv('CAR_START', '-1'))
CAR_END = int(os.getenv('CAR_END', '20'))
MIN_CAR_FRAC = float(os.getenv('MIN_CAR_FRAC', '0.4'))

PLACEBO_SHIFTS = os.getenv('PLACEBO_SHIFTS', '-6,6')
SHORT_WINDOW = os.getenv('PLACEBO_SHORT', '1') == '1'

METRICS_SUFFIX = os.getenv('METRICS_SUFFIX', 'phaseout')


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
            vwretd[f'{date[:4]}-{date[4:6]}-{date[6:]}'] = (mktrf_val + rf_val) / 100.0
    return vwretd


def shift_date(date_str, months):
    # date_str: YYYY-MM-DD
    y = int(date_str[:4])
    m = int(date_str[5:7])
    d = int(date_str[8:10])
    m += months
    while m > 12:
        y += 1
        m -= 12
    while m < 1:
        y -= 1
        m += 12
    # clamp day to 28 for simplicity
    d = min(d, 28)
    return f'{y:04d}-{m:02d}-{d:02d}'


def apply_transform(vals, transform):
    if transform == 'base':
        return vals
    if transform == 'log1p':
        return [math.log1p(v) if v > 0 else 0.0 for v in vals]
    if transform == 'zscore':
        mean = sum(vals) / len(vals) if vals else 0.0
        var = sum((v - mean) ** 2 for v in vals) / len(vals) if vals else 0.0
        std = math.sqrt(var) if var > 1e-12 else 1.0
        return [(v - mean) / std for v in vals]
    return vals


def ols_cluster(y, x, clusters):
    n = len(y)
    if n == 0:
        return None
    # x is list of [1, w]
    XtX = [[0.0, 0.0], [0.0, 0.0]]
    Xty = [0.0, 0.0]
    for i in range(n):
        xi0, xi1 = x[i]
        XtX[0][0] += xi0 * xi0
        XtX[0][1] += xi0 * xi1
        XtX[1][0] += xi1 * xi0
        XtX[1][1] += xi1 * xi1
        Xty[0] += xi0 * y[i]
        Xty[1] += xi1 * y[i]
    det = XtX[0][0] * XtX[1][1] - XtX[0][1] * XtX[1][0]
    if abs(det) < 1e-12:
        return None
    inv = [[XtX[1][1] / det, -XtX[0][1] / det],
           [-XtX[1][0] / det, XtX[0][0] / det]]
    beta0 = inv[0][0] * Xty[0] + inv[0][1] * Xty[1]
    beta1 = inv[1][0] * Xty[0] + inv[1][1] * Xty[1]
    resid = [y[i] - beta0 * x[i][0] - beta1 * x[i][1] for i in range(n)]
    clus = defaultdict(list)
    for i, cid in enumerate(clusters):
        clus[cid].append(i)
    S = [[0.0, 0.0], [0.0, 0.0]]
    for _, idxs in clus.items():
        xu0 = 0.0
        xu1 = 0.0
        for i in idxs:
            xu0 += x[i][0] * resid[i]
            xu1 += x[i][1] * resid[i]
        S[0][0] += xu0 * xu0
        S[0][1] += xu0 * xu1
        S[1][0] += xu1 * xu0
        S[1][1] += xu1 * xu1
    cov = [
        [inv[0][0] * S[0][0] * inv[0][0] + inv[0][1] * S[1][0] * inv[0][0] + inv[0][0] * S[0][1] * inv[0][1] + inv[0][1] * S[1][1] * inv[0][1],
         inv[0][0] * S[0][0] * inv[1][0] + inv[0][1] * S[1][0] * inv[1][0] + inv[0][0] * S[0][1] * inv[1][1] + inv[0][1] * S[1][1] * inv[1][1]],
        [inv[1][0] * S[0][0] * inv[0][0] + inv[1][1] * S[1][0] * inv[0][0] + inv[1][0] * S[0][1] * inv[0][1] + inv[1][1] * S[1][1] * inv[0][1],
         inv[1][0] * S[0][0] * inv[1][0] + inv[1][1] * S[1][0] * inv[1][0] + inv[1][0] * S[0][1] * inv[1][1] + inv[1][1] * S[1][1] * inv[1][1]]
    ]
    se1 = math.sqrt(cov[1][1]) if cov[1][1] > 0 else float('nan')
    t1 = beta1 / se1 if se1 and se1 > 0 else float('nan')
    return beta1, se1, t1, len(clus), n


def build_dataset(events, daily_ar, W, fundamentals, event_dates, car_start, car_end):
    obs = []
    for event_id, e in enumerate(events):
        event_date = event_dates[event_id]
        event_gvkeys = set(e['gvkeys'])
        for fm_gk in event_gvkeys:
            if fm_gk not in W:
                continue
            neighbors = W[fm_gk]
            neighbor_gks = set(neighbors.keys()) - event_gvkeys
            non_connected = [gk for gk in fundamentals if gk not in event_gvkeys and gk not in neighbors]
            n_ctrl = min(len(non_connected), max(CONTROL_MULT * len(neighbor_gks), 20))
            ctrl_sample = non_connected[:n_ctrl]
            candidate_firms = list(neighbor_gks) + ctrl_sample
            for gk in candidate_firms:
                w = neighbors.get(gk, 0.0)
                if gk not in daily_ar:
                    continue
                dates, ar = daily_ar[gk]
                # find event index
                idx = bisect.bisect_left(dates, event_date)
                if idx is None or idx >= len(dates):
                    continue
                # window
                vals = []
                for offset in range(car_start, car_end + 1):
                    j = idx + offset
                    if 0 <= j < len(dates):
                        vals.append(ar[j])
                if len(vals) < (car_end - car_start + 1) * MIN_CAR_FRAC:
                    continue
                car = sum(vals)
                obs.append({
                    'event_id': event_id,
                    'w': w,
                    'car': car,
                })
    return obs


# Load daily returns and compute AR
market = load_ff_factors_daily(raw_path('factors', 'F-F_Research_Data_Factors_daily.csv'))
if not market:
    raise RuntimeError('Missing F-F daily factors for vwretd.')

daily_ret = defaultdict(dict)
with open(derived_path('returns', 'daily_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        date = row['datadate']
        try:
            daily_ret[gk][date] = float(row['ret_daily'])
        except (ValueError, TypeError):
            continue

daily_ar = {}
for gk, dct in daily_ret.items():
    dates = [d for d in dct.keys() if d in market]
    if not dates:
        continue
    dates.sort()
    ar = [dct[d] - market[d] for d in dates]
    daily_ar[gk] = (dates, ar)

# Load W
W = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

# Fundamentals for controls
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

# Events
tiers = set([t.strip() for t in TIER_FILTER.split(',') if t.strip()]) if TIER_FILTER else None
events = []
event_dates = []
with open(EVENTS_PATH, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        if tiers and row.get('exogeneity_tier', '') not in tiers:
            continue
        if BINDING_ONLY and row.get('binding', '').strip().lower() != 'yes':
            continue
        event_date = row.get('event_date', '')
        if not event_date or len(event_date) < 10:
            continue
        events.append({
            'gvkeys': [g for g in row['matched_gvkeys'].split(';') if g],
        })
        event_dates.append(event_date)

# Build base dataset
obs = build_dataset(events, daily_ar, W, fundamentals, event_dates, CAR_START, CAR_END)

# Exposure transform
w_raw = [o['w'] for o in obs]
w_t = apply_transform(w_raw, TRANSFORM_SET)
for i, o in enumerate(obs):
    o['w_t'] = w_t[i]

# Main regression
X = [[1.0, o['w_t']] for o in obs]
y = [o['car'] for o in obs]
clusters = [o['event_id'] for o in obs]
main = ols_cluster(y, X, clusters)

# Leave-one-out by event
loo = []
for eid in range(len(events)):
    idx = [i for i, o in enumerate(obs) if o['event_id'] != eid]
    if len(idx) < 10:
        continue
    X_ = [X[i] for i in idx]
    y_ = [y[i] for i in idx]
    c_ = [clusters[i] for i in idx]
    res = ols_cluster(y_, X_, c_)
    if res:
        loo.append(res[0:3])  # beta, se, t

# Placebos ±6 months
placebo_results = []
shifts = [int(x.strip()) for x in PLACEBO_SHIFTS.split(',') if x.strip()]
for shift in shifts:
    shifted = [shift_date(d, shift) for d in event_dates]
    car_start = CAR_START
    car_end = CAR_END
    if SHORT_WINDOW:
        car_start = -1
        car_end = 5
    obs_p = build_dataset(events, daily_ar, W, fundamentals, shifted, car_start, car_end)
    if not obs_p:
        placebo_results.append((shift, None))
        continue
    w_raw_p = [o['w'] for o in obs_p]
    w_t_p = apply_transform(w_raw_p, TRANSFORM_SET)
    Xp = [[1.0, w_t_p[i]] for i in range(len(w_t_p))]
    yp = [o['car'] for o in obs_p]
    cp = [o['event_id'] for o in obs_p]
    res = ols_cluster(yp, Xp, cp)
    placebo_results.append((shift, res))

# Write summary
out_path = results_path('metrics', f'strategy3_phaseout_car_robustness_{METRICS_SUFFIX}.md')
lines = [
    '# Phase-out CAR Robustness',
    '',
    f'- events: {len(events)} (tier_filter={TIER_FILTER}, binding_only={BINDING_ONLY})',
    f'- transform: {TRANSFORM_SET}',
    f'- CAR window: [{CAR_START},{CAR_END}]',
    f'- N obs: {len(obs)}',
    '',
]
if main:
    lines += [
        '## Main regression (event-clustered)',
        f'beta={main[0]:+.4f}, se={main[1]:.4f}, t={main[2]:.2f}, clusters={main[3]}, N={main[4]}',
        '',
    ]

if loo:
    betas = [b for b, _, _ in loo]
    tstats = [t for _, _, t in loo]
    lines += [
        '## Leave-one-out (event)',
        f'- beta min={min(betas):+.4f}, max={max(betas):+.4f}, mean={sum(betas)/len(betas):+.4f}',
        f'- t min={min(tstats):+.2f}, max={max(tstats):+.2f}, mean={sum(tstats)/len(tstats):+.2f}',
        '',
    ]

lines.append('## Placebo timing (shifted months)')
for shift, res in placebo_results:
    if res:
        lines.append(f'- shift {shift:+d}: beta={res[0]:+.4f}, se={res[1]:.4f}, t={res[2]:.2f}, clusters={res[3]}, N={res[4]}')
    else:
        lines.append(f'- shift {shift:+d}: no data')

os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))
print(f'Wrote: {out_path}')
