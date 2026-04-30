"""Event-time path of the fuel-similarity coefficient.

For each event-time month tau in [-12, +6], runs the same cross-sectional
regression as robust_inference.py but with the dependent variable
replaced by the *single-month* abnormal return at event-month + tau. This
generates the event-time path of beta_fuel(tau), which is the input to:

  1. The event-time plot figure (results/figures/fig_event_time_path.pdf).
  2. The Honest DID sensitivity analysis (honest_did.py).

A clean exposure design predicts beta_fuel(tau) ~= 0 for tau in [-12, -2]
(pre-period) and beta_fuel(tau) < 0 for tau >= 0 (post-period). The pre-event
balance test reported in Section 3.5 of the paper -- t = -1.87 in window
[-5, -2] -- is a coarser version of this exhibit.

Output: results/metrics/event_time_path.md
        results/summaries/event_time_betas.csv
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


# ── Inlined OLS helpers (mirror robust_inference.py) ──────

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


# ── Load data (same as robust_inference.py) ───────────────

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

_print('Loading weight matrices...')
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


_print('Loading events...')
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
_print(f'  First-mover events: {len(all_events)}')


# ── Event-time AR helper ─────────────────────────────────────────────

def compute_event_month(event_month_str, tau):
    """Add tau months to event_month string 'YYYY-MM'."""
    y = int(event_month_str[:4])
    m = int(event_month_str[5:7])
    total = y * 12 + (m - 1) + tau
    return f'{total // 12:04d}-{total % 12 + 1:02d}'


def compute_ar_at_tau(gvkey, event_month_str, tau):
    """Single-month abnormal return at event_month + tau."""
    target = compute_event_month(event_month_str, tau)
    if gvkey not in monthly_ret or target not in monthly_ret[gvkey]:
        return None
    if target not in market_ret_monthly:
        return None
    return monthly_ret[gvkey][target] - market_ret_monthly[target]


# ── Build per-event candidate firm sets (one set per event, reused) ──

_print('\nBuilding per-event candidate firm sets...')

MIN_OBS_PER_EVENT = 20
TAU_MIN = -12
TAU_MAX = 6
TAUS = list(range(TAU_MIN, TAU_MAX + 1))

event_candidates = {}  # event_id -> {'event_month', 'rows': [obs dicts without AR]}

for event_id, event in enumerate(all_events):
    event_gvkeys = set(event['gvkeys'])
    event_month_str = event['event_month']

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

    event_candidates[event_id] = {'event_month': event_month_str, 'rows': rows}

_print(f'  Events with candidate firm sets: {len(event_candidates)}')


# ── Run cross-sectional OLS for each (event, tau) ────────────────────

_print('\nRunning cross-sectional regressions across event-time...')

SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']

# tau -> list of (event_id, beta_fuel, beta_geo, n)
tau_results = {tau: [] for tau in TAUS}

for event_id, info in event_candidates.items():
    event_month_str = info['event_month']
    base_rows = info['rows']

    for tau in TAUS:
        obs = []
        for r in base_rows:
            ar = compute_ar_at_tau(r['gvkey'], event_month_str, tau)
            if ar is None:
                continue
            obs.append({
                'ar': ar,
                'w_geo': r['w_geo'],
                'w_fuel': r['w_fuel'],
                'w_reg': r['w_reg'],
                'same_sector': r['same_sector'],
            })
        if len(obs) < MIN_OBS_PER_EVENT:
            continue

        ss_vals = set(o['same_sector'] for o in obs)
        use_vars = SPEC_VARS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']

        y = [o['ar'] for o in obs]
        X = [[1.0] + [o[v] for v in use_vars] for o in obs]
        result = ols_simple(y, X)
        if result is None:
            continue

        names = ['intercept'] + use_vars
        beta_dict = dict(zip(names, result['beta']))
        tau_results[tau].append({
            'event_id': event_id,
            'beta_fuel': beta_dict.get('w_fuel', float('nan')),
            'beta_geo': beta_dict.get('w_geo', float('nan')),
            'n': result['n'],
        })


# ── Fama-MacBeth + Newey-West aggregation across events ──────────────

def newey_west_se(series, lag=4):
    """Newey-West HAC SE for the mean of a time series."""
    series = [x for x in series if not (isinstance(x, float) and math.isnan(x))]
    T = len(series)
    if T < 2:
        return float('nan'), float('nan')
    mean = sum(series) / T
    dev = [x - mean for x in series]
    gamma0 = sum(d * d for d in dev) / T
    var_nw = gamma0
    for L in range(1, min(lag, T - 1) + 1):
        weight = 1.0 - L / (lag + 1)
        cov_L = sum(dev[t] * dev[t - L] for t in range(L, T)) / T
        var_nw += 2.0 * weight * cov_L
    if var_nw <= 0:
        return mean, float('nan')
    return mean, math.sqrt(var_nw / T)


_print('\nFama-MacBeth aggregation across events at each tau...')

path_summary = []  # list of dicts per tau
for tau in TAUS:
    rs = tau_results[tau]
    if not rs:
        path_summary.append({
            'tau': tau, 'n_events': 0, 'mean_n': 0,
            'beta_fuel': float('nan'), 'se_fuel': float('nan'),
            'beta_geo': float('nan'), 'se_geo': float('nan'),
        })
        continue
    fuel_series = [r['beta_fuel'] for r in rs]
    geo_series = [r['beta_geo'] for r in rs]
    mean_fuel, se_fuel = newey_west_se(fuel_series, lag=4)
    mean_geo, se_geo = newey_west_se(geo_series, lag=4)
    path_summary.append({
        'tau': tau,
        'n_events': len(rs),
        'mean_n': sum(r['n'] for r in rs) / len(rs),
        'beta_fuel': mean_fuel, 'se_fuel': se_fuel,
        'beta_geo': mean_geo, 'se_geo': se_geo,
    })


# ── Output ───────────────────────────────────────────────────────────

# CSV for figure generation
out_csv = results_path('summaries', 'event_time_betas.csv')
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['tau', 'n_events', 'mean_n_per_event',
                'beta_fuel', 'se_fuel', 't_fuel',
                'beta_geo', 'se_geo', 't_geo'])
    for row in path_summary:
        t_fuel = row['beta_fuel'] / row['se_fuel'] if (
            row['se_fuel'] and not math.isnan(row['se_fuel']) and row['se_fuel'] > 0
        ) else float('nan')
        t_geo = row['beta_geo'] / row['se_geo'] if (
            row['se_geo'] and not math.isnan(row['se_geo']) and row['se_geo'] > 0
        ) else float('nan')
        w.writerow([
            row['tau'], row['n_events'], f'{row["mean_n"]:.1f}',
            f'{row["beta_fuel"]:.6f}', f'{row["se_fuel"]:.6f}', f'{t_fuel:.3f}',
            f'{row["beta_geo"]:.6f}', f'{row["se_geo"]:.6f}', f'{t_geo:.3f}',
        ])
_print(f'\n  Wrote {out_csv}')

# Markdown summary
out_md = results_path('metrics', 'event_time_path.md')
os.makedirs(os.path.dirname(out_md), exist_ok=True)
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Event-Time Path of the Fuel-Similarity Coefficient\n\n')
    f.write('Cross-sectional regression of single-month abnormal returns on '
            'network weights (`w_geo`, `w_fuel`, `w_reg`, same-sector indicator), '
            'run separately for each event-time month tau in [-12, +6]. '
            'Beta(tau) is averaged across events using Fama-MacBeth with '
            'Newey-West (lag=4) standard errors.\n\n')
    f.write(f'- Tau range: [{TAU_MIN}, {TAU_MAX}]\n')
    f.write(f'- Events with valid event-month data: {len(event_candidates)}\n')
    f.write(f'- Min obs per event regression: {MIN_OBS_PER_EVENT}\n\n')

    f.write('## Beta_fuel(tau) path\n\n')
    f.write('| tau | n_events | beta_fuel | SE | t-stat | beta_geo | SE | t-stat |\n')
    f.write('|---|---|---|---|---|---|---|---|\n')
    for row in path_summary:
        bf, sf = row['beta_fuel'], row['se_fuel']
        bg, sg = row['beta_geo'], row['se_geo']
        tf = bf / sf if sf and not math.isnan(sf) and sf > 0 else float('nan')
        tg = bg / sg if sg and not math.isnan(sg) and sg > 0 else float('nan')
        f.write(
            f'| {row["tau"]:+d} | {row["n_events"]} | '
            f'{bf:+.4f} | {sf:.4f} | {tf:+.2f} | '
            f'{bg:+.4f} | {sg:.4f} | {tg:+.2f} |\n'
        )

    # Pre/post analysis
    pre_taus = [r for r in path_summary if -12 <= r['tau'] <= -2 and r['n_events'] > 0]
    post_taus = [r for r in path_summary if -1 <= r['tau'] <= 3 and r['n_events'] > 0]
    if pre_taus and post_taus:
        max_pre_abs = max(abs(r['beta_fuel']) for r in pre_taus
                          if not math.isnan(r['beta_fuel']))
        mean_post = sum(r['beta_fuel'] for r in post_taus) / len(post_taus)
        f.write('\n## Pre/post summary\n\n')
        f.write(f'- Pre-period [-12, -2] max |beta_fuel|: {max_pre_abs:.4f}\n')
        f.write(f'- Post-period [-1, +3] mean beta_fuel: {mean_post:+.4f}\n')
        f.write(f'- Ratio (post mean / pre max abs): {mean_post / max_pre_abs:.2f}\n')

_print(f'  Wrote {out_md}')
_print('\nDone.')
