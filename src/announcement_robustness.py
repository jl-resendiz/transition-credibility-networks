"""Announcement-date vs physical-retirement robustness.

The headline pipeline sets t=0 at the announcement date when available
(via EIA-860 disclosures and other public announcements stored in the
`announcement_date` column of `coal_retirement_events.csv`), falling back
to the physical retirement date when no announcement is documented.
Quant-finance referee concern: physical retirement may be a stale
information event, since rational markets price the announcement.

This script reports two contrasting specifications:

  (A) Announcement-when-available (HEADLINE REFERENCE).
      Identical to the existing pipeline: t=0 at announcement_date if
      present, retirement_date otherwise.

  (B) Physical-retirement-only.
      Force t=0 at retirement_date even when announcement_date is
      available. Compare gamma_fuel against (A).

  (C) Announcement-only subsample.
      Restrict to events with announcement_date populated. Tests whether
      the channel is detectable purely from announcement-driven repricing.

The comparison illuminates whether the headline coefficient is driven by
announcement (information) or physical retirement (drift / rebalancing).

Inputs:  data/derived/returns/monthly_returns.csv
         data/derived/events/coal_retirement_events.csv (has both date columns)
         data/derived/networks/weight_matrix_W_*.csv
         data/raw/factors/F-F_Research_Data_Factors.csv

Outputs: results/metrics/announcement_robustness.md
         results/summaries/announcement_robustness.csv
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


# ── OLS helpers ─────────────────────────────────────────────────────

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


def load_market(path):
    out = {}
    with open(path, 'r', encoding='utf-8') as f:
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
            out[f'{parts[0][:4]}-{parts[0][4:6]}'] = (mktrf + rf) / 100.0
    return out


market_ret = load_market(raw_path('factors', 'F-F_Research_Data_Factors.csv'))

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


# Load events with BOTH dates exposed
events_full = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        if row.get('is_first_mover') != 'True':
            continue
        ann = row.get('announcement_date', '').strip()
        ret = row.get('event_date', '').strip()
        events_full.append({
            'plant': row['plant_name'],
            'announcement_date': ann if (ann and len(ann) >= 7) else None,
            'retirement_date': ret if (ret and len(ret) >= 7) else None,
            'gvkeys': row['matched_gvkeys'].split(';'),
        })

n_with_announcement = sum(1 for e in events_full if e['announcement_date'])
n_with_retirement = sum(1 for e in events_full if e['retirement_date'])
_print(f'Events: {len(events_full)} total, {n_with_announcement} with announcement date, '
       f'{n_with_retirement} with retirement date')


# ── Single-factor CAR ────────────────────────────────────────────────

PRE_MONTHS = 24
POST_MONTHS = 3


def compute_car(gvkey, event_month_str):
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
    for offset in range(-1, POST_MONTHS + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months) and months[idx] in monthly_ret[gvkey]:
            m = months[idx]
            if m in market_ret:
                car += (monthly_ret[gvkey][m] - market_ret[m]) - pre_mean
                n += 1
    if n < 3:
        return None
    return car


# ── Build datasets and run FM under each spec ───────────────────────

MIN_OBS = 20
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']


def fm_run(events_to_use, label):
    """Run FM cross-sectional regressions on the given event list."""
    fuel_betas = []
    geo_betas = []
    n_events_used = 0
    for event_id, event in enumerate(events_to_use):
        em = event.get('event_month')
        if not em:
            continue
        event_gvkeys = set(event['gvkeys'])
        fm_sic4 = None
        for gk in event_gvkeys:
            fm_sic4 = get_sic4(gk)
            if fm_sic4:
                break
        obs = []
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
                car = compute_car(gk, em)
                if car is None:
                    continue
                j_sic4 = get_sic4(gk)
                obs.append({
                    'car': car,
                    'w_geo': neighbors.get(gk, 0.0),
                    'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                    'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                    'same_sector': 1.0 if (fm_sic4 and j_sic4 and fm_sic4 == j_sic4) else 0.0,
                })
        if len(obs) < MIN_OBS:
            continue
        ss_vals = set(o['same_sector'] for o in obs)
        use_vars = SPEC_VARS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']
        y = [o['car'] for o in obs]
        X = [[1.0] + [o[v] for v in use_vars] for o in obs]
        result = ols_simple(y, X)
        if result is None:
            continue
        names = ['intercept'] + use_vars
        bd = dict(zip(names, result['beta']))
        fuel_betas.append(bd.get('w_fuel', float('nan')))
        geo_betas.append(bd.get('w_geo', float('nan')))
        n_events_used += 1
    return fuel_betas, geo_betas, n_events_used


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
    return mean, math.sqrt(var / n) if var > 0 else float('nan')


# Build three event lists

events_A = []  # announcement when available
events_B = []  # forced physical
events_C = []  # announcement only subsample
for e in events_full:
    if e['announcement_date'] or e['retirement_date']:
        em_A = (e['announcement_date'] or e['retirement_date'])[:7]
        events_A.append({**e, 'event_month': em_A})
    if e['retirement_date']:
        events_B.append({**e, 'event_month': e['retirement_date'][:7]})
    if e['announcement_date']:
        events_C.append({**e, 'event_month': e['announcement_date'][:7]})

specs = [
    ('A_announce_when_available', events_A,
     'Announcement when available, retirement otherwise (HEADLINE)'),
    ('B_force_physical_retirement', events_B,
     'Force physical retirement date for all events'),
    ('C_announcement_only', events_C,
     'Restrict to events with announcement date'),
]

results = []
for code, evs, desc in specs:
    _print(f'\n[{code}] {desc} (N={len(evs)})')
    fuel, geo, n_used = fm_run(evs, code)
    mean_f, se_f = newey_west_se(fuel, lag=4)
    mean_g, se_g = newey_west_se(geo, lag=4)
    t_f = mean_f / se_f if (not math.isnan(se_f) and se_f > 0) else float('nan')
    t_g = mean_g / se_g if (not math.isnan(se_g) and se_g > 0) else float('nan')
    _print(f'  fuel: beta = {mean_f:+.4f}, SE = {se_f:.4f}, t = {t_f:+.2f}')
    _print(f'  geo:  beta = {mean_g:+.4f}, SE = {se_g:.4f}, t = {t_g:+.2f}')
    _print(f'  n_events_FM: {n_used}')
    results.append({
        'code': code, 'desc': desc, 'n_events_total': len(evs),
        'n_events_FM': n_used,
        'fuel_mean': mean_f, 'fuel_se': se_f, 'fuel_t': t_f,
        'geo_mean': mean_g, 'geo_se': se_g, 'geo_t': t_g,
    })


# ── Output ──────────────────────────────────────────────────────────

out_csv = results_path('summaries', 'announcement_robustness.csv')
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['spec', 'description', 'n_events_total', 'n_events_FM',
                'fuel_mean', 'fuel_se', 'fuel_t',
                'geo_mean', 'geo_se', 'geo_t'])
    for r in results:
        w.writerow([
            r['code'], r['desc'], r['n_events_total'], r['n_events_FM'],
            f'{r["fuel_mean"]:+.6f}',
            f'{r["fuel_se"]:.6f}' if not math.isnan(r['fuel_se']) else 'NA',
            f'{r["fuel_t"]:+.4f}' if not math.isnan(r['fuel_t']) else 'NA',
            f'{r["geo_mean"]:+.6f}',
            f'{r["geo_se"]:.6f}' if not math.isnan(r['geo_se']) else 'NA',
            f'{r["geo_t"]:+.4f}' if not math.isnan(r['geo_t']) else 'NA',
        ])
_print(f'\nWrote {out_csv}')

out_md = results_path('metrics', 'announcement_robustness.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Announcement-Date vs Physical-Retirement Robustness\n\n')
    f.write('Tests whether the headline channel is driven by announcement-date '
            'repricing (information) or by physical retirement (drift / '
            'rebalancing). Three specifications: (A) headline default '
            '(announcement when available, retirement otherwise); (B) forced '
            'physical retirement; (C) announcement-only subsample.\n\n')
    f.write(f'Among {len(events_full)} first-mover events, '
            f'{n_with_announcement} have a populated announcement date, '
            f'{n_with_retirement} have a populated retirement date. '
            'For the headline spec (A), announcement is used whenever present.\n\n')

    f.write('## Specifications\n\n')
    f.write('| Spec | Description | N events (total) | N events (FM) | gamma_fuel | SE | t | gamma_geo | SE | t |\n')
    f.write('|---|---|---|---|---|---|---|---|---|---|\n')
    for r in results:
        f.write(
            f'| {r["code"]} | {r["desc"]} | {r["n_events_total"]} | {r["n_events_FM"]} | '
            f'{r["fuel_mean"]:+.4f} | {r["fuel_se"]:.4f} | {r["fuel_t"]:+.2f} | '
            f'{r["geo_mean"]:+.4f} | {r["geo_se"]:.4f} | {r["geo_t"]:+.2f} |\n'
        )

    f.write('\n## Interpretation\n\n')
    head = results[0]
    forced = results[1]
    ann_only = results[2]
    f.write(f'Headline (A): $\\hat\\gamma_{{\\mathrm{{fuel}}}} = {head["fuel_mean"]:+.2f}$ '
            f'($t={head["fuel_t"]:+.2f}$).\n\n')
    f.write(f'Forced physical retirement (B): $\\hat\\gamma_{{\\mathrm{{fuel}}}} = '
            f'{forced["fuel_mean"]:+.2f}$ ($t={forced["fuel_t"]:+.2f}$).\n\n')
    f.write(f'Announcement-only (C): $\\hat\\gamma_{{\\mathrm{{fuel}}}} = '
            f'{ann_only["fuel_mean"]:+.2f}$ ($t={ann_only["fuel_t"]:+.2f}$).\n\n')

    if (not math.isnan(ann_only['fuel_t']) and abs(ann_only['fuel_t']) >= 1.96
            and ann_only['fuel_mean'] < 0):
        f.write('The announcement-only subsample shows a negative and statistically '
                'significant fuel coefficient, supporting the interpretation that '
                'the channel is driven by announcement-date information.\n')
    else:
        f.write('The announcement-only subsample is small; coefficient direction '
                'and significance should be read with caution.\n')

_print(f'Wrote {out_md}')
_print('\nDone.')
