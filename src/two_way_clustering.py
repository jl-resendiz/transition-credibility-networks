"""Cameron-Gelbach-Miller (2011) two-way clustered standard errors for the
channel decomposition (Table 2 column 3).

Why this script exists
----------------------
The headline pooled-OLS coefficients are computed by `joint_tests.jl` (Julia,
fast for the F-test bootstrap). The Julia version reports event-clustered SEs
only. The manuscript also reports two-way (event + firm) clustered SEs in
Table 2 column 3 and in §3.5 ("two-way clustered $t = 4.32$"). Those SEs are
implemented in `joint_tests.py` but the Python script is only run as a
fallback when Julia is unavailable, so the two-way numbers fall out of the
DAG when Julia is on PATH.

This script writes the two-way clustered SEs to its own metric file
(`results/metrics/two_way_clustering.md`) regardless of which version of
`joint_tests` ran upstream, restoring full reproducibility for Table 2.

Specification (identical to joint_tests):
  CAR_j = alpha + beta_geo*w^geo + beta_fuel*w^fuel
        + beta_reg*w^reg + beta_s*SameSector + eps_j
"""
import csv
import json
import math
import os
import random
import hashlib
import sys
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

# ── Configuration (matches joint_tests.py and joint_tests.jl) ─────────

POST_MONTHS = 3
PRE_MONTHS = 24
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']
CHANNEL_VARS = ['w_geo', 'w_fuel', 'w_reg']


def _print(msg=''):
    print(msg)
    sys.stdout.flush()


# ── Linear-algebra helpers ────────────────────────────────────────────

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


def mat_mul(a, b):
    rows_a, mid, cols_b = len(a), len(b), len(b[0])
    out = [[0.0] * cols_b for _ in range(rows_a)]
    for i in range(rows_a):
        for k in range(mid):
            aik = a[i][k]
            if aik == 0.0:
                continue
            for j in range(cols_b):
                out[i][j] += aik * b[k][j]
    return out


# ── OLS with full vcov ────────────────────────────────────────────────

def ols(data, y_var, x_vars):
    n = len(data)
    k = len(x_vars) + 1
    if n <= k + 1:
        return None
    y = [d[y_var] for d in data]
    X = [[1.0] + [d[xv] for xv in x_vars] for d in data]
    XtX = [[sum(X[i][a] * X[i][b] for i in range(n))
            for b in range(k)] for a in range(k)]
    Xty = [sum(X[i][a] * y[i] for i in range(n)) for a in range(k)]
    inv_XtX = invert_matrix(XtX)
    if inv_XtX is None:
        return None
    beta = [sum(inv_XtX[a][b] * Xty[b] for b in range(k)) for a in range(k)]
    y_hat = [sum(X[i][a] * beta[a] for a in range(k)) for i in range(n)]
    resid = [y[i] - y_hat[i] for i in range(n)]
    return {'beta': beta, 'X': X, 'resid': resid, 'inv_XtX': inv_XtX,
            'n': n, 'k': k, 'names': ['intercept'] + x_vars}


# ── Cameron-Gelbach-Miller two-way clustering ────────────────────────

def _cluster_meat(X, resid, cluster_keys, k):
    """Σ_g (X_g' u_g)(u_g' X_g) for a single clustering dimension."""
    cluster_map = defaultdict(list)
    for i, c in enumerate(cluster_keys):
        cluster_map[c].append(i)
    S = [[0.0] * k for _ in range(k)]
    for idxs in cluster_map.values():
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
    return S, len(cluster_map)


def twoway_cgm_vcov(X, resid, event_ids, firm_ids, inv_XtX, k, n):
    """V_twoway = V_event + V_firm − V_(event×firm)."""
    S_event, G_event = _cluster_meat(X, resid, event_ids, k)
    S_firm, G_firm = _cluster_meat(X, resid, firm_ids, k)
    S_ef, G_ef = _cluster_meat(X, resid, list(zip(event_ids, firm_ids)), k)

    def vmat(S, G):
        if G < 2:
            return [[0.0] * k for _ in range(k)]
        V = mat_mul(mat_mul(inv_XtX, S), inv_XtX)
        scale = (G / (G - 1)) * ((n - 1) / (n - k))
        for a in range(k):
            for b in range(k):
                V[a][b] *= scale
        return V

    V_event = vmat(S_event, G_event)
    V_firm = vmat(S_firm, G_firm)
    V_ef = vmat(S_ef, G_ef)
    V_tw = [[V_event[a][b] + V_firm[a][b] - V_ef[a][b]
             for b in range(k)] for a in range(k)]
    # Fall back to max(V_event, V_firm) if a diagonal entry is non-positive
    for a in range(k):
        if V_tw[a][a] <= 0:
            V_tw[a][a] = max(V_event[a][a], V_firm[a][a])
    return V_tw, V_event, V_firm, G_event, G_firm


# ── Two-sided p-value via normal CDF ─────────────────────────────────

def normal_cdf(x):
    if x < -8:
        return 0.0
    if x > 8:
        return 1.0
    ax = abs(x)
    b0, b1, b2, b3, b4, b5 = (0.2316419, 0.319381530, -0.356563782,
                              1.781477937, -1.821255978, 1.330274429)
    t = 1.0 / (1.0 + b0 * ax)
    phi = (1.0 / math.sqrt(2.0 * math.pi)) * math.exp(-0.5 * ax * ax)
    cdf = 1.0 - phi * (b1*t + b2*t**2 + b3*t**3 + b4*t**4 + b5*t**5)
    return cdf if x >= 0 else 1.0 - cdf


def p_two_sided(t):
    return 2.0 * (1.0 - normal_cdf(abs(t)))


def stars(p):
    if p < 0.01:
        return '***'
    if p < 0.05:
        return '**'
    if p < 0.10:
        return '*'
    return ''


# ── Data loading (matches joint_tests.py) ────────────────────────────

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
_print(f'  Monthly: {len(monthly_ret)} firms')

_print('Loading Fama-French factors...')
market_ret_monthly = {}
ff_path = raw_path('factors', 'F-F_Research_Data_Factors.csv')
with open(ff_path, 'r', encoding='utf-8') as f:
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
        market_ret_monthly[f'{date[:4]}-{date[4:6]}'] = (mktrf_val + rf_val) / 100.0
_print(f'  Market months: {len(market_ret_monthly)}')

_print('Loading weight matrices...')
W_geo = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        W_geo[row['gvkey_i']][row['gvkey_j']] = float(row['w_ij'])

W_fuel = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_fuel.csv'), 'r', encoding='utf-8') as f:
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
    if f and f.get('sic'):
        return f['sic'][:4]
    return None


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
        if effective_date and len(effective_date) >= 4 and effective_date[:4].isdigit():
            event_year = int(effective_date[:4])
        else:
            event_year = int(row['ret_year']) if row.get('ret_year') else None
        all_events.append({
            'plant': row['plant_name'],
            'year': event_year,
            'event_date': effective_date,
            'gvkeys': row['matched_gvkeys'].split(';'),
        })
_print(f'  First-mover events: {len(all_events)}')


# ── CAR with pre-event demeaning (matches joint_tests.py) ────────────

def compute_monthly_car(gvkey, event_month, post=POST_MONTHS):
    if gvkey not in monthly_ret:
        return None
    months = sorted(monthly_ret[gvkey].keys())
    event_idx = next((i for i, m in enumerate(months) if m >= event_month), None)
    if event_idx is None:
        return None
    pre_rets = [monthly_ret[gvkey][months[i]]
                for i in range(max(0, event_idx - PRE_MONTHS), event_idx)
                if months[i] in monthly_ret[gvkey]]
    if len(pre_rets) < 12:
        return None
    ar_list = []
    for i in range(max(0, event_idx - PRE_MONTHS), event_idx):
        m = months[i]
        if m in monthly_ret[gvkey] and m in market_ret_monthly:
            ar_list.append(monthly_ret[gvkey][m] - market_ret_monthly[m])
    pre_mean_ar = (sum(ar_list) / len(ar_list)) if ar_list else 0.0
    car = 0.0
    for offset in range(-1, post + 1):
        idx = event_idx + offset
        if 0 <= idx < len(months) and months[idx] in monthly_ret[gvkey]:
            m = months[idx]
            if m in market_ret_monthly:
                ar = monthly_ret[gvkey][m] - market_ret_monthly[m]
                car += ar - pre_mean_ar
    return car


# ── Build observation panel ──────────────────────────────────────────

_print('\nBuilding observation panel...')
obs = []
for event_id, event in enumerate(all_events):
    event_gvkeys = set(event['gvkeys'])
    event_date = event.get('event_date', '')
    event_month = (event_date[:7] if event_date and len(event_date) >= 7
                   else (f'{event["year"]}-07' if event['year'] else None))
    if not event_month:
        continue

    fm_sic4 = None
    for gk in event_gvkeys:
        fm_sic4 = get_sic4(gk)
        if fm_sic4:
            break

    for fm_gk in event_gvkeys:
        if fm_gk not in W_geo:
            continue
        neighbors = W_geo[fm_gk]
        neighbor_gks = set(neighbors.keys()) - event_gvkeys
        non_connected = [gk for gk in fundamentals
                         if gk not in event_gvkeys and gk not in neighbors]
        # Stable per-firm seed for control sampling — matches joint_tests.py
        stable_seed = int(hashlib.md5(
            str(fm_gk).encode('utf-8')).hexdigest()[:8], 16)
        random.seed(stable_seed)
        n_ctrl = min(len(non_connected),
                     max(5 * len(neighbor_gks), 20))
        ctrl_sample = (random.sample(non_connected, n_ctrl)
                       if len(non_connected) > n_ctrl else non_connected)
        candidates = list(neighbor_gks) + ctrl_sample

        for gk in candidates:
            j_sic4 = get_sic4(gk)
            same_sector = 1.0 if (fm_sic4 and j_sic4 and fm_sic4 == j_sic4) else 0.0
            car = compute_monthly_car(gk, event_month)
            if car is None:
                continue
            obs.append({
                'car': car,
                'w_geo': neighbors.get(gk, 0.0),
                'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                'same_sector': same_sector,
                'event_id': event_id,
                'gvkey': gk,
            })

_print(f'  N observations: {len(obs)}')
_print(f'  Events used: {len(set(o["event_id"] for o in obs))}')
_print(f'  Unique firms in panel: {len(set(o["gvkey"] for o in obs))}')


# ── Estimate pooled OLS and CGM SEs ───────────────────────────────────

ss_vals = set(o['same_sector'] for o in obs)
spec_vars = SPEC_VARS if len(ss_vals) > 1 else CHANNEL_VARS
res = ols(obs, 'car', spec_vars)
if res is None:
    sys.exit('OLS failed.')

event_ids = [o['event_id'] for o in obs]
firm_ids = [o['gvkey'] for o in obs]

_print('\nComputing CGM two-way clustered SEs...')
V_tw, V_event, V_firm, G_event, G_firm = twoway_cgm_vcov(
    res['X'], res['resid'], event_ids, firm_ids,
    res['inv_XtX'], res['k'], res['n']
)
_print(f'  Event clusters: {G_event}')
_print(f'  Firm clusters:  {G_firm}')

names = res['names']
betas = dict(zip(names, res['beta']))
se_event = {names[a]: math.sqrt(V_event[a][a]) if V_event[a][a] > 0 else 0.0
            for a in range(res['k'])}
se_tw = {names[a]: math.sqrt(V_tw[a][a]) if V_tw[a][a] > 0 else 0.0
         for a in range(res['k'])}
t_event = {n: (betas[n] / se_event[n] if se_event[n] > 1e-15 else 0.0)
           for n in names}
t_tw = {n: (betas[n] / se_tw[n] if se_tw[n] > 1e-15 else 0.0) for n in names}

# Channel difference (γ_geo − γ_fuel)
geo_idx = names.index('w_geo')
fuel_idx = names.index('w_fuel')
diff = betas['w_geo'] - betas['w_fuel']
var_diff_event = (V_event[geo_idx][geo_idx] + V_event[fuel_idx][fuel_idx]
                  - 2 * V_event[geo_idx][fuel_idx])
var_diff_tw = (V_tw[geo_idx][geo_idx] + V_tw[fuel_idx][fuel_idx]
               - 2 * V_tw[geo_idx][fuel_idx])
se_diff_event = math.sqrt(max(var_diff_event, 0.0))
se_diff_tw = math.sqrt(max(var_diff_tw, 0.0))
t_diff_event = diff / se_diff_event if se_diff_event > 1e-15 else 0.0
t_diff_tw = diff / se_diff_tw if se_diff_tw > 1e-15 else 0.0
p_diff_tw = p_two_sided(t_diff_tw)


# ── Write metric file ─────────────────────────────────────────────────

out_path = results_path('metrics', 'two_way_clustering.md')
lines = [
    '# Two-Way Clustered Standard Errors (Cameron, Gelbach & Miller 2011)',
    '',
    'Pooled OLS of CAR on the four channel regressors, with two-way',
    'clustered SEs on (event × firm) following the CGM (2011) formula:',
    '',
    '    V_twoway = V_event + V_firm − V_(event × firm)',
    '',
    'This complements the event-clustered SEs in `joint_tests.md`. The',
    'pooled-OLS coefficients are identical to those in `joint_tests.md`',
    'by construction; only the SEs change.',
    '',
    f'Window: [-1, +{POST_MONTHS}] months (monthly CARs, vwretd).',
    f'Spec: CAR = α + γ_geo·w^geo + γ_fuel·w^fuel + γ_reg·w^reg + γ_s·SameSector + ε.',
    f'N observations: {res["n"]:,}',
    f'Event clusters: {G_event}',
    f'Firm clusters:  {G_firm}',
    '',
    '## Headline coefficients with two-way clustered SEs',
    '',
    '| Variable | β | SE (event) | t (event) | SE (two-way CGM) | t (two-way) | p (two-way) |',
    '|---|---:|---:|---:|---:|---:|---:|',
]
for nm in names:
    p_tw = p_two_sided(t_tw[nm])
    lines.append(
        f'| {nm} | {betas[nm]:+.6f} | {se_event[nm]:.6f} | {t_event[nm]:+.3f} '
        f'| {se_tw[nm]:.6f} | {t_tw[nm]:+.3f} | {p_tw:.4f}{stars(p_tw)} |'
    )

lines += [
    '',
    '## Channel difference test (γ_geo − γ_fuel)',
    '',
    f'- Difference: {diff:+.6f}',
    f'- SE (event-clustered):  {se_diff_event:.6f}, t = {t_diff_event:+.3f}',
    f'- SE (two-way CGM):      {se_diff_tw:.6f}, t = {t_diff_tw:+.3f}, p = {p_diff_tw:.4f}',
    '',
    '## Notes',
    '',
    '- Two-way clustering on (event × firm) is the appropriate variance estimator',
    '  when within-event correlation (events as clusters of observations) and',
    '  within-firm correlation (each firm appears in many events) are both present.',
    '- The Fama-MacBeth estimator avoids this concern by construction (separate',
    '  cross-sectional regression per event); the two-way clustered pooled OLS is',
    '  reported alongside FM in Table 2 to triangulate inference under different',
    '  dependence structures.',
    '- The CGM finite-sample correction `G/(G-1) * (N-1)/(N-K)` is applied to each',
    '  of V_event, V_firm, and V_(event × firm) before combination.',
    '',
]

with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

# ── Emit panel facts JSON for downstream consumers (summary_statistics) ──
panel_firms = set(o['gvkey'] for o in obs)
panel_events = set(o['event_id'] for o in obs)
facts = {
    'n_observations': len(obs),
    'n_event_clusters': len(panel_events),
    'n_firm_clusters': len(panel_firms),
    'n_first_mover_events_total': len(all_events),
    'panel_firms_gvkeys': sorted(panel_firms),
}
facts_path = results_path('summaries', 'panel_facts.json')
os.makedirs(os.path.dirname(facts_path), exist_ok=True)
with open(facts_path, 'w', encoding='utf-8') as f:
    json.dump(facts, f, indent=2)

_print(f'\nWrote: {out_path}')
_print(f'Wrote: {facts_path}')
_print('Done.')
