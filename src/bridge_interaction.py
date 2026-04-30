"""Empirical test of the bridge equation in Section 2.5.

The decomposition in Section 2.5 of the manuscript predicts that the firm-
specific load -A*phi*alpha_i, when projected onto a single scalar
gamma_fuel via cross-sectional OLS, leaves a refutable interaction:

    CAR_i = a + gamma_fuel * w_fuel_i
            + gamma_het * w_fuel_i * (alpha_i^pre2014 - alpha_bar)
            + epsilon_i

with predicted gamma_het < 0. Pre-2014 coal share serves as a pre-determined
proxy for alpha_i, avoiding endogeneity to the contemporaneous retirement
event. Single-factor CARs are used to match the headline reference; a
multi-factor robustness version is straightforward and noted in the docs.

Inputs:  data/derived/returns/monthly_returns.csv
         data/derived/fundamentals/firm_alpha_panel.csv (2010-)
         data/derived/networks/weight_matrix_W_*.csv
         data/derived/events/coal_retirement_events.csv

Outputs: results/metrics/bridge_interaction.md
         results/summaries/bridge_interaction_betas.csv
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


def load_market_monthly(path):
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


market_ret = load_market_monthly(raw_path('factors', 'F-F_Research_Data_Factors.csv'))

# Pre-2014 alpha: mean coal share over 2010-2013 for each firm
_print('Loading pre-2014 alpha (coal share 2010-2013 mean)...')
firm_pre_alpha = {}
alpha_year_records = defaultdict(list)
with open(derived_path('fundamentals', 'firm_alpha_panel.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        try:
            yr = int(row['year'])
            cs = float(row['coal_share'])
        except (ValueError, KeyError):
            continue
        if 2010 <= yr <= 2013:
            alpha_year_records[gk].append(cs)
for gk, cs_list in alpha_year_records.items():
    if cs_list:
        firm_pre_alpha[gk] = sum(cs_list) / len(cs_list)
_print(f'  Firms with pre-2014 coal share: {len(firm_pre_alpha)}')

# Cross-sectional mean coal share (used to demean the interaction)
alpha_bar = (sum(firm_pre_alpha.values()) / len(firm_pre_alpha)
             if firm_pre_alpha else 0.0)
_print(f'  Cross-sectional mean alpha (pre-2014): {alpha_bar:.4f}')

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


# ── Single-factor CAR (matches headline) ─────────────────────────────

PRE_MONTHS = 24
POST_MONTHS = 3


def compute_car_single(gvkey, event_month_str):
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


# ── Build per-event datasets with pre-2014 alpha attached ───────────

_print('\nBuilding per-event datasets with pre-2014 alpha...')
MIN_OBS_PER_EVENT = 20
event_datasets = {}

for event_id, event in enumerate(all_events):
    event_gvkeys = set(event['gvkeys'])
    em = event['event_month']
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
        candidate_firms = list(neighbor_gks) + ctrl_sample
        for gk in candidate_firms:
            j_sic4 = get_sic4(gk)
            same_sector = 1.0 if (fm_sic4 and j_sic4 and fm_sic4 == j_sic4) else 0.0
            car = compute_car_single(gk, em)
            if car is None:
                continue
            alpha_pre = firm_pre_alpha.get(gk)
            if alpha_pre is None:
                continue
            obs.append({
                'car': car,
                'w_geo': neighbors.get(gk, 0.0),
                'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                'same_sector': same_sector,
                'alpha_dev': alpha_pre - alpha_bar,
            })
    if len(obs) >= MIN_OBS_PER_EVENT:
        event_datasets[event_id] = obs

n_valid = len(event_datasets)
total_obs = sum(len(v) for v in event_datasets.values())
_print(f'  Valid events: {n_valid}, total firm-event obs: {total_obs}')


# ── Augmented FM regression ─────────────────────────────────────────

def newey_west_se(series, lag=4):
    series = [x for x in series if not (isinstance(x, float) and math.isnan(x))]
    n = len(series)
    if n < 2:
        return float('nan')
    mean = sum(series) / n
    dev = [x - mean for x in series]
    gamma0 = sum(d * d for d in dev) / n
    var_nw = gamma0
    for L in range(1, min(lag, n - 1) + 1):
        weight = 1.0 - L / (lag + 1)
        cov_L = sum(dev[t] * dev[t - L] for t in range(L, n)) / n
        var_nw += 2.0 * weight * cov_L
    if var_nw <= 0:
        return float('nan')
    return math.sqrt(var_nw / n)


SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector', 'w_fuel_x_alpha']

event_betas = defaultdict(list)
event_ns = []

for event_id in sorted(event_datasets.keys()):
    obs = event_datasets[event_id]
    # Build w_fuel * (alpha - alpha_bar) interaction
    for o in obs:
        o['w_fuel_x_alpha'] = o['w_fuel'] * o['alpha_dev']
    ss_vals = set(o['same_sector'] for o in obs)
    use_vars = SPEC_VARS if len(ss_vals) > 1 else [
        'w_geo', 'w_fuel', 'w_reg', 'w_fuel_x_alpha']
    y = [o['car'] for o in obs]
    X = [[1.0] + [o[v] for v in use_vars] for o in obs]
    result = ols_simple(y, X)
    if result is None:
        continue
    names = ['intercept'] + use_vars
    for i, name in enumerate(names):
        event_betas[name].append(result['beta'][i])
    for v in SPEC_VARS:
        if v not in use_vars:
            event_betas[v].append(float('nan'))
    event_ns.append(result['n'])

T = len(event_ns)
_print(f'  Events with valid augmented FM regressions: {T}')


# Aggregate
fm_results = {}
for name in ['w_geo', 'w_fuel', 'w_reg', 'w_fuel_x_alpha']:
    series = [b for b in event_betas[name] if not math.isnan(b)]
    if not series:
        continue
    n = len(series)
    mean = sum(series) / n
    se_nw = newey_west_se(series, lag=4)
    fm_results[name] = {
        'mean': mean,
        'se_nw': se_nw,
        't': mean / se_nw if (not math.isnan(se_nw) and se_nw > 0) else float('nan'),
        'n_events': n,
    }


# ── Output ──────────────────────────────────────────────────────────

out_csv = results_path('summaries', 'bridge_interaction_betas.csv')
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['variable', 'mean', 'se_nw_lag4', 't', 'n_events',
                'predicted_sign'])
    for v, sign in [('w_fuel', 'negative'),
                    ('w_geo', 'attenuated to zero'),
                    ('w_reg', 'positive'),
                    ('w_fuel_x_alpha', 'negative (theory)')]:
        if v in fm_results:
            r = fm_results[v]
            w.writerow([v, f'{r["mean"]:+.6f}',
                        f'{r["se_nw"]:.6f}' if not math.isnan(r['se_nw']) else 'NA',
                        f'{r["t"]:+.4f}' if not math.isnan(r['t']) else 'NA',
                        r['n_events'], sign])
_print(f'\nWrote {out_csv}')

out_md = results_path('metrics', 'bridge_interaction.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Bridge Interaction: Empirical Test of Section 2.5\n\n')
    f.write('Augmented Fama-MacBeth specification:\n\n')
    f.write('$$\\mathrm{CAR}_i = a + \\gamma_{\\mathrm{fuel}} w^{\\mathrm{fuel}}_i '
            '+ \\gamma_{\\mathrm{het}} w^{\\mathrm{fuel}}_i (\\alpha_i^{pre-2014} - \\bar\\alpha) '
            '+ \\gamma_{\\mathrm{geo}} w^{\\mathrm{geo}}_i + \\gamma_{\\mathrm{reg}} w^{\\mathrm{reg}}_i '
            '+ \\gamma_{ss} \\mathbb{1}_{ss} + \\varepsilon_i$$\n\n')
    f.write('with $\\alpha_i^{pre-2014}$ = mean firm coal share 2010-2013, '
            f'$\\bar\\alpha = {alpha_bar:.4f}$ (cross-sectional mean of pre-2014 alpha). '
            'Single-factor (market-adjusted) CARs over [-1, +3] are used to match the '
            'headline reference table.\n\n')

    f.write('## Augmented FM coefficients\n\n')
    f.write('| Variable | Mean | SE (NW lag=4) | t-stat | N events | Predicted sign |\n')
    f.write('|---|---|---|---|---|---|\n')
    for v, lbl, sign in [
        ('w_fuel', '$\\gamma_{\\mathrm{fuel}}$', 'negative'),
        ('w_fuel_x_alpha', '$\\gamma_{\\mathrm{het}}$ ($w^{\\mathrm{fuel}} \\times (\\alpha-\\bar\\alpha)$)', 'negative (theory)'),
        ('w_geo', '$\\gamma_{\\mathrm{geo}}$', 'attenuated to zero'),
        ('w_reg', '$\\gamma_{\\mathrm{reg}}$', 'positive'),
    ]:
        if v in fm_results:
            r = fm_results[v]
            f.write(f'| {lbl} | {r["mean"]:+.4f} | '
                    f'{r["se_nw"]:.4f} | {r["t"]:+.2f} | {r["n_events"]} | {sign} |\n')

    f.write('\n## Interpretation\n\n')
    if 'w_fuel_x_alpha' in fm_results:
        r = fm_results['w_fuel_x_alpha']
        if not math.isnan(r['t']):
            t_abs = abs(r['t'])
            sign_str = 'negative' if r['mean'] < 0 else 'positive'
            f.write(f'The interaction coefficient $\\gamma_{{\\mathrm{{het}}}}$ is '
                    f'**{sign_str}** ($\\hat\\gamma_{{\\mathrm{{het}}}} = {r["mean"]:+.4f}$, '
                    f'$|t| = {t_abs:.2f}$). ')
            if r['mean'] < 0 and t_abs >= 1.96:
                f.write('The sign is consistent with the theory prediction in '
                        'Section 2.5, and the coefficient is statistically significant '
                        'at the 5% level. This corroborates the firm-specific '
                        'heterogeneity implied by the load $-A\\phi\\alpha_i$.\n')
            elif r['mean'] < 0 and t_abs >= 1.65:
                f.write('The sign is consistent with theory; significance at the 10% '
                        'level provides directional support without precise '
                        'identification of the magnitude.\n')
            elif r['mean'] < 0:
                f.write('The sign is consistent with theory but not statistically '
                        'distinguishable from zero. Power is limited by the within-'
                        'event compression of $\\alpha_i$ variation.\n')
            else:
                f.write('The sign is opposite to the theory prediction. '
                        'This warrants a reassessment of the linearity assumption '
                        'in Assumption 2 of the model.\n')

    f.write('\n## Notes\n\n')
    f.write(f'- Events with pre-2014 alpha coverage and >= {MIN_OBS_PER_EVENT} firms: {T}\n')
    f.write(f'- Total firm-event observations: {total_obs}\n')
    f.write('- Pre-2014 alpha used as a pre-determined proxy avoids endogeneity to '
            'the contemporaneous retirement event.\n')
    f.write('- A multi-factor robustness version (using FF3 + utility CARs) is '
            'a natural extension; see multifactor_inference.py for the '
            'multi-factor CAR construction.\n')

_print(f'Wrote {out_md}')
_print('\nDone.')
