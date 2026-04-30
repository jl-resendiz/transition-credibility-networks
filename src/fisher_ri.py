"""Fisher randomization-inference (RI) p-value for the FM fuel coefficient.

The Honest DID placebo distribution in `honest_did.md` reports
the maximum |beta_fuel| across pre-event windows alongside the headline.
A more rigorous identification defense is randomization inference: under
the sharp null of no treatment effect, randomly permute the firm-level
fuel-similarity weights within each event and recompute the FM coefficient.
The fraction of permutations more extreme than the observed headline is
the Fisher RI p-value.

This is the cleanest possible defense against the "1.72x placebo gap"
critique: it does not assume a parametric model and is immune to spatial
or serial dependence in residuals.

Inputs:  monthly_returns, weight matrices, events (same as headline)
Outputs: results/metrics/fisher_ri.md
         results/summaries/fisher_ri_distribution.csv
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


# ── CAR computation ─────────────────────────────────────────────────

PRE_MONTHS = 24


def compute_car(gvkey, em):
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
    for offset in range(-1, 4):
        idx = event_idx + offset
        if 0 <= idx < len(months) and months[idx] in monthly_ret[gvkey]:
            m = months[idx]
            if m in market_ret:
                car += (monthly_ret[gvkey][m] - market_ret[m]) - pre_mean
                n += 1
    if n < 3:
        return None
    return car


# ── Build per-event observations (CARs + weights) ───────────────────

MIN_OBS = 20
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']


def get_sic4(gvkey):
    f = fundamentals.get(gvkey)
    return f['sic'][:4] if (f and f.get('sic')) else None


_print('Building per-event observation panels...')
event_data = {}
for event_id, event in enumerate(all_events):
    em = event['event_month']
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
            j_sic4 = get_sic4(gk)
            car = compute_car(gk, em)
            if car is None:
                continue
            obs.append({
                'car': car,
                'w_geo': neighbors.get(gk, 0.0),
                'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                'same_sector': 1.0 if (fm_sic4 and j_sic4 and fm_sic4 == j_sic4) else 0.0,
            })
    if len(obs) >= MIN_OBS:
        event_data[event_id] = obs


# ── Compute observed FM beta_fuel ───────────────────────────────────

def fm_beta_fuel_from(event_data_dict, perm_within_event=False, seed=0):
    """Compute FM beta_fuel. If perm_within_event, randomly permute w_fuel
    within each event before estimation."""
    rng = random.Random(seed)
    fuel_betas = []
    for event_id, obs_list in event_data_dict.items():
        if perm_within_event:
            w_fuels = [o['w_fuel'] for o in obs_list]
            permuted = w_fuels[:]
            rng.shuffle(permuted)
            obs_list = [{**o, 'w_fuel': permuted[i]} for i, o in enumerate(obs_list)]
        ss_vals = set(o['same_sector'] for o in obs_list)
        use_vars = SPEC_VARS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']
        y = [o['car'] for o in obs_list]
        X = [[1.0] + [o[v] for v in use_vars] for o in obs_list]
        result = ols_simple(y, X)
        if result is None:
            continue
        bd = dict(zip(['intercept'] + use_vars, result['beta']))
        fuel_betas.append(bd.get('w_fuel', float('nan')))
    if not fuel_betas:
        return float('nan')
    return sum(fuel_betas) / len(fuel_betas)


_print('Computing observed FM beta_fuel...')
obs_beta = fm_beta_fuel_from(event_data, perm_within_event=False)
_print(f'  Observed FM beta_fuel = {obs_beta:+.4f}')


# ── Fisher RI: permute w_fuel within event and recompute ─────────────

N_PERM = 999
_print(f'\nRunning {N_PERM} within-event permutations...')

perm_betas = []
report_every = max(1, N_PERM // 10)
for p in range(1, N_PERM + 1):
    pb = fm_beta_fuel_from(event_data, perm_within_event=True, seed=p)
    perm_betas.append(pb)
    if p % report_every == 0:
        _print(f'  Permutation {p}/{N_PERM} (last beta = {pb:+.4f})')


# ── Compute one-sided RI p-value ────────────────────────────────────

# H0: w_fuel has no effect on CARs (sharp null)
# Alternative: beta_fuel < observed (one-sided in direction of observed sign)
n_more_extreme = sum(1 for pb in perm_betas
                     if not math.isnan(pb) and pb <= obs_beta)
ri_p_one_sided = (n_more_extreme + 1) / (len(perm_betas) + 1)

# Two-sided
n_two_sided = sum(1 for pb in perm_betas
                  if not math.isnan(pb) and abs(pb) >= abs(obs_beta))
ri_p_two_sided = (n_two_sided + 1) / (len(perm_betas) + 1)


# ── Output ──────────────────────────────────────────────────────────

# CSV: distribution of permuted betas
out_csv = results_path('summaries', 'fisher_ri_distribution.csv')
os.makedirs(os.path.dirname(out_csv), exist_ok=True)
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['permutation_id', 'fm_beta_fuel'])
    w.writerow([0, f'{obs_beta:+.6f}'])  # observed
    for i, pb in enumerate(perm_betas, start=1):
        w.writerow([i, f'{pb:+.6f}' if not math.isnan(pb) else 'NA'])
_print(f'\nWrote {out_csv}')

# Markdown
out_md = results_path('metrics', 'fisher_ri.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Fisher Randomization Inference for the FM Fuel Coefficient\n\n')
    f.write(f'Within-event permutation of $w^{{\\mathrm{{fuel}}}}_i$ values, '
            f'recomputing the Fama-MacBeth average. {N_PERM} permutations. '
            'Tests the sharp null that fuel-mix similarity carries no '
            'cross-sectional information about CARs at retirement events.\n\n')

    f.write('## Distribution\n\n')
    valid = [pb for pb in perm_betas if not math.isnan(pb)]
    if valid:
        valid_sorted = sorted(valid)
        n_v = len(valid)
        f.write(f'- Permutation mean: {sum(valid)/n_v:+.4f}\n')
        f.write(f'- Permutation 1st percentile: {valid_sorted[max(0, n_v//100)]:+.4f}\n')
        f.write(f'- Permutation 5th percentile: {valid_sorted[max(0, n_v//20)]:+.4f}\n')
        f.write(f'- Permutation 50th percentile (median): {valid_sorted[n_v//2]:+.4f}\n')
        f.write(f'- Permutation 95th percentile: {valid_sorted[min(n_v-1, 19*n_v//20)]:+.4f}\n')
        f.write(f'- Permutation 99th percentile: {valid_sorted[min(n_v-1, 99*n_v//100)]:+.4f}\n')
        f.write(f'- Range: [{valid_sorted[0]:+.4f}, {valid_sorted[-1]:+.4f}]\n\n')

    f.write('## Test\n\n')
    f.write(f'- **Observed FM beta_fuel: {obs_beta:+.4f}**\n')
    f.write(f'- One-sided RI p-value (P[perm <= observed]): {ri_p_one_sided:.4f}\n')
    f.write(f'- Two-sided RI p-value (P[|perm| >= |observed|]): {ri_p_two_sided:.4f}\n\n')

    if ri_p_one_sided < 0.001:
        f.write(f'**Reject** the sharp null at $p < 0.001$. The observed '
                f'FM coefficient is more extreme than {N_PERM} of {N_PERM} '
                f'within-event permutations.\n')
    elif ri_p_one_sided < 0.01:
        f.write(f'**Reject** the sharp null at $p < 0.01$. RI p = {ri_p_one_sided:.4f}.\n')
    elif ri_p_one_sided < 0.05:
        f.write(f'**Reject** the sharp null at $p < 0.05$. RI p = {ri_p_one_sided:.4f}.\n')
    else:
        f.write(f'**Fail to reject** the sharp null at conventional levels. '
                f'RI p = {ri_p_one_sided:.4f}.\n')

    f.write('\n## Interpretation\n\n')
    f.write('Randomization inference tests the sharp null that the firm-level '
            '$w^{\\mathrm{fuel}}$ values carry no information about CARs '
            'within an event. Under this null, randomly relabelling the firms '
            'should produce a distribution of FM coefficients centred at zero. '
            'The observed coefficient, if it falls in the tail of the '
            'permutation distribution, is unlikely to have arisen by chance. '
            'This test is robust to spatial dependence, serial correlation, '
            'and any parametric assumptions about residual structure.\n')

_print(f'Wrote {out_md}')
_print(f'\nFisher RI p-value (one-sided): {ri_p_one_sided:.4f}')
_print(f'Fisher RI p-value (two-sided):   {ri_p_two_sided:.4f}')
_print('Done.')
