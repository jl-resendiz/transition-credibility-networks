"""ESG horse race under Fama-MacBeth inference (FM joint test).

The paper text claims: "Under Fama-MacBeth inference, both ESG and fuel-mix
similarity survive a joint test." The existing esg_horse_race.py
implements this as pooled OLS with event-clustered standard errors. This
script adds the FM specification: cross-sectional regression event-by-event
within the ESG-coverage subsample, then aggregate via Newey-West.

The joint test is Wald(H0: gamma_ESG = 0 AND gamma_fuel = 0) computed under
the FM aggregation, treating the two coefficient time series as bivariate
with covariance estimated from the event-level beta panel.

Inputs:  data/raw/refinitiv/refinitiv_esg.csv (env_score)
         (everything else as in robust_inference.py)

Outputs: results/metrics/esg_fm_joint.md
         results/summaries/esg_fm_joint.csv
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


# ESG scores
esg_scores = {}
esg_path = raw_path('refinitiv', 'refinitiv_esg.csv')
if os.path.exists(esg_path):
    with open(esg_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row.get('gvkey', '').strip()
            env = row.get('env_score', '').strip()
            if gk and env:
                try:
                    esg_scores[gk] = float(env) / 100.0
                except ValueError:
                    pass
_print(f'ESG scores loaded: {len(esg_scores)} firms')

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


# ── CAR ─────────────────────────────────────────────────────────────

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


# ── Build per-event datasets restricted to ESG firms ────────────────

MIN_OBS_PER_EVENT = 10
event_datasets = {}
for event_id, event in enumerate(all_events):
    em = event['event_month']
    event_gvkeys = set(event['gvkeys'])
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
            if gk not in esg_scores:
                continue
            car = compute_car(gk, em)
            if car is None:
                continue
            obs.append({
                'car': car,
                'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                'w_geo': neighbors.get(gk, 0.0),
                'esg_score': esg_scores[gk],
            })
    if len(obs) >= MIN_OBS_PER_EVENT:
        event_datasets[event_id] = obs

n_valid = len(event_datasets)
total_obs = sum(len(v) for v in event_datasets.values())
_print(f'Events with >= {MIN_OBS_PER_EVENT} ESG-covered firms: {n_valid}')
_print(f'Total firm-event obs in ESG subsample: {total_obs}')


# ── Run FM regressions: 3 specifications ────────────────────────────

specs = {
    'ESG only': ['esg_score'],
    'Spatial only (fuel + geo)': ['w_fuel', 'w_geo'],
    'Both (ESG + spatial)': ['esg_score', 'w_fuel', 'w_geo'],
}

per_spec_betas = {name: defaultdict(list) for name in specs}
for event_id in sorted(event_datasets.keys()):
    obs = event_datasets[event_id]
    if len(obs) < MIN_OBS_PER_EVENT:
        continue
    for spec_name, vars_list in specs.items():
        y = [o['car'] for o in obs]
        X = [[1.0] + [o[v] for v in vars_list] for o in obs]
        result = ols_simple(y, X)
        if result is None:
            continue
        names = ['intercept'] + vars_list
        for i, name in enumerate(names):
            per_spec_betas[spec_name][name].append(result['beta'][i])


# ── Aggregate FM with NW SEs ────────────────────────────────────────

def newey_west_se(series, lag=4):
    series = [x for x in series if not (isinstance(x, float) and math.isnan(x))]
    n = len(series)
    if n < 2:
        return float('nan'), float('nan'), 0
    mean = sum(series) / n
    dev = [x - mean for x in series]
    gamma0 = sum(d * d for d in dev) / n
    var = gamma0
    for L in range(1, min(lag, n - 1) + 1):
        weight = 1.0 - L / (lag + 1)
        cov_L = sum(dev[t] * dev[t - L] for t in range(L, n)) / n
        var += 2.0 * weight * cov_L
    if var <= 0:
        return mean, float('nan'), n
    return mean, math.sqrt(var / n), n


def fm_covariance(s1, s2, lag=4):
    """Newey-West covariance of two FM time series."""
    n = min(len(s1), len(s2))
    s1, s2 = s1[:n], s2[:n]
    if n < 2:
        return float('nan')
    m1, m2 = sum(s1) / n, sum(s2) / n
    d1 = [x - m1 for x in s1]
    d2 = [x - m2 for x in s2]
    gamma0 = sum(d1[t] * d2[t] for t in range(n)) / n
    cov = gamma0
    for L in range(1, min(lag, n - 1) + 1):
        weight = 1.0 - L / (lag + 1)
        ct1 = sum(d1[t] * d2[t - L] for t in range(L, n)) / n
        ct2 = sum(d2[t] * d1[t - L] for t in range(L, n)) / n
        cov += weight * (ct1 + ct2)
    return cov / n


# Collect FM results per spec
fm_results = {}
for spec_name, betas_dict in per_spec_betas.items():
    spec_results = {}
    for var, series in betas_dict.items():
        mean, se, n = newey_west_se(series, lag=4)
        spec_results[var] = {'mean': mean, 'se': se, 'n': n,
                             't': mean / se if (not math.isnan(se) and se > 0) else float('nan')}
    fm_results[spec_name] = spec_results


# ── Joint Wald test for "Both" specification ────────────────────────

both_betas = per_spec_betas['Both (ESG + spatial)']
esg_series = both_betas['esg_score']
fuel_series = both_betas['w_fuel']
n_common = min(len(esg_series), len(fuel_series))
mean_esg = sum(esg_series[:n_common]) / n_common if n_common else float('nan')
mean_fuel = sum(fuel_series[:n_common]) / n_common if n_common else float('nan')
var_esg = newey_west_se(esg_series[:n_common], lag=4)[1] ** 2 if n_common else float('nan')
var_fuel = newey_west_se(fuel_series[:n_common], lag=4)[1] ** 2 if n_common else float('nan')
cov_esg_fuel = fm_covariance(esg_series[:n_common], fuel_series[:n_common], lag=4)

# 2x2 covariance matrix of (mean_esg, mean_fuel)
V = [[var_esg, cov_esg_fuel],
     [cov_esg_fuel, var_fuel]]
V_inv = invert_matrix(V)
if V_inv is not None:
    diff = [mean_esg, mean_fuel]
    wald = sum(diff[i] * V_inv[i][j] * diff[j] for i in range(2) for j in range(2))
else:
    wald = float('nan')

# p-value: chi-squared with 2 df. Approximation via gamma.
# Using regularized: P(X>w) for chi-sq_2 = exp(-w/2)
joint_p = math.exp(-wald / 2) if not math.isnan(wald) and wald >= 0 else float('nan')


# ── Output ──────────────────────────────────────────────────────────

out_csv = results_path('summaries', 'esg_fm_joint.csv')
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['spec', 'variable', 'mean', 'se_nw_lag4', 't', 'n_events'])
    for spec_name, results_dict in fm_results.items():
        for var, r in results_dict.items():
            w.writerow([
                spec_name, var,
                f'{r["mean"]:+.6f}',
                f'{r["se"]:.6f}' if not math.isnan(r['se']) else 'NA',
                f'{r["t"]:+.4f}' if not math.isnan(r['t']) else 'NA',
                r['n'],
            ])
_print(f'Wrote {out_csv}')

out_md = results_path('metrics', 'esg_fm_joint.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# ESG Horse Race under Fama-MacBeth Inference\n\n')
    f.write(f'Sample restricted to firms with Refinitiv Environmental Score '
            f'({total_obs} firm-event observations across {n_valid} events). '
            f'For each event, runs cross-sectional OLS within the ESG-covered '
            f'subsample; aggregates with Newey-West (lag 4) HAC standard errors.\n\n')

    for spec_name, results_dict in fm_results.items():
        f.write(f'## Spec: {spec_name}\n\n')
        f.write('| Variable | Mean | SE (NW lag=4) | t | N events |\n')
        f.write('|---|---|---|---|---|\n')
        for var, r in results_dict.items():
            f.write(f'| {var} | {r["mean"]:+.4f} | '
                    f'{r["se"]:.4f} | {r["t"]:+.2f} | {r["n"]} |\n')
        f.write('\n')

    f.write('## Joint test (Both spec): H0: gamma_ESG = 0 AND gamma_fuel = 0\n\n')
    f.write(f'- Estimated mean ESG coefficient: {mean_esg:+.4f}\n')
    f.write(f'- Estimated mean fuel coefficient: {mean_fuel:+.4f}\n')
    f.write(f'- Wald statistic: {wald:.3f} (chi-squared with 2 df)\n')
    f.write(f'- Approximate p-value: {joint_p:.4f}\n\n')
    if not math.isnan(joint_p) and joint_p < 0.05:
        f.write('**Reject the joint null** at the 5% level: at least one of the '
                'two coefficients is non-zero in the FM cross-section.\n\n')
        # Check whether each survives marginally
        esg_t = fm_results['Both (ESG + spatial)']['esg_score']['t']
        fuel_t = fm_results['Both (ESG + spatial)']['w_fuel']['t']
        if abs(esg_t) >= 1.96 and abs(fuel_t) >= 1.96:
            f.write('Both individually survive marginal significance at the 5% level. '
                    'The paper claim that "both ESG and fuel-mix similarity survive '
                    'a joint test" under FM inference is supported.\n')
        elif abs(esg_t) >= 1.96 and abs(fuel_t) < 1.96:
            f.write('ESG individually survives (|t| >= 1.96); fuel does not. The '
                    'paper claim that both survive a joint FM test is partially '
                    'supported (joint significance holds, but only ESG carries '
                    'the marginal weight).\n')
        elif abs(fuel_t) >= 1.96 and abs(esg_t) < 1.96:
            f.write('Fuel individually survives; ESG does not.\n')
        else:
            f.write('Neither coefficient survives individual marginal significance, '
                    'though the joint test rejects.\n')
    else:
        f.write('Joint null is not rejected at 5%.\n')

_print(f'Wrote {out_md}')
_print('\nDone.')
