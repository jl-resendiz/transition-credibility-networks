"""Multi-factor abnormal returns and cross-sectional inference.

The headline inference in robust_inference.py uses single-factor
(market-adjusted) abnormal returns. The original quant referee critique:
single-factor adjustment is inadequate for 703 utilities across 80 countries,
because (i) utility-industry beta heterogeneity is substantial and (ii) the
fuel coefficient may load on the utility-industry factor.

This script re-estimates the headline regression using a 4-factor model:

    (r_it - rf_t) = alpha_i + b_M (Mkt-RF)_t + b_S SMB_t + b_V HML_t
                    + b_U (UTL_excess)_t + epsilon_it

where UTL is the equal-weighted utility-sector portfolio constructed from
the paper's own sample (all firms with monthly returns at month t). Betas
are estimated on a 24-month pre-event estimation window for each (firm,
event) pair. The abnormal return is the within-window prediction error
under these firm-specific betas.

Currency factor: omitted. Constructing a USD trade-weighted index requires
external data (FRED/DXY) not in the repo. The market factor (FF global
Mkt-RF) absorbs much aggregate dollar exposure for non-US firms; remaining
currency idiosyncrasies are absorbed into the firm-level alpha. This is
documented as a limitation in the manuscript.

Inputs:  data/derived/returns/monthly_returns.csv
         data/raw/factors/F-F_Research_Data_Factors.csv
         data/derived/networks/weight_matrix_W_*.csv
         data/derived/events/coal_retirement_events.csv

Outputs: results/metrics/multifactor_inference.md
         results/summaries/event_level_betas_mf.csv
         results/summaries/multifactor_comparison.csv
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
    ss_tot = sum((yi - sum(y) / n) ** 2 for yi in y)
    ss_res = sum(r ** 2 for r in resid)
    return {'beta': beta, 'resid': resid, 'n': n,
            'r2': 1 - ss_res / ss_tot if ss_tot > 1e-15 else 0.0}


# ── Load monthly returns ────────────────────────────────────────────

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


# ── Load FF3 factors (Mkt-RF, SMB, HML, RF) at monthly frequency ────

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
if not ff3:
    raise RuntimeError('Missing FF3 factors.')
_print(f'  FF3 months: {len(ff3)}')


# ── Construct utility industry factor (equal-weighted sample mean) ──

_print('Constructing utility industry factor from sample...')
month_returns = defaultdict(list)
for gk, dates in monthly_ret.items():
    for m, r in dates.items():
        month_returns[m].append(r)

util_factor = {}  # month -> mean return across sample (in decimal)
for m, rs in month_returns.items():
    if len(rs) >= 30:  # require at least 30 firms for stable sector mean
        util_factor[m] = sum(rs) / len(rs)

# Combine with FF3 to produce a single factors dict per month
factors_panel = {}
for m, ff in ff3.items():
    if m in util_factor:
        factors_panel[m] = {
            **ff,
            'utl_excess': util_factor[m] - ff['rf'],
        }
_print(f'  Months with all factors: {len(factors_panel)}')


# ── Load weight matrices, fundamentals, events ──────────────────────

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


# ── Multi-factor CAR ─────────────────────────────────────────────────

PRE_MONTHS = 24
PRE_MIN = 12
POST_MONTHS = 3


def compute_car_multifactor(gvkey, event_month_str):
    """4-factor abnormal return CAR over [-1, +POST_MONTHS]."""
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

    # Event-window CAR using firm-specific multi-factor model
    car = 0.0
    n_months = 0
    for offset in range(-1, POST_MONTHS + 1):
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
    if n_months < (POST_MONTHS + 1) // 2 + 1:
        return None
    return car


# ── Build per-event candidate firm sets and compute multi-factor CARs ──

_print('\nComputing multi-factor CARs per (event, firm)...')
MIN_OBS_PER_EVENT = 20
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']

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
            car_mf = compute_car_multifactor(gk, em)
            if car_mf is None:
                continue
            obs.append({
                'car': car_mf,
                'w_geo': neighbors.get(gk, 0.0),
                'w_fuel': W_fuel.get(fm_gk, {}).get(gk, 0.0),
                'w_reg': W_reg.get(fm_gk, {}).get(gk, 0.0),
                'same_sector': same_sector,
                'gvkey': gk,
            })
    if len(obs) >= MIN_OBS_PER_EVENT:
        event_datasets[event_id] = obs

n_valid = len(event_datasets)
total_obs = sum(len(v) for v in event_datasets.values())
_print(f'  Valid events: {n_valid}, total obs: {total_obs}')


# ── Fama-MacBeth cross-sectional regression ──────────────────────────

_print('\nRunning FM cross-sectional regressions...')

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


event_betas = defaultdict(list)
event_r2s, event_ns, event_ids_used = [], [], []
for event_id in sorted(event_datasets.keys()):
    obs = event_datasets[event_id]
    ss_vals = set(o['same_sector'] for o in obs)
    use_vars = SPEC_VARS if len(ss_vals) > 1 else ['w_geo', 'w_fuel', 'w_reg']
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
    event_r2s.append(result['r2'])
    event_ns.append(result['n'])
    event_ids_used.append(event_id)

T = len(event_ids_used)
_print(f'  Events with valid FM regressions: {T}')


# Aggregate FM betas with NW lag=4
fm_results = {}
for name in ['w_geo', 'w_fuel', 'w_reg']:
    series = [b for b in event_betas[name] if not math.isnan(b)]
    if not series:
        continue
    n = len(series)
    mean = sum(series) / n
    se_nw = newey_west_se(series, lag=4)
    fm_results[name] = {
        'mean': mean, 'se_nw': se_nw,
        't': mean / se_nw if (not math.isnan(se_nw) and se_nw > 0) else float('nan'),
        'n_events': n,
    }

# Difference geo - fuel
diff_series = [event_betas['w_geo'][i] - event_betas['w_fuel'][i]
               for i in range(T)
               if not math.isnan(event_betas['w_geo'][i])
               and not math.isnan(event_betas['w_fuel'][i])]
if diff_series:
    n = len(diff_series)
    mean = sum(diff_series) / n
    se_nw = newey_west_se(diff_series, lag=4)
    fm_results['diff_geo_fuel'] = {
        'mean': mean, 'se_nw': se_nw,
        't': mean / se_nw if (not math.isnan(se_nw) and se_nw > 0) else float('nan'),
        'n_events': n,
    }


# ── Output ───────────────────────────────────────────────────────────

out_csv_betas = results_path('summaries', 'event_level_betas_mf.csv')
os.makedirs(os.path.dirname(out_csv_betas), exist_ok=True)
with open(out_csv_betas, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['event_id', 'event_n', 'event_r2', 'beta_fuel', 'beta_geo',
                'beta_reg', 'beta_same_sector'])
    for i in range(T):
        w.writerow([
            event_ids_used[i], event_ns[i], f'{event_r2s[i]:.6f}',
            f'{event_betas["w_fuel"][i]:.6f}',
            f'{event_betas["w_geo"][i]:.6f}',
            f'{event_betas["w_reg"][i]:.6f}',
            f'{event_betas["same_sector"][i]:.6f}'
            if not math.isnan(event_betas['same_sector'][i]) else 'NA',
        ])

# Comparison CSV (single-factor headline numbers from previous run)
out_csv_cmp = results_path('summaries', 'multifactor_comparison.csv')
with open(out_csv_cmp, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['channel', 'spec', 'mean', 'se_nw_lag4', 't', 'n_events'])
    # Single-factor numbers from existing output (reference; not recomputed here)
    sf_ref = {
        'w_fuel': (-4.7656, 0.6508, -7.32, 117),
        'w_geo': (-0.5427, 0.3090, -1.76, 117),
        'w_reg': (2.6975, 0.9518, 2.83, 117),
        'diff_geo_fuel': (4.2229, 0.7076, 5.97, 117),
    }
    for ch, lbl in [('w_fuel', 'fuel'), ('w_geo', 'geo'),
                    ('w_reg', 'reg'), ('diff_geo_fuel', 'diff_geo_fuel')]:
        if ch in sf_ref:
            m, s, t, n = sf_ref[ch]
            w.writerow([lbl, 'single_factor',
                        f'{m:+.6f}', f'{s:.6f}', f'{t:+.4f}', n])
        if ch in fm_results:
            r = fm_results[ch]
            w.writerow([lbl, 'multi_factor',
                        f'{r["mean"]:+.6f}',
                        f'{r["se_nw"]:.6f}' if not math.isnan(r['se_nw']) else 'NA',
                        f'{r["t"]:+.4f}' if not math.isnan(r['t']) else 'NA',
                        r['n_events']])

# Markdown
out_md = results_path('metrics', 'multifactor_inference.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Multi-Factor Abnormal Returns: FM Cross-Sectional Inference\n\n')
    f.write('Replaces single-factor (market-adjusted) abnormal returns with a 4-factor '
            'model (FF3 Mkt-RF, SMB, HML + sample-constructed utility industry excess '
            'return) estimated firm-by-firm on a 24-month pre-event window. CAR is '
            'the within-window prediction error over the headline window [-1, +3].\n\n')

    f.write('## Headline comparison (FM + Newey-West, lag=4)\n\n')
    f.write('| Channel | Single-factor | Multi-factor | Shrinkage |\n')
    f.write('|---|---|---|---|\n')
    sf_ref = {
        'w_fuel': (-4.7656, 0.6508, -7.32),
        'w_geo': (-0.5427, 0.3090, -1.76),
        'w_reg': (2.6975, 0.9518, 2.83),
        'diff_geo_fuel': (4.2229, 0.7076, 5.97),
    }
    labels = {
        'w_fuel': '$\\gamma_{\\text{fuel}}$',
        'w_geo': '$\\gamma_{\\text{geo}}$',
        'w_reg': '$\\gamma_{\\text{reg}}$',
        'diff_geo_fuel': '$\\gamma_{\\text{geo}} - \\gamma_{\\text{fuel}}$',
    }
    for ch in ['w_fuel', 'w_geo', 'w_reg', 'diff_geo_fuel']:
        sf_m, sf_s, sf_t = sf_ref[ch]
        if ch not in fm_results:
            continue
        r = fm_results[ch]
        mf_m, mf_s, mf_t = r['mean'], r['se_nw'], r['t']
        shrinkage_pct = (1 - abs(mf_m / sf_m)) * 100 if sf_m != 0 else float('nan')
        f.write(
            f'| {labels[ch]} | '
            f'{sf_m:+.4f} ({sf_s:.4f}) [t={sf_t:+.2f}] | '
            f'{mf_m:+.4f} ({mf_s:.4f}) [t={mf_t:+.2f}] | '
            f'{shrinkage_pct:+.1f}% |\n'
        )

    f.write('\n## Specification\n\n')
    f.write('Pre-event regression for each (firm, event):\n\n')
    f.write('$$(r_{it} - rf_t) = \\alpha_i + \\beta_M (Mkt-RF)_t + \\beta_S SMB_t '
            '+ \\beta_V HML_t + \\beta_U (UTL_t - rf_t) + \\epsilon_{it}$$\n\n')
    f.write('estimated on the 24-month window pre-event (minimum 12 obs). '
            'Abnormal return at month $t$ is $AR_{it} = (r_{it}-rf_t) - '
            '\\hat{\\alpha}_i - \\hat{\\beta} \\cdot factors_t$.\n\n')

    f.write('## Notes\n\n')
    f.write(f'- Events with valid FM regressions: {T}\n')
    f.write(f'- Total firm-event observations: {total_obs}\n')
    f.write('- Utility factor: equal-weighted mean of all sample firm returns at each '
            'month (require >= 30 firms).\n')
    f.write('- Currency factor: omitted. The FF Mkt-RF is a US factor; for non-US '
            'firms, a USD trade-weighted index would be appropriate but is not in '
            'the repo. Limitation documented.\n')
    f.write('- Honest DID and lag sensitivity should be re-run on these multi-factor '
            'CARs in a future revision.\n')

_print(f'\nWrote {out_md}')
_print(f'Wrote {out_csv_betas}')
_print(f'Wrote {out_csv_cmp}')
_print('\nDone.')
