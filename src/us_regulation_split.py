"""US restructured-vs-regulated split for the saturation hypothesis.

Section 4.4 of the manuscript reports that the fuel coefficient is near-zero
in the United States. The paper currently lists four candidate explanations
(saturation, regulation, informational efficiency, ownership) without
selecting among them. This script implements the regulation test: split the
US event sample by whether the retiring plant operates in a restructured
retail market (electricity is sold competitively, ISO/RTO wholesale) or in
a vertically integrated state with rate-of-return regulation (PUC).

Hypothesis: under rate-of-return regulation, utilities can pass stranded
costs through to ratepayers, dampening the equity-price reaction. In
restructured markets, the cost falls on the equity holder, so the channel
should be detectable.

Restructured states/jurisdictions (academic standard): CA, CT, DE, DC, IL,
ME, MD, MA, NH, NJ, NY, OH, PA, RI, TX. All other US states are classified
as regulated.

Inputs:  data/derived/gem/gem_coal.csv (state of retiring plant)
         data/derived/events/coal_retirement_events.csv
         (everything else as in robust_inference.py)

Outputs: results/metrics/us_regulation_split.md
         results/summaries/us_regulation_split.csv
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


# ── State classification ────────────────────────────────────────────

RESTRUCTURED_STATES = {
    'California', 'Connecticut', 'Delaware', 'District of Columbia',
    'Illinois', 'Maine', 'Maryland', 'Massachusetts', 'New Hampshire',
    'New Jersey', 'New York', 'Ohio', 'Pennsylvania', 'Rhode Island',
    'Texas',
}

# State abbreviation forms sometimes appear in GEM data
RESTRUCTURED_ABBREVS = {
    'CA', 'CT', 'DE', 'DC', 'IL', 'ME', 'MD', 'MA', 'NH', 'NJ', 'NY',
    'OH', 'PA', 'RI', 'TX',
}


def classify_state(state_str):
    """Return 'restructured', 'regulated', or 'unknown'."""
    if not state_str:
        return 'unknown'
    s = state_str.strip()
    if s in RESTRUCTURED_STATES or s.upper() in RESTRUCTURED_ABBREVS:
        return 'restructured'
    # Match common variants
    for full, abbrev in zip(['California', 'Texas', 'Pennsylvania'],
                             ['Calif', 'Tex', 'Penn']):
        if abbrev.lower() in s.lower():
            return 'restructured'
    return 'regulated'


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


# Load GEM coal plant -> state mapping (by plant name)
plant_state = {}
with open(derived_path('gem', 'gem_coal.csv'), 'r', encoding='utf-8') as f:
    reader = csv.DictReader(f)
    for row in reader:
        country = (row.get('Country/Area') or '').strip()
        if country not in ('United States', 'United States of America', 'USA', 'US'):
            continue
        plant_name = (row.get('Plant name') or '').strip()
        state = (row.get('Subnational unit (province, state)') or '').strip()
        if plant_name and state:
            plant_state[plant_name] = state
_print(f'GEM US coal plants with state: {len(plant_state)}')


# Load events with retirement plant state
events_full = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        if row.get('is_first_mover') != 'True':
            continue
        ann = row.get('announcement_date', '').strip()
        ret = row.get('event_date', '').strip()
        effective_date = ann if (ann and len(ann) >= 7) else ret
        if not effective_date or len(effective_date) < 7:
            continue
        country = (row.get('country') or '').strip()
        plant_name = (row.get('plant_name') or '').strip()
        state = plant_state.get(plant_name, '')
        if country in ('United States', 'United States of America', 'USA', 'US'):
            regime = classify_state(state)
        else:
            regime = 'non_US'
        events_full.append({
            'plant': plant_name,
            'event_month': effective_date[:7],
            'gvkeys': row['matched_gvkeys'].split(';'),
            'country': country,
            'state': state,
            'regime': regime,
        })

n_us = sum(1 for e in events_full if e['regime'] in ('restructured', 'regulated'))
n_restruct = sum(1 for e in events_full if e['regime'] == 'restructured')
n_reg = sum(1 for e in events_full if e['regime'] == 'regulated')
n_unknown_us = sum(1 for e in events_full
                   if (e['country'] in ('United States', 'United States of America', 'USA', 'US')
                       and e['regime'] == 'unknown'))
n_non_us = sum(1 for e in events_full if e['regime'] == 'non_US')
_print(f'Events: {len(events_full)} total; '
       f'US restructured = {n_restruct}, US regulated = {n_reg}, '
       f'US unknown state = {n_unknown_us}, non-US = {n_non_us}')


# ── CAR helper (single-factor, matches headline) ────────────────────

PRE_MONTHS = 24
POST_MONTHS = 3


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


# ── Per-event candidate sets (computed once, reused per regime) ────

MIN_OBS = 20
SPEC_VARS = ['w_geo', 'w_fuel', 'w_reg', 'same_sector']

per_event_obs = {}
for event_id, event in enumerate(events_full):
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
        candidate_firms = list(neighbor_gks) + ctrl_sample
        for gk in candidate_firms:
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
    per_event_obs[event_id] = obs


# ── FM aggregation per regime subset ───────────────────────────────

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


def fm_subset(event_filter, label):
    fuel_betas = []
    geo_betas = []
    n_used = 0
    for event_id, event in enumerate(events_full):
        if not event_filter(event):
            continue
        obs = per_event_obs.get(event_id, [])
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
        geo_betas.append(bd.get('w_geo', float('nan')))
        n_used += 1
    mean_f, se_f, _ = newey_west_se(fuel_betas, lag=4)
    mean_g, se_g, _ = newey_west_se(geo_betas, lag=4)
    return {
        'label': label, 'n_events_FM': n_used,
        'fuel_mean': mean_f, 'fuel_se': se_f,
        'fuel_t': mean_f / se_f if (not math.isnan(se_f) and se_f > 0) else float('nan'),
        'geo_mean': mean_g, 'geo_se': se_g,
        'geo_t': mean_g / se_g if (not math.isnan(se_g) and se_g > 0) else float('nan'),
    }


splits = [
    ('US: All', lambda e: e['regime'] in ('restructured', 'regulated')),
    ('US: Restructured', lambda e: e['regime'] == 'restructured'),
    ('US: Regulated', lambda e: e['regime'] == 'regulated'),
    ('Non-US', lambda e: e['regime'] == 'non_US'),
]

results = [fm_subset(filt, lbl) for lbl, filt in splits]


# ── Output ─────────────────────────────────────────────────────────

out_csv = results_path('summaries', 'us_regulation_split.csv')
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['split', 'n_events_FM',
                'fuel_mean', 'fuel_se', 'fuel_t',
                'geo_mean', 'geo_se', 'geo_t'])
    for r in results:
        w.writerow([
            r['label'], r['n_events_FM'],
            f'{r["fuel_mean"]:+.6f}',
            f'{r["fuel_se"]:.6f}' if not math.isnan(r['fuel_se']) else 'NA',
            f'{r["fuel_t"]:+.4f}' if not math.isnan(r['fuel_t']) else 'NA',
            f'{r["geo_mean"]:+.6f}',
            f'{r["geo_se"]:.6f}' if not math.isnan(r['geo_se']) else 'NA',
            f'{r["geo_t"]:+.4f}' if not math.isnan(r['geo_t']) else 'NA',
        ])
_print(f'Wrote {out_csv}')

out_md = results_path('metrics', 'us_regulation_split.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# US Restructured vs Regulated Market Split\n\n')
    f.write('Tests the regulation hypothesis for the US-null finding in '
            'Section 4.4. The headline US fuel coefficient is approximately zero. '
            'If the channel is dampened by rate-of-return regulation (which '
            'allows utilities to pass stranded costs through to ratepayers), the '
            'channel should be detectable in restructured retail markets where '
            'equity holders bear the cost.\n\n')
    f.write('Restructured states (academic standard, ~15 states): '
            + ', '.join(sorted(RESTRUCTURED_STATES)) + '. '
            'All other US states are classified as regulated.\n\n')

    f.write('## Subsample results (FM + NW lag 4, single-factor CARs)\n\n')
    f.write('| Split | N events (FM) | gamma_fuel | SE | t | gamma_geo | SE | t |\n')
    f.write('|---|---|---|---|---|---|---|---|\n')
    for r in results:
        f.write(
            f'| {r["label"]} | {r["n_events_FM"]} | '
            f'{r["fuel_mean"]:+.4f} | {r["fuel_se"]:.4f} | {r["fuel_t"]:+.2f} | '
            f'{r["geo_mean"]:+.4f} | {r["geo_se"]:.4f} | {r["geo_t"]:+.2f} |\n'
        )

    f.write('\n## Interpretation\n\n')
    res_dict = {r['label']: r for r in results}
    us_all = res_dict.get('US: All', {})
    us_restruct = res_dict.get('US: Restructured', {})
    us_reg = res_dict.get('US: Regulated', {})
    non_us = res_dict.get('Non-US', {})

    if (us_restruct and us_reg
            and not math.isnan(us_restruct.get('fuel_t', float('nan')))
            and not math.isnan(us_reg.get('fuel_t', float('nan')))):
        f.write(f'- Non-US: $\\hat\\gamma_{{\\mathrm{{fuel}}}} = '
                f'{non_us["fuel_mean"]:+.2f}$ ($t={non_us["fuel_t"]:+.2f}$).\n')
        f.write(f'- US (all): $\\hat\\gamma_{{\\mathrm{{fuel}}}} = '
                f'{us_all["fuel_mean"]:+.2f}$ ($t={us_all["fuel_t"]:+.2f}$).\n')
        f.write(f'- US restructured: $\\hat\\gamma_{{\\mathrm{{fuel}}}} = '
                f'{us_restruct["fuel_mean"]:+.2f}$ ($t={us_restruct["fuel_t"]:+.2f}$).\n')
        f.write(f'- US regulated: $\\hat\\gamma_{{\\mathrm{{fuel}}}} = '
                f'{us_reg["fuel_mean"]:+.2f}$ ($t={us_reg["fuel_t"]:+.2f}$).\n\n')

        if (abs(us_restruct['fuel_t']) > abs(us_reg['fuel_t'])
                and us_restruct['fuel_mean'] < us_reg['fuel_mean']):
            f.write('The channel is more pronounced in restructured states than in '
                    'regulated states, consistent with the regulation hypothesis: '
                    'rate-of-return ratemaking attenuates the equity-price reaction '
                    'because regulated utilities can recover stranded costs through '
                    'PUC-approved rates.\n')
        elif abs(us_reg['fuel_t']) > abs(us_restruct['fuel_t']):
            f.write('The channel is at least as strong in regulated states as in '
                    'restructured ones; the regulation hypothesis is not confirmed. '
                    'The US-null pattern likely reflects other mechanisms '
                    '(saturation, informational efficiency, ownership structure).\n')
        else:
            f.write('The two US subsamples produce coefficients of similar '
                    'magnitude and significance; the regulation hypothesis is not '
                    'sharply supported by this split.\n')

_print(f'Wrote {out_md}')
_print('\nDone.')
