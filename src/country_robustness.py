"""Country-level robustness for the non-US fuel coefficient.

Two diagnostics:

  (A) Country event counts (descriptive). Non-US identifying variation may be
      concentrated in a few countries; this table shows where the FM mass
      comes from.

  (B) Leave-one-country-out for non-US. Drop one country at a time, re-run
      the FM regression, and report the resulting coefficient and t-stat.
      A robust non-US result should survive removal of any single country.

  (C) Developed-ex-US vs emerging-market split. Discriminates the
      "informational efficiency" and "ownership structure" explanations
      flagged in Section 4.4 (developed-ex-US markets are nearly as efficient
      as the US, so an effect in this subsample weakens the efficiency
      hypothesis).

Inputs:  data/derived/events/coal_retirement_events.csv (country field)
         (everything else as in robust_inference.py)

Outputs: results/metrics/country_robustness.md
         results/summaries/country_robustness_loo.csv
         results/summaries/country_robustness_split.csv
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


# ── Country classification (MSCI standard) ──────────────────────────

EMERGING_MARKETS = {
    'Brazil', 'Chile', 'China', 'Colombia', 'Czech Republic', 'Egypt',
    'Greece', 'Hungary', 'India', 'Indonesia', 'Korea', 'South Korea',
    'Kuwait', 'Malaysia', 'Mexico', 'Peru', 'Philippines', 'Poland',
    'Qatar', 'Saudi Arabia', 'South Africa', 'Taiwan', 'Thailand',
    'Turkey', 'United Arab Emirates', 'UAE',
}

DEVELOPED_EX_US = {
    'Australia', 'Austria', 'Belgium', 'Canada', 'Denmark', 'Finland',
    'France', 'Germany', 'Hong Kong', 'Ireland', 'Israel', 'Italy',
    'Japan', 'Netherlands', 'New Zealand', 'Norway', 'Portugal',
    'Singapore', 'Spain', 'Sweden', 'Switzerland', 'United Kingdom', 'UK',
}


def classify(country):
    if not country:
        return 'unknown'
    c = country.strip()
    if c in ('United States', 'USA', 'US', 'United States of America'):
        return 'US'
    if c in DEVELOPED_EX_US:
        return 'developed_ex_US'
    if c in EMERGING_MARKETS:
        return 'emerging'
    return 'frontier_or_other'


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


def get_sic4(gvkey):
    f = fundamentals.get(gvkey)
    return f['sic'][:4] if (f and f.get('sic')) else None


events_full = []
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
        country = (row.get('country') or '').strip()
        events_full.append({
            'event_month': ed[:7],
            'gvkeys': row['matched_gvkeys'].split(';'),
            'country': country,
            'classification': classify(country),
        })


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


# ── Per-event obs (compute once, reuse) ─────────────────────────────

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
    per_event_obs[event_id] = obs


# ── FM aggregator ───────────────────────────────────────────────────

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


def fm_subset(event_filter):
    fuel_betas, geo_betas = [], []
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
        'n_events': n_used,
        'fuel_mean': mean_f, 'fuel_se': se_f,
        'fuel_t': mean_f / se_f if (not math.isnan(se_f) and se_f > 0) else float('nan'),
        'geo_mean': mean_g, 'geo_se': se_g,
        'geo_t': mean_g / se_g if (not math.isnan(se_g) and se_g > 0) else float('nan'),
    }


# ── (A) Country event counts (non-US) ───────────────────────────────

country_counts = defaultdict(int)
for e in events_full:
    if e['classification'] != 'US':
        country_counts[e['country'] or 'unknown'] += 1
country_counts_sorted = sorted(country_counts.items(),
                                key=lambda x: -x[1])

# ── (B) Leave-one-country-out on non-US ─────────────────────────────

non_us_full = fm_subset(lambda e: e['classification'] != 'US')
loo_results = [{'dropped': '(none)', **non_us_full}]
for country, count in country_counts_sorted[:10]:  # top 10 countries by event count
    res = fm_subset(lambda e: e['classification'] != 'US' and e['country'] != country)
    loo_results.append({'dropped': country, 'n_dropped': count, **res})

# ── (C) Developed-ex-US vs emerging split ───────────────────────────

dev_ex_us = fm_subset(lambda e: e['classification'] == 'developed_ex_US')
emerging = fm_subset(lambda e: e['classification'] == 'emerging')
frontier = fm_subset(lambda e: e['classification'] == 'frontier_or_other')


# ── Output ──────────────────────────────────────────────────────────

out_loo = results_path('summaries', 'country_robustness_loo.csv')
with open(out_loo, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['dropped_country', 'n_events_dropped', 'n_events_FM',
                'fuel_mean', 'fuel_se', 'fuel_t'])
    for r in loo_results:
        w.writerow([
            r['dropped'], r.get('n_dropped', 0), r['n_events'],
            f'{r["fuel_mean"]:+.6f}',
            f'{r["fuel_se"]:.6f}' if not math.isnan(r['fuel_se']) else 'NA',
            f'{r["fuel_t"]:+.4f}' if not math.isnan(r['fuel_t']) else 'NA',
        ])

out_split = results_path('summaries', 'country_robustness_split.csv')
with open(out_split, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['split', 'n_events',
                'fuel_mean', 'fuel_se', 'fuel_t',
                'geo_mean', 'geo_se', 'geo_t'])
    for label, r in [('developed_ex_US', dev_ex_us),
                     ('emerging', emerging),
                     ('frontier_or_other', frontier)]:
        w.writerow([
            label, r['n_events'],
            f'{r["fuel_mean"]:+.6f}',
            f'{r["fuel_se"]:.6f}' if not math.isnan(r['fuel_se']) else 'NA',
            f'{r["fuel_t"]:+.4f}' if not math.isnan(r['fuel_t']) else 'NA',
            f'{r["geo_mean"]:+.6f}',
            f'{r["geo_se"]:.6f}' if not math.isnan(r['geo_se']) else 'NA',
            f'{r["geo_t"]:+.4f}' if not math.isnan(r['geo_t']) else 'NA',
        ])

# Markdown
out_md = results_path('metrics', 'country_robustness.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Country-Level Robustness for the Non-US Coefficient\n\n')

    f.write('## (A) Non-US event counts by country\n\n')
    f.write('| Country | Events |\n|---|---|\n')
    for c, n in country_counts_sorted[:20]:
        f.write(f'| {c} | {n} |\n')
    if len(country_counts_sorted) > 20:
        f.write(f'| (other countries) | {sum(n for _, n in country_counts_sorted[20:])} |\n')
    f.write(f'\nTotal non-US events: {sum(country_counts.values())}\n\n')

    f.write('## (B) Leave-one-country-out\n\n')
    f.write('| Dropped | N events removed | N FM events | gamma_fuel | SE | t |\n')
    f.write('|---|---|---|---|---|---|\n')
    for r in loo_results:
        nd = r.get('n_dropped', 0)
        f.write(
            f'| {r["dropped"]} | {nd} | {r["n_events"]} | '
            f'{r["fuel_mean"]:+.4f} | {r["fuel_se"]:.4f} | {r["fuel_t"]:+.2f} |\n'
        )

    f.write('\n## (C) Developed-ex-US vs Emerging vs Frontier\n\n')
    f.write('| Split | N events | gamma_fuel | SE | t | gamma_geo | SE | t |\n')
    f.write('|---|---|---|---|---|---|---|---|\n')
    for label, r in [('Developed (ex-US)', dev_ex_us),
                     ('Emerging', emerging),
                     ('Frontier / Other', frontier)]:
        f.write(
            f'| {label} | {r["n_events"]} | '
            f'{r["fuel_mean"]:+.4f} | {r["fuel_se"]:.4f} | {r["fuel_t"]:+.2f} | '
            f'{r["geo_mean"]:+.4f} | {r["geo_se"]:.4f} | {r["geo_t"]:+.2f} |\n'
        )

_print(f'Wrote {out_loo}')
_print(f'Wrote {out_split}')
_print(f'Wrote {out_md}')
_print('\nDone.')
