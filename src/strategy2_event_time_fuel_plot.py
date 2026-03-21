"""Event-time plot: High vs Low fuel-similarity neighbours.

Compute demeaned monthly abnormal returns by event time (tau) for
high vs low fuel-similarity neighbours around first-mover events.
Outputs:
  - results/summaries/strategy2_event_time_fuel.csv
  - finance_draft/figures/event_time_fuel_similarity.pdf
"""
import csv
import os
import math
from collections import defaultdict

from _paths import derived_path, raw_path, results_path

EVENT_SCOPE = 'first_mover'
EXACT_ONLY = False
TAU_START = -6
TAU_END = 12
MIN_PRE = 3  # minimum pre-event months for demeaning


def add_months(ym, delta):
    y, m = ym.split('-')
    y = int(y)
    m = int(m)
    m += delta
    while m > 12:
        y += 1
        m -= 12
    while m < 1:
        y -= 1
        m += 12
    return f"{y:04d}-{m:02d}"


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
            vwretd_dec = (mktrf_val + rf_val) / 100.0
            date_fmt = f'{date[:4]}-{date[4:6]}'
            vwretd[date_fmt] = vwretd_dec
    return vwretd


def is_exact_source(src):
    if not src:
        return False
    s = src.lower()
    if 'proxy' in s or 'approx' in s or 'mid' in s or 'month' in s:
        return False
    return True


# Load monthly returns
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            monthly_ret[row['gvkey']][row['datadate'][:7]] = float(row['ret_monthly'])
        except (ValueError, TypeError):
            continue
print(f'Monthly returns: {len(monthly_ret)} firms')

vwretd = load_ff_factors_monthly(raw_path('factors', 'F-F_Research_Data_Factors.csv'))
if not vwretd:
    raise SystemExit('Missing Fama-French monthly factors.')
print(f'Market monthly months: {len(vwretd)}')

# Load fuel-similarity weights (row-normalized)
W_fuel = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_fuel.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gi = row['gvkey_i']
        gj = row['gvkey_j']
        try:
            W_fuel[gi][gj] = float(row['w_ij'])
        except (ValueError, TypeError):
            continue
print(f'Fuel W rows: {len(W_fuel)}')

# Load events
all_events = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        ann_date = row.get('announcement_date', '').strip()
        ret_date = row.get('event_date', '').strip()
        ann_src = row.get('announcement_source', '').strip()
        if EXACT_ONLY:
            if not ann_date:
                continue
            if ann_src and not is_exact_source(ann_src):
                continue
        effective_date = ann_date if ann_date else ret_date
        event_year = None
        if effective_date and len(effective_date) >= 4 and effective_date[:4].isdigit():
            event_year = int(effective_date[:4])
        else:
            event_year = int(row['ret_year']) if row.get('ret_year') else None
        all_events.append({
            'event_date': effective_date,
            'year': event_year,
            'gvkeys': row['matched_gvkeys'].split(';'),
            'is_first_mover': row.get('is_first_mover') == 'True',
        })

events = [e for e in all_events if (e['is_first_mover'] if EVENT_SCOPE == 'first_mover' else True)]
print(f'Events used: {len(events)}')

# Aggregate demeaned ARs by tau and group
group_tau = {'high': defaultdict(list), 'low': defaultdict(list)}

for event in events:
    event_gvkeys = set(event['gvkeys'])
    event_date = event.get('event_date', '')
    if event_date and len(event_date) >= 7:
        event_month = event_date[:7]
    else:
        year = event.get('year')
        event_month = f'{year}-07' if year else None
    if not event_month:
        continue

    for fm_gk in event_gvkeys:
        if fm_gk not in W_fuel:
            continue
        neighbors = W_fuel[fm_gk]
        neighbor_gks = [gk for gk in neighbors.keys() if gk not in event_gvkeys]
        if not neighbor_gks:
            continue

        weights = [neighbors[gk] for gk in neighbor_gks if neighbors.get(gk, 0.0) > 0]
        if not weights:
            continue
        weights.sort()
        median_w = weights[len(weights) // 2]

        for gk in neighbor_gks:
            w_f = neighbors.get(gk, 0.0)
            if w_f <= 0:
                continue
            group = 'high' if w_f >= median_w else 'low'
            if gk not in monthly_ret:
                continue

            # Pre-event demeaning
            pre_vals = []
            for tau in range(TAU_START, 0):
                ym = add_months(event_month, tau)
                if ym in monthly_ret[gk] and ym in vwretd:
                    ar = monthly_ret[gk][ym] - vwretd[ym]
                    pre_vals.append(ar)
            if len(pre_vals) < MIN_PRE:
                continue
            pre_mean = sum(pre_vals) / len(pre_vals)

            for tau in range(TAU_START, TAU_END + 1):
                ym = add_months(event_month, tau)
                if ym not in monthly_ret[gk] or ym not in vwretd:
                    continue
                ar = monthly_ret[gk][ym] - vwretd[ym]
                ar_demean = ar - pre_mean
                group_tau[group][tau].append(ar_demean)

# Write summary CSV
out_csv = results_path('summaries', 'strategy2_event_time_fuel.csv')
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['group', 'tau', 'mean_ar', 'se', 'n'])
    for group in ['low', 'high']:
        for tau in range(TAU_START, TAU_END + 1):
            vals = group_tau[group].get(tau, [])
            n = len(vals)
            if n == 0:
                continue
            mean = sum(vals) / n
            if n > 1:
                var = sum((v - mean) ** 2 for v in vals) / (n - 1)
                se = math.sqrt(var / n)
            else:
                se = 0.0
            w.writerow([group, tau, f'{mean:.6f}', f'{se:.6f}', n])

print(f'Wrote event-time summary: {out_csv}')

# Plot
try:
    import matplotlib.pyplot as plt
    import numpy as np
except Exception as e:
    print(f'Plotting skipped (matplotlib not available): {e}')
    raise SystemExit(0)

def series_from_group(group):
    taus = []
    means = []
    ses = []
    for tau in range(TAU_START, TAU_END + 1):
        vals = group_tau[group].get(tau, [])
        if not vals:
            continue
        n = len(vals)
        mean = sum(vals) / n
        var = sum((v - mean) ** 2 for v in vals) / (n - 1) if n > 1 else 0.0
        se = math.sqrt(var / n) if n > 1 else 0.0
        taus.append(tau)
        means.append(mean)
        ses.append(se)
    return np.array(taus), np.array(means), np.array(ses)

tau_l, mean_l, se_l = series_from_group('low')
tau_h, mean_h, se_h = series_from_group('high')

fig, ax = plt.subplots(figsize=(6.5, 4))
ax.axhline(0, color='black', linewidth=0.8)
ax.axvline(0, color='gray', linestyle='--', linewidth=0.8)

ax.plot(tau_l, mean_l, label='Low fuel similarity', color='#1f77b4')
ax.fill_between(tau_l, mean_l - 1.96 * se_l, mean_l + 1.96 * se_l, color='#1f77b4', alpha=0.2)

ax.plot(tau_h, mean_h, label='High fuel similarity', color='#d62728')
ax.fill_between(tau_h, mean_h - 1.96 * se_h, mean_h + 1.96 * se_h, color='#d62728', alpha=0.2)

ax.set_xlabel('Event time (months)')
ax.set_ylabel('Demeaned abnormal return')
ax.set_title('Event-time ARs by Fuel Similarity')
ax.legend(frameon=False)
ax.grid(True, alpha=0.2)

fig_dir = os.path.join(os.path.dirname(derived_path('dummy')), '..', 'finance_draft', 'figures')
fig_dir = os.path.abspath(fig_dir)
os.makedirs(fig_dir, exist_ok=True)
out_pdf = os.path.join(fig_dir, 'event_time_fuel_similarity.pdf')
out_png = os.path.join(fig_dir, 'event_time_fuel_similarity.png')
fig.tight_layout()
fig.savefig(out_pdf)
fig.savefig(out_png, dpi=200)
print(f'Saved figure: {out_pdf}')
