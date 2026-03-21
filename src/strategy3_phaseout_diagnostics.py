"""Diagnostics for coal phase-out regulatory shocks.

Outputs: finance_data/results/metrics/strategy3_phaseout_diagnostics.md
"""
import csv
import math
import os
from collections import Counter, defaultdict

from _paths import derived_path, results_path

EVENTS_PATH = derived_path('events', 'coal_phaseout_shocks_events.csv')
TAU_START = int(os.getenv('TAU_START', '-6'))
TAU_END = int(os.getenv('TAU_END', '12'))
CONTROL_MULT = int(os.getenv('CONTROL_MULT', '5'))


def stats(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return {}
    vals.sort()
    n = len(vals)
    mean = sum(vals) / n
    var = sum((v - mean) ** 2 for v in vals) / (n - 1) if n > 1 else 0.0
    sd = math.sqrt(var)
    return {
        'n': n,
        'mean': mean,
        'sd': sd,
        'p25': vals[n // 4],
        'p50': vals[n // 2],
        'p75': vals[(3 * n) // 4],
        'p95': vals[int(0.95 * (n - 1))],
    }


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


# Load returns for coverage
monthly_ret = defaultdict(dict)
with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        ym = row['datadate'][:7]
        try:
            monthly_ret[gk][ym] = float(row['ret_monthly'])
        except (ValueError, TypeError):
            continue

# Load W
W = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gi = row['gvkey_i']
        gj = row['gvkey_j']
        try:
            W[gi][gj] = float(row['w_ij'])
        except (ValueError, TypeError):
            continue

# Load fundamentals for control pool
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

events = []
with open(EVENTS_PATH, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        events.append(row)

years = Counter()
tiers = Counter()
bindings = Counter()
treated_counts = []
coverage_ok = 0
obs_pairs = 0
neighbor_counts = []
exposure_vals = []

for e in events:
    event_date = e['event_date']
    if len(event_date) >= 7:
        event_month = event_date[:7]
    else:
        continue
    if event_date[:4].isdigit():
        years[int(event_date[:4])] += 1
    tiers[e.get('exogeneity_tier', '')] += 1
    bindings[e.get('binding', '')] += 1

    gvkeys = [g for g in e['matched_gvkeys'].split(';') if g]
    treated_counts.append(len(gvkeys))

    for fm_gk in gvkeys:
        if fm_gk not in W:
            continue
        neighbors = W[fm_gk]
        neighbor_gks = set(neighbors.keys()) - set(gvkeys)
        neighbor_counts.append(len(neighbor_gks))
        non_connected = [gk for gk in fundamentals if gk not in gvkeys and gk not in neighbors]
        n_ctrl = min(len(non_connected), max(CONTROL_MULT * len(neighbor_gks), 20))
        ctrl_sample = non_connected[:n_ctrl]
        candidate_firms = list(neighbor_gks) + ctrl_sample
        for gk in candidate_firms:
            w_ij = neighbors.get(gk, 0.0)
            exposure_vals.append(w_ij)
            obs_pairs += 1
            if gk in monthly_ret:
                need = 0
                have = 0
                for tau in range(TAU_START, TAU_END + 1):
                    ym = add_months(event_month, tau)
                    need += 1
                    if ym in monthly_ret[gk]:
                        have += 1
                if need > 0 and have / need >= 0.6:
                    coverage_ok += 1

treated_stats = stats(treated_counts)
neighbor_stats = stats(neighbor_counts)
exposure_stats = stats(exposure_vals)

lines = [
    '# Coal Phase-Out Shock Diagnostics',
    '',
    f'- events_total: {len(events)}',
    f'- unique_tiers: {dict(tiers)}',
    f'- binding: {dict(bindings)}',
    '',
    '## Event years (top 10)',
]
for yr, cnt in years.most_common(10):
    lines.append(f'- {yr}: {cnt}')

lines += [
    '',
    '## Treated firms per shock',
]
if treated_stats:
    lines += [
        f"- mean={treated_stats['mean']:.2f}, sd={treated_stats['sd']:.2f}, p25={treated_stats['p25']}, p50={treated_stats['p50']}, p75={treated_stats['p75']}, p95={treated_stats['p95']}",
    ]

lines += [
    '',
    '## Exposure mapping',
    f'- firm-event pairs (neighbors + controls): {obs_pairs}',
    f'- avg neighbors per event: {neighbor_stats["mean"]:.2f}' if neighbor_stats else '- avg neighbors per event: NA',
    f'- monthly return coverage (>=60% window): {coverage_ok} ({(coverage_ok/obs_pairs*100):.1f}%)' if obs_pairs else '- monthly return coverage: NA',
    '',
    '## Exposure distribution (W_geo)',
]
if exposure_stats:
    lines.append(f"- mean={exposure_stats['mean']:.4f}, sd={exposure_stats['sd']:.4f}, p25={exposure_stats['p25']:.4f}, p50={exposure_stats['p50']:.4f}, p75={exposure_stats['p75']:.4f}, p95={exposure_stats['p95']:.4f}")

out_path = results_path('metrics', 'strategy3_phaseout_diagnostics.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

print(f'Wrote: {out_path}')
