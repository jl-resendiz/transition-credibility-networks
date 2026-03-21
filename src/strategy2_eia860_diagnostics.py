"""Diagnostics for EIA-860-style announcement shocks.

Outputs a compact summary of:
  - Event coverage (counts, exact vs month-level)
  - Exposure mapping (neighbor counts, exposure distribution)
  - Return coverage (monthly availability around events)
  - Exposure transformations (base, log1p, zscore)

Results: finance_data/results/metrics/strategy2_eia860_diagnostics.md
"""
import csv
import math
import os
from collections import defaultdict, Counter

from _paths import derived_path, results_path, raw_path

EVENTS_PATH = os.getenv('EVENTS_PATH', '')
EVENT_SCOPE = os.getenv('EVENT_SCOPE', 'all_matched')
CONTROL_MULT = int(os.getenv('CONTROL_MULT', '5'))
TAU_START = int(os.getenv('TAU_START', '-6'))
TAU_END = int(os.getenv('TAU_END', '12'))


def pct(v):
    return f"{100.0 * v:.1f}%"


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


def load_monthly_returns():
    monthly = defaultdict(dict)
    with open(derived_path('returns', 'monthly_returns.csv'), 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            ym = row['datadate'][:7]
            try:
                monthly[gk][ym] = float(row['ret_monthly'])
            except (ValueError, TypeError):
                continue
    return monthly


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


# Load W (geographic exposure)
W = defaultdict(dict)
with open(derived_path('networks', 'weight_matrix_W_geo.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gi = row['gvkey_i']
        gj = row['gvkey_j']
        try:
            W[gi][gj] = float(row['w_ij'])
        except (ValueError, TypeError):
            continue

# Load fundamentals for control sampling
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

# Load events
events_path = EVENTS_PATH if EVENTS_PATH else derived_path('events', 'eia860_announcement_events.csv')
all_events = []
with open(events_path, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        ann_date = row.get('announcement_date', '').strip()
        if not ann_date:
            continue
        effective_date = ann_date
        if len(effective_date) >= 7:
            event_month = effective_date[:7]
        else:
            continue
        all_events.append({
            'event_month': event_month,
            'announcement_date': ann_date,
            'announcement_source': row.get('announcement_source', ''),
            'gvkeys': [g.strip() for g in row['matched_gvkeys'].split(';') if g.strip()],
            'is_first_mover': row.get('is_first_mover') == 'True',
        })

if EVENT_SCOPE == 'first_mover':
    events = [e for e in all_events if e['is_first_mover']]
else:
    events = list(all_events)

# Coverage diagnostics
exact_cnt = sum(1 for e in events if len(e['announcement_date']) >= 10)
month_cnt = sum(1 for e in events if len(e['announcement_date']) == 7)
event_years = Counter(int(e['announcement_date'][:4]) for e in events if e['announcement_date'][:4].isdigit())
unique_gvkeys = set(gk for e in events for gk in e['gvkeys'])

# Exposure mapping
neighbor_counts = []
exposure_vals = []
obs_pairs = 0
monthly_ret = load_monthly_returns()
coverage_ok = 0

for event_id, e in enumerate(events):
    event_month = e['event_month']
    event_gvkeys = set(e['gvkeys'])
    for fm_gk in event_gvkeys:
        if fm_gk not in W:
            continue
        neighbors = W[fm_gk]
        neighbor_gks = set(neighbors.keys()) - event_gvkeys
        neighbor_counts.append(len(neighbor_gks))
        non_connected = [gk for gk in fundamentals if gk not in event_gvkeys and gk not in neighbors]
        n_ctrl = min(len(non_connected), max(CONTROL_MULT * len(neighbor_gks), 20))
        ctrl_sample = non_connected[:n_ctrl]
        candidate_firms = list(neighbor_gks) + ctrl_sample
        for gk in candidate_firms:
            w_ij = neighbors.get(gk, 0.0)
            exposure_vals.append(w_ij)
            obs_pairs += 1
            # return coverage check: at least 60% of months in window
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

# Exposure transforms
base_stats = stats(exposure_vals)
log1p_stats = stats([math.log1p(v) if v > 0 else 0.0 for v in exposure_vals])
if exposure_vals:
    mean = base_stats['mean']
    std = base_stats['sd'] if base_stats['sd'] > 1e-12 else 1.0
    z_stats = stats([(v - mean) / std for v in exposure_vals])
else:
    z_stats = {}

lines = [
    '# EIA-860 Announcement Shock Diagnostics',
    '',
    f'- events_path: {events_path}',
    f'- event_scope: {EVENT_SCOPE}',
    f'- events_total: {len(events)}',
    f'- unique_gvkeys: {len(unique_gvkeys)}',
    f'- exact_day: {exact_cnt}',
    f'- month_only: {month_cnt}',
    '',
    '## Event year counts (top 10)',
]
for yr, cnt in event_years.most_common(10):
    lines.append(f'- {yr}: {cnt}')

lines += [
    '',
    '## Exposure mapping',
    f'- firm-event pairs (neighbors + controls): {obs_pairs}',
    f'- avg neighbors per event: {sum(neighbor_counts)/len(neighbor_counts):.2f}' if neighbor_counts else '- avg neighbors per event: NA',
    f'- median neighbors per event: {sorted(neighbor_counts)[len(neighbor_counts)//2]}' if neighbor_counts else '- median neighbors per event: NA',
    f'- monthly return coverage (>=60% window): {coverage_ok} ({pct(coverage_ok / obs_pairs) if obs_pairs else "NA"})',
    '',
    '## Exposure distribution',
]
if base_stats:
    lines += [
        f'- base mean={base_stats["mean"]:.4f}, sd={base_stats["sd"]:.4f}, p25={base_stats["p25"]:.4f}, p50={base_stats["p50"]:.4f}, p75={base_stats["p75"]:.4f}, p95={base_stats["p95"]:.4f}',
        f'- log1p mean={log1p_stats["mean"]:.4f}, sd={log1p_stats["sd"]:.4f}, p25={log1p_stats["p25"]:.4f}, p50={log1p_stats["p50"]:.4f}, p75={log1p_stats["p75"]:.4f}, p95={log1p_stats["p95"]:.4f}',
        f'- zscore mean={z_stats["mean"]:.4f}, sd={z_stats["sd"]:.4f}, p25={z_stats["p25"]:.4f}, p50={z_stats["p50"]:.4f}, p75={z_stats["p75"]:.4f}, p95={z_stats["p95"]:.4f}',
    ]

out_path = results_path('metrics', 'strategy2_eia860_diagnostics.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

print(f'Wrote: {out_path}')
