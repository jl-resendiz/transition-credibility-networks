"""Diagnostics for coal-share treatment intensity around phase-out shocks.

Uses firm-level coal_share (coal_mw / total_mw) from firm_alpha_panel.csv.
Treated firms are those matched to each shock (country/state).

Outputs: results/metrics/strategy3_phaseout_coalshare_diagnostics.md
"""
import csv
import math
import os
from collections import Counter, defaultdict

from _paths import derived_path, results_path

EVENTS_PATH = derived_path('events', 'coal_phaseout_shocks_events.csv')
ALPHA_PANEL = derived_path('fundamentals', 'firm_alpha_panel.csv')
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


def load_coal_share(panel_path):
    years_by_gvkey = defaultdict(list)
    coal_by_year = defaultdict(dict)
    with open(panel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            yr = row['year']
            cs = row.get('coal_share', '')
            if not cs:
                continue
            try:
                cs_val = float(cs)
            except (ValueError, TypeError):
                continue
            coal_by_year[gk][int(yr)] = cs_val
            years_by_gvkey[gk].append(int(yr))
    for gk in years_by_gvkey:
        years_by_gvkey[gk] = sorted(set(years_by_gvkey[gk]))

    def get_share(gk, year):
        if gk not in coal_by_year or not years_by_gvkey[gk]:
            return None
        if year in coal_by_year[gk]:
            return coal_by_year[gk][year]
        years = years_by_gvkey[gk]
        # nearest year (prefer <= year)
        prior = [y for y in years if y <= year]
        if prior:
            return coal_by_year[gk][max(prior)]
        # otherwise nearest later year
        return coal_by_year[gk][years[0]]

    return get_share


# Load coal_share lookup
get_coal_share = load_coal_share(ALPHA_PANEL)

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
        if not row.get('event_date'):
            continue
        events.append(row)

years = Counter()
tiers = Counter()
bindings = Counter()
treated_counts = []
coverage_ok = 0
obs_pairs = 0
coal_treated = []
coal_controls = []

for e in events:
    event_date = e['event_date']
    if len(event_date) >= 7:
        event_month = event_date[:7]
        event_year = int(event_date[:4])
    else:
        continue
    years[event_year] += 1
    tiers[e.get('exogeneity_tier', '')] += 1
    bindings[e.get('binding', '')] += 1

    gvkeys = [g for g in e['matched_gvkeys'].split(';') if g]
    treated_counts.append(len(gvkeys))
    treated_set = set(gvkeys)

    # Build candidate firms: treated + controls
    non_treated = [gk for gk in fundamentals if gk not in treated_set]
    n_ctrl = min(len(non_treated), max(CONTROL_MULT * len(treated_set), 20))
    ctrl_sample = non_treated[:n_ctrl]
    candidate_firms = list(treated_set) + ctrl_sample

    for gk in candidate_firms:
        cs = get_coal_share(gk, event_year)
        if cs is None:
            continue
        if gk in treated_set:
            coal_treated.append(cs)
        else:
            coal_controls.append(cs)
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
coal_t_stats = stats(coal_treated)
coal_c_stats = stats(coal_controls)

lines = [
    '# Coal Phase-Out Coal-Share Diagnostics',
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
        f"- mean={treated_stats['mean']:.2f}, sd={treated_stats['sd']:.2f}, "
        f"p25={treated_stats['p25']}, p50={treated_stats['p50']}, "
        f"p75={treated_stats['p75']}, p95={treated_stats['p95']}",
    ]

lines += [
    '',
    '## Coal-share distribution (treated vs controls)',
]
if coal_t_stats:
    lines.append(
        f"- treated: mean={coal_t_stats['mean']:.4f}, sd={coal_t_stats['sd']:.4f}, "
        f"p25={coal_t_stats['p25']:.4f}, p50={coal_t_stats['p50']:.4f}, "
        f"p75={coal_t_stats['p75']:.4f}, p95={coal_t_stats['p95']:.4f}"
    )
if coal_c_stats:
    lines.append(
        f"- control: mean={coal_c_stats['mean']:.4f}, sd={coal_c_stats['sd']:.4f}, "
        f"p25={coal_c_stats['p25']:.4f}, p50={coal_c_stats['p50']:.4f}, "
        f"p75={coal_c_stats['p75']:.4f}, p95={coal_c_stats['p95']:.4f}"
    )

lines += [
    '',
    '## Coverage',
    f'- firm-event pairs (treated + controls with coal_share): {obs_pairs}',
    f'- monthly return coverage (>=60% window): {coverage_ok} ({(coverage_ok/obs_pairs*100):.1f}%)' if obs_pairs else '- monthly return coverage: NA',
]

out_path = results_path('metrics', 'strategy3_phaseout_coalshare_diagnostics.md')
os.makedirs(os.path.dirname(out_path), exist_ok=True)
with open(out_path, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines))

print(f'Wrote: {out_path}')
