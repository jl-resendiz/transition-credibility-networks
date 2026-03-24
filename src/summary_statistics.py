"""Generate summary statistics table for the paper.

Writes both to stdout and to results/summaries/summary_statistics.md
so that all manuscript claims are traceable to a pipeline output.
"""
import csv, os, math
from collections import defaultdict, Counter

from _paths import derived_path, results_path, raw_path

# ── Load firm fundamentals (latest year per firm) ──
firms = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in firms or fy > firms[gk]['fyear']:
            firms[gk] = row

# ── Load density (raw weight sum preferred) ──
density = {}
with open(derived_path('networks', 'firm_centroids.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        wsum = row.get('w_sum', '')
        if wsum not in ('', None):
            try:
                density[gk] = float(wsum)
            except (ValueError, TypeError):
                density[gk] = float(row.get('n_neighbors', 0))
        else:
            density[gk] = float(row.get('n_neighbors', 0))

# ── Load GEM capacity totals per firm ──
gem_capacity = defaultdict(lambda: {'fossil_mw': 0, 'total_mw': 0})
gem_match_path = derived_path('mappings', 'gem_compustat_matches.csv')
if os.path.exists(gem_match_path):
    with open(gem_match_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            coal = float(row.get('coal_mw', 0) or 0)
            gas = float(row.get('gas_mw', 0) or 0)
            solar = float(row.get('solar_mw', 0) or 0)
            wind = float(row.get('wind_mw', 0) or 0)
            gem_capacity[gk]['fossil_mw'] += coal + gas
            gem_capacity[gk]['total_mw'] += coal + gas + solar + wind

# ── Build sample ──
sample = []
for gk, f in firms.items():
    d = {}
    d['gvkey'] = gk
    d['conm'] = f['conm']
    d['fic'] = f['fic']

    for var in ['at', 'sale', 'ppent']:
        try:
            d[var] = float(f[var]) if f.get(var) else None
        except (ValueError, TypeError):
            d[var] = None

    for var in ['alpha', 'lambda', 'rho', 'kappa', 'delta']:
        try:
            d[var] = float(f[var]) if f.get(var) else None
        except (ValueError, TypeError):
            d[var] = None

    d['density'] = density.get(gk)
    d['has_gem'] = d['alpha'] is not None
    d['total_mw'] = gem_capacity.get(gk, {}).get('total_mw', 0)
    sample.append(d)

# ── Summary stats function ──
def stats(vals):
    vals = [v for v in vals if v is not None]
    if not vals:
        return {}
    n = len(vals)
    vals.sort()
    mean = sum(vals) / n
    var = sum((v - mean)**2 for v in vals) / (n - 1) if n > 1 else 0
    std = math.sqrt(var)
    return {
        'N': n, 'Mean': mean, 'SD': std, 'Min': vals[0], 'Max': vals[-1],
        'P25': vals[n // 4], 'Median': vals[n // 2], 'P75': vals[3 * n // 4],
    }

# ── Dual output: stdout + lines list ──
lines = []
def out(s=''):
    print(s)
    lines.append(s)

# ── Analysis sample (complete theta) ──
analysis_sample = [d for d in sample if all(
    d.get(v) is not None for v in ['alpha', 'lambda', 'rho', 'kappa'])]

out('# Summary Statistics')
out()

# ── Key sample counts (traceable provenance for manuscript claims) ──
all_countries = set(d['fic'] for d in sample if d['fic'])
analysis_countries = set(d['fic'] for d in analysis_sample if d['fic'])
total_capacity_mw = sum(d['total_mw'] for d in analysis_sample)
total_capacity_tw = total_capacity_mw / 1e6

out('## Sample Overview')
out()
out(f'- Total Compustat utility firms: {len(sample)}')
out(f'- Analysis sample (complete theta): {len(analysis_sample)} firms')
out(f'- Countries (all): {len(all_countries)}')
out(f'- Countries (analysis sample): {len(analysis_countries)}')
out(f'- Total installed capacity (analysis sample): {total_capacity_mw:,.0f} MW = {total_capacity_tw:.1f} TW')
out()

# ── Panel A: Analysis sample ──
out(f'## Panel A: Analysis Sample (N = {len(analysis_sample)} firms, latest fiscal year)')
out()
fmt_hdr = '| {:<25} | {:>6} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} |'
fmt_row = '| {:<25} | {:>6} | {:>10} | {:>10} | {:>10} | {:>10} | {:>10} |'
out(fmt_hdr.format('Variable', 'N', 'Mean', 'SD', 'Min', 'Median', 'Max'))
out('|' + '---|' * 7)

for label, var in [
    ('Leverage (lambda)', 'lambda'),
    ('Return spread (rho)', 'rho'),
    ('Cash flow adequacy (kappa)', 'kappa'),
    ('Legacy intensity (alpha)', 'alpha'),
]:
    s = stats([d[var] for d in analysis_sample])
    if s:
        out(fmt_row.format(label, s['N'],
            f"{s['Mean']:.3f}", f"{s['SD']:.3f}",
            f"{s['Min']:.3f}", f"{s['Median']:.3f}", f"{s['Max']:.3f}"))
out()

# ── Panel B: GEM-matched subsample ──
gem_sample = [d for d in sample if d['has_gem']]
out(f'## Panel B: GEM-Matched Subsample (N = {len(gem_sample)} firms)')
out()
out(fmt_hdr.format('Variable', 'N', 'Mean', 'SD', 'Min', 'Median', 'Max'))
out('|' + '---|' * 7)

for label, var in [
    ('Legacy Intensity (alpha)', 'alpha'),
    ('Leverage (lambda)', 'lambda'),
    ('Operating ROA (rho)', 'rho'),
    ('Interest Coverage (kappa)', 'kappa'),
    ('Obligation Rigidity (delta)', 'delta'),
    ('Network Density', 'density'),
    ('Total Assets ($M)', 'at'),
]:
    s = stats([d[var] for d in gem_sample])
    if s:
        out(fmt_row.format(label, s['N'],
            f"{s['Mean']:.3f}", f"{s['SD']:.3f}",
            f"{s['Min']:.3f}", f"{s['Median']:.3f}", f"{s['Max']:.3f}"))
out()

# ── Panel C: Country distribution ──
out('## Panel C: Geographic Distribution')
out()
countries_all = Counter(d['fic'] for d in analysis_sample if d['fic'])
out(f'Countries in analysis sample: {len(countries_all)}')
out()
out('| Country | Firms |')
out('|---|---|')
for c, n in countries_all.most_common():
    out(f'| {c} | {n} |')
out()

# ── Panel D: Retirement events ──
out('## Panel D: Coal Retirement Events')
out()
events = []
evt_path = derived_path('events', 'coal_retirement_events.csv')
if os.path.exists(evt_path):
    with open(evt_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            events.append(row)

    total = len(events)
    matched = sum(1 for e in events if e.get('is_matched') == 'True')
    first_mover = sum(1 for e in events if e.get('is_first_mover') == 'True')
    fm_matched = sum(1 for e in events if e.get('is_first_mover') == 'True'
                     and e.get('is_matched') == 'True')
    total_mw = sum(float(e.get('capacity_mw', 0)) for e in events)
    matched_mw = sum(float(e.get('capacity_mw', 0)) for e in events
                     if e.get('is_matched') == 'True')

    # Country distribution of first-mover events
    fm_events = [e for e in events if e.get('is_first_mover') == 'True']
    fm_countries = Counter(e.get('country', '') for e in fm_events if e.get('country'))
    us_events = fm_countries.get('United States', 0) + fm_countries.get('US', 0) + fm_countries.get('USA', 0)

    out(f'- Total retirements: {total}')
    out(f'- Matched to Compustat: {matched}')
    out(f'- First-mover events: {first_mover}')
    out(f'- First-mover + matched: {fm_matched}')
    out(f'- Total retired MW: {total_mw:,.0f}')
    out(f'- Matched retired MW: {matched_mw:,.0f}')
    out(f'- Countries with first-mover events: {len(fm_countries)}')
    out(f'- US first-mover events: {us_events}')
    out()

    out('First-mover events by country:')
    out()
    out('| Country | Events |')
    out('|---|---|')
    for c, n in fm_countries.most_common():
        out(f'| {c} | {n} |')
else:
    out('(coal_retirement_events.csv not found)')
out()

# ── Panel E: EIA-860 and phase-out events ──
eia_path = derived_path('events', 'eia860_announcement_events.csv')
if os.path.exists(eia_path):
    eia_events = []
    with open(eia_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            eia_events.append(row)
    out(f'## Panel E: EIA-860 Announcement Events')
    out()
    out(f'- Total EIA-860 events: {len(eia_events)}')
    out()

phaseout_path = derived_path('events', 'coal_phaseout_shocks_events.csv')
if os.path.exists(phaseout_path):
    phaseout_events = []
    with open(phaseout_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            phaseout_events.append(row)
    out(f'## Panel F: Coal Phase-out Events')
    out()
    out(f'- Total phase-out events: {len(phaseout_events)}')
    out()

# ── Write to file ──
os.makedirs(results_path('summaries'), exist_ok=True)
outpath = results_path('summaries', 'summary_statistics.md')
with open(outpath, 'w', encoding='utf-8') as f:
    f.write('\n'.join(lines) + '\n')

print(f'\nSaved to {outpath}')
