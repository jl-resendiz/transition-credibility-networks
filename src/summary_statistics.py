"""Generate summary statistics table for the paper."""
import csv, os, math
from collections import defaultdict

from _paths import derived_path

# Load firm fundamentals (latest year per firm, with complete theta)
firms = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in firms or fy > firms[gk]['fyear']:
            firms[gk] = row

# Load density (raw weight sum preferred)
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

# Build sample
sample = []
for gk, f in firms.items():
    d = {}
    d['gvkey'] = gk
    d['conm'] = f['conm']
    d['fic'] = f['fic']

    # Financials
    for var in ['at', 'sale', 'ppent']:
        try:
            d[var] = float(f[var]) if f.get(var) else None
        except (ValueError, TypeError):
            d[var] = None

    # Theta
    for var in ['alpha', 'lambda', 'rho', 'kappa', 'delta']:
        try:
            d[var] = float(f[var]) if f.get(var) else None
        except (ValueError, TypeError):
            d[var] = None

    d['density'] = density.get(gk)
    d['has_gem'] = d['alpha'] is not None
    sample.append(d)

# Summary stats function
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
        'N': n,
        'Mean': mean,
        'SD': std,
        'P25': vals[n // 4],
        'Median': vals[n // 2],
        'P75': vals[3 * n // 4],
    }

# Panel A: Full sample
print('=' * 80)
print('TABLE 1: SUMMARY STATISTICS')
print('=' * 80)

print(f'\nPanel A: Full Compustat Sample (N = {len(sample)} firms, latest fiscal year)')
fmt = '{:<25} {:>6} {:>10} {:>10} {:>10} {:>10} {:>10}'
print(fmt.format('Variable', 'N', 'Mean', 'SD', 'P25', 'Median', 'P75'))
print('-' * 85)

variables = [
    ('Total Assets ($M)', 'at'),
    ('Revenue ($M)', 'sale'),
    ('PP&E ($M)', 'ppent'),
    ('Leverage (lambda)', 'lambda'),
    ('Operating ROA (rho)', 'rho'),
    ('Interest Coverage (kappa)', 'kappa'),
    ('Obligation Rigidity (delta)', 'delta'),
]

for label, var in variables:
    s = stats([d[var] for d in sample])
    if s:
        print(fmt.format(label, s['N'], f"{s['Mean']:,.2f}", f"{s['SD']:,.2f}",
                          f"{s['P25']:,.2f}", f"{s['Median']:,.2f}", f"{s['P75']:,.2f}"))

# Panel B: GEM-matched subsample
gem_sample = [d for d in sample if d['has_gem']]
print(f'\nPanel B: GEM-Matched Subsample (N = {len(gem_sample)} firms)')
print(fmt.format('Variable', 'N', 'Mean', 'SD', 'P25', 'Median', 'P75'))
print('-' * 85)

variables_gem = [
    ('Legacy Intensity (alpha)', 'alpha'),
    ('Delivery (1-alpha)', None),
    ('Leverage (lambda)', 'lambda'),
    ('Operating ROA (rho)', 'rho'),
    ('Interest Coverage (kappa)', 'kappa'),
    ('Obligation Rigidity (delta)', 'delta'),
    ('Network Density', 'density'),
    ('Total Assets ($M)', 'at'),
]

for label, var in variables_gem:
    if var is None:
        # Delivery
        s = stats([1.0 - d['alpha'] for d in gem_sample if d['alpha'] is not None])
    else:
        s = stats([d[var] for d in gem_sample])
    if s:
        print(fmt.format(label, s['N'], f"{s['Mean']:,.3f}", f"{s['SD']:,.3f}",
                          f"{s['P25']:,.3f}", f"{s['Median']:,.3f}", f"{s['P75']:,.3f}"))

# Panel C: Country distribution
print(f'\nPanel C: Geographic Distribution')
from collections import Counter
countries = Counter(d['fic'] for d in sample if d['fic'])
gem_countries = Counter(d['fic'] for d in gem_sample if d['fic'])
print(f'{"Country":<8} {"All":>6} {"GEM-matched":>12}')
print('-' * 30)
for c, n in countries.most_common(20):
    gn = gem_countries.get(c, 0)
    print(f'{c:<8} {n:>6} {gn:>12}')

# Panel D: Retirement events
print(f'\nPanel D: Coal Retirement Events (2015-2025)')
events = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        events.append(row)

total = len(events)
matched = sum(1 for e in events if e['is_matched'] == 'True')
first_mover = sum(1 for e in events if e['is_first_mover'] == 'True')
fm_matched = sum(1 for e in events if e['is_first_mover'] == 'True' and e['is_matched'] == 'True')
print(f'  Total retirements: {total}')
print(f'  Matched to Compustat: {matched}')
print(f'  First-mover events: {first_mover}')
print(f'  First-mover + matched: {fm_matched}')
total_mw = sum(float(e['capacity_mw']) for e in events)
matched_mw = sum(float(e['capacity_mw']) for e in events if e['is_matched'] == 'True')
print(f'  Total retired MW: {total_mw:,.0f}')
print(f'  Matched retired MW: {matched_mw:,.0f}')
