"""Parse all 4 GEM trackers, compute alpha (fossil intensity) per parent entity."""
import re, csv
from collections import defaultdict
from _paths import derived_path

def parse_parents(field):
    """Parse 'AES Corp [100.0%]; Enel SpA [50.0%]' into [(name, pct), ...]"""
    if not field or str(field).strip() == '':
        return []
    parts = str(field).split(';')
    results = []
    for p in parts:
        p = p.strip()
        match = re.match(r'^(.+?)\s*\[(\d+\.?\d*)%\]$', p)
        if match:
            results.append((match.group(1).strip(), float(match.group(2))))
        elif p:
            results.append((p.strip(), None))
    return results

all_parents = defaultdict(lambda: {'coal_mw': 0, 'gas_mw': 0, 'solar_mw': 0, 'wind_mw': 0, 'units': 0})

trackers = [
    ('gem_coal.csv', 'Parent', 'coal_mw'),
    ('gem_gas.csv', 'Parent(s)', 'gas_mw'),
    ('gem_solar.csv', 'Owner', 'solar_mw'),
    ('gem_wind.csv', 'Owner', 'wind_mw'),
]

for fname, parent_col, mw_key in trackers:
    fpath = derived_path('gem', fname)
    print(f'Processing {fname}...')
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            status = row.get('Status', '')
            if status != 'operating':
                continue
            cap = row.get('Capacity (MW)', '')
            if cap == '':
                continue
            try:
                cap = float(cap)
            except (ValueError, TypeError):
                continue
            parsed = parse_parents(row.get(parent_col, ''))
            if not parsed:
                continue
            for name, pct in parsed:
                share = (pct / 100.0) if pct else 1.0 / len(parsed) if len(parsed) > 1 else 1.0
                all_parents[name][mw_key] += cap * share
                all_parents[name]['units'] += 1

# Compute totals
for p in all_parents:
    all_parents[p]['fossil_mw'] = all_parents[p]['coal_mw'] + all_parents[p]['gas_mw']
    all_parents[p]['clean_mw'] = all_parents[p]['solar_mw'] + all_parents[p]['wind_mw']
    all_parents[p]['total_mw'] = all_parents[p]['fossil_mw'] + all_parents[p]['clean_mw']

sorted_parents = sorted(all_parents.items(), key=lambda x: x[1]['total_mw'], reverse=True)

print(f'\nTotal unique parent/owner entities: {len(all_parents)}')
print(f'\nTop 30 by total MW:')
fmt = '{:<45} {:>8} {:>8} {:>8} {:>8} {:>8} {:>6}'
print(fmt.format('Name', 'Coal', 'Gas', 'Solar', 'Wind', 'Total', 'Alpha'))
for name, d in sorted_parents[:30]:
    alpha = d['fossil_mw'] / d['total_mw'] if d['total_mw'] > 0 else 0
    print(fmt.format(name[:45], int(d['coal_mw']), int(d['gas_mw']), int(d['solar_mw']), int(d['wind_mw']), int(d['total_mw']), f'{alpha:.2f}'))

# Save full parent list
outpath = derived_path('mappings', 'gem_parents_parsed.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['parent_name', 'coal_mw', 'gas_mw', 'solar_mw', 'wind_mw', 'fossil_mw', 'clean_mw', 'total_mw', 'alpha', 'units'])
    for name, d in sorted_parents:
        alpha = d['fossil_mw'] / d['total_mw'] if d['total_mw'] > 0 else 0
        w.writerow([name, f"{d['coal_mw']:.1f}", f"{d['gas_mw']:.1f}", f"{d['solar_mw']:.1f}", f"{d['wind_mw']:.1f}",
                     f"{d['fossil_mw']:.1f}", f"{d['clean_mw']:.1f}", f"{d['total_mw']:.1f}", f'{alpha:.4f}', d['units']])

print(f'\nSaved {outpath} ({len(all_parents)} entities)')
