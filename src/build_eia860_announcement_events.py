"""Build EIA-860-style announcement events from existing retirement announcements.

This is a proxy builder that filters to US firms (fic=USA) and retains
only events with announcement dates. Output is compatible with
panel_did.py via EVENTS_PATH.
"""
import csv
from _paths import derived_path

# Input: derived retirement events
in_path = derived_path('events', 'coal_retirement_events.csv')
# Output: announcement event file
out_path = derived_path('events', 'eia860_announcement_events.csv')

# Load latest firm country (fic)
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

rows = []
with open(in_path, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if not row.get('matched_gvkeys'):
            continue
        ann_date = row.get('announcement_date', '').strip()
        if not ann_date:
            continue
        gvkeys = [g.strip() for g in row['matched_gvkeys'].split(';') if g.strip()]
        # Keep US only to mimic EIA 860 scope
        us_gvkeys = []
        for gk in gvkeys:
            fic = fundamentals.get(gk, {}).get('fic')
            if fic == 'USA':
                us_gvkeys.append(gk)
        if not us_gvkeys:
            continue
        rows.append({
            'event_id': row.get('gem_id'),
            'announcement_date': ann_date,
            'announcement_source': row.get('announcement_source', ''),
            'event_date': row.get('event_date', ''),
            'ret_year': row.get('ret_year', ''),
            'matched_gvkeys': ';'.join(us_gvkeys),
            'is_first_mover': row.get('is_first_mover', ''),
        })

with open(out_path, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
    if rows:
        w.writeheader()
        w.writerows(rows)

print(f'Wrote {len(rows)} events to {out_path}')
