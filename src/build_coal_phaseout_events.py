"""Build derived coal phase-out shock events with matched gvkeys.

Outputs:
  - derived/events/coal_phaseout_shocks_events.csv
  - derived/events/coal_phaseout_shocks_state_mapping.csv (US state coal exposure)
"""
import csv
import os
import re
from collections import defaultdict

from _paths import raw_path, derived_path

RAW_SHOCKS = raw_path('policy', 'coal_phaseout_shocks.csv')
OUT_EVENTS = derived_path('events', 'coal_phaseout_shocks_events.csv')
OUT_STATE = derived_path('events', 'coal_phaseout_shocks_state_mapping.csv')

GEM_COAL = derived_path('gem', 'gem_coal.csv')
GEM_MATCH = derived_path('mappings', 'gem_compustat_matches.csv')


EU_ISO3 = {
    'AUT','BEL','BGR','HRV','CYP','CZE','DNK','EST','FIN','FRA','DEU','GRC','HUN','IRL','ITA','LVA','LTU','LUX',
    'MLT','NLD','POL','PRT','ROU','SVK','SVN','ESP','SWE'
}

STATE_ALIASES = {
    'CA': 'California',
    'WA': 'Washington',
    'OR': 'Oregon',
    'NM': 'New Mexico',
    'CO': 'Colorado',
    'NY': 'New York',
    'VA': 'Virginia',
    'IL': 'Illinois',
    'NC': 'North Carolina',
    'MN': 'Minnesota',
}


def normalize_name(s):
    if s is None:
        return ''
    s = s.lower().strip()
    s = re.sub(r'[^a-z0-9 ]+', ' ', s)
    s = re.sub(r'\s+', ' ', s).strip()
    return s


def parse_parents(field):
    if not field:
        return []
    parts = str(field).split(';')
    out = []
    for p in parts:
        p = p.strip()
        m = re.match(r'^(.+?)\s*\[(\d+\.?\d*)%\]$', p)
        if m:
            out.append((m.group(1).strip(), float(m.group(2))))
        elif p:
            out.append((p.strip(), None))
    return out


# Load gvkey -> fic (latest)
fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

gvkey_fic = {gk: row.get('fic') for gk, row in fundamentals.items()}

# Load GEM parent -> gvkeys map
parent_to_gvkeys = defaultdict(set)
with open(GEM_MATCH, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gp = normalize_name(row.get('gem_parent', ''))
        gk = row.get('gvkey')
        if gp and gk:
            parent_to_gvkeys[gp].add(gk)


# Build gvkey -> state coal MW (US only)
gvkey_state_mw = defaultdict(lambda: defaultdict(float))
if os.path.exists(GEM_COAL):
    with open(GEM_COAL, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        headers = reader.fieldnames
        needed = ['Country/Area', 'Subnational unit (province, state)', 'Parent', 'Capacity (MW)', 'Status']
        if all(k in headers for k in needed):
            for row in reader:
                country = row['Country/Area']
                state = row['Subnational unit (province, state)']
                parent = row['Parent']
                status = row['Status']
                cap = row['Capacity (MW)']
                if not country or country.strip() != 'United States':
                    continue
                if not state:
                    continue
                if status and status.strip().lower() != 'operating':
                    continue
                try:
                    cap = float(cap)
                except (ValueError, TypeError):
                    continue
                parents = parse_parents(parent)
                if not parents:
                    continue
                for name, pct in parents:
                    share = (pct / 100.0) if pct else 1.0 / len(parents)
                    gp = normalize_name(name)
                    for gk in parent_to_gvkeys.get(gp, []):
                        gvkey_state_mw[gk][state.strip()] += cap * share

# Write state mapping output
os.makedirs(os.path.dirname(OUT_STATE), exist_ok=True)
with open(OUT_STATE, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gvkey', 'state', 'coal_mw'])
    for gk, states in gvkey_state_mw.items():
        for st, mw in states.items():
            w.writerow([gk, st, f'{mw:.2f}'])


def parse_affected_countries(val):
    if not val:
        return []
    v = str(val).strip()
    if v.lower().startswith('eu'):
        return list(EU_ISO3)
    if 'EU ETS' in v or 'EU-wide' in v or 'EU-wide' in v:
        return list(EU_ISO3)
    # split by comma or semicolon
    parts = re.split(r'[;,]', v)
    out = []
    for p in parts:
        p = p.strip()
        if not p:
            continue
        out.append(p)
    return out


# Build events
rows = []
with open(RAW_SHOCKS, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        event_date = row.get('event_date', '').strip()
        if not event_date:
            continue
        tier = row.get('exogeneity_tier', '').strip()
        binding = row.get('binding', '').strip().lower()
        jurisdiction = row.get('jurisdiction', '').strip()
        affected = row.get('affected_countries', '').strip()

        treated = set()
        # State-level US laws (jurisdiction name)
        if row.get('shock_id', '').startswith('US_') and jurisdiction:
            state_name = jurisdiction
            # allow abbreviations
            if jurisdiction in STATE_ALIASES:
                state_name = STATE_ALIASES[jurisdiction]
            for gk, states in gvkey_state_mw.items():
                if state_name in states and states[state_name] > 0:
                    treated.add(gk)
        else:
            countries = parse_affected_countries(affected if affected else jurisdiction)
            for gk, fic in gvkey_fic.items():
                if fic in countries:
                    treated.add(gk)

        if not treated:
            continue

        rows.append({
            'shock_id': row.get('shock_id'),
            'shock_name': row.get('shock_name'),
            'jurisdiction': jurisdiction,
            'event_date': event_date,
            'legal_date': row.get('legal_date', ''),
            'instrument_type': row.get('instrument_type', ''),
            'binding': row.get('binding', ''),
            'exogeneity_tier': tier,
            'affected_countries': affected,
            'matched_gvkeys': ';'.join(sorted(treated)),
            'is_first_mover': 'True' if tier == '1' else 'False',
            'source': row.get('source', ''),
        })

os.makedirs(os.path.dirname(OUT_EVENTS), exist_ok=True)
with open(OUT_EVENTS, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys()) if rows else [])
    if rows:
        w.writeheader()
        w.writerows(rows)

print(f'Wrote {len(rows)} events to {OUT_EVENTS}')
print(f'State mapping rows: {sum(len(s) for s in gvkey_state_mw.values())}')
