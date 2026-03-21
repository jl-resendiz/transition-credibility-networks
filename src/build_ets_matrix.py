"""Build ETS regulatory component of weight matrix W.

Uses World Bank Carbon Pricing Dashboard to identify which firms share an ETS.
W_reg is row-normalized: w_ij = 1/deg(i) if firms i and j share an ETS.
"""
import openpyxl, csv, os
from collections import defaultdict

from _paths import raw_path, derived_path

# 1. Extract operating ETS jurisdictions from World Bank data
wb = openpyxl.load_workbook(raw_path('policy', 'carbon_pricing_worldbank.xlsx'), read_only=True)
ws = wb['Compliance_Gen Info']
rows = list(ws.iter_rows(values_only=True))

ets_instruments = []
for row in rows[5:]:
    if row[0] is None:
        continue
    if row[2] == 'ETS' and row[3] == 'Implemented':
        ets_instruments.append({
            'uid': row[0],
            'name': row[1],
            'jurisdiction': row[4],
        })
wb.close()

print(f'Operating ETS instruments: {len(ets_instruments)}')

# 2. Map ETS jurisdictions to ISO country codes (matching Compustat fic codes)
# Multi-country systems map to all member countries
ETS_TO_COUNTRIES = {
    # EU ETS covers EU27 + EEA
    'EU ETS': ['AUT', 'BEL', 'BGR', 'HRV', 'CYP', 'CZE', 'DNK', 'EST', 'FIN', 'FRA',
               'DEU', 'GRC', 'HUN', 'IRL', 'ITA', 'LVA', 'LTU', 'LUX', 'MLT', 'NLD',
               'POL', 'PRT', 'ROU', 'SVK', 'SVN', 'ESP', 'SWE', 'ISL', 'LIE', 'NOR'],
    'UK ETS': ['GBR'],
    'Switzerland ETS': ['CHE'],
    'China national ETS': ['CHN'],
    'Korea ETS': ['KOR'],
    'New Zealand ETS': ['NZL'],
    'Kazakhstan ETS': ['KAZ'],
    'Indonesia ETS': ['IDN'],
    'Montenegro ETS': ['MNE'],
    'Mexico ETS': ['MEX'],
    'Australia Safeguard Mechanism': ['AUS'],
    # Subnational systems map to their country
    'California CaT': ['USA'],
    'Quebec CaT': ['CAN'],
    'RGGI': ['USA'],
    'Alberta TIER': ['CAN'],
    'BC OBPS': ['CAN'],
    'Canada federal OBPS': ['CAN'],
    'Massachusetts ETS': ['USA'],
    'Oregon ETS': ['USA'],
    'Washington CCA': ['USA'],
    'Colorado GHG crediting trading system': ['USA'],
    'Ontario EPS': ['CAN'],
    'New Brunswick OBPS': ['CAN'],
    'Newfoundland and Labrador PSS': ['CAN'],
    'Nova Scotia OBPS': ['CAN'],
    'Germany ETS': ['DEU'],  # national fuel ETS, separate from EU ETS
    'Austria ETS': ['AUT'],  # national fuel ETS, separate from EU ETS
    # Chinese pilots
    'Beijing pilot ETS': ['CHN'],
    'Shanghai pilot ETS': ['CHN'],
    'Shenzhen pilot ETS': ['CHN'],
    'Guangdong pilot ETS': ['CHN'],
    'Hubei pilot ETS': ['CHN'],
    'Chongqing pilot ETS': ['CHN'],
    'Fujian pilot ETS': ['CHN'],
    'Tianjin pilot ETS': ['CHN'],
    # Japanese subnational
    'Tokyo CaT': ['JPN'],
    'Saitama ETS': ['JPN'],
}

# Compustat uses 3-letter ISO but some use 2-letter. Map to Compustat fic codes.
ISO3_TO_FIC = {
    'AUT': 'AUT', 'BEL': 'BEL', 'BGR': 'BGR', 'HRV': 'HRV', 'CYP': 'CYP',
    'CZE': 'CZE', 'DNK': 'DNK', 'EST': 'EST', 'FIN': 'FIN', 'FRA': 'FRA',
    'DEU': 'DEU', 'GRC': 'GRC', 'HUN': 'HUN', 'IRL': 'IRL', 'ITA': 'ITA',
    'LVA': 'LVA', 'LTU': 'LTU', 'LUX': 'LUX', 'MLT': 'MLT', 'NLD': 'NLD',
    'POL': 'POL', 'PRT': 'PRT', 'ROU': 'ROU', 'SVK': 'SVK', 'SVN': 'SVN',
    'ESP': 'ESP', 'SWE': 'SWE', 'ISL': 'ISL', 'LIE': 'LIE', 'NOR': 'NOR',
    'GBR': 'GBR', 'CHE': 'CHE', 'CHN': 'CHN', 'KOR': 'KOR', 'NZL': 'NZL',
    'KAZ': 'KAZ', 'IDN': 'IDN', 'MNE': 'MNE', 'MEX': 'MEX', 'AUS': 'AUS',
    'USA': 'USA', 'CAN': 'CAN', 'JPN': 'JPN',
}

# Build: country_fic -> set of ETS names it belongs to
country_ets = defaultdict(set)
for ets in ets_instruments:
    name = ets['name']
    countries = ETS_TO_COUNTRIES.get(name, [])
    for c in countries:
        fic = ISO3_TO_FIC.get(c, c)
        country_ets[fic].add(name)

print(f'Countries with ETS coverage: {len(country_ets)}')

# 3. Load firm fundamentals to get fic per gvkey
firm_fic = {}
for src in ['compustat_global_utilities.csv', 'compustat_na_utilities.csv']:
    fpath = raw_path('compustat', src)
    if not os.path.exists(fpath):
        continue
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            fy = row['fyear']
            if gk not in firm_fic or fy > firm_fic[gk][1]:
                firm_fic[gk] = (row['fic'], fy)

firm_fic = {gk: v[0] for gk, v in firm_fic.items()}

# 4. Load GEM-Compustat matches (firms with plant data)
gvkeys_in_w = set()
with open(derived_path('mappings', 'gem_compustat_matches.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gvkeys_in_w.add(row['gvkey'])

gvkeys_in_w = sorted(gvkeys_in_w)
print(f'Firms with GEM matches: {len(gvkeys_in_w)}')

# 5. Build regulatory W: same_ets(i,j) = 1 if they share at least one ETS
def shares_ets(gk1, gk2):
    fic1 = firm_fic.get(gk1, '')
    fic2 = firm_fic.get(gk2, '')
    ets1 = country_ets.get(fic1, set())
    ets2 = country_ets.get(fic2, set())
    return len(ets1 & ets2) > 0

# Count firms under ETS
n_ets_firms = sum(1 for gk in gvkeys_in_w if firm_fic.get(gk, '') in country_ets)
print(f'Firms under an ETS: {n_ets_firms} / {len(gvkeys_in_w)}')

# Build sparse regulatory adjacency
reg_adj = defaultdict(set)
for i, gi in enumerate(gvkeys_in_w):
    for j, gj in enumerate(gvkeys_in_w):
        if i >= j:
            continue
        if shares_ets(gi, gj):
            reg_adj[gi].add(gj)
            reg_adj[gj].add(gi)

undirected_edges = sum(len(v) for v in reg_adj.values()) // 2
print(f'Regulatory edges (same ETS pairs): {undirected_edges}')

# 6. Save regulatory component (row-normalized)
outpath = derived_path('networks', 'weight_matrix_W_regulatory.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gvkey_i', 'gvkey_j', 'w_reg'])
    directed_edges = 0
    for gi, nbrs in reg_adj.items():
        deg = len(nbrs)
        if deg == 0:
            continue
        w_ij = 1.0 / deg
        for gj in nbrs:
            w.writerow([gi, gj, f'{w_ij:.6f}'])
            directed_edges += 1
print(f'Saved {outpath} ({directed_edges} directed edges)')

# 7. Save firm ETS membership
outpath2 = derived_path('networks', 'firm_ets_membership.csv')
with open(outpath2, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gvkey', 'fic', 'ets_names', 'has_ets'])
    for gk in gvkeys_in_w:
        fic = firm_fic.get(gk, '')
        ets_names = ';'.join(sorted(country_ets.get(fic, set())))
        has_ets = 1 if ets_names else 0
        w.writerow([gk, fic, ets_names, has_ets])
print(f'Saved {outpath2}')

# Summary by ETS
from collections import Counter
ets_firm_count = Counter()
for gk in gvkeys_in_w:
    fic = firm_fic.get(gk, '')
    for ets_name in country_ets.get(fic, set()):
        ets_firm_count[ets_name] += 1

print(f'\nFirms per ETS (top 10):')
for ets_name, count in ets_firm_count.most_common(10):
    print(f'  {ets_name:<40} {count} firms')
