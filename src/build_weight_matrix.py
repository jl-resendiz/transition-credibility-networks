"""Build spatial weight matrix W for matched Compustat firms.

Geographic component:
  - Compute centroid of GEM plants (MW-weighted).
  - Weight w_ij = exp(-d_ij / DECAY_KM) / d_ij (no hard cutoff).
Regulatory component (optional):
  - If weight_matrix_W_regulatory.csv exists, add REG_WEIGHT * w_reg_ij.
Row-normalize so each row sums to 1.
Outputs:
  - firm_centroids.csv
  - weight_matrix_W_geo.csv (geographic only)
  - weight_matrix_W.csv (composite)
"""
import csv, os, math, re
from collections import defaultdict
from _paths import derived_path

def parse_parents(field):
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

def haversine(lat1, lon1, lat2, lon2):
    """Distance in km between two GPS points."""
    R = 6371
    dlat = math.radians(lat2 - lat1)
    dlon = math.radians(lon2 - lon1)
    a = math.sin(dlat/2)**2 + math.cos(math.radians(lat1)) * math.cos(math.radians(lat2)) * math.sin(dlon/2)**2
    return R * 2 * math.asin(math.sqrt(a))

# 1. Load GEM->Compustat matches to know which GEM parents map to which gvkeys
parent_to_gvkeys = defaultdict(set)
with open(derived_path('mappings', 'gem_compustat_matches.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        parent_to_gvkeys[row['gem_parent']].add(row['gvkey'])

print(f'Matched parents: {len(parent_to_gvkeys)}')

# 2. Read plant-level GPS from GEM trackers, aggregate to firm centroid
gvkey_plants = defaultdict(list)  # gvkey -> [(lat, lon, mw), ...]

trackers = [
    ('gem_coal.csv', 'Parent'),
    ('gem_gas.csv', 'Parent(s)'),
    ('gem_solar.csv', 'Owner'),
    ('gem_wind.csv', 'Owner'),
]

for fname, parent_col in trackers:
    fpath = derived_path('gem', fname)
    print(f'Reading {fname}...')
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            status = row.get('Status', '')
            if status != 'operating':
                continue
            try:
                cap = float(row['Capacity (MW)'])
                lat = float(row['Latitude'])
                lon = float(row['Longitude'])
            except (ValueError, TypeError, KeyError):
                continue

            parsed = parse_parents(row.get(parent_col, ''))
            for name, pct in parsed:
                if name in parent_to_gvkeys:
                    share = (pct / 100.0) if pct else 1.0 / len(parsed) if len(parsed) > 1 else 1.0
                    for gvkey in parent_to_gvkeys[name]:
                        gvkey_plants[gvkey].append((lat, lon, cap * share))

print(f'\nFirms with plant GPS data: {len(gvkey_plants)}')

# 3. Compute MW-weighted centroid per firm
centroids = {}
for gvkey, plants in gvkey_plants.items():
    total_mw = sum(mw for _, _, mw in plants)
    if total_mw <= 0:
        continue
    wlat = sum(lat * mw for lat, lon, mw in plants) / total_mw
    wlon = sum(lon * mw for lat, lon, mw in plants) / total_mw
    centroids[gvkey] = (wlat, wlon, total_mw, len(plants))

print(f'Firms with valid centroids: {len(centroids)}')

# 4. Build distance matrix and geographic weight matrix
gvkeys = sorted(centroids.keys())
n = len(gvkeys)
# Exponential decay scale. We set the half-life to 500 km to avoid a hard cutoff.
DECAY_KM = 500 / math.log(2)
REG_WEIGHT = 1.0    # weight on regulatory component (if available)

# Compute pairwise distances and inverse-distance weights
W_geo = {}
row_sums = {}
row_avg_dist = {}
for i, gi in enumerate(gvkeys):
    lat_i, lon_i = centroids[gi][:2]
    row_sum = 0
    neighbors = {}
    dist_weight_sum = 0.0
    for j, gj in enumerate(gvkeys):
        if i == j:
            continue
        lat_j, lon_j = centroids[gj][:2]
        d = haversine(lat_i, lon_i, lat_j, lon_j)
        if d > 0:
            w = math.exp(-d / DECAY_KM) / d
            neighbors[gj] = w
            row_sum += w
            dist_weight_sum += d * w
    # Row-normalize
    if row_sum > 0:
        W_geo[gi] = {gj: w / row_sum for gj, w in neighbors.items()}
    else:
        W_geo[gi] = {}
    row_sums[gi] = row_sum
    row_avg_dist[gi] = (dist_weight_sum / row_sum) if row_sum > 0 else None

# 4b. Load regulatory component (if available) and build composite W
reg_path = derived_path('networks', 'weight_matrix_W_regulatory.csv')
reg_edges = defaultdict(dict)
if os.path.exists(reg_path):
    with open(reg_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gi = row['gvkey_i']
            gj = row['gvkey_j']
            try:
                reg_edges[gi][gj] = float(row.get('w_reg', 1.0))
            except (ValueError, TypeError):
                reg_edges[gi][gj] = 1.0
    print(f'Loaded regulatory edges: {sum(len(v) for v in reg_edges.values())}')
else:
    print('No regulatory matrix found; using geographic component only.')

W = {}
for gi in gvkeys:
    combined = {}
    for gj, w in W_geo.get(gi, {}).items():
        combined[gj] = combined.get(gj, 0.0) + w
    for gj, w in reg_edges.get(gi, {}).items():
        combined[gj] = combined.get(gj, 0.0) + REG_WEIGHT * w
    row_sum = sum(combined.values())
    if row_sum > 0:
        W[gi] = {gj: w / row_sum for gj, w in combined.items()}
    else:
        W[gi] = {}

# Stats
n_connected = sum(1 for gi in gvkeys if len(W.get(gi, {})) > 0)
avg_neighbors = sum(len(W.get(gi, {})) for gi in gvkeys) / n if n > 0 else 0
max_neighbors = max(len(W.get(gi, {})) for gi in gvkeys) if gvkeys else 0

print(f'\n=== WEIGHT MATRIX W ===')
print(f'Firms: {n}')
print(f'Decay scale: {DECAY_KM:.1f} km (half-life = 500 km)')
print(f'Regulatory weight: {REG_WEIGHT}')
print(f'Connected firms: {n_connected} ({100*n_connected/n:.1f}%)')
print(f'Avg neighbors: {avg_neighbors:.1f}')
print(f'Max neighbors: {max_neighbors}')

# Save centroid data
outpath = derived_path('networks', 'firm_centroids.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gvkey', 'centroid_lat', 'centroid_lon', 'total_mw', 'n_plants', 'n_neighbors', 'w_sum', 'avg_dist_km'])
    for gk in gvkeys:
        lat, lon, mw, nplants = centroids[gk]
        nn = len(W.get(gk, {}))
        wsum = row_sums.get(gk, 0.0)
        avgd = row_avg_dist.get(gk)
        avgd_str = f'{avgd:.2f}' if avgd is not None else ''
        w.writerow([gk, f'{lat:.4f}', f'{lon:.4f}', f'{mw:.1f}', nplants, nn, f'{wsum:.6f}', avgd_str])
print(f'Saved firm_centroids.csv')

# Save geographic W matrix (only non-zero entries)
outpath = derived_path('networks', 'weight_matrix_W_geo.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gvkey_i', 'gvkey_j', 'w_ij'])
    for gi in gvkeys:
        for gj, wij in sorted(W_geo.get(gi, {}).items()):
            w.writerow([gi, gj, f'{wij:.6f}'])
n_edges_geo = sum(len(v) for v in W_geo.values())
print(f'Saved weight_matrix_W_geo.csv ({n_edges_geo} edges)')

# Save composite W matrix (only non-zero entries)
outpath = derived_path('networks', 'weight_matrix_W.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gvkey_i', 'gvkey_j', 'w_ij'])
    for gi in gvkeys:
        for gj, wij in sorted(W.get(gi, {}).items()):
            w.writerow([gi, gj, f'{wij:.6f}'])
n_edges = sum(len(v) for v in W.values())
print(f'Saved weight_matrix_W.csv ({n_edges} edges)')
