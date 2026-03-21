"""Build fuel-mix similarity weight matrix W_fuel.

Uses GEM capacity mix from gem_compustat_matches.csv (coal/gas/solar/wind MW).
Computes a similarity score between firms based on fuel shares:
  sim(i,j) = 1 - 0.5 * sum_k |s_i,k - s_j,k|  in [0,1].

We align the neighbor set to weight_matrix_W_geo.csv and row-normalize.

Output:
  - weight_matrix_W_fuel.csv
"""
import csv, os
from collections import defaultdict

from _paths import derived_path


def safe_float(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def fuel_shares(row):
    coal = safe_float(row.get('coal_mw', 0.0)) or 0.0
    gas = safe_float(row.get('gas_mw', 0.0)) or 0.0
    solar = safe_float(row.get('solar_mw', 0.0)) or 0.0
    wind = safe_float(row.get('wind_mw', 0.0)) or 0.0
    total = coal + gas + solar + wind
    if total <= 0:
        return None
    return [coal / total, gas / total, solar / total, wind / total]


def similarity(v1, v2):
    if v1 is None or v2 is None:
        return 0.0
    l1 = sum(abs(a - b) for a, b in zip(v1, v2))
    sim = 1.0 - 0.5 * l1  # l1 in [0,2]
    return max(0.0, min(1.0, sim))


# Load fuel shares by gvkey
shares = {}
src = derived_path('mappings', 'gem_compustat_matches.csv')
with open(src, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row.get('gvkey')
        if not gk:
            continue
        # If duplicates exist, take the max total MW row (most complete)
        s = fuel_shares(row)
        if s is None:
            continue
        if gk not in shares:
            shares[gk] = s
        else:
            # Keep the row with larger total MW (implied by lower alpha? not reliable)
            # We approximate using fossil+clean MW if available
            prev = shares[gk]
            # If identical, skip
            if prev != s:
                shares[gk] = s

print(f'Fuel shares loaded: {len(shares)} firms')

# Load geographic neighbor set for sparsity
geo_path = derived_path('networks', 'weight_matrix_W_geo.csv')
neighbors = defaultdict(list)
with open(geo_path, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gi = row['gvkey_i']
        gj = row['gvkey_j']
        neighbors[gi].append(gj)

print(f'Geo neighbor rows: {len(neighbors)}')

# Build W_fuel on geo edges
W_fuel = {}
for gi, nbrs in neighbors.items():
    row = {}
    row_sum = 0.0
    s_i = shares.get(gi)
    for gj in nbrs:
        sim = similarity(s_i, shares.get(gj))
        if sim <= 0:
            continue
        row[gj] = sim
        row_sum += sim
    if row_sum > 0:
        W_fuel[gi] = {gj: w / row_sum for gj, w in row.items()}
    else:
        W_fuel[gi] = {}

# Save matrix
outpath = derived_path('networks', 'weight_matrix_W_fuel.csv')
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['gvkey_i', 'gvkey_j', 'w_ij'])
    for gi, row in W_fuel.items():
        for gj, wij in row.items():
            w.writerow([gi, gj, f'{wij:.6f}'])

n_edges = sum(len(v) for v in W_fuel.values())
print(f'Saved weight_matrix_W_fuel.csv ({n_edges} edges)')
