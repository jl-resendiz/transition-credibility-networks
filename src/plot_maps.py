"""Generate publication-quality maps for the paper.

Figure 1: Global sample — firm centroids colored by alpha (fossil intensity)
Figure 2: Network structure — W edges with density-colored nodes
Figure 3: Coal retirement first-mover events
"""
import csv, os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import numpy as np

from _paths import derived_path, DATA_ROOT
OUT = os.path.join(os.path.dirname(DATA_ROOT), 'finance_draft', 'figures')
os.makedirs(OUT, exist_ok=True)

# Style
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 9,
    'axes.linewidth': 0.5,
    'figure.dpi': 300,
})

# Load data
centroids = {}
with open(derived_path('networks', 'firm_centroids.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        wsum = row.get('w_sum', '')
        if wsum not in ('', None):
            try:
                dens = float(wsum)
            except (ValueError, TypeError):
                dens = float(row.get('n_neighbors', 0))
        else:
            dens = float(row.get('n_neighbors', 0))
        centroids[row['gvkey']] = {
            'lat': float(row['centroid_lat']),
            'lon': float(row['centroid_lon']),
            'mw': float(row['total_mw']),
            'n_neighbors': int(row['n_neighbors']),
            'density': dens,
        }

fundamentals = {}
with open(derived_path('fundamentals', 'firm_fundamentals.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        fy = row['fyear']
        if gk not in fundamentals or fy > fundamentals[gk]['fyear']:
            fundamentals[gk] = row

W_edges = []
with open(derived_path('networks', 'weight_matrix_W.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gi, gj = row['gvkey_i'], row['gvkey_j']
        if gi < gj:  # avoid duplicates
            W_edges.append((gi, gj, float(row['w_ij'])))

retirements = []
with open(derived_path('events', 'coal_retirement_events.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        if row['is_first_mover'] == 'True':
            retirements.append({
                'lat': float(row['lat']),
                'lon': float(row['lon']),
                'mw': float(row['capacity_mw']),
                'year': int(row['ret_year']),
                'matched': row['is_matched'] == 'True',
                'country': row['country'],
            })


# =========================================================================
# FIGURE 1: Global sample with alpha gradient
# =========================================================================
print('Generating Figure 1: Global sample with alpha...')
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.add_feature(cfeature.LAND, facecolor='#f0f0f0', edgecolor='none')
ax.add_feature(cfeature.BORDERS, linewidth=0.3, edgecolor='#cccccc')
ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#999999')

# Colormap: green (clean) to red (fossil)
cmap = plt.cm.RdYlGn_r
norm = mcolors.Normalize(vmin=0, vmax=1)

lats, lons, alphas, sizes = [], [], [], []
for gk, c in centroids.items():
    f = fundamentals.get(gk)
    if f and f.get('alpha', '') != '':
        alpha = float(f['alpha'])
        lats.append(c['lat'])
        lons.append(c['lon'])
        alphas.append(alpha)
        sizes.append(max(8, min(60, c['mw'] / 500)))

sc = ax.scatter(lons, lats, c=alphas, cmap=cmap, norm=norm,
                s=sizes, alpha=0.7, edgecolors='#333333', linewidths=0.3,
                transform=ccrs.PlateCarree(), zorder=5)

cbar = plt.colorbar(sc, ax=ax, orientation='horizontal', pad=0.05,
                     fraction=0.04, shrink=0.6)
cbar.set_label(r'Legacy intensity ($\alpha$: fossil MW / total MW)', fontsize=9)
cbar.set_ticks([0, 0.25, 0.5, 0.75, 1.0])

ax.set_title('Panel A: Global Sample of Listed Power Utilities', fontsize=11, pad=10)

fig.savefig(os.path.join(OUT, 'fig1_global_sample.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'fig1_global_sample.png'), bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'  Saved fig1_global_sample.pdf/png ({len(lats)} firms)')


# =========================================================================
# FIGURE 2: Network structure with W edges
# =========================================================================
print('Generating Figure 2: Network structure...')

# Panel A: Global network
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.add_feature(cfeature.LAND, facecolor='#f5f5f5', edgecolor='none')
ax.add_feature(cfeature.BORDERS, linewidth=0.2, edgecolor='#dddddd')
ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#999999')

# Draw edges first (light)
for gi, gj, wij in W_edges:
    if gi in centroids and gj in centroids:
        ci, cj = centroids[gi], centroids[gj]
        ax.plot([ci['lon'], cj['lon']], [ci['lat'], cj['lat']],
                color='#4488cc', alpha=0.08, linewidth=0.3,
                transform=ccrs.PlateCarree(), zorder=3)

# Nodes colored by density quartile
density_cmap = plt.cm.YlOrRd
max_d = max(c['density'] for c in centroids.values())
norm_d = mcolors.Normalize(vmin=0, vmax=max_d)

for gk, c in centroids.items():
    color = density_cmap(norm_d(c['density']))
    size = max(6, min(40, c['mw'] / 800))
    ax.scatter(c['lon'], c['lat'], c=[color], s=size,
               alpha=0.8, edgecolors='#333333', linewidths=0.2,
               transform=ccrs.PlateCarree(), zorder=5)

sm = plt.cm.ScalarMappable(cmap=density_cmap, norm=norm_d)
cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', pad=0.05,
                     fraction=0.04, shrink=0.6)
cbar.set_label('Network density (exponential kernel, half-life 500 km)', fontsize=9)

ax.set_title('Panel B: Spatial Network Structure', fontsize=11, pad=10)

fig.savefig(os.path.join(OUT, 'fig2_network_structure.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'fig2_network_structure.png'), bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'  Saved fig2_network_structure.pdf/png ({len(W_edges)} edges)')


# =========================================================================
# FIGURE 2B: Regional zoom-ins (Europe, China, India)
# =========================================================================
print('Generating Figure 2B: Regional zoom-ins...')

regions = [
    ('Europe', [-12, 35, 35, 72]),
    ('East & South Asia', [65, 145, -10, 50]),
    ('North America', [-130, -60, 22, 55]),
]

fig, axes = plt.subplots(1, 3, figsize=(14, 4),
                          subplot_kw={'projection': ccrs.PlateCarree()})

for idx, (title, extent) in enumerate(regions):
    ax = axes[idx]
    ax.set_extent(extent, crs=ccrs.PlateCarree())
    ax.add_feature(cfeature.LAND, facecolor='#f0f0f0', edgecolor='none')
    ax.add_feature(cfeature.BORDERS, linewidth=0.4, edgecolor='#cccccc')
    ax.add_feature(cfeature.COASTLINE, linewidth=0.4, edgecolor='#999999')

    # Edges
    for gi, gj, wij in W_edges:
        if gi in centroids and gj in centroids:
            ci, cj = centroids[gi], centroids[gj]
            if (extent[0] <= ci['lon'] <= extent[1] and extent[2] <= ci['lat'] <= extent[3]) or \
               (extent[0] <= cj['lon'] <= extent[1] and extent[2] <= cj['lat'] <= extent[3]):
                ax.plot([ci['lon'], cj['lon']], [ci['lat'], cj['lat']],
                        color='#4488cc', alpha=0.15, linewidth=0.5, zorder=3)

    # Nodes
    for gk, c in centroids.items():
        if extent[0] <= c['lon'] <= extent[1] and extent[2] <= c['lat'] <= extent[3]:
            f = fundamentals.get(gk)
            alpha_val = float(f['alpha']) if f and f.get('alpha', '') != '' else 0.5
            color = cmap(norm(alpha_val))
            size = max(10, min(80, c['mw'] / 300))
            ax.scatter(c['lon'], c['lat'], c=[color], s=size,
                       alpha=0.8, edgecolors='#333333', linewidths=0.3, zorder=5)

    ax.set_title(title, fontsize=10)

fig.suptitle('Panel C: Regional Network Detail', fontsize=11, y=1.02)
fig.tight_layout()
fig.savefig(os.path.join(OUT, 'fig2b_regional_networks.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'fig2b_regional_networks.png'), bbox_inches='tight', dpi=300)
plt.close(fig)
print('  Saved fig2b_regional_networks.pdf/png')


# =========================================================================
# FIGURE 3: Coal retirement first-mover events
# =========================================================================
print('Generating Figure 3: Coal retirement events...')
fig = plt.figure(figsize=(10, 5))
ax = fig.add_subplot(1, 1, 1, projection=ccrs.Robinson())
ax.set_global()
ax.add_feature(cfeature.LAND, facecolor='#f5f5f5', edgecolor='none')
ax.add_feature(cfeature.BORDERS, linewidth=0.2, edgecolor='#dddddd')
ax.add_feature(cfeature.COASTLINE, linewidth=0.3, edgecolor='#999999')

# Color by year
year_cmap = plt.cm.plasma
year_norm = mcolors.Normalize(vmin=2015, vmax=2025)

for r in retirements:
    color = year_cmap(year_norm(r['year']))
    size = max(5, min(80, r['mw'] / 50))
    marker = 'o' if r['matched'] else 'x'
    alpha = 0.8 if r['matched'] else 0.4
    ax.scatter(r['lon'], r['lat'], c=[color], s=size,
               marker=marker, alpha=alpha, edgecolors='none',
               transform=ccrs.PlateCarree(), zorder=5)

sm = plt.cm.ScalarMappable(cmap=year_cmap, norm=year_norm)
cbar = plt.colorbar(sm, ax=ax, orientation='horizontal', pad=0.05,
                     fraction=0.04, shrink=0.6)
cbar.set_label('Retirement year', fontsize=9)

# Legend for matched vs unmatched
from matplotlib.lines import Line2D
legend_elements = [
    Line2D([0], [0], marker='o', color='w', markerfacecolor='#888888',
           markersize=6, label=f'Matched to listed firm (n={sum(1 for r in retirements if r["matched"])})'),
    Line2D([0], [0], marker='x', color='#888888', markerfacecolor='none',
           markersize=6, label=f'Unmatched (n={sum(1 for r in retirements if not r["matched"])})', linestyle='None'),
]
ax.legend(handles=legend_elements, loc='lower left', fontsize=7, framealpha=0.9)

ax.set_title('Panel D: First-Mover Coal Retirement Events (2015\u20132025)', fontsize=11, pad=10)

fig.savefig(os.path.join(OUT, 'fig3_retirement_events.pdf'), bbox_inches='tight')
fig.savefig(os.path.join(OUT, 'fig3_retirement_events.png'), bbox_inches='tight', dpi=300)
plt.close(fig)
print(f'  Saved fig3_retirement_events.pdf/png ({len(retirements)} events)')

print(f'\nAll figures saved to {OUT}')
