"""Generate publication-quality figures for the manuscript.

Figure 1: Portfolio-sort quintile CARs (fuel similarity and geographic proximity).
Figure 2: Fuel-similarity Q5-Q1 spread by ETS status.

Data sources: strategy2_portfolio_sorts.md, strategy2_credibility_interaction.md.
"""

import os
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from _paths import results_path

# ---------------------------------------------------------------------------
# Global style
# ---------------------------------------------------------------------------
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 11,
    'axes.linewidth': 0.5,
    'axes.spines.top': False,
    'axes.spines.right': False,
    'figure.dpi': 300,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
    'savefig.pad_inches': 0.1,
})

# Output directory
ROOT_DIR = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
FIG_DIR = os.path.join(ROOT_DIR, "manuscript", "figures")
os.makedirs(FIG_DIR, exist_ok=True)

# Colours that reproduce well in greyscale
COLOUR_FUEL = '#C44E52'   # muted red/brick
COLOUR_GEO  = '#4C72B0'   # muted blue


# ===================================================================
# Figure 1 — Portfolio-sort quintile CARs
# ===================================================================
def figure_portfolio_sorts():
    """Side-by-side vertical bar charts for fuel and geographic quintile CARs."""

    quintiles = ['Q1', 'Q2', 'Q3', 'Q4', 'Q5']
    fuel_car  = [0.84, 1.19, -4.68, -4.14, -4.61]
    geo_car   = [-5.05, -0.23, -8.14, 1.15, 0.98]

    fuel_spread = fuel_car[4] - fuel_car[0]   # Q5 - Q1
    geo_spread  = geo_car[4] - geo_car[0]

    fig, (ax_fuel, ax_geo) = plt.subplots(1, 2, figsize=(10, 4.5),
                                           sharey=True)

    # --- Panel A: Fuel Similarity ---
    bars_f = ax_fuel.bar(quintiles, fuel_car, color=COLOUR_FUEL,
                         edgecolor='white', width=0.6, zorder=3)
    ax_fuel.axhline(0, color='grey', linewidth=0.5, zorder=1)
    ax_fuel.set_title('Panel A: Fuel Similarity', fontsize=13,
                      fontweight='bold', pad=10)
    ax_fuel.set_ylabel('Mean CAR (%)', fontsize=12)
    ax_fuel.set_xlabel('Quintile (low to high similarity)', fontsize=12)

    # Label each bar
    for bar, val in zip(bars_f, fuel_car):
        y_offset = 0.25 if val >= 0 else -0.45
        ax_fuel.text(bar.get_x() + bar.get_width() / 2, val + y_offset,
                     f'{val:+.2f}', ha='center', va='bottom' if val >= 0 else 'top',
                     fontsize=10, color='black')

    # Annotate spread
    ax_fuel.annotate(
        f'Q5\u2013Q1 = {fuel_spread:+.2f} pp',
        xy=(0.97, 0.03), xycoords='axes fraction',
        ha='right', va='bottom', fontsize=11,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0',
                  edgecolor='grey', linewidth=0.5))

    # --- Panel B: Geographic Proximity ---
    bars_g = ax_geo.bar(quintiles, geo_car, color=COLOUR_GEO,
                        edgecolor='white', width=0.6, zorder=3)
    ax_geo.axhline(0, color='grey', linewidth=0.5, zorder=1)
    ax_geo.set_title('Panel B: Geographic Proximity', fontsize=13,
                     fontweight='bold', pad=10)
    ax_geo.set_xlabel('Quintile (low to high proximity)', fontsize=12)

    for bar, val in zip(bars_g, geo_car):
        y_offset = 0.25 if val >= 0 else -0.45
        ax_geo.text(bar.get_x() + bar.get_width() / 2, val + y_offset,
                    f'{val:+.2f}', ha='center', va='bottom' if val >= 0 else 'top',
                    fontsize=10, color='black')

    ax_geo.annotate(
        f'Q5\u2013Q1 = {geo_spread:+.2f} pp',
        xy=(0.97, 0.03), xycoords='axes fraction',
        ha='right', va='bottom', fontsize=11,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0',
                  edgecolor='grey', linewidth=0.5))

    # Shared formatting
    for ax in (ax_fuel, ax_geo):
        ax.tick_params(axis='both', labelsize=11)
        ax.set_ylim(-10, 3)

    fig.tight_layout(w_pad=3)

    pdf_path = os.path.join(FIG_DIR, 'fig_portfolio_sorts.pdf')
    png_path = os.path.join(FIG_DIR, 'fig_portfolio_sorts.png')
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    print(f'Saved  {pdf_path}')
    print(f'Saved  {png_path}')


# ===================================================================
# Figure 2 — ETS interaction: fuel-similarity spread
# ===================================================================
def figure_ets_interaction():
    """Bar chart comparing fuel Q5-Q1 spread for ETS vs non-ETS firms."""

    labels  = ['ETS firms', 'Non-ETS firms']
    spreads = [-2.54, 1.15]
    colours = [COLOUR_FUEL, COLOUR_GEO]
    diff_t  = -3.66

    fig, ax = plt.subplots(figsize=(5, 4.5))

    bars = ax.bar(labels, spreads, color=colours, edgecolor='white',
                  width=0.5, zorder=3)
    ax.axhline(0, color='grey', linewidth=0.5, zorder=1)

    # Label bars
    for bar, val in zip(bars, spreads):
        y_offset = 0.15 if val >= 0 else -0.15
        ax.text(bar.get_x() + bar.get_width() / 2, val + y_offset,
                f'{val:+.2f}%', ha='center',
                va='bottom' if val >= 0 else 'top',
                fontsize=12, fontweight='bold', color='black')

    ax.set_ylabel('Fuel-Similarity Q5\u2013Q1 Spread (%)', fontsize=12)
    ax.set_title('Fuel-Similarity Quintile Spread by ETS Status',
                 fontsize=13, fontweight='bold', pad=12)
    ax.tick_params(axis='both', labelsize=12)

    # Annotate difference
    diff_val = spreads[0] - spreads[1]
    ax.annotate(
        f'Difference = {diff_val:+.2f} pp\n$t$ = {diff_t:.2f}',
        xy=(0.97, 0.97), xycoords='axes fraction',
        ha='right', va='top', fontsize=11,
        bbox=dict(boxstyle='round,pad=0.3', facecolor='#f0f0f0',
                  edgecolor='grey', linewidth=0.5))

    ax.set_ylim(-4.5, 2.5)

    fig.tight_layout()

    pdf_path = os.path.join(FIG_DIR, 'fig_ets_interaction.pdf')
    png_path = os.path.join(FIG_DIR, 'fig_ets_interaction.png')
    fig.savefig(pdf_path)
    fig.savefig(png_path, dpi=300)
    plt.close(fig)
    print(f'Saved  {pdf_path}')
    print(f'Saved  {png_path}')


# ===================================================================
# Main
# ===================================================================
if __name__ == '__main__':
    print('Generating publication figures ...')
    figure_portfolio_sorts()
    figure_ets_interaction()
    print('Done.')
