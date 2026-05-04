"""Generate Figure 5: cumulative gamma_fuel(tau) over the daily event window.

Visualizes the four phases identified in §5 daily event-study:
  (i)   Pre-event drift (taus -21 to -2)
  (ii)  Announcement reaction (-1 to +1)
  (iii) Post-announcement recovery (+2 to +10)
  (iv)  Long-horizon resumption (+11 to +21)

Inputs:
  - results/metrics/daily_event_time_path.md (parsed)

Outputs:
  - results/figures/fig5_daily_path.pdf
  - manuscript/figures/fig5_daily_path.pdf
"""
import csv
import os
import re
import sys
from collections import OrderedDict

import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import results_path

INPUT = os.path.join(results_path('metrics'), 'daily_event_time_path.md')
OUT_RES = os.path.join(results_path('figures'), 'fig5_daily_path.pdf')
OUT_MS = os.path.join('manuscript', 'figures', 'fig5_daily_path.pdf')


def parse_md_table(path):
    """Parse the markdown table from daily_event_time_path.md."""
    rows = []
    with open(path, 'r', encoding='utf-8') as f:
        for line in f:
            line = line.strip()
            # Match data rows: | tau | T | gamma_fuel | t | gamma_geo | t |
            if not line.startswith('|') or '---' in line:
                continue
            parts = [p.strip() for p in line.strip('|').split('|')]
            if len(parts) < 6:
                continue
            try:
                tau = int(parts[0])
                T = int(parts[1])
                g_fuel = float(parts[2].replace('+', ''))
                t_fuel = float(parts[3].replace('+', ''))
                g_geo = float(parts[4].replace('+', ''))
                t_geo = float(parts[5].replace('+', ''))
                rows.append((tau, T, g_fuel, t_fuel, g_geo, t_geo))
            except (ValueError, IndexError):
                continue
    return sorted(rows, key=lambda r: r[0])


def main():
    rows = parse_md_table(INPUT)
    if not rows:
        sys.exit(f'No data parsed from {INPUT}')

    taus = [r[0] for r in rows]
    g_fuel = [r[2] for r in rows]
    t_fuel = [r[3] for r in rows]

    # Cumulative
    cum = []
    s = 0.0
    for g in g_fuel:
        s += g
        cum.append(s)

    # SE envelope: back out per-day SE from t-stat, then cumulative SE
    # (treats days as independent — slight understatement, but OK for visualization)
    se_daily = [abs(g_fuel[i] / t_fuel[i]) if t_fuel[i] != 0 else 0
                for i in range(len(g_fuel))]
    cum_var = []
    v = 0.0
    for s_ in se_daily:
        v += s_ ** 2
        cum_var.append(v)
    cum_se = [v ** 0.5 for v in cum_var]
    cum_lo = [cum[i] - 1.96 * cum_se[i] for i in range(len(cum))]
    cum_hi = [cum[i] + 1.96 * cum_se[i] for i in range(len(cum))]

    # Plot
    fig, ax = plt.subplots(figsize=(6.5, 3.6), dpi=120)

    # Phase shading
    phases = [
        (-21.5, -1.5, '#e8e8e8', 'Pre-event drift'),
        (-1.5, +1.5, '#d4a76a', 'Announcement'),
        (+1.5, +10.5, '#cce8d4', 'Recovery'),
        (+10.5, +21.5, '#e8d8d8', 'Long horizon'),
    ]
    y_top = max(cum_hi) * 1.05
    y_bot = min(cum_lo) * 1.05
    for x0, x1, color, _ in phases:
        ax.axvspan(x0, x1, alpha=0.45, color=color, zorder=0)

    # CI ribbon
    ax.fill_between(taus, cum_lo, cum_hi, color='gray', alpha=0.25,
                    label='95% CI', zorder=1)

    # Cumulative line
    ax.plot(taus, cum, color='#1a1a1a', linewidth=1.3, zorder=3)
    ax.scatter(taus, cum, color='#1a1a1a', s=14, zorder=4)

    # Reference lines
    ax.axhline(0, color='gray', linewidth=0.4, zorder=2)
    ax.axvline(0, color='gray', linewidth=0.5, linestyle='--', zorder=2)

    # Phase labels
    label_y = y_top * 0.92
    for x0, x1, _, lab in phases:
        x_mid = (x0 + x1) / 2
        ax.text(x_mid, label_y, lab, ha='center', va='top',
                fontsize=8, color='#404040', fontstyle='italic')

    ax.set_xlabel(r'Trading-day offset from announcement ($\tau$)', fontsize=10)
    ax.set_ylabel(r'Cumulative $\hat\gamma_{\mathrm{fuel}}(\tau)$', fontsize=10)
    ax.set_xlim(-21.5, 21.5)
    ax.set_xticks(range(-21, 22, 7))
    ax.tick_params(axis='both', labelsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.grid(axis='y', linewidth=0.3, color='#dddddd', zorder=0)

    plt.tight_layout()

    os.makedirs(os.path.dirname(OUT_RES), exist_ok=True)
    os.makedirs(os.path.dirname(OUT_MS), exist_ok=True)
    fig.savefig(OUT_RES)
    fig.savefig(OUT_MS)
    print(f'Wrote {OUT_RES}')
    print(f'Wrote {OUT_MS}')

    # Print summary so manuscript text matches
    pre_sum = sum(g for tau, _, g, *_ in rows if -21 <= tau <= -2)
    ann_sum = sum(g for tau, _, g, *_ in rows if -1 <= tau <= 1)
    rec_sum = sum(g for tau, _, g, *_ in rows if 2 <= tau <= 10)
    long_sum = sum(g for tau, _, g, *_ in rows if 11 <= tau <= 21)
    print(f'\nPhase sums:')
    print(f'  Pre-event [-21,-2]: {pre_sum:+.2f}')
    print(f'  Announcement [-1,+1]: {ann_sum:+.2f}')
    print(f'  Recovery [+2,+10]: {rec_sum:+.2f}')
    print(f'  Long horizon [+11,+21]: {long_sum:+.2f}')


if __name__ == '__main__':
    main()
