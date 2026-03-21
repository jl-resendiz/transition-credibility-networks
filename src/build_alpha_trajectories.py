"""Build pre-event alpha (CO2/Revenue) trajectories from Refinitiv panel.

Outputs: alpha_trajectory_panel.csv
Columns:
  gvkey, year, n_points, slope_log, trend_group
Where slope_log is the OLS slope of log(CO2/Revenue) over years [year-3, year-1].
"""
import csv
import math
import os
from collections import defaultdict

from _paths import raw_path, derived_path

PANEL_PATH = raw_path('refinitiv', 'refinitiv_panel.csv')
OUT_PATH = derived_path('trajectories', 'alpha_trajectory_panel.csv')

WINDOW = 3  # years back from t-1
EPS = 0.02  # threshold for trend classification in log-units per year


def safe_float(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def ols_slope(xs, ys):
    n = len(xs)
    if n < 2:
        return None
    xbar = sum(xs) / n
    ybar = sum(ys) / n
    num = sum((x - xbar) * (y - ybar) for x, y in zip(xs, ys))
    den = sum((x - xbar) ** 2 for x in xs)
    if den == 0:
        return None
    return num / den


def classify_trend(slope):
    if slope is None:
        return ''
    if slope <= -EPS:
        return 'declining'
    if slope >= EPS:
        return 'rising'
    return 'stable'


# Load panel CO2/Revenue
panel = defaultdict(dict)  # gvkey -> year -> co2_to_revenue
years_all = set()
if not os.path.exists(PANEL_PATH):
    raise SystemExit('refinitiv_panel.csv not found')

with open(PANEL_PATH, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row.get('gvkey')
        if not gk:
            continue
        try:
            year = int(row.get('year'))
        except (ValueError, TypeError):
            continue
        co2_rev = safe_float(row.get('co2_to_revenue'))
        if co2_rev is None or co2_rev <= 0:
            continue
        panel[gk][year] = co2_rev
        years_all.add(year)

years_all = sorted(years_all)
if not years_all:
    raise SystemExit('No CO2/Revenue values found in refinitiv_panel.csv')

min_year = min(years_all)
max_year = max(years_all)

rows_out = []

for gk, yr_map in panel.items():
    years = sorted(yr_map.keys())
    for year in range(min_year + 1, max_year + 1):
        # pre-event window: [year-3, year-1]
        window_years = [y for y in years if (year - WINDOW) <= y <= (year - 1)]
        if len(window_years) < 2:
            continue
        xs = window_years
        ys = []
        for y in window_years:
            val = yr_map.get(y)
            if val is None or val <= 0:
                continue
            ys.append(math.log(val))
        if len(ys) != len(xs) or len(ys) < 2:
            continue
        slope = ols_slope(xs, ys)
        rows_out.append({
            'gvkey': gk,
            'year': year,
            'n_points': len(xs),
            'slope_log': f'{slope:.6f}' if slope is not None else '',
            'trend_group': classify_trend(slope),
        })

with open(OUT_PATH, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=['gvkey', 'year', 'n_points', 'slope_log', 'trend_group'])
    w.writeheader()
    w.writerows(rows_out)

print(f'Wrote {len(rows_out)} rows to {OUT_PATH}')
