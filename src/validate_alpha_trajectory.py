"""Validate alpha trajectory against ESG scores.

Outputs:
  alpha_trajectory_validation.csv
"""
import csv
import math
import os
from collections import defaultdict

from _paths import raw_path, derived_path

TRAJ_PATH = derived_path('trajectories', 'alpha_trajectory_panel.csv')
PANEL_PATH = raw_path('refinitiv', 'refinitiv_panel.csv')
EXTRA_PATH = raw_path('refinitiv', 'refinitiv_extra.csv')
OUT_PATH = derived_path('validation', 'alpha_trajectory_validation.csv')


def safe_float(x):
    try:
        return float(x)
    except (ValueError, TypeError):
        return None


def corr(xs, ys):
    n = len(xs)
    if n < 3:
        return None
    mx = sum(xs) / n
    my = sum(ys) / n
    num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    denx = sum((x - mx) ** 2 for x in xs)
    deny = sum((y - my) ** 2 for y in ys)
    if denx == 0 or deny == 0:
        return None
    return num / math.sqrt(denx * deny)


# Load alpha slopes
alpha_slope = defaultdict(dict)  # gvkey -> year -> slope
with open(TRAJ_PATH, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        year = int(row['year'])
        slope = safe_float(row.get('slope_log'))
        if slope is None:
            continue
        alpha_slope[gk][year] = slope

# Load panel env_score
env_panel = defaultdict(dict)
with open(PANEL_PATH, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row.get('gvkey')
        if not gk:
            continue
        try:
            year = int(row.get('year'))
        except (ValueError, TypeError):
            continue
        env = safe_float(row.get('env_score'))
        if env is None:
            continue
        env_panel[gk][year] = env

# Build env_score slopes over same window
env_slope = defaultdict(dict)
for gk, yrs in env_panel.items():
    years = sorted(yrs.keys())
    for year in range(min(years) + 1, max(years) + 1):
        window_years = [y for y in years if (year - 3) <= y <= (year - 1)]
        if len(window_years) < 2:
            continue
        xs = window_years
        ys = [env_panel[gk][y] for y in window_years]
        # simple slope
        mx = sum(xs) / len(xs)
        my = sum(ys) / len(ys)
        num = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
        den = sum((x - mx) ** 2 for x in xs)
        if den == 0:
            continue
        env_slope[gk][year] = num / den

# Cross-sectional ESG subscores
extra = {}
if os.path.exists(EXTRA_PATH):
    with open(EXTRA_PATH, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row.get('gvkey')
            if not gk:
                continue
            extra[gk] = {
                'emissions_score': safe_float(row.get('emissions_score')),
                'env_innovation_score': safe_float(row.get('env_innovation_score')),
                'resource_use_score': safe_float(row.get('resource_use_score')),
            }

# Correlate alpha_slope with env_score_slope (panel)
xs, ys = [], []
for gk, years in alpha_slope.items():
    for year, slope in years.items():
        if gk in env_slope and year in env_slope[gk]:
            xs.append(slope)
            ys.append(env_slope[gk][year])

panel_corr = corr(xs, ys)

# Cross-sectional: use latest slope per gvkey
latest_slope = {}
for gk, years in alpha_slope.items():
    yr = max(years.keys())
    latest_slope[gk] = years[yr]

xs_em, ys_em = [], []
xs_inno, ys_inno = [], []
xs_res, ys_res = [], []
for gk, slope in latest_slope.items():
    e = extra.get(gk, {})
    if e.get('emissions_score') is not None:
        xs_em.append(slope)
        ys_em.append(e['emissions_score'])
    if e.get('env_innovation_score') is not None:
        xs_inno.append(slope)
        ys_inno.append(e['env_innovation_score'])
    if e.get('resource_use_score') is not None:
        xs_res.append(slope)
        ys_res.append(e['resource_use_score'])

em_corr = corr(xs_em, ys_em)
inno_corr = corr(xs_inno, ys_inno)
res_corr = corr(xs_res, ys_res)

with open(OUT_PATH, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['metric', 'n', 'corr'])
    w.writerow(['alpha_slope_vs_env_score_slope', len(xs), panel_corr])
    w.writerow(['alpha_slope_vs_emissions_score', len(xs_em), em_corr])
    w.writerow(['alpha_slope_vs_env_innovation_score', len(xs_inno), inno_corr])
    w.writerow(['alpha_slope_vs_resource_use_score', len(xs_res), res_corr])

print(f'Panel slope correlation (alpha vs env_score): {panel_corr}')
print(f'Cross-sectional correlations (latest slope):')
print(f'  emissions_score: {em_corr}')
print(f'  env_innovation_score: {inno_corr}')
print(f'  resource_use_score: {res_corr}')
print(f'Wrote {OUT_PATH}')
