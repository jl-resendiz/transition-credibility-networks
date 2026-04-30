"""Newey-West lag sensitivity for the Fama-MacBeth coefficients.

The headline inference uses NW lag=4 on the time series of event-level FM
coefficients. With 5-month CAR windows and 63% of calendar months hosting
multiple active event windows (`robust_inference.py`, Approach 8),
serial correlation in the FM beta series can extend beyond lag 4.

Standard lag rules of thumb for T=117 events:
  - Newey-West (1994): L = 4 * (T/100)^(2/9) approx 4.2
  - Stock-Watson (2008): L = 0.75 * T^(1/3) approx 3.7
  - T^(1/4) rule: approx 3.3

These rules ignore the overlapping-window structure. Because each window is
5 months wide and many windows overlap, the effective serial correlation in
the FM time series can extend to roughly the window width (5) plus the
typical inter-event gap. We therefore report SEs at lags {4, 8, 12, 18} for
transparency.

Driscoll-Kraay (1998) standard errors are conceptually for panel data with
contemporaneous cross-sectional dependence. Applied to a Fama-MacBeth time
series of scalar coefficients, the cross-section has already been averaged
out, so Driscoll-Kraay reduces algebraically to Newey-West with the same
lag. We therefore report only the Newey-West lag sensitivity.

Inputs:  results/summaries/event_level_betas.csv
Outputs: results/metrics/lag_sensitivity.md
         results/summaries/lag_sensitivity.csv
"""
import csv
import math
import os
import sys

from _paths import results_path


def _print(msg=''):
    print(msg)
    sys.stdout.flush()


# ── Inputs ──────────────────────────────────────────────────────────

in_csv = results_path('summaries', 'event_level_betas.csv')
if not os.path.exists(in_csv):
    raise FileNotFoundError(
        f'{in_csv} not found. Run robust_inference.py first.')

events = []
with open(in_csv, 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        try:
            events.append({
                'event_id': int(row['event_id']),
                'beta_fuel': float(row['beta_fuel']),
                'beta_geo': float(row['beta_geo']),
                'beta_reg': float(row['beta_reg']),
                'beta_same_sector': float(row['beta_same_sector'])
                if row['beta_same_sector'] not in ('', 'nan') else float('nan'),
            })
        except (ValueError, KeyError):
            continue

T = len(events)
_print(f'Loaded {T} event-level FM coefficients from {in_csv}')


# ── Newey-West HAC SE for the mean of a time series ─────────────────

def newey_west_se(series, lag):
    """Newey-West (1987) HAC SE for the sample mean of a stationary series.

    Uses the Bartlett kernel: weight(L) = 1 - L/(lag + 1).
    """
    series = [x for x in series if not (isinstance(x, float) and math.isnan(x))]
    n = len(series)
    if n < 2:
        return float('nan')
    mean = sum(series) / n
    dev = [x - mean for x in series]
    gamma0 = sum(d * d for d in dev) / n
    var_nw = gamma0
    for L in range(1, min(lag, n - 1) + 1):
        weight = 1.0 - L / (lag + 1)
        cov_L = sum(dev[t] * dev[t - L] for t in range(L, n)) / n
        var_nw += 2.0 * weight * cov_L
    if var_nw <= 0:
        return float('nan')
    return math.sqrt(var_nw / n)


def fm_summary(series, lags):
    """Mean and NW SE at each lag, plus simple SE (lag=0)."""
    series = [x for x in series if not (isinstance(x, float) and math.isnan(x))]
    n = len(series)
    if n < 2:
        return None
    mean = sum(series) / n
    sd = math.sqrt(sum((x - mean) ** 2 for x in series) / (n - 1))
    se_simple = sd / math.sqrt(n)
    out = {'mean': mean, 'n': n, 'se_simple': se_simple}
    for L in lags:
        out[f'se_nw{L}'] = newey_west_se(series, L)
    return out


# ── Compute per channel and difference ──────────────────────────────

LAGS = [4, 8, 12, 18]

channels = {
    'fuel': [e['beta_fuel'] for e in events],
    'geo': [e['beta_geo'] for e in events],
    'reg': [e['beta_reg'] for e in events],
    'diff_geo_minus_fuel': [e['beta_geo'] - e['beta_fuel'] for e in events],
}

results = {name: fm_summary(series, LAGS) for name, series in channels.items()}


# ── Output ──────────────────────────────────────────────────────────

# CSV
out_csv = results_path('summaries', 'lag_sensitivity.csv')
with open(out_csv, 'w', newline='', encoding='utf-8') as f:
    w = csv.writer(f)
    w.writerow(['channel', 'mean', 'n', 'se_iid'] +
               [f'se_nw_lag{L}' for L in LAGS] +
               [f't_nw_lag{L}' for L in LAGS])
    for name in ['fuel', 'geo', 'reg', 'diff_geo_minus_fuel']:
        r = results[name]
        if r is None:
            continue
        row = [name, f'{r["mean"]:+.6f}', r['n'], f'{r["se_simple"]:.6f}']
        for L in LAGS:
            se = r[f'se_nw{L}']
            row.append(f'{se:.6f}' if not math.isnan(se) else 'NA')
        for L in LAGS:
            se = r[f'se_nw{L}']
            t = r['mean'] / se if (not math.isnan(se) and se > 0) else float('nan')
            row.append(f'{t:+.4f}' if not math.isnan(t) else 'NA')
        w.writerow(row)
_print(f'Wrote {out_csv}')

# Markdown
out_md = results_path('metrics', 'lag_sensitivity.md')
with open(out_md, 'w', encoding='utf-8') as f:
    f.write('# Newey-West Lag Sensitivity for Fama-MacBeth Coefficients\n\n')
    f.write(f'Time series of event-level FM coefficients: T = {T} events.\n\n')
    f.write('Standard lag rules of thumb suggest lags of 4 to 5 for this T. '
            'However, the event windows are 5 months wide and 63% of months '
            'host multiple overlapping windows, which can extend serial '
            'correlation in the FM time series beyond rule-of-thumb lags. '
            'We therefore report lags {4, 8, 12, 18} for transparency.\n\n')

    f.write('Driscoll-Kraay (1998) HAC standard errors are designed for panel '
            'data with contemporaneous cross-sectional dependence. The Fama-'
            'MacBeth time series has already collapsed the cross-section into '
            'a scalar coefficient per event, so Driscoll-Kraay applied to this '
            'series reduces algebraically to Newey-West with the same lag. '
            'Only Newey-West is reported.\n\n')

    f.write('## Lag sensitivity table\n\n')
    header = ('| Channel | Mean | SE iid | SE NW(4) | SE NW(8) | SE NW(12) | '
              'SE NW(18) | t NW(4) | t NW(8) | t NW(12) | t NW(18) |\n')
    f.write(header)
    f.write('|---|---|---|---|---|---|---|---|---|---|---|\n')
    for name, label in [
        ('fuel', '$\\gamma_{\\text{fuel}}$'),
        ('geo', '$\\gamma_{\\text{geo}}$'),
        ('reg', '$\\gamma_{\\text{reg}}$'),
        ('diff_geo_minus_fuel', '$\\gamma_{\\text{geo}} - \\gamma_{\\text{fuel}}$'),
    ]:
        r = results[name]
        if r is None:
            continue
        f.write(f'| {label} | {r["mean"]:+.4f} | {r["se_simple"]:.4f}')
        for L in LAGS:
            se = r[f'se_nw{L}']
            f.write(f' | {se:.4f}' if not math.isnan(se) else ' | NA')
        for L in LAGS:
            se = r[f'se_nw{L}']
            t = r['mean'] / se if (not math.isnan(se) and se > 0) else float('nan')
            f.write(f' | {t:+.2f}' if not math.isnan(t) else ' | NA')
        f.write(' |\n')

    f.write('\n## Interpretation\n\n')
    f.write('A coefficient that loses statistical significance only at long '
            'lags suggests serial correlation in the FM time series that the '
            'baseline lag may understate. A coefficient whose t-statistic is '
            'stable across lags indicates a result that does not depend on '
            'lag choice.\n\n')

    # Auto-interpretation for the fuel channel
    r_fuel = results['fuel']
    t_fuel_lags = [r_fuel['mean'] / r_fuel[f'se_nw{L}']
                   if (not math.isnan(r_fuel[f'se_nw{L}'])
                       and r_fuel[f'se_nw{L}'] > 0) else float('nan')
                   for L in LAGS]
    valid_t = [abs(t) for t in t_fuel_lags if not math.isnan(t)]
    if valid_t:
        min_t = min(valid_t)
        max_t = max(valid_t)
        f.write(f'For the headline fuel coefficient, |t| ranges from '
                f'{min_t:.2f} to {max_t:.2f} across lags {LAGS}. ')
        if min_t > 1.96:
            f.write('All lags reject the null at the 5% level, so the result '
                    'does not depend on the lag choice within this range.\n')
        elif min_t > 1.65:
            f.write('All lags reject at the 10% level; the strictest rejection '
                    'at 5% is borderline at the longest lag.\n')
        else:
            f.write('At least one lag fails to reject at conventional levels; '
                    'the result is sensitive to lag choice.\n')

_print(f'Wrote {out_md}')
_print('Done.')
