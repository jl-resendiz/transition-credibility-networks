"""Construct firm fundamental variables (theta vector) from Compustat + GEM matches.

Variables per the paper:
  alpha (legacy intensity): CO2/Revenue from Refinitiv panel (preferred),
        Refinitiv cross-section (fallback), or fossil_MW/total_MW from GEM
  lambda (leverage): (dltt + dlc) / at
  rho (return spread): oibdp / at  (operating ROA proxy)
  kappa (cash flow adequacy): oancf / xint  (interest coverage)
  delta (obligation rigidity): dltt / (dltt + dlc)  (long-term debt share)

Alpha source priority:
  1. Refinitiv panel (refinitiv_panel.csv) — exact (gvkey, fyear) match
  2. Refinitiv panel — nearest year for same firm
  3. Refinitiv cross-section (refinitiv_esg.csv) — latest observation
  4. GEM capacity (gem_compustat_matches.csv) — fossil_MW / total_MW
"""
import csv, os
from collections import defaultdict

from _paths import raw_path, derived_path


def safe_float(val):
    try:
        return float(val)
    except (ValueError, TypeError):
        return None


# ---------- 1. GEM capacity-based alpha (lowest priority fallback) ----------
gvkey_mw = defaultdict(lambda: {'coal_mw': 0, 'gas_mw': 0, 'solar_mw': 0, 'wind_mw': 0})
with open(derived_path('mappings', 'gem_compustat_matches.csv'), 'r', encoding='utf-8') as f:
    for row in csv.DictReader(f):
        gk = row['gvkey']
        gvkey_mw[gk]['coal_mw'] += float(row['coal_mw'])
        gvkey_mw[gk]['gas_mw'] += float(row['gas_mw'])
        gvkey_mw[gk]['solar_mw'] += float(row['solar_mw'])
        gvkey_mw[gk]['wind_mw'] += float(row['wind_mw'])

gvkey_alpha_gem = {}
for gk, mw in gvkey_mw.items():
    fossil = mw['coal_mw'] + mw['gas_mw']
    total = fossil + mw['solar_mw'] + mw['wind_mw']
    gvkey_alpha_gem[gk] = fossil / total if total > 0 else None

print(f'GEM-matched gvkeys with alpha: {len(gvkey_alpha_gem)}')


# ---------- 2. Refinitiv panel alpha (highest priority) ----------
# Build (gvkey, year) -> scaled CO2/Revenue from annual panel
panel_path = raw_path('refinitiv', 'refinitiv_panel.csv')
panel_alpha = {}        # (gvkey, year_str) -> scaled alpha
panel_years = {}        # gvkey -> sorted list of available years
_panel_raw = {}         # (gvkey, year_str) -> raw co2_to_revenue (for winsorization)

if os.path.exists(panel_path):
    with open(panel_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row.get('gvkey')
            if not gk:
                continue
            co2_rev = safe_float(row.get('co2_to_revenue'))
            if co2_rev is None:
                continue
            year = row.get('year', '')
            if not year:
                continue
            _panel_raw[(gk, str(year))] = co2_rev

    # Winsorize and scale panel values to [0,1]
    if _panel_raw:
        all_vals = sorted(_panel_raw.values())
        n = len(all_vals)
        p05 = all_vals[int(0.05 * (n - 1))]  # 5th/95th for CO2-to-revenue scaling (not Fama-French)
        p95 = all_vals[int(0.95 * (n - 1))]
        denom = (p95 - p05) if p95 > p05 else None
        if denom:
            for (gk, yr), x in _panel_raw.items():
                scaled = max(0.0, min(1.0, (x - p05) / denom))
                panel_alpha[(gk, yr)] = scaled
                panel_years.setdefault(gk, []).append(yr)
            for gk in panel_years:
                panel_years[gk] = sorted(panel_years[gk])

    print(f'Refinitiv panel alpha: {len(panel_alpha)} firm-year obs '
          f'({len(panel_years)} firms)')
else:
    print('No refinitiv_panel.csv found (panel alpha not available)')


# ---------- 3. Refinitiv cross-sectional alpha (second fallback) ----------
ref_path = raw_path('refinitiv', 'refinitiv_esg.csv')
gvkey_alpha_ref = {}

if os.path.exists(ref_path):
    latest = {}
    with open(ref_path, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row.get('gvkey')
            if not gk:
                continue
            co2_rev = safe_float(row.get('co2_to_revenue'))
            if co2_rev is None:
                continue
            # Keep latest observation per gvkey
            year = None
            yr_str = row.get('refinitiv_period') or ''
            if yr_str and yr_str[:4].isdigit():
                year = int(yr_str[:4])
            if gk not in latest or (year and (latest[gk]['year'] is None or year > latest[gk]['year'])):
                latest[gk] = {'year': year, 'co2_to_revenue': co2_rev}

    vals = sorted(v['co2_to_revenue'] for v in latest.values())
    if vals:
        n = len(vals)
        p05 = vals[int(0.05 * (n - 1))]  # 5th/95th for CO2-to-revenue scaling (not Fama-French)
        p95 = vals[int(0.95 * (n - 1))]
        denom = (p95 - p05) if p95 > p05 else None
        if denom:
            for gk, v in latest.items():
                scaled = max(0.0, min(1.0, (v['co2_to_revenue'] - p05) / denom))
                gvkey_alpha_ref[gk] = scaled

    print(f'Refinitiv cross-sectional alpha: {len(gvkey_alpha_ref)} firms')
else:
    print('No refinitiv_esg.csv found (cross-sectional alpha not available)')


# ---------- Alpha lookup function ----------
def get_alpha(gk, fyear):
    """Return (alpha, source_tag) using priority chain: panel > xsec > GEM."""
    yr = str(fyear)

    # Priority 1: Panel exact year match
    if (gk, yr) in panel_alpha:
        return panel_alpha[(gk, yr)], 'panel'

    # Priority 2: Panel nearest year for same firm
    if gk in panel_years:
        years = panel_years[gk]
        closest = min(years, key=lambda y: abs(int(y) - int(yr)))
        if abs(int(closest) - int(yr)) <= 3:  # within 3 years
            return panel_alpha[(gk, closest)], 'panel_near'

    # Priority 3: Refinitiv cross-sectional (latest observation)
    if gk in gvkey_alpha_ref:
        return gvkey_alpha_ref[gk], 'refinitiv'

    # Priority 4: GEM capacity-based
    if gk in gvkey_alpha_gem and gvkey_alpha_gem[gk] is not None:
        return gvkey_alpha_gem[gk], 'gem'

    return None, ''


# ---------- 4. Build fundamentals from Compustat ----------
output_rows = []
alpha_source_counts = defaultdict(int)

for src_file in ['compustat_global_utilities.csv', 'compustat_na_utilities.csv']:
    fpath = raw_path('compustat', src_file)
    if not os.path.exists(fpath):
        continue
    with open(fpath, 'r', encoding='utf-8') as f:
        for row in csv.DictReader(f):
            gk = row['gvkey']
            fyear = row['fyear']
            at = safe_float(row.get('at'))
            dltt = safe_float(row.get('dltt'))
            dlc = safe_float(row.get('dlc'))
            oibdp = safe_float(row.get('oibdp'))
            xint = safe_float(row.get('xint'))
            oancf = safe_float(row.get('oancf'))
            sale = safe_float(row.get('sale'))
            ppent = safe_float(row.get('ppent'))
            capx = safe_float(row.get('capx'))
            ceq = safe_float(row.get('ceq'))

            # lambda: leverage
            lam = None
            if dltt is not None and dlc is not None and at and at > 0:
                lam = (dltt + dlc) / at

            # rho: operating return on assets
            rho = None
            if oibdp is not None and at and at > 0:
                rho = oibdp / at

            # kappa: cash flow adequacy (interest coverage)
            kappa = None
            if oancf is not None and xint and xint > 0:
                kappa = oancf / xint

            # delta: obligation rigidity (long-term debt share)
            delta = None
            if dltt is not None and dlc is not None and (dltt + dlc) > 0:
                delta = dltt / (dltt + dlc)

            # alpha with priority chain
            alpha, alpha_src = get_alpha(gk, fyear)
            if alpha_src:
                alpha_source_counts[alpha_src] += 1

            output_rows.append({
                'gvkey': gk,
                'fyear': fyear,
                'conm': row['conm'],
                'fic': row['fic'],
                'sic': row['sic'],
                'at': at,
                'sale': sale,
                'ppent': ppent,
                'capx': capx,
                'ceq': ceq,
                'dltt': dltt,
                'dlc': dlc,
                'oibdp': oibdp,
                'xint': xint,
                'oancf': oancf,
                'alpha': f'{alpha:.4f}' if alpha is not None else '',
                'lambda': f'{lam:.4f}' if lam is not None else '',
                'rho': f'{rho:.4f}' if rho is not None else '',
                'kappa': f'{kappa:.4f}' if kappa is not None else '',
                'delta': f'{delta:.4f}' if delta is not None else '',
                'source': src_file,
            })

# Save
outpath = derived_path('fundamentals', 'firm_fundamentals.csv')
fieldnames = ['gvkey', 'fyear', 'conm', 'fic', 'sic', 'at', 'sale', 'ppent', 'capx', 'ceq',
              'dltt', 'dlc', 'oibdp', 'xint', 'oancf', 'alpha', 'lambda', 'rho', 'kappa', 'delta', 'source']
with open(outpath, 'w', newline='', encoding='utf-8') as f:
    w = csv.DictWriter(f, fieldnames=fieldnames)
    w.writeheader()
    w.writerows(output_rows)

# Summary stats
n_rows = len(output_rows)
n_firms = len(set(r['gvkey'] for r in output_rows))
n_with_alpha = sum(1 for r in output_rows if r['alpha'] != '')
n_with_all = sum(1 for r in output_rows if all(r[v] != '' for v in ['alpha', 'lambda', 'rho', 'kappa', 'delta']))

print(f'\nOutput: {n_rows} firm-years, {n_firms} unique firms')
print(f'With alpha: {n_with_alpha} firm-years')
print(f'With complete theta vector: {n_with_all} firm-years')

print(f'\nAlpha source breakdown:')
for src, count in sorted(alpha_source_counts.items(), key=lambda x: -x[1]):
    label = {
        'panel': 'Refinitiv panel (exact year)',
        'panel_near': 'Refinitiv panel (nearest year)',
        'refinitiv': 'Refinitiv cross-section',
        'gem': 'GEM capacity',
    }.get(src, src)
    print(f'  {label}: {count} firm-years')

# Distribution of alpha for most recent year
latest = {}
for r in output_rows:
    if r['alpha'] != '':
        gk = r['gvkey']
        if gk not in latest or r['fyear'] > latest[gk]['fyear']:
            latest[gk] = r
alphas = [float(r['alpha']) for r in latest.values()]
alphas.sort()
n = len(alphas)
if n > 0:
    print(f'\nAlpha distribution (latest year, {n} firms):')
    print(f'  Min:    {alphas[0]:.3f}')
    print(f'  P25:    {alphas[n//4]:.3f}')
    print(f'  Median: {alphas[n//2]:.3f}')
    print(f'  P75:    {alphas[3*n//4]:.3f}')
    print(f'  Max:    {alphas[-1]:.3f}')
    print(f'  Mean:   {sum(alphas)/n:.3f}')

print(f'\nSaved to {outpath}')
