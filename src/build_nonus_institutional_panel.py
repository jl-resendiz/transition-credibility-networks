"""Build non-US institutional concentration panel from Refinitiv free-float data.

PRE-STAGED for Phase 6 of EXECUTION_PLAN.md. Run AFTER pull_refinitiv_extra.py.

Inputs:
  - data/raw/refinitiv/refinitiv_extra.csv  (must exist; produced by pull_refinitiv_extra.py)
    Expected columns: gvkey, ric, free_float_pct, shares_outstanding, ...

Output:
  - data/derived/institutional/nonus_concentration.csv
    Columns: gvkey, ric, free_float_pct, concentrated_ownership_pct, n_shares_outstanding

Concept: Refinitiv free-float is firm-level cross-sectional (one observation per
firm). Use the most recent snapshot. The concentration metric is:
    concentrated_ownership_pct = 1 - free_float_pct/100
where free_float_pct is in [0, 100]. Higher concentrated_ownership_pct = more
ownership held in non-public hands (insiders, sovereign, strategic blocks),
which is a proxy for institutional concentration in the non-US setting where
13F filings are unavailable.
"""
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import raw_path, derived_path

IN_PATH = raw_path('refinitiv', 'refinitiv_extra.csv')
OUT_PATH = os.path.join(derived_path('institutional'), 'nonus_concentration.csv')


def main():
    if not os.path.exists(IN_PATH):
        print(f'SKIP: {IN_PATH} not found.')
        print('  Run pull_refinitiv_extra.py first (requires Eikon access).')
        print('  This is the optional Phase 6 non-US institutional extension.')
        sys.exit(0)

    print(f'Reading {IN_PATH}...')
    rows = []
    with open(IN_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                gk = str(row['gvkey']).strip().split('.')[0].zfill(6)
                ric = row.get('ric', '').strip()
                ff = row.get('free_float_pct', '').strip()
                so = row.get('shares_outstanding', '').strip()
            except KeyError:
                continue
            if not gk or gk == 'NAN':
                continue

            # Free-float can be in [0, 100] (percent) or [0, 1] (fraction)
            ff_val = None
            if ff and ff not in ('', 'NA', 'nan', 'None'):
                try:
                    ff_val = float(ff)
                    # Normalize to [0, 100] if Refinitiv returned a fraction
                    if ff_val <= 1.0 and ff_val > 0:
                        ff_val = ff_val * 100.0
                    if ff_val < 0 or ff_val > 100:
                        ff_val = None
                except ValueError:
                    pass

            so_val = None
            if so and so not in ('', 'NA', 'nan', 'None'):
                try:
                    so_val = float(so)
                except ValueError:
                    pass

            if ff_val is None:
                continue  # Cannot compute concentration without free-float

            rows.append({
                'gvkey': gk,
                'ric': ric,
                'free_float_pct': ff_val,
                'concentrated_ownership_pct': max(0.0, 100.0 - ff_val),
                'shares_outstanding': so_val if so_val is not None else '',
            })

    print(f'  {len(rows)} firms with valid free-float data')
    if not rows:
        sys.exit('No usable rows. Check refinitiv_extra.csv content.')

    # Cross-check with the existing 13F US sub-sample. The concentrated_ownership
    # measure should NOT collide with 13F-covered (US-listed) firms; we use this
    # only for non-US firms.
    us_path = os.path.join(derived_path('institutional'), 'institutional_ownership.csv')
    us_gvkeys = set()
    if os.path.exists(us_path):
        with open(us_path, 'r', encoding='utf-8') as f:
            for r in csv.DictReader(f):
                us_gvkeys.add(str(r['gvkey']).split('.')[0].zfill(6))
    print(f'  US-13F-covered firms (to flag): {len(us_gvkeys)}')

    n_nonus = sum(1 for r in rows if r['gvkey'] not in us_gvkeys)
    print(f'  Non-US firms with free-float: {n_nonus}')

    # Distribution diagnostics
    ff_values = sorted(r['free_float_pct'] for r in rows if r['gvkey'] not in us_gvkeys)
    if ff_values:
        n = len(ff_values)
        print('\n  Free-float distribution (non-US sample):')
        print(f'    p10:    {ff_values[n//10]:.1f}%')
        print(f'    p25:    {ff_values[n//4]:.1f}%')
        print(f'    median: {ff_values[n//2]:.1f}%')
        print(f'    p75:    {ff_values[3*n//4]:.1f}%')
        print(f'    p90:    {ff_values[9*n//10]:.1f}%')
        print(f'    range:  {ff_values[0]:.1f}% – {ff_values[-1]:.1f}%')

    # Tag US vs non-US (mainly for downstream filtering)
    for r in rows:
        r['us_listed'] = '1' if r['gvkey'] in us_gvkeys else '0'

    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
    with open(OUT_PATH, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=['gvkey', 'ric', 'free_float_pct',
                                          'concentrated_ownership_pct',
                                          'shares_outstanding', 'us_listed'])
        w.writeheader()
        w.writerows(rows)
    print(f'\nWrote {OUT_PATH}')

    # Decision trigger note for the user
    if ff_values:
        n = len(ff_values)
        clustered_high = sum(1 for v in ff_values if v >= 80)
        if clustered_high / n > 0.85:
            print('\nDECISION TRIGGER: free-float distribution is clustered at high values')
            print('  (>85% of firms have free_float >= 80%). Tercile split may have low power.')
            print('  Spawn /quant-finance agent to discuss alternative cuts (e.g., quintile,')
            print('  or use complementary insider-holdings concentration).')


if __name__ == '__main__':
    main()
