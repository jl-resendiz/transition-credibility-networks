"""Build institutional ownership panel from raw 13F holdings.

Inputs:
  - data/raw/wrds/holdings_13f.csv (raw 13F holdings, one row per
    manager-firm-quarter)

Output:
  - data/derived/institutional_ownership.csv
    cols: gvkey, fdate, total_inst_shares, n_managers, hhi_managers,
          top5_concentration, shrout_thousands, inst_pct, inst_pct_capped

Units:
  - shares (13F): thousands of shares held by each manager
  - shrout2 (13F): thousands of total shares outstanding (firm-level, per fdate)
  - inst_pct = SUM(shares) / shrout2  (cap at 100%)
"""
import csv
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _paths import raw_path, derived_path

IN_PATH = os.path.join(raw_path('wrds'), 'holdings_13f.csv')
OUT_PATH = os.path.join(derived_path('institutional'), 'institutional_ownership.csv')


def main():
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

    # Aggregate using stdlib (file is 131 MB, but only 1.7M rows)
    panel = defaultdict(lambda: {
        'shares_per_mgr': defaultdict(float),
        'shrout': set(),
    })

    print(f'Reading {IN_PATH}...')
    with open(IN_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for i, row in enumerate(reader):
            if i % 200_000 == 0:
                print(f'  rows: {i:,}')
            gvkey = row.get('gvkey', '').strip()
            if not gvkey or gvkey == 'nan':
                continue
            try:
                fdate = row['fdate']
                mgrno = row['mgrno']
                shares = float(row['shares'])
                shrout2 = row.get('shrout2', '').strip()
            except (KeyError, ValueError):
                continue
            key = (gvkey, fdate)
            panel[key]['shares_per_mgr'][mgrno] += shares
            if shrout2 and shrout2 not in ('', 'nan'):
                try:
                    panel[key]['shrout'].add(float(shrout2))
                except ValueError:
                    pass

    print(f'\nAggregated: {len(panel):,} (gvkey, fdate) rows')

    # Compute summary stats
    rows = []
    for (gvkey, fdate), v in panel.items():
        shares = v['shares_per_mgr']
        if not shares:
            continue
        total = sum(shares.values())
        n_mgr = len(shares)
        hhi = sum((s / total) ** 2 for s in shares.values()) if total > 0 else float('nan')
        top5 = sum(sorted(shares.values(), reverse=True)[:5]) / total if total > 0 else float('nan')
        # Use median shrout across reporting managers (more robust to outliers)
        shrout_vals = sorted(v['shrout'])
        shrout = shrout_vals[len(shrout_vals) // 2] if shrout_vals else float('nan')
        inst_pct = total / shrout if shrout and shrout > 0 else float('nan')
        inst_pct_capped = min(inst_pct, 1.0) if inst_pct == inst_pct else float('nan')
        rows.append({
            'gvkey': gvkey,
            'fdate': fdate,
            'total_inst_shares': total,
            'n_managers': n_mgr,
            'hhi_managers': hhi,
            'top5_concentration': top5,
            'shrout_thousands': shrout,
            'inst_pct': inst_pct,
            'inst_pct_capped': inst_pct_capped,
        })

    rows.sort(key=lambda r: (r['gvkey'], r['fdate']))
    with open(OUT_PATH, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    print(f'\nWrote {OUT_PATH}')

    # ── Quick distribution diagnostics ──
    print('\n=== Distribution diagnostics ===')
    print(f'Total firm-quarters: {len(rows):,}')
    inst_pct_vals = [r['inst_pct_capped'] for r in rows
                     if r['inst_pct_capped'] == r['inst_pct_capped']]
    inst_pct_vals.sort()
    n = len(inst_pct_vals)
    print(f'inst_pct distribution (capped at 1.0):')
    print(f'  N      : {n:,}')
    if n > 0:
        print(f'  mean   : {sum(inst_pct_vals)/n:.3f}')
        print(f'  p10    : {inst_pct_vals[n//10]:.3f}')
        print(f'  p25    : {inst_pct_vals[n//4]:.3f}')
        print(f'  p50    : {inst_pct_vals[n//2]:.3f}')
        print(f'  p75    : {inst_pct_vals[3*n//4]:.3f}')
        print(f'  p90    : {inst_pct_vals[9*n//10]:.3f}')

    n_mgr_vals = sorted(r['n_managers'] for r in rows)
    if n_mgr_vals:
        print(f'\nn_managers distribution:')
        print(f'  median : {n_mgr_vals[len(n_mgr_vals)//2]}')
        print(f'  p10/p90: {n_mgr_vals[len(n_mgr_vals)//10]} / {n_mgr_vals[9*len(n_mgr_vals)//10]}')

    hhi_vals = sorted(r['hhi_managers'] for r in rows
                      if r['hhi_managers'] == r['hhi_managers'])
    if hhi_vals:
        print(f'\nhhi_managers distribution:')
        print(f'  p10/p50/p90: {hhi_vals[len(hhi_vals)//10]:.4f} / '
              f'{hhi_vals[len(hhi_vals)//2]:.4f} / '
              f'{hhi_vals[9*len(hhi_vals)//10]:.4f}')

    # Bimodality check on inst_pct
    if n > 100:
        n_below_50 = sum(1 for v in inst_pct_vals if v < 0.50)
        n_above_80 = sum(1 for v in inst_pct_vals if v >= 0.80)
        print(f'\nBimodality check on inst_pct:')
        print(f'  inst_pct < 0.50 : {n_below_50} ({100*n_below_50/n:.1f}%)')
        print(f'  inst_pct >= 0.80: {n_above_80} ({100*n_above_80/n:.1f}%)')


if __name__ == '__main__':
    main()
