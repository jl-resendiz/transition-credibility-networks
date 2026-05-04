"""Pull Thomson Reuters 13F institutional holdings for US-linked panel firms.

Inputs:
  - data/raw/wrds/ccm_link_us.csv  (from pull_wrds_crsp_daily.py)
  - data/raw/wrds/crsp_daily_us.csv (for cusip lookup)

Outputs:
  - data/raw/wrds/holdings_13f.csv (raw quarterly holdings, one row per
    manager-firm-quarter)
  - data/raw/wrds/holdings_13f_aggregated.csv (one row per firm-quarter,
    aggregated across managers)

The aggregation collapses managers and reports:
  - total_inst_shares
  - n_managers
  - hhi_managers   (Herfindahl of share-of-shares-held by manager)
  - top5_concentration

Date range: 2008-01-01 to 2025-09-30 (matches s34 coverage end).
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _credentials

import pandas as pd
import wrds

from _paths import raw_path

OUT_DIR = raw_path('wrds')
START_DATE = '2008-01-01'
END_DATE = '2025-09-30'


def main():
    # Read panel firm cusips (from CRSP daily pull)
    crsp_path = os.path.join(OUT_DIR, 'crsp_daily_us.csv')
    if not os.path.exists(crsp_path):
        sys.exit(f'Run pull_wrds_crsp_daily.py first: {crsp_path} missing.')
    df_crsp = pd.read_csv(crsp_path, usecols=['gvkey', 'permno', 'cusip', 'date'])
    df_crsp = df_crsp.drop_duplicates(['gvkey', 'cusip']).dropna(subset=['cusip'])

    # Get unique cusip-permno-gvkey triplets seen in CRSP daily
    cusip_map = df_crsp.groupby('cusip')['gvkey'].first().reset_index()
    cusip_map['gvkey'] = cusip_map['gvkey'].astype(str).str.zfill(6)
    cusips = cusip_map['cusip'].astype(str).str.zfill(8).tolist()
    print(f'CUSIPs to query: {len(cusips)} (from {df_crsp["gvkey"].nunique()} gvkeys)')

    cusip_list = ','.join(f"'{c}'" for c in cusips)

    user = os.environ['WRDS_USERNAME']
    pwd = os.environ['WRDS_PASSWORD']
    print(f'Connecting to WRDS as {user!r}...')
    db = wrds.Connection(wrds_username=user, wrds_password=pwd)

    print('\n[1/2] Pulling tr_13f.s34 holdings...')
    holdings = db.raw_sql(
        f"""SELECT fdate, cusip, mgrno, shares, sole, shared, prc, shrout1, shrout2
            FROM tr_13f.s34
            WHERE cusip IN ({cusip_list})
              AND fdate BETWEEN '{START_DATE}' AND '{END_DATE}'
              AND shares > 0""",
    )
    print(f'  {len(holdings):,} raw holding rows.')

    # Map cusip → gvkey
    holdings['cusip'] = holdings['cusip'].astype(str).str.zfill(8)
    holdings = holdings.merge(cusip_map, on='cusip', how='left')

    print(f'  Unique fdates: {holdings["fdate"].nunique()}')
    print(f'  Unique gvkeys covered: {holdings["gvkey"].nunique()}')
    print(f'  Date range: {holdings["fdate"].min()} to {holdings["fdate"].max()}')

    # Save raw
    raw_out = os.path.join(OUT_DIR, 'holdings_13f.csv')
    holdings.to_csv(raw_out, index=False)
    print(f'  Wrote {raw_out} ({os.path.getsize(raw_out)/1e6:.1f} MB)')

    # ── Aggregate by (gvkey, fdate) ──
    print('\n[2/2] Aggregating by (gvkey, fdate)...')
    holdings = holdings.dropna(subset=['gvkey'])

    def agg(grp):
        total = grp['shares'].sum()
        n_mgr = grp['mgrno'].nunique()
        # HHI of manager share-of-shares-held
        if total > 0:
            shares_per_mgr = grp.groupby('mgrno')['shares'].sum()
            wts = shares_per_mgr / total
            hhi = float((wts ** 2).sum())
            top5 = float(shares_per_mgr.nlargest(5).sum() / total)
        else:
            hhi = float('nan')
            top5 = float('nan')
        return pd.Series({
            'total_inst_shares': total,
            'n_managers': n_mgr,
            'hhi_managers': hhi,
            'top5_concentration': top5,
        })

    agg_df = holdings.groupby(['gvkey', 'fdate']).apply(
        agg, include_groups=False).reset_index()
    print(f'  {len(agg_df):,} (gvkey × quarter) rows.')

    agg_out = os.path.join(OUT_DIR, 'holdings_13f_aggregated.csv')
    agg_df.to_csv(agg_out, index=False)
    print(f'  Wrote {agg_out}')

    # Quick sanity check: ownership % requires total shares outstanding from CRSP
    # We compute it in build_institutional_panel.py (Phase 2.1).
    db.close()
    print('\nDone.')


if __name__ == '__main__':
    main()
