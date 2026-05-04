"""Pull CRSP daily returns + market index for US-linked firms in the panel.

Inputs (from existing repo state):
  - results/summaries/panel_facts.json (gvkey universe)
  - WRDS credentials in .env

Outputs:
  - data/raw/wrds/crsp_daily_us.csv
      cols: gvkey, permno, date, ret, retx, prc, shrout, vol, cusip, cfacshr
  - data/raw/wrds/crsp_dsi.csv
      cols: date, vwretd, ewretd, sprtrn   (CRSP US market index)
  - data/raw/wrds/ccm_link_us.csv
      cols: gvkey, permno, linkdt, linkenddt, linktype, linkprim
      (panel firms that map cleanly to a US permno)

The CCM link uses LU/LC/LS link types (canonical filter for CRSP-Compustat
research) and primary-link rows only (linkprim ∈ {P, C}). Date range covers
2008-2026 to give estimation windows for events 2010+.
"""
import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _credentials  # auto-loads .env

import wrds
import pandas as pd

from _paths import raw_path, results_path

OUT_DIR = raw_path('wrds')
os.makedirs(OUT_DIR, exist_ok=True)

START_DATE = '2008-01-01'
END_DATE = '2026-02-28'


def main():
    # Read panel firm gvkeys
    facts_path = results_path('summaries', 'panel_facts.json')
    with open(facts_path) as f:
        facts = json.load(f)
    panel_gvkeys = [g.zfill(6) for g in facts['panel_firms_gvkeys']]
    print(f'Panel gvkeys: {len(panel_gvkeys)}')

    # Connect
    user = os.environ['WRDS_USERNAME']
    pwd = os.environ['WRDS_PASSWORD']
    print(f'Connecting to WRDS as {user!r}...')
    db = wrds.Connection(wrds_username=user, wrds_password=pwd)

    # ── 1) CCM link table → resolve gvkey to (permno, valid date range) ──
    gvkey_list = ','.join(f"'{g}'" for g in panel_gvkeys)
    print('\n[1/3] Loading CCM link table...')
    ccm = db.raw_sql(
        f"""SELECT gvkey, lpermno AS permno, linkdt, linkenddt, linktype, linkprim
            FROM crsp.ccmxpf_lnkhist
            WHERE gvkey IN ({gvkey_list})
              AND linktype IN ('LU', 'LC', 'LS')
              AND linkprim IN ('P', 'C')""",
    )
    ccm['linkdt'] = pd.to_datetime(ccm['linkdt'], errors='coerce')
    ccm['linkenddt'] = pd.to_datetime(ccm['linkenddt'], errors='coerce')
    ccm['linkenddt'] = ccm['linkenddt'].fillna(pd.Timestamp(END_DATE))
    ccm = ccm.dropna(subset=['permno']).reset_index(drop=True)
    ccm['permno'] = ccm['permno'].astype(int)
    print(f'  {len(ccm)} link rows for {ccm["gvkey"].nunique()} unique gvkeys / '
          f'{ccm["permno"].nunique()} unique permnos.')
    ccm.to_csv(os.path.join(OUT_DIR, 'ccm_link_us.csv'), index=False)

    permnos = sorted(ccm['permno'].unique().tolist())
    permno_list_sql = ','.join(str(p) for p in permnos)
    print(f'  Querying CRSP for {len(permnos)} permnos.')

    # ── 2) Daily stock file ──
    print('\n[2/3] Loading crsp.dsf (daily stock file)...')
    dsf = db.raw_sql(
        f"""SELECT permno, date, ret, retx, prc, vol, shrout, cfacshr, cusip
            FROM crsp.dsf
            WHERE permno IN ({permno_list_sql})
              AND date BETWEEN '{START_DATE}' AND '{END_DATE}'""",
    )
    print(f'  {len(dsf):,} daily rows pulled.')

    # Map permno → gvkey (using time-aware merge)
    dsf['date'] = pd.to_datetime(dsf['date'])
    dsf = dsf.merge(ccm[['permno', 'gvkey', 'linkdt', 'linkenddt']],
                    on='permno', how='left')
    in_range = (dsf['date'] >= dsf['linkdt']) & (dsf['date'] <= dsf['linkenddt'])
    dsf = dsf[in_range].drop(columns=['linkdt', 'linkenddt'])
    dsf = dsf.sort_values(['gvkey', 'date']).reset_index(drop=True)
    print(f'  After CCM date-range filter: {len(dsf):,} rows.')

    out_path = os.path.join(OUT_DIR, 'crsp_daily_us.csv')
    dsf.to_csv(out_path, index=False)
    print(f'  Wrote {out_path}')

    # ── 3) Daily index (vwretd, ewretd) for US market-adjusted CAR ──
    print('\n[3/3] Loading crsp.dsi (daily index)...')
    dsi = db.raw_sql(
        f"""SELECT date, vwretd, ewretd, vwretx, sprtrn
            FROM crsp.dsi
            WHERE date BETWEEN '{START_DATE}' AND '{END_DATE}'""",
    )
    print(f'  {len(dsi):,} index rows.')
    dsi.to_csv(os.path.join(OUT_DIR, 'crsp_dsi.csv'), index=False)

    db.close()
    print('\nDone.')


if __name__ == '__main__':
    main()
