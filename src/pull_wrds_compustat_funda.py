"""Pull Compustat annual fundamentals for US-linked panel firms.

Used to compute book-equity for the DGTW-style B/M sort:
  book_equity = ceq + txditc - pstk
  (Davis-Fama-French 2000 / Daniel-Titman 1997 standard definition)

Output: data/raw/wrds/compustat_funda_us.csv
"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import _credentials

import pandas as pd
import wrds

from _paths import raw_path

OUT_DIR = raw_path('wrds')


def main():
    ccm = pd.read_csv(os.path.join(OUT_DIR, 'ccm_link_us.csv'),
                      dtype={'gvkey': str})
    panel_gvkeys = ccm['gvkey'].astype(str).str.zfill(6).unique().tolist()
    print(f'US-linked gvkeys: {len(panel_gvkeys)}')

    user = os.environ['WRDS_USERNAME']
    pwd = os.environ['WRDS_PASSWORD']
    db = wrds.Connection(wrds_username=user, wrds_password=pwd)

    gvkey_list = ','.join(f"'{g}'" for g in panel_gvkeys)
    print('Pulling comp.funda...')
    df = db.raw_sql(
        f"""SELECT gvkey, datadate, fyear, fyr, ceq, txditc, pstk, pstkrv,
                   pstkl, at, seq, ni, dvc, dvp
            FROM comp.funda
            WHERE gvkey IN ({gvkey_list})
              AND datafmt = 'STD' AND popsrc = 'D' AND consol = 'C'
              AND indfmt = 'INDL'
              AND fyear BETWEEN 2007 AND 2024""",
    )
    print(f'  {len(df):,} firm-year rows ({df["gvkey"].nunique()} unique gvkeys)')

    # Compute book equity (Davis-Fama-French)
    # BE = CEQ + TXDITC - PSTK_pref
    # where PSTK_pref = PSTKRV (preferred stock at redemption value) if available,
    # else PSTKL (liquidation), else PSTK (par value).
    df['pstk_pref'] = df['pstkrv'].fillna(df['pstkl']).fillna(df['pstk']).fillna(0)
    df['txditc'] = df['txditc'].fillna(0)
    df['book_equity'] = df['ceq'] + df['txditc'] - df['pstk_pref']

    out_path = os.path.join(OUT_DIR, 'compustat_funda_us.csv')
    df.to_csv(out_path, index=False)
    print(f'Wrote {out_path}')

    # Coverage check
    has_be = df['book_equity'].notna().sum()
    print(f'\nBook-equity populated: {has_be} of {len(df)} ({100*has_be/len(df):.1f}%)')
    print(f'fyear range: {df["fyear"].min():.0f} to {df["fyear"].max():.0f}')

    db.close()
    print('Done.')


if __name__ == '__main__':
    main()
