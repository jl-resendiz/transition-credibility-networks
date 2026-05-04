"""Build DGTW-style size × B/M × momentum characteristic deciles for US panel firms.

Following Daniel-Grinblatt-Titman-Wermers (1997) intent, but adapted to a
small (87-firm) cross-section. We use TERCILES instead of quintiles to avoid
thin-cell pathologies, and rank firms WITHIN their event-month (not vs the
full CRSP universe). The within-event ranking is appropriate because the paper's
identifying variation is cross-firm dispersion within an event, not absolute
factor exposure.

Inputs:
  - data/raw/wrds/crsp_daily_us.csv  (gvkey, date, prc, shrout, ret, retx)
  - data/raw/wrds/compustat_funda_us.csv  (gvkey, fyear, datadate, book_equity)

Output:
  - data/derived/dgtw_chars_us.csv
    cols: gvkey, month_end, log_size, log_bm, mom12_2, size_t, bm_t, mom_t

  Where:
    size_t   ∈ {1,2,3} = tercile of log_size at month_end across the US sample
    bm_t     ∈ {1,2,3} = tercile of log_bm   at month_end across the US sample
    mom_t    ∈ {1,2,3} = tercile of mom12_2  at month_end across the US sample
    log_size = log(prc × shrout) at month_end (market cap, $thousands)
    log_bm   = log(book_equity / market_value)
              with book_equity from most recent fiscal year ending >= 6mo before month_end
    mom12_2  = cumulative return over months [t-12, t-2]  (skip month t-1)
"""
import csv
import os
import sys
import math
from collections import defaultdict
from datetime import datetime

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _paths import raw_path, derived_path

CRSP_PATH = os.path.join(raw_path('wrds'), 'crsp_daily_us.csv')
FUNDA_PATH = os.path.join(raw_path('wrds'), 'compustat_funda_us.csv')
OUT_PATH = os.path.join(derived_path('dgtw'), 'dgtw_chars_us.csv')


def parse_date(s):
    return datetime.strptime(s[:10], '%Y-%m-%d').date()


def month_end(d):
    """Return last day of d's month as a date."""
    if d.month == 12:
        nxt = d.replace(year=d.year + 1, month=1, day=1)
    else:
        nxt = d.replace(month=d.month + 1, day=1)
    from datetime import timedelta
    return nxt - timedelta(days=1)


def main():
    os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)

    # ── 1) Read CRSP daily, build month-end (prc, shrout, monthly_ret) per firm ──
    print(f'Reading {CRSP_PATH}...')
    monthly = defaultdict(dict)  # (gvkey, ym) -> {prc, shrout, monthly_ret}
    daily_returns = defaultdict(list)  # (gvkey, ym) -> list of daily rets

    with open(CRSP_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        n = 0
        for row in reader:
            n += 1
            try:
                gvkey = str(row['gvkey']).strip().split('.')[0].zfill(6)
                if gvkey == 'NAN' or gvkey == '':
                    continue
                date = parse_date(row['date'])
                prc = float(row['prc']) if row['prc'] not in ('', 'nan') else None
                shrout = float(row['shrout']) if row['shrout'] not in ('', 'nan') else None
                ret = float(row['ret']) if row['ret'] not in ('', 'nan', 'B', 'C') else None
            except (ValueError, KeyError):
                continue
            ym = date.strftime('%Y-%m')
            key = (gvkey, ym)
            if ret is not None:
                daily_returns[key].append(ret)
            # Track most recent observation in the month for prc/shrout
            cur = monthly[key]
            if 'date' not in cur or date > cur['date']:
                cur['date'] = date
                cur['prc'] = abs(prc) if prc is not None else None  # CRSP signs prc neg if bid/ask midpoint
                cur['shrout'] = shrout
        print(f'  {n:,} daily rows scanned, {len(monthly):,} (firm, month) cells.')

    # Compute monthly returns from daily compounds
    for key, rets in daily_returns.items():
        compound = 1.0
        for r in rets:
            compound *= (1.0 + r)
        monthly[key]['monthly_ret'] = compound - 1.0

    # ── 2) Read Compustat fundamentals → build {gvkey: [(datadate, book_equity)]} ──
    print(f'\nReading {FUNDA_PATH}...')
    funda = defaultdict(list)
    with open(FUNDA_PATH, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            try:
                gvkey = str(row['gvkey']).strip().split('.')[0].zfill(6)
                datadate = parse_date(row['datadate'])
                be = float(row['book_equity']) if row['book_equity'] not in ('', 'nan') else None
            except (ValueError, KeyError):
                continue
            if be is not None and be > 0:
                funda[gvkey].append((datadate, be))
    for k in funda:
        funda[k].sort()
    print(f'  Compustat: {len(funda)} gvkeys with valid book equity.')

    def get_book_equity_at(gvkey, target_date):
        """Most recent fyear-end book equity at least 6 months before target."""
        history = funda.get(gvkey, [])
        from datetime import timedelta
        cutoff = target_date - timedelta(days=180)
        eligible = [(d, be) for (d, be) in history if d <= cutoff]
        return eligible[-1][1] if eligible else None

    # ── 3) For each (gvkey, month_end), compute log_size, log_bm, mom12_2 ──
    print('\nComputing characteristics per (firm, month)...')
    keys_sorted = sorted(monthly.keys())  # (gvkey, ym) tuples

    # Build per-firm chronologically-sorted list of (ym, monthly_ret) for momentum
    firm_returns = defaultdict(list)
    for key in keys_sorted:
        gvkey, ym = key
        if 'monthly_ret' in monthly[key]:
            firm_returns[gvkey].append((ym, monthly[key]['monthly_ret']))
    for k in firm_returns:
        firm_returns[k].sort()

    rows = []
    for (gvkey, ym), rec in monthly.items():
        if 'prc' not in rec or rec['prc'] is None:
            continue
        if 'shrout' not in rec or rec['shrout'] is None or rec['shrout'] <= 0:
            continue
        prc = rec['prc']
        shrout = rec['shrout']  # in thousands
        market_value = prc * shrout  # in $thousands

        if market_value <= 0:
            continue

        # Book-equity at this month-end
        me_date = month_end(rec['date'])
        be = get_book_equity_at(gvkey, me_date)

        # B/M (book equity in $millions, market value in $thousands)
        if be is None or be <= 0:
            log_bm = None
        else:
            # be is millions, market_value is thousands → ratio: be * 1000 / market_value
            log_bm = math.log((be * 1000.0) / market_value)

        log_size = math.log(market_value)

        # Momentum: cumulative ret over months [t-12, t-2]
        history = firm_returns[gvkey]
        idx = next((i for i, (y, _) in enumerate(history) if y == ym), None)
        if idx is None or idx < 12:
            mom = None
        else:
            window = history[idx - 12: idx - 1]  # months [t-12, t-2]
            if len(window) >= 11:  # require at least 11 of 11 months
                compound = 1.0
                for _, r in window:
                    compound *= (1.0 + r)
                mom = compound - 1.0
            else:
                mom = None

        rows.append({
            'gvkey': gvkey,
            'month_end': me_date.strftime('%Y-%m-%d'),
            'log_size': log_size,
            'log_bm': log_bm,
            'mom12_2': mom,
        })

    print(f'  {len(rows):,} (firm, month-end) characteristic rows.')

    # ── 4) Within-month tercile assignment for size, B/M, momentum ──
    by_month = defaultdict(list)
    for r in rows:
        by_month[r['month_end']].append(r)

    def assign_terciles(month_rows, key):
        vals = [(r[key], r) for r in month_rows if r[key] is not None]
        if len(vals) < 3:
            return  # too few firms; leave as None
        vals.sort()
        n = len(vals)
        for i, (v, r) in enumerate(vals):
            if i < n / 3:
                r[key + '_t'] = 1
            elif i < 2 * n / 3:
                r[key + '_t'] = 2
            else:
                r[key + '_t'] = 3

    for me, mrows in by_month.items():
        for k in ('log_size', 'log_bm', 'mom12_2'):
            assign_terciles(mrows, k)

    # Rename tercile keys for output clarity
    for r in rows:
        r['size_t'] = r.pop('log_size_t', None)
        r['bm_t'] = r.pop('log_bm_t', None)
        r['mom_t'] = r.pop('mom12_2_t', None)

    # ── 5) Write output ──
    rows.sort(key=lambda r: (r['gvkey'], r['month_end']))
    fieldnames = ['gvkey', 'month_end', 'log_size', 'log_bm', 'mom12_2',
                  'size_t', 'bm_t', 'mom_t']
    with open(OUT_PATH, 'w', encoding='utf-8', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        w.writerows(rows)
    print(f'\nWrote {OUT_PATH}')

    # ── 6) Diagnostics ──
    print('\n=== Diagnostics ===')
    print(f'Total (firm, month) rows: {len(rows):,}')
    print(f'Unique firms: {len(set(r["gvkey"] for r in rows))}')
    print(f'Unique months: {len(by_month)}')
    n_size = sum(1 for r in rows if r['size_t'] is not None)
    n_bm = sum(1 for r in rows if r['bm_t'] is not None)
    n_mom = sum(1 for r in rows if r['mom_t'] is not None)
    n_all = sum(1 for r in rows if r['size_t'] is not None and r['bm_t'] is not None and r['mom_t'] is not None)
    print(f'Firms with size_t : {n_size:,} ({100*n_size/len(rows):.1f}%)')
    print(f'Firms with bm_t   : {n_bm:,} ({100*n_bm/len(rows):.1f}%)')
    print(f'Firms with mom_t  : {n_mom:,} ({100*n_mom/len(rows):.1f}%)')
    print(f'Firms with all 3  : {n_all:,} ({100*n_all/len(rows):.1f}%)')


if __name__ == '__main__':
    main()
