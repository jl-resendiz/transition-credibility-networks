"""Convert GEM xlsx tracker files to CSV for faster loading.

Reads the four GEM xlsx files from data/raw/gem/ and writes CSV versions
to data/derived/gem/, preserving all columns and rows.  This is a
pre-processing step: downstream scripts can then read lightweight CSVs
instead of parsing multi-MB Excel workbooks on every run.

Usage:
    python src/build_gem_csv.py
"""

import csv
import os
import time

import openpyxl

from _paths import raw_path, derived_path

# ── Configuration ────────────────────────────────────────────────────

TRACKERS = [
    # (xlsx filename, sheet name, output csv filename)
    ('Global-Coal-Plant-Tracker-January-2026.xlsx',
     'Units',
     'gem_coal.csv'),
    ('Global-Oil-and-Gas-Plant-Tracker-GOGPT-January-2026.xlsx',
     'Gas & Oil Units',
     'gem_gas.csv'),
    ('Global-Solar-Power-Tracker-February-2026.xlsx',
     'Utility-Scale (1 MW+)',
     'gem_solar.csv'),
    ('Global-Wind-Power-Tracker-February-2026.xlsx',
     'Data',
     'gem_wind.csv'),
]

OUT_DIR = derived_path('gem')

# ── Main ─────────────────────────────────────────────────────────────


def convert_one(xlsx_fname, sheet_name, csv_fname):
    """Read *sheet_name* from an xlsx file and write all rows to CSV."""
    fpath = raw_path('gem', xlsx_fname)
    if not os.path.exists(fpath):
        print(f'  SKIP: {xlsx_fname} not found')
        return

    t0 = time.perf_counter()
    wb = openpyxl.load_workbook(fpath, read_only=True)
    ws = wb[sheet_name]

    rows = ws.iter_rows(values_only=True)
    headers = [str(h) if h is not None else '' for h in next(rows)]

    out_path = os.path.join(OUT_DIR, csv_fname)
    n_rows = 0
    with open(out_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(headers)
        for row in rows:
            writer.writerow([str(v) if v is not None else '' for v in row])
            n_rows += 1

    wb.close()
    elapsed = time.perf_counter() - t0
    print(f'  {xlsx_fname}')
    print(f'    -> {csv_fname}  ({n_rows:,} rows, {len(headers)} cols, '
          f'{elapsed:.1f}s)')


def main():
    os.makedirs(OUT_DIR, exist_ok=True)
    print('Converting GEM xlsx -> CSV ...')
    t_total = time.perf_counter()

    for xlsx_fname, sheet_name, csv_fname in TRACKERS:
        convert_one(xlsx_fname, sheet_name, csv_fname)

    elapsed = time.perf_counter() - t_total
    print(f'\nDone.  Total time: {elapsed:.1f}s')
    print(f'Output directory: {OUT_DIR}')


if __name__ == '__main__':
    main()
