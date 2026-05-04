"""Download Ken French momentum (UMD) factor: monthly + daily.

Source: https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/data_library.html

Outputs:
  - data/raw/factors/F-F_Momentum_Factor.csv
  - data/raw/factors/F-F_Momentum_Factor_daily.csv

Both stored as plain CSV with columns: date, Mom (in percent, matching FF3 format).
"""
import io
import os
import sys
import urllib.request
import zipfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from _paths import raw_path

URLS = {
    'F-F_Momentum_Factor.csv': (
        'https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/'
        'ftp/F-F_Momentum_Factor_CSV.zip'
    ),
    'F-F_Momentum_Factor_daily.csv': (
        'https://mba.tuck.dartmouth.edu/pages/faculty/ken.french/'
        'ftp/F-F_Momentum_Factor_daily_CSV.zip'
    ),
}

OUT_DIR = raw_path('factors')
os.makedirs(OUT_DIR, exist_ok=True)


def download_one(filename, url):
    out_path = os.path.join(OUT_DIR, filename)
    if os.path.exists(out_path):
        print(f'{filename}: already exists, skipping.')
        return
    print(f'Downloading {url}')
    req = urllib.request.Request(url, headers={'User-Agent': 'Mozilla/5.0'})
    with urllib.request.urlopen(req, timeout=60) as r:
        data = r.read()
    print(f'  Got {len(data):,} bytes (zipped).')
    z = zipfile.ZipFile(io.BytesIO(data))
    names = z.namelist()
    csv_name = next((n for n in names if n.endswith('.CSV') or n.endswith('.csv')),
                    None)
    if csv_name is None:
        raise RuntimeError(f'No CSV inside zip: {names}')
    extracted = z.read(csv_name).decode('latin-1')
    with open(out_path, 'w', encoding='utf-8', newline='') as f:
        f.write(extracted)
    print(f'  Wrote {out_path} ({os.path.getsize(out_path):,} bytes)')


def main():
    for filename, url in URLS.items():
        download_one(filename, url)
    print('\nDone.')


if __name__ == '__main__':
    main()
