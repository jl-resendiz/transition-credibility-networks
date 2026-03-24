"""Pull monthly total returns from Eikon for non-US firms.

Uses TR.TotalReturn (includes dividends, splits, corporate actions)
via ek.get_data() with monthly frequency and params dict.

Requires:
  - Refinitiv Workspace Desktop running locally
  - EIKON_APP_KEY environment variable set
  - pip install eikon pandas

Usage:
  python src/pull_eikon_returns.py

Outputs:
  data/raw/eikon/eikon_monthly_returns.csv
  Columns: gvkey, datadate, ret_monthly
"""
import csv
import os
import sys
import time

from _paths import raw_path, derived_path

try:
    import eikon as ek
    import pandas as pd
except ImportError:
    print("Missing packages. Install with: pip install eikon pandas")
    sys.exit(1)

# ── Configuration ──
APP_KEY = os.getenv("EIKON_APP_KEY") or os.getenv("REFINITIV_APP_KEY")
if not APP_KEY:
    print("Set EIKON_APP_KEY (or REFINITIV_APP_KEY) in your environment.")
    sys.exit(1)
ek.set_app_key(APP_KEY)

START_DATE = "2013-01-01"
END_DATE = "2026-01-31"
CHUNK_SIZE = 20
OUT_PATH = raw_path("eikon", "eikon_monthly_returns.csv")

# ── Load RIC mapping ──
MAP_PATH = derived_path("mappings", "gvkey_ric_map.csv")

ric_to_gvkey = {}
with open(MAP_PATH, "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        gk = row.get("gvkey", "")
        ric = row.get("ric", "")
        if gk and ric:
            ric_to_gvkey[ric] = gk

print(f"RICs to pull: {len(ric_to_gvkey)}")


def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


# ── Pull monthly total returns ──
all_results = []
rics = list(ric_to_gvkey.keys())
n_success = 0
n_fail = 0

print(f"\nPulling monthly total returns for {len(rics)} RICs in chunks of {CHUNK_SIZE}...")

for i, chunk in enumerate(chunked(rics, CHUNK_SIZE)):
    try:
        df, err = ek.get_data(
            chunk,
            ['TR.TotalReturn.date', 'TR.TotalReturn'],
            {'SDate': START_DATE, 'EDate': END_DATE, 'Frq': 'M'},
        )
    except Exception as e:
        n_fail += len(chunk)
        time.sleep(2)
        continue

    if df is None or df.empty:
        n_fail += len(chunk)
        time.sleep(0.5)
        continue

    chunk_ok = set()
    for _, row in df.iterrows():
        ric = row.get("Instrument")
        if not ric or ric not in ric_to_gvkey:
            continue

        date_val = row.get("Date")
        ret_val = row.get("Total Return")

        if pd.isna(date_val) or pd.isna(ret_val):
            continue

        # Parse date
        try:
            dt = pd.to_datetime(date_val)
            datadate = dt.strftime("%Y%m%d")
        except Exception:
            continue

        # Convert from percentage to decimal, cap at +/- 50%
        ret = float(ret_val) / 100.0
        ret = max(-0.5, min(0.5, ret))

        all_results.append({
            "gvkey": ric_to_gvkey[ric],
            "datadate": datadate,
            "ret_monthly": f"{ret:.6f}",
        })
        chunk_ok.add(ric)

    n_success += len(chunk_ok)
    n_fail += len(chunk) - len(chunk_ok)

    if (i + 1) % 3 == 0:
        print(f"  Progress: {min((i+1)*CHUNK_SIZE, len(rics))}/{len(rics)} RICs "
              f"({n_success} ok, {n_fail} failed, {len(all_results):,} obs)")

    time.sleep(0.5)

# ── Write output ──
os.makedirs(os.path.dirname(OUT_PATH), exist_ok=True)
with open(OUT_PATH, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=["gvkey", "datadate", "ret_monthly"])
    w.writeheader()
    w.writerows(all_results)

n_firms = len(set(r["gvkey"] for r in all_results))
print(f"\nTotal return observations: {len(all_results):,}")
print(f"Firms with data: {n_firms}")
print(f"Saved to {OUT_PATH}")
if all_results:
    dates = sorted(r["datadate"] for r in all_results)
    print(f"Date range: {dates[0]} to {dates[-1]}")
print("\nNext step: re-run compute_returns.py to incorporate these returns.")
