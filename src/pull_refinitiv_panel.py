"""Pull historical annual ESG/emissions panel from Refinitiv (2014-2024).

Pulls one year at a time for all mapped firms, ensuring reliable year labelling.
This is the key data upgrade: alpha varies by firm-year rather than being a
static cross-sectional snapshot.

Outputs: refinitiv_panel.csv
  (gvkey, isin, ric, year, esg_score, env_score, scope1, scope2, scope3,
   co2_total, co2_to_revenue)
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
    print("Missing packages. Install with `pip install eikon pandas`.")
    sys.exit(1)

MAP_PATH = derived_path("mappings", "gvkey_ric_map.csv")
OUT_PATH = raw_path("refinitiv", "refinitiv_panel.csv")

APP_KEY = os.getenv("EIKON_APP_KEY") or os.getenv("REFINITIV_APP_KEY")
if not APP_KEY:
    print("Set EIKON_APP_KEY (or REFINITIV_APP_KEY) in your environment.")
    sys.exit(1)
ek.set_app_key(APP_KEY)

YEARS = list(range(2014, 2025))  # 2014-2024

# Fields to pull (concept -> display name returned by Eikon)
FIELDS = {
    "esg_score":      "ESG Score",
    "env_score":      "Environmental Pillar Score",
    "scope1":         "CO2 Equivalent Emissions Direct, Scope 1",
    "scope2":         "CO2 Equivalent Emissions Indirect, Scope 2",
    "scope3":         "CO2 Equivalent Emissions Indirect, Scope 3",
    "co2_total":      "CO2 Emission Total",
    "co2_to_revenue": "Total CO2 Equivalent Emissions To Revenues USD in million",
}

# Field codes to request
FIELD_CODES = [
    "TR.TRESGScore",
    "TR.EnvironmentPillarScore",
    "TR.CO2DirectScope1",
    "TR.CO2IndirectScope2",
    "TR.CO2IndirectScope3",
    "TR.CO2EmissionTotal",
    "TR.AnalyticCO2",
]

# display name -> concept
DISPLAY_TO_CONCEPT = {v: k for k, v in FIELDS.items()}


def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


def load_mapping():
    rows = []
    with open(MAP_PATH, "r", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            ric = row.get("ric") or row.get("RIC") or ""
            if not ric:
                continue
            rows.append({
                "gvkey": row.get("gvkey"),
                "isin": row.get("isin"),
                "ric": ric,
            })
    return rows


def main():
    rows = load_mapping()
    if not rows:
        print("No RICs found in gvkey_ric_map.csv")
        sys.exit(1)

    print(f"Pulling annual panel ({YEARS[0]}-{YEARS[-1]}) for {len(rows)} firms")
    print(f"Fields: {list(FIELDS.keys())}")

    rics = [r["ric"] for r in rows]
    ric_to_meta = {r["ric"]: r for r in rows}
    out_rows = []
    total_calls = 0

    for year in YEARS:
        year_params = {
            'SDate': f'{year}-01-01',
            'EDate': f'{year}-12-31',
            'Frq': 'FY',
        }
        year_rows = 0

        for chunk in chunked(rics, 50):
            try:
                df, err = ek.get_data(chunk, FIELD_CODES,
                                      parameters=year_params)
            except Exception as e:
                print(f"  {year} exception: {e}")
                time.sleep(2)
                continue
            total_calls += 1

            if df is None or df.empty:
                continue

            for _, row in df.iterrows():
                ric = row.get("Instrument")
                meta = ric_to_meta.get(ric, {})

                # Check if row has any actual data
                has_data = False
                out = {
                    "gvkey": meta.get("gvkey"),
                    "isin": meta.get("isin"),
                    "ric": ric,
                    "year": year,
                }
                for col in row.index:
                    concept = DISPLAY_TO_CONCEPT.get(col)
                    if concept:
                        val = row[col]
                        try:
                            if pd.isna(val):
                                val = ""
                            else:
                                has_data = True
                        except (TypeError, ValueError):
                            if val not in (None, ""):
                                has_data = True
                        out[concept] = val

                # Only keep rows with at least one populated field
                if has_data:
                    out_rows.append(out)
                    year_rows += 1

            time.sleep(0.3)

        print(f"  {year}: {year_rows} firms with data ({total_calls} API calls so far)")

    # Write output
    fieldnames = ["gvkey", "isin", "ric", "year"] + list(FIELDS.keys())
    with open(OUT_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    # Summary
    n_firms = len(set(r.get("gvkey") for r in out_rows if r.get("gvkey")))
    years_covered = sorted(set(r.get("year") for r in out_rows if r.get("year")))
    has_co2rev = sum(1 for r in out_rows if r.get("co2_to_revenue") not in ("", None))
    has_scope1 = sum(1 for r in out_rows if r.get("scope1") not in ("", None))
    has_esg = sum(1 for r in out_rows if r.get("esg_score") not in ("", None))
    print(f"\nSaved {len(out_rows)} rows to {OUT_PATH}")
    print(f"  Firms: {n_firms}")
    if years_covered:
        print(f"  Years: {years_covered[0]}-{years_covered[-1]} ({len(years_covered)} unique)")
    print(f"  With ESG Score: {has_esg}")
    print(f"  With CO2/Revenue: {has_co2rev}")
    print(f"  With Scope 1: {has_scope1}")


if __name__ == "__main__":
    main()
