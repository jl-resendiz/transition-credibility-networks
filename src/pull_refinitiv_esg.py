"""Pull Refinitiv ESG / emissions data for mapped firms via Eikon API.

Outputs: refinitiv_esg.csv (gvkey, ric, isin, scope1/2/3, env_score, co2_to_revenue, trbc, date/period)
"""
import csv
import os
import sys
import time
from datetime import datetime

try:
    import eikon as ek
except ImportError:
    print("Missing eikon package. Install with `pip install eikon`.")
    sys.exit(1)

from _paths import raw_path, derived_path

MAP_PATH = derived_path("mappings", "gvkey_ric_map.csv")
OUT_PATH = raw_path("refinitiv", "refinitiv_esg.csv")

APP_KEY = os.getenv("EIKON_APP_KEY") or os.getenv("REFINITIV_APP_KEY")
if not APP_KEY:
    print("Set EIKON_APP_KEY (or REFINITIV_APP_KEY) in your environment.")
    sys.exit(1)
ek.set_app_key(APP_KEY)


def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


# Field codes that are confirmed to work with the Eikon Data API.
# Column names returned by get_data are display names, not field codes,
# so we map concept -> (field_code, display_name).
FIELDS = {
    "scope1":         ("TR.CO2DirectScope1",        "CO2 Equivalent Emissions Direct, Scope 1"),
    "scope2":         ("TR.CO2IndirectScope2",       "CO2 Equivalent Emissions Indirect, Scope 2"),
    "scope3":         ("TR.CO2IndirectScope3",       "CO2 Equivalent Emissions Indirect, Scope 3"),
    "co2_to_revenue": ("TR.AnalyticCO2",            "Total CO2 Equivalent Emissions To Revenues USD in million"),
    "env_score":      ("TR.EnvironmentPillarScore",  "Environmental Pillar Score"),
    "trbc":           ("TR.TRBCIndustry",            "TRBC Industry Name"),
}


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

    field_codes = [v[0] for v in FIELDS.values()]
    display_to_concept = {v[1]: k for k, v in FIELDS.items()}

    print("Fields to pull:")
    for concept, (code, display) in FIELDS.items():
        print(f"  {concept}: {code} -> {display}")

    out_rows = []
    rics = [r["ric"] for r in rows]
    ric_to_meta = {r["ric"]: r for r in rows}

    # Pull in chunks to respect API limits
    for i, chunk in enumerate(chunked(rics, 50)):
        try:
            df, err = ek.get_data(chunk, field_codes)
        except Exception as e:
            print(f"  Chunk {i+1} exception: {e}")
            time.sleep(1)
            continue
        if df is None or df.empty:
            continue

        for _, row in df.iterrows():
            ric = row.get("Instrument")
            meta = ric_to_meta.get(ric, {})
            out = {
                "gvkey": meta.get("gvkey"),
                "isin": meta.get("isin"),
                "ric": ric,
            }
            # Map display-name columns back to standardized concept names
            for col in row.index:
                if col in display_to_concept:
                    val = row[col]
                    try:
                        import pandas as pd
                        if pd.isna(val):
                            val = ""
                    except (TypeError, ValueError):
                        pass
                    out[display_to_concept[col]] = val
            out_rows.append(out)
        time.sleep(0.3)
        if (i + 1) % 5 == 0:
            print(f"  Progress: {min((i+1)*50, len(rics))}/{len(rics)} RICs, {len(out_rows)} rows")

    # Write output
    fieldnames = ["gvkey", "isin", "ric"] + list(FIELDS.keys())
    with open(OUT_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    # Summary
    has_scope1 = sum(1 for r in out_rows if r.get("scope1"))
    has_co2rev = sum(1 for r in out_rows if r.get("co2_to_revenue"))
    print(f"\nSaved {len(out_rows)} rows to {OUT_PATH}")
    print(f"  With Scope 1 data: {has_scope1}")
    print(f"  With CO2/Revenue: {has_co2rev}")


if __name__ == "__main__":
    main()
