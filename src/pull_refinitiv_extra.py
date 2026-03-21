"""Pull additional cross-sectional ESG/governance fields from Refinitiv.

Supplements the core ESG/emissions pull with sub-pillar scores, energy use,
ownership structure, and policy indicators.

Pulls in two batches because certain ESG sub-scores (Emissions, Resource Use,
Environmental Innovation) use different field codes (TR.TRESG* prefix) and
get silently dropped if mixed with other field families in one API call.

Outputs: refinitiv_extra.csv
  (gvkey, isin, ric, esg_score, social_score, governance_score,
   emissions_score, resource_use_score, env_innovation_score, energy_use,
   free_float_pct, shares_outstanding, policy_emissions, policy_energy_eff)
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
OUT_PATH = raw_path("refinitiv", "refinitiv_extra.csv")

APP_KEY = os.getenv("EIKON_APP_KEY") or os.getenv("REFINITIV_APP_KEY")
if not APP_KEY:
    print("Set EIKON_APP_KEY (or REFINITIV_APP_KEY) in your environment.")
    sys.exit(1)
ek.set_app_key(APP_KEY)


def chunked(seq, n):
    for i in range(0, len(seq), n):
        yield seq[i:i + n]


# Pull in two batches to avoid field-code conflicts
# Batch A: ESG sub-scores (TR.TRESG* prefix required)
BATCH_A = {
    "emissions_score":      ("TR.TRESGEmissionsScore",     "Emissions Score"),
    "resource_use_score":   ("TR.TRESGResourceUseScore",   "Resource Use Score"),
    "env_innovation_score": ("TR.TRESGInnovationScore",    "Environmental Innovation Score"),
}

# Batch B: pillar scores + operational/ownership
BATCH_B = {
    "esg_score":            ("TR.TRESGScore",              "ESG Score"),
    "social_score":         ("TR.SocialPillarScore",       "Social Pillar Score"),
    "governance_score":     ("TR.GovernancePillarScore",   "Governance Pillar Score"),
    "energy_use":           ("TR.EnergyUseTotal",          "Energy Use Total"),
    "free_float_pct":       ("TR.FreeFloatPct",            "Free Float (Percent)"),
    "shares_outstanding":   ("TR.SharesOutstanding",       "Outstanding Shares"),
    "policy_emissions":     ("TR.PolicyEmissions",         "Policy Emissions"),
    "policy_energy_eff":    ("TR.PolicyEnergyEfficiency",  "Policy Energy Efficiency"),
}

ALL_CONCEPTS = list(BATCH_A.keys()) + list(BATCH_B.keys())


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


def pull_batch(rics, ric_to_meta, fields_dict, batch_label):
    """Pull one batch of fields for all RICs. Returns dict keyed by RIC."""
    field_codes = [v[0] for v in fields_dict.values()]
    display_to_concept = {v[1]: k for k, v in fields_dict.items()}

    print(f"\n  {batch_label}: {list(fields_dict.keys())}")
    results = {}  # ric -> {concept: value}
    n_chunks = (len(rics) + 49) // 50

    for i, chunk in enumerate(chunked(rics, 50)):
        try:
            df, err = ek.get_data(chunk, field_codes)
        except Exception as e:
            print(f"    Chunk {i+1}/{n_chunks} exception: {e}")
            time.sleep(1)
            continue
        if df is None or df.empty:
            continue

        for _, row in df.iterrows():
            ric = row.get("Instrument")
            if not ric:
                continue
            data = {}
            for col in row.index:
                concept = display_to_concept.get(col)
                if concept:
                    val = row[col]
                    try:
                        if pd.isna(val):
                            val = ""
                    except (TypeError, ValueError):
                        pass
                    data[concept] = val
            results[ric] = data

        time.sleep(0.3)
        if (i + 1) % 5 == 0:
            print(f"    Progress: {min((i+1)*50, len(rics))}/{len(rics)} RICs")

    return results


def main():
    rows = load_mapping()
    if not rows:
        print("No RICs found in gvkey_ric_map.csv")
        sys.exit(1)

    rics = [r["ric"] for r in rows]
    ric_to_meta = {r["ric"]: r for r in rows}

    print(f"Pulling extra fields for {len(rows)} firms")

    # Pull each batch
    results_a = pull_batch(rics, ric_to_meta, BATCH_A, "Batch A (ESG sub-scores)")
    results_b = pull_batch(rics, ric_to_meta, BATCH_B, "Batch B (pillar + operational)")

    # Merge into output rows
    out_rows = []
    for r in rows:
        ric = r["ric"]
        out = {
            "gvkey": r["gvkey"],
            "isin": r["isin"],
            "ric": ric,
        }
        out.update(results_a.get(ric, {}))
        out.update(results_b.get(ric, {}))
        out_rows.append(out)

    # Write output
    fieldnames = ["gvkey", "isin", "ric"] + ALL_CONCEPTS
    with open(OUT_PATH, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    # Summary
    populated = {}
    for concept in ALL_CONCEPTS:
        populated[concept] = sum(1 for r in out_rows if r.get(concept) not in ("", None))
    print(f"\nSaved {len(out_rows)} rows to {OUT_PATH}")
    for concept, count in populated.items():
        print(f"  {concept}: {count}/{len(out_rows)} populated")


if __name__ == "__main__":
    main()
