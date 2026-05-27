"""
02 — Product Dictionary Build: Match DAEN records to pre-specified product sets

Reads the merged DAEN dataset and the pre-specified product dictionary.
For each study drug, identifies matching DAEN records using case-insensitive
substring matching against the search terms defined in the pre-registration.

Outputs:
  - data/processed/product_dictionary.csv        Final matched dictionary
  - data/processed/daen_study_drugs.parquet       DAEN records matching study drugs
  - data/processed/product_match_report.txt       QC gate report

Pre-registration QC gate (Section 9.2.1):
  - Match rate >= 95% for study drugs
  - Unmatched fraction reported

Usage:
    python scripts/03_build_product_dict.py
"""

import re
import sys
from pathlib import Path

import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
DAEN_PATH = PROCESSED_DIR / "daen_merged.parquet"
PRESPEC_PATH = PROCESSED_DIR / "product_dictionary_prespec.csv"


def load_prespec():
    df = pd.read_csv(PRESPEC_PATH)
    products = []
    for _, row in df.iterrows():
        search_terms = [
            t.strip().lower()
            for t in str(row["search_terms"]).split(";")
            if t.strip()
        ]
        products.append({
            "product_set": row["product_set"],
            "product_id": row["product_id"],
            "canonical_name": row["canonical_name"],
            "active_ingredient": row["active_ingredient"],
            "atc_code": row["atc_code"],
            "search_terms": search_terms,
            "indication_expected_pts": str(row.get("indication_expected_pts", "")),
        })
    return products


def match_records(daen: pd.DataFrame, products: list[dict]) -> pd.DataFrame:
    """Match DAEN records to study drugs using pre-specified search terms."""
    # Build a searchable text column from medicines and active_ingredient fields
    search_cols = []
    for col in ("medicines", "active_ingredient"):
        if col in daen.columns:
            search_cols.append(col)

    if not search_cols:
        print("ERROR: No 'medicines' or 'active_ingredient' column found.")
        sys.exit(1)

    daen["_search_text"] = ""
    for col in search_cols:
        daen["_search_text"] += " " + daen[col].fillna("").astype(str).str.lower()

    matches = []
    for prod in products:
        pattern = "|".join(re.escape(t) for t in prod["search_terms"] if t)
        if not pattern:
            continue

        mask = daen["_search_text"].str.contains(pattern, regex=True, na=False)
        matched = daen[mask].copy()
        matched["product_id"] = prod["product_id"]
        matched["product_set"] = prod["product_set"]
        matched["canonical_name"] = prod["canonical_name"]
        matched["active_ingredient_canonical"] = prod["active_ingredient"]
        matched["atc_code"] = prod["atc_code"]
        matches.append(matched)

        print(f"  {prod['canonical_name']:35s}  {mask.sum():>6,} reports matched")

    daen.drop(columns=["_search_text"], inplace=True)

    if not matches:
        return pd.DataFrame()
    return pd.concat(matches, ignore_index=True)


def main():
    print("=" * 70)
    print("  Product Dictionary Build — DAEN Tropical PV")
    print("=" * 70)

    if not DAEN_PATH.exists():
        print(f"\nError: {DAEN_PATH} not found.")
        print("Run src/ingest/parse_daen.py first.")
        sys.exit(1)

    if not PRESPEC_PATH.exists():
        print(f"\nError: {PRESPEC_PATH} not found.")
        sys.exit(1)

    daen = pd.read_parquet(DAEN_PATH)
    print(f"\nLoaded DAEN: {len(daen):,} records")

    products = load_prespec()
    print(f"Pre-specified products: {len(products)}")

    print(f"\nMatching records:\n")
    matched = match_records(daen, products)

    if matched.empty:
        print("\nWARNING: No matches found. Check search terms and DAEN column names.")
        sys.exit(1)

    total_matched = len(matched)
    unique_cases = matched["case_number"].nunique() if "case_number" in matched.columns else total_matched

    # QC report
    report_lines = [
        "PRODUCT DICTIONARY MATCH REPORT",
        "=" * 50,
        f"Total DAEN records: {len(daen):,}",
        f"Total matched records (study drugs): {total_matched:,}",
        f"Unique matched cases: {unique_cases:,}",
        "",
        "Per-product breakdown:",
    ]

    for prod in products:
        prod_matches = matched[matched["product_id"] == prod["product_id"]]
        report_lines.append(
            f"  {prod['product_set']}/{prod['product_id']}  "
            f"{prod['canonical_name']:35s}  "
            f"{len(prod_matches):>6,} records"
        )

    # Per-set summary
    for s in ("A", "B"):
        set_matches = matched[matched["product_set"] == s]
        set_name = "Antivenoms" if s == "A" else "Tropical medicines"
        report_lines.append(f"\n  Set {s} ({set_name}): {len(set_matches):,} total records")

    report_lines.extend([
        "",
        f"Study drug records as % of total DAEN: {100 * total_matched / len(daen):.2f}%",
        "",
        "QC GATE: Match rate assessment",
        "  (Match rate is assessed post-data-cleaning, not here.)",
        "  (This step identifies raw matches; cleaning may refine counts.)",
    ])

    report_text = "\n".join(report_lines)
    print(f"\n{report_text}")

    # Save
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    matched.to_parquet(PROCESSED_DIR / "daen_study_drugs.parquet", index=False)
    print(f"\nSaved: {PROCESSED_DIR / 'daen_study_drugs.parquet'}")

    # Save final product dictionary with match counts
    dict_rows = []
    for prod in products:
        prod_matches = matched[matched["product_id"] == prod["product_id"]]
        dict_rows.append({
            "product_set": prod["product_set"],
            "product_id": prod["product_id"],
            "canonical_name": prod["canonical_name"],
            "active_ingredient": prod["active_ingredient"],
            "atc_code": prod["atc_code"],
            "search_terms": ";".join(prod["search_terms"]),
            "n_matched_records": len(prod_matches),
            "indication_expected_pts": prod["indication_expected_pts"],
        })
    pd.DataFrame(dict_rows).to_csv(PROCESSED_DIR / "product_dictionary.csv", index=False)
    print(f"Saved: {PROCESSED_DIR / 'product_dictionary.csv'}")

    report_path = PROCESSED_DIR / "product_match_report.txt"
    report_path.write_text(report_text)
    print(f"Saved: {report_path}")

    print(f"\n{'=' * 70}")
    print(f"  Next: python scripts/04_build_reaction_dict.py")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
