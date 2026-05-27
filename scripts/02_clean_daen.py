"""
02 — DAEN Data Cleaning: Parse multi-value fields into normalised tables

Reads the merged DAEN parquet and:
  1. Parses the medicines field → one row per case × drug (with role, brand, ingredient)
  2. Parses the reactions field → one row per case × reaction PT
  3. Builds the drug–reaction pair table for disproportionality analysis
  4. Matches study drugs via the product dictionary
  5. Reports descriptive statistics and QC gate metrics

Outputs (in data/processed/):
  - daen_cases.parquet              One row per case
  - daen_case_drugs.parquet         One row per case × drug
  - daen_case_reactions.parquet     One row per case × reaction
  - daen_drug_reaction_pairs.parquet  All (drug, reaction) pairs for 2×2 tables
  - cleaning_report.txt            QC gate report

Usage:
    python scripts/02_clean_daen.py
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from tqdm import tqdm

PROJECT_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
MERGED_PATH = PROCESSED_DIR / "daen_merged.parquet"
PRODUCT_DICT_PATH = PROCESSED_DIR / "product_dictionary_prespec.csv"


def parse_medicines_field(text: str) -> list[dict]:
    """Parse a DAEN medicines field into structured drug entries.

    Format: • BRAND (active_ingredient) - Role
    """
    if pd.isna(text) or not str(text).strip():
        return []

    entries = []
    for line in str(text).split("\n"):
        line = line.strip().lstrip("•").strip()
        if not line:
            continue

        role = "Unknown"
        for r in ("Suspected", "Not suspected", "Interaction"):
            if line.endswith(f"- {r}"):
                role = r
                line = line[: -(len(r) + 2)].strip()
                break

        brand = line
        ingredient = ""
        paren_match = re.search(r"\(([^)]+)\)\s*$", line)
        if paren_match:
            ingredient = paren_match.group(1).strip()
            brand = line[: paren_match.start()].strip()

        entries.append({
            "brand_name": brand,
            "active_ingredient": ingredient.lower() if ingredient else "",
            "drug_role": role,
            "raw_text": line,
        })

    return entries


def parse_reactions_field(text: str) -> list[str]:
    """Parse a DAEN reactions field into individual PT strings."""
    if pd.isna(text) or not str(text).strip():
        return []

    terms = []
    for line in str(text).split("\n"):
        line = line.strip().lstrip("•").strip()
        if line:
            terms.append(line)
    return terms


def build_case_drugs(df: pd.DataFrame) -> pd.DataFrame:
    """Expand medicines field into one row per case × drug."""
    rows = []
    for case_number, medicines in tqdm(
        zip(df["case_number"], df["medicines"]),
        total=len(df),
        desc="Parsing medicines",
    ):
        for entry in parse_medicines_field(medicines):
            entry["case_number"] = case_number
            rows.append(entry)

    return pd.DataFrame(rows)


def build_case_reactions(df: pd.DataFrame) -> pd.DataFrame:
    """Expand reactions field into one row per case × reaction."""
    rows = []
    for case_number, reactions in tqdm(
        zip(df["case_number"], df["reactions"]),
        total=len(df),
        desc="Parsing reactions",
    ):
        for pt in parse_reactions_field(reactions):
            rows.append({"case_number": case_number, "reaction_pt": pt})

    return pd.DataFrame(rows)


def match_study_drugs(case_drugs: pd.DataFrame) -> pd.DataFrame:
    """Tag case_drugs rows that match study products.

    Polyvalent antivenom is matched first: if the active_ingredient field
    contains >=3 individual antivenom names (semicolon-delimited constituents),
    it is classified as polyvalent per pre-registration Section 9.1.1.
    """
    prespec = pd.read_csv(PRODUCT_DICT_PATH)

    case_drugs["product_id"] = None
    case_drugs["product_set"] = None
    case_drugs["canonical_name"] = None

    search_col = (
        case_drugs["brand_name"].str.lower()
        + " "
        + case_drugs["active_ingredient"]
    )

    monovalent_av_names = [
        "brown snake", "tiger snake", "black snake",
        "taipan", "death adder",
    ]

    # Step 1: Match polyvalent first — ingredient field lists multiple antivenoms
    poly_count = search_col.apply(
        lambda x: sum(1 for av in monovalent_av_names if av in str(x))
    )
    poly_mask = poly_count >= 3
    polyvalent_row = prespec[prespec["product_id"] == "A06"]
    if not polyvalent_row.empty:
        case_drugs.loc[poly_mask, "product_id"] = "A06"
        case_drugs.loc[poly_mask, "product_set"] = "A"
        case_drugs.loc[poly_mask, "canonical_name"] = "Polyvalent snake antivenom"
        print(f"  Polyvalent (multi-constituent match): {poly_mask.sum()} rows")

    # Also match explicit "polyvalent" keyword for cases not caught above
    poly_kw_mask = search_col.str.contains("polyvalent", na=False) & case_drugs["product_id"].isna()
    if poly_kw_mask.any():
        case_drugs.loc[poly_kw_mask, "product_id"] = "A06"
        case_drugs.loc[poly_kw_mask, "product_set"] = "A"
        case_drugs.loc[poly_kw_mask, "canonical_name"] = "Polyvalent snake antivenom"
        print(f"  Polyvalent (keyword match): {poly_kw_mask.sum()} rows")

    # Step 2: Match all other products (skip polyvalent in the loop)
    for _, prod in prespec.iterrows():
        if prod["product_id"] == "A06":
            continue

        search_terms = [
            t.strip().lower()
            for t in str(prod["search_terms"]).split(";")
            if t.strip()
        ]
        if not search_terms:
            continue

        pattern = "|".join(re.escape(t) for t in search_terms)
        mask = search_col.str.contains(pattern, regex=True, na=False)
        mask = mask & case_drugs["product_id"].isna()

        case_drugs.loc[mask, "product_id"] = prod["product_id"]
        case_drugs.loc[mask, "product_set"] = prod["product_set"]
        case_drugs.loc[mask, "canonical_name"] = prod["canonical_name"]

    return case_drugs


def main():
    print("=" * 70)
    print("  DAEN Data Cleaning — DAEN Tropical PV")
    print("=" * 70)

    if not MERGED_PATH.exists():
        print(f"Error: {MERGED_PATH} not found. Run parse_daen.py first.")
        sys.exit(1)

    df = pd.read_parquet(MERGED_PATH)
    print(f"\nLoaded: {len(df):,} records")

    # Build cases table
    cases = df[["case_number", "report_date", "age", "sex"]].copy()
    cases["age_numeric"] = pd.to_numeric(cases["age"], errors="coerce")
    cases["age_group"] = pd.cut(
        cases["age_numeric"],
        bins=[0, 18, 65, 200],
        labels=["<18", "18-64", ">=65"],
        right=False,
    )
    print(f"Cases: {len(cases):,}")

    # Parse medicines → case_drugs
    print()
    case_drugs = build_case_drugs(df)
    print(f"\nCase-drug rows: {len(case_drugs):,}")
    print(f"Unique drugs: {case_drugs['brand_name'].nunique():,}")

    print(f"\nDrug role distribution:")
    for role, count in case_drugs["drug_role"].value_counts().items():
        pct = 100 * count / len(case_drugs)
        print(f"  {role:20s}  {count:>10,} ({pct:5.1f}%)")

    # Parse reactions → case_reactions
    print()
    case_reactions = build_case_reactions(df)
    print(f"\nCase-reaction rows: {len(case_reactions):,}")
    print(f"Unique reaction PTs: {case_reactions['reaction_pt'].nunique():,}")

    # Match study drugs
    print("\nMatching study drugs...")
    case_drugs = match_study_drugs(case_drugs)

    study_drug_rows = case_drugs[case_drugs["product_id"].notna()]
    print(f"Study drug case-drug rows: {len(study_drug_rows):,}")
    print(f"\nPer-product matches:")
    for name, group in study_drug_rows.groupby("canonical_name"):
        n_suspected = (group["drug_role"] == "Suspected").sum()
        n_total = len(group)
        print(f"  {name:35s}  {n_total:>5,} total  ({n_suspected:>5,} Suspected)")

    # Build drug-reaction pairs for 2×2 tables
    print("\nBuilding drug-reaction pairs...")
    dr_pairs = case_drugs.merge(case_reactions, on="case_number")
    print(f"Total drug-reaction pairs: {len(dr_pairs):,}")

    study_dr_pairs = dr_pairs[dr_pairs["product_id"].notna()]
    print(f"Study drug-reaction pairs: {len(study_dr_pairs):,}")

    # QC report
    report_lines = [
        "DAEN CLEANING REPORT",
        "=" * 50,
        f"Total cases: {len(cases):,}",
        f"Total case-drug rows: {len(case_drugs):,}",
        f"Total case-reaction rows: {len(case_reactions):,}",
        f"Total drug-reaction pairs: {len(dr_pairs):,}",
        "",
        f"Study drug case-drug rows: {len(study_drug_rows):,}",
        f"Study drug-reaction pairs: {len(study_dr_pairs):,}",
        "",
        "Drug role distribution:",
    ]
    for role, count in case_drugs["drug_role"].value_counts().items():
        report_lines.append(f"  {role:20s}  {count:>10,}")

    report_lines.extend([
        "",
        "Age completeness:",
        f"  Numeric age available: {cases['age_numeric'].notna().sum():,} / {len(cases):,} "
        f"({100 * cases['age_numeric'].notna().mean():.1f}%)",
        "",
        "Sex distribution:",
    ])
    for sex, count in cases["sex"].value_counts().items():
        report_lines.append(f"  {sex:20s}  {count:>10,}")

    report_lines.extend([
        "",
        "Top 10 reaction PTs (all drugs):",
    ])
    for pt, count in case_reactions["reaction_pt"].value_counts().head(10).items():
        report_lines.append(f"  {pt:40s}  {count:>8,}")

    report_text = "\n".join(report_lines)
    print(f"\n{report_text}")

    # Save
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    cases.to_parquet(PROCESSED_DIR / "daen_cases.parquet", index=False)
    case_drugs.to_parquet(PROCESSED_DIR / "daen_case_drugs.parquet", index=False)
    case_reactions.to_parquet(PROCESSED_DIR / "daen_case_reactions.parquet", index=False)
    dr_pairs.to_parquet(PROCESSED_DIR / "daen_drug_reaction_pairs.parquet", index=False)

    (PROCESSED_DIR / "cleaning_report.txt").write_text(report_text)

    for name in ("daen_cases", "daen_case_drugs", "daen_case_reactions", "daen_drug_reaction_pairs"):
        p = PROCESSED_DIR / f"{name}.parquet"
        mb = p.stat().st_size / (1024 * 1024)
        print(f"  Saved: {p.name} ({mb:.1f} MB)")

    print(f"\n{'=' * 70}")
    print(f"  Next: python scripts/04_build_reaction_dict.py")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
