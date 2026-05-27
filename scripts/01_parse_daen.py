"""
01 — DAEN Data Acquisition: Merge Batched Exports

Scans data/raw/ for all Excel (.xlsx) and CSV (.csv) files downloaded from
the TGA DAEN web interface. Merges them into a single dataset, removes
duplicate cases, and saves to data/processed/daen_merged.parquet.

Adapted from the P14 TGA DAEN pipeline, with modifications for the
Antivenom/tropical-medicine pre-registration QC gates.

Usage:
    python scripts/01_parse_daen.py
"""

import hashlib
import json
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = PROJECT_DIR / "data" / "raw"
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"

COLUMN_MAP = {
    "case number": "case_number",
    "case no": "case_number",
    "case no.": "case_number",
    "case_no": "case_number",
    "case id": "case_number",
    "report entry date": "report_date",
    "report date": "report_date",
    "date of report": "report_date",
    "date": "report_date",
    "age": "age",
    "age (years)": "age",
    "age group": "age",
    "patient age": "age",
    "sex": "sex",
    "gender": "sex",
    "patient sex": "sex",
    "patient gender": "sex",
    "medicines": "medicines",
    "medicine": "medicines",
    "medicines reported as being taken": "medicines",
    "drug": "medicines",
    "drugs": "medicines",
    "drug name": "medicines",
    "product name": "medicines",
    "active ingredient": "active_ingredient",
    "active ingredients": "active_ingredient",
    "ingredient": "active_ingredient",
    "reactions": "reactions",
    "reaction": "reactions",
    "reaction term": "reactions",
    "meddra reaction terms": "reactions",
    "meddra reaction term": "reactions",
    "adverse event": "reactions",
    "adverse events": "reactions",
    "adverse reaction": "reactions",
    "preferred term": "reactions",
    "meddra preferred term": "reactions",
    "reporter type": "reporter_type",
    "reporter": "reporter_type",
    "source": "reporter_type",
    "report source": "reporter_type",
    "outcome": "outcome",
    "seriousness": "seriousness",
}


def print_download_guide():
    print("""
╔══════════════════════════════════════════════════════════════════════╗
║         DAEN DOWNLOAD GUIDE — DAEN Tropical PV             ║
╠══════════════════════════════════════════════════════════════════════╣
║                                                                    ║
║  FULL DATABASE EXPORT REQUIRED for 2×2 table denominators.         ║
║  Download ALL reports (blank product name), not just study drugs.   ║
║                                                                    ║
╚══════════════════════════════════════════════════════════════════════╝

STEPS:

  1. Go to: https://daen.tga.gov.au/medicines-search/
     (redirects to https://aems.tga.gov.au)
  2. Accept the terms and conditions
  3. Leave the product name field BLANK (to get ALL reports)
  4. Set a date range (see batches below)
  5. Click "Search"
  6. Click the "List of Reports" tab
  7. Click the three-dot menu icon (⋮) → "Data with current layout"
  8. Save the Excel file to: data/raw/
  9. Repeat for each date range batch

SUGGESTED DATE RANGE BATCHES:

  Each batch must stay under 150,000 rows. Use these ranges as a
  starting point — if a batch exceeds the limit, split it further.

  Batch 01:  01/01/1971  →  31/12/2005   (early sparse years)
  Batch 02:  01/01/2006  →  31/12/2010
  Batch 03:  01/01/2011  →  31/12/2013
  Batch 04:  01/01/2014  →  31/12/2016
  Batch 05:  01/01/2017  →  31/12/2018
  Batch 06:  01/01/2019  →  31/12/2019
  Batch 07:  01/01/2020  →  31/12/2020
  Batch 08:  01/01/2021  →  31/12/2021
  Batch 09:  01/01/2022  →  31/12/2022
  Batch 10:  01/01/2023  →  31/12/2023
  Batch 11:  01/01/2024  →  31/12/2024
  Batch 12:  01/01/2025  →  today

  Date format on the DAEN site: DD/MM/YYYY.

NAMING CONVENTION:

  daen_batch_01_1971_2005.xlsx
  daen_batch_02_2006_2010.xlsx
  ... etc.

After downloading, re-run this script to merge all batches.
""")


def find_raw_files():
    if not RAW_DIR.exists():
        return []
    files = []
    for ext in ("*.xlsx", "*.xls", "*.csv"):
        files.extend(RAW_DIR.glob(ext))
    return sorted(f for f in files if not f.name.startswith(("~", ".")))


def file_sha256(filepath: Path) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def read_file(filepath: Path) -> pd.DataFrame | None:
    suffix = filepath.suffix.lower()
    try:
        if suffix == ".csv":
            return pd.read_csv(filepath, low_memory=False)
        if suffix in (".xlsx", ".xls"):
            engine = "openpyxl" if suffix == ".xlsx" else None
            return pd.read_excel(filepath, engine=engine)
    except Exception as e:
        print(f"  Error reading {filepath.name}: {e}")
    return None


def standardise_columns(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = df.columns.str.strip().str.lower()
    rename = {c: COLUMN_MAP[c] for c in df.columns if c in COLUMN_MAP}
    return df.rename(columns=rename)


def main():
    print("=" * 70)
    print("  DAEN Data Acquisition — DAEN Tropical PV")
    print("=" * 70)

    files = find_raw_files()
    if not files:
        print_download_guide()
        sys.exit(1)

    print(f"\nFound {len(files)} file(s) in {RAW_DIR}/\n")

    snapshot_log = {
        "export_date": datetime.now().strftime("%Y-%m-%d"),
        "files": [],
    }

    dfs = []
    for f in files:
        size_mb = f.stat().st_size / (1024 * 1024)
        print(f"  Reading {f.name} ({size_mb:.1f} MB) ...")

        snapshot_log["files"].append({
            "name": f.name,
            "size_bytes": f.stat().st_size,
            "sha256": file_sha256(f),
        })

        df = read_file(f)
        if df is None:
            continue

        df = standardise_columns(df)
        print(f"    Rows: {len(df):,}  |  Columns: {list(df.columns)}")
        df["_source_file"] = f.name
        dfs.append(df)

    if not dfs:
        print("\nNo valid data read. Check files and try again.")
        sys.exit(1)

    merged = pd.concat(dfs, ignore_index=True)
    print(f"\nConcatenated total: {len(merged):,} rows")

    # Deduplicate on case_number (within-channel dedup per pre-registration 9.1.2)
    if "case_number" in merged.columns:
        before = len(merged)
        merged = merged.drop_duplicates(subset="case_number", keep="last")
        dupes = before - len(merged)
        print(f"Deduplicated on case_number: removed {dupes:,} duplicates (kept most recent)")
    else:
        before = len(merged)
        dedup_cols = [c for c in merged.columns if c != "_source_file"]
        merged = merged.drop_duplicates(subset=dedup_cols, keep="first")
        dupes = before - len(merged)
        print(f"Warning: No 'case_number' column. Full-row dedup removed {dupes:,}")

    # Drop metadata/filter rows that DAEN exports append
    if "case_number" in merged.columns:
        non_numeric = pd.to_numeric(merged["case_number"], errors="coerce").isna()
        n_meta = non_numeric.sum()
        if n_meta > 0:
            sample = merged.loc[non_numeric, "case_number"].head(3).tolist()
            print(f"Dropped {n_meta:,} metadata/filter rows (e.g., {sample[0][:60]}...)")
            merged = merged[~non_numeric].copy()
            merged["case_number"] = merged["case_number"].astype(int)

    # Parse dates
    if "report_date" in merged.columns:
        merged["report_date"] = pd.to_datetime(
            merged["report_date"], errors="coerce", dayfirst=True
        )
        valid_dates = merged["report_date"].notna().sum()
        if valid_dates > 0:
            print(f"Date range: {merged['report_date'].min().date()} → "
                  f"{merged['report_date'].max().date()}")
            print(f"Valid dates: {valid_dates:,} / {len(merged):,}")

    # Summary
    print(f"\n{'=' * 70}")
    print(f"  MERGED DATASET SUMMARY")
    print(f"{'=' * 70}")
    print(f"\n  Total unique records: {len(merged):,}")

    print(f"\n  Column completeness:")
    for col in merged.columns:
        if col == "_source_file":
            continue
        non_null = merged[col].notna().sum()
        pct = 100 * non_null / len(merged)
        print(f"    {col:25s}  {non_null:>10,} non-null ({pct:5.1f}%)")

    if "reporter_type" in merged.columns:
        print(f"\n  Reporter type distribution:")
        for rt, count in merged["reporter_type"].value_counts().items():
            pct = 100 * count / len(merged)
            print(f"    {rt:25s}  {count:>10,} ({pct:5.1f}%)")

    # Save
    PROCESSED_DIR.mkdir(parents=True, exist_ok=True)

    parquet_path = PROCESSED_DIR / "daen_merged.parquet"
    merged.to_parquet(parquet_path, index=False)
    pq_mb = parquet_path.stat().st_size / (1024 * 1024)
    print(f"\n  Saved: {parquet_path} ({pq_mb:.1f} MB)")

    snapshot_log["total_records"] = len(merged)
    snapshot_log["total_records_before_dedup"] = before
    snapshot_log["duplicates_removed"] = dupes
    snapshot_log["merged_sha256"] = file_sha256(parquet_path)

    log_path = PROCESSED_DIR / "data_snapshot_log.json"
    with open(log_path, "w") as f:
        json.dump(snapshot_log, f, indent=2, default=str)
    print(f"  Snapshot log: {log_path}")

    print(f"\n{'=' * 70}")
    print(f"  Next: python scripts/03_build_product_dict.py")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
