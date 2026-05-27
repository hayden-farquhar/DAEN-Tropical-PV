"""
03 — FAERS Comparator Data: Fetch tropical-basket drug AE reports from openFDA

Queries the openFDA FAERS API for each Set B (tropical-medicine) drug to
build the FAERS comparator dataset for triangulation (pre-registration 9.12).

For Set A (antivenoms), FAERS is expected to have near-zero reports —
this is documented as comparator-absence per the pre-registration.

Outputs:
  - data/raw/faers_tropical_basket.parquet    FAERS reports for Set B drugs
  - data/raw/faers_antivenom_check.json       Confirms near-zero AV reports

Usage:
    python scripts/05_parse_faers_comparator.py
"""

import json
import sys
import time
from pathlib import Path

import pandas as pd
import requests

PROJECT_DIR = Path(__file__).resolve().parent.parent
RAW_DIR = PROJECT_DIR / "data" / "raw"

FAERS_BASE = "https://api.fda.gov/drug/event.json"

TROPICAL_QUERIES = {
    "doxycycline": 'patient.drug.openfda.generic_name:"doxycycline"',
    "primaquine": 'patient.drug.openfda.generic_name:"primaquine"',
    "tafenoquine": 'patient.drug.openfda.generic_name:"tafenoquine"',
    "artemether_lumefantrine": (
        'patient.drug.openfda.generic_name:"artemether"+'
        'patient.drug.openfda.generic_name:"lumefantrine"'
    ),
    "ivermectin": 'patient.drug.openfda.generic_name:"ivermectin"',
    "praziquantel": 'patient.drug.openfda.generic_name:"praziquantel"',
}

ANTIVENOM_QUERIES = {
    "antivenom": 'patient.drug.openfda.generic_name:"antivenom"',
    "antivenin": 'patient.drug.openfda.generic_name:"antivenin"',
    "antivenene": 'patient.drug.openfda.generic_name:"antivenene"',
}


def query_faers_count(search: str) -> int:
    params = {"search": search, "limit": 1}
    try:
        r = requests.get(FAERS_BASE, params=params, timeout=30)
        if r.status_code == 200:
            return r.json().get("meta", {}).get("results", {}).get("total", 0)
        if r.status_code == 404:
            return 0
    except requests.RequestException as e:
        print(f"    Request error: {e}")
    return -1


def query_faers_reactions(search: str, limit: int = 100) -> list[dict]:
    params = {
        "search": search,
        "count": "patient.reaction.reactionmeddrapt.exact",
        "limit": limit,
    }
    try:
        r = requests.get(FAERS_BASE, params=params, timeout=30)
        if r.status_code == 200:
            return r.json().get("results", [])
    except requests.RequestException:
        pass
    return []


def main():
    print("=" * 70)
    print("  FAERS Comparator Data — DAEN Tropical PV")
    print("=" * 70)

    RAW_DIR.mkdir(parents=True, exist_ok=True)

    # Check antivenom comparator-absence
    print("\nAntivenom comparator-absence check:")
    av_results = {}
    for name, query in ANTIVENOM_QUERIES.items():
        count = query_faers_count(query)
        av_results[name] = count
        print(f"  {name:20s}  {count:>8,} reports")
        time.sleep(0.5)

    av_path = RAW_DIR / "faers_antivenom_check.json"
    with open(av_path, "w") as f:
        json.dump(av_results, f, indent=2)
    print(f"\n  Saved: {av_path}")

    # Fetch tropical basket reaction profiles
    print("\nTropical-basket reaction profiles:")
    all_rows = []
    for drug_name, query in TROPICAL_QUERIES.items():
        count = query_faers_count(query)
        print(f"\n  {drug_name:25s}  {count:>8,} total FAERS reports")

        reactions = query_faers_reactions(query, limit=100)
        for rx in reactions:
            all_rows.append({
                "drug": drug_name,
                "reaction_pt": rx["term"],
                "faers_count": rx["count"],
                "faers_total_for_drug": count,
            })

        if reactions:
            top3 = reactions[:3]
            for rx in top3:
                print(f"    {rx['term']:40s}  {rx['count']:>8,}")

        time.sleep(1)

    if all_rows:
        df = pd.DataFrame(all_rows)
        out_path = RAW_DIR / "faers_tropical_basket.parquet"
        df.to_parquet(out_path, index=False)
        print(f"\n  Saved: {out_path} ({len(df)} drug-reaction pairs)")
    else:
        print("\n  WARNING: No FAERS reaction data retrieved.")

    print(f"\n{'=' * 70}")
    print(f"  FAERS comparator data complete.")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
