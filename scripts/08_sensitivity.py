"""
06 — Sensitivity Analyses S1–S8 (Phase 6, Section 9.9)

Pre-registered sensitivity analyses:
  S1: Serious reports only          → INFEASIBLE (no seriousness field in DAEN export)
  S2: HCP reports only              → INFEASIBLE (no reporter-type field; see Amendment 1)
  S3: Raised minimum n ≥ 10         → Feasible
  S4: Raised IC025 > 0.5 threshold  → Feasible
  S5: Drop most recent calendar year → Feasible
  S6: Restrict to 2010–present      → Feasible
  S7: Remove anaphylaxis PTs (Set A) → Feasible
  S8: Remove indication-expected PTs (Set A) → Feasible

Outputs (in outputs/tables/):
  - sensitivity_S3_n10.csv
  - sensitivity_S4_ic05.csv
  - sensitivity_S5_drop_recent.csv
  - sensitivity_S6_post2010.csv
  - sensitivity_S7_no_anaphylaxis.csv
  - sensitivity_S8_no_indication.csv
  - sensitivity_comparison.csv      Signal-level comparison across all analyses
  - sensitivity_summary.txt

Usage:
    python scripts/08_sensitivity.py
"""

from pathlib import Path

import numpy as np
import pandas as pd

import sys as _sys
_sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module as _import
_disp = _import("06_disproportionality")
run_disproportionality = _disp.run_disproportionality
MIN_N = _disp.MIN_N
PRR_THRESHOLD = _disp.PRR_THRESHOLD
CHI2_THRESHOLD = _disp.CHI2_THRESHOLD
IC025_THRESHOLD = _disp.IC025_THRESHOLD
EB05_THRESHOLD = _disp.EB05_THRESHOLD

PROJECT_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
RESULTS_DIR = PROJECT_DIR / "outputs" / "tables"

ANAPHYLAXIS_PTS = {
    "Anaphylactic reaction",
    "Anaphylactic shock",
    "Anaphylactoid reaction",
    "Anaphylactoid shock",
    "Hypersensitivity",
    "Type I hypersensitivity",
    "Drug hypersensitivity",
    "Allergic reaction",
}

INDICATION_EXPECTED_PTS = {
    # Envenomation syndrome PTs (expected given the indication for antivenoms)
    "Venom poisoning",
    "Snake bite",
    "Snakebite",
    "Envenomation",
    "Coagulopathy",
    "Venom-induced consumptive coagulopathy",
    "Neurotoxicity",
    "Myotoxicity",
    "Rhabdomyolysis",
    "Paralysis",
    "Pain in extremity",
    "Swelling",
    "Oedema peripheral",
    "Limb discomfort",
    "Bite site reaction",
    "Sting",
    "Jellyfish sting",
}


def load_data():
    cases = pd.read_parquet(PROCESSED_DIR / "daen_cases.parquet")
    case_drugs = pd.read_parquet(PROCESSED_DIR / "daen_case_drugs.parquet")
    case_reactions = pd.read_parquet(PROCESSED_DIR / "daen_case_reactions.parquet")
    product_dict = pd.read_csv(PROCESSED_DIR / "product_dictionary_prespec.csv")
    return cases, case_drugs, case_reactions, product_dict


def run_sensitivity(
    case_drugs: pd.DataFrame,
    case_reactions: pd.DataFrame,
    product_dict: pd.DataFrame,
    label: str,
) -> pd.DataFrame:
    """Wrapper: run disproportionality and return signals."""
    result = run_disproportionality(
        case_drugs, case_reactions, product_dict, label=label
    )
    return result


def s3_raised_n(primary_df: pd.DataFrame) -> pd.DataFrame:
    """S3: Re-threshold the primary results at n >= 10."""
    df = primary_df.copy()
    df["signal_s3"] = (
        (df["n"] >= 10)
        & (df["prr"] >= PRR_THRESHOLD)
        & (df["chi2"] >= CHI2_THRESHOLD)
        & (df["ic025"] > IC025_THRESHOLD)
        & (df["eb05"] >= EB05_THRESHOLD)
    )
    return df


def s4_raised_ic(primary_df: pd.DataFrame) -> pd.DataFrame:
    """S4: Re-threshold the primary results at IC025 > 0.5."""
    df = primary_df.copy()
    df["signal_s4"] = (
        (df["n"] >= MIN_N)
        & (df["prr"] >= PRR_THRESHOLD)
        & (df["chi2"] >= CHI2_THRESHOLD)
        & (df["ic025"] > 0.5)
        & (df["eb05"] >= EB05_THRESHOLD)
    )
    return df


def s5_drop_recent(
    cases: pd.DataFrame,
    case_drugs: pd.DataFrame,
    case_reactions: pd.DataFrame,
    product_dict: pd.DataFrame,
) -> pd.DataFrame:
    """S5: Drop the most recent calendar year."""
    max_year = cases["report_date"].dt.year.max()
    keep_cases = set(cases[cases["report_date"].dt.year < max_year]["case_number"])
    cd = case_drugs[case_drugs["case_number"].isin(keep_cases)]
    cr = case_reactions[case_reactions["case_number"].isin(keep_cases)]
    return run_sensitivity(cd, cr, product_dict, label=f"S5_drop_{max_year}")


def s6_post2010(
    cases: pd.DataFrame,
    case_drugs: pd.DataFrame,
    case_reactions: pd.DataFrame,
    product_dict: pd.DataFrame,
) -> pd.DataFrame:
    """S6: Restrict to 2010-present."""
    keep_cases = set(cases[cases["report_date"].dt.year >= 2010]["case_number"])
    cd = case_drugs[case_drugs["case_number"].isin(keep_cases)]
    cr = case_reactions[case_reactions["case_number"].isin(keep_cases)]
    return run_sensitivity(cd, cr, product_dict, label="S6_post2010")


def s7_no_anaphylaxis(
    case_drugs: pd.DataFrame,
    case_reactions: pd.DataFrame,
    product_dict: pd.DataFrame,
) -> pd.DataFrame:
    """S7: Remove anaphylaxis/hypersensitivity PTs for Set A only."""
    set_a = product_dict[product_dict["product_set"] == "A"]
    cr = case_reactions[~case_reactions["reaction_pt"].isin(ANAPHYLAXIS_PTS)]
    return run_sensitivity(case_drugs, cr, set_a, label="S7_no_anaphylaxis")


def s8_no_indication(
    case_drugs: pd.DataFrame,
    case_reactions: pd.DataFrame,
    product_dict: pd.DataFrame,
) -> pd.DataFrame:
    """S8: Remove indication-expected PTs for Set A only."""
    set_a = product_dict[product_dict["product_set"] == "A"]
    cr = case_reactions[~case_reactions["reaction_pt"].isin(INDICATION_EXPECTED_PTS)]
    return run_sensitivity(case_drugs, cr, set_a, label="S8_no_indication")


def compare_signals(primary: pd.DataFrame, analyses: dict) -> pd.DataFrame:
    """Build a signal-level comparison across all sensitivity analyses.

    For each primary signal, check whether it remains a signal in each sensitivity.
    """
    primary_signals = primary[primary["signal"]].copy()
    if primary_signals.empty:
        return pd.DataFrame()

    comparison = primary_signals[
        ["product_id", "product_set", "canonical_name", "reaction_pt", "n", "prr", "ic025"]
    ].copy()
    comparison = comparison.rename(columns={"n": "n_primary", "prr": "prr_primary", "ic025": "ic025_primary"})
    comparison["signal_primary"] = True

    for sa_label, sa_df in analyses.items():
        if sa_df.empty:
            comparison[f"signal_{sa_label}"] = np.nan
            comparison[f"n_{sa_label}"] = np.nan
            continue

        if "signal" in sa_df.columns:
            signal_col = "signal"
        elif f"signal_{sa_label.lower()}" in sa_df.columns:
            signal_col = f"signal_{sa_label.lower()}"
        else:
            comparison[f"signal_{sa_label}"] = np.nan
            comparison[f"n_{sa_label}"] = np.nan
            continue

        sa_lookup = sa_df.set_index(["product_id", "reaction_pt"])
        for idx, row in comparison.iterrows():
            key = (row["product_id"], row["reaction_pt"])
            if key in sa_lookup.index:
                sa_row = sa_lookup.loc[key]
                if isinstance(sa_row, pd.DataFrame):
                    sa_row = sa_row.iloc[0]
                comparison.at[idx, f"signal_{sa_label}"] = bool(sa_row[signal_col])
                comparison.at[idx, f"n_{sa_label}"] = sa_row.get("n", np.nan)
            else:
                comparison.at[idx, f"signal_{sa_label}"] = False
                comparison.at[idx, f"n_{sa_label}"] = 0

    return comparison


def main():
    print("=" * 70)
    print("  Sensitivity Analyses (Phase 6) — DAEN Tropical PV")
    print("=" * 70)

    cases, case_drugs, case_reactions, product_dict = load_data()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # Load primary results for comparison
    primary = pd.read_csv(RESULTS_DIR / "signals_full_scan.csv")
    primary_signals = primary[primary["signal"]].copy()
    n_primary = primary_signals.shape[0]
    print(f"\nPrimary signals: {n_primary}")

    analyses = {}

    # ── S1 & S2: Infeasible ───────────────────────────────────────────
    print("\n--- S1: Serious reports only ---")
    print("  INFEASIBLE — DAEN public export lacks seriousness/outcome field.")
    print("  Documented as pre-registered but data-infeasible.")

    print("\n--- S2: HCP reports only ---")
    print("  INFEASIBLE — DAEN public export lacks reporter-type field.")
    print("  See Amendment 1 (drug-role substitution for H2).")

    # ── S3: Raised minimum n ≥ 10 ─────────────────────────────────────
    print("\n--- S3: Raised minimum n ≥ 10 ---")
    s3 = s3_raised_n(primary)
    n_s3 = s3["signal_s3"].sum()
    print(f"  Signals at n ≥ 10: {n_s3} (vs {n_primary} at n ≥ 3)")
    s3_signals = s3[s3["signal_s3"]].copy()
    s3_signals["signal"] = s3_signals["signal_s3"]
    s3_signals.to_csv(RESULTS_DIR / "sensitivity_S3_n10.csv", index=False)
    analyses["S3"] = s3_signals

    # ── S4: Raised IC025 > 0.5 ────────────────────────────────────────
    print("\n--- S4: Raised IC025 > 0.5 ---")
    s4 = s4_raised_ic(primary)
    n_s4 = s4["signal_s4"].sum()
    print(f"  Signals at IC025 > 0.5: {n_s4} (vs {n_primary} at IC025 > 0)")
    s4_signals = s4[s4["signal_s4"]].copy()
    s4_signals["signal"] = s4_signals["signal_s4"]
    s4_signals.to_csv(RESULTS_DIR / "sensitivity_S4_ic05.csv", index=False)
    analyses["S4"] = s4_signals

    # ── S5: Drop most recent calendar year ─────────────────────────────
    print("\n--- S5: Drop most recent calendar year ---")
    s5 = s5_drop_recent(cases, case_drugs, case_reactions, product_dict)
    n_s5 = s5["signal"].sum()
    max_year = cases["report_date"].dt.year.max()
    print(f"  Dropped year: {max_year}")
    print(f"  Signals: {n_s5} (vs {n_primary} primary)")
    s5[s5["signal"]].to_csv(RESULTS_DIR / "sensitivity_S5_drop_recent.csv", index=False)
    analyses["S5"] = s5

    # ── S6: Restrict to 2010–present ──────────────────────────────────
    print("\n--- S6: Restrict to 2010–present ---")
    s6 = s6_post2010(cases, case_drugs, case_reactions, product_dict)
    n_s6 = s6["signal"].sum()
    print(f"  Signals (2010+): {n_s6} (vs {n_primary} primary)")
    s6[s6["signal"]].to_csv(RESULTS_DIR / "sensitivity_S6_post2010.csv", index=False)
    analyses["S6"] = s6

    # ── S7: Remove anaphylaxis PTs (Set A) ─────────────────────────────
    print("\n--- S7: Remove anaphylaxis/hypersensitivity PTs (Set A) ---")
    s7 = s7_no_anaphylaxis(case_drugs, case_reactions, product_dict)
    n_s7 = s7["signal"].sum()
    set_a_primary = primary_signals[primary_signals["product_set"] == "A"].shape[0]
    print(f"  Set A signals without anaphylaxis PTs: {n_s7} (vs {set_a_primary} Set A primary)")
    removed_pts = ANAPHYLAXIS_PTS & set(primary_signals[primary_signals["product_set"] == "A"]["reaction_pt"])
    if removed_pts:
        print(f"  Anaphylaxis PTs that were signals: {removed_pts}")
    s7[s7["signal"]].to_csv(RESULTS_DIR / "sensitivity_S7_no_anaphylaxis.csv", index=False)
    analyses["S7"] = s7

    # ── S8: Remove indication-expected PTs (Set A) ─────────────────────
    print("\n--- S8: Remove indication-expected PTs (Set A) ---")
    s8 = s8_no_indication(case_drugs, case_reactions, product_dict)
    n_s8 = s8["signal"].sum()
    print(f"  Set A signals without indication PTs: {n_s8} (vs {set_a_primary} Set A primary)")
    removed_ind = INDICATION_EXPECTED_PTS & set(primary_signals[primary_signals["product_set"] == "A"]["reaction_pt"])
    if removed_ind:
        print(f"  Indication PTs that were signals: {removed_ind}")
    s8[s8["signal"]].to_csv(RESULTS_DIR / "sensitivity_S8_no_indication.csv", index=False)
    analyses["S8"] = s8

    # ── Signal comparison across all analyses ──────────────────────────
    print("\n--- Signal comparison across analyses ---")
    comp = compare_signals(primary, analyses)
    if not comp.empty:
        signal_cols = [c for c in comp.columns if c.startswith("signal_") and c != "signal_primary"]
        comp["retained_in"] = comp[signal_cols].apply(
            lambda row: sum(1 for v in row if v is True), axis=1
        )
        comp["retained_pct"] = comp["retained_in"] / len(signal_cols)

        robust = comp[comp["retained_pct"] >= 0.8]
        fragile = comp[comp["retained_pct"] < 0.5]

        print(f"  Total primary signals compared: {len(comp)}")
        print(f"  Robust (retained in ≥80% of analyses): {len(robust)}")
        print(f"  Fragile (retained in <50% of analyses): {len(fragile)}")

        if len(fragile) > 0:
            print("\n  Fragile signals:")
            for _, row in fragile.iterrows():
                print(f"    {row['canonical_name']:35s} × {row['reaction_pt']:35s} "
                      f"retained={row['retained_in']}/{len(signal_cols)}")

        comp.to_csv(RESULTS_DIR / "sensitivity_comparison.csv", index=False)

    # ── Summary ────────────────────────────────────────────────────────
    summary_lines = [
        "SENSITIVITY ANALYSIS SUMMARY (Phase 6, Section 9.9)",
        "=" * 55,
        "",
        f"Primary signals: {n_primary}",
        "",
        "S1  Serious only:          INFEASIBLE (no seriousness field)",
        "S2  HCP only:              INFEASIBLE (no reporter-type field)",
        f"S3  n ≥ 10:                {n_s3} signals ({n_s3/n_primary*100:.0f}% retained)",
        f"S4  IC025 > 0.5:           {n_s4} signals ({n_s4/n_primary*100:.0f}% retained)",
        f"S5  Drop {max_year}:            {n_s5} signals ({n_s5/n_primary*100:.0f}% retained)",
        f"S6  2010+ only:            {n_s6} signals ({n_s6/n_primary*100:.0f}% retained)",
        f"S7  No anaphylaxis (A):    {n_s7} signals (Set A: {n_s7}/{set_a_primary})",
        f"S8  No indication (A):     {n_s8} signals (Set A: {n_s8}/{set_a_primary})",
    ]

    if not comp.empty:
        summary_lines.extend([
            "",
            "Signal robustness:",
            f"  Robust (≥80% retention): {len(robust)}",
            f"  Fragile (<50% retention): {len(fragile)}",
        ])

    summary_text = "\n".join(summary_lines)
    (RESULTS_DIR / "sensitivity_summary.txt").write_text(summary_text)

    print(f"\n{summary_text}")
    print(f"\nSaved to: {RESULTS_DIR}/")
    print(f"\n{'=' * 70}")
    print("  Sensitivity analyses complete (Phase 6).")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
