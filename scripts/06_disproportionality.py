"""
04 — Disproportionality Analysis: PRR, ROR, BCPNN (IC), MGPS (EBGM)

Computes all four disproportionality measures for every (study drug, reaction PT)
pair with n >= 3 reports. Applies the pre-registered conjunction criterion:
  PRR >= 2 AND chi2 >= 4 AND n >= 3 AND IC025 > 0 AND EB05 >= 2

Also runs:
  - Drug-role stratification (H2, Amendment 1): Suspected-only vs all-role
  - Validation controls (positive and negative)

Outputs:
  - outputs/tables/signals_primary.csv         All pairs meeting conjunction criterion
  - outputs/tables/signals_full_scan.csv       All pairs with n >= 3 (for near-miss review)
  - outputs/tables/drug_role_stratification.csv  H2 masking candidates
  - outputs/tables/validation_controls.csv     Positive/negative control results
  - outputs/tables/disproportionality_summary.txt

Usage:
    python scripts/06_disproportionality.py
"""

import sys
from math import log2, sqrt
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

PROJECT_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
RESULTS_DIR = PROJECT_DIR / "outputs" / "tables"

MIN_N = 3
PRR_THRESHOLD = 2.0
CHI2_THRESHOLD = 4.0
IC025_THRESHOLD = 0.0
EB05_THRESHOLD = 2.0


def build_2x2(drug_cases: set, event_cases: set, total_cases: int) -> tuple:
    """Build the 2x2 contingency table."""
    a = len(drug_cases & event_cases)
    b = len(drug_cases - event_cases)
    c = len(event_cases - drug_cases)
    d = total_cases - len(drug_cases | event_cases)
    return a, b, c, d


def compute_prr(a: int, b: int, c: int, d: int) -> dict:
    """Proportional Reporting Ratio with chi-squared."""
    if a == 0 or (a + b) == 0 or (c + d) == 0 or c == 0:
        return {"prr": np.nan, "prr_ci_lo": np.nan, "prr_ci_hi": np.nan, "chi2": np.nan}

    prr = (a / (a + b)) / (c / (c + d))

    # Chi-squared from 2x2 table
    observed = np.array([[a, b], [c, d]])
    if observed.min() >= 0:
        chi2_val = stats.chi2_contingency(observed, correction=False)[0]
    else:
        chi2_val = np.nan

    # 95% CI for log(PRR)
    se_log = sqrt(1/a - 1/(a+b) + 1/c - 1/(c+d)) if a > 0 and c > 0 else np.nan
    if not np.isnan(se_log):
        log_prr = np.log(prr)
        ci_lo = np.exp(log_prr - 1.96 * se_log)
        ci_hi = np.exp(log_prr + 1.96 * se_log)
    else:
        ci_lo = ci_hi = np.nan

    return {"prr": prr, "prr_ci_lo": ci_lo, "prr_ci_hi": ci_hi, "chi2": chi2_val}


def compute_ror(a: int, b: int, c: int, d: int) -> dict:
    """Reporting Odds Ratio with 95% CI."""
    if a == 0 or b == 0 or c == 0 or d == 0:
        return {"ror": np.nan, "ror_ci_lo": np.nan, "ror_ci_hi": np.nan}

    ror = (a * d) / (b * c)
    se_log = sqrt(1/a + 1/b + 1/c + 1/d)
    log_ror = np.log(ror)
    ci_lo = np.exp(log_ror - 1.96 * se_log)
    ci_hi = np.exp(log_ror + 1.96 * se_log)

    return {"ror": ror, "ror_ci_lo": ci_lo, "ror_ci_hi": ci_hi}


def compute_ic(a: int, b: int, c: int, d: int) -> dict:
    """Information Component (BCPNN) with IC025 (lower 2.5% credible bound).

    Uses the simplified BCPNN approximation (Bate et al. 1998).
    """
    n_total = a + b + c + d
    if n_total == 0 or a == 0:
        return {"ic": np.nan, "ic025": np.nan}

    n_drug = a + b
    n_event = a + c
    expected = (n_drug * n_event) / n_total

    if expected == 0:
        return {"ic": np.nan, "ic025": np.nan}

    # IC = log2(observed / expected)
    ic = log2(a / expected)

    # Approximate variance of IC (Norén et al. 2006)
    # V(IC) ≈ 1 / (a * ln(2)^2) for the simplified version
    # More robust: use the posterior variance from the BCPNN
    # For the simplified version: SE ≈ sqrt(1/a) / ln(2)
    se_ic = sqrt(1 / a) / np.log(2) if a > 0 else np.nan
    ic025 = ic - 1.96 * se_ic if not np.isnan(se_ic) else np.nan

    return {"ic": ic, "ic025": ic025}


def compute_ebgm_simple(a: int, b: int, c: int, d: int) -> dict:
    """Simplified EBGM approximation.

    For the full MGPS with 3-component mixture hyperpriors, use the openEBGM
    R package (pre-registered). This Python implementation provides a simplified
    shrinkage estimate for initial screening.

    The simplified EBGM uses a single-component Gamma prior with parameters
    estimated from the marginals, following DuMouchel (1999).
    """
    n_total = a + b + c + d
    if n_total == 0 or a == 0:
        return {"ebgm": np.nan, "eb05": np.nan}

    n_drug = a + b
    n_event = a + c
    expected = (n_drug * n_event) / n_total

    if expected == 0:
        return {"ebgm": np.nan, "eb05": np.nan}

    # Simplified: EBGM ≈ (a + 0.5) / (expected + 0.5) as a shrunken estimate
    # This is a rough approximation; the full MGPS uses EM-estimated hyperpriors
    alpha_prior = 0.5
    beta_prior = 0.5
    ebgm = (a + alpha_prior) / (expected + beta_prior)

    # EB05: lower 5th percentile of the posterior
    # Approximate using Gamma posterior: Gamma(a + alpha, 1/(expected + beta))
    shape = a + alpha_prior
    rate = expected + beta_prior
    eb05 = stats.gamma.ppf(0.05, a=shape, scale=1/rate)

    return {"ebgm": ebgm, "eb05": eb05}


def compute_all_measures(a: int, b: int, c: int, d: int) -> dict:
    """Compute all four disproportionality measures."""
    result = {"n": a}
    result.update(compute_prr(a, b, c, d))
    result.update(compute_ror(a, b, c, d))
    result.update(compute_ic(a, b, c, d))
    result.update(compute_ebgm_simple(a, b, c, d))
    return result


def meets_conjunction(row: dict) -> bool:
    """Check if a drug-reaction pair meets the pre-registered conjunction criterion."""
    return (
        row["n"] >= MIN_N
        and not np.isnan(row.get("prr", np.nan))
        and row["prr"] >= PRR_THRESHOLD
        and row["chi2"] >= CHI2_THRESHOLD
        and row["ic025"] > IC025_THRESHOLD
        and row["eb05"] >= EB05_THRESHOLD
    )


def run_disproportionality(
    case_drugs: pd.DataFrame,
    case_reactions: pd.DataFrame,
    product_filter: pd.DataFrame | None = None,
    role_filter: str | None = None,
    label: str = "primary",
) -> pd.DataFrame:
    """Run disproportionality for all (study drug, PT) pairs.

    Args:
        case_drugs: case-drug table
        case_reactions: case-reaction table
        product_filter: if provided, restrict to these products
        role_filter: if "Suspected", restrict to Suspected-role only
        label: analysis label for output
    """
    if role_filter:
        case_drugs = case_drugs[case_drugs["drug_role"] == role_filter].copy()

    study_drugs = case_drugs[case_drugs["product_id"].notna()]
    if product_filter is not None:
        study_drugs = study_drugs[
            study_drugs["product_id"].isin(product_filter["product_id"])
        ]

    total_cases = case_drugs["case_number"].nunique()

    # Build case sets per drug and per reaction
    drug_case_sets = {}
    for prod_id, group in study_drugs.groupby("product_id"):
        drug_case_sets[prod_id] = set(group["case_number"])

    event_case_sets = {}
    for pt, group in case_reactions.groupby("reaction_pt"):
        event_case_sets[pt] = set(group["case_number"])

    # Product metadata
    prod_meta = (
        study_drugs[["product_id", "product_set", "canonical_name"]]
        .drop_duplicates()
        .set_index("product_id")
    )

    results = []
    products = sorted(drug_case_sets.keys())
    for prod_id in products:
        drug_cases = drug_case_sets[prod_id]
        meta = prod_meta.loc[prod_id]

        # Get reaction PTs observed for this drug's cases
        drug_reactions = case_reactions[
            case_reactions["case_number"].isin(drug_cases)
        ]["reaction_pt"].unique()

        for pt in drug_reactions:
            event_cases = event_case_sets.get(pt, set())
            a, b, c, d = build_2x2(drug_cases, event_cases, total_cases)

            if a < MIN_N:
                continue

            measures = compute_all_measures(a, b, c, d)
            measures["product_id"] = prod_id
            measures["product_set"] = meta["product_set"]
            measures["canonical_name"] = meta["canonical_name"]
            measures["reaction_pt"] = pt
            measures["signal"] = meets_conjunction(measures)
            measures["analysis"] = label
            results.append(measures)

    return pd.DataFrame(results)


def main():
    print("=" * 70)
    print("  Disproportionality Analysis — DAEN Tropical PV")
    print("=" * 70)

    case_drugs = pd.read_parquet(PROCESSED_DIR / "daen_case_drugs.parquet")
    case_reactions = pd.read_parquet(PROCESSED_DIR / "daen_case_reactions.parquet")
    product_dict = pd.read_csv(PROCESSED_DIR / "product_dictionary_prespec.csv")

    print(f"\nCase-drugs: {len(case_drugs):,}")
    print(f"Case-reactions: {len(case_reactions):,}")
    print(f"Study products: {len(product_dict)}")

    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # ── Primary analysis (all-role) ─────────────────────────────────────
    print("\n--- Primary analysis (all-role) ---")
    primary = run_disproportionality(
        case_drugs, case_reactions, product_dict, role_filter=None, label="primary_all_role"
    )
    n_signals = primary["signal"].sum()
    print(f"Pairs with n >= {MIN_N}: {len(primary):,}")
    print(f"Signals (conjunction criterion): {n_signals}")

    if n_signals > 0:
        print("\nSignals:")
        signals = primary[primary["signal"]].sort_values("prr", ascending=False)
        for _, row in signals.iterrows():
            print(
                f"  {row['canonical_name']:35s}  {row['reaction_pt']:40s}  "
                f"n={row['n']:>3}  PRR={row['prr']:.1f}  IC025={row['ic025']:.2f}  "
                f"EB05={row['eb05']:.2f}"
            )

    # ── Drug-role stratification (H2, Amendment 1) ──────────────────────
    print("\n--- Drug-role stratification (Suspected-only) ---")
    suspected = run_disproportionality(
        case_drugs, case_reactions, product_dict,
        role_filter="Suspected", label="suspected_only"
    )
    n_signals_susp = suspected["signal"].sum()
    print(f"Pairs with n >= {MIN_N}: {len(suspected):,}")
    print(f"Signals (conjunction criterion): {n_signals_susp}")

    # Identify masking candidates: Suspected-only signal=True AND all-role signal=False
    merged_h2 = primary.merge(
        suspected,
        on=["product_id", "reaction_pt"],
        suffixes=("_all", "_susp"),
        how="outer",
    )
    masking = merged_h2[
        (merged_h2["signal_susp"] == True) &
        (merged_h2["signal_all"] != True)
    ]
    print(f"\nDrug-role masking candidates (H2): {len(masking)}")
    if len(masking) > 0:
        for _, row in masking.iterrows():
            print(
                f"  {row['canonical_name_susp']:35s}  {row['reaction_pt']:40s}  "
                f"Susp IC025={row['ic025_susp']:.2f}  All IC025={row.get('ic025_all', np.nan):.2f}"
            )

    # Also check IC025-based masking (pre-registered definition)
    ic_masking = merged_h2[
        (merged_h2["ic025_susp"] > 0) &
        (merged_h2["ic025_all"].fillna(0) <= 0)
    ]
    print(f"IC025-based masking candidates: {len(ic_masking)}")

    # ── Validation controls ─────────────────────────────────────────────
    print("\n--- Validation controls ---")
    controls = {
        "positive": [
            ("Polyvalent snake antivenom", "Anaphylactic reaction"),
            ("Doxycycline", "Photosensitivity reaction"),
            ("Primaquine", "Methaemoglobinaemia"),
        ],
        "negative": [
            ("Box jellyfish antivenom", "Alopecia"),
            ("Praziquantel", "Myocardial infarction"),
        ],
    }

    control_results = []
    for control_type, pairs in controls.items():
        for drug_name, pt in pairs:
            match = primary[
                (primary["canonical_name"] == drug_name) &
                (primary["reaction_pt"] == pt)
            ]
            if len(match) == 1:
                row = match.iloc[0]
                result = {
                    "control_type": control_type,
                    "drug": drug_name,
                    "reaction_pt": pt,
                    "n": row["n"],
                    "prr": row["prr"],
                    "ic025": row["ic025"],
                    "eb05": row["eb05"],
                    "signal": row["signal"],
                    "status": "PASS" if (control_type == "positive" and row["signal"])
                             or (control_type == "negative" and not row["signal"])
                             else "FAIL",
                }
            else:
                result = {
                    "control_type": control_type,
                    "drug": drug_name,
                    "reaction_pt": pt,
                    "n": 0,
                    "prr": np.nan,
                    "ic025": np.nan,
                    "eb05": np.nan,
                    "signal": False,
                    "status": "N<3" if control_type == "positive" else "PASS (absent)",
                }
            control_results.append(result)
            status = result["status"]
            n_val = result["n"]
            print(f"  [{control_type:8s}] {drug_name:35s} × {pt:30s}  n={n_val:>3}  {status}")

    # ── Save outputs ────────────────────────────────────────────────────
    signals_df = primary[primary["signal"]].sort_values(
        ["product_set", "canonical_name", "prr"], ascending=[True, True, False]
    )
    signals_df.to_csv(RESULTS_DIR / "signals_primary.csv", index=False)

    primary.sort_values(
        ["product_set", "canonical_name", "n"], ascending=[True, True, False]
    ).to_csv(RESULTS_DIR / "signals_full_scan.csv", index=False)

    if len(masking) > 0:
        masking.to_csv(RESULTS_DIR / "drug_role_stratification.csv", index=False)
    else:
        ic_masking.to_csv(RESULTS_DIR / "drug_role_stratification.csv", index=False)

    pd.DataFrame(control_results).to_csv(RESULTS_DIR / "validation_controls.csv", index=False)

    # Summary
    summary = [
        "DISPROPORTIONALITY ANALYSIS SUMMARY",
        "=" * 50,
        f"Total DAEN cases (denominator): {case_drugs['case_number'].nunique():,}",
        f"Study drug-reaction pairs analysed (n >= {MIN_N}): {len(primary):,}",
        f"Signals (conjunction criterion): {n_signals}",
        f"",
        f"Drug-role stratification (Suspected-only):",
        f"  Pairs analysed: {len(suspected):,}",
        f"  Signals: {n_signals_susp}",
        f"  Masking candidates (conjunction): {len(masking)}",
        f"  Masking candidates (IC025): {len(ic_masking)}",
        f"",
        f"Validation controls: {sum(1 for r in control_results if 'PASS' in r['status'])}/{len(control_results)} passed",
    ]
    summary_text = "\n".join(summary)
    (RESULTS_DIR / "disproportionality_summary.txt").write_text(summary_text)

    print(f"\n{summary_text}")
    print(f"\nSaved to: {RESULTS_DIR}/")
    print(f"\n{'=' * 70}")
    print(f"  Disproportionality analysis complete.")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
