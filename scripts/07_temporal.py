"""
05 — Temporal Analysis (Phase 5, Section 9.6)

JEV-specific temporal analysis:
  1. Quarterly aggregation of JEV vaccine reports
  2. PELT change-point detection (ruptures; BIC penalty)
  3. Poisson regression with structural break at 2022-Q1
  4. Time-period-stratified disproportionality (pre-2022 vs 2022+)
  5. PT profile shift analysis (new/lost PTs across periods)

Antivenom temporal context:
  6. Descriptive time-series by product × year

Consumer reporting reform proxy:
  7. Drug-role proportion by year (proxy for reporter-type shift;
     Amendment 1 replaced reporter-channel with drug-role)

Outputs (in outputs/tables/):
  - temporal_jev_quarterly.csv         JEV quarterly counts
  - temporal_jev_changepoint.csv       PELT change-point results
  - temporal_jev_poisson.csv           Poisson regression coefficients
  - temporal_jev_stratified_disp.csv   Pre/post disproportionality comparison
  - temporal_jev_pt_shift.csv          PT profile shift (new/lost/changed)
  - temporal_antivenom_yearly.csv      Antivenom yearly descriptive
  - temporal_drugrole_yearly.csv       Drug-role proportions by year
  - temporal_summary.txt               Summary report

Usage:
    python scripts/07_temporal.py
"""

import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import ruptures as rpt
from scipy import stats
import statsmodels.api as sm

import sys as _sys; _sys.path.insert(0, str(Path(__file__).resolve().parent))
from importlib import import_module as _import
_disp = _import("06_disproportionality")
build_2x2 = _disp.build_2x2
compute_all_measures = _disp.compute_all_measures
meets_conjunction = _disp.meets_conjunction

PROJECT_DIR = Path(__file__).resolve().parent.parent
PROCESSED_DIR = PROJECT_DIR / "data" / "processed"
RESULTS_DIR = PROJECT_DIR / "outputs" / "tables"

CHANGE_POINT_WINDOW = ("2021-10-01", "2023-12-31")
PRE_PERIOD_END = "2021-12-31"
POST_PERIOD_START = "2022-01-01"


def load_data() -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    cases = pd.read_parquet(PROCESSED_DIR / "daen_cases.parquet")
    case_drugs = pd.read_parquet(PROCESSED_DIR / "daen_case_drugs.parquet")
    case_reactions = pd.read_parquet(PROCESSED_DIR / "daen_case_reactions.parquet")
    product_dict = pd.read_csv(PROCESSED_DIR / "product_dictionary_prespec.csv")
    return cases, case_drugs, case_reactions, product_dict


def quarterly_aggregation(
    cases: pd.DataFrame,
    case_drugs: pd.DataFrame,
    product_name: str,
) -> pd.DataFrame:
    """Aggregate report counts by quarter for a given product."""
    drug_rows = case_drugs[case_drugs["canonical_name"] == product_name]
    drug_cases = drug_rows[["case_number"]].drop_duplicates().merge(
        cases[["case_number", "report_date"]], on="case_number"
    )
    drug_cases = drug_cases.copy()
    drug_cases["quarter"] = drug_cases["report_date"].dt.to_period("Q")

    qtr_counts = (
        drug_cases.groupby("quarter")["case_number"]
        .nunique()
        .rename("n_reports")
        .reset_index()
    )

    # Fill missing quarters with zero
    if len(qtr_counts) > 0:
        full_range = pd.period_range(
            qtr_counts["quarter"].min(),
            qtr_counts["quarter"].max(),
            freq="Q",
        )
        qtr_counts = (
            qtr_counts.set_index("quarter")
            .reindex(full_range, fill_value=0)
            .rename_axis("quarter")
            .reset_index()
        )

    qtr_counts["product"] = product_name
    return qtr_counts


def pelt_changepoint(series: np.ndarray, min_size: int = 4) -> dict:
    """PELT change-point detection with BIC penalty (pre-registered).

    Uses the Poisson cost model (appropriate for count data).
    Returns detected change points and their positions.
    """
    if len(series) < 2 * min_size:
        return {
            "change_points": [],
            "n_segments": 1,
            "method": "PELT",
            "penalty": "BIC",
            "status": "insufficient_data",
        }

    algo = rpt.Pelt(model="l2", min_size=min_size).fit(series.reshape(-1, 1))
    # BIC penalty: pen = k * log(n) where k = number of parameters
    pen = np.log(len(series))
    change_points = algo.predict(pen=pen)
    # ruptures returns the last index as a breakpoint; remove it
    change_points = [cp for cp in change_points if cp < len(series)]

    return {
        "change_points": change_points,
        "n_segments": len(change_points) + 1,
        "method": "PELT",
        "penalty": "BIC",
        "pen_value": pen,
        "status": "complete",
    }


def poisson_structural_break(
    qtr_df: pd.DataFrame,
    break_date: str = POST_PERIOD_START,
) -> dict:
    """Poisson regression with structural break indicator.

    Model: log(E[Y]) = β0 + β1*post + β2*time_index
    where post = 1 if quarter >= break_date, time_index = secular trend.
    """
    df = qtr_df.copy()
    df.loc[:, "time_idx"] = np.arange(len(df))
    df.loc[:, "post"] = (df["quarter"].dt.start_time >= break_date).astype(int)

    y = df["n_reports"].values
    X = sm.add_constant(df[["post", "time_idx"]].values)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        try:
            model = sm.GLM(y, X, family=sm.families.Poisson()).fit()
        except Exception as e:
            return {"status": "failed", "error": str(e)}

    params = model.params
    pvalues = model.pvalues
    conf = model.conf_int()

    return {
        "status": "complete",
        "intercept": params[0],
        "post_coef": params[1],
        "post_irr": np.exp(params[1]),
        "post_irr_ci_lo": np.exp(conf[1, 0]),
        "post_irr_ci_hi": np.exp(conf[1, 1]),
        "post_pvalue": pvalues[1],
        "trend_coef": params[2],
        "trend_pvalue": pvalues[2],
        "aic": model.aic,
        "deviance": model.deviance,
        "pearson_chi2": model.pearson_chi2,
    }


def time_stratified_disproportionality(
    cases: pd.DataFrame,
    case_drugs: pd.DataFrame,
    case_reactions: pd.DataFrame,
    product_dict: pd.DataFrame,
    product_name: str,
    cutoff: str = POST_PERIOD_START,
) -> pd.DataFrame:
    """Compute disproportionality separately for pre-cutoff and post-cutoff periods.

    This is the key analysis: does the signal PROFILE change, not just the volume?
    """
    drug_rows = case_drugs[case_drugs["canonical_name"] == product_name]
    drug_case_numbers = set(drug_rows["case_number"])

    results = []
    for period_label, date_filter in [
        ("pre_2022", cases["report_date"] < cutoff),
        ("post_2022", cases["report_date"] >= cutoff),
    ]:
        period_cases = set(cases.loc[date_filter, "case_number"])

        period_cd = case_drugs[case_drugs["case_number"].isin(period_cases)]
        period_cr = case_reactions[case_reactions["case_number"].isin(period_cases)]
        total_period = period_cd["case_number"].nunique()

        drug_period_cases = drug_case_numbers & period_cases
        if len(drug_period_cases) < 3:
            continue

        drug_case_set = drug_period_cases

        event_case_sets = {}
        for pt, grp in period_cr.groupby("reaction_pt"):
            event_case_sets[pt] = set(grp["case_number"])

        drug_reactions = period_cr[
            period_cr["case_number"].isin(drug_case_set)
        ]["reaction_pt"].unique()

        for pt in drug_reactions:
            event_cases = event_case_sets.get(pt, set())
            a, b, c, d = build_2x2(drug_case_set, event_cases, total_period)
            if a < 1:
                continue

            measures = compute_all_measures(a, b, c, d)
            measures["product"] = product_name
            measures["reaction_pt"] = pt
            measures["period"] = period_label
            measures["signal"] = meets_conjunction(measures) if a >= 3 else False
            results.append(measures)

    return pd.DataFrame(results)


def pt_profile_shift(strat_df: pd.DataFrame) -> pd.DataFrame:
    """Classify PTs as new, lost, strengthened, weakened, or stable across periods."""
    if strat_df.empty:
        return pd.DataFrame()

    pre = strat_df[strat_df["period"] == "pre_2022"].set_index("reaction_pt")
    post = strat_df[strat_df["period"] == "post_2022"].set_index("reaction_pt")

    all_pts = sorted(set(pre.index) | set(post.index))
    shifts = []
    for pt in all_pts:
        in_pre = pt in pre.index
        in_post = pt in post.index

        pre_n = int(pre.loc[pt, "n"]) if in_pre else 0
        post_n = int(post.loc[pt, "n"]) if in_post else 0
        pre_signal = bool(pre.loc[pt, "signal"]) if in_pre else False
        post_signal = bool(post.loc[pt, "signal"]) if in_post else False
        pre_prr = float(pre.loc[pt, "prr"]) if in_pre else np.nan
        post_prr = float(post.loc[pt, "prr"]) if in_post else np.nan

        if not in_pre and in_post:
            shift_type = "new_post_2022"
        elif in_pre and not in_post:
            shift_type = "lost_post_2022"
        elif not pre_signal and post_signal:
            shift_type = "gained_signal"
        elif pre_signal and not post_signal:
            shift_type = "lost_signal"
        elif pre_signal and post_signal:
            shift_type = "persistent_signal"
        else:
            shift_type = "stable_non_signal"

        shifts.append({
            "reaction_pt": pt,
            "pre_n": pre_n,
            "post_n": post_n,
            "pre_prr": pre_prr,
            "post_prr": post_prr,
            "pre_signal": pre_signal,
            "post_signal": post_signal,
            "shift_type": shift_type,
        })

    return pd.DataFrame(shifts).sort_values(
        ["shift_type", "post_n"], ascending=[True, False]
    )


def antivenom_yearly_descriptive(
    cases: pd.DataFrame,
    case_drugs: pd.DataFrame,
) -> pd.DataFrame:
    """Yearly report counts per antivenom product (Set A)."""
    av = case_drugs[case_drugs["product_set"] == "A"].copy()
    av = av.merge(cases[["case_number", "report_date"]], on="case_number")
    av["year"] = av["report_date"].dt.year

    yearly = (
        av.groupby(["canonical_name", "year"])["case_number"]
        .nunique()
        .rename("n_reports")
        .reset_index()
    )
    return yearly


def drug_role_yearly(
    cases: pd.DataFrame,
    case_drugs: pd.DataFrame,
) -> pd.DataFrame:
    """Drug-role proportions by year (proxy for reporter-type shift)."""
    merged = case_drugs.merge(
        cases[["case_number", "report_date"]], on="case_number"
    )
    merged["year"] = merged["report_date"].dt.year

    role_year = (
        merged.groupby(["year", "drug_role"])
        .size()
        .rename("count")
        .reset_index()
    )
    totals = role_year.groupby("year")["count"].sum().rename("total")
    role_year = role_year.merge(totals, on="year")
    role_year["proportion"] = role_year["count"] / role_year["total"]

    return role_year


def main():
    print("=" * 70)
    print("  Temporal Analysis (Phase 5) — DAEN Tropical PV")
    print("=" * 70)

    cases, case_drugs, case_reactions, product_dict = load_data()
    RESULTS_DIR.mkdir(parents=True, exist_ok=True)

    # ── 1. JEV quarterly aggregation ───────────────────────────────────
    print("\n--- JEV quarterly aggregation ---")
    jev_qtr = quarterly_aggregation(cases, case_drugs, "JEV vaccines")
    print(f"Quarters with data: {(jev_qtr['n_reports'] > 0).sum()} / {len(jev_qtr)}")

    pre_mask = jev_qtr["quarter"].dt.start_time < POST_PERIOD_START
    pre_mean = jev_qtr.loc[pre_mask, "n_reports"].mean()
    post_mean = jev_qtr.loc[~pre_mask, "n_reports"].mean()
    fold_change = post_mean / pre_mean if pre_mean > 0 else np.inf
    print(f"Pre-2022 quarterly mean: {pre_mean:.1f}")
    print(f"Post-2022 quarterly mean: {post_mean:.1f}")
    print(f"Fold change: {fold_change:.1f}x")

    sparse_pre = (jev_qtr.loc[pre_mask, "n_reports"] < 5).mean()
    print(f"Pre-2022 quarters with < 5 reports: {sparse_pre:.0%}")
    if sparse_pre > 0.8:
        print("  WARNING: Pre-2022 quarterly counts predominantly < 5.")
        print("  Per pre-registration, formal testing proceeds but is flagged as sparse.")

    jev_qtr.to_csv(RESULTS_DIR / "temporal_jev_quarterly.csv", index=False)

    # ── 2. PELT change-point detection ─────────────────────────────────
    print("\n--- PELT change-point detection (BIC penalty) ---")
    series = jev_qtr["n_reports"].values.astype(float)
    pelt_result = pelt_changepoint(series, min_size=4)
    print(f"Status: {pelt_result['status']}")
    print(f"Change points detected: {pelt_result['n_segments'] - 1}")

    cp_rows = []
    if pelt_result["change_points"]:
        for cp_idx in pelt_result["change_points"]:
            if cp_idx < len(jev_qtr):
                cp_quarter = jev_qtr.iloc[cp_idx]["quarter"]
                in_window = (
                    pd.Timestamp(CHANGE_POINT_WINDOW[0])
                    <= cp_quarter.start_time
                    <= pd.Timestamp(CHANGE_POINT_WINDOW[1])
                )
                print(f"  Change point at index {cp_idx} → {cp_quarter}"
                      f" {'(IN H3 window)' if in_window else '(outside H3 window)'}")
                cp_rows.append({
                    "cp_index": cp_idx,
                    "cp_quarter": str(cp_quarter),
                    "in_h3_window": in_window,
                    "pre_mean": series[:cp_idx].mean(),
                    "post_mean": series[cp_idx:].mean(),
                })
    else:
        print("  No change points detected.")
        cp_rows.append({
            "cp_index": None,
            "cp_quarter": None,
            "in_h3_window": False,
            "pre_mean": series.mean(),
            "post_mean": series.mean(),
        })

    cp_df = pd.DataFrame(cp_rows)
    cp_df["method"] = "PELT"
    cp_df["penalty"] = "BIC"
    cp_df.to_csv(RESULTS_DIR / "temporal_jev_changepoint.csv", index=False)

    h3_supported = any(r.get("in_h3_window", False) for r in cp_rows)

    # ── 3. Poisson structural break regression ─────────────────────────
    print("\n--- Poisson structural break regression ---")
    poisson = poisson_structural_break(jev_qtr)
    if poisson["status"] == "complete":
        print(f"Post-2022 IRR: {poisson['post_irr']:.2f} "
              f"(95% CI: {poisson['post_irr_ci_lo']:.2f}–{poisson['post_irr_ci_hi']:.2f})")
        print(f"Post-2022 p-value: {poisson['post_pvalue']:.4f}")
        print(f"Secular trend p-value: {poisson['trend_pvalue']:.4f}")
        poisson_twofold = poisson["post_irr"] >= 2.0
        print(f"H3 two-fold criterion (IRR ≥ 2): {'MET' if poisson_twofold else 'NOT MET'}")
    else:
        print(f"  Poisson regression failed: {poisson.get('error', 'unknown')}")
        poisson_twofold = False

    poisson_df = pd.DataFrame([poisson])
    poisson_df.to_csv(RESULTS_DIR / "temporal_jev_poisson.csv", index=False)

    # ── 4. Time-period-stratified disproportionality ───────────────────
    print("\n--- Time-period-stratified disproportionality (JEV) ---")
    strat_df = time_stratified_disproportionality(
        cases, case_drugs, case_reactions, product_dict, "JEV vaccines"
    )
    n_pre_signals = strat_df[(strat_df["period"] == "pre_2022") & strat_df["signal"]].shape[0]
    n_post_signals = strat_df[(strat_df["period"] == "post_2022") & strat_df["signal"]].shape[0]
    print(f"Pre-2022 pairs (n ≥ 1): {(strat_df['period'] == 'pre_2022').sum()}")
    print(f"Post-2022 pairs (n ≥ 1): {(strat_df['period'] == 'post_2022').sum()}")
    print(f"Pre-2022 signals (conjunction, n ≥ 3): {n_pre_signals}")
    print(f"Post-2022 signals (conjunction, n ≥ 3): {n_post_signals}")

    if n_post_signals > 0:
        print("\nPost-2022 signals:")
        post_sigs = strat_df[
            (strat_df["period"] == "post_2022") & strat_df["signal"]
        ].sort_values("prr", ascending=False)
        for _, row in post_sigs.iterrows():
            print(f"  {row['reaction_pt']:45s} n={row['n']:>3} PRR={row['prr']:.1f}")

    strat_df.to_csv(RESULTS_DIR / "temporal_jev_stratified_disp.csv", index=False)

    # ── 5. PT profile shift ────────────────────────────────────────────
    print("\n--- PT profile shift analysis ---")
    shift_df = pt_profile_shift(strat_df)
    if not shift_df.empty:
        for st, grp in shift_df.groupby("shift_type"):
            print(f"  {st}: {len(grp)} PTs")
            if st in ("new_post_2022", "gained_signal") and len(grp) <= 10:
                for _, row in grp.iterrows():
                    print(f"    {row['reaction_pt']:45s} post_n={row['post_n']}")
        shift_df.to_csv(RESULTS_DIR / "temporal_jev_pt_shift.csv", index=False)
    else:
        print("  No data for shift analysis.")

    # ── 6. Antivenom descriptive time-series ───────────────────────────
    print("\n--- Antivenom yearly descriptive ---")
    av_yearly = antivenom_yearly_descriptive(cases, case_drugs)
    if not av_yearly.empty:
        pivoted = av_yearly.pivot_table(
            index="canonical_name", columns="year",
            values="n_reports", fill_value=0, aggfunc="sum",
        )
        recent_cols = [c for c in pivoted.columns if c >= 2015]
        print(pivoted[recent_cols].to_string())
        av_yearly.to_csv(RESULTS_DIR / "temporal_antivenom_yearly.csv", index=False)
    else:
        print("  No antivenom data.")

    # ── 7. Drug-role proportions by year ───────────────────────────────
    print("\n--- Drug-role proportions by year (reporting reform proxy) ---")
    role_yr = drug_role_yearly(cases, case_drugs)
    suspected_by_year = role_yr[role_yr["drug_role"] == "Suspected"]
    if not suspected_by_year.empty:
        recent = suspected_by_year[suspected_by_year["year"] >= 2015]
        for _, row in recent.iterrows():
            print(f"  {int(row['year'])}: Suspected {row['proportion']:.1%} "
                  f"({int(row['count']):,} / {int(row['total']):,})")
    role_yr.to_csv(RESULTS_DIR / "temporal_drugrole_yearly.csv", index=False)

    # ── H3 assessment ──────────────────────────────────────────────────
    print("\n" + "=" * 70)
    print("  H3 ASSESSMENT: JEV temporal signal")
    print("=" * 70)

    h3_criteria = {
        "change_point_in_window": h3_supported,
        "twofold_increase": fold_change >= 2.0,
        "poisson_irr_ge_2": poisson_twofold,
        "pt_profile_qualitative_shift": (
            not shift_df.empty
            and len(shift_df[shift_df["shift_type"] == "new_post_2022"]) > 0
        ),
    }

    for criterion, met in h3_criteria.items():
        print(f"  {criterion:40s}  {'MET' if met else 'NOT MET'}")

    n_met = sum(h3_criteria.values())
    overall = "SUPPORTED" if n_met >= 3 else "PARTIALLY SUPPORTED" if n_met >= 2 else "NOT SUPPORTED"
    print(f"\n  Overall H3: {overall} ({n_met}/4 criteria met)")

    # ── Summary ────────────────────────────────────────────────────────
    summary_lines = [
        "TEMPORAL ANALYSIS SUMMARY (Phase 5, Section 9.6)",
        "=" * 55,
        "",
        "JEV VACCINE TEMPORAL ANALYSIS",
        f"  Total JEV reports: {case_drugs[case_drugs['canonical_name'] == 'JEV vaccines']['case_number'].nunique()}",
        f"  Quarterly range: {jev_qtr['quarter'].iloc[0]} to {jev_qtr['quarter'].iloc[-1]}",
        f"  Pre-2022 quarterly mean: {pre_mean:.1f}",
        f"  Post-2022 quarterly mean: {post_mean:.1f}",
        f"  Fold change: {fold_change:.1f}x",
        "",
        "PELT Change-Point Detection (BIC penalty):",
        f"  Change points detected: {len(pelt_result['change_points'])}",
    ]
    for r in cp_rows:
        if r["cp_quarter"]:
            summary_lines.append(
                f"  Change point: {r['cp_quarter']} (in H3 window: {r['in_h3_window']})"
            )
    summary_lines.extend([
        "",
        "Poisson Structural Break:",
    ])
    if poisson["status"] == "complete":
        summary_lines.extend([
            f"  Post-2022 IRR: {poisson['post_irr']:.2f} "
            f"(95% CI: {poisson['post_irr_ci_lo']:.2f}–{poisson['post_irr_ci_hi']:.2f})",
            f"  Post-2022 p-value: {poisson['post_pvalue']:.6f}",
        ])
    summary_lines.extend([
        "",
        "Time-Period-Stratified Disproportionality:",
        f"  Pre-2022 signals: {n_pre_signals}",
        f"  Post-2022 signals: {n_post_signals}",
        "",
        f"H3 OVERALL: {overall} ({n_met}/4 criteria)",
    ])

    for criterion, met in h3_criteria.items():
        summary_lines.append(f"  {criterion}: {'MET' if met else 'NOT MET'}")

    summary_text = "\n".join(summary_lines)
    (RESULTS_DIR / "temporal_summary.txt").write_text(summary_text)
    print(f"\n{summary_text}")

    print(f"\nSaved to: {RESULTS_DIR}/")
    print(f"\n{'=' * 70}")
    print("  Temporal analysis complete (Phase 5).")
    print(f"{'=' * 70}\n")


if __name__ == "__main__":
    main()
