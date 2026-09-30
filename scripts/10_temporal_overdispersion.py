"""
10 — JEV temporal analysis: overdispersion refits and per-quarter rates

Computes:
  1. Overdispersion: dispersion statistics + negative-binomial and
     quasi-Poisson refits of the identical structural-break model.
  2. Per-quarter rates: classical-reaction PTs do not decline per unit time;
     administrative-error PTs emerge on top.
  3. Exact pre/post quarterly means and fold change (arithmetic self-consistency).

Reuses the quarterly series and model specification of 07_temporal.py.

Inputs (outputs/tables/): temporal_jev_quarterly.csv, temporal_jev_pt_shift.csv
Outputs (outputs/tables/):
  - temporal_jev_poisson_overdispersion.csv   Supplementary Table S3
  - temporal_jev_per_quarter_rates.csv        Table 1 (per-quarter rates)

Usage:
    python scripts/10_temporal_overdispersion.py
"""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm

TABLES = Path(__file__).resolve().parent.parent / "outputs" / "tables"
BREAK = pd.Timestamp("2022-01-01")


def quarter_start(q: str) -> pd.Timestamp:
    year, qn = int(q[:4]), int(q[-1])
    return pd.Timestamp(year=year, month={1: 1, 2: 4, 3: 7, 4: 10}[qn], day=1)


def main():
    qtr = pd.read_csv(TABLES / "temporal_jev_quarterly.csv")
    qtr = qtr.assign(start=qtr["quarter"].apply(quarter_start))
    qtr = qtr.sort_values("start").reset_index(drop=True)

    y = qtr["n_reports"].values.astype(float)
    post = (qtr["start"] >= BREAK).astype(int).values
    time_idx = np.arange(len(qtr))
    X = sm.add_constant(np.column_stack([post, time_idx]))

    # --- Poisson (reproduce) ---
    pois = sm.GLM(y, X, family=sm.families.Poisson()).fit()
    irr_p = np.exp(pois.params[1])
    ci_p = np.exp(pois.conf_int()[1])
    disp_pearson = pois.pearson_chi2 / pois.df_resid
    disp_dev = pois.deviance / pois.df_resid

    # --- Quasi-Poisson (Pearson-scaled SEs) ---
    quasi = sm.GLM(y, X, family=sm.families.Poisson()).fit(scale="X2")
    irr_q = np.exp(quasi.params[1])
    ci_q = np.exp(quasi.conf_int()[1])

    # --- Negative binomial (NB2; estimate alpha via discrete NB) ---
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        nb = sm.NegativeBinomial(y, X).fit(disp=0)
    irr_nb = np.exp(nb.params[1])
    ci_nb = np.exp(nb.conf_int()[1])
    p_nb = nb.pvalues[1]

    # --- Exact pre/post quarterly means + fold change ---
    pre_mask = qtr["start"] < BREAK
    n_pre_q = int(pre_mask.sum())
    n_post_q = int((~pre_mask).sum())
    pre_mean = y[pre_mask.values].mean()
    post_mean = y[~pre_mask.values].mean()
    fold = post_mean / pre_mean

    print("=" * 64)
    print("OVERDISPERSION + REFITS (structural break: post + secular trend)")
    print("=" * 64)
    print(f"n quarters: pre={n_pre_q}  post={n_post_q}  total={len(qtr)}")
    print(f"Pre-2022 quarterly mean : {pre_mean:.3f}")
    print(f"Post-2022 quarterly mean: {post_mean:.3f}")
    print(f"Fold change (post/pre)  : {fold:.2f}x")
    print("-" * 64)
    print(f"Dispersion Pearson chi2/df : {disp_pearson:.3f}")
    print(f"Dispersion deviance/df     : {disp_dev:.3f}")
    print("-" * 64)
    print(f"Poisson       IRR {irr_p:.2f}  95% CI {ci_p[0]:.2f}-{ci_p[1]:.2f}  p={pois.pvalues[1]:.2e}")
    print(f"Quasi-Poisson IRR {irr_q:.2f}  95% CI {ci_q[0]:.2f}-{ci_q[1]:.2f}")
    print(f"Neg-binomial  IRR {irr_nb:.2f}  95% CI {ci_nb[0]:.2f}-{ci_nb[1]:.2f}  p={p_nb:.2e}  alpha={np.exp(nb.params[-1]):.3f}")

    # --- Per-quarter rates for selected PTs ---
    shift = pd.read_csv(TABLES / "temporal_jev_pt_shift.csv")
    sel = {
        "Headache": "classical",
        "Urticaria": "classical",
        "Injection site pain": "classical",
        "Incorrect route of product administration": "error",
        "Contraindication to vaccination": "error",
        "Vaccination error": "error",
    }
    print("-" * 64)
    print("PER-QUARTER RATES (reports / quarter)")
    print(f"{'PT':<45}{'pre/q':>8}{'post/q':>8}")
    rows = []
    for _, r in shift.iterrows():
        pt = r["reaction_pt"]
        if pt in sel:
            pre_rate = r["pre_n"] / n_pre_q
            post_rate = r["post_n"] / n_post_q
            print(f"{pt[:44]:<45}{pre_rate:>8.3f}{post_rate:>8.3f}")
            rows.append({
                "reaction_pt": pt, "category": sel[pt],
                "pre_n": int(r["pre_n"]), "post_n": int(r["post_n"]),
                "pre_per_quarter": round(pre_rate, 4),
                "post_per_quarter": round(post_rate, 4),
            })

    # --- Persist results tables ---
    out = pd.DataFrame([{
        "model": "Poisson", "post_irr": irr_p, "ci_lo": ci_p[0], "ci_hi": ci_p[1],
        "pvalue": pois.pvalues[1], "dispersion_pearson": disp_pearson,
        "dispersion_deviance": disp_dev, "n_pre_quarters": n_pre_q,
        "n_post_quarters": n_post_q, "pre_mean": pre_mean, "post_mean": post_mean,
        "fold_change": fold,
    }, {
        "model": "Quasi-Poisson", "post_irr": irr_q, "ci_lo": ci_q[0], "ci_hi": ci_q[1],
        "pvalue": np.nan, "dispersion_pearson": disp_pearson,
        "dispersion_deviance": disp_dev, "n_pre_quarters": n_pre_q,
        "n_post_quarters": n_post_q, "pre_mean": pre_mean, "post_mean": post_mean,
        "fold_change": fold,
    }, {
        "model": "NegativeBinomial", "post_irr": irr_nb, "ci_lo": ci_nb[0],
        "ci_hi": ci_nb[1], "pvalue": p_nb, "dispersion_pearson": disp_pearson,
        "dispersion_deviance": disp_dev, "n_pre_quarters": n_pre_q,
        "n_post_quarters": n_post_q, "pre_mean": pre_mean, "post_mean": post_mean,
        "fold_change": fold,
    }])
    TABLES.mkdir(parents=True, exist_ok=True)
    out.to_csv(TABLES / "temporal_jev_poisson_overdispersion.csv", index=False)
    pd.DataFrame(rows).to_csv(TABLES / "temporal_jev_per_quarter_rates.csv", index=False)
    print("-" * 64)
    print(f"Wrote {TABLES / 'temporal_jev_poisson_overdispersion.csv'}")
    print(f"Wrote {TABLES / 'temporal_jev_per_quarter_rates.csv'}")


if __name__ == "__main__":
    main()
