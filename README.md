# Disproportionality and Temporal Analysis of Australian Antivenoms, Tropical Medicines and Japanese Encephalitis Vaccine in the TGA DAEN

Code repository for two papers built on the same pre-registered DAEN pipeline:

1. **Drug-Role-Stratified Disproportionality Analysis of Australian Antivenoms and Tropical Medicines in the TGA Database of Adverse Event Notifications: A Pre-Registered READUS-PV-Compliant Framework** (scripts 01–09)
2. **Vaccination errors surface during outbreak-driven mass Japanese encephalitis vaccination: a within-product temporal analysis of Australian adverse event reports** (scripts 01–02, 07, 10–13; see [JEV vaccine temporal analysis](#jev-vaccine-temporal-analysis))

Hayden Farquhar MBBS MPHTM, Independent researcher, Finley, NSW, Australia

[![ORCID](https://img.shields.io/badge/ORCID-0009--0002--6226--440X-green)](https://orcid.org/0009-0002-6226-440X)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20404265.svg)](https://doi.org/10.5281/zenodo.20404265)

- **Pre-registration:** https://doi.org/10.17605/OSF.IO/TWGX4
- **Code archive:** https://doi.org/10.5281/zenodo.20404265
- **Preprint (antivenom and tropical-medicine paper):** https://doi.org/10.5281/zenodo.20411526

## Overview

This repository contains the analysis code to reproduce a READUS-PV 2024-compliant disproportionality analysis of the TGA Database of Adverse Event Notifications (DAEN). The study examines adverse event signals for Australian antivenoms (8 products) and a tropical-medicine basket (7 products) using PRR, ROR, IC (BCPNN), and EBGM measures with a pre-specified conjunction criterion. Analyses include drug-role stratification for masked signal detection, PELT change-point detection for the 2022 JEV vaccine rollout, eight pre-registered sensitivity analyses, and blinded clinical face-validity review.

## Data Sources

| Source | URL | Access | License |
|--------|-----|--------|---------|
| TGA DAEN (Medicines) | https://aems.tga.gov.au | Free, no registration | Crown copyright, attribution required |
| FDA FAERS (openFDA) | https://api.fda.gov/drug/event.json | Free API, no key required | Public domain |

Raw DAEN data cannot be redistributed. See [`data/raw/README.md`](data/raw/README.md) for step-by-step acquisition instructions.

## Requirements

- Python ≥ 3.11
- All dependencies pinned in `requirements.txt`

```bash
pip install -r requirements.txt
```

## Reproduction

Run all scripts from the repository root in order. Total runtime is approximately 15–20 minutes on a modern laptop (CPU only).

```bash
# Step 1: Parse raw DAEN exports into merged parquet
#         (requires DAEN xlsx files in data/raw/ — see data/raw/README.md)
python scripts/01_parse_daen.py

# Step 2: Clean data — parse medicines/reactions fields, match study drugs
python scripts/02_clean_daen.py

# Step 3: Build product dictionary matches
python scripts/03_build_product_dict.py

# Step 4: Validate reaction dictionary
python scripts/04_build_reaction_dict.py

# Step 5: Fetch FAERS comparator data (requires internet)
python scripts/05_parse_faers_comparator.py

# Step 6: Primary disproportionality analysis + drug-role stratification + validation controls
python scripts/06_disproportionality.py

# Step 7: Temporal analysis (JEV change-point, antivenom descriptive)
python scripts/07_temporal.py

# Step 8: Sensitivity analyses S1–S8
python scripts/08_sensitivity.py

# Step 9: Clinical face-validity review
python scripts/09_clinical_review.py
```

## JEV vaccine temporal analysis

The Japanese encephalitis (JEV) vaccine paper uses the shared parsing and cleaning steps plus the temporal module, followed by four JEV-specific scripts:

```bash
python scripts/01_parse_daen.py              # requires DAEN xlsx files in data/raw/
python scripts/02_clean_daen.py
python scripts/07_temporal.py                # quarterly series, PELT, Poisson break, PT shift
python scripts/10_temporal_overdispersion.py # dispersion + quasi-Poisson / negative-binomial refits; per-quarter rates
python scripts/11_temporal_subperiod.py      # post-2022 sub-period composition (errors vs reactions)
python scripts/12_figure_jev_temporal.py     # Figure 1
python scripts/13_figure_jev_pt_shift.py     # Figure 2
```

The aggregate output tables these steps produce are committed under `outputs/tables/temporal_jev_*.csv`, and the two figures under `outputs/figures/`. Scripts 10, 12 and 13 read only those aggregate tables, so the model refits and both figures can be regenerated without the raw DAEN exports:

```bash
python scripts/10_temporal_overdispersion.py
python scripts/12_figure_jev_temporal.py
python scripts/13_figure_jev_pt_shift.py
```

Script 11 reads the case-level processed parquet files and therefore needs scripts 01–02 to have been run first. Script 13 checks every count it plots against `temporal_jev_pt_shift.csv` and the 136 pre-2022 / 18 post-2022 quarter split against `temporal_jev_quarterly.csv`, and stops if either differs. Running the full sequence from the 25 May 2026 DAEN export reproduces the committed tables byte-for-byte and the committed figures pixel-for-pixel.

## Script Descriptions

| Script | Description | Inputs | Outputs |
|--------|-------------|--------|---------|
| `01_parse_daen.py` | Parse raw DAEN xlsx exports, merge, deduplicate, compute file hashes | `data/raw/*.xlsx` | `data/processed/daen_merged.parquet`, `data_snapshot_log.json` |
| `02_clean_daen.py` | Parse multi-value medicines/reactions fields into normalised tables, match study drugs | `daen_merged.parquet`, `product_dictionary_prespec.csv` | `daen_cases.parquet`, `daen_case_drugs.parquet`, `daen_case_reactions.parquet` |
| `03_build_product_dict.py` | Match DAEN records to pre-specified product sets via string search | `daen_merged.parquet`, `product_dictionary_prespec.csv` | `product_dictionary.csv`, `daen_study_drugs.parquet` |
| `04_build_reaction_dict.py` | Validate reaction PT vocabulary and build frequency dictionary | `daen_case_reactions.parquet` | `reaction_dictionary.csv`, `reaction_validation.txt` |
| `05_parse_faers_comparator.py` | Fetch FAERS reaction profiles for Set B drugs via openFDA API; confirm antivenom comparator-absence | (API) | `faers_tropical_basket.parquet`, `faers_antivenom_check.json` |
| `06_disproportionality.py` | Compute PRR, ROR, IC, EBGM for all study drug–PT pairs; apply conjunction criterion; drug-role stratification; validation controls | `daen_case_drugs.parquet`, `daen_case_reactions.parquet`, `product_dictionary_prespec.csv` | `signals_primary.csv`, `signals_full_scan.csv`, `drug_role_stratification.csv`, `validation_controls.csv` |
| `07_temporal.py` | JEV quarterly aggregation, PELT change-point detection, Poisson structural break, time-period-stratified disproportionality, antivenom descriptive time-series | `daen_cases.parquet`, `daen_case_drugs.parquet`, `daen_case_reactions.parquet` | `temporal_jev_*.csv`, `temporal_antivenom_yearly.csv`, `temporal_drugrole_yearly.csv` |
| `08_sensitivity.py` | Six feasible pre-registered sensitivity analyses (S3–S8); signal robustness comparison | `signals_full_scan.csv`, `daen_*.parquet` | `sensitivity_S*.csv`, `sensitivity_comparison.csv` |
| `09_clinical_review.py` | Blinded pharmacological plausibility rating (Established / Plausible novel / Noise) for all 68 signals | `signals_primary.csv` | `clinical_review_log.csv`, `clinical_review_summary.txt` |
| `10_temporal_overdispersion.py` | JEV structural-break model: dispersion statistics, Poisson / quasi-Poisson / negative-binomial refits, pre/post quarterly means and fold change, per-quarter rates for selected PTs | `temporal_jev_quarterly.csv`, `temporal_jev_pt_shift.csv` | `temporal_jev_poisson_overdispersion.csv`, `temporal_jev_per_quarter_rates.csv` |
| `11_temporal_subperiod.py` | Split of the post-2022 JEV period (2022–2023, 2024–2026, 2025 alone); share of PT occurrences that are administrative/programmatic errors vs classical reactions | `daen_cases.parquet`, `daen_case_drugs.parquet`, `daen_case_reactions.parquet` | `temporal_jev_subperiod_composition.csv` |
| `12_figure_jev_temporal.py` | Quarterly JEV vaccine report counts with all PELT change points | `temporal_jev_quarterly.csv`, `temporal_jev_changepoint.csv` | `outputs/figures/figure_jev_temporal.pdf`, `.png` |
| `13_figure_jev_pt_shift.py` | Per-quarter JEV reporting rates for classical reactions, administrative errors and a population-shift indicator, pre-2022 vs 2022 onwards | `temporal_jev_pt_shift.csv`, `temporal_jev_quarterly.csv` | `outputs/figures/figure_jev_pt_shift.pdf`, `.png` |

## Outputs

Antivenom and tropical-medicine paper (regenerate by running scripts 01–09):

| File | Paper reference |
|------|----------------|
| `outputs/tables/signals_primary.csv` | Table 1 (antivenom signals); Supplementary Table S1 (tropical signals) |
| `outputs/tables/validation_controls.csv` | Table 2 |
| `outputs/tables/drug_role_stratification.csv` | Supplementary Table S2 |
| `outputs/tables/sensitivity_comparison.csv` | Table 3 |
| `outputs/tables/temporal_jev_quarterly.csv` | Supplementary Figure S1 data |
| `outputs/tables/temporal_jev_changepoint.csv` | Section 3.6 (H3 results) |
| `outputs/tables/temporal_jev_poisson.csv` | Section 3.6 (Poisson IRR) |
| `outputs/tables/temporal_jev_pt_shift.csv` | Supplementary Table S5 |
| `outputs/tables/clinical_review_log.csv` | Section 3.7; Supplementary Table S1 (rating column) |

JEV vaccine paper (committed in this repository):

| File | Paper reference |
|------|----------------|
| `outputs/figures/figure_jev_temporal.pdf` | Figure 1 |
| `outputs/figures/figure_jev_pt_shift.pdf` | Figure 2 |
| `outputs/tables/temporal_jev_pt_shift.csv`, `temporal_jev_stratified_disp.csv`, `temporal_jev_per_quarter_rates.csv` | Table 1 (counts, per-quarter rates, post-2022 PRR and IC025) |
| `outputs/tables/temporal_jev_quarterly.csv` | Supplementary Table S1 (full quarterly series, 1988Q1–2026Q2) |
| `outputs/tables/temporal_jev_changepoint.csv` | Supplementary Table S2 |
| `outputs/tables/temporal_jev_poisson_overdispersion.csv` (Poisson row also in `temporal_jev_poisson.csv`) | Supplementary Table S3 |
| `outputs/tables/temporal_jev_pt_shift.csv` | Supplementary Table S4 (full PT-level table, all 187 PTs) |
| `outputs/tables/temporal_jev_stratified_disp.csv` | Supplementary Table S5 |
| `outputs/tables/temporal_jev_subperiod_composition.csv` | Supplementary Table S7 |

The committed JEV tables contain aggregate report counts and derived statistics only, computed from publicly released DAEN data (source: Therapeutic Goods Administration, https://aems.tga.gov.au).

## Citation

If you use this code, please cite:

```
Farquhar H. Drug-Role-Stratified Disproportionality Analysis of Australian
Antivenoms and Tropical Medicines in the TGA Database of Adverse Event
Notifications: A Pre-Registered READUS-PV-Compliant Framework. 2026.
Code: https://doi.org/10.5281/zenodo.20404265
Pre-registration: https://doi.org/10.17605/OSF.IO/TWGX4
```

For the JEV vaccine analysis:

```
Farquhar H. Vaccination errors surface during outbreak-driven mass Japanese
encephalitis vaccination: a within-product temporal analysis of Australian
adverse event reports. 2026.
Code: https://doi.org/10.5281/zenodo.20404265
```

## License

Code: [MIT License](LICENSE)

The pre-specified product dictionary (`data/processed/product_dictionary_prespec.csv`) is released under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).
