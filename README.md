# Drug-Role-Stratified Disproportionality Analysis of Australian Antivenoms and Tropical Medicines in the TGA DAEN

Code repository for: **Drug-Role-Stratified Disproportionality Analysis of Australian Antivenoms and Tropical Medicines in the TGA Database of Adverse Event Notifications: A Pre-Registered READUS-PV-Compliant Framework**

Hayden Farquhar MBBS MPHTM, Independent researcher, Finley, NSW, Australia

[![ORCID](https://img.shields.io/badge/ORCID-0009--0002--6226--440X-green)](https://orcid.org/0009-0002-6226-440X)

- **Pre-registration:** https://doi.org/10.17605/OSF.IO/TWGX4
- **Preprint:** to be posted

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

## Outputs

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

## Citation

If you use this code, please cite:

```
Farquhar H. Drug-Role-Stratified Disproportionality Analysis of Australian
Antivenoms and Tropical Medicines in the TGA Database of Adverse Event
Notifications: A Pre-Registered READUS-PV-Compliant Framework. 2026.
Pre-registration: https://doi.org/10.17605/OSF.IO/TWGX4
```

## License

Code: [MIT License](LICENSE)

The pre-specified product dictionary (`data/processed/product_dictionary_prespec.csv`) is released under [CC-BY 4.0](https://creativecommons.org/licenses/by/4.0/).
