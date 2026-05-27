# Raw Data Acquisition

The raw DAEN export files are not included in this repository due to size (~300 MB total) and TGA copyright terms. Follow these steps to obtain the data:

## TGA DAEN Export

1. Navigate to https://aems.tga.gov.au (or https://daen.tga.gov.au, which redirects)
2. Accept the terms and conditions
3. Leave the product name field **blank** (to export ALL reports — required for 2×2 table denominators)
4. Set a date range (see batches below)
5. Click "Search"
6. Click the **"List of Reports"** tab
7. Click the three-dot menu icon (⋮) → **"Data with current layout"**
8. Save the Excel file to this directory (`data/raw/`)
9. Repeat for each date range batch

### Suggested batches

Each batch should stay under 150,000 rows. Adjust if the DAEN interface imposes limits:

| Batch | Date range |
|---|---|
| 01 | 01/01/1971 – 31/12/2005 |
| 02 | 01/01/2006 – 31/12/2010 |
| 03 | 01/01/2011 – 31/12/2013 |
| 04 | 01/01/2014 – 31/12/2016 |
| 05 | 01/01/2017 – 31/12/2018 |
| 06 | 01/01/2019 – 31/12/2019 |
| 07 | 01/01/2020 – 31/12/2020 |
| 08 | 01/01/2021 – 31/12/2021 |
| 09 | 01/01/2022 – 31/12/2022 |
| 10 | 01/01/2023 – 31/12/2023 |
| 11 | 01/01/2024 – 31/12/2024 |
| 12 | 01/01/2025 – present |

The date format on the DAEN site is DD/MM/YYYY.

Any file naming convention is acceptable — `01_parse_daen.py` auto-detects `.xlsx`, `.xls`, and `.csv` files in this directory.

## FAERS Comparator

FAERS data is fetched automatically by `05_parse_faers_comparator.py` via the openFDA API. No manual download required.

## Expected output

After running `01_parse_daen.py`, you should have a merged parquet file at `data/processed/daen_merged.parquet` containing approximately 634,000 unique records (exact count depends on download date, as the DAEN is a living database).
