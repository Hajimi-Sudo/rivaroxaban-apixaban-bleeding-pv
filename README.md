# Rivaroxaban–Apixaban Bleeding Pharmacovigilance

Reproducible code and aggregate outputs for a multi-database pharmacovigilance study of outcome-definition-dependent uterine, menstrual, and vaginal bleeding reporting with rivaroxaban versus apixaban.

## Repository Layout

```text
.
├── src/
│   ├── faers_pilot.py
│   ├── faers_full.py
│   ├── external_validation.py
│   ├── dechallenge_followup.py
│   └── pv_stats.py
├── scripts/
│   ├── download_faers_ascii.py
│   ├── render_figures.py
│   ├── render_numeric_flow.py
│   ├── country_strata_postprocess.py
│   └── validate_revision_outputs.py
├── tests/
│   ├── test_faers_pilot.py
│   ├── test_faers_full.py
│   ├── test_external_validation.py
│   └── test_pv_stats.py
├── data/
│   ├── README.md
│   ├── raw/README.md
│   └── aggregate/
├── figures/
├── results/
├── .github/workflows/tests.yml
├── run_pipeline.sh
├── verify_release.py
├── requirements.txt
├── pyproject.toml
├── CITATION.cff
├── LICENSE
└── LICENSE-DATA.md
```

- `src/faers_pilot.py`: schema audit, FAERS deduplication, cohort construction, and pilot analysis.
- `src/faers_full.py`: complete FAERS primary, sensitivity, outcome-definition decomposition, temporal, and descriptive analyses.
- `src/external_validation.py`: separate Canada Vigilance and JADER analyses.
- `src/dechallenge_followup.py`: report-level dechallenge and rechallenge code summaries.
- `src/pv_stats.py`: ROR, PRR, Fisher test, confidence intervals, BH-FDR, and minimum detectable odds ratio.
- `scripts/download_faers_ascii.py`: discovery, download, ZIP validation, and SHA-256 manifest generation for FDA quarterly ASCII archives.
- `scripts/render_figures.py`: regeneration of manuscript Figures 2–4 from frozen aggregate tables.
- `scripts/render_numeric_flow.py`: regeneration of the numerical FAERS report-flow Figure 1.
- `scripts/country_strata_postprocess.py`: aggregate-only US/non-US stratified analysis and heterogeneity test.
- `scripts/validate_revision_outputs.py`: frozen-cell, family-size, missingness, and exact-interval validation gates.
- `data/aggregate/`: small, non-patient-level tables underlying the manuscript results and figures.
- `figures/`: publication figures in PNG format.
- `tests/`: unit tests that do not require raw pharmacovigilance data.

## Overview

The study compares spontaneous-report patterns for rivaroxaban and apixaban in FAERS, with Canada Vigilance and JADER analyzed separately as external reporting environments.

FAERS is the prespecified primary database. Canada Vigilance and JADER are not pooled with FAERS because the databases differ in reporting systems, schemas, age representation, and eligible populations.

The primary outcome is a frozen five-term coding-based definition:

1. abnormal uterine bleeding;
2. heavy menstrual bleeding;
3. intermenstrual bleeding;
4. uterine haemorrhage; and
5. vaginal haemorrhage.

The principal estimand is a reporting odds ratio, not incidence, prevalence, comparative clinical risk, or causality.

The broad FAERS contrast was positive, but narrower and term-level analyses showed pronounced outcome-definition heterogeneity. In particular, the narrow abnormal-uterine/heavy-menstrual definition reversed direction, whereas uterine-haemorrhage and vaginal-haemorrhage coding produced strong positive contrasts.

No predictive or causal machine-learning model is claimed in the associated manuscript.

## Dataset

### Sources

The analysis uses publicly released spontaneous-report datasets:

- FDA FAERS/AEMS quarterly ASCII extracts:  
  https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html
- FDA latest quarterly data page:  
  https://www.fda.gov/drugs/fda-adverse-event-monitoring-system-aems/fda-adverse-event-monitoring-system-aems-latest-quarterly-data-files
- Health Canada Canada Vigilance extract:  
  https://www.canada.ca/en/health-canada/services/drugs-health-products/medeffect-canada/adverse-reaction-database/canada-vigilance-online-database-data-extract.html
- PMDA JADER download page:  
  https://www.pmda.go.jp/safety/info-services/drugs/adr-info/suspected-adr/0004.html

### Study releases

- FAERS: 2015 Q1 through 2026 Q1.
- Canada Vigilance: full extract current through 31 March 2026.
- JADER: July 2026 public extract.

### Raw-data policy

Raw source archives are not committed to this repository.

The raw-data directories are ignored by Git because:

- the combined source data are several gigabytes;
- the databases remain available from their official custodians;
- redistribution may be governed by source-specific terms;
- public Git history is unsuitable for large mutable archives.

Place downloaded data under the structure documented in `data/raw/README.md`.

### Aggregate release

The `data/aggregate/` directory contains frozen non-patient-level outputs:

- primary (2×2) cells and effect estimates;
- cohort-flow counts;
- sensitivity-analysis cells;
- external-database cells;
- post-result outcome-definition decomposition;
- full-database background, symmetric active-comparator, reporter-source, country, and calendar-period analyses;
- conditional exact intervals for sparse tables;
- quarterly counts;
- time-to-onset summaries;
- dechallenge/rechallenge completeness and code distributions;
- serious-outcome descriptor counts.

These aggregate files contain no individual case safety reports and no patient-level data.

## Feature Engineering

### FAERS case-version resolution

Repeated primary-report identifiers are checked for agreement on case identifier and case version. The latest FDA date is retained, with archive period used only as a tie breaker. Longitudinal case identifiers are then reduced to their latest retained version.

### Age harmonization

Reported ages are converted to years from supported year, month, week, day, decade, and hour units. Unsupported units and implausible values outside 0–130 years are treated as unavailable.

### Exposure construction

Drug and active-ingredient fields are normalized to uppercase alphanumeric tokens. Prespecified generic and brand aliases are used for rivaroxaban and apixaban. The primary analysis requires mutually exclusive primary-suspect exposure.

### Outcome construction

A report contributes once to an outcome even if multiple qualifying preferred terms are present. The broad outcome, narrow sensitivity outcome, individual preferred terms, and transparent clinical groupings are encoded separately.

### External-database harmonization

Canada Vigilance uses official report identifiers, drug involvement, age, sex, and English MedDRA preferred-term fields. Explicit duplicates are collapsed; linked reports are not automatically treated as duplicates.

JADER uses database case identifiers and a categorical 20s–40s age proxy because continuous age is not consistently available.

### Descriptive fields

Time to onset is calculated only from valid exact dates. Negative and implausibly long intervals are excluded. Dechallenge, rechallenge, hospitalization, and death fields are descriptive and are not treated as causal evidence or comparative rates.

## Model Baselines

This repository does not present a predictive model.

The statistical analysis consists of:

- reporting odds ratio (ROR);
- proportional reporting ratio (PRR);
- Fisher exact test;
- log-scale 95% confidence intervals;
- Jeffreys 0.5 correction for zero cells;
- Benjamini–Hochberg false-discovery-rate adjustment;
- minimum detectable odds ratio at 80% power;
- prespecified sensitivity analyses;
- post-result explanatory outcome-definition decomposition;
- database-specific external analyses without pooling.

A previously explored temporal model failed its development gate and is intentionally excluded because it does not support any claim in the manuscript.

## Requirements

Validated scientific-analysis environment:

- Python 3.12.3
- pandas 3.0.3
- NumPy 2.5.3
- SciPy 1.17.1

Figure and test tooling:

- matplotlib 3.10.7
- pytest 8.4.2

The exact package list is in `requirements.txt`.

## Installation

Clone the repository:

```bash
git clone git@github.com:Hajimi-Sudo/rivaroxaban-apixaban-bleeding-pv.git
cd rivaroxaban-apixaban-bleeding-pv
```

Create an isolated environment:

```bash
python -m venv .venv
```

Activate it on Linux or macOS:

```bash
source .venv/bin/activate
```

Activate it on Windows PowerShell:

```powershell
.venv\Scripts\Activate.ps1
```

Install dependencies:

```bash
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
```

Run the unit tests:

```bash
python -m pytest
```

Verify the frozen aggregate release:

```bash
python verify_release.py
```

## Usage

### 1. Download FAERS/AEMS archives

From the repository root:

```bash
python scripts/download_faers_ascii.py
```

The downloader writes official archives and manifests under:

```text
data/raw/faers/archives/
data/raw/faers/manifest/
```

To run from another working directory, set `PV_PROJECT_ROOT` to the repository root.

### 2. Add Canada Vigilance and JADER

Download the official releases manually and arrange them as described in `data/raw/README.md`.

The repository does not redistribute the raw Canada Vigilance or JADER files.

### 3. Run a pilot FAERS analysis

```bash
python src/faers_pilot.py \
  --project-root . \
  --periods 2025Q1 2025Q2 2025Q3 2025Q4 \
  --output-dir results/faers_pilot
```

### 4. Run the complete FAERS analysis

```bash
python src/faers_full.py \
  --project-root . \
  --output-dir results/faers_primary
```

If `--periods` is omitted, all periods listed in the verified manifest are analyzed.

### 5. Run Canada Vigilance and JADER analyses

```bash
python src/external_validation.py \
  --project-root . \
  --output-dir results/external_validation
```

The databases are processed separately. The script does not pool records or estimates.

### 6. Run dechallenge/rechallenge summaries

```bash
python src/dechallenge_followup.py \
  --project-root . \
  --output-dir results/dechallenge_followup \
  --workers 8
```

### 7. Run the complete pipeline

After all raw datasets are placed correctly:

```bash
bash run_pipeline.sh
```

### 8. Regenerate paper figures

Figures 2–4 can be regenerated directly from the frozen aggregate tables:

```bash
python scripts/render_figures.py
```

Figure 1 can be regenerated with `python scripts/render_numeric_flow.py`; it contains only audited aggregate cohort counts.

## Key Results Summary

Primary FAERS active-comparator analysis:

- rivaroxaban: 665 event and 3,528 non-event reports;
- apixaban: 90 event and 1,895 non-event reports;
- ROR 3.969;
- 95% CI 3.162–4.981;
- PRR 3.498;
- Fisher exact (p=1.89×10^{-42});
- minimum detectable OR at 80% power: 1.410.

External reporting environments:

- Canada Vigilance ROR 4.472, 95% CI 1.674–11.948;
- JADER ROR 1.857, 95% CI 0.845–4.082.

Interpretive boundary:

- the estimates describe reporting disproportionality;
- they do not estimate incidence or prevalence;
- they do not establish treatment risk or causality;
- external databases are not pooled;
- post-result outcome-definition decomposition is explanatory.

## Output Files

### FAERS primary output

Expected files include:

| File | Purpose |
|---|---|
| `primary_contrast.json` | Primary (2×2) cells, ROR, PRR, confidence intervals, Fisher test, and MDE |
| `cohort_flow.csv` | Cohort construction and reconciliation |
| `family_A_drug_contrasts.csv` | Secondary active comparators |
| `family_C_sensitivities.csv` | Prespecified sensitivity family |
| `post_result_phenotype_decomposition.csv` | Preferred-term and fixed-group outcome-definition decomposition |
| `quarterly_active_comparator_counts.csv` | Quarterly report and outcome counts |
| `time_to_onset_summary.csv` | Exact-date time-to-onset completeness |
| `dechallenge_rechallenge_completeness.csv` | Descriptive follow-up-field completeness |
| `serious_outcome_descriptors.csv` | Non-exclusive serious-outcome counts |
| `run_metadata.json` | Runtime versions and analysis metadata |

### External output

| File | Purpose |
|---|---|
| `canada_primary_contrast.json` | Canada Vigilance primary contrast |
| `jader_primary_contrast.json` | JADER primary contrast |
| `family_B_external_replication.csv` | Database-specific external estimates |

### Released aggregate snapshot

The files under `data/aggregate/` are the frozen values underlying the manuscript tables and figures.

## Citation

The associated manuscript is:

> Phenotype-dependent reporting of uterine and vaginal bleeding with rivaroxaban versus apixaban: a cross-national pharmacovigilance study.

Citation metadata are provided in `CITATION.cff`.

Until a DOI is assigned, cite this repository by URL and commit hash:

```text
https://github.com/Hajimi-Sudo/rivaroxaban-apixaban-bleeding-pv
```

## License

Source code is released under the MIT License; see `LICENSE`.

Repository-created aggregate tables are released under CC BY 4.0; see `LICENSE-DATA.md`.

The original FAERS/AEMS, Canada Vigilance, and JADER datasets are not redistributed. Their use remains subject to the terms and guidance of the respective data custodians.
