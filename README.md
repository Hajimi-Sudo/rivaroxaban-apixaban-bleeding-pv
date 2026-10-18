# Rivaroxaban–apixaban bleeding pharmacovigilance code

This repository contains the reproducible analysis code for the current pharmacovigilance study:

> Separating menstrual and anatomical bleeding phenotypes in anticoagulant safety surveillance: a multi-system pharmacovigilance and clinical-evidence triangulation study

The public repository intentionally contains code and methodological documentation only. Raw FAERS, Canada Vigilance, and JADER files, case identifiers, narratives, and patient-level outputs are not redistributed.

## Contents

- `experiment/code/`: core FAERS parsing, statistics, full-analysis, external-validation, and follow-up modules.
- `upgrade_code/`: current endpoint-taxonomy audit, low-memory FAERS pipeline, aggregate case-level complement, figure generation, and PubMed clinical-evidence search.
- `src/` and `scripts/`: reusable modules and acquisition/validation utilities retained from the analysis workflow.
- `tests/`: unit tests for the reusable analysis modules.
- `docs/`: clinical-evidence protocol and READUS-PV mapping.

## Data and execution

The scripts accept an external project root and manifest so that source archives remain outside the repository. For the 45-quarter reference audit:

```bash
python upgrade_code/low_memory_faers_counts.py \
  --project-root /path/to/runtime_project \
  --manifest /path/to/runtime_project/upgrade_2026-10-14/runtime_manifest/download_manifest_sha256.csv \
  --periods 2015Q1 2015Q2 ... 2026Q1 \
  --output-dir /path/to/results/faers_2015Q1_2026Q1
```

For the aggregate-only case-level complement:

```bash
python upgrade_code/case_level_faers_analysis.py \
  --project-root /path/to/runtime_project \
  --manifest /path/to/runtime_project/upgrade_2026-10-14/runtime_manifest/download_manifest_sha256.csv \
  --periods 2015Q1 2015Q2 ... 2026Q1 \
  --output-dir /path/to/results/case_level_faers_2015Q1_2026Q1
```

The case-level script exports aggregate descriptors only; no PRIMARYID, CASEID, narrative text, or other direct report identifier is written. All estimates are reporting-level disproportionality measures and must not be interpreted as incidence, comparative clinical risk, or causality.

The PubMed search is run with `upgrade_code/clinical_evidence_search.py`. The search export and screening protocol are documented in `docs/`.

## Verification

```bash
python -m py_compile experiment/code/*.py upgrade_code/*.py
pytest
```

The reference analysis is versioned by the commit used for reproduction; the raw-data manifest and archive checksums are kept in the private runtime workspace.
