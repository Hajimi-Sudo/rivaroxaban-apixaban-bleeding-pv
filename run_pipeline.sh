#!/usr/bin/env bash
set -euo pipefail

python src/faers_full.py   --project-root .   --output-dir results/faers_primary

python src/external_validation.py   --project-root .   --output-dir results/external_validation

python src/dechallenge_followup.py   --project-root .   --output-dir results/dechallenge_followup   --workers "${PV_WORKERS:-8}"

python scripts/render_figures.py
python verify_release.py
