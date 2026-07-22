"""Verify the compact aggregate release without accessing raw case reports."""

from __future__ import annotations

import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent
AGGREGATE = ROOT / "data" / "aggregate"

REQUIRED = {
    "primary_contrast.json",
    "cohort_flow.csv",
    "family_A_drug_contrasts.csv",
    "family_B_external_replication.csv",
    "family_C_sensitivities.csv",
    "post_result_phenotype_decomposition.csv",
    "quarterly_active_comparator_counts.csv",
    "time_to_onset_summary.csv",
    "dechallenge_rechallenge_completeness.csv",
    "dechallenge_rechallenge_code_distribution.csv",
    "serious_outcome_descriptors.csv",
    "README.md",
}


def main() -> None:
    missing = sorted(name for name in REQUIRED if not (AGGREGATE / name).is_file())
    if missing:
        raise SystemExit(f"Missing aggregate release files: {missing}")

    primary = json.loads((AGGREGATE / "primary_contrast.json").read_text(encoding="utf-8"))
    expected_cells = {
        "a_exposed_event": 665,
        "b_exposed_non_event": 3528,
        "c_comparator_event": 90,
        "d_comparator_non_event": 1895,
    }
    for key, expected in expected_cells.items():
        if primary.get(key) != expected:
            raise SystemExit(f"Primary cell mismatch for {key}: {primary.get(key)} != {expected}")

    if not primary.get("reconciliation_pass"):
        raise SystemExit("Primary reconciliation flag is false")
    if not math.isclose(primary["ror"], 3.9688051146384478, rel_tol=0, abs_tol=1e-12):
        raise SystemExit(f"Primary ROR mismatch: {primary['ror']}")

    print("Aggregate release verification passed.")
    print("Primary cells: 665/3528 versus 90/1895")
    print(f"Primary ROR: {primary['ror']:.6f}")


if __name__ == "__main__":
    main()
