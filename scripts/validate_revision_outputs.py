"""Validate revision aggregate outputs before manuscript fill-in."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import pandas as pd


EXPECTED_PRIMARY = {
    "a_exposed_event": 665,
    "b_exposed_non_event": 3528,
    "c_comparator_event": 90,
    "d_comparator_non_event": 1895,
}


def require(condition: bool, message: str, checks: list[dict[str, object]]) -> None:
    checks.append({"check": message, "pass": bool(condition)})
    if not condition:
        raise AssertionError(message)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--faers-dir", type=Path, required=True)
    parser.add_argument("--external-dir", type=Path, required=True)
    parser.add_argument("--json-out", type=Path, required=True)
    args = parser.parse_args()
    checks: list[dict[str, object]] = []

    primary = json.loads((args.faers_dir / "primary_contrast.json").read_text(encoding="utf-8"))
    for key, expected in EXPECTED_PRIMARY.items():
        require(primary[key] == expected, f"primary {key} equals frozen value {expected}", checks)
    require(primary["reconciliation_pass"] is True, "primary cohort reconciliation passes", checks)

    background = pd.read_csv(args.faers_dir / "background_disproportionality.csv")
    secondary = pd.read_csv(args.faers_dir / "symmetric_active_comparators.csv")
    sensitivity = pd.read_csv(args.faers_dir / "family_C_sensitivities.csv")
    decomposition = pd.read_csv(args.faers_dir / "post_result_phenotype_decomposition.csv")
    periods = pd.read_csv(args.faers_dir / "calendar_period_strata.csv")
    period_het = pd.read_csv(args.faers_dir / "calendar_period_heterogeneity.csv")
    reporter = pd.read_csv(args.faers_dir / "reporter_source_interaction.csv")
    missingness = pd.read_csv(args.faers_dir / "missingness_by_drug.csv")
    characteristics = pd.read_csv(args.faers_dir / "cohort_characteristics.csv")
    sparse = pd.read_csv(args.faers_dir / "sparse_exact_intervals.csv")
    jader = pd.read_csv(args.external_dir / "jader_age_rule_sensitivity.csv")

    require(len(background) == 2, "background family contains two tests", checks)
    require(len(secondary) == 6, "symmetric comparator family contains six tests", checks)
    require(len(sensitivity) == 8, "prespecified sensitivity family contains eight tests", checks)
    require(len(decomposition) == 8, "outcome decomposition family contains eight unique tests", checks)
    require(decomposition["analysis"].nunique() == 8, "outcome decomposition has no duplicate analysis", checks)
    require(len(periods) == 3, "calendar analysis contains pre/during/post strata", checks)
    require(len(reporter) == 2, "reporter interaction contains HCP and consumer strata", checks)
    require(len(jader) == 3, "JADER analysis contains three age rules", checks)
    require((missingness["denominator_n"] > 0).all(), "all missingness denominators are positive", checks)
    require((missingness["missing_fraction"].between(0, 1)).all(), "all missingness fractions are valid", checks)
    require({"age_years", "reporter_group", "reporter_country_group", "calendar_period", "indication_group"}.issubset(set(characteristics["variable"])), "cohort characteristics cover required variables", checks)
    require((sparse["conditional_exact_ci95_lower"].notna()).all(), "all sparse tables have exact lower intervals", checks)
    require((sparse["conditional_exact_ci95_upper"].notna()).all(), "all sparse tables have exact upper intervals", checks)
    require(math.isfinite(float(period_het.loc[0, "p_heterogeneity"])), "calendar heterogeneity p value is finite", checks)
    require(math.isfinite(float(reporter.loc[0, "heterogeneity_p"])), "reporter heterogeneity p value is finite", checks)

    output = {
        "status": "pass",
        "checks": checks,
        "files_checked": sorted(str(path) for path in list(args.faers_dir.glob("*")) + list(args.external_dir.glob("*"))),
    }
    args.json_out.parent.mkdir(parents=True, exist_ok=True)
    args.json_out.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(json.dumps({"status": "pass", "check_count": len(checks)}, indent=2))


if __name__ == "__main__":
    main()
