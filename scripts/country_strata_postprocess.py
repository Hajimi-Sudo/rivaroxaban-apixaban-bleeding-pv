"""Compute US/non-US reporting-disproportionality strata from frozen aggregate counts."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from pv_stats import log_or_heterogeneity, two_by_two_statistics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--faers-dir", type=Path, required=True)
    args = parser.parse_args()

    source = pd.read_csv(args.faers_dir / "reporter_source_associations.csv")
    source = source[source["reporter_country_group"].isin(["United States", "non-US"])]
    grouped = source.groupby(
        ["drug", "event_status", "reporter_country_group"], dropna=False
    )["n"].sum()

    tables: list[tuple[str, int, int, int, int]] = []
    result_rows: list[dict[str, object]] = []
    for country in ["United States", "non-US"]:
        def get(drug: str, status: str) -> int:
            return int(grouped.get((drug, status, country), 0))

        a = get("rivaroxaban", "broad_outcome")
        b = get("rivaroxaban", "no_broad_outcome")
        c = get("apixaban", "broad_outcome")
        d = get("apixaban", "no_broad_outcome")
        row = two_by_two_statistics(a, b, c, d)
        row["stratum"] = country
        result_rows.append(row)
        tables.append((country, a, b, c, d))

    pd.DataFrame(result_rows).to_csv(args.faers_dir / "reporter_country_strata.csv", index=False)
    heterogeneity = log_or_heterogeneity(tables)
    (args.faers_dir / "reporter_country_heterogeneity.json").write_text(
        json.dumps(heterogeneity, indent=2), encoding="utf-8"
    )
    print(json.dumps({"strata": result_rows, "heterogeneity": heterogeneity}, indent=2))


if __name__ == "__main__":
    main()
