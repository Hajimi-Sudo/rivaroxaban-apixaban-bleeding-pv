"""Low-memory FAERS active-comparator/outcome-group counts.

The script avoids the previous all-DEMO/all-REAC pandas concatenation. It first
indexes candidate report versions from DRUG, then retains only candidate DEMO
rows needed for latest-case resolution, and finally scans reactions for the
eligible primary IDs. It is intended for the upgrade baseline and validation,
not for unreviewed clinical inference.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pandas as pd


TERMS = {
    "menstrual_specific": {"HEAVY MENSTRUAL BLEEDING", "INTERMENSTRUAL BLEEDING"},
    "abnormal_uterine": {"ABNORMAL UTERINE BLEEDING"},
    "uterine_haemorrhage": {"UTERINE HAEMORRHAGE"},
    "vaginal_haemorrhage": {"VAGINAL HAEMORRHAGE"},
    "historical_broad_composite": {
        "ABNORMAL UTERINE BLEEDING",
        "HEAVY MENSTRUAL BLEEDING",
        "INTERMENSTRUAL BLEEDING",
        "UTERINE HAEMORRHAGE",
        "VAGINAL HAEMORRHAGE",
    },
}


def latest_row(old: dict[str, object] | None, row: dict[str, object]) -> dict[str, object]:
    if old is None:
        return row
    old_key = (int(old["fda_dt_num"]), int(old["primaryid_num"]), str(old["period"]))
    new_key = (int(row["fda_dt_num"]), int(row["primaryid_num"]), str(row["period"]))
    return row if new_key >= old_key else old


def run(project_root: Path, manifest_path: Path, periods: list[str], output_dir: Path) -> dict[str, object]:
    sys.path.insert(0, str(project_root / "experiment" / "code"))
    from faers_pilot import DRUG_ALIASES, drug_mask, load_manifest, read_table_chunks
    from pv_stats import two_by_two_statistics

    archive_dir = project_root / ".aris" / "data" / "faers" / "archives"
    records = load_manifest(manifest_path, archive_dir, periods)
    for record in records:
        if not record.path.exists():
            raise FileNotFoundError(record.path)

    # Pass 1: identify all report versions carrying a primary-suspect target drug.
    candidate_ids: set[str] = set()
    for record in records:
        usecols = ["primaryid", "role_cod", "drugname", "prod_ai"]
        for chunk in read_table_chunks(record.path, "DRUG", usecols):
            ps = chunk["role_cod"].str.strip().str.upper().eq("PS")
            chunk = chunk[ps]
            if chunk.empty:
                continue
            any_match = pd.Series(False, index=chunk.index)
            for aliases in DRUG_ALIASES.values():
                any_match |= drug_mask(chunk["prod_ai"], chunk["drugname"], aliases)
            candidate_ids.update(chunk.loc[any_match, "primaryid"].astype(str).tolist())

    # Pass 2: resolve corrected PRIMARYID rows and latest CASEID rows only among candidates.
    best_by_primary: dict[str, dict[str, object]] = {}
    usecols = ["primaryid", "caseid", "caseversion", "fda_dt", "age", "age_cod", "sex"]
    for record in records:
        for chunk in read_table_chunks(record.path, "DEMO", usecols):
            chunk = chunk[chunk["primaryid"].isin(candidate_ids)].copy()
            if chunk.empty:
                continue
            chunk["primaryid_num"] = pd.to_numeric(chunk["primaryid"], errors="coerce")
            chunk["fda_dt_num"] = pd.to_numeric(chunk["fda_dt"], errors="coerce")
            chunk["age_num"] = pd.to_numeric(chunk["age"].replace("", pd.NA), errors="coerce")
            code = chunk["age_cod"].str.strip().str.upper()
            factors = code.map({"YR": 1.0, "YEAR": 1.0, "YEARS": 1.0, "DEC": 10.0, "DECADE": 10.0, "DECADES": 10.0, "MON": 1 / 12, "MONTH": 1 / 12, "MONTHS": 1 / 12})
            chunk["age_years"] = chunk["age_num"] * factors
            if chunk[["primaryid_num", "fda_dt_num"]].isna().any().any():
                raise ValueError(f"Non-numeric primaryid/fda_dt in {record.period}")
            for row in chunk.to_dict("records"):
                key = str(row["primaryid"])
                normalized = {**row, "period": record.period}
                best_by_primary[key] = latest_row(best_by_primary.get(key), normalized)

    best_by_case: dict[str, dict[str, object]] = {}
    for row in best_by_primary.values():
        key = str(row["caseid"])
        best_by_case[key] = latest_row(best_by_case.get(key), row)
    eligible = {
        str(row["primaryid"]): row
        for row in best_by_case.values()
        if str(row["sex"]).strip().upper() == "F"
        and pd.notna(row["age_years"])
        and 15 <= float(row["age_years"]) <= 49
    }
    eligible_ids = set(eligible)
    if not eligible_ids:
        raise RuntimeError("No eligible reports after streaming latest-case resolution")

    # Pass 3: classify the mutually exclusive primary-suspect drug cohorts.
    exposure_ids = {drug: set() for drug in DRUG_ALIASES}
    for record in records:
        for chunk in read_table_chunks(record.path, "DRUG", ["primaryid", "role_cod", "drugname", "prod_ai"]):
            chunk = chunk[chunk["primaryid"].isin(eligible_ids)].copy()
            chunk = chunk[chunk["role_cod"].str.strip().str.upper().eq("PS")]
            if chunk.empty:
                continue
            for drug, aliases in DRUG_ALIASES.items():
                matched = drug_mask(chunk["prod_ai"], chunk["drugname"], aliases)
                exposure_ids[drug].update(chunk.loc[matched, "primaryid"].astype(str).tolist())

    riv = exposure_ids["rivaroxaban"]
    api = exposure_ids["apixaban"]
    riv_only, api_only = riv - api, api - riv
    if not riv_only or not api_only:
        raise RuntimeError("One active-comparator cohort is empty")

    # Pass 4: collect event IDs by term and by quarter.
    term_ids = {term: set() for group in TERMS.values() for term in group}
    for record in records:
        for chunk in read_table_chunks(record.path, "REAC", ["primaryid", "pt"]):
            chunk = chunk[chunk["primaryid"].isin(eligible_ids)].copy()
            if chunk.empty:
                continue
            pt = chunk["pt"].str.strip().str.upper()
            for term in term_ids:
                term_ids[term].update(chunk.loc[pt.eq(term), "primaryid"].astype(str).tolist())

    def cells(event_ids: set[str]) -> tuple[int, int, int, int]:
        return (len(riv_only & event_ids), len(riv_only - event_ids), len(api_only & event_ids), len(api_only - event_ids))

    rows = []
    for group, terms in TERMS.items():
        event_ids = set().union(*(term_ids[t] for t in terms))
        a, b, c, d = cells(event_ids)
        stats = two_by_two_statistics(a, b, c, d)
        rows.append({"outcome_group": group, "a_rivaroxaban_event": a, "b_rivaroxaban_non_event": b, "c_apixaban_event": c, "d_apixaban_non_event": d, **stats})

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(output_dir / "outcome_group_counts.csv", index=False)
    summary = {
        "periods": periods,
        "candidate_primary_ids": len(candidate_ids),
        "latest_candidate_cases": len(best_by_case),
        "eligible_female_15_49_reports": len(eligible_ids),
        "rivaroxaban_only_reports": len(riv_only),
        "apixaban_only_reports": len(api_only),
        "excluded_dual_exposure_reports": len(riv & api),
        "outcome_groups": rows,
    }
    (output_dir / "streaming_run_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    ap.add_argument("--periods", nargs="+", required=True)
    args = ap.parse_args()
    print(json.dumps(run(args.project_root.resolve(), args.manifest.resolve(), args.periods, args.output_dir.resolve()), indent=2))


if __name__ == "__main__":
    main()
