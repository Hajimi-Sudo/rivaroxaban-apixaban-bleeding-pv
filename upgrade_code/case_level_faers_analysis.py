"""Aggregate case-level FAERS descriptors for the current manuscript release.

This analysis reuses the current low-memory latest-case resolution and the
frozen female 15--49 primary-suspect rivaroxaban/apixaban comparator rule. It
does not export report identifiers or narratives. Outputs are aggregate only
and are intended to complement the disproportionality analysis, not to estimate
incidence, causality, or comparative treatment risk.
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd


def latest_row(old: dict[str, object] | None, row: dict[str, object]) -> dict[str, object]:
    if old is None:
        return row
    old_key = (int(old["fda_dt_num"]), int(old["primaryid_num"]), str(old["period"]))
    new_key = (int(row["fda_dt_num"]), int(row["primaryid_num"]), str(row["period"]))
    return row if new_key >= old_key else old


def normalize(values: pd.Series) -> pd.Series:
    return values.astype("string").fillna("").str.strip().str.upper()


def collect_latest_candidate_demo(records, candidate_ids: set[str], read_table_chunks):
    fields = [
        "primaryid", "caseid", "caseversion", "fda_dt", "event_dt", "age", "age_cod",
        "sex", "rept_cod", "occp_cod", "reporter_country", "occr_country",
    ]
    best_by_primary: dict[str, dict[str, object]] = {}
    for record in records:
        for chunk in read_table_chunks(record.path, "DEMO", fields):
            chunk = chunk[chunk["primaryid"].isin(candidate_ids)].copy()
            if chunk.empty:
                continue
            chunk["primaryid_num"] = pd.to_numeric(chunk["primaryid"], errors="coerce")
            chunk["fda_dt_num"] = pd.to_numeric(chunk["fda_dt"], errors="coerce")
            chunk["age_num"] = pd.to_numeric(chunk["age"].replace("", pd.NA), errors="coerce")
            age_code = normalize(chunk["age_cod"])
            factors = age_code.map({
                "YR": 1.0, "YEAR": 1.0, "YEARS": 1.0,
                "DEC": 10.0, "DECADE": 10.0, "DECADES": 10.0,
                "MON": 1 / 12, "MONTH": 1 / 12, "MONTHS": 1 / 12,
            })
            chunk["age_years"] = chunk["age_num"] * factors
            if chunk[["primaryid_num", "fda_dt_num"]].isna().any().any():
                raise ValueError(f"Non-numeric PRIMARYID/FDA_DT in {record.period}")
            for row in chunk.to_dict("records"):
                key = str(row["primaryid"])
                best_by_primary[key] = latest_row(
                    best_by_primary.get(key), {**row, "period": record.period}
                )

    best_by_case: dict[str, dict[str, object]] = {}
    for row in best_by_primary.values():
        key = str(row["caseid"])
        best_by_case[key] = latest_row(best_by_case.get(key), row)
    eligible = {
        str(row["primaryid"]): row
        for row in best_by_case.values()
        if normalize(pd.Series([row["sex"]])).iloc[0] == "F"
        and pd.notna(row["age_years"])
        and 15 <= float(row["age_years"]) <= 49
    }
    return eligible


def add_descriptor_rows(rows: list[dict[str, object]], subset: pd.DataFrame, drug: str, cohort: str, fields: list[str]) -> None:
    for field in fields:
        values = normalize(subset[field]).replace("", "MISSING")
        for value, count in values.value_counts(dropna=False).items():
            rows.append({
                "drug": drug,
                "cohort": cohort,
                "field": field,
                "value": str(value),
                "n": int(count),
                "fraction": float(count / len(subset)) if len(subset) else None,
                "reports_in_cohort": int(len(subset)),
            })


def run(project_root: Path, manifest_path: Path, periods: list[str], output_dir: Path) -> dict[str, object]:
    sys.path.insert(0, str(project_root / "experiment" / "code"))
    from faers_pilot import DRUG_ALIASES, drug_mask, load_manifest, read_table_chunks

    archive_dir = project_root / ".aris" / "data" / "faers" / "archives"
    records = load_manifest(manifest_path, archive_dir, periods)

    candidate_ids: set[str] = set()
    for record in records:
        for chunk in read_table_chunks(record.path, "DRUG", ["primaryid", "role_cod", "drugname", "prod_ai"]):
            ps = normalize(chunk["role_cod"]).eq("PS")
            chunk = chunk[ps]
            if chunk.empty:
                continue
            match = pd.Series(False, index=chunk.index)
            for drug in ("rivaroxaban", "apixaban"):
                match |= drug_mask(chunk["prod_ai"], chunk["drugname"], DRUG_ALIASES[drug])
            candidate_ids.update(chunk.loc[match, "primaryid"].astype(str).tolist())

    eligible = collect_latest_candidate_demo(records, candidate_ids, read_table_chunks)
    eligible_ids = set(eligible)
    if not eligible_ids:
        raise RuntimeError("No eligible female 15--49 reports were retained")

    exposure_ids = {"rivaroxaban": set(), "apixaban": set()}
    ps_keys = {"rivaroxaban": set(), "apixaban": set()}
    followup = {drug: {"dechal": {}, "rechal": {}} for drug in exposure_ids}
    for record in records:
        cols = ["primaryid", "drug_seq", "role_cod", "drugname", "prod_ai", "dechal", "rechal"]
        for chunk in read_table_chunks(record.path, "DRUG", cols):
            chunk = chunk[chunk["primaryid"].isin(eligible_ids)].copy()
            if chunk.empty:
                continue
            ps = normalize(chunk["role_cod"]).eq("PS")
            for drug in exposure_ids:
                matched = ps & drug_mask(chunk["prod_ai"], chunk["drugname"], DRUG_ALIASES[drug])
                ids = chunk.loc[matched, "primaryid"].astype(str)
                exposure_ids[drug].update(ids.tolist())
                ps_keys[drug].update(zip(ids.tolist(), chunk.loc[matched, "drug_seq"].tolist()))
                for field in ("dechal", "rechal"):
                    values = normalize(chunk.loc[matched, field])
                    for pid, code in zip(ids.tolist(), values.tolist()):
                        if code:
                            followup[drug][field].setdefault(pid, set()).add(code)

    riv_only = exposure_ids["rivaroxaban"] - exposure_ids["apixaban"]
    api_only = exposure_ids["apixaban"] - exposure_ids["rivaroxaban"]
    cohort_by_drug = {"rivaroxaban": riv_only, "apixaban": api_only}
    cohort_ids = riv_only | api_only

    aub_terms = {
        "ABNORMAL UTERINE BLEEDING", "HEAVY MENSTRUAL BLEEDING",
        "INTERMENSTRUAL BLEEDING", "UTERINE HAEMORRHAGE", "VAGINAL HAEMORRHAGE",
    }
    term_ids = {term: set() for term in aub_terms}
    for record in records:
        for chunk in read_table_chunks(record.path, "REAC", ["primaryid", "pt"]):
            chunk = chunk[chunk["primaryid"].isin(cohort_ids)].copy()
            if chunk.empty:
                continue
            pt = normalize(chunk["pt"])
            for term in aub_terms:
                term_ids[term].update(chunk.loc[pt.eq(term), "primaryid"].astype(str).tolist())
    broad_events = set().union(*term_ids.values())
    event_by_drug = {drug: ids & broad_events for drug, ids in cohort_by_drug.items()}

    demo = pd.DataFrame.from_dict(eligible, orient="index")
    demo.index = demo.index.astype(str)
    descriptor_fields = ["sex", "age_cod", "rept_cod", "occp_cod", "reporter_country", "occr_country"]
    descriptors: list[dict[str, object]] = []
    for drug, ids in cohort_by_drug.items():
        cohort_subset = demo.loc[sorted(ids)]
        event_subset = demo.loc[sorted(event_by_drug[drug])]
        add_descriptor_rows(descriptors, cohort_subset, drug, "all_active_comparator_reports", descriptor_fields)
        add_descriptor_rows(descriptors, event_subset, drug, "historical_broad_composite_events", descriptor_fields)
    pd.DataFrame(descriptors).to_csv(output_dir / "case_level_descriptors.csv", index=False)

    outcome_rows: list[dict[str, object]] = []
    event_union = set().union(*event_by_drug.values())
    outc_codes: dict[str, set[str]] = {}
    for record in records:
        for chunk in read_table_chunks(record.path, "OUTC", ["primaryid", "outc_cod"]):
            chunk = chunk[chunk["primaryid"].isin(event_union)].copy()
            if chunk.empty:
                continue
            codes = normalize(chunk["outc_cod"]).replace("", "MISSING")
            for code, ids in chunk.assign(code=codes).groupby("code")["primaryid"]:
                outc_codes.setdefault(str(code), set()).update(ids.astype(str).tolist())
    for drug, ids in event_by_drug.items():
        for code, coded_ids in sorted(outc_codes.items()):
            n = len(ids & coded_ids)
            outcome_rows.append({"drug": drug, "outcome_code": code, "n": n, "fraction_of_event_reports": n / len(ids) if ids else None, "event_reports": len(ids)})
    pd.DataFrame(outcome_rows).to_csv(output_dir / "case_level_serious_outcomes.csv", index=False)

    followup_rows: list[dict[str, object]] = []
    for drug, ids in event_by_drug.items():
        for field in ("dechal", "rechal"):
            any_recorded = set().union(*followup[drug][field].values()) if followup[drug][field] else set()
            recorded_ids = set(followup[drug][field]) & ids
            code_counts = {}
            for pid in recorded_ids:
                for code in followup[drug][field][pid]:
                    code_counts[code] = code_counts.get(code, 0) + 1
            for code, n in sorted(code_counts.items()):
                followup_rows.append({"drug": drug, "field": field, "code": code, "n": n, "fraction_of_event_reports": n / len(ids) if ids else None, "event_reports": len(ids)})
            followup_rows.append({"drug": drug, "field": field, "code": "ANY_RECORDED", "n": len(recorded_ids), "fraction_of_event_reports": len(recorded_ids) / len(ids) if ids else None, "event_reports": len(ids)})
    pd.DataFrame(followup_rows).to_csv(output_dir / "case_level_followup.csv", index=False)

    indication_rows: list[dict[str, object]] = []
    indication_patterns = {
        "venous_thromboembolism": re.compile(r"DEEP VEIN THROMBOSIS|PULMONARY EMBOLISM|VENOUS THROMBO|\bVTE\b", re.I),
        "atrial_fibrillation": re.compile(r"ATRIAL FIBRILLATION|ATRIAL FLUTTER", re.I),
    }
    all_keys = set().union(*ps_keys.values())
    indication_ids = {drug: {name: set() for name in indication_patterns} for drug in cohort_by_drug}
    for record in records:
        for chunk in read_table_chunks(record.path, "INDI", ["primaryid", "indi_drug_seq", "indi_pt"]):
            keys = list(zip(chunk["primaryid"].astype(str), chunk["indi_drug_seq"].tolist()))
            keep = [key in all_keys for key in keys]
            chunk = chunk.loc[keep].copy()
            if chunk.empty:
                continue
            text = chunk["indi_pt"].astype("string").fillna("")
            for drug, drug_keys in ps_keys.items():
                dmask = [key in drug_keys for key in zip(chunk["primaryid"].astype(str), chunk["indi_drug_seq"].tolist())]
                for name, pattern in indication_patterns.items():
                    indication_ids[drug][name].update(chunk.loc[[a and bool(pattern.search(str(b))) for a, b in zip(dmask, text.tolist())], "primaryid"].astype(str).tolist())
    for drug, ids in cohort_by_drug.items():
        for name, matched in indication_ids[drug].items():
            n_all = len(ids & matched)
            n_event = len(event_by_drug[drug] & matched)
            indication_rows.append({"drug": drug, "indication_group": name, "all_cohort_reports": n_all, "all_cohort_fraction": n_all / len(ids) if ids else None, "event_reports": n_event, "event_fraction": n_event / len(event_by_drug[drug]) if event_by_drug[drug] else None})
    pd.DataFrame(indication_rows).to_csv(output_dir / "case_level_indications.csv", index=False)

    summary = {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "periods": periods,
        "candidate_primary_ids": len(candidate_ids),
        "eligible_female_15_49_reports": len(eligible_ids),
        "rivaroxaban_only_reports": len(riv_only),
        "apixaban_only_reports": len(api_only),
        "excluded_dual_target_reports": len(exposure_ids["rivaroxaban"] & exposure_ids["apixaban"]),
        "historical_broad_composite_event_reports": {drug: len(ids) for drug, ids in event_by_drug.items()},
        "outcome_term_event_reports": {term: {drug: len(ids & cohort_by_drug[drug]) for drug in cohort_by_drug} for term, ids in term_ids.items()},
        "outputs_are_aggregate_only": True,
        "interpretation": "case-level descriptive analysis of report characteristics; no incidence, causality, or comparative treatment risk",
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "case_level_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    return summary


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--periods", nargs="+", required=True)
    args = parser.parse_args()
    args.output_dir.mkdir(parents=True, exist_ok=True)
    print(json.dumps(run(args.project_root.resolve(), args.manifest.resolve(), args.periods, args.output_dir.resolve()), indent=2))


if __name__ == "__main__":
    main()
