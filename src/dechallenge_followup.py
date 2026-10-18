"""Targeted FAERS dechallenge/rechallenge coded-value follow-up.

This script preserves the frozen cohort, exposure, and broad outcome rules but
parallelizes independent quarterly DRUG/REAC scans. It reports non-exclusive
report-level code distributions only; it performs no causality assessment.
"""

from __future__ import annotations

import argparse
import json
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from faers_pilot import (
    AUB_TERMS,
    DRUG_ALIASES,
    drug_mask,
    load_and_deduplicate_demo,
    load_manifest,
    read_table_chunks,
)
from faers_full import period_range_from_manifest


TARGET_DRUGS = ("rivaroxaban", "apixaban")
FOLLOWUP_FIELDS = ("dechal", "rechal")
_WORKER_UNIVERSE_IDS: set[str] = set()


def initialize_worker(universe_ids: set[str]) -> None:
    """Copy the frozen cohort once per worker, not once per submitted quarter."""
    global _WORKER_UNIVERSE_IDS
    _WORKER_UNIVERSE_IDS = universe_ids


def scan_record_worker(record) -> dict[str, object]:
    return scan_record(record, _WORKER_UNIVERSE_IDS)


def normalize_code(series: pd.Series) -> pd.Series:
    return series.astype("string").fillna("").str.strip().str.upper()


def scan_record(record, universe_ids: set[str]) -> dict[str, object]:
    exposures = {drug: set() for drug in TARGET_DRUGS}
    codes = {
        drug: {field: {} for field in FOLLOWUP_FIELDS}
        for drug in TARGET_DRUGS
    }
    outcomes: set[str] = set()

    drug_cols = [
        "primaryid", "caseid", "drug_seq", "role_cod", "drugname", "prod_ai",
        "dechal", "rechal",
    ]
    for chunk in read_table_chunks(record.path, "DRUG", drug_cols):
        chunk = chunk[chunk["primaryid"].isin(universe_ids)].copy()
        if chunk.empty:
            continue
        ps = normalize_code(chunk["role_cod"]).eq("PS")
        for drug in TARGET_DRUGS:
            matched = ps & drug_mask(chunk["prod_ai"], chunk["drugname"], DRUG_ALIASES[drug])
            exposures[drug].update(chunk.loc[matched, "primaryid"].tolist())
            for field in FOLLOWUP_FIELDS:
                normalized = normalize_code(chunk[field])
                for code, ids in chunk.loc[matched & normalized.ne("")].assign(
                    normalized_code=normalized[matched & normalized.ne("")]
                ).groupby("normalized_code")["primaryid"]:
                    codes[drug][field].setdefault(str(code), set()).update(ids.tolist())

    reac_cols = ["primaryid", "caseid", "pt"]
    for chunk in read_table_chunks(record.path, "REAC", reac_cols):
        chunk = chunk[chunk["primaryid"].isin(universe_ids)].copy()
        pt = normalize_code(chunk["pt"])
        outcomes.update(chunk.loc[pt.isin(AUB_TERMS), "primaryid"].tolist())

    return {"exposures": exposures, "codes": codes, "outcomes": outcomes}


def merge_scan_results(results: list[dict[str, object]]) -> tuple[dict, dict, set[str]]:
    exposures = {drug: set() for drug in TARGET_DRUGS}
    codes = {drug: {field: {} for field in FOLLOWUP_FIELDS} for drug in TARGET_DRUGS}
    outcomes: set[str] = set()
    for result in results:
        outcomes.update(result["outcomes"])
        for drug in TARGET_DRUGS:
            exposures[drug].update(result["exposures"][drug])
            for field in FOLLOWUP_FIELDS:
                for code, ids in result["codes"][drug][field].items():
                    codes[drug][field].setdefault(code, set()).update(ids)
    return exposures, codes, outcomes


def summarize_code_sets(
    drug: str, field: str, event_ids: set[str], code_sets: dict[str, set[str]]
) -> tuple[list[dict[str, object]], dict[str, object]]:
    rows = []
    recorded_union: set[str] = set()
    membership_count: dict[str, int] = {}
    for code, ids in sorted(code_sets.items()):
        matched = event_ids & ids
        recorded_union.update(matched)
        for pid in matched:
            membership_count[pid] = membership_count.get(pid, 0) + 1
        rows.append(
            {
                "drug": drug,
                "field": field,
                "code": code,
                "event_reports": len(event_ids),
                "reports_with_code": len(matched),
                "fraction_of_event_reports": len(matched) / len(event_ids) if event_ids else None,
                "codes_nonexclusive": True,
            }
        )
    audit = {
        "drug": drug,
        "field": field,
        "event_reports": len(event_ids),
        "reports_with_any_recorded_code": len(recorded_union),
        "recorded_fraction": len(recorded_union) / len(event_ids) if event_ids else None,
        "reports_with_conflicting_multiple_codes": sum(value > 1 for value in membership_count.values()),
    }
    return rows, audit


def run(project_root: Path, output_dir: Path, workers: int = 8) -> dict[str, object]:
    archives_dir = project_root / "data" / "raw" / "faers" / "archives"
    manifest_path = project_root / "data" / "raw" / "faers" / "manifest" / "download_manifest_sha256.csv"
    periods = period_range_from_manifest(manifest_path)
    records = load_manifest(manifest_path, archives_dir, periods)
    output_dir.mkdir(parents=True, exist_ok=True)

    demo, _ = load_and_deduplicate_demo(records)
    female = demo[demo["sex"].str.strip().str.upper().eq("F") & demo["age_years"].notna()]
    cohort = set(female.loc[female["age_years"].between(15, 49), "primaryid"])

    with ProcessPoolExecutor(
        max_workers=workers,
        initializer=initialize_worker,
        initargs=(cohort,),
    ) as pool:
        results = list(pool.map(scan_record_worker, records))
    exposures, codes, broad_aub = merge_scan_results(results)

    riv_only = (exposures["rivaroxaban"] & cohort) - exposures["apixaban"]
    api_only = (exposures["apixaban"] & cohort) - exposures["rivaroxaban"]
    event_ids = {
        "rivaroxaban": riv_only & broad_aub,
        "apixaban": api_only & broad_aub,
    }
    if len(event_ids["rivaroxaban"]) != 665 or len(event_ids["apixaban"]) != 90:
        raise RuntimeError(f"Frozen event-cell mismatch: { {k: len(v) for k, v in event_ids.items()} }")

    rows = []
    audits = []
    for drug in TARGET_DRUGS:
        for field in FOLLOWUP_FIELDS:
            field_rows, audit = summarize_code_sets(drug, field, event_ids[drug], codes[drug][field])
            rows.extend(field_rows)
            audits.append(audit)
    pd.DataFrame(rows).to_csv(output_dir / "dechallenge_rechallenge_code_distribution.csv", index=False)
    pd.DataFrame(audits).to_csv(output_dir / "dechallenge_rechallenge_code_audit.csv", index=False)

    metadata = {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "periods": periods,
        "workers": workers,
        "frozen_primary_event_cells": {key: len(value) for key, value in event_ids.items()},
        "interpretation": "descriptive non-exclusive report-level code distributions; no causality assessment",
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return metadata


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=8)
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    print(json.dumps(run(args.project_root.resolve(), args.output_dir, args.workers), indent=2))
