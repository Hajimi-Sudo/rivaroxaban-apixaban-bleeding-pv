"""Reproducible FAERS pilot for oral-anticoagulant-associated AUB.

Reads official quarterly ASCII ZIP archives without modifying or extracting
the raw files. The pilot covers 2025Q1-Q4 and implements cross-quarter case
deduplication before cohort, exposure, and outcome filtering.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import re
import sys
import zipfile
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


TABLE_PATTERNS = {
    "DEMO": re.compile(r"(?i)(?:^|/)DEMO\d{2}Q[1-4](?:_new)?\.txt$"),
    "DRUG": re.compile(r"(?i)(?:^|/)DRUG\d{2}Q[1-4]\.txt$"),
    "REAC": re.compile(r"(?i)(?:^|/)REAC\d{2}Q[1-4]\.txt$"),
    "INDI": re.compile(r"(?i)(?:^|/)INDI\d{2}Q[1-4]\.txt$"),
    "THER": re.compile(r"(?i)(?:^|/)THER\d{2}Q[1-4]\.txt$"),
    "OUTC": re.compile(r"(?i)(?:^|/)OUTC\d{2}Q[1-4]\.txt$"),
}

REQUIRED_COLUMNS = {
    "DEMO": {"primaryid", "caseid", "caseversion", "fda_dt", "age", "age_cod", "sex"},
    "DRUG": {"primaryid", "caseid", "drug_seq", "role_cod", "drugname", "prod_ai"},
    "REAC": {"primaryid", "caseid", "pt"},
    "INDI": {"primaryid", "caseid", "indi_drug_seq", "indi_pt"},
    "THER": {"primaryid", "caseid", "dsg_drug_seq", "start_dt", "end_dt"},
}

AUB_TERMS = {
    "ABNORMAL UTERINE BLEEDING",
    "HEAVY MENSTRUAL BLEEDING",
    "INTERMENSTRUAL BLEEDING",
    "UTERINE HAEMORRHAGE",
    "VAGINAL HAEMORRHAGE",
}

DRUG_ALIASES = {
    "rivaroxaban": ("RIVAROXABAN", "XARELTO"),
    "apixaban": ("APIXABAN", "ELIQUIS"),
    "edoxaban": ("EDOXABAN", "SAVAYSA", "LIXIANA"),
    "dabigatran": ("DABIGATRAN", "PRADAXA"),
    "warfarin": ("WARFARIN", "COUMADIN", "JANTOVEN"),
}

AGE_TO_YEARS = {
    "YR": 1.0,
    "YEAR": 1.0,
    "YEARS": 1.0,
    "MON": 1.0 / 12.0,
    "MONTH": 1.0 / 12.0,
    "MONTHS": 1.0 / 12.0,
    "WK": 7.0 / 365.25,
    "WEEK": 7.0 / 365.25,
    "WEEKS": 7.0 / 365.25,
    "DY": 1.0 / 365.25,
    "DAY": 1.0 / 365.25,
    "DAYS": 1.0 / 365.25,
    "DEC": 10.0,
    "DECADE": 10.0,
    "DECADES": 10.0,
    "HR": 1.0 / (365.25 * 24.0),
    "HOUR": 1.0 / (365.25 * 24.0),
    "HOURS": 1.0 / (365.25 * 24.0),
}


@dataclass(frozen=True)
class QuarterArchive:
    period: str
    path: Path
    expected_sha256: str


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def find_member(archive: zipfile.ZipFile, table: str) -> str:
    matches = [name for name in archive.namelist() if TABLE_PATTERNS[table].search(name)]
    if len(matches) != 1:
        raise ValueError(f"Expected one {table} TXT member, found {matches}")
    return matches[0]


def header_columns(archive_path: Path, table: str) -> list[str]:
    with zipfile.ZipFile(archive_path) as archive:
        member = find_member(archive, table)
        with archive.open(member) as handle:
            header = handle.readline().decode("latin-1").strip("\r\n")
    return [column.strip().lower() for column in header.split("$")]


def read_table_chunks(
    archive_path: Path,
    table: str,
    usecols: list[str],
    chunksize: int = 250_000,
):
    with zipfile.ZipFile(archive_path) as archive:
        member = find_member(archive, table)
        with archive.open(member) as handle:
            yield from pd.read_csv(
                handle,
                sep="$",
                usecols=usecols,
                dtype="string",
                encoding="latin-1",
                chunksize=chunksize,
                keep_default_na=False,
            )


def age_in_years(age: pd.Series, age_code: pd.Series) -> pd.Series:
    numeric_age = pd.to_numeric(age.replace("", pd.NA), errors="coerce")
    normalized_code = age_code.astype("string").str.strip().str.upper()
    factors = normalized_code.map(AGE_TO_YEARS)
    result = numeric_age * factors
    return result.where((result >= 0) & (result <= 130))


def normalize_text(series: pd.Series) -> pd.Series:
    return (
        series.astype("string")
        .str.upper()
        .str.replace(r"[^A-Z0-9]+", " ", regex=True)
        .str.replace(r"\s+", " ", regex=True)
        .str.strip()
    )


def drug_mask(prod_ai: pd.Series, drugname: pd.Series, aliases: tuple[str, ...]) -> pd.Series:
    combined = normalize_text(prod_ai.fillna("") + " " + drugname.fillna(""))
    pattern = r"(?:^|\s)(?:" + "|".join(re.escape(alias) for alias in aliases) + r")(?:\s|$)"
    return combined.str.contains(pattern, regex=True, na=False)


def reporting_odds_ratio(a: int, b: int, c: int, d: int) -> dict[str, float | int | bool]:
    cells = np.array([a, b, c, d], dtype=float)
    corrected = bool(np.any(cells == 0))
    if corrected:
        cells += 0.5
    aa, bb, cc, dd = cells
    ror = (aa * dd) / (bb * cc)
    se = math.sqrt(1.0 / aa + 1.0 / bb + 1.0 / cc + 1.0 / dd)
    log_ror = math.log(ror)
    return {
        "a_rivaroxaban_aub": a,
        "b_rivaroxaban_non_aub": b,
        "c_apixaban_aub": c,
        "d_apixaban_non_aub": d,
        "ror": float(ror),
        "ci95_lower": float(math.exp(log_ror - 1.96 * se)),
        "ci95_upper": float(math.exp(log_ror + 1.96 * se)),
        "jeffreys_0_5_correction": corrected,
    }


def load_manifest(
    manifest_path: Path, archives_dir: Path, periods: list[str]
) -> list[QuarterArchive]:
    manifest = pd.read_csv(manifest_path, dtype="string")
    selected = manifest[manifest["Period"].isin(periods)].copy()
    if len(selected) != len(periods):
        found = set(selected["Period"].tolist())
        raise ValueError(f"Missing requested periods: {sorted(set(periods) - found)}")
    records = []
    for row in selected.sort_values("Period").itertuples(index=False):
        records.append(
            QuarterArchive(
                period=str(row.Period),
                path=archives_dir / str(row.FileName),
                expected_sha256=str(row.SHA256).upper(),
            )
        )
    return records


def schema_audit(records: list[QuarterArchive]) -> dict[str, object]:
    audit: dict[str, object] = {
        "pilot_periods": [record.period for record in records],
        "archives": {},
        "required_variable_coverage": None,
    }
    required_total = 0
    present_total = 0
    for record in records:
        if not record.path.exists():
            raise FileNotFoundError(record.path)
        actual_hash = sha256(record.path)
        if actual_hash != record.expected_sha256:
            raise ValueError(f"SHA-256 mismatch for {record.path.name}")
        with zipfile.ZipFile(record.path) as archive:
            bad_member = archive.testzip()
        if bad_member is not None:
            raise zipfile.BadZipFile(f"CRC failure in {record.path.name}: {bad_member}")

        table_results: dict[str, object] = {}
        for table, required in REQUIRED_COLUMNS.items():
            columns = header_columns(record.path, table)
            missing = sorted(required - set(columns))
            required_total += len(required)
            present_total += len(required) - len(missing)
            table_results[table] = {
                "columns": columns,
                "missing_required": missing,
                "readable": not missing,
            }
        audit["archives"][record.period] = {
            "file": record.path.name,
            "bytes": record.path.stat().st_size,
            "sha256": actual_hash,
            "tables": table_results,
        }
    audit["required_variable_coverage"] = present_total / required_total
    return audit


def deduplicate_demo_frame(demo: pd.DataFrame) -> pd.DataFrame:
    """Apply PRIMARYID correction handling, then the prespecified CASEID rule."""
    demo["primaryid_num"] = pd.to_numeric(demo["primaryid"], errors="coerce")
    demo["caseversion_num"] = pd.to_numeric(demo["caseversion"], errors="coerce")
    demo["fda_dt_num"] = pd.to_numeric(demo["fda_dt"], errors="coerce")
    if demo[["primaryid_num", "caseversion_num", "fda_dt_num"]].isna().any().any():
        raise ValueError("Non-numeric primaryid, caseversion, or fda_dt encountered in DEMO")

    # The same report version can recur in adjacent extracts. CASEID and
    # CASEVERSION must agree. FDA_DT can be corrected in a later official
    # extract without a new PRIMARYID; retain the greatest FDA_DT and audit it.
    duplicated = demo[demo.duplicated("primaryid", keep=False)]
    fda_dt_conflicts = pd.DataFrame()
    if not duplicated.empty:
        consistency = duplicated.groupby("primaryid").agg(
            caseid_n=("caseid", "nunique"),
            caseversion_n=("caseversion_num", "nunique"),
            fda_dt_n=("fda_dt_num", "nunique"),
        )
        structural_conflict = (consistency[["caseid_n", "caseversion_n"]] > 1).any(axis=1)
        if structural_conflict.any():
            bad_ids = consistency[structural_conflict].index.tolist()[:10]
            raise ValueError(f"Conflicting duplicate primaryid records: {bad_ids}")
        conflict_ids = set(consistency.index[consistency["fda_dt_n"] > 1])
        if conflict_ids:
            fda_dt_conflicts = duplicated[duplicated["primaryid"].isin(conflict_ids)][
                ["primaryid", "caseid", "caseversion", "fda_dt", "period"]
            ].sort_values(["primaryid", "fda_dt", "period"])
        demo = demo.sort_values(["primaryid", "fda_dt_num", "period"], kind="mergesort")
        demo = demo.drop_duplicates("primaryid", keep="last")
    demo = demo.sort_values(["caseid", "fda_dt_num", "primaryid_num"], kind="mergesort")
    dedup = demo.drop_duplicates("caseid", keep="last").copy()
    dedup["age_years"] = age_in_years(dedup["age"], dedup["age_cod"])
    dedup.attrs["primaryid_fda_dt_conflicts"] = fda_dt_conflicts.to_dict("records")
    return dedup


def load_and_deduplicate_demo(records: list[QuarterArchive]) -> tuple[pd.DataFrame, list[dict[str, object]]]:
    frames = []
    raw_by_quarter = []
    usecols = [
        "primaryid", "caseid", "caseversion", "fda_dt", "event_dt",
        "age", "age_cod", "age_grp", "sex", "rept_cod", "occp_cod",
        "reporter_country", "occr_country",
    ]
    for record in records:
        quarter_frames = list(read_table_chunks(record.path, "DEMO", usecols))
        quarter = pd.concat(quarter_frames, ignore_index=True)
        quarter["period"] = record.period
        frames.append(quarter)
        raw_by_quarter.append({"period": record.period, "raw_demo_rows": len(quarter)})
    demo = deduplicate_demo_frame(pd.concat(frames, ignore_index=True))
    return demo, raw_by_quarter


def collect_exposures(records: list[QuarterArchive], eligible_ids: set[str]) -> dict[str, set[str]]:
    exposures = {drug: set() for drug in DRUG_ALIASES}
    usecols = ["primaryid", "caseid", "drug_seq", "role_cod", "drugname", "prod_ai"]
    for record in records:
        for chunk in read_table_chunks(record.path, "DRUG", usecols):
            chunk = chunk[chunk["primaryid"].isin(eligible_ids)]
            chunk = chunk[chunk["role_cod"].str.strip().str.upper().eq("PS")]
            if chunk.empty:
                continue
            for drug, aliases in DRUG_ALIASES.items():
                matched = drug_mask(chunk["prod_ai"], chunk["drugname"], aliases)
                exposures[drug].update(chunk.loc[matched, "primaryid"].tolist())
    return exposures


def collect_aub_cases(records: list[QuarterArchive], eligible_ids: set[str]) -> tuple[set[str], pd.DataFrame]:
    cases: set[str] = set()
    term_rows = []
    usecols = ["primaryid", "caseid", "pt"]
    for record in records:
        quarter_counts = {term: 0 for term in sorted(AUB_TERMS)}
        for chunk in read_table_chunks(record.path, "REAC", usecols):
            chunk = chunk[chunk["primaryid"].isin(eligible_ids)].copy()
            normalized_pt = chunk["pt"].str.strip().str.upper()
            matched = chunk[normalized_pt.isin(AUB_TERMS)].copy()
            matched["pt_normalized"] = normalized_pt[normalized_pt.isin(AUB_TERMS)]
            cases.update(matched["primaryid"].tolist())
            for term, count in matched["pt_normalized"].value_counts().items():
                quarter_counts[str(term)] += int(count)
        term_rows.extend(
            {"period": record.period, "pt": term, "reaction_rows": count}
            for term, count in quarter_counts.items()
        )
    return cases, pd.DataFrame(term_rows)


def run(project_root: Path, periods: list[str], output_dir: Path) -> dict[str, object]:
    archives_dir = project_root / "data" / "raw" / "faers" / "archives"
    manifest_path = project_root / "data" / "raw" / "faers" / "manifest" / "download_manifest_sha256.csv"
    if not output_dir.is_absolute():
        output_dir = project_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    records = load_manifest(manifest_path, archives_dir, periods)
    audit = schema_audit(records)
    (output_dir / "schema_audit.json").write_text(
        json.dumps(audit, indent=2, ensure_ascii=False), encoding="utf-8"
    )
    if audit["required_variable_coverage"] != 1.0:
        raise RuntimeError("Required-variable coverage is below 100%")

    demo, raw_by_quarter = load_and_deduplicate_demo(records)
    female_demo = demo[demo["sex"].str.strip().str.upper().eq("F")].copy()
    female_convertible = female_demo[female_demo["age_years"].notna()].copy()
    eligible = female_convertible[
        female_convertible["age_years"].between(15, 49, inclusive="both")
    ].copy()
    eligible_ids = set(eligible["primaryid"].tolist())

    period_label = f"{periods[0]}-{periods[-1]}"
    flow_rows = raw_by_quarter + [
        {
            "period": period_label,
            "stage": "raw_demo_rows",
            "n": sum(int(row["raw_demo_rows"]) for row in raw_by_quarter),
        },
        {"period": period_label, "stage": "deduplicated_cases", "n": len(demo)},
        {"period": period_label, "stage": "female_deduplicated_cases", "n": len(female_demo)},
        {"period": period_label, "stage": "female_convertible_age_cases", "n": len(female_convertible)},
        {"period": period_label, "stage": "female_age_15_49_cases", "n": len(eligible)},
    ]
    # Normalize the quarter rows to the same schema as the aggregate flow rows.
    normalized_flow = []
    for row in flow_rows:
        if "raw_demo_rows" in row:
            normalized_flow.append({"period": row["period"], "stage": "raw_demo_rows", "n": row["raw_demo_rows"]})
        else:
            normalized_flow.append(row)
    pd.DataFrame(normalized_flow).to_csv(output_dir / "cohort_flow.csv", index=False)

    exposures = collect_exposures(records, eligible_ids)
    if any(len(ids) == 0 for ids in exposures.values()):
        missing = [drug for drug, ids in exposures.items() if not ids]
        raise RuntimeError(f"Zero eligible reports after normalization for: {missing}")
    aub_cases, term_counts = collect_aub_cases(records, eligible_ids)
    term_counts.to_csv(output_dir / "aub_term_counts_by_quarter.csv", index=False)

    count_rows = []
    for drug, ids in exposures.items():
        event_ids = ids & aub_cases
        count_rows.append(
            {
                "drug": drug,
                "eligible_primary_suspect_reports": len(ids),
                "aub_reports": len(event_ids),
                "non_aub_reports": len(ids - aub_cases),
                "stable_20_event_flag": len(event_ids) >= 20,
            }
        )
    counts = pd.DataFrame(count_rows).sort_values("drug")
    counts.to_csv(output_dir / "drug_event_counts.csv", index=False)

    riv_only = exposures["rivaroxaban"] - exposures["apixaban"]
    api_only = exposures["apixaban"] - exposures["rivaroxaban"]
    overlap = exposures["rivaroxaban"] & exposures["apixaban"]
    a = len(riv_only & aub_cases)
    b = len(riv_only - aub_cases)
    c = len(api_only & aub_cases)
    d = len(api_only - aub_cases)
    contrast = reporting_odds_ratio(a, b, c, d)
    reconciliation = {
        "rivaroxaban_cells_equal_exclusive_reports": a + b == len(riv_only),
        "apixaban_cells_equal_exclusive_reports": c + d == len(api_only),
        "four_cells_equal_pairwise_union": a + b + c + d == len(riv_only | api_only),
        "event_cells_disjoint": (riv_only & aub_cases).isdisjoint(api_only & aub_cases),
        "non_event_cells_disjoint": (riv_only - aub_cases).isdisjoint(api_only - aub_cases),
    }
    if not all(reconciliation.values()):
        raise RuntimeError(f"2x2 reconciliation failure: {reconciliation}")
    pd.DataFrame(
        [
            {"exposure": "rivaroxaban_only", "outcome": "AUB", "n": a},
            {"exposure": "rivaroxaban_only", "outcome": "non_AUB", "n": b},
            {"exposure": "apixaban_only", "outcome": "AUB", "n": c},
            {"exposure": "apixaban_only", "outcome": "non_AUB", "n": d},
        ]
    ).to_csv(output_dir / "pairwise_2x2_cells.csv", index=False)
    contrast.update(
        {
            "pairwise_rule": "mutually exclusive rivaroxaban-only versus apixaban-only reports",
            "excluded_dual_exposure_reports": len(overlap),
            "pilot_period": period_label,
            "reconciliation": reconciliation,
            "thresholds": {
                "rivaroxaban_aub_cases_gte_20": bool(a >= 20),
                "apixaban_aub_cases_gte_1": bool(c >= 1),
                "ror_gt_1": bool(contrast["ror"] > 1.0),
            },
        }
    )
    contrast["pilot_gate_pass"] = bool(all(contrast["thresholds"].values()))
    (output_dir / "primary_contrast.json").write_text(
        json.dumps(contrast, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    run_metadata = {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "input_periods": [record.period for record in records],
        "input_files": [record.path.name for record in records],
        "aub_terms": sorted(AUB_TERMS),
        "drug_aliases": {drug: list(aliases) for drug, aliases in DRUG_ALIASES.items()},
    }
    (output_dir / "run_metadata.json").write_text(
        json.dumps(run_metadata, indent=2, ensure_ascii=False), encoding="utf-8"
    )

    riv_row = counts.loc[counts["drug"].eq("rivaroxaban")].iloc[0]
    api_row = counts.loc[counts["drug"].eq("apixaban")].iloc[0]
    summary = f"""# FAERS {period_label} Pilot Summary

- Required-variable coverage: {audit['required_variable_coverage']:.1%}
- Deduplicated cases: {len(demo):,}
- Eligible female reports aged 15–49 years: {len(eligible):,}
- Rivaroxaban primary-suspect reports / AUB reports: {int(riv_row['eligible_primary_suspect_reports']):,} / {int(riv_row['aub_reports']):,}
- Apixaban primary-suspect reports / AUB reports: {int(api_row['eligible_primary_suspect_reports']):,} / {int(api_row['aub_reports']):,}
- Exclusive active-comparator ROR: {contrast['ror']:.3f} (95% CI {contrast['ci95_lower']:.3f}–{contrast['ci95_upper']:.3f})
- Dual rivaroxaban/apixaban reports excluded from the pairwise contrast: {len(overlap):,}
- Pilot gate: {'PASS' if contrast['pilot_gate_pass'] else 'FAIL'}

This is a data-and-direction gate, not the final confirmatory estimate. The full analysis must repeat cross-quarter deduplication across all 45 quarters.
"""
    (output_dir / "pilot_summary.md").write_text(summary, encoding="utf-8")
    return {"output_dir": str(output_dir), "contrast": contrast, "eligible_cases": len(eligible)}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--project-root",
        type=Path,
        default=Path(__file__).resolve().parents[1],
    )
    parser.add_argument(
        "--periods",
        nargs="+",
        default=[f"2025Q{quarter}" for quarter in range(1, 5)],
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("results/faers_pilot"),
    )
    return parser.parse_args()


if __name__ == "__main__":
    arguments = parse_args()
    result = run(arguments.project_root.resolve(), arguments.periods, arguments.output_dir)
    print(json.dumps(result, indent=2, ensure_ascii=False))
