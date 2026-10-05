"""Full prespecified FAERS analysis for oral-anticoagulant-associated AUB."""

from __future__ import annotations

import argparse
import json
import platform
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from faers_pilot import (
    AUB_TERMS,
    DRUG_ALIASES,
    drug_mask,
    load_and_deduplicate_demo,
    load_manifest,
    read_table_chunks,
    schema_audit,
)
from pv_stats import (
    benjamini_hochberg,
    log_or_heterogeneity,
    minimum_detectable_or,
    two_by_two_statistics,
)


NARROW_AUB_TERMS = {"ABNORMAL UTERINE BLEEDING", "HEAVY MENSTRUAL BLEEDING"}
EXPLANATORY_OUTCOME_GROUPS = {
    "menstrual_specific": {"HEAVY MENSTRUAL BLEEDING", "INTERMENSTRUAL BLEEDING"},
    "uterine_specific": {"ABNORMAL UTERINE BLEEDING", "UTERINE HAEMORRHAGE"},
    "broad_excluding_vaginal_haemorrhage": set(AUB_TERMS) - {"VAGINAL HAEMORRHAGE"},
}
ALL_SUSPECT_ROLES = {"PS", "SS"}
VTE_PATTERN = re.compile(r"\b(?:DEEP VEIN THROMBOSIS|PULMONARY EMBOLISM|VENOUS THROMBO|VTE)\b", re.I)
AF_PATTERN = re.compile(r"\b(?:ATRIAL FIBRILLATION|ATRIAL FLUTTER)\b", re.I)
# FDA AEMS identified menorrhagia as a potential class-wide oral-anticoagulant
# signal in 2021Q1; the relevant labels were updated in April 2021.
SAFETY_COMMUNICATION_PERIODS = {"2021Q1", "2021Q2"}
HEALTHCARE_OCCUPATIONS = {"MD", "PH", "HP", "OT"}
CONSUMER_OCCUPATIONS = {"CN"}
LAWYER_OCCUPATIONS = {"LW"}
MAX_PLAUSIBLE_TTO_DAYS = 3650


def collect_drugs(
    records, universe_ids: set[str]
) -> tuple[
    dict[str, dict[str, set[str]]],
    dict[str, set[tuple[str, str]]],
    dict[str, dict[str, set[str]]],
]:
    exposures = {
        "ps": {drug: set() for drug in DRUG_ALIASES},
        "all_suspect": {drug: set() for drug in DRUG_ALIASES},
    }
    ps_keys = {drug: set() for drug in DRUG_ALIASES}
    followup = {
        "dechallenge_recorded": {drug: set() for drug in DRUG_ALIASES},
        "rechallenge_recorded": {drug: set() for drug in DRUG_ALIASES},
    }
    usecols = [
        "primaryid", "caseid", "drug_seq", "role_cod", "drugname", "prod_ai",
        "dechal", "rechal",
    ]
    for record in records:
        for chunk in read_table_chunks(record.path, "DRUG", usecols):
            chunk = chunk[chunk["primaryid"].isin(universe_ids)].copy()
            if chunk.empty:
                continue
            roles = chunk["role_cod"].str.strip().str.upper()
            for drug, aliases in DRUG_ALIASES.items():
                matched = drug_mask(chunk["prod_ai"], chunk["drugname"], aliases)
                ps = matched & roles.eq("PS")
                suspect = matched & roles.isin(ALL_SUSPECT_ROLES)
                exposures["ps"][drug].update(chunk.loc[ps, "primaryid"].tolist())
                exposures["all_suspect"][drug].update(chunk.loc[suspect, "primaryid"].tolist())
                ps_keys[drug].update(
                    zip(chunk.loc[ps, "primaryid"].tolist(), chunk.loc[ps, "drug_seq"].tolist())
                )
                followup["dechallenge_recorded"][drug].update(
                    chunk.loc[ps & chunk["dechal"].fillna("").str.strip().ne(""), "primaryid"].tolist()
                )
                followup["rechallenge_recorded"][drug].update(
                    chunk.loc[ps & chunk["rechal"].fillna("").str.strip().ne(""), "primaryid"].tolist()
                )
    return exposures, ps_keys, followup


def collect_outcomes(
    records, universe_ids: set[str]
) -> tuple[set[str], set[str], dict[str, set[str]], pd.DataFrame]:
    broad: set[str] = set()
    narrow: set[str] = set()
    reports_by_term = {term: set() for term in sorted(AUB_TERMS)}
    rows = []
    usecols = ["primaryid", "caseid", "pt"]
    for record in records:
        quarter_seen = {term: set() for term in sorted(AUB_TERMS)}
        for chunk in read_table_chunks(record.path, "REAC", usecols):
            chunk = chunk[chunk["primaryid"].isin(universe_ids)].copy()
            pt = chunk["pt"].str.strip().str.upper()
            for term in AUB_TERMS:
                ids = set(chunk.loc[pt.eq(term), "primaryid"].tolist())
                quarter_seen[term].update(ids)
                reports_by_term[term].update(ids)
                broad.update(ids)
                if term in NARROW_AUB_TERMS:
                    narrow.update(ids)
        rows.extend(
            {"period": record.period, "pt": term, "unique_reports": len(ids)}
            for term, ids in quarter_seen.items()
        )
    return broad, narrow, reports_by_term, pd.DataFrame(rows)


def build_explanatory_outcomes(reports_by_term: dict[str, set[str]]) -> dict[str, set[str]]:
    """Build frozen post-result PT decompositions without changing the primary phenotype."""
    missing = set(AUB_TERMS) - set(reports_by_term)
    if missing:
        raise ValueError(f"Missing term-level report sets: {sorted(missing)}")
    outcomes = {f"pt_{term.lower().replace(' ', '_')}": set(ids) for term, ids in reports_by_term.items()}
    for name, terms in EXPLANATORY_OUTCOME_GROUPS.items():
        outcomes[name] = set().union(*(reports_by_term[term] for term in terms))
    return outcomes


def parse_exact_faers_date(values: pd.Series) -> pd.Series:
    """Parse only complete YYYYMMDD dates; incomplete dates remain missing."""
    text = values.astype("string").str.strip().str.replace(r"\.0$", "", regex=True)
    exact = text.where(text.str.fullmatch(r"\d{8}", na=False))
    return pd.to_datetime(exact, format="%Y%m%d", errors="coerce")


def collect_therapy_starts(
    records, keys_by_drug: dict[str, set[tuple[str, str]]]
) -> dict[str, dict[str, pd.Timestamp]]:
    starts: dict[str, dict[str, pd.Timestamp]] = {drug: {} for drug in keys_by_drug}
    all_keys = set().union(*keys_by_drug.values())
    usecols = ["primaryid", "caseid", "dsg_drug_seq", "start_dt"]
    for record in records:
        for chunk in read_table_chunks(record.path, "THER", usecols):
            keys = list(zip(chunk["primaryid"].tolist(), chunk["dsg_drug_seq"].tolist()))
            keep = pd.Series([key in all_keys for key in keys], index=chunk.index)
            chunk = chunk.loc[keep].copy()
            if chunk.empty:
                continue
            chunk["start_date"] = parse_exact_faers_date(chunk["start_dt"])
            chunk = chunk[chunk["start_date"].notna()]
            for drug, drug_keys in keys_by_drug.items():
                dmask = pd.Series(
                    [key in drug_keys for key in zip(chunk["primaryid"], chunk["dsg_drug_seq"])],
                    index=chunk.index,
                )
                for pid, date in chunk.loc[dmask, ["primaryid", "start_date"]].itertuples(index=False):
                    old = starts[drug].get(pid)
                    if old is None or date < old:
                        starts[drug][pid] = date
    return starts


def summarize_time_to_onset(
    drug: str,
    event_ids: set[str],
    event_dates: dict[str, pd.Timestamp],
    therapy_starts: dict[str, pd.Timestamp],
) -> dict[str, object]:
    exact_event_ids = event_ids & set(event_dates)
    linked_ids = exact_event_ids & set(therapy_starts)
    intervals = pd.Series(
        [(event_dates[pid] - therapy_starts[pid]).days for pid in linked_ids], dtype="float64"
    )
    valid = intervals[(intervals >= 0) & (intervals <= MAX_PLAUSIBLE_TTO_DAYS)]
    return {
        "drug": drug,
        "event_reports": len(event_ids),
        "exact_event_date_reports": len(exact_event_ids),
        "linked_exact_start_and_event_reports": len(linked_ids),
        "valid_tto_reports": int(len(valid)),
        "valid_tto_fraction": float(len(valid) / len(event_ids)) if event_ids else None,
        "median_days": float(valid.median()) if len(valid) else None,
        "q1_days": float(valid.quantile(0.25)) if len(valid) else None,
        "q3_days": float(valid.quantile(0.75)) if len(valid) else None,
        "excluded_negative_or_gt_3650": int(len(intervals) - len(valid)),
    }


def collect_serious_outcomes(records, event_universe: set[str]) -> dict[str, set[str]]:
    by_code: dict[str, set[str]] = {}
    usecols = ["primaryid", "caseid", "outc_cod"]
    for record in records:
        for chunk in read_table_chunks(record.path, "OUTC", usecols):
            chunk = chunk[chunk["primaryid"].isin(event_universe)].copy()
            codes = chunk["outc_cod"].fillna("MISSING").str.strip().str.upper().replace("", "MISSING")
            for code, ids in chunk.assign(outc_normalized=codes).groupby("outc_normalized")["primaryid"]:
                by_code.setdefault(str(code), set()).update(ids.tolist())
    return by_code


def collect_indication_reports(records, ps_keys: dict[str, set[tuple[str, str]]]) -> dict[str, dict[str, set[str]]]:
    result = {kind: {drug: set() for drug in DRUG_ALIASES} for kind in ("any", "vte", "af")}
    all_keys = set().union(*ps_keys.values())
    usecols = ["primaryid", "caseid", "indi_drug_seq", "indi_pt"]
    for record in records:
        for chunk in read_table_chunks(record.path, "INDI", usecols):
            keys = list(zip(chunk["primaryid"].tolist(), chunk["indi_drug_seq"].tolist()))
            keep = pd.Series([key in all_keys for key in keys], index=chunk.index)
            chunk = chunk.loc[keep].copy()
            if chunk.empty:
                continue
            normalized = chunk["indi_pt"].str.upper().str.strip()
            for drug, drug_keys in ps_keys.items():
                dmask = pd.Series(
                    [key in drug_keys for key in zip(chunk["primaryid"], chunk["indi_drug_seq"])],
                    index=chunk.index,
                )
                result["any"][drug].update(
                    chunk.loc[dmask & normalized.notna() & normalized.ne(""), "primaryid"]
                )
                result["vte"][drug].update(chunk.loc[dmask & normalized.str.contains(VTE_PATTERN, na=False), "primaryid"])
                result["af"][drug].update(chunk.loc[dmask & normalized.str.contains(AF_PATTERN, na=False), "primaryid"])
    return result


def contrast(
    exposed: set[str], comparator: set[str], outcome: set[str], cohort: set[str],
    exposed_name: str, comparator_name: str, analysis: str,
) -> tuple[dict[str, object], set[str], set[str]]:
    e = (exposed & cohort) - comparator
    c = (comparator & cohort) - exposed
    overlap = exposed & comparator & cohort
    a, b = len(e & outcome), len(e - outcome)
    cc, d = len(c & outcome), len(c - outcome)
    stats = two_by_two_statistics(a, b, cc, d)
    stats.update(
        {
            "analysis": analysis,
            "exposed": exposed_name,
            "comparator": comparator_name,
            "excluded_dual_exposure_reports": len(overlap),
            "reconciliation_pass": bool(a + b == len(e) and cc + d == len(c) and e.isdisjoint(c)),
        }
    )
    if not stats["reconciliation_pass"]:
        raise RuntimeError(f"2x2 reconciliation failed for {analysis}")
    return stats, e, c


def period_range_from_manifest(manifest_path: Path) -> list[str]:
    manifest = pd.read_csv(manifest_path, dtype="string")
    return sorted(manifest["Period"].tolist())


def fda_date_to_period(value: int | float | str) -> str:
    digits = str(int(float(value)))
    if len(digits) != 8:
        raise ValueError(f"Invalid FDA_DT value: {value}")
    year, month = int(digits[:4]), int(digits[4:6])
    if month < 1 or month > 12:
        raise ValueError(f"Invalid FDA_DT month: {value}")
    return f"{year}Q{(month - 1) // 3 + 1}"


def run(project_root: Path, output_dir: Path, periods: list[str] | None = None) -> dict[str, object]:
    archives_dir = project_root / ".aris" / "data" / "faers" / "archives"
    manifest_path = project_root / ".aris" / "data" / "faers" / "manifest" / "download_manifest_sha256.csv"
    periods = periods or period_range_from_manifest(manifest_path)
    if not output_dir.is_absolute():
        output_dir = project_root / output_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    records = load_manifest(manifest_path, archives_dir, periods)
    audit = schema_audit(records)
    (output_dir / "schema_audit.json").write_text(json.dumps(audit, indent=2), encoding="utf-8")
    if audit["required_variable_coverage"] != 1.0:
        raise RuntimeError("Required-variable coverage is below 100%")

    demo, raw_by_quarter = load_and_deduplicate_demo(records)
    pd.DataFrame(demo.attrs.get("primaryid_fda_dt_conflicts", [])).to_csv(
        output_dir / "primaryid_fda_dt_correction_audit.csv", index=False
    )
    female = demo[demo["sex"].str.strip().str.upper().eq("F")].copy()
    age_ok = female[female["age_years"].notna()].copy()
    cohort_primary = set(age_ok.loc[age_ok["age_years"].between(15, 49), "primaryid"])
    cohort_age_12_55 = set(age_ok.loc[age_ok["age_years"].between(12, 55), "primaryid"])
    cohort_all_female = set(female["primaryid"])
    universe = cohort_all_female
    # Identical PRIMARYIDs can recur in adjacent published ZIPs. Assign the
    # retained report version by its FDA receipt/update date, not by whichever
    # quarterly archive happened to contain the repeated row last.
    if demo["fda_dt_num"].isna().any():
        raise ValueError("Retained FAERS cases contain missing FDA_DT; temporal assignment is undefined")
    period_by_id = dict(zip(demo["primaryid"], demo["fda_dt_num"].map(fda_date_to_period)))

    exposures, ps_keys, followup_fields = collect_drugs(records, universe)
    broad_aub, narrow_aub, reports_by_term, term_counts = collect_outcomes(records, universe)
    indications = collect_indication_reports(records, ps_keys)

    primary, riv_only, api_only = contrast(
        exposures["ps"]["rivaroxaban"], exposures["ps"]["apixaban"], broad_aub,
        cohort_primary, "rivaroxaban", "apixaban", "confirmatory_primary",
    )
    primary["minimum_detectable_or_80pct"] = minimum_detectable_or(
        primary["a_exposed_event"] + primary["b_exposed_non_event"],
        primary["c_comparator_event"] + primary["d_comparator_non_event"],
        primary["c_comparator_event"],
    )
    primary["confirmatory_estimable"] = bool(
        primary["a_exposed_event"] >= 20
        and primary["c_comparator_event"] >= 20
        and primary["minimum_detectable_or_80pct"] <= 2.0
    )
    (output_dir / "primary_contrast.json").write_text(json.dumps(primary, indent=2), encoding="utf-8")

    background_rows = []
    for drug in ("rivaroxaban", "apixaban"):
        exposed_ids = exposures["ps"][drug] & cohort_primary
        comparator_ids = cohort_primary - exposed_ids
        background = two_by_two_statistics(
            len(exposed_ids & broad_aub),
            len(exposed_ids - broad_aub),
            len(comparator_ids & broad_aub),
            len(comparator_ids - broad_aub),
        )
        background.update(
            {
                "analysis": f"background_{drug}_vs_all_other_reports",
                "exposed": drug,
                "comparator": "all_other_reports_in_eligible_faers_background",
                "eligible_background_reports": len(cohort_primary),
            }
        )
        background_rows.append(background)
    qvals = benjamini_hochberg(row["fisher_exact_p"] for row in background_rows)
    for row, q in zip(background_rows, qvals):
        row["bh_fdr_q_background_family"] = q
    pd.DataFrame(background_rows).to_csv(
        output_dir / "background_disproportionality.csv", index=False
    )

    demo_by_id = demo.set_index("primaryid", drop=False)
    event_dates = {
        pid: date for pid, date in zip(
            demo_by_id.index,
            parse_exact_faers_date(demo_by_id["event_dt"]),
        ) if pd.notna(date)
    }
    therapy_starts = collect_therapy_starts(records, ps_keys)
    event_ids_by_drug = {
        "rivaroxaban": riv_only & broad_aub,
        "apixaban": api_only & broad_aub,
    }

    occupation_codes = demo_by_id["occp_cod"].fillna("").str.strip().str.upper()

    def reporter_group(code: str) -> str:
        if code in HEALTHCARE_OCCUPATIONS:
            return "healthcare_professional"
        if code in CONSUMER_OCCUPATIONS:
            return "consumer"
        if code in LAWYER_OCCUPATIONS:
            return "lawyer"
        if not code:
            return "missing"
        return "other"

    reporter_groups = occupation_codes.map(reporter_group)

    def period_group(period: str) -> str:
        if period < "2021Q1":
            return "pre_2021q1"
        if period in SAFETY_COMMUNICATION_PERIODS:
            return "2021q1_q2"
        return "post_2021q2"

    characteristic_rows = []
    missingness_rows = []
    association_rows = []
    for drug, ids in (("rivaroxaban", riv_only), ("apixaban", api_only)):
        subset = demo_by_id.loc[sorted(ids)].copy()
        ages = subset["age_years"].dropna().astype(float)
        for statistic, value in (
            ("mean", ages.mean()),
            ("sd", ages.std(ddof=1)),
            ("median", ages.median()),
            ("q1", ages.quantile(0.25)),
            ("q3", ages.quantile(0.75)),
        ):
            characteristic_rows.append(
                {
                    "drug": drug,
                    "variable": "age_years",
                    "level_or_statistic": statistic,
                    "n": len(ages),
                    "fraction": None,
                    "value": float(value),
                }
            )
        drug_reporter = reporter_groups.loc[subset.index]
        country = subset["reporter_country"].fillna("").str.strip().str.upper()
        country_group = country.map(
            lambda value: "missing" if not value else ("United States" if value in {"US", "USA", "UNITED STATES"} else "non-US")
        )
        drug_period = pd.Series(
            [period_group(period_by_id[pid]) for pid in subset.index], index=subset.index
        )
        any_ind = indications["any"][drug]
        vte_ind = indications["vte"][drug]
        af_ind = indications["af"][drug]

        def indication_group(pid: str) -> str:
            if pid in vte_ind and pid in af_ind:
                return "VTE_and_AF"
            if pid in vte_ind:
                return "VTE"
            if pid in af_ind:
                return "AF"
            if pid in any_ind:
                return "other_recorded_indication"
            return "no_recorded_indication"

        indication_groups = pd.Series(
            [indication_group(pid) for pid in subset.index], index=subset.index
        )
        for variable, values in (
            ("reporter_group", drug_reporter),
            ("reporter_country_group", country_group),
            ("calendar_period", drug_period),
            ("indication_group", indication_groups),
        ):
            for level, count in values.value_counts(dropna=False).items():
                characteristic_rows.append(
                    {
                        "drug": drug,
                        "variable": variable,
                        "level_or_statistic": str(level),
                        "n": int(count),
                        "fraction": float(count / len(subset)),
                        "value": None,
                    }
                )
        frame = pd.DataFrame(
            {
                "drug": drug,
                "event_status": ["broad_outcome" if pid in broad_aub else "no_broad_outcome" for pid in subset.index],
                "reporter_group": drug_reporter.values,
                "reporter_country_group": country_group.values,
                "calendar_period": drug_period.values,
            }
        )
        association_rows.append(frame)

        female_drug_ids = (exposures["ps"][drug] & cohort_all_female)
        other = "apixaban" if drug == "rivaroxaban" else "rivaroxaban"
        female_drug_ids = female_drug_ids - exposures["ps"][other]
        female_drug = demo_by_id.loc[sorted(female_drug_ids)]
        missingness_rows.extend(
            [
                {
                    "population": f"exclusive_{drug}_primary_suspect_female_reports",
                    "variable": "age_missing_or_unconvertible",
                    "missing_n": int(female_drug["age_years"].isna().sum()),
                    "denominator_n": int(len(female_drug)),
                },
                {
                    "population": f"eligible_{drug}_primary_suspect_reports",
                    "variable": "reporter_occupation_missing",
                    "missing_n": int(drug_reporter.eq("missing").sum()),
                    "denominator_n": int(len(subset)),
                },
                {
                    "population": f"eligible_{drug}_primary_suspect_reports",
                    "variable": "reporter_country_missing",
                    "missing_n": int(country.eq("").sum()),
                    "denominator_n": int(len(subset)),
                },
                {
                    "population": f"eligible_{drug}_primary_suspect_reports",
                    "variable": "indication_not_recorded",
                    "missing_n": int(sum(pid not in any_ind for pid in subset.index)),
                    "denominator_n": int(len(subset)),
                },
                {
                    "population": f"eligible_{drug}_primary_suspect_reports",
                    "variable": "exact_event_date_missing",
                    "missing_n": int(sum(pid not in event_dates for pid in subset.index)),
                    "denominator_n": int(len(subset)),
                },
            ]
        )

    sex_code = demo["sex"].fillna("").str.strip().str.upper()
    missingness_rows.extend(
        [
            {
                "population": "all_deduplicated_faers_cases",
                "variable": "sex_missing",
                "missing_n": int(sex_code.eq("").sum()),
                "denominator_n": int(len(demo)),
            },
            {
                "population": "female_deduplicated_faers_cases",
                "variable": "age_missing_or_unconvertible",
                "missing_n": int(female["age_years"].isna().sum()),
                "denominator_n": int(len(female)),
            },
        ]
    )
    characteristics = pd.DataFrame(characteristic_rows)
    characteristics.to_csv(output_dir / "cohort_characteristics.csv", index=False)
    characteristics.loc[
        characteristics["variable"].eq("indication_group")
    ].to_csv(output_dir / "indication_distribution.csv", index=False)
    associations = pd.concat(association_rows, ignore_index=True)
    associations.groupby(
        ["drug", "event_status", "reporter_group", "reporter_country_group", "calendar_period"],
        dropna=False,
    ).size().reset_index(name="n").to_csv(
        output_dir / "reporter_source_associations.csv", index=False
    )

    for drug, event_ids in event_ids_by_drug.items():
        missingness_rows.extend(
            [
                {
                    "population": f"eligible_{drug}_broad_outcome_reports",
                    "variable": "exact_event_date_missing",
                    "missing_n": int(sum(pid not in event_dates for pid in event_ids)),
                    "denominator_n": int(len(event_ids)),
                },
                {
                    "population": f"eligible_{drug}_broad_outcome_reports",
                    "variable": "exact_drug_start_date_missing",
                    "missing_n": int(sum(pid not in therapy_starts[drug] for pid in event_ids)),
                    "denominator_n": int(len(event_ids)),
                },
            ]
        )
    missingness = pd.DataFrame(missingness_rows)
    missingness["missing_fraction"] = missingness["missing_n"] / missingness["denominator_n"]
    missingness.to_csv(output_dir / "missingness_by_drug.csv", index=False)

    tto_rows = [
        summarize_time_to_onset(drug, event_ids_by_drug[drug], event_dates, therapy_starts[drug])
        for drug in ("rivaroxaban", "apixaban")
    ]
    pd.DataFrame(tto_rows).to_csv(output_dir / "time_to_onset_summary.csv", index=False)

    followup_rows = []
    for drug, event_ids in event_ids_by_drug.items():
        for field in ("dechallenge_recorded", "rechallenge_recorded"):
            n = len(event_ids & followup_fields[field][drug])
            followup_rows.append({
                "drug": drug,
                "field": field,
                "event_reports": len(event_ids),
                "recorded_reports": n,
                "completeness": n / len(event_ids) if event_ids else None,
                "analysis_allowed_at_20pct": bool(event_ids and n / len(event_ids) >= 0.20),
            })
    pd.DataFrame(followup_rows).to_csv(output_dir / "dechallenge_rechallenge_completeness.csv", index=False)

    descriptor_rows = []
    descriptor_fields = ["rept_cod", "occp_cod", "reporter_country", "occr_country"]
    for drug, event_ids in event_ids_by_drug.items():
        subset = demo_by_id.loc[sorted(event_ids)]
        for field in descriptor_fields:
            values = subset[field].fillna("MISSING").str.strip().replace("", "MISSING")
            for value, count in values.value_counts(dropna=False).items():
                descriptor_rows.append({
                    "drug": drug, "field": field, "value": value,
                    "n": int(count), "fraction": float(count / len(subset)) if len(subset) else None,
                })
    pd.DataFrame(descriptor_rows).to_csv(output_dir / "event_report_descriptors.csv", index=False)

    serious_by_code = collect_serious_outcomes(records, set().union(*event_ids_by_drug.values()))
    serious_rows = []
    for drug, event_ids in event_ids_by_drug.items():
        for code, ids in sorted(serious_by_code.items()):
            serious_rows.append({
                "drug": drug, "outc_cod": code, "n": len(event_ids & ids),
                "fraction": len(event_ids & ids) / len(event_ids) if event_ids else None,
            })
    pd.DataFrame(serious_rows).to_csv(output_dir / "serious_outcome_descriptors.csv", index=False)

    quarter_rows = []
    for period in periods:
        r = {pid for pid in riv_only if period_by_id.get(pid) == period}
        a = {pid for pid in api_only if period_by_id.get(pid) == period}
        quarter_rows.append(
            {
                "period": period,
                "rivaroxaban_reports": len(r),
                "rivaroxaban_aub": len(r & broad_aub),
                "apixaban_reports": len(a),
                "apixaban_aub": len(a & broad_aub),
            }
        )
    pd.DataFrame(quarter_rows).to_csv(output_dir / "quarterly_active_comparator_counts.csv", index=False)
    term_counts.to_csv(output_dir / "aub_term_counts_by_quarter.csv", index=False)

    family_a_specs = [
        ("rivaroxaban", "edoxaban"),
        ("rivaroxaban", "dabigatran"),
        ("rivaroxaban", "warfarin"),
        ("apixaban", "edoxaban"),
        ("apixaban", "dabigatran"),
        ("apixaban", "warfarin"),
    ]
    family_a = []
    for exposed_name, comparator_name in family_a_specs:
        result, _, _ = contrast(
            exposures["ps"][exposed_name], exposures["ps"][comparator_name], broad_aub,
            cohort_primary, exposed_name, comparator_name, f"family_A_{exposed_name}_vs_{comparator_name}",
        )
        family_a.append(result)
    qvals = benjamini_hochberg(row["fisher_exact_p"] for row in family_a)
    for row, q in zip(family_a, qvals):
        row["bh_fdr_q"] = q
    family_a_frame = pd.DataFrame(family_a)
    family_a_frame.to_csv(output_dir / "family_A_drug_contrasts.csv", index=False)
    family_a_frame.to_csv(output_dir / "symmetric_active_comparators.csv", index=False)

    base_riv, base_api = exposures["ps"]["rivaroxaban"], exposures["ps"]["apixaban"]
    occupations = demo_by_id["occp_cod"].fillna("").str.strip().str.upper()
    healthcare_cohort = cohort_primary & set(demo_by_id.index[occupations.isin(HEALTHCARE_OCCUPATIONS)])
    consumer_cohort = cohort_primary & set(demo_by_id.index[occupations.isin(CONSUMER_OCCUPATIONS)])

    reporter_tables = []
    reporter_results = []
    for name, stratum in (
        ("healthcare_professional", healthcare_cohort),
        ("consumer", consumer_cohort),
    ):
        result, _, _ = contrast(
            base_riv, base_api, broad_aub, stratum,
            "rivaroxaban", "apixaban", f"reporter_{name}",
        )
        reporter_results.append({"stratum": name, **result})
        reporter_tables.append(
            (
                name,
                result["a_exposed_event"],
                result["b_exposed_non_event"],
                result["c_comparator_event"],
                result["d_comparator_non_event"],
            )
        )
    reporter_heterogeneity = log_or_heterogeneity(reporter_tables)
    for row in reporter_results:
        row["heterogeneity_method"] = reporter_heterogeneity["method"]
        row["heterogeneity_q"] = reporter_heterogeneity["q_statistic"]
        row["heterogeneity_df"] = reporter_heterogeneity["degrees_freedom"]
        row["heterogeneity_p"] = reporter_heterogeneity["p_heterogeneity"]
    pd.DataFrame(reporter_results).to_csv(
        output_dir / "reporter_source_interaction.csv", index=False
    )

    calendar_specs = [
        ("pre_2021q1", {pid for pid in cohort_primary if period_by_id.get(pid) < "2021Q1"}),
        ("2021q1_q2", {pid for pid in cohort_primary if period_by_id.get(pid) in SAFETY_COMMUNICATION_PERIODS}),
        ("post_2021q2", {pid for pid in cohort_primary if period_by_id.get(pid) > "2021Q2"}),
    ]
    calendar_rows = []
    calendar_tables = []
    for name, stratum in calendar_specs:
        result, _, _ = contrast(
            base_riv, base_api, broad_aub, stratum,
            "rivaroxaban", "apixaban", f"calendar_{name}",
        )
        calendar_rows.append({"period_group": name, **result})
        calendar_tables.append(
            (
                name,
                result["a_exposed_event"],
                result["b_exposed_non_event"],
                result["c_comparator_event"],
                result["d_comparator_non_event"],
            )
        )
    calendar_heterogeneity = log_or_heterogeneity(calendar_tables)
    pd.DataFrame(calendar_rows).to_csv(output_dir / "calendar_period_strata.csv", index=False)
    pd.DataFrame(
        [
            {
                "method": calendar_heterogeneity["method"],
                "q_statistic": calendar_heterogeneity["q_statistic"],
                "degrees_freedom": calendar_heterogeneity["degrees_freedom"],
                "p_heterogeneity": calendar_heterogeneity["p_heterogeneity"],
            }
        ]
    ).to_csv(output_dir / "calendar_period_heterogeneity.csv", index=False)

    sensitivity_inputs = [
        ("all_suspect", exposures["all_suspect"]["rivaroxaban"], exposures["all_suspect"]["apixaban"], broad_aub, cohort_primary),
        ("age_12_55", base_riv, base_api, broad_aub, cohort_age_12_55),
        ("all_females", base_riv, base_api, broad_aub, cohort_all_female),
        ("vte_only", indications["vte"]["rivaroxaban"], indications["vte"]["apixaban"], broad_aub, cohort_primary),
        ("atrial_fibrillation_only", indications["af"]["rivaroxaban"], indications["af"]["apixaban"], broad_aub, cohort_primary),
        ("narrow_pt", base_riv, base_api, narrow_aub, cohort_primary),
        ("healthcare_professional_reporter", base_riv, base_api, broad_aub, healthcare_cohort),
        ("consumer_reporter", base_riv, base_api, broad_aub, consumer_cohort),
    ]
    family_c = []
    for name, riv, api, outcome, cohort in sensitivity_inputs:
        result, _, _ = contrast(riv, api, outcome, cohort, "rivaroxaban", "apixaban", name)
        family_c.append(result)
    qvals = benjamini_hochberg(row["fisher_exact_p"] for row in family_c)
    for row, q in zip(family_c, qvals):
        row["bh_fdr_q"] = q
    pd.DataFrame(family_c).to_csv(output_dir / "family_C_sensitivities.csv", index=False)

    # These decompositions were added only after the narrow-phenotype sensitivity
    # reversed direction. They explain heterogeneity and never replace the frozen
    # broad primary contrast. Multiplicity is controlled across the entire family.
    phenotype_rows = []
    for name, outcome in build_explanatory_outcomes(reports_by_term).items():
        result, _, _ = contrast(
            base_riv, base_api, outcome, cohort_primary,
            "rivaroxaban", "apixaban", f"post_result_{name}",
        )
        result["analysis_status"] = "post_result_explanatory"
        phenotype_rows.append(result)
    qvals = benjamini_hochberg(row["fisher_exact_p"] for row in phenotype_rows)
    for row, q in zip(phenotype_rows, qvals):
        row["bh_fdr_q"] = q
    pd.DataFrame(phenotype_rows).to_csv(
        output_dir / "post_result_phenotype_decomposition.csv", index=False
    )

    sparse_rows = []
    for family_name, rows in (
        ("primary", [primary]),
        ("background", background_rows),
        ("secondary_active_comparator", family_a),
        ("prespecified_sensitivity", family_c),
        ("post_result_outcome_decomposition", phenotype_rows),
        ("calendar_period", calendar_rows),
    ):
        for row in rows:
            cells = [
                row["a_exposed_event"], row["b_exposed_non_event"],
                row["c_comparator_event"], row["d_comparator_non_event"],
            ]
            if min(cells) < 20:
                sparse_rows.append({"analysis_family": family_name, **row})
    pd.DataFrame(sparse_rows).to_csv(output_dir / "sparse_exact_intervals.csv", index=False)

    flow = pd.DataFrame(
        [
            {"stage": "raw_demo_rows", "n": sum(int(x["raw_demo_rows"]) for x in raw_by_quarter)},
            {"stage": "deduplicated_cases", "n": len(demo)},
            {"stage": "female_cases", "n": len(female)},
            {"stage": "female_convertible_age", "n": len(age_ok)},
            {"stage": "female_age_15_49", "n": len(cohort_primary)},
            {"stage": "exclusive_rivaroxaban_primary", "n": len(riv_only)},
            {"stage": "exclusive_apixaban_primary", "n": len(api_only)},
        ]
    )
    flow.to_csv(output_dir / "cohort_flow.csv", index=False)

    metadata = {
        "run_utc": datetime.now(timezone.utc).isoformat(),
        "python": sys.version,
        "platform": platform.platform(),
        "pandas": pd.__version__,
        "numpy": np.__version__,
        "periods": periods,
        "aub_terms": sorted(AUB_TERMS),
        "narrow_aub_terms": sorted(NARROW_AUB_TERMS),
        "post_result_explanatory_groups": {
            key: sorted(value) for key, value in EXPLANATORY_OUTCOME_GROUPS.items()
        },
        "drug_aliases": {key: list(value) for key, value in DRUG_ALIASES.items()},
        "safety_communication_exclusion": sorted(SAFETY_COMMUNICATION_PERIODS),
        "reporter_occupation_groups": {
            "healthcare_professional": sorted(HEALTHCARE_OCCUPATIONS),
            "consumer": sorted(CONSUMER_OCCUPATIONS),
        },
        "time_to_onset_rule": "exact YYYYMMDD dates; retain 0 through 3650 days inclusive",
    }
    (output_dir / "run_metadata.json").write_text(json.dumps(metadata, indent=2), encoding="utf-8")
    return {"output_dir": str(output_dir), "period_count": len(periods), "primary": primary}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--output-dir", type=Path, default=Path("experiment/main/results/faers_primary"))
    parser.add_argument("--periods", nargs="+")
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    print(json.dumps(run(args.project_root.resolve(), args.output_dir, args.periods), indent=2))
