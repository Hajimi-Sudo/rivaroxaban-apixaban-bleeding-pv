"""Database-specific Canada Vigilance replication and JADER directional validation."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from faers_pilot import AUB_TERMS, DRUG_ALIASES, normalize_text
from pv_stats import benjamini_hochberg, two_by_two_statistics


CANADA_REPORT_COLUMNS = [
    "report_id", "report_no", "version_no", "date_received", "date_initial_received", "mah_no",
    "report_type_code", "report_type_eng", "report_type_fr", "gender_code", "gender_eng", "gender_fr",
    "age", "age_y", "age_unit_eng", "age_unit_fr", "outcome_code", "outcome_eng", "outcome_fr",
    "weight", "weight_unit_eng", "weight_unit_fr", "height", "height_unit_eng", "height_unit_fr",
    "seriousness_code", "seriousness_eng", "seriousness_fr", "death", "disability", "congenital_anomaly",
    "life_threatening", "hospital_required", "other_medically_important", "reporter_type_eng",
    "reporter_type_fr", "source_code", "source_eng", "source_fr", "e2b_id", "authority_number", "company_number",
]
CANADA_DRUG_COLUMNS = [
    "report_drug_id", "report_id", "drug_product_id", "drugname", "role_eng", "role_fr", "route_eng", "route_fr",
    "dose", "dose_unit_eng", "dose_unit_fr", "frequency", "freq_time", "frequency_time_eng", "frequency_time_fr",
    "freq_unit_eng", "freq_unit_fr", "therapy_duration", "therapy_duration_unit_eng", "therapy_duration_unit_fr",
    "dosageform_eng", "dosageform_fr",
]
CANADA_INGREDIENT_COLUMNS = ["ingredient_row_id", "drug_product_id", "drugname", "ingredient_id", "ingredient_name"]
CANADA_REACTION_COLUMNS = [
    "reaction_id", "report_id", "duration", "duration_unit_eng", "duration_unit_fr", "pt_name_eng", "pt_name_fr",
    "soc_name_eng", "soc_name_fr", "meddra_version",
]
CANADA_LINK_COLUMNS = ["report_link_id", "report_id", "record_type_eng", "record_type_fr", "report_link_no"]

JADER_DRUG_ALIASES = {
    "rivaroxaban": ("RIVAROXABAN", "リバーロキサバン"),
    "apixaban": ("APIXABAN", "アピキサバン"),
    "edoxaban": ("EDOXABAN", "エドキサバン"),
    "dabigatran": ("DABIGATRAN", "ダビガトラン"),
    "warfarin": ("WARFARIN", "ワルファリン"),
}
JADER_AUB_TERMS = {"重度月経出血", "月経中間期出血", "子宮出血", "腟出血", "異常子宮出血"}


class UnionFind:
    def __init__(self) -> None:
        self.parent: dict[str, str] = {}

    def find(self, item: str) -> str:
        self.parent.setdefault(item, item)
        if self.parent[item] != item:
            self.parent[item] = self.find(self.parent[item])
        return self.parent[item]

    def union(self, left: str, right: str) -> None:
        a, b = self.find(left), self.find(right)
        if a != b:
            self.parent[max(a, b)] = min(a, b)


def read_canada(path: Path, names: list[str], chunksize: int | None = None):
    return pd.read_csv(
        path, sep="$", quotechar='"', header=None, names=names, dtype="string",
        keep_default_na=False, encoding="utf-8", chunksize=chunksize,
    )


def canada_canonical_map(data_dir: Path, reports: pd.DataFrame) -> tuple[dict[str, str], pd.DataFrame]:
    links = read_canada(data_dir / "report_links.txt", CANADA_LINK_COLUMNS)
    uf = UnionFind()
    # "Linked" indicates related reports and is not synonymous with a duplicate.
    # Only rows explicitly coded Duplicate are collapsed.
    duplicates = links[
        links["record_type_eng"].str.strip().str.upper().eq("DUPLICATE")
        & links["report_id"].str.fullmatch(r"\d+", na=False)
        & links["report_link_no"].str.fullmatch(r"\d+", na=False)
    ]
    duplicates = duplicates.copy()
    duplicates["left_num"] = pd.to_numeric(duplicates["report_id"], errors="coerce")
    duplicates["right_num"] = pd.to_numeric(duplicates["report_link_no"], errors="coerce")
    duplicates = duplicates.dropna(subset=["left_num", "right_num"])
    for row in duplicates.itertuples(index=False):
        left_num, right_num = int(row.left_num), int(row.right_num)
        if left_num <= 0 or right_num <= 0:
            continue
        left, right = str(left_num), str(right_num)
        uf.union(left, right)
    reports = reports.copy()
    reports["report_id"] = reports["report_id"].str.lstrip("0").replace("", "0")
    reports["group"] = reports["report_id"].map(lambda x: uf.find(str(x)) if str(x) in uf.parent else str(x))
    reports["date_received_parsed"] = pd.to_datetime(reports["date_received"], format="%d-%b-%y", errors="coerce")
    reports["report_id_num"] = pd.to_numeric(reports["report_id"], errors="coerce")
    ordered = reports.sort_values(["group", "date_received_parsed", "report_id_num"], kind="mergesort")
    canonical_rows = ordered.drop_duplicates("group", keep="last").copy()
    canonical_by_group = dict(zip(canonical_rows["group"], canonical_rows["report_id"]))
    mapping = {row.report_id: canonical_by_group[row.group] for row in reports.itertuples(index=False)}
    return mapping, canonical_rows


def active_comparator_result(
    exposures: dict[str, set[str]], outcome: set[str], exposed: str = "rivaroxaban", comparator: str = "apixaban"
) -> dict[str, object]:
    e = exposures[exposed] - exposures[comparator]
    c = exposures[comparator] - exposures[exposed]
    a, b = len(e & outcome), len(e - outcome)
    cc, d = len(c & outcome), len(c - outcome)
    result = two_by_two_statistics(a, b, cc, d)
    result.update(
        {
            "exposed": exposed,
            "comparator": comparator,
            "excluded_dual_exposure_reports": len(exposures[exposed] & exposures[comparator]),
            "reconciliation_pass": bool(a + b == len(e) and cc + d == len(c) and e.isdisjoint(c)),
        }
    )
    if not result["reconciliation_pass"]:
        raise RuntimeError("External active-comparator reconciliation failure")
    return result


def run_canada(project_root: Path, output_dir: Path) -> dict[str, object]:
    data_dir = project_root / "data" / "raw" / "canada" / "cvponline_extract_20260331"
    reports = read_canada(data_dir / "reports.txt", CANADA_REPORT_COLUMNS)
    mapping, canonical = canada_canonical_map(data_dir, reports)
    canonical["age_y_num"] = pd.to_numeric(canonical["age_y"], errors="coerce")
    eligible = set(
        canonical.loc[
            canonical["gender_eng"].str.strip().str.upper().eq("FEMALE")
            & canonical["age_y_num"].between(15, 49),
            "report_id",
        ]
    )

    ingredients = read_canada(data_dir / "drug_product_ingredients.txt", CANADA_INGREDIENT_COLUMNS)
    ingredient_norm = normalize_text(ingredients["ingredient_name"])
    product_ids = {}
    for drug, aliases in DRUG_ALIASES.items():
        pattern = "|".join(aliases)
        product_ids[drug] = set(ingredients.loc[ingredient_norm.str.contains(pattern, regex=True, na=False), "drug_product_id"])

    exposures = {drug: set() for drug in DRUG_ALIASES}
    for chunk in read_canada(data_dir / "report_drug.txt", CANADA_DRUG_COLUMNS, chunksize=500_000):
        chunk = chunk[chunk["role_eng"].str.strip().str.upper().eq("SUSPECT")].copy()
        chunk["source_id"] = chunk["report_id"].str.lstrip("0").replace("", "0")
        chunk["canonical_id"] = chunk["source_id"].map(mapping).fillna(chunk["source_id"])
        chunk = chunk[chunk["canonical_id"].isin(eligible)]
        name_norm = normalize_text(chunk["drugname"])
        for drug, aliases in DRUG_ALIASES.items():
            pattern = "|".join(aliases)
            matched = chunk["drug_product_id"].isin(product_ids[drug]) | name_norm.str.contains(pattern, regex=True, na=False)
            exposures[drug].update(chunk.loc[matched, "canonical_id"])

    outcome: set[str] = set()
    for chunk in read_canada(data_dir / "reactions.txt", CANADA_REACTION_COLUMNS, chunksize=500_000):
        chunk["source_id"] = chunk["report_id"].str.lstrip("0").replace("", "0")
        chunk["canonical_id"] = chunk["source_id"].map(mapping).fillna(chunk["source_id"])
        matched = chunk["canonical_id"].isin(eligible) & chunk["pt_name_eng"].str.strip().str.upper().isin(AUB_TERMS)
        outcome.update(chunk.loc[matched, "canonical_id"])

    result = active_comparator_result(exposures, outcome)
    result.update(
        {
            "database": "Canada Vigilance",
            "extract_date": "2026-03-31",
            "eligible_reports": len(eligible),
            "duplicate_source_reports_collapsed": int(len(reports) - len(canonical)),
        }
    )
    counts = [
        {"drug": drug, "eligible_suspect_reports": len(ids), "aub_reports": len(ids & outcome), "stable_20": len(ids & outcome) >= 20}
        for drug, ids in exposures.items()
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(counts).to_csv(output_dir / "canada_drug_event_counts.csv", index=False)
    (output_dir / "canada_primary_contrast.json").write_text(json.dumps(result, indent=2), encoding="utf-8")
    return result


def jader_age_20_49(age: pd.Series) -> pd.Series:
    normalized = age.astype("string").str.strip()
    return normalized.isin({"20歳代", "30歳代", "40歳代"})


def run_jader(project_root: Path, output_dir: Path) -> dict[str, object]:
    data_dir = project_root / "data" / "raw" / "jader" / "extracted"
    demo = pd.read_csv(data_dir / "demo202607.csv", dtype="string", encoding="cp932", keep_default_na=False)
    demo["report_n"] = pd.to_numeric(demo["報告回数"], errors="raise")
    demo = demo.sort_values(["識別番号", "report_n"], kind="mergesort").drop_duplicates("識別番号", keep="last")
    eligible_demo = demo[demo["性別"].eq("女性") & jader_age_20_49(demo["年齢"])].copy()
    female_age_missing_or_unparseable = int(
        (demo["性別"].eq("女性") & ~jader_age_20_49(demo["年齢"]) & ~demo["年齢"].str.match(r"^(?:0|10|50|60|70|80|90|100)歳代$")).sum()
    )
    retained = set(zip(eligible_demo["識別番号"], eligible_demo["報告回数"]))
    eligible_ids = set(eligible_demo["識別番号"])

    drug_table = pd.read_csv(data_dir / "drug202607.csv", dtype="string", encoding="cp932", keep_default_na=False)
    keys = pd.Series(list(zip(drug_table["識別番号"], drug_table["報告回数"])), index=drug_table.index)
    drug_table = drug_table[keys.isin(retained) & drug_table["医薬品の関与"].eq("被疑薬")].copy()
    combined = drug_table["医薬品（一般名）"].str.upper() + " " + drug_table["医薬品（販売名）"].str.upper()
    exposures = {drug: set() for drug in JADER_DRUG_ALIASES}
    for drug, aliases in JADER_DRUG_ALIASES.items():
        matched = combined.apply(lambda value: any(alias in value for alias in aliases))
        exposures[drug].update(drug_table.loc[matched, "識別番号"])

    reactions = pd.read_csv(data_dir / "reac202607.csv", dtype="string", encoding="cp932", keep_default_na=False)
    reaction_keys = pd.Series(list(zip(reactions["識別番号"], reactions["報告回数"])), index=reactions.index)
    reactions = reactions[reaction_keys.isin(retained)]
    outcome = set(reactions.loc[reactions["有害事象"].isin(JADER_AUB_TERMS), "識別番号"])

    result = active_comparator_result(exposures, outcome)
    result.update(
        {
            "database": "JADER",
            "extract": "2026-07",
            "eligible_reports": len(eligible_ids),
            "age_proxy": "female 20s, 30s, or 40s",
            "directional_only": True,
            "female_age_missing_or_unparseable": female_age_missing_or_unparseable,
        }
    )
    counts = [
        {"drug": drug, "eligible_suspect_reports": len(ids), "aub_reports": len(ids & outcome), "stable_20": len(ids & outcome) >= 20}
        for drug, ids in exposures.items()
    ]
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(counts).to_csv(output_dir / "jader_drug_event_counts.csv", index=False)
    (output_dir / "jader_primary_contrast.json").write_text(json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8")
    return result


def run(project_root: Path, output_dir: Path) -> dict[str, object]:
    if not output_dir.is_absolute():
        output_dir = project_root / output_dir
    canada = run_canada(project_root, output_dir)
    jader = run_jader(project_root, output_dir)
    qvals = benjamini_hochberg([canada["fisher_exact_p"], jader["fisher_exact_p"]])
    comparison = pd.DataFrame(
        [
            {**canada, "bh_fdr_q_family_B": qvals[0]},
            {**jader, "bh_fdr_q_family_B": qvals[1]},
        ]
    )
    comparison.to_csv(output_dir / "family_B_external_replication.csv", index=False)
    return {"canada": canada, "jader": jader, "family_B_q": qvals}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--output-dir", type=Path, default=Path("results/external_validation"))
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    print(json.dumps(run(args.project_root.resolve(), args.output_dir), indent=2, ensure_ascii=False))
