import pandas as pd

from faers_full import (
    AUB_TERMS,
    EXPLANATORY_OUTCOME_GROUPS,
    build_explanatory_outcomes,
    parse_exact_faers_date,
    summarize_time_to_onset,
)
from dechallenge_followup import summarize_code_sets


def test_explanatory_outcomes_preserve_individual_terms_and_fixed_groups():
    reports_by_term = {
        term: {f"{term}-shared", f"{term}-only"} for term in AUB_TERMS
    }

    outcomes = build_explanatory_outcomes(reports_by_term)

    assert len([name for name in outcomes if name.startswith("pt_")]) == len(AUB_TERMS)
    for name, terms in EXPLANATORY_OUTCOME_GROUPS.items():
        expected = set().union(*(reports_by_term[term] for term in terms))
        assert outcomes[name] == expected


def test_explanatory_outcomes_reject_missing_term_set():
    reports_by_term = {term: set() for term in AUB_TERMS}
    reports_by_term.pop(next(iter(AUB_TERMS)))

    try:
        build_explanatory_outcomes(reports_by_term)
    except ValueError as exc:
        assert "Missing term-level report sets" in str(exc)
    else:
        raise AssertionError("missing terms must be rejected")


def test_parse_exact_faers_date_rejects_partial_and_invalid_dates():
    parsed = parse_exact_faers_date(pd.Series(["20200102", "202001", "bad", None]))
    assert parsed.iloc[0] == pd.Timestamp("2020-01-02")
    assert parsed.iloc[1:].isna().all()


def test_tto_summary_excludes_negative_and_implausibly_long_intervals():
    event_ids = {"a", "b", "c", "d"}
    event_dates = {
        "a": pd.Timestamp("2020-01-11"),
        "b": pd.Timestamp("2020-01-01"),
        "c": pd.Timestamp("2035-01-01"),
    }
    therapy_starts = {
        "a": pd.Timestamp("2020-01-01"),
        "b": pd.Timestamp("2020-01-02"),
        "c": pd.Timestamp("2020-01-01"),
    }
    result = summarize_time_to_onset("rivaroxaban", event_ids, event_dates, therapy_starts)
    assert result["valid_tto_reports"] == 1
    assert result["median_days"] == 10
    assert result["excluded_negative_or_gt_3650"] == 2


def test_dechallenge_code_summary_is_report_level_and_flags_conflicts():
    rows, audit = summarize_code_sets(
        "rivaroxaban",
        "dechal",
        {"a", "b", "c"},
        {"Y": {"a", "b"}, "N": {"b", "outside"}},
    )
    by_code = {row["code"]: row for row in rows}
    assert by_code["Y"]["reports_with_code"] == 2
    assert by_code["N"]["reports_with_code"] == 1
    assert audit["reports_with_any_recorded_code"] == 2
    assert audit["reports_with_conflicting_multiple_codes"] == 1
