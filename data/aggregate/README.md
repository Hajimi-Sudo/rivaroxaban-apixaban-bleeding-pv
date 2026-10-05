# Aggregate data tables

These files are the frozen aggregate outputs underlying the manuscript tables and figures. They contain no individual case safety reports and no patient-level data.

- `primary_contrast.json`: primary FAERS 2x2 cells, ROR, PRR, confidence intervals, Fisher test and minimum detectable odds ratio.
- `cohort_flow.csv`: FAERS cohort funnel.
- `background_disproportionality.csv`: conventional full-FAERS-background analyses for both index drugs.
- `family_A_drug_contrasts.csv` and `symmetric_active_comparators.csv`: six symmetric secondary active-comparator analyses.
- `family_B_external_replication.csv`: separate Canada Vigilance and JADER estimates.
- `family_C_sensitivities.csv`: prespecified FAERS sensitivity analyses.
- `post_result_phenotype_decomposition.csv`: eight-test post-result explanatory outcome-definition family and BH-adjusted q values.
- `cohort_characteristics.csv` and `missingness_by_drug.csv`: aggregate cohort characteristics and missingness denominators.
- `reporter_source_interaction.csv`, `reporter_country_strata.csv`, and `calendar_period_strata.csv`: reporter-source, country, and calendar-period contrasts with heterogeneity records.
- `jader_age_rule_sensitivity.csv`: JADER age-rule sensitivity analyses.
- `sparse_exact_intervals.csv`: conditional maximum-likelihood odds ratios and exact confidence intervals for sparse tables.
- `quarterly_active_comparator_counts.csv`: quarterly eligible-report and broad-outcome counts used in Figure 2.
- `time_to_onset_summary.csv`: exact-date time-to-onset completeness and descriptive summaries.
- `dechallenge_rechallenge_completeness.csv`: completeness by drug and field.
- `dechallenge_rechallenge_code_distribution.csv`: raw Y/N/U/D code distributions.
- `serious_outcome_descriptors.csv`: non-exclusive serious-outcome code counts.

All effect estimates describe reporting disproportionality and must not be interpreted as incidence, prevalence, comparative treatment risk, or causality.

The release intentionally excludes raw individual case safety reports, individual report identifiers, and the identifier-level deduplication audit.
