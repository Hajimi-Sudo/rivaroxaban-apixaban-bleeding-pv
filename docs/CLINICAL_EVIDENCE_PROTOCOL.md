# Clinical evidence layer — protocol v1.0 (locked 2026-10-18)

This protocol defines a focused, reproducible evidence synthesis that extends rather than duplicates the 2026 JTH systematic review/meta-analysis. That review included 21 studies (10 randomized trials, 11 cohorts) searched through June 2024 and found that outcome definition strongly affects apparent AUB risk. The present update prioritizes post-June-2024 records, direct rivaroxaban/apixaban comparisons, endpoint-by-design separation, and a transparent case-level FAERS complement.

## Focused question

In people of reproductive potential receiving rivaroxaban or apixaban for VTE or AF, how do **direct comparative clinical estimates** differ across: (a) patient-centered heavy menstrual/AUB outcomes; (b) clinician-coded uterine bleeding; and (c) severe uterine bleeding requiring care, transfusion, or surgery? Secondary goal: identify whether published data permit a fair direct comparison within the same indication and outcome definition.

## Sources and search

Search MEDLINE/PubMed, Embase (if institutionally available), CENTRAL, ClinicalTrials.gov, and forward/backward citations from the 2026 review, FDA Sentinel study, and HEMBLED registry. The reproducible PubMed search was run on 2026-10-18 through NCBI E-utilities and is archived in `clinical_review_2026-10-18/search_protocol.json`, with all 46 retrieved records and a screening log. The exact PubMed query was:

```text
(rivaroxaban[Title/Abstract] OR apixaban[Title/Abstract] OR "direct oral anticoagulant"[Title/Abstract]) AND ("abnormal uterine bleeding"[Title/Abstract] OR "heavy menstrual bleeding"[Title/Abstract] OR menorrhagia[Title/Abstract] OR "uterine bleeding"[Title/Abstract]) AND (female[Title/Abstract] OR women[Title/Abstract] OR menstru*[Title/Abstract])
```

Record exact search strings, run dates, database coverage, deduplication, and all full-text exclusions. Search from inception to the final freeze date, emphasizing studies published after the June 2024 search cutoff of the existing review. Use PRISMA 2020 for reporting and its flow diagram. If Embase, CENTRAL, or ClinicalTrials.gov access is unavailable, state that limitation rather than claiming those sources were searched.

## Eligibility and extraction

- Include randomized, prospective, or retrospective studies with exposure denominators and extractable rivaroxaban/apixaban uterine or menstrual bleeding outcomes; mark whether comparisons are direct or indirect.
- Separate VTE from AF and avoid pooling mixed indications without stratified data. Record age, menstrual status, pregnancy exclusion, baseline bleeding, hormone therapy, dose, follow-up, endpoint ascertainment, and covariate control.
- Extract numerator, exposed N/person-time, adjusted and unadjusted estimates, uncertainty, and exact outcome definition. Cross-link multiple reports of the same cohort to avoid double counting.
- Appraise randomized and observational designs with appropriate risk-of-bias tools, with decisions and justifications retained. At least two human-author screeners and conflict resolution should be documented for a formal systematic review; AI can assist retrieval/extraction but cannot be represented as independent human screening.

## Synthesis rule

Start with an endpoint-by-design evidence map. A meta-analysis is **conditional**, not automatic: pool only directly comparable populations, exposure contrasts, outcomes, estimands, and nonoverlapping cohorts. If these conditions fail, use structured narrative synthesis with interval plots, not a pooled headline ranking. Do not combine RORs with clinical HRs/RRs or treat the spontaneous-report systems as trials/cohorts.

## Case-level FAERS complement

The primary FAERS analysis will be complemented by an aggregate-only case-level description using the same latest-case, female 15--49 years, primary-suspect rivaroxaban-only and apixaban-only cohorts. The analysis will summarize reporter type, occupation, reporter and occurrence country, coded serious outcomes, indication groups, and dechallenge/rechallenge fields separately for the active-comparator cohorts and historical broad-composite event reports. No PRIMARYID, CASEID, narrative text, or re-identification key will be released. These summaries describe the composition and completeness of reported cases; they do not estimate incidence, establish causality, or rank clinical treatment risk.

## Reporting and supplementary materials

The final manuscript will include a PRISMA-style search-flow diagram, a complete study-characteristics/evidence-boundary table, the search strategy and screening log as supplementary material, and a completed READUS-PV checklist. The checklist will be mapped to the title, abstract, methods, results, discussion, declarations, data/code availability, and protocol-registration fields where applicable.

## Verified anchor sources (not yet extracted or appraised)

- FDA Sentinel active-comparator new-user cohort on severe AUB: https://sentinelinitiative.org/news-events/publications-presentations/risk-severe-abnormal-uterine-bleeding-associated-rivaroxaban
- 2026 JTH systematic review/meta-analysis, PMID 41791665: https://pubmed.ncbi.nlm.nih.gov/41791665/
- 2025 prospective HEMBLED registry, PMID 41224211: https://pubmed.ncbi.nlm.nih.gov/41224211/
- 2026 MarketScan case-control study, PMID 41505758: https://pubmed.ncbi.nlm.nih.gov/41505758/
- PRISMA 2020: https://www.prisma-statement.org/prisma-2020-statement

These sources establish a starting inventory, not the direction or certainty of a final conclusion.
