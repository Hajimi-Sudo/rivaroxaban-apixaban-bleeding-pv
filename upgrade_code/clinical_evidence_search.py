"""Reproducible PubMed search and screening scaffold for the clinical layer."""

from __future__ import annotations

import argparse
import csv
import json
import re
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from datetime import datetime, timezone
from pathlib import Path


QUERY = (
    '(rivaroxaban[Title/Abstract] OR apixaban[Title/Abstract] OR '
    '"direct oral anticoagulant"[Title/Abstract]) AND '
    '("abnormal uterine bleeding"[Title/Abstract] OR '
    '"heavy menstrual bleeding"[Title/Abstract] OR menorrhagia[Title/Abstract] OR '
    '"uterine bleeding"[Title/Abstract]) AND '
    '(female[Title/Abstract] OR women[Title/Abstract] OR menstru*[Title/Abstract])'
)


def fetch_json(url: str) -> dict:
    with urllib.request.urlopen(url, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_xml(url: str) -> ET.Element:
    with urllib.request.urlopen(url, timeout=60) as response:
        return ET.fromstring(response.read())


def text(node: ET.Element | None) -> str:
    return " ".join("".join(node.itertext()).split()) if node is not None else ""


def parse_article(article: ET.Element) -> dict[str, str]:
    pmid = text(article.find(".//PMID"))
    title = text(article.find(".//ArticleTitle"))
    abstract = " ".join(text(n) for n in article.findall(".//Abstract/AbstractText"))
    journal = text(article.find(".//Journal/Title"))
    year = text(article.find(".//PubDate/Year")) or text(article.find(".//PubDate/MedlineDate"))
    doi = ""
    for aid in article.findall(".//ArticleId"):
        if aid.attrib.get("IdType") == "doi":
            doi = text(aid)
    types = "; ".join(text(n) for n in article.findall(".//PublicationType"))
    keywords = "; ".join(text(n) for n in article.findall(".//Keyword"))
    return {"pmid": pmid, "title": title, "abstract": abstract, "journal": journal, "pub_year": year, "doi": doi, "publication_types": types, "keywords": keywords}


def run(output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    run_date = datetime.now(timezone.utc).date().isoformat()
    params = urllib.parse.urlencode({"db": "pubmed", "retmode": "json", "retmax": "500", "sort": "pub date", "term": QUERY})
    search_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?{params}"
    search = fetch_json(search_url)
    result = search["esearchresult"]
    ids = result.get("idlist", [])
    fetch_params = urllib.parse.urlencode({"db": "pubmed", "retmode": "xml", "id": ",".join(ids)})
    fetch_url = f"https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?{fetch_params}"
    root = fetch_xml(fetch_url)
    records = [parse_article(node) for node in root.findall(".//PubmedArticle")]
    records.sort(key=lambda row: (row["pub_year"], row["pmid"]), reverse=True)

    with (output_dir / "pubmed_hits.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]) if records else ["pmid"])
        writer.writeheader()
        writer.writerows(records)

    post_cutoff = [row for row in records if re.search(r"20(25|26)", row["pub_year"])]
    with (output_dir / "pubmed_hits_post_2024.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(records[0]) if records else ["pmid"])
        writer.writeheader()
        writer.writerows(post_cutoff)

    screening = []
    for row in records:
        screening.append({"pmid": row["pmid"], "title": row["title"], "year": row["pub_year"], "screening_stage": "title_abstract_pending", "include": "", "exclusion_reason": "", "direct_rivaroxaban_apixaban_comparison": "", "exposure_denominator": "", "endpoint_class": "", "notes": ""})
    with (output_dir / "screening_log.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(screening[0]) if screening else ["pmid"])
        writer.writeheader()
        writer.writerows(screening)

    protocol = {
        "database": "PubMed/MEDLINE via NCBI E-utilities",
        "run_date_utc": datetime.now(timezone.utc).isoformat(),
        "publication_date_window": "inception through 2026-10-18",
        "query": QUERY,
        "search_url": search_url,
        "record_count": int(result.get("count", 0)),
        "retrieved_count": len(records),
        "post_2024_count": len(post_cutoff),
        "cutoff_rationale": "The existing 2026 JTH review searched through June 2024; post-cutoff records are prioritized for incremental evidence while all retrieved records remain in the screening log.",
        "eligibility": "Human studies with an exposure denominator and extractable uterine, vaginal, menstrual, or abnormal-uterine bleeding outcomes among women receiving rivaroxaban, apixaban, or another direct oral anticoagulant.",
        "exclusions": "Non-human, male-only, pediatric-only, non-bleeding, pharmacokinetic-only, editorial/commentary, duplicate-cohort, or no-denominator reports; exclusions are recorded after title/abstract and full-text screening.",
        "synthesis": "Endpoint-by-design evidence map; no pooling with spontaneous-report RORs and no pooled clinical estimate unless population, exposure, endpoint, estimand, and cohort overlap are directly comparable.",
    }
    (output_dir / "search_protocol.json").write_text(json.dumps(protocol, indent=2), encoding="utf-8")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    run(args.output_dir.resolve())


if __name__ == "__main__":
    main()
