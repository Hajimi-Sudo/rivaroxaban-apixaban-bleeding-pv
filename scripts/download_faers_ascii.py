"""Download and verify official FDA FAERS/AEMS quarterly ASCII archives."""

from __future__ import annotations

import csv
import hashlib
import http.cookiejar
import re
import shutil
import socket
import subprocess
import time
import urllib.error
import urllib.request
import zipfile
from concurrent.futures import ThreadPoolExecutor, as_completed
import os
from pathlib import Path


PROJECT_ROOT = Path(os.environ.get("PV_PROJECT_ROOT", Path.cwd())).resolve()
INDEX_URL = "https://fis.fda.gov/extensions/FPD-QDE-FAERS/FPD-QDE-FAERS.html"
START_YEAR = 2015
END_YEAR = 2026
WORKERS = 4
ATTEMPTS = 10

FAERS_ROOT = PROJECT_ROOT / "data" / "raw" / "faers"
ARCHIVES = FAERS_ROOT / "archives"
MANIFEST = FAERS_ROOT / "manifest"
ARCHIVES.mkdir(parents=True, exist_ok=True)
MANIFEST.mkdir(parents=True, exist_ok=True)


def make_opener() -> urllib.request.OpenerDirector:
    jar = http.cookiejar.CookieJar()
    return urllib.request.build_opener(urllib.request.HTTPCookieProcessor(jar))


def establish_session() -> tuple[urllib.request.OpenerDirector, str]:
    opener = make_opener()
    request = urllib.request.Request(INDEX_URL, headers={"User-Agent": "Mozilla/5.0"})
    with opener.open(request, timeout=30) as response:
        html = response.read().decode("utf-8", errors="replace")
    return opener, html


def discover_records() -> list[dict[str, object]]:
    _, html = establish_session()
    pattern = re.compile(
        r'href=["\'](https://fis\.fda\.gov/content/Exports/faers_ascii_(20\d{2})[qQ]([1-4])\.zip)["\']',
        re.IGNORECASE,
    )
    unique: dict[str, dict[str, object]] = {}
    for url, year_text, quarter_text in pattern.findall(html):
        year, quarter = int(year_text), int(quarter_text)
        if START_YEAR <= year <= END_YEAR:
            period = f"{year}Q{quarter}"
            unique[period] = {
                "period": period,
                "year": year,
                "quarter": quarter,
                "url": url,
                "filename": url.rsplit("/", 1)[-1],
            }
    records = sorted(unique.values(), key=lambda item: (item["year"], item["quarter"]))
    if len(records) != 45:
        raise RuntimeError(f"Expected 45 quarters from 2015Q1 through 2026Q1; found {len(records)}")
    return records


def valid_zip(path: Path) -> bool:
    if not path.exists() or path.stat().st_size == 0:
        return False
    try:
        with zipfile.ZipFile(path) as archive:
            return bool(archive.infolist()) and archive.testzip() is None
    except (OSError, zipfile.BadZipFile):
        return False


def download_one(record: dict[str, object]) -> tuple[str, Path, str]:
    period = str(record["period"])
    target = ARCHIVES / str(record["filename"])
    partial = Path(f"{target}.part")
    if valid_zip(target):
        return period, target, "existing"
    if target.exists():
        target.unlink()

    cookie_file = MANIFEST / f"cookie_{period}.txt"
    for attempt in range(1, ATTEMPTS + 1):
        try:
            subprocess.run(
                [
                    "curl.exe", "--silent", "--show-error", "--fail", "--location",
                    "--cookie-jar", str(cookie_file), INDEX_URL,
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                check=True,
                timeout=45,
            )
            subprocess.run(
                [
                    "curl.exe", "--silent", "--show-error", "--fail", "--location",
                    "--connect-timeout", "30", "--max-time", "180",
                    "--speed-limit", "1024", "--speed-time", "30",
                    "--cookie", str(cookie_file), "--continue-at", "-",
                    "--output", str(partial), str(record["url"]),
                ],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.PIPE,
                check=True,
                timeout=200,
            )

            if valid_zip(partial):
                partial.replace(target)
                cookie_file.unlink(missing_ok=True)
                return period, target, "downloaded"
            if partial.exists() and partial.stat().st_size == 0:
                partial.unlink()
            raise zipfile.BadZipFile("response did not form a valid ZIP archive")
        except (
            OSError,
            socket.timeout,
            subprocess.CalledProcessError,
            subprocess.TimeoutExpired,
            urllib.error.URLError,
            zipfile.BadZipFile,
        ) as exc:
            if attempt == ATTEMPTS:
                raise RuntimeError(f"{period} failed after {ATTEMPTS} attempts: {exc}") from exc
            time.sleep(min(3 * attempt, 15))

    raise AssertionError("unreachable")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(4 * 1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest().upper()


def main() -> None:
    records = discover_records()
    with (MANIFEST / "official_source_urls.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=["period", "year", "quarter", "url", "filename"])
        writer.writeheader()
        writer.writerows(records)

    completed = 0
    with ThreadPoolExecutor(max_workers=WORKERS) as pool:
        futures = {pool.submit(download_one, record): record for record in records}
        for future in as_completed(futures):
            period, path, state = future.result()
            completed += 1
            print(f"[{completed:02d}/45] {period} {state}: {path.name} ({path.stat().st_size / 1024**2:.1f} MB)", flush=True)

    rows = []
    for record in records:
        path = ARCHIVES / str(record["filename"])
        if not valid_zip(path):
            raise RuntimeError(f"Final validation failed: {path}")
        with zipfile.ZipFile(path) as archive:
            entries = len(archive.infolist())
        rows.append(
            {
                "Period": record["period"],
                "FileName": path.name,
                "Bytes": path.stat().st_size,
                "SHA256": sha256(path),
                "ZipEntries": entries,
                "SourceUrl": record["url"],
            }
        )

    with (MANIFEST / "download_manifest_sha256.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)
    print("All 45 quarterly archives downloaded and verified.", flush=True)


if __name__ == "__main__":
    main()
