"""Run the existing audited FAERS engine against an isolated 46-quarter runtime.

This wrapper is intentionally a baseline update, not the final upgraded analysis.
It appends the verified 2026 Q2 archive to an upgrade-only manifest, uses a
runtime root with symlinked read-only archives, and writes outputs only below
the upgrade workspace. It never modifies the old project manifest or results.
"""

from __future__ import annotations

import argparse
import hashlib
import os
import sys
from pathlib import Path

import pandas as pd


Q2_NAME = "faers_ascii_2026q2.zip"
Q2_SHA256 = "90213500721ecaf976ee45d8ad04aa2dbb80861d93c87fef67ccaf515919bf2a"
Q2_URL = "https://fis.fda.gov/content/Exports/faers_ascii_2026q2.zip"


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(4 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def build_runtime(project_root: Path, upgrade_root: Path) -> tuple[Path, list[str]]:
    archive_dir = project_root / ".aris" / "data" / "faers" / "archives"
    source_manifest = project_root / ".aris" / "data" / "faers" / "manifest" / "download_manifest_sha256.csv"
    q2_archive = archive_dir / Q2_NAME
    if not q2_archive.exists():
        raise FileNotFoundError(q2_archive)
    if sha256(q2_archive) != Q2_SHA256:
        raise ValueError("2026 Q2 SHA-256 mismatch")

    upgrade_manifest_dir = upgrade_root / "runtime_manifest"
    upgrade_manifest_dir.mkdir(parents=True, exist_ok=True)
    manifest = pd.read_csv(source_manifest, dtype="string")
    if "2026Q2" in set(manifest["Period"].tolist()):
        raise ValueError("source manifest unexpectedly already contains 2026Q2")
    q2_row = {
        "Period": "2026Q2",
        "FileName": Q2_NAME,
        "Bytes": q2_archive.stat().st_size,
        "SHA256": Q2_SHA256.upper(),
        "ZipEntries": 18,
        "SourceUrl": Q2_URL,
    }
    manifest = pd.concat([manifest, pd.DataFrame([q2_row])], ignore_index=True)
    upgrade_manifest = upgrade_manifest_dir / "download_manifest_sha256.csv"
    manifest.to_csv(upgrade_manifest, index=False, encoding="utf-8-sig")

    runtime = upgrade_root / "runtime_project"
    runtime_manifest = runtime / ".aris" / "data" / "faers" / "manifest"
    runtime_archives = runtime / ".aris" / "data" / "faers" / "archives"
    runtime_manifest.mkdir(parents=True, exist_ok=True)
    runtime_archives.parent.mkdir(parents=True, exist_ok=True)
    link = runtime_archives
    if link.exists() or link.is_symlink():
        if link.is_symlink() or link.is_file():
            link.unlink()
        else:
            raise FileExistsError(link)
    link.symlink_to(archive_dir, target_is_directory=True)
    (runtime_manifest / "download_manifest_sha256.csv").write_bytes(upgrade_manifest.read_bytes())
    periods = sorted(manifest["Period"].tolist())
    return runtime, periods


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--project-root", type=Path, required=True)
    ap.add_argument("--upgrade-root", type=Path, required=True)
    ap.add_argument("--periods", nargs="+", default=None)
    args = ap.parse_args()
    project_root = args.project_root.resolve()
    upgrade_root = args.upgrade_root.resolve()
    runtime, periods = build_runtime(project_root, upgrade_root)
    sys.path.insert(0, str(project_root / "experiment" / "code"))
    from faers_full import run  # imported only after the isolated runtime is ready

    selected = args.periods or periods
    out = upgrade_root / "results" / ("faers_" + selected[0] + "_" + selected[-1] + "_baseline")
    result = run(runtime, out, periods=selected)
    print(result)


if __name__ == "__main__":
    main()
