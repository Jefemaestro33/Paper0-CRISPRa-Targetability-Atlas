#!/usr/bin/env python3
"""Download and archive compressed public reference inputs for release.

The primary Snakemake workflow downloads decompressed mm39 FASTA and GENCODE
vM33 GTF resources. This helper keeps the final archival release lightweight and
auditable by packaging the original compressed public reference files plus a
tracked URL/size/SHA-256 audit table.
"""
from __future__ import annotations

import argparse
import hashlib
import shutil
import tarfile
import urllib.request
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_CACHE = ROOT / "release_archives" / "reference_inputs"
DEFAULT_OUTDIR = ROOT / "release_archives"
AUDIT_PATH = ROOT / "reference" / "source_input_audit.tsv"
REFERENCE_INPUTS = [
    {
        "name": "mm39.fa.gz",
        "resource": "GRCm39/mm39 FASTA",
        "url": "https://hgdownload.soe.ucsc.edu/goldenPath/mm39/bigZips/mm39.fa.gz",
        "workflow_output": "workflow/resources/mm39.fa",
    },
    {
        "name": "gencode.vM33.annotation.gtf.gz",
        "resource": "GENCODE mouse vM33 annotation",
        "url": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M33/gencode.vM33.annotation.gtf.gz",
        "workflow_output": "workflow/resources/gencode.vM33.annotation.gtf",
    },
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def download(url: str, path: Path) -> None:
    if path.exists() and path.stat().st_size > 0:
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".part")
    with urllib.request.urlopen(url) as response, tmp.open("wb") as handle:
        shutil.copyfileobj(response, handle)
    tmp.replace(path)


def write_audit(rows: list[dict[str, str]]) -> None:
    fields = [
        "resource",
        "archive_file",
        "source_url",
        "bytes",
        "sha256",
        "workflow_output",
    ]
    lines = ["\t".join(fields)]
    for row in rows:
        lines.append("\t".join(row[field] for field in fields))
    AUDIT_PATH.write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--prefix", default="paper0_reference_inputs")
    args = parser.parse_args()

    args.cache_dir.mkdir(parents=True, exist_ok=True)
    args.outdir.mkdir(parents=True, exist_ok=True)

    rows = []
    files = []
    for item in REFERENCE_INPUTS:
        path = args.cache_dir / item["name"]
        download(item["url"], path)
        files.append(path)
        rows.append(
            {
                "resource": item["resource"],
                "archive_file": path.name,
                "source_url": item["url"],
                "bytes": str(path.stat().st_size),
                "sha256": sha256(path),
                "workflow_output": item["workflow_output"],
            }
        )

    write_audit(rows)

    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    package_name = f"{args.prefix}_{stamp}"
    tar_path = args.outdir / f"{package_name}.tar.gz"
    sha_path = args.outdir / f"{tar_path.name}.sha256"
    readme = [
        f"Package: {package_name}",
        f"Created UTC: {stamp}",
        "",
        "Contents:",
        "- compressed public reference inputs used by the Snakemake workflow",
        "- SOURCE_INPUT_AUDIT.tsv with source URLs, byte counts, and SHA-256 hashes",
        "",
        "The workflow decompresses these files to the paths listed in the audit table.",
        "Raw FASTQs are not included; their ENA URLs and MD5 checksums are in config/samples.tsv.",
    ]

    with tarfile.open(tar_path, "w:gz") as archive:
        for path in files:
            archive.add(path, arcname=f"{package_name}/{path.name}")
        archive.add(AUDIT_PATH, arcname=f"{package_name}/SOURCE_INPUT_AUDIT.tsv")
        readme_bytes = ("\n".join(readme) + "\n").encode()
        info = tarfile.TarInfo(f"{package_name}/README.txt")
        info.size = len(readme_bytes)
        archive.addfile(info, fileobj=__import__("io").BytesIO(readme_bytes))

    digest = sha256(tar_path)
    sha_path.write_text(f"{digest}  {tar_path.name}\n")
    print(f"Wrote {AUDIT_PATH}")
    print(f"Wrote {tar_path}")
    print(f"Wrote {sha_path}")
    print(f"sha256={digest}")


if __name__ == "__main__":
    main()
