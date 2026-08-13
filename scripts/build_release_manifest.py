#!/usr/bin/env python3
"""Hash release files for an auditable source/output manifest.

The manifest includes tracked source/output files plus selected ignored derived
workflow products when they are present locally. The ignored products are not
committed to git because of size, but they are part of the archival evidence for
the BMC revision.
"""
from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "release_manifest.sha256"
OPTIONAL_RELEASE_GLOBS = [
    "figures/output/*.png",
    "workflow/logs/*.vm.log",
    "workflow/results/bigwig/*.bw",
    "workflow/results/peaks/primary/*.narrowPeak",
    "workflow/results/peaks/replicate/*.narrowPeak",
    "workflow/results/peaks/pooled/*.narrowPeak",
    "workflow/results/matched_depth/*.narrowPeak",
    "workflow/results/matched_depth/*.tsv",
    "workflow/results/sensitivity/*.tsv",
    "workflow/results/sensitivity/*.csv",
    "workflow/results/sensitivity/tss/*.tsv",
    "workflow/results/shared_peak_signal.tsv",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    listing = subprocess.check_output(
        ["git", "ls-files", "--cached", "--others", "--exclude-standard", "-z"],
        cwd=ROOT,
    )
    paths = []
    for raw in listing.split(b"\0"):
        if not raw:
            continue
        relative = Path(raw.decode())
        path = ROOT / relative
        if path == OUTPUT or not path.is_file():
            continue
        paths.append(relative)
    for pattern in OPTIONAL_RELEASE_GLOBS:
        for path in ROOT.glob(pattern):
            if path.is_file() and path != OUTPUT:
                paths.append(path.relative_to(ROOT))
    paths = sorted(set(paths))
    lines = [f"{sha256(ROOT / relative)}  {relative.as_posix()}" for relative in paths]
    OUTPUT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUTPUT} with {len(lines)} file hashes")


if __name__ == "__main__":
    main()
