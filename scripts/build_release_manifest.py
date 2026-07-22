#!/usr/bin/env python3
"""Hash every non-ignored release file for an auditable source/output manifest."""
from __future__ import annotations

import hashlib
import subprocess
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "release_manifest.sha256"


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
    lines = [f"{sha256(ROOT / relative)}  {relative.as_posix()}" for relative in sorted(paths)]
    OUTPUT.write_text("\n".join(lines) + "\n")
    print(f"Wrote {OUTPUT} with {len(lines)} file hashes")


if __name__ == "__main__":
    main()
