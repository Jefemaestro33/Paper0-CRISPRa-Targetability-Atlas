#!/usr/bin/env python3
"""Create a local GitHub/Zenodo-style archival package.

The package contains the tracked repository state plus optional derived signal
tracks, peak files, sensitivity outputs, PNG renderings, and VM logs if they are
present locally. It intentionally avoids private editorial folders and does not
attempt to create a remote DOI. The output tarball can be uploaded to a
persistent archive after final approval.
"""
from __future__ import annotations

import argparse
import hashlib
import subprocess
import tarfile
from datetime import datetime, timezone
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUTDIR = ROOT / "release_archives"
BIGWIGS = [
    ROOT / "workflow/results/bigwig/homeostatic.bw",
    ROOT / "workflow/results/bigwig/PP_control.bw",
    ROOT / "workflow/results/bigwig/PL_acute_LPS.bw",
    ROOT / "workflow/results/bigwig/LL_tolerized.bw",
    ROOT / "workflow/results/bigwig/sham_WT.bw",
    ROOT / "workflow/results/bigwig/stroke_WT.bw",
]
OPTIONAL_RELEASE_GLOBS = [
    "figures/output/*.png",
    "workflow/logs/*.vm.log",
    "workflow/results/bigwig/*.bw",
    "workflow/results/peaks/primary/*.narrowPeak",
    "workflow/results/peaks/primary/*.narrowPeak.gz",
    "workflow/results/peaks/replicate/*.narrowPeak",
    "workflow/results/peaks/replicate/*.narrowPeak.gz",
    "workflow/results/peaks/pooled/*.narrowPeak",
    "workflow/results/peaks/pooled/*.narrowPeak.gz",
    "workflow/results/matched_depth/*.narrowPeak",
    "workflow/results/matched_depth/*.narrowPeak.gz",
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


def git_output(args: list[str]) -> str:
    return subprocess.check_output(["git", *args], cwd=ROOT, text=True).strip()


def tracked_files() -> list[Path]:
    raw = subprocess.check_output(["git", "ls-files", "-z"], cwd=ROOT)
    return [Path(item.decode()) for item in raw.split(b"\0") if item]


def optional_release_files() -> list[Path]:
    paths: set[Path] = set()
    for pattern in OPTIONAL_RELEASE_GLOBS:
        for path in ROOT.glob(pattern):
            if path.is_file():
                paths.add(path)
    return sorted(paths)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, default=DEFAULT_OUTDIR)
    parser.add_argument("--prefix", default="paper0_bmc_revision")
    args = parser.parse_args()

    commit = git_output(["rev-parse", "--short", "HEAD"])
    dirty = bool(git_output(["status", "--porcelain"]))
    stamp = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    package_name = f"{args.prefix}_{commit}_{stamp}"
    outdir = args.outdir
    outdir.mkdir(parents=True, exist_ok=True)
    tar_path = outdir / f"{package_name}.tar.gz"
    sha_path = outdir / f"{tar_path.name}.sha256"

    files = tracked_files()
    optional_files = optional_release_files()
    present_bigwigs = [path for path in BIGWIGS if path.exists()]
    missing_bigwigs = [path for path in BIGWIGS if not path.exists()]

    readme_lines = [
        f"Package: {package_name}",
        f"Git commit: {commit}",
        f"Working tree dirty at package time: {dirty}",
        f"Created UTC: {stamp}",
        "",
        "Contents:",
        "- tracked repository files from git ls-files",
        "- optional ignored workflow/figure/log files if present locally",
        "",
        "BigWig status:",
    ]
    if present_bigwigs:
        for path in present_bigwigs:
            readme_lines.append(f"- included: {path.relative_to(ROOT)} ({path.stat().st_size} bytes)")
    if missing_bigwigs:
        for path in missing_bigwigs:
            readme_lines.append(f"- missing: {path.relative_to(ROOT)}")
        readme_lines.append("")
        readme_lines.append(
            "Missing bigWigs are derived visualization intermediates needed to rebuild Figure 6. "
            "They must be recreated by the Snakemake workflow or restored before a final archival release "
            "if the archive is intended to contain signal tracks."
        )

    with tarfile.open(tar_path, "w:gz") as archive:
        for relative in files:
            archive.add(ROOT / relative, arcname=f"{package_name}/{relative.as_posix()}")
        for path in optional_files:
            archive.add(path, arcname=f"{package_name}/{path.relative_to(ROOT).as_posix()}")
        readme_bytes = ("\n".join(readme_lines) + "\n").encode()
        info = tarfile.TarInfo(f"{package_name}/ARCHIVAL_PACKAGE_README.txt")
        info.size = len(readme_bytes)
        archive.addfile(info, fileobj=__import__("io").BytesIO(readme_bytes))

    digest = sha256(tar_path)
    sha_path.write_text(f"{digest}  {tar_path.name}\n")
    print(f"Wrote {tar_path}")
    print(f"Wrote {sha_path}")
    print(f"sha256={digest}")
    if missing_bigwigs:
        print("WARNING: bigWigs missing; see ARCHIVAL_PACKAGE_README.txt inside archive.")


if __name__ == "__main__":
    main()
