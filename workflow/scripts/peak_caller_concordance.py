#!/usr/bin/env python3
"""Summarize promoter-level concordance between matched-depth peak callers.

The Jaccard index is calculated on promoter sets with any peak overlap, not on
base-pair or peak-interval sets.
"""
from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def genes_overlapping(promoters: Path, peaks: Path) -> set[str]:
    output = subprocess.check_output(
        ["bedtools", "intersect", "-a", str(promoters), "-b", str(peaks), "-u"], text=True
    )
    return {line.split("\t")[3] for line in output.splitlines() if line}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--promoters", type=Path, required=True)
    parser.add_argument("--matched-depth-dir", type=Path, required=True)
    parser.add_argument("--conditions", nargs="+", required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    total = sum(1 for _ in args.promoters.open())
    rows = []
    for condition in args.conditions:
        genrich = genes_overlapping(args.promoters, args.matched_depth_dir / f"{condition}.genrich.narrowPeak")
        macs3 = genes_overlapping(args.promoters, args.matched_depth_dir / f"{condition}.macs3.narrowPeak")
        both = genrich & macs3
        either = genrich | macs3
        rows.append({
            "condition": condition,
            "n_promoters": total,
            "genrich_open": len(genrich),
            "macs3_open": len(macs3),
            "both_open": len(both),
            "either_open": len(either),
            "jaccard": f"{len(both) / max(len(either), 1):.6f}",
            "binary_concordance": f"{(len(both) + total - len(either)) / total:.6f}",
        })
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
