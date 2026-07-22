#!/usr/bin/env python3
from __future__ import annotations

import argparse
import csv
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.sample_manifest.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    condition_layout = {}
    for row in rows:
        condition_layout.setdefault(row["condition"], row["layout"])
        if condition_layout[row["condition"]] != row["layout"]:
            raise SystemExit(f"Mixed layouts in condition {row['condition']}")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["condition", "layout", "bam"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for condition in sorted(condition_layout):
            writer.writerow({
                "condition": condition,
                "layout": condition_layout[condition],
                "bam": str((args.bam_dir / f"{condition}.bam").resolve()),
            })


if __name__ == "__main__":
    main()
