#!/usr/bin/env python3
"""Lift the official ENCODE mm10 blacklist to mm39 with full provenance."""
from __future__ import annotations

import argparse
import csv
import hashlib
import subprocess
import tempfile
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def chromosome_key(chromosome: str) -> tuple[int, int | str]:
    value = chromosome.removeprefix("chr")
    if value.isdigit():
        return 0, int(value)
    order = {"X": 20, "Y": 21, "M": 22}
    return (0, order[value]) if value in order else (1, value)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-mm10", type=Path, required=True)
    parser.add_argument("--chain", type=Path, required=True)
    parser.add_argument("--liftover", type=Path, required=True)
    parser.add_argument("--output-mm39", type=Path, required=True)
    parser.add_argument("--provenance", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    source = []
    with args.input_mm10.open() as handle:
        for index, line in enumerate(handle, start=1):
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            start, end = int(fields[1]), int(fields[2])
            if end <= start:
                raise ValueError(f"Invalid source interval on line {index}: {line.rstrip()}")
            source.append({
                "id": f"BL{index:06d}", "source_chrom": fields[0],
                "source_start": start, "source_end": end,
                "label": fields[3] if len(fields) > 3 else "ENCODE blacklist v2",
            })

    by_id = {row["id"]: row for row in source}
    with tempfile.TemporaryDirectory(prefix="paper0_blacklist_") as directory:
        directory = Path(directory)
        safe_input = directory / "mm10.safe.bed"
        mapped_path = directory / "mm39.mapped.bed"
        unmapped_path = directory / "mm39.unmapped.bed"
        with safe_input.open("w") as handle:
            for row in source:
                handle.write(
                    f"{row['source_chrom']}\t{row['source_start']}\t{row['source_end']}\t{row['id']}\n"
                )
        subprocess.run([
            str(args.liftover), str(safe_input), str(args.chain),
            str(mapped_path), str(unmapped_path),
        ], check=True)

        mapped = []
        with mapped_path.open() as handle:
            for line in handle:
                chrom, start_text, end_text, identifier = line.rstrip("\n").split("\t")[:4]
                start, end = int(start_text), int(end_text)
                if identifier not in by_id:
                    raise ValueError(f"Unknown liftOver identifier: {identifier}")
                if end <= start:
                    raise ValueError(f"Invalid lifted interval: {line.rstrip()}")
                mapped.append({
                    **by_id[identifier], "mapped_chrom": chrom,
                    "mapped_start": start, "mapped_end": end, "status": "mapped",
                })

    mapped_ids = {row["id"] for row in mapped}
    provenance = mapped + [
        {**row, "mapped_chrom": "", "mapped_start": "", "mapped_end": "", "status": "unmapped"}
        for row in source if row["id"] not in mapped_ids
    ]
    mapped.sort(key=lambda row: (chromosome_key(row["mapped_chrom"]), row["mapped_start"], row["mapped_end"]))

    args.output_mm39.parent.mkdir(parents=True, exist_ok=True)
    with args.output_mm39.open("w") as handle:
        for row in mapped:
            handle.write(
                f"{row['mapped_chrom']}\t{row['mapped_start']}\t{row['mapped_end']}\t{row['label']}\n"
            )

    fields = [
        "id", "source_chrom", "source_start", "source_end", "label",
        "mapped_chrom", "mapped_start", "mapped_end", "status",
    ]
    args.provenance.parent.mkdir(parents=True, exist_ok=True)
    with args.provenance.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(provenance, key=lambda row: row["id"]))

    print(
        f"Official source SHA-256: {sha256(args.input_mm10)}\n"
        f"Mapped {len(mapped):,}/{len(source):,}; unmapped {len(source)-len(mapped):,}; "
        f"mm39 SHA-256: {sha256(args.output_mm39)}"
    )


if __name__ == "__main__":
    main()
