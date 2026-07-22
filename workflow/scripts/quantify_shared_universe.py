#!/usr/bin/env python3
"""Quantify normalized ATAC signal on a shared union of primary peaks."""
from __future__ import annotations

import argparse
import csv
import subprocess
import tempfile
from pathlib import Path


def fragment_bed(bam: Path, layout: str, directory: Path, condition: str) -> tuple[Path, int]:
    """Write one BED interval per paired fragment or single-end alignment."""
    output = directory / f"{condition}.fragments.bed"
    total = 0
    if layout == "PE":
        name_sorted = directory / f"{condition}.name_sorted.bam"
        subprocess.run(["samtools", "sort", "-n", "-o", str(name_sorted), str(bam)], check=True)
        process = subprocess.Popen(
            ["bedtools", "bamtobed", "-bedpe", "-i", str(name_sorted)],
            stdout=subprocess.PIPE, text=True,
        )
        assert process.stdout is not None
        with output.open("w") as handle:
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6 or fields[0] != fields[3]:
                    continue
                start = min(int(fields[1]), int(fields[4]))
                end = max(int(fields[2]), int(fields[5]))
                if end > start:
                    handle.write(f"{fields[0]}\t{start}\t{end}\n")
                    total += 1
        if process.wait() != 0:
            raise subprocess.CalledProcessError(process.returncode, process.args)
    else:
        process = subprocess.Popen(
            ["bedtools", "bamtobed", "-i", str(bam)],
            stdout=subprocess.PIPE, text=True,
        )
        assert process.stdout is not None
        with output.open("w") as handle:
            for line in process.stdout:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 3:
                    handle.write("\t".join(fields[:3]) + "\n")
                    total += 1
        if process.wait() != 0:
            raise subprocess.CalledProcessError(process.returncode, process.args)
    return output, total


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--primary-peaks", nargs="+", type=Path, required=True)
    parser.add_argument("--bam-manifest", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()

    with tempfile.TemporaryDirectory() as tmpdir:
        combined = Path(tmpdir) / "combined.bed"
        with combined.open("w") as out:
            for peak_path in args.primary_peaks:
                with peak_path.open() as handle:
                    for line in handle:
                        if line.strip() and not line.startswith(("#", "track")):
                            fields = line.split("\t")
                            out.write("\t".join(fields[:3]) + "\n")
        universe = Path(tmpdir) / "universe.bed"
        sort_proc = subprocess.Popen(["bedtools", "sort", "-i", str(combined)], stdout=subprocess.PIPE)
        with universe.open("w") as out:
            subprocess.run(["bedtools", "merge", "-i", "stdin"], stdin=sort_proc.stdout, stdout=out, check=True)
        sort_proc.wait()

        regions = []
        with universe.open() as handle:
            for index, line in enumerate(handle, start=1):
                chrom, start, end = line.rstrip("\n").split("\t")[:3]
                regions.append({"peak_universe_id": f"U{index}", "chrom": chrom, "start": start, "end": end})

        with args.bam_manifest.open(newline="") as handle:
            bams = list(csv.DictReader(handle, delimiter="\t"))

        for row in bams:
            condition, bam = row["condition"], Path(row["bam"])
            fragments, total = fragment_bed(bam, row["layout"], Path(tmpdir), condition)
            coverage = subprocess.check_output(
                ["bedtools", "coverage", "-a", str(universe), "-b", str(fragments), "-counts"], text=True
            )
            for region, line in zip(regions, coverage.splitlines(), strict=True):
                count = int(line.rsplit("\t", 1)[1])
                region[f"fragments_{condition}"] = count
                region[f"cpm_{condition}"] = f"{count * 1_000_000 / max(total, 1):.6f}"

    args.output.parent.mkdir(parents=True, exist_ok=True)
    fields = list(regions[0].keys())
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(regions)
    print(f"Quantified {len(regions):,} shared peak intervals across {len(bams)} conditions")


if __name__ == "__main__":
    main()
