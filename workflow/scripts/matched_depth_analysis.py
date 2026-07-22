#!/usr/bin/env python3
"""Generate replicate-aware deterministic matched-depth peak sensitivities."""
from __future__ import annotations

import argparse
import csv
import shutil
import subprocess
from pathlib import Path


def run(command: list[str], stdout=None) -> None:
    subprocess.run(command, check=True, stdout=stdout)


def alignment_count(bam: Path, paired: bool) -> int:
    command = ["samtools", "view", "-c", "-F", "3852"]
    if paired:
        command.extend(["-f", "64"])
    command.append(str(bam))
    return int(subprocess.check_output(command, text=True).strip())


def allocate_target(capacities: list[int], target: int) -> list[int]:
    """Allocate an exact requested total as evenly as capacities permit."""
    if target > sum(capacities):
        raise ValueError("Requested target exceeds available fragments")
    allocation = [0] * len(capacities)
    remaining = target
    active = list(range(len(capacities)))
    while remaining and active:
        share, extra = divmod(remaining, len(active))
        saturated = [index for index in active if capacities[index] - allocation[index] <= share]
        if saturated:
            for index in saturated:
                available = capacities[index] - allocation[index]
                allocation[index] += available
                remaining -= available
                active.remove(index)
            continue
        for position, index in enumerate(active):
            amount = share + (position < extra)
            allocation[index] += amount
            remaining -= amount
        break
    assert sum(allocation) == target
    assert all(value <= capacity for value, capacity in zip(allocation, capacities, strict=True))
    return allocation


def subsample_bam(source: Path, output: Path, requested: int, available: int, seed: int) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    output.unlink(missing_ok=True)
    Path(str(output) + ".bai").unlink(missing_ok=True)
    fraction = min(1.0, requested / available)
    command = ["samtools", "view", "-b"]
    if fraction < 1:
        digits = f"{fraction:.12f}".split(".", 1)[1]
        command.extend(["-s", f"{seed}.{digits}"])
    command.extend(["-o", str(output), str(source)])
    run(command)
    run(["samtools", "index", str(output)])


def call_genrich(bam: Path, output: Path, layout: str) -> None:
    output.unlink(missing_ok=True)
    name_sorted = output.with_suffix(".name_sorted.bam")
    name_sorted.unlink(missing_ok=True)
    run(["samtools", "sort", "-n", "-o", str(name_sorted), str(bam)])
    command = ["Genrich", "-t", str(name_sorted), "-o", str(output), "-j"]
    if layout == "SE":
        command.append("-y")
    run(command)
    name_sorted.unlink()


def call_macs3(bam: Path, output: Path, layout: str) -> None:
    output.unlink(missing_ok=True)
    prefix = output.name.removesuffix(".narrowPeak")
    command = [
        "macs3", "callpeak", "-t", str(bam),
        "-f", "BAMPE" if layout == "PE" else "BAM",
        "-g", "mm", "-n", prefix, "--outdir", str(output.parent),
        "--keep-dup", "all", "-q", "0.01",
    ]
    if layout == "SE":
        command.extend(["--nomodel", "--shift", "-100", "--extsize", "200"])
    run(command)
    produced = output.parent / f"{prefix}_peaks.narrowPeak"
    shutil.move(produced, output)


def consensus(inputs: list[Path], output: Path, script: Path) -> None:
    output.unlink(missing_ok=True)
    run([
        "python", str(script), "--inputs", *map(str, inputs),
        "--output", str(output), "--min-replicates", "2",
        "--reciprocal-overlap", "0.50",
    ])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--condition-manifest", type=Path, required=True)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--run-bam-dir", type=Path, required=True)
    parser.add_argument("--consensus-script", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--seed", type=int, default=1729)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    with args.condition_manifest.open(newline="") as handle:
        condition_rows = list(csv.DictReader(handle, delimiter="\t"))
    with args.sample_manifest.open(newline="") as handle:
        sample_rows = list(csv.DictReader(handle, delimiter="\t"))

    condition_counts = {
        row["condition"]: alignment_count(Path(row["bam"]), row["layout"] == "PE")
        for row in condition_rows
    }
    target = min(condition_counts.values())
    summary = []
    for condition_offset, condition_row in enumerate(condition_rows):
        condition = condition_row["condition"]
        layout = condition_row["layout"]
        relevant = [row for row in sample_rows if row["condition"] == condition]
        biological = all(row["replicate_type"] == "biological" for row in relevant)
        if biological:
            source_bams = [args.run_bam_dir / f"{row['run_accession']}.bam" for row in relevant]
            source_labels = [row["run_accession"] for row in relevant]
            capacities = [alignment_count(path, layout == "PE") for path in source_bams]
            requested = allocate_target(capacities, target)
            matched_rule = "replicate_level_downsampling_then_replicate_consensus"
        else:
            assert all(row["replicate_type"] == "technical" for row in relevant)
            source_bams = [Path(condition_row["bam"])]
            source_labels = [f"{condition}_technical_pool"]
            capacities = [condition_counts[condition]]
            requested = [target]
            matched_rule = "cross_lane_deduplicated_technical_pool_downsampling"

        genrich_replicates: list[Path] = []
        macs3_replicates: list[Path] = []
        actual_total = 0
        seed_values = []
        for source_offset, (source, label, capacity, requested_n) in enumerate(
            zip(source_bams, source_labels, capacities, requested, strict=True)
        ):
            seed = args.seed + condition_offset * 100 + source_offset
            seed_values.append(seed)
            downsampled = args.outdir / f"{condition}.{label}.matched.bam"
            subsample_bam(source, downsampled, requested_n, capacity, seed)
            actual_total += alignment_count(downsampled, layout == "PE")
            genrich_peak = args.outdir / f"{condition}.{label}.genrich.replicate.narrowPeak"
            macs3_peak = args.outdir / f"{condition}.{label}.macs3.replicate.narrowPeak"
            call_genrich(downsampled, genrich_peak, layout)
            call_macs3(downsampled, macs3_peak, layout)
            genrich_replicates.append(genrich_peak)
            macs3_replicates.append(macs3_peak)

        genrich_final = args.outdir / f"{condition}.genrich.narrowPeak"
        macs3_final = args.outdir / f"{condition}.macs3.narrowPeak"
        if biological:
            consensus(genrich_replicates, genrich_final, args.consensus_script)
            consensus(macs3_replicates, macs3_final, args.consensus_script)
        else:
            shutil.copy2(genrich_replicates[0], genrich_final)
            shutil.copy2(macs3_replicates[0], macs3_final)

        summary.append({
            "condition": condition,
            "original_fragments_or_reads": condition_counts[condition],
            "requested_target_fragments_or_reads": target,
            "actual_matched_fragments_or_reads": actual_total,
            "sampling_fraction": f"{target / condition_counts[condition]:.8f}",
            "seeds": ";".join(map(str, seed_values)),
            "matched_depth_rule": matched_rule,
            "n_source_units": len(source_bams),
            "genrich_peaks": sum(1 for _ in genrich_final.open()),
            "macs3_peaks": sum(1 for _ in macs3_final.open()),
        })

    with (args.outdir / "matched_depth_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary[0].keys(), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(summary)


if __name__ == "__main__":
    main()
