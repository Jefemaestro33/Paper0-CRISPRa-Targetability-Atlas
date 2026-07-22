#!/usr/bin/env python3
"""Compute per-run ATAC-seq QC and primary-peak provenance metrics."""
from __future__ import annotations

import argparse
import bisect
import csv
import re
from collections import defaultdict
from pathlib import Path

import numpy as np
import pysam


class PeakIndex:
    def __init__(self, path: Path):
        intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
        with path.open() as handle:
            for line in handle:
                if line.strip() and not line.startswith(("#", "track")):
                    fields = line.split("\t")
                    intervals[fields[0]].append((int(fields[1]), int(fields[2])))
        self.intervals = {chrom: sorted(values) for chrom, values in intervals.items()}
        self.starts = {chrom: [start for start, _ in values] for chrom, values in self.intervals.items()}
        self.n = sum(len(values) for values in self.intervals.values())

    def overlaps(self, chrom: str, start: int, end: int) -> bool:
        values = self.intervals.get(chrom, [])
        starts = self.starts.get(chrom, [])
        if not values:
            return False
        index = bisect.bisect_right(starts, start) - 1
        if index >= 0 and values[index][1] > start:
            return True
        index += 1
        return index < len(values) and values[index][0] < end


def load_tss(path: Path) -> tuple[dict[str, list[int]], dict[tuple[str, int], str]]:
    positions: dict[str, list[int]] = defaultdict(list)
    strands = {}
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["definition"] == "ensembl_canonical":
                chrom, tss = row["chrom"], int(row["tss"])
                positions[chrom].append(tss)
                strands[(chrom, tss)] = row["strand"]
    return {chrom: sorted(values) for chrom, values in positions.items()}, strands


def alignment_rate(path: Path) -> float | None:
    if not path.exists():
        return None
    matches = re.findall(r"([0-9.]+)% overall alignment rate", path.read_text(errors="replace"))
    return float(matches[-1]) if matches else None


def duplication_rate(path: Path) -> float | None:
    if not path.exists():
        return None
    lines = path.read_text(errors="replace").splitlines()
    for index, line in enumerate(lines[:-1]):
        if line.startswith("LIBRARY\t") and "PERCENT_DUPLICATION" in line:
            fields = line.split("\t")
            values = lines[index + 1].split("\t")
            if len(values) == len(fields):
                return float(dict(zip(fields, values))["PERCENT_DUPLICATION"])
    return None


def run_qc(bam_path: Path, peak_path: Path, layout: str, tss_by_chrom, tss_strands) -> dict:
    peaks = PeakIndex(peak_path)
    bam = pysam.AlignmentFile(str(bam_path), "rb")
    usable_units = units_in_peaks = 0
    fragment_sizes = []
    tss_profile = np.zeros(4001, dtype=np.int64)
    insertion_count = 0
    for read in bam.fetch(until_eof=True):
        if read.is_unmapped or read.is_secondary or read.is_supplementary:
            continue
        chrom = bam.get_reference_name(read.reference_id)
        # Count fragments for paired-end runs and reads for single-end runs.
        if layout == "PE":
            if read.is_read1 and read.is_proper_pair:
                start = min(read.reference_start, read.next_reference_start)
                end = start + abs(read.template_length)
                if end <= start:
                    end = read.reference_end
                usable_units += 1
                fragment_sizes.append(end - start)
                units_in_peaks += peaks.overlaps(chrom, start, end)
        else:
            usable_units += 1
            units_in_peaks += peaks.overlaps(chrom, read.reference_start, read.reference_end)

        cut = read.reference_start + 4 if not read.is_reverse else read.reference_end - 5
        positions = tss_by_chrom.get(chrom, [])
        left = bisect.bisect_left(positions, cut - 2000)
        right = bisect.bisect_right(positions, cut + 2000)
        for tss in positions[left:right]:
            offset = cut - tss
            if tss_strands[(chrom, tss)] == "-":
                offset = -offset
            tss_profile[offset + 2000] += 1
        insertion_count += 1
    bam.close()
    flank = np.concatenate([tss_profile[:100], tss_profile[-100:]])
    background = max(float(flank.mean()), 1e-12)
    tss_enrichment = float(tss_profile[1950:2051].max() / background)
    sizes = np.asarray(fragment_sizes, dtype=float)
    return {
        "usable_fragments_or_reads": usable_units,
        "frip": units_in_peaks / usable_units if usable_units else float("nan"),
        "replicate_peak_count": peaks.n,
        "tss_enrichment_max": tss_enrichment,
        "insertions_for_tss_profile": insertion_count,
        "median_fragment_size": float(np.median(sizes)) if len(sizes) else "NA",
        "fraction_subnucleosomal_lt150": float((sizes < 150).mean()) if len(sizes) else "NA",
        "fraction_mononucleosomal_150_300": float(((sizes >= 150) & (sizes < 300)).mean()) if len(sizes) else "NA",
    }


def count_peaks(path: Path) -> int:
    with path.open() as handle:
        return sum(1 for line in handle if line.strip() and not line.startswith(("#", "track")))


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--bam-dir", type=Path, required=True)
    parser.add_argument("--replicate-peak-dir", type=Path, required=True)
    parser.add_argument("--primary-peak-dir", type=Path, required=True)
    parser.add_argument("--bowtie-log-dir", type=Path, required=True)
    parser.add_argument("--duplication-metrics-dir", type=Path)
    parser.add_argument("--tss-selection", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    duplication_dir = args.duplication_metrics_dir or args.bowtie_log_dir
    with args.sample_manifest.open(newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    tss, strands = load_tss(args.tss_selection)
    condition_counts = defaultdict(int)
    for row in samples:
        condition_counts[row["condition"]] += 1
    rows = []
    for sample in samples:
        run = sample["run_accession"]
        metrics = run_qc(
            args.bam_dir / f"{run}.bam",
            args.replicate_peak_dir / f"{run}.narrowPeak",
            sample["layout"], tss, strands,
        )
        rows.append({
            "level": "run", "run_accession": run, "condition": sample["condition"],
            "study": sample["study"], "replicate_type": sample["replicate_type"],
            "biological_replicate": sample["biological_replicate"],
            "technical_replicate": sample.get("technical_replicate", sample["biological_replicate"] if sample["replicate_type"] == "technical" else ""), "layout": sample["layout"],
            "bowtie2_overall_alignment_rate_pct": alignment_rate(args.bowtie_log_dir / f"{run}.bowtie2.log"),
            "picard_percent_duplication": duplication_rate(duplication_dir / f"{run}.duplication_metrics.txt"),
            "primary_peak_rule": "biological_replicate_consensus" if sample["replicate_type"] == "biological" else "technical_runs_pooled_no_biological_reproducibility_claim",
            **metrics,
        })
    for condition, count in condition_counts.items():
        technical_pool_metrics = duplication_dir / f"{condition}.technical_pool_duplication_metrics.txt"
        rows.append({
            "level": "condition", "run_accession": "", "condition": condition,
            "study": next(row["study"] for row in samples if row["condition"] == condition),
            "replicate_type": next(row["replicate_type"] for row in samples if row["condition"] == condition),
            "biological_replicate": "", "technical_replicate": "", "layout": "",
            "bowtie2_overall_alignment_rate_pct": "", "usable_fragments_or_reads": "",
            "picard_percent_duplication": duplication_rate(technical_pool_metrics) if technical_pool_metrics.exists() else "",
            "frip": "", "replicate_peak_count": "", "tss_enrichment_max": "",
            "insertions_for_tss_profile": "", "median_fragment_size": "",
            "fraction_subnucleosomal_lt150": "", "fraction_mononucleosomal_150_300": "",
            "primary_peak_rule": "support_in_at_least_2_biological_replicates_reciprocal_overlap_0.5" if condition in {"homeostatic", "sham_WT", "stroke_WT"} else "pooled_technical_runs_cross_run_deduplicated_no_biological_reproducibility_claim",
            "primary_peak_count": count_peaks(args.primary_peak_dir / f"{condition}.narrowPeak"),
            "n_runs": count,
        })
    fields = sorted({key for row in rows for key in row})
    preferred = ["level", "run_accession", "condition", "study", "replicate_type", "biological_replicate", "technical_replicate", "layout"]
    fields = preferred + [field for field in fields if field not in preferred]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    print(f"Wrote QC for {len(samples)} runs and {len(condition_counts)} conditions")


if __name__ == "__main__":
    main()
