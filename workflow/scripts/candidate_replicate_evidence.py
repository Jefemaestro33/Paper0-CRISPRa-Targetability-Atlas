#!/usr/bin/env python3
"""Attach run-level ATAC peak evidence to every reported candidate site."""
from __future__ import annotations

import argparse
import bisect
import csv
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Peak:
    start: int
    end: int
    peak_id: str
    signal: float
    summit: int


class PeakIndex:
    def __init__(self, path: Path):
        by_chrom: dict[str, list[Peak]] = defaultdict(list)
        with path.open() as handle:
            for line_number, line in enumerate(handle, start=1):
                if not line.strip() or line.startswith(("#", "track")):
                    continue
                fields = line.rstrip("\n").split("\t")
                start, end = int(fields[1]), int(fields[2])
                peak_id = fields[3] if len(fields) > 3 and fields[3] not in {"", "."} else f"{path.stem}:{line_number}"
                signal = float(fields[6]) if len(fields) > 6 and fields[6] not in {"", ".", "-1"} else 0.0
                offset = int(float(fields[9])) if len(fields) > 9 and fields[9] not in {"", ".", "-1"} else (end - start) // 2
                by_chrom[fields[0]].append(Peak(start, end, peak_id, signal, start + offset))
        self.peaks = {chrom: sorted(values, key=lambda peak: (peak.start, peak.end)) for chrom, values in by_chrom.items()}
        self.starts = {chrom: [peak.start for peak in values] for chrom, values in self.peaks.items()}

    def evidence(self, chrom: str, start: int, end: int) -> dict:
        values = self.peaks.get(chrom, [])
        starts = self.starts.get(chrom, [])
        index = max(0, bisect.bisect_right(starts, start) - 1)
        while index > 0 and values[index - 1].end > start:
            index -= 1
        overlaps = []
        while index < len(values) and values[index].start < end:
            peak = values[index]
            overlap = min(end, peak.end) - max(start, peak.start)
            if overlap > 0:
                overlaps.append((peak, overlap))
            index += 1
        if not overlaps:
            return {
                "guide_fully_in_replicate_peak": False,
                "guide_any_replicate_peak_overlap": False,
                "overlap_bp": 0, "replicate_peak_id": "NA",
                "replicate_peak_signal": "NA", "distance_to_summit": "NA",
            }
        complete = [(peak, overlap) for peak, overlap in overlaps if peak.start <= start and peak.end >= end]
        ranked = sorted(
            complete or overlaps,
            key=lambda item: (-item[1], -item[0].signal, abs((start + end) // 2 - item[0].summit), item[0].peak_id),
        )
        peak, overlap = ranked[0]
        return {
            "guide_fully_in_replicate_peak": bool(complete),
            "guide_any_replicate_peak_overlap": True,
            "overlap_bp": overlap, "replicate_peak_id": peak.peak_id,
            "replicate_peak_signal": f"{peak.signal:.6f}",
            "distance_to_summit": (start + end) // 2 - peak.summit,
        }


def parse_interval(value: str) -> tuple[str, int, int]:
    chrom, coordinates = value.split(":", 1)
    start, end = map(int, coordinates.split("-", 1))
    return chrom, start, end


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--sample-manifest", type=Path, required=True)
    parser.add_argument("--replicate-peak-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    with args.candidates.open(newline="") as handle:
        candidates = list(csv.DictReader(handle))
    with args.sample_manifest.open(newline="") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    indexes = {
        row["run_accession"]: PeakIndex(args.replicate_peak_dir / f"{row['run_accession']}.narrowPeak")
        for row in samples
    }
    fields = [
        "candidate_row", "gene_symbol", "nuclease_pam_class", "rank",
        "target_interval", "run_accession", "study", "condition",
        "replicate_type", "biological_replicate", "technical_replicate",
        "guide_fully_in_replicate_peak", "guide_any_replicate_peak_overlap",
        "overlap_bp", "replicate_peak_id", "replicate_peak_signal",
        "distance_to_summit",
    ]
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for candidate_index, candidate in enumerate(candidates, start=2):
            chrom, start, end = parse_interval(candidate["target_interval"])
            for sample in samples:
                writer.writerow({
                    "candidate_row": candidate_index,
                    "gene_symbol": candidate["gene_symbol"],
                    "nuclease_pam_class": candidate["nuclease_pam_class"],
                    "rank": candidate["rank"],
                    "target_interval": candidate["target_interval"],
                    "run_accession": sample["run_accession"],
                    "study": sample["study"], "condition": sample["condition"],
                    "replicate_type": sample["replicate_type"],
                    "biological_replicate": sample["biological_replicate"],
                    "technical_replicate": sample["technical_replicate"],
                    **indexes[sample["run_accession"]].evidence(chrom, start, end),
                })
    print(f"Wrote run-level evidence for {len(candidates):,} candidates x {len(samples)} runs")


if __name__ == "__main__":
    main()
