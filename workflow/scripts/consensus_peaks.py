#!/usr/bin/env python3
"""Create conservative replicate-consensus narrowPeaks.

For the biological datasets used here, the primary rule is support in at least
two independent replicates.  Every pair of peaks from different replicates is
required to pass the requested reciprocal-overlap threshold.  The pairwise
intersections are then merged, so every retained base is supported by at least
one qualifying biological-replicate pair.  This avoids the transitive-chain
failure of connected-component approaches (A overlaps B and B overlaps C even
when A does not overlap C).
"""
from __future__ import annotations

import argparse
import itertools
import statistics
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class Peak:
    replicate: int
    chrom: str
    start: int
    end: int
    signal: float
    summit: int


def read_peaks(path: Path, replicate: int) -> list[Peak]:
    peaks = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            start, end = int(fields[1]), int(fields[2])
            signal = float(fields[6]) if len(fields) > 6 and fields[6] not in {".", "-1"} else 0.0
            offset = int(float(fields[9])) if len(fields) > 9 and fields[9] not in {".", "-1"} else (end - start) // 2
            peaks.append(Peak(replicate, fields[0], start, end, signal, start + offset))
    return peaks


def reciprocal_overlap(a: Peak, b: Peak) -> tuple[float, float]:
    overlap = max(0, min(a.end, b.end) - max(a.start, b.start))
    return overlap / (a.end - a.start), overlap / (b.end - b.start)


def qualifying_pair_intersections(
    peaks: list[Peak], threshold: float
) -> list[tuple[str, int, int, float, int, frozenset[int]]]:
    by_chrom: dict[str, list[Peak]] = defaultdict(list)
    for peak in peaks:
        by_chrom[peak.chrom].append(peak)
    intersections: list[tuple[str, int, int, float, int, frozenset[int]]] = []
    for chrom, chrom_peaks in by_chrom.items():
        ordered = sorted(chrom_peaks, key=lambda peak: (peak.start, peak.end, peak.replicate))
        for i, left in enumerate(ordered):
            for j in range(i + 1, len(ordered)):
                right = ordered[j]
                if right.start >= left.end:
                    break
                if left.replicate == right.replicate:
                    continue
                frac_left, frac_right = reciprocal_overlap(left, right)
                if frac_left >= threshold and frac_right >= threshold:
                    start, end = max(left.start, right.start), min(left.end, right.end)
                    summit = round(statistics.median([left.summit, right.summit]))
                    summit = min(max(summit, start), end - 1)
                    intersections.append((
                        chrom, start, end, statistics.mean([left.signal, right.signal]),
                        summit, frozenset([left.replicate, right.replicate]),
                    ))
    return intersections


def merge_intersections(
    intersections: list[tuple[str, int, int, float, int, frozenset[int]]]
) -> list[tuple[str, int, int, float, int, int]]:
    """Union overlapping pairwise intersections without bridging unsupported gaps."""
    rows: list[tuple[str, int, int, float, int, int]] = []
    ordered = sorted(intersections, key=lambda item: (item[0], item[1], item[2], item[4]))
    for chrom, group_iter in itertools.groupby(ordered, key=lambda item: item[0]):
        group = list(group_iter)
        current = [group[0]]
        current_end = group[0][2]
        for item in group[1:]:
            if item[1] <= current_end:
                current.append(item)
                current_end = max(current_end, item[2])
                continue
            start = min(value[1] for value in current)
            end = max(value[2] for value in current)
            signal = statistics.mean(value[3] for value in current)
            summit = round(statistics.median(value[4] for value in current))
            support = len(set().union(*(value[5] for value in current)))
            rows.append((chrom, start, end, signal, min(max(summit, start), end - 1), support))
            current = [item]
            current_end = item[2]
        start = min(value[1] for value in current)
        end = max(value[2] for value in current)
        signal = statistics.mean(value[3] for value in current)
        summit = round(statistics.median(value[4] for value in current))
        support = len(set().union(*(value[5] for value in current)))
        rows.append((chrom, start, end, signal, min(max(summit, start), end - 1), support))
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", nargs="+", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--min-replicates", type=int, required=True)
    parser.add_argument("--reciprocal-overlap", type=float, default=0.50)
    args = parser.parse_args()

    if args.min_replicates != 2:
        raise ValueError("This implementation is validated for the prespecified >=2-replicate rule")
    peaks = [peak for replicate, path in enumerate(args.inputs, start=1) for peak in read_peaks(path, replicate)]
    intersections = qualifying_pair_intersections(peaks, args.reciprocal_overlap)
    rows = merge_intersections(intersections) if intersections else []

    rows.sort(key=lambda row: (row[0], row[1], row[2]))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w") as handle:
        for index, (chrom, start, end, signal, summit, support) in enumerate(rows, start=1):
            score = min(1000, round(signal))
            handle.write(
                f"{chrom}\t{start}\t{end}\tconsensus_{index}\t{score}\t.\t"
                f"{signal:.6f}\t-1\t-1\t{summit-start}\t{support}\n"
            )
    print(f"Wrote {len(rows):,} consensus peaks supported by >= {args.min_replicates} replicates")


if __name__ == "__main__":
    main()
