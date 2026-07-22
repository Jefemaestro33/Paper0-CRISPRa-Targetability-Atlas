#!/usr/bin/env python3
"""Add a PAM-aware Bowtie-1 off-target screen for Un1Cas12f1 candidates.

This remains a computational pre-screen.  It enumerates GRCm39 protospacer
alignments with up to three mismatches, then retains only alignments with a
correctly oriented TTTR PAM.  It does not model bulges, chromatin, cleavage,
or experimental genotoxicity.
"""
from __future__ import annotations

import argparse
import bisect
import csv
import re
import subprocess
import tempfile
from collections import defaultdict
from pathlib import Path

import pysam


PAM_RE = re.compile(r"^TTT[AG]$")


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1].upper()


class IntervalIndex:
    def __init__(self, intervals: dict[str, list[tuple[int, int]]]):
        self.intervals = {}
        self.starts = {}
        for chrom, values in intervals.items():
            merged: list[list[int]] = []
            for start, end in sorted(values):
                if merged and start <= merged[-1][1]:
                    merged[-1][1] = max(merged[-1][1], end)
                else:
                    merged.append([start, end])
            self.intervals[chrom] = [(start, end) for start, end in merged]
            self.starts[chrom] = [start for start, _ in merged]

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


def exon_index(gtf: Path) -> IntervalIndex:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with gtf.open() as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 5 and fields[2] == "exon":
                intervals[fields[0]].append((int(fields[3]) - 1, int(fields[4])))
    return IntervalIndex(intervals)


def promoter_index(selection: Path) -> IntervalIndex:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with selection.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["definition"] == "ensembl_canonical":
                tss = int(row["tss"])
                intervals[row["chrom"]].append((max(0, tss - 1000), tss + 1001))
    return IntervalIndex(intervals)


def mismatch_count(field: str) -> int:
    return 0 if not field else len([item for item in field.split(",") if item])


def parse_interval(value: str) -> tuple[str, int, int]:
    chrom, coordinates = value.split(":", 1)
    start, end = coordinates.split("-", 1)
    return chrom, int(start), int(end)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--fasta", type=Path, required=True)
    parser.add_argument("--bowtie-index", type=Path, required=True)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--tss-selection", type=Path, required=True)
    parser.add_argument("--annotated-out", type=Path, required=True)
    parser.add_argument("--alignments-out", type=Path, required=True)
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    with args.candidates.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
        original_fields = list(rows[0]) if rows else []
    evaluated = [index for index, row in enumerate(rows) if row["nuclease_pam_class"] == "Un1Cas12f1_TTTR"]
    if not Path(str(args.fasta) + ".fai").exists():
        subprocess.run(["samtools", "faidx", str(args.fasta)], check=True)
    genome = pysam.FastaFile(str(args.fasta))
    exons = exon_index(args.gtf)
    promoters = promoter_index(args.tss_selection)

    with tempfile.TemporaryDirectory(prefix="paper0_offtargets_") as directory:
        query_path = Path(directory) / "queries.fa"
        output_path = Path(directory) / "bowtie.tsv"
        with query_path.open("w") as handle:
            for index in evaluated:
                handle.write(f">candidate_{index}\n{rows[index]['protospacer_sequence']}\n")
        with output_path.open("w") as handle:
            subprocess.run([
                "bowtie", "-f", "-v", "3", "-a", "--best", "--quiet",
                str(args.bowtie_index), str(query_path),
            ], stdout=handle, check=True)

        annotations: dict[int, list[dict]] = defaultdict(list)
        with output_path.open() as handle:
            for line in handle:
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 7:
                    continue
                index = int(fields[0].removeprefix("candidate_"))
                strand, chrom, start = fields[1], fields[2], int(fields[3])
                end = start + len(rows[index]["protospacer_sequence"])
                try:
                    if strand == "+":
                        pam_start, pam_end = start - 4, start
                        pam = genome.fetch(chrom, max(0, pam_start), max(0, pam_end)).upper()
                    else:
                        pam_start, pam_end = end, end + 4
                        pam = reverse_complement(genome.fetch(chrom, pam_start, pam_end))
                except (KeyError, ValueError):
                    continue
                if not PAM_RE.fullmatch(pam):
                    continue
                intended_chrom, intended_start, intended_end = parse_interval(rows[index]["protospacer_interval"])
                on_target = chrom == intended_chrom and start == intended_start and end == intended_end and strand == rows[index]["strand"]
                mismatches = mismatch_count(fields[7] if len(fields) > 7 else "")
                annotations[index].append({
                    "candidate_row": index + 2, "gene_symbol": rows[index]["gene_symbol"],
                    "protospacer_sequence": rows[index]["protospacer_sequence"],
                    "alignment_chrom": chrom, "alignment_start": start, "alignment_end": end,
                    "alignment_strand": strand, "mismatches": mismatches, "adjacent_pam": pam,
                    "on_target": on_target, "overlaps_coding_or_noncoding_exon": exons.overlaps(chrom, start, end),
                    "within_1kb_of_canonical_tss": promoters.overlaps(chrom, start, end),
                })

    alignment_fields = [
        "candidate_row", "gene_symbol", "protospacer_sequence", "alignment_chrom",
        "alignment_start", "alignment_end", "alignment_strand", "mismatches",
        "adjacent_pam", "on_target", "overlaps_coding_or_noncoding_exon",
        "within_1kb_of_canonical_tss",
    ]
    args.alignments_out.parent.mkdir(parents=True, exist_ok=True)
    with args.alignments_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=alignment_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index in evaluated:
            writer.writerows(annotations[index])

    added_fields = [
        "pam_valid_offtargets_0mm", "pam_valid_offtargets_1mm", "pam_valid_offtargets_2mm",
        "pam_valid_offtargets_3mm", "pam_valid_offtargets_total_le3mm",
        "pam_valid_offtargets_exonic_le3mm", "pam_valid_offtargets_promoter_le3mm",
    ]
    fields = original_fields + [field for field in added_fields if field not in original_fields]
    for index, row in enumerate(rows):
        if index not in evaluated:
            row["off_target_status"] = "not_evaluated_requires_nuclease_specific_model"
            row["off_target_summary"] = "Only Un1Cas12f1/TTTR candidates received the PAM-aware Bowtie screen."
            for field in added_fields:
                row[field] = ""
            continue
        off_targets = [hit for hit in annotations[index] if not hit["on_target"]]
        counts = {mismatch: sum(hit["mismatches"] == mismatch for hit in off_targets) for mismatch in range(4)}
        row["pam_valid_offtargets_0mm"] = counts[0]
        row["pam_valid_offtargets_1mm"] = counts[1]
        row["pam_valid_offtargets_2mm"] = counts[2]
        row["pam_valid_offtargets_3mm"] = counts[3]
        row["pam_valid_offtargets_total_le3mm"] = len(off_targets)
        row["pam_valid_offtargets_exonic_le3mm"] = sum(hit["overlaps_coding_or_noncoding_exon"] for hit in off_targets)
        row["pam_valid_offtargets_promoter_le3mm"] = sum(hit["within_1kb_of_canonical_tss"] for hit in off_targets)
        row["off_target_status"] = "preliminary_complete_alignment_PAM_aware_to_3_mismatches"
        row["off_target_summary"] = (
            f"{len(off_targets)} PAM-valid genomic alignments excluding intended target at <=3 mismatches; "
            f"{row['pam_valid_offtargets_exonic_le3mm']} overlap exons and "
            f"{row['pam_valid_offtargets_promoter_le3mm']} lie within 1 kb of a canonical TSS. "
            "Bulges and experimental genotoxicity were not evaluated."
        )
    args.annotated_out.parent.mkdir(parents=True, exist_ok=True)
    with args.annotated_out.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)
    genome.close()
    print(f"Annotated {len(evaluated)} Un1Cas12f1 candidates")


if __name__ == "__main__":
    main()
