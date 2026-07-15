#!/usr/bin/env python3
"""
Rebuild the promoter-level CRISPRa targetability atlas.

This script is the canonical matrix-level rebuild used for the manuscript
analyses. It scans both DNA strands with proper IUPAC reverse-complement
handling, integrates PAM hits with per-state ATAC-seq peaks, and writes:

  * supplementary/table_S2_targetability_full.tsv
  * supplementary/table_S3_sgrna_recommendations.csv

Raw FASTA and upstream ATAC files are intentionally supplied as command-line
inputs because the mouse genome FASTA is not stored in this repository.
"""
from __future__ import annotations

import argparse
import csv
import re
from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

STATES = [
    "homeostatic",
    "PP_naive",
    "PL_acute_LPS",
    "LL_tolerized",
    "sham_WT",
    "stroke_WT",
]

CAS_ORDER = [
    "HEAL_Un1Cas12f1",
    "SminiCRa_Un1Cas12f1",
    "SaCas9",
    "SpCas9",
    "CjCas9_MiniCAFE",
    "Nme2Cas9",
]

CAS_ORTHOLOGS = {
    "HEAL_Un1Cas12f1": {"pam": "TTTR", "pam_side": "5prime", "spacer_length": 20},
    "SminiCRa_Un1Cas12f1": {"pam": "TTTR", "pam_side": "5prime", "spacer_length": 20},
    "SaCas9": {"pam": "NNGRRT", "pam_side": "3prime", "spacer_length": 21},
    "SpCas9": {"pam": "NGG", "pam_side": "3prime", "spacer_length": 20},
    "CjCas9_MiniCAFE": {"pam": "NNNVRYM", "pam_side": "3prime", "spacer_length": 22},
    "Nme2Cas9": {"pam": "NNNNCC", "pam_side": "3prime", "spacer_length": 24},
}

IUPAC_REGEX = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "R": "[AG]",
    "Y": "[CT]",
    "M": "[AC]",
    "K": "[GT]",
    "S": "[GC]",
    "W": "[AT]",
    "H": "[ACT]",
    "B": "[CGT]",
    "V": "[ACG]",
    "D": "[AGT]",
    "N": "[ACGT]",
}

DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
IUPAC_COMP = str.maketrans("ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn")
RESTRICTION_SITES = ("CGTCTC", "GAGACG", "GGTCTC", "GAGACC")


@dataclass(frozen=True)
class Promoter:
    chrom: str
    start: int
    end: int
    name: str
    strand: str

    @property
    def midpoint(self) -> int:
        return (self.start + self.end) // 2


@dataclass(frozen=True)
class PamHit:
    gene: str
    cas: str
    chrom: str
    pam_strand: str
    pam_start: int
    pam_end: int
    pam_seq: str
    spacer_start: int
    spacer_end: int
    protospacer: str
    gc_content: float
    heuristic_score: str
    pams_in_peak_by_state: dict[str, bool]


class PeakIndex:
    def __init__(self, intervals_by_chrom: dict[str, list[tuple[int, int]]]):
        self.starts: dict[str, list[int]] = {}
        self.ends: dict[str, list[int]] = {}
        for chrom, intervals in intervals_by_chrom.items():
            merged = merge_intervals(intervals)
            self.starts[chrom] = [start for start, _ in merged]
            self.ends[chrom] = [end for _, end in merged]

    def contains(self, chrom: str, pos: int) -> bool:
        starts = self.starts.get(chrom)
        if not starts:
            return False
        idx = bisect_right(starts, pos) - 1
        return idx >= 0 and pos < self.ends[chrom][idx]


def merge_intervals(intervals: list[tuple[int, int]]) -> list[tuple[int, int]]:
    if not intervals:
        return []
    intervals = sorted(intervals)
    merged = [intervals[0]]
    for start, end in intervals[1:]:
        last_start, last_end = merged[-1]
        if start <= last_end:
            merged[-1] = (last_start, max(last_end, end))
        else:
            merged.append((start, end))
    return merged


def reverse_complement_dna(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1].upper()


def reverse_complement_iupac(seq: str) -> str:
    return seq.translate(IUPAC_COMP)[::-1].upper()


def iupac_to_regex(seq: str) -> str:
    return "".join(IUPAC_REGEX[base] for base in seq.upper())


def parse_fasta(path: Path) -> dict[str, str]:
    genome: dict[str, list[str]] = {}
    chrom = None
    with path.open() as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                chrom = line[1:].split()[0]
                genome[chrom] = []
            elif chrom is not None:
                genome[chrom].append(line.upper())
    return {chrom: "".join(parts) for chrom, parts in genome.items()}


def parse_promoters(path: Path) -> list[Promoter]:
    promoters: list[Promoter] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            strand = fields[5] if len(fields) > 5 else "+"
            promoters.append(Promoter(fields[0], int(fields[1]), int(fields[2]), fields[3], strand))
    return promoters


def load_peaks(path: Path) -> PeakIndex:
    intervals: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 3:
                intervals[fields[0]].append((int(fields[1]), int(fields[2])))
    return PeakIndex(intervals)


def load_therapeutic_genes(path: Path) -> list[str]:
    delimiter = "\t" if path.suffix == ".tsv" else ","
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        return [row["gene_symbol"] for row in reader]


def compute_gc(seq: str) -> float:
    return sum(1 for base in seq.upper() if base in "GC") / len(seq)


def has_homopolymer(seq: str, max_len: int = 4) -> bool:
    seq = seq.upper()
    return any(base * (max_len + 1) in seq for base in "ACGT")


def has_restriction_site(seq: str) -> bool:
    seq = seq.upper()
    rc = reverse_complement_dna(seq)
    return any(site in seq or site in rc for site in RESTRICTION_SITES)


def gc_score(gc: float) -> float:
    if 0.40 <= gc <= 0.60:
        return 1.0
    if 0.30 <= gc < 0.40:
        return (gc - 0.30) / 0.10
    if 0.60 < gc <= 0.70:
        return (0.70 - gc) / 0.10
    return 0.0


def longest_homopolymer(seq: str) -> int:
    longest = 0
    current_base = None
    current_len = 0
    for base in seq.upper():
        if base == current_base:
            current_len += 1
        else:
            current_base = base
            current_len = 1
        longest = max(longest, current_len)
    return longest


def heuristic_score(seq: str, cas: str) -> str:
    if cas not in {"HEAL_Un1Cas12f1", "SminiCRa_Un1Cas12f1"}:
        return "NA"
    gc = compute_gc(seq)
    hp = longest_homopolymer(seq)
    homopolymer_score = 1.0 if hp <= 3 else 0.5 if hp == 4 else 0.0
    poly_t_score = 0.0 if "TTTTT" in seq.upper() else 1.0
    score = 0.5 * gc_score(gc) + 0.3 * homopolymer_score + 0.2 * poly_t_score
    return f"{score:.3f}"


def scan_promoter(promoter: Promoter, seq: str, cas: str, peak_by_state: dict[str, PeakIndex]) -> list[PamHit]:
    info = CAS_ORTHOLOGS[cas]
    pam = info["pam"]
    pam_side = info["pam_side"]
    spacer_len = info["spacer_length"]
    pam_len = len(pam)
    seq_len = len(seq)
    hits: list[PamHit] = []

    patterns = [
        ("+", re.compile(f"(?=({iupac_to_regex(pam)}))")),
        ("-", re.compile(f"(?=({iupac_to_regex(reverse_complement_iupac(pam))}))")),
    ]

    for pam_strand, pattern in patterns:
        for match in pattern.finditer(seq):
            pam_start = match.start()
            pam_end = pam_start + pam_len
            if pam_side == "5prime":
                if pam_strand == "+":
                    spacer_start = pam_end
                    spacer_end = spacer_start + spacer_len
                else:
                    spacer_end = pam_start
                    spacer_start = spacer_end - spacer_len
            else:
                if pam_strand == "+":
                    spacer_end = pam_start
                    spacer_start = spacer_end - spacer_len
                else:
                    spacer_start = pam_end
                    spacer_end = spacer_start + spacer_len
            if spacer_start < 0 or spacer_end > seq_len:
                continue

            protospacer_sense = seq[spacer_start:spacer_end]
            protospacer = protospacer_sense if pam_strand == "+" else reverse_complement_dna(protospacer_sense)
            if "N" in protospacer:
                continue

            gc = compute_gc(protospacer)
            if not (0.30 <= gc <= 0.70):
                continue
            if has_homopolymer(protospacer) or has_restriction_site(protospacer):
                continue

            genomic_pam_start = promoter.start + pam_start
            genomic_pam_end = promoter.start + pam_end
            pams_in_peak = {
                state: peak_by_state[state].contains(promoter.chrom, genomic_pam_start)
                for state in STATES
            }
            pam_seq = seq[pam_start:pam_end] if pam_strand == "+" else reverse_complement_dna(seq[pam_start:pam_end])
            hits.append(PamHit(
                gene=promoter.name,
                cas=cas,
                chrom=promoter.chrom,
                pam_strand=pam_strand,
                pam_start=genomic_pam_start,
                pam_end=genomic_pam_end,
                pam_seq=pam_seq,
                spacer_start=promoter.start + spacer_start,
                spacer_end=promoter.start + spacer_end,
                protospacer=protospacer,
                gc_content=gc,
                heuristic_score=heuristic_score(protospacer, cas),
                pams_in_peak_by_state=pams_in_peak,
            ))
    return hits


def recommendation_class(targetable: dict[str, bool]) -> str:
    n_states = sum(targetable.values())
    if n_states == 6:
        return "constitutive_atlas_candidate"
    if n_states == 0:
        return "not_targetable_under_atlas"
    if targetable["sham_WT"] or targetable["stroke_WT"]:
        return "state_conditional_atlas_candidate"
    return "limited_state_atlas_candidate"


def write_table_s2(
    path: Path,
    promoters: list[Promoter],
    hits_by_gene_cas: dict[tuple[str, str], list[PamHit]],
    accessible_by_gene_state: dict[tuple[str, str], bool],
    therapeutic: set[str],
) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow([
            "gene",
            "cas",
            "state",
            "promoter_accessible",
            "pams_total_passing",
            "pams_in_peak",
            "targetable",
            "is_therapeutic",
        ])
        for promoter in sorted(promoters, key=lambda item: item.name):
            for cas in CAS_ORDER:
                hits = hits_by_gene_cas.get((promoter.name, cas), [])
                for state in STATES:
                    accessible = accessible_by_gene_state[(promoter.name, state)]
                    in_peak = sum(1 for hit in hits if hit.pams_in_peak_by_state[state])
                    targetable = accessible and in_peak > 0
                    writer.writerow([
                        promoter.name,
                        cas,
                        state,
                        str(accessible),
                        len(hits),
                        in_peak,
                        str(targetable),
                        str(promoter.name in therapeutic),
                    ])


def write_table_s3(
    path: Path,
    therapeutic_genes: list[str],
    hits_by_gene_cas: dict[tuple[str, str], list[PamHit]],
    targetable_by_gene_cas_state: dict[tuple[str, str, str], bool],
) -> None:
    fields = [
        "gene_symbol",
        "cas_ortholog",
        "rank",
        "protospacer_sequence",
        "pam_sequence",
        "strand",
        "genomic_position",
        "gc_content",
        "heuristic_score",
        "atlas_targetable_homeostatic",
        "atlas_targetable_PP_naive",
        "atlas_targetable_PL_acute_LPS",
        "atlas_targetable_LL_tolerized",
        "atlas_targetable_sham_WT",
        "atlas_targetable_stroke_WT",
        "atlas_n_targetable_states",
        "recommendation_class",
        "note",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for gene in therapeutic_genes:
            for cas in CAS_ORDER:
                hits = hits_by_gene_cas.get((gene, cas), [])
                if not hits:
                    continue
                midpoint = sum(1 for state in STATES if targetable_by_gene_cas_state[(gene, cas, state)])
                ranked = sorted(
                    hits,
                    key=lambda hit: (
                        -sum(hit.pams_in_peak_by_state.values()),
                        -float(hit.heuristic_score) if hit.heuristic_score != "NA" else 0.0,
                        abs(hit.gc_content - 0.50),
                        hit.pam_start,
                        hit.pam_strand,
                    ),
                )[:3]
                targetable = {
                    state: targetable_by_gene_cas_state[(gene, cas, state)]
                    for state in STATES
                }
                for rank, hit in enumerate(ranked, start=1):
                    row = {
                        "gene_symbol": gene,
                        "cas_ortholog": cas,
                        "rank": rank,
                        "protospacer_sequence": hit.protospacer,
                        "pam_sequence": hit.pam_seq,
                        "strand": hit.pam_strand,
                        "genomic_position": f"{hit.chrom}:{hit.pam_start}-{hit.pam_end}",
                        "gc_content": f"{hit.gc_content:.3f}",
                        "heuristic_score": hit.heuristic_score,
                        "atlas_targetable_homeostatic": str(targetable["homeostatic"]),
                        "atlas_targetable_PP_naive": str(targetable["PP_naive"]),
                        "atlas_targetable_PL_acute_LPS": str(targetable["PL_acute_LPS"]),
                        "atlas_targetable_LL_tolerized": str(targetable["LL_tolerized"]),
                        "atlas_targetable_sham_WT": str(targetable["sham_WT"]),
                        "atlas_targetable_stroke_WT": str(targetable["stroke_WT"]),
                        "atlas_n_targetable_states": midpoint,
                        "recommendation_class": recommendation_class(targetable),
                        "note": (
                            "Predictive atlas candidate. State flags indicate gene/Cas-level "
                            "PAM+chromatin targetability; guide efficacy and off-target risk "
                            "require experimental/orthogonal validation."
                        ),
                    }
                    writer.writerow(row)


def default_peak_paths(root: Path) -> dict[str, Path]:
    return {
        "homeostatic": root / "data/phase2_results/peaks/gosselin_2017/homeostatic_peaks.narrowPeak",
        "PP_naive": root / "data/phase2_results/peaks/holtman_wendeln_2022/PP_peaks.narrowPeak",
        "PL_acute_LPS": root / "data/phase2_results/peaks/holtman_wendeln_2022/PL_peaks.narrowPeak",
        "LL_tolerized": root / "data/phase2_results/peaks/holtman_wendeln_2022/LL_peaks.narrowPeak",
        "sham_WT": root / "data/phase2_results/peaks_zhang/sham_WT_peaks.narrowPeak",
        "stroke_WT": root / "data/phase2_results/peaks_zhang/stroke_WT_peaks.narrowPeak",
    }


def parse_args() -> argparse.Namespace:
    defaults = default_peak_paths(ROOT)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=ROOT / "data/references/mm39.fa")
    parser.add_argument("--promoters", type=Path, default=ROOT / "data/processed/promoters_crispra_optimal.bed")
    parser.add_argument("--therapeutic", type=Path, default=ROOT / "supplementary/table_S1_therapeutic_genes.csv")
    parser.add_argument("--table-s2-out", type=Path, default=ROOT / "supplementary/table_S2_targetability_full.tsv")
    parser.add_argument("--table-s3-out", type=Path, default=ROOT / "supplementary/table_S3_sgrna_recommendations.csv")
    for state in STATES:
        parser.add_argument(f"--{state.replace('_', '-')}-peaks", type=Path, default=defaults[state])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    peak_paths = {state: getattr(args, f"{state}_peaks") for state in STATES}

    print("Loading genome FASTA...")
    genome = parse_fasta(args.fasta)
    print(f"  chromosomes: {len(genome)}")

    print("Loading promoters and peaks...")
    promoters = parse_promoters(args.promoters)
    therapeutic_genes = load_therapeutic_genes(args.therapeutic)
    therapeutic = set(therapeutic_genes)
    peaks = {state: load_peaks(path) for state, path in peak_paths.items()}
    print(f"  promoters: {len(promoters)}")
    print(f"  therapeutic genes: {len(therapeutic_genes)}")

    accessible_by_gene_state: dict[tuple[str, str], bool] = {}
    hits_by_gene_cas: dict[tuple[str, str], list[PamHit]] = {}
    targetable_by_gene_cas_state: dict[tuple[str, str, str], bool] = {}

    print("Scanning promoters...")
    for index, promoter in enumerate(promoters, start=1):
        if promoter.chrom not in genome:
            continue
        seq = genome[promoter.chrom][promoter.start:promoter.end].upper()
        for state in STATES:
            accessible_by_gene_state[(promoter.name, state)] = peaks[state].contains(promoter.chrom, promoter.midpoint)
        for cas in CAS_ORDER:
            hits = scan_promoter(promoter, seq, cas, peaks)
            hits_by_gene_cas[(promoter.name, cas)] = hits
            for state in STATES:
                targetable_by_gene_cas_state[(promoter.name, cas, state)] = (
                    accessible_by_gene_state[(promoter.name, state)]
                    and any(hit.pams_in_peak_by_state[state] for hit in hits)
                )
        if index % 5000 == 0:
            print(f"  scanned {index:,}/{len(promoters):,} promoters")

    print("Writing supplementary tables...")
    write_table_s2(args.table_s2_out, promoters, hits_by_gene_cas, accessible_by_gene_state, therapeutic)
    write_table_s3(args.table_s3_out, therapeutic_genes, hits_by_gene_cas, targetable_by_gene_cas_state)

    print("\nSummary: therapeutic PAM availability")
    for cas in CAS_ORDER:
        n = sum(1 for gene in therapeutic if hits_by_gene_cas.get((gene, cas)))
        print(f"  {cas}: {n}/{len(therapeutic)} ({100*n/len(therapeutic):.1f}%)")

    print("\nSummary: therapeutic targetability")
    for cas in CAS_ORDER:
        values = []
        for state in STATES:
            n = sum(1 for gene in therapeutic if targetable_by_gene_cas_state[(gene, cas, state)])
            values.append(f"{state}={n}/{len(therapeutic)} ({100*n/len(therapeutic):.1f}%)")
        print(f"  {cas}: " + "; ".join(values))

    for cas in ("HEAL_Un1Cas12f1", "SminiCRa_Un1Cas12f1"):
        invalid = sorted({
            hit.pam_seq
            for (gene, hit_cas), hits in hits_by_gene_cas.items()
            if hit_cas == cas
            for hit in hits
            if hit.pam_seq not in {"TTTA", "TTTG"}
        })
        if invalid:
            raise SystemExit(f"Invalid {cas} PAM sequences after strict rebuild: {invalid}")


if __name__ == "__main__":
    main()
