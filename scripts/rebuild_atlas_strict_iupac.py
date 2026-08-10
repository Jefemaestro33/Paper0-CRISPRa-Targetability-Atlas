#!/usr/bin/env python3
"""Rebuild the guide-site-aware CRISPRa targetability matrix.

Primary targetability requires the complete protospacer-plus-PAM interval to be
contained in the context-specific primary ATAC-seq peak set.  Biological
contexts use replicate consensus; technical-only contexts remain labelled.
Promoter midpoint
coverage and any guide/peak overlap are retained only as prespecified
sensitivities.  The output distinguishes five nuclease/PAM classes; HEAL and
SminiCRa are recorded as two activation architectures represented by the same
Un1Cas12f1/TTTR targeting class.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import re
from bisect import bisect_right
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import pysam


ROOT = Path(__file__).resolve().parents[1]

STATES = [
    "homeostatic",
    "PP_control",
    "PL_acute_LPS",
    "LL_tolerized",
    "sham_WT",
    "stroke_WT",
]

CAS_ORDER = [
    "Un1Cas12f1_TTTR",
    "SaCas9_NNGRRT",
    "SpCas9_NGG",
    "CjCas9_NNNVRYM",
    "Nme2Cas9_NNNNCC",
]

CAS_CLASSES = {
    "Un1Cas12f1_TTTR": {
        "pam": "TTTR", "pam_side": "5prime", "spacer_length": 20,
        "systems": "HEAL;SminiCRa",
    },
    "SaCas9_NNGRRT": {
        "pam": "NNGRRT", "pam_side": "3prime", "spacer_length": 21,
        "systems": "SaCas9_CRISPRa",
    },
    "SpCas9_NGG": {
        "pam": "NGG", "pam_side": "3prime", "spacer_length": 20,
        "systems": "SpCas9_CRISPRa",
    },
    "CjCas9_NNNVRYM": {
        "pam": "NNNVRYM", "pam_side": "3prime", "spacer_length": 22,
        "systems": "MiniCAFE",
    },
    "Nme2Cas9_NNNNCC": {
        "pam": "NNNNCC", "pam_side": "3prime", "spacer_length": 24,
        "systems": "proposed_dNme2Cas9_activator",
    },
}

IUPAC_REGEX = {
    "A": "A", "C": "C", "G": "G", "T": "T", "R": "[AG]", "Y": "[CT]",
    "M": "[AC]", "K": "[GT]", "S": "[GC]", "W": "[AT]", "H": "[ACT]",
    "B": "[CGT]", "V": "[ACG]", "D": "[AGT]", "N": "[ACGT]",
}
DNA_COMP = str.maketrans("ACGTacgt", "TGCAtgca")
IUPAC_COMP = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn", "TGCAYRSWMKVHDBNtgcayrswmkvhdbn"
)
RESTRICTION_SITES = ("CGTCTC", "GAGACG", "GGTCTC", "GAGACC")


@dataclass(frozen=True)
class Promoter:
    chrom: str
    start: int
    end: int
    name: str
    strand: str
    tss: int
    gene_id: str
    transcript_id: str
    tss_definition: str
    selection_source: str

    @property
    def midpoint(self) -> int:
        # Preserve the historical transcription-oriented -225 coordinate on
        # both strands.  The arithmetic midpoint of the negative-strand
        # half-open BED interval is one base farther upstream because the
        # -50 endpoint is excluded.
        return self.tss - 225 if self.strand == "+" else self.tss + 225


@dataclass(frozen=True)
class PeakRecord:
    chrom: str
    start: int
    end: int
    peak_id: str
    signal: float
    summit: int


@dataclass(frozen=True)
class PeakMatch:
    fully_contained: bool
    any_overlap: bool
    overlap_bp: int
    peak_id: str
    peak_signal: float | None
    distance_to_summit: int | None


@dataclass(frozen=True)
class GuideHit:
    gene: str
    cas_class: str
    chrom: str
    target_strand: str
    pam_start: int
    pam_end: int
    pam_seq: str
    spacer_start: int
    spacer_end: int
    protospacer: str
    gc_content: float
    heuristic_score: str
    distance_to_tss: int
    peak_matches: dict[str, PeakMatch]

    @property
    def target_start(self) -> int:
        return min(self.pam_start, self.spacer_start)

    @property
    def target_end(self) -> int:
        return max(self.pam_end, self.spacer_end)


class PeakIndex:
    def __init__(self, records_by_chrom: dict[str, list[PeakRecord]]):
        self.records: dict[str, list[PeakRecord]] = {}
        self.starts: dict[str, list[int]] = {}
        for chrom, records in records_by_chrom.items():
            ordered = sorted(records, key=lambda record: (record.start, record.end, record.peak_id))
            self.records[chrom] = ordered
            self.starts[chrom] = [record.start for record in ordered]

    def overlaps(self, chrom: str, start: int, end: int) -> list[PeakRecord]:
        records = self.records.get(chrom, [])
        starts = self.starts.get(chrom, [])
        if not records or end <= start:
            return []
        index = max(0, bisect_right(starts, start) - 1)
        while index > 0 and records[index - 1].end > start:
            index -= 1
        found: list[PeakRecord] = []
        while index < len(records) and records[index].start < end:
            record = records[index]
            if record.end > start:
                found.append(record)
            index += 1
        return found

    def contains_point(self, chrom: str, position: int) -> bool:
        return bool(self.overlaps(chrom, position, position + 1))

    def match(self, chrom: str, start: int, end: int) -> PeakMatch:
        overlaps = self.overlaps(chrom, start, end)
        if not overlaps:
            return PeakMatch(False, False, 0, "", None, None)
        center = (start + end) // 2
        ranked = sorted(
            overlaps,
            key=lambda peak: (
                -(min(end, peak.end) - max(start, peak.start)),
                -peak.signal,
                abs(center - peak.summit),
                peak.peak_id,
            ),
        )
        full = [peak for peak in ranked if peak.start <= start and peak.end >= end]
        best = full[0] if full else ranked[0]
        return PeakMatch(
            fully_contained=bool(full),
            any_overlap=True,
            overlap_bp=min(end, best.end) - max(start, best.start),
            peak_id=best.peak_id,
            peak_signal=best.signal,
            distance_to_summit=center - best.summit,
        )


def reverse_complement_dna(seq: str) -> str:
    return seq.translate(DNA_COMP)[::-1].upper()


def reverse_complement_iupac(seq: str) -> str:
    return seq.translate(IUPAC_COMP)[::-1].upper()


def iupac_to_regex(seq: str) -> str:
    return "".join(IUPAC_REGEX[base] for base in seq.upper())


def parse_promoters(path: Path) -> list[Promoter]:
    promoters: list[Promoter] = []
    with path.open() as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 11:
                raise ValueError(
                    f"Promoter BED must contain 11 columns from prepare_reference.py: {line[:120]}"
                )
            promoters.append(Promoter(
                chrom=fields[0], start=int(fields[1]), end=int(fields[2]), name=fields[3],
                strand=fields[5], tss=int(fields[6]), gene_id=fields[7], transcript_id=fields[8],
                tss_definition=fields[9], selection_source=fields[10],
            ))
    return promoters


def load_peaks(path: Path) -> PeakIndex:
    records: dict[str, list[PeakRecord]] = defaultdict(list)
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            start, end = int(fields[1]), int(fields[2])
            peak_id = fields[3] if len(fields) > 3 and fields[3] not in {"", "."} else f"{path.stem}:{line_number}"
            signal = float(fields[6]) if len(fields) > 6 and fields[6] not in {"", ".", "-1"} else 0.0
            offset = int(float(fields[9])) if len(fields) > 9 and fields[9] not in {"", ".", "-1"} else (end - start) // 2
            records[fields[0]].append(PeakRecord(fields[0], start, end, peak_id, signal, start + offset))
    return PeakIndex(records)


def load_therapeutic_genes(path: Path) -> list[str]:
    delimiter = "\t" if path.suffix == ".tsv" else ","
    with path.open(newline="") as handle:
        return [row["gene_symbol"] for row in csv.DictReader(handle, delimiter=delimiter)]


def compute_gc(seq: str) -> float:
    return sum(base in "GCgc" for base in seq) / len(seq)


def longest_homopolymer(seq: str) -> int:
    longest = current = 0
    previous = ""
    for base in seq.upper():
        current = current + 1 if base == previous else 1
        previous = base
        longest = max(longest, current)
    return longest


def passes_sequence_filters(seq: str) -> bool:
    gc = compute_gc(seq)
    if not 0.30 <= gc <= 0.70 or longest_homopolymer(seq) > 4:
        return False
    sequence_and_rc = seq.upper() + reverse_complement_dna(seq)
    return not any(site in sequence_and_rc for site in RESTRICTION_SITES)


def gc_score(gc: float) -> float:
    if 0.40 <= gc <= 0.60:
        return 1.0
    if 0.30 <= gc < 0.40:
        return (gc - 0.30) / 0.10
    if 0.60 < gc <= 0.70:
        return (0.70 - gc) / 0.10
    return 0.0


def heuristic_score(seq: str, cas_class: str) -> str:
    if cas_class != "Un1Cas12f1_TTTR":
        return "NA"
    hp = longest_homopolymer(seq)
    hp_score = 1.0 if hp <= 3 else 0.5 if hp == 4 else 0.0
    poly_t_score = 0.0 if "TTTTT" in seq.upper() else 1.0
    score = 0.5 * gc_score(compute_gc(seq)) + 0.3 * hp_score + 0.2 * poly_t_score
    return f"{score:.3f}"


def oriented_distance(promoter: Promoter, genomic_position: int) -> int:
    return genomic_position - promoter.tss if promoter.strand == "+" else promoter.tss - genomic_position


def scan_promoter(
    promoter: Promoter,
    seq: str,
    cas_class: str,
    peaks: dict[str, PeakIndex],
) -> list[GuideHit]:
    info = CAS_CLASSES[cas_class]
    pam = info["pam"]
    pam_len = len(pam)
    spacer_len = info["spacer_length"]
    patterns = (
        ("+", re.compile(f"(?=({iupac_to_regex(pam)}))")),
        ("-", re.compile(f"(?=({iupac_to_regex(reverse_complement_iupac(pam))}))")),
    )
    hits: list[GuideHit] = []
    for target_strand, pattern in patterns:
        for match in pattern.finditer(seq):
            pam_start = match.start()
            pam_end = pam_start + pam_len
            if info["pam_side"] == "5prime":
                spacer_start, spacer_end = (
                    (pam_end, pam_end + spacer_len) if target_strand == "+"
                    else (pam_start - spacer_len, pam_start)
                )
            else:
                spacer_start, spacer_end = (
                    (pam_start - spacer_len, pam_start) if target_strand == "+"
                    else (pam_end, pam_end + spacer_len)
                )
            if spacer_start < 0 or spacer_end > len(seq):
                continue
            spacer_sense = seq[spacer_start:spacer_end]
            protospacer = spacer_sense if target_strand == "+" else reverse_complement_dna(spacer_sense)
            if "N" in protospacer or not passes_sequence_filters(protospacer):
                continue
            genomic_pam_start = promoter.start + pam_start
            genomic_pam_end = promoter.start + pam_end
            genomic_spacer_start = promoter.start + spacer_start
            genomic_spacer_end = promoter.start + spacer_end
            target_start = min(genomic_pam_start, genomic_spacer_start)
            target_end = max(genomic_pam_end, genomic_spacer_end)
            target_center = (target_start + target_end) // 2
            matches = {
                state: peaks[state].match(promoter.chrom, target_start, target_end)
                for state in STATES
            }
            pam_seq = seq[pam_start:pam_end]
            if target_strand == "-":
                pam_seq = reverse_complement_dna(pam_seq)
            hits.append(GuideHit(
                gene=promoter.name,
                cas_class=cas_class,
                chrom=promoter.chrom,
                target_strand=target_strand,
                pam_start=genomic_pam_start,
                pam_end=genomic_pam_end,
                pam_seq=pam_seq,
                spacer_start=genomic_spacer_start,
                spacer_end=genomic_spacer_end,
                protospacer=protospacer,
                gc_content=compute_gc(protospacer),
                heuristic_score=heuristic_score(protospacer, cas_class),
                distance_to_tss=oriented_distance(promoter, target_center),
                peak_matches=matches,
            ))
    return hits


def candidate_class(accessible: dict[str, bool]) -> str:
    n_states = sum(accessible.values())
    if n_states == len(STATES):
        return "constitutive_guide_site_candidate"
    if n_states == 0:
        return "no_primary_peak_support"
    if accessible["sham_WT"] or accessible["stroke_WT"]:
        return "dataset_context_conditional_candidate"
    return "limited_context_candidate"


def minimum_primary_signal(hit: GuideHit) -> float:
    """Minimum signal among contexts that fully contain the complete guide."""
    values = [
        match.peak_signal for match in hit.peak_matches.values()
        if match.fully_contained and match.peak_signal is not None
    ]
    return min(values) if values else 0.0


def write_table_s2(
    path: Path,
    promoters: list[Promoter],
    hits_by_gene_class: dict[tuple[str, str], list[GuideHit]],
    peaks: dict[str, PeakIndex],
    therapeutic: set[str],
) -> None:
    fields = [
        "gene", "gene_id", "transcript_id", "tss", "strand", "tss_definition",
        "selection_source", "cas", "represented_systems", "state",
        "promoter_midpoint_accessible", "promoter_any_peak_overlap",
        "protospacers_total_passing", "protospacers_fully_in_peak",
        "protospacers_any_peak_overlap", "targetable", "is_therapeutic",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for promoter in sorted(promoters, key=lambda item: item.name):
            for cas_class in CAS_ORDER:
                hits = hits_by_gene_class.get((promoter.name, cas_class), [])
                for state in STATES:
                    full = sum(hit.peak_matches[state].fully_contained for hit in hits)
                    any_overlap = sum(hit.peak_matches[state].any_overlap for hit in hits)
                    writer.writerow({
                        "gene": promoter.name,
                        "gene_id": promoter.gene_id,
                        "transcript_id": promoter.transcript_id,
                        "tss": promoter.tss,
                        "strand": promoter.strand,
                        "tss_definition": promoter.tss_definition,
                        "selection_source": promoter.selection_source,
                        "cas": cas_class,
                        "represented_systems": CAS_CLASSES[cas_class]["systems"],
                        "state": state,
                        "promoter_midpoint_accessible": peaks[state].contains_point(promoter.chrom, promoter.midpoint),
                        "promoter_any_peak_overlap": bool(peaks[state].overlaps(promoter.chrom, promoter.start, promoter.end)),
                        "protospacers_total_passing": len(hits),
                        "protospacers_fully_in_peak": full,
                        "protospacers_any_peak_overlap": any_overlap,
                        "targetable": full > 0,
                        "is_therapeutic": promoter.name in therapeutic,
                    })


def state_fields() -> list[str]:
    fields: list[str] = []
    for state in STATES:
        fields.extend([
            f"guide_fully_in_peak_{state}",
            f"guide_any_peak_overlap_{state}",
            f"peak_id_{state}",
            f"peak_signal_{state}",
            f"distance_to_summit_{state}",
        ])
    return fields


def write_table_s3(
    path: Path,
    promoter_by_gene: dict[str, Promoter],
    therapeutic_genes: list[str],
    hits_by_gene_class: dict[tuple[str, str], list[GuideHit]],
) -> None:
    fields = [
        "gene_symbol", "gene_id", "transcript_id", "tss", "tss_definition",
        "nuclease_pam_class", "represented_systems", "rank", "protospacer_sequence",
        "pam_sequence", "strand", "target_interval", "protospacer_interval", "pam_interval",
        "distance_to_tss", "gc_content", "heuristic_score", "n_primary_states",
        "candidate_class",
    ] + state_fields() + ["off_target_status", "off_target_summary", "note"]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for gene in therapeutic_genes:
            promoter = promoter_by_gene.get(gene)
            if promoter is None:
                continue
            for cas_class in CAS_ORDER:
                hits = hits_by_gene_class.get((gene, cas_class), [])
                ranked = sorted(
                    hits,
                    key=lambda hit: (
                        -sum(hit.peak_matches[state].fully_contained for state in STATES),
                        -sum(hit.peak_matches[state].any_overlap for state in STATES),
                        -minimum_primary_signal(hit),
                        -float(hit.heuristic_score) if hit.heuristic_score != "NA" else 0.0,
                        abs(hit.gc_content - 0.50),
                        abs(hit.distance_to_tss + 200),
                        hit.target_start,
                    ),
                )[:5]
                for rank, hit in enumerate(ranked, start=1):
                    primary = {state: hit.peak_matches[state].fully_contained for state in STATES}
                    row = {
                        "gene_symbol": gene,
                        "gene_id": promoter.gene_id,
                        "transcript_id": promoter.transcript_id,
                        "tss": promoter.tss,
                        "tss_definition": promoter.tss_definition,
                        "nuclease_pam_class": cas_class,
                        "represented_systems": CAS_CLASSES[cas_class]["systems"],
                        "rank": rank,
                        "protospacer_sequence": hit.protospacer,
                        "pam_sequence": hit.pam_seq,
                        "strand": hit.target_strand,
                        "target_interval": f"{hit.chrom}:{hit.target_start}-{hit.target_end}",
                        "protospacer_interval": f"{hit.chrom}:{hit.spacer_start}-{hit.spacer_end}",
                        "pam_interval": f"{hit.chrom}:{hit.pam_start}-{hit.pam_end}",
                        "distance_to_tss": hit.distance_to_tss,
                        "gc_content": f"{hit.gc_content:.3f}",
                        "heuristic_score": hit.heuristic_score,
                        "n_primary_states": sum(primary.values()),
                        "candidate_class": candidate_class(primary),
                        "off_target_status": "not_evaluated_by_targetability_rebuild",
                        "off_target_summary": "See nuclease-aware external validation before use.",
                        "note": (
                            "Guide-specific predictive candidate. Primary flags require the complete "
                            "protospacer+PAM interval inside a primary ATAC peak; replicate status "
                            "is recorded separately, and the flags do not "
                            "demonstrate CRISPRa efficacy or safety."
                        ),
                    }
                    for state in STATES:
                        match = hit.peak_matches[state]
                        row[f"guide_fully_in_peak_{state}"] = match.fully_contained
                        row[f"guide_any_peak_overlap_{state}"] = match.any_overlap
                        row[f"peak_id_{state}"] = match.peak_id
                        row[f"peak_signal_{state}"] = "" if match.peak_signal is None else f"{match.peak_signal:.4f}"
                        row[f"distance_to_summit_{state}"] = "" if match.distance_to_summit is None else match.distance_to_summit
                    writer.writerow(row)


def default_peak_paths(root: Path) -> dict[str, Path]:
    base = root / "workflow/results/peaks/primary"
    return {state: base / f"{state}.narrowPeak" for state in STATES}


def parse_args() -> argparse.Namespace:
    defaults = default_peak_paths(ROOT)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--fasta", type=Path, default=ROOT / "workflow/resources/mm39.fa")
    parser.add_argument("--promoters", type=Path, default=ROOT / "reference/promoters_ensembl_canonical.bed")
    parser.add_argument("--therapeutic", type=Path, default=ROOT / "config/therapeutic_genes_locked.csv")
    parser.add_argument("--table-s2-out", type=Path, default=ROOT / "supplementary/table_S2_targetability_full.tsv.gz")
    parser.add_argument("--table-s3-out", type=Path, default=ROOT / "supplementary/table_S3_candidate_protospacers.csv")
    for state in STATES:
        parser.add_argument(f"--{state.replace('_', '-')}-peaks", type=Path, default=defaults[state])
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    peak_paths = {state: getattr(args, f"{state}_peaks") for state in STATES}
    missing = [path for path in [args.fasta, args.promoters, *peak_paths.values()] if not path.exists()]
    if missing:
        raise SystemExit("Missing required inputs:\n  " + "\n  ".join(map(str, missing)))

    print("Loading genome, promoters, and peak sets...")
    if not Path(str(args.fasta) + ".fai").exists():
        pysam.faidx(str(args.fasta))
    genome = pysam.FastaFile(str(args.fasta))
    promoters = parse_promoters(args.promoters)
    promoter_by_gene = {promoter.name: promoter for promoter in promoters}
    therapeutic_genes = load_therapeutic_genes(args.therapeutic)
    therapeutic = set(therapeutic_genes)
    peaks = {state: load_peaks(path) for state, path in peak_paths.items()}
    print(f"  chromosomes={len(genome.references)} promoters={len(promoters):,} therapeutic={len(therapeutic)}")

    s2_fields = [
        "gene", "gene_id", "transcript_id", "tss", "strand", "tss_definition",
        "selection_source", "cas", "represented_systems", "state",
        "promoter_midpoint_accessible", "promoter_any_peak_overlap",
        "protospacers_total_passing", "protospacers_fully_in_peak",
        "protospacers_any_peak_overlap", "targetable", "is_therapeutic",
    ]
    therapeutic_hits: dict[tuple[str, str], list[GuideHit]] = {}
    invalid_tttr: set[str] = set()
    args.table_s2_out.parent.mkdir(parents=True, exist_ok=True)
    table_s2_open = gzip.open if args.table_s2_out.suffix == ".gz" else open
    with table_s2_open(args.table_s2_out, "wt", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=s2_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for index, promoter in enumerate(promoters, start=1):
            try:
                sequence = genome.fetch(promoter.chrom, promoter.start, promoter.end).upper()
            except (KeyError, ValueError):
                sequence = ""
            if len(sequence) != promoter.end - promoter.start:
                continue
            for cas_class in CAS_ORDER:
                hits = scan_promoter(promoter, sequence, cas_class, peaks)
                if promoter.name in therapeutic:
                    therapeutic_hits[(promoter.name, cas_class)] = hits
                if cas_class == "Un1Cas12f1_TTTR":
                    invalid_tttr.update(hit.pam_seq for hit in hits if hit.pam_seq not in {"TTTA", "TTTG"})
                for state in STATES:
                    full = sum(hit.peak_matches[state].fully_contained for hit in hits)
                    any_overlap = sum(hit.peak_matches[state].any_overlap for hit in hits)
                    writer.writerow({
                        "gene": promoter.name,
                        "gene_id": promoter.gene_id,
                        "transcript_id": promoter.transcript_id,
                        "tss": promoter.tss,
                        "strand": promoter.strand,
                        "tss_definition": promoter.tss_definition,
                        "selection_source": promoter.selection_source,
                        "cas": cas_class,
                        "represented_systems": CAS_CLASSES[cas_class]["systems"],
                        "state": state,
                        "promoter_midpoint_accessible": peaks[state].contains_point(promoter.chrom, promoter.midpoint),
                        "promoter_any_peak_overlap": bool(peaks[state].overlaps(promoter.chrom, promoter.start, promoter.end)),
                        "protospacers_total_passing": len(hits),
                        "protospacers_fully_in_peak": full,
                        "protospacers_any_peak_overlap": any_overlap,
                        "targetable": full > 0,
                        "is_therapeutic": promoter.name in therapeutic,
                    })
            if index % 5000 == 0:
                print(f"  scanned {index:,}/{len(promoters):,}")

    write_table_s3(args.table_s3_out, promoter_by_gene, therapeutic_genes, therapeutic_hits)
    genome.close()

    print("Therapeutic-gene summary")
    for cas_class in CAS_ORDER:
        pam_n = sum(bool(therapeutic_hits.get((gene, cas_class))) for gene in therapeutic)
        values = []
        for state in STATES:
            n = sum(
                any(hit.peak_matches[state].fully_contained for hit in therapeutic_hits.get((gene, cas_class), []))
                for gene in therapeutic
            )
            values.append(f"{state}={n}/{len(therapeutic)}")
        print(f"  {cas_class}: PAM={pam_n}/{len(therapeutic)}; " + "; ".join(values))

    invalid = sorted(invalid_tttr)
    if invalid:
        raise SystemExit(f"Invalid TTTR PAM sequences after strict rebuild: {invalid}")


if __name__ == "__main__":
    main()
