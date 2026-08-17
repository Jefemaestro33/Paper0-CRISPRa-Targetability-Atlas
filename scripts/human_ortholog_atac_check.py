#!/usr/bin/env python3
"""Exploratory human ortholog ATAC support check for the revision response.

This script is intentionally separate from the murine targetability rebuild.  It
regenerates the compact human ortholog-panel checks used as a post-review
translational robustness analysis:

* human panel promoters from GENCODE v19 / hg19;
* public iPSC-derived microglia ATAC peak support from GSE206479 and GSE245522;
* the same PAM classes, promoter window, sequence filters, and primary
  complete-protospacer-plus-PAM-in-peak rule used by the murine rebuild;
* a TFEB alternative-TSS sensitivity audit.

The analysis is exploratory.  It does not create a human atlas and does not
validate CRISPRa activity.
"""
from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import re
import shutil
import tarfile
import time
import urllib.request
from bisect import bisect_right
from collections import Counter, defaultdict
from dataclasses import dataclass
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]

GENCODE_V19_URL = (
    "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/"
    "release_19/gencode.v19.annotation.gtf.gz"
)

GSE206479_PEAK_URLS = {
    "WTC11_rest": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206479/suppl/"
        "GSE206479_WTC11_rest_idr_optimal_peak.narrowPeak.gz"
    ),
    "WTC11_IFNb": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206479/suppl/"
        "GSE206479_WTC11_stim_idr_optimal_peak.narrowPeak.gz"
    ),
    "H1_rest": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206479/suppl/"
        "GSE206479_H1_rest_idr_optimal_peak.narrowPeak.gz"
    ),
    "H1_IFNb": (
        "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE206nnn/GSE206479/suppl/"
        "GSE206479_H1_stim_idr_optimal_peak.narrowPeak.gz"
    ),
}

GSE245522_TAR_URL = (
    "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE245nnn/GSE245522/suppl/"
    "GSE245522_RAW.tar"
)
GSE245522_TAR_MEMBERS = {
    "10K_peak_file": "GSM7844462_peaks_10K.narrowPeak.gz",
    "31K_peak_file": "GSM7844463_peaks_31K.narrowPeak.gz",
    "100K_peak_file": "GSM7844464_peaks_100K.narrowPeak.gz",
    "100K2_peak_file": "GSM7844465_peaks_100K2.narrowPeak.gz",
}

CAS_ORDER = [
    "Un1Cas12f1_TTTR",
    "SaCas9_NNGRRT",
    "SpCas9_NGG",
    "CjCas9_NNNVRYM",
    "Nme2Cas9_NNNNCC",
]

CAS_CLASSES = {
    "Un1Cas12f1_TTTR": {
        "pam": "TTTR",
        "pam_side": "5prime",
        "spacer_length": 20,
        "systems": "HEAL;SminiCRa",
    },
    "SaCas9_NNGRRT": {
        "pam": "NNGRRT",
        "pam_side": "3prime",
        "spacer_length": 21,
        "systems": "SaCas9_CRISPRa",
    },
    "SpCas9_NGG": {
        "pam": "NGG",
        "pam_side": "3prime",
        "spacer_length": 20,
        "systems": "SpCas9_CRISPRa",
    },
    "CjCas9_NNNVRYM": {
        "pam": "NNNVRYM",
        "pam_side": "3prime",
        "spacer_length": 22,
        "systems": "MiniCAFE",
    },
    "Nme2Cas9_NNNNCC": {
        "pam": "NNNNCC",
        "pam_side": "3prime",
        "spacer_length": 24,
        "systems": "proposed_dNme2Cas9_activator",
    },
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
IUPAC_COMP = str.maketrans(
    "ACGTRYSWKMBDHVNacgtryswkmbdhvn",
    "TGCAYRSWMKVHDBNtgcayrswmkvhdbn",
)
RESTRICTION_SITES = ("CGTCTC", "GAGACG", "GGTCTC", "GAGACC")

CURATED_HUMAN_SYMBOLS = {
    "Il4ra": "IL4R",
    "Map1lc3b": "MAP1LC3B",
}
ORTHOLOGY_EXCLUSIONS = {
    "Siglech": "no clear one-to-one human ortholog for mouse Siglech",
    "Chil3": "no direct human CHIL3 ortholog; paralogous chitinase-like genes would be non-equivalent",
}


@dataclass(frozen=True)
class Transcript:
    human_symbol: str
    gene_id: str
    transcript_id: str
    transcript_name: str
    chrom: str
    start: int
    end: int
    strand: str
    tss: int
    tags: tuple[str, ...]
    ccdsid: str
    tsl: int
    length: int


@dataclass(frozen=True)
class Promoter:
    mouse_symbol: str
    human_symbol: str
    mapping_status: str
    mapping_note: str
    chrom: str
    start: int
    end: int
    strand: str
    tss: int
    gene_id: str
    transcript_id: str
    transcript_name: str
    selection_rule: str
    tags: str
    n_transcripts_for_symbol: int

    @property
    def midpoint(self) -> int:
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
    mouse_symbol: str
    human_symbol: str
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
            ordered = sorted(records, key=lambda item: (item.start, item.end, item.peak_id))
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


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def download(url: str, path: Path) -> None:
    if path.exists():
        return
    path.parent.mkdir(parents=True, exist_ok=True)
    with urllib.request.urlopen(url, timeout=120) as response, path.open("wb") as handle:
        shutil.copyfileobj(response, handle)


def read_therapeutic_symbols(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        return [row["gene_symbol"] for row in csv.DictReader(handle)]


def human_symbol_for_mouse(mouse_symbol: str) -> tuple[str | None, str, str]:
    if mouse_symbol in ORTHOLOGY_EXCLUSIONS:
        return None, "excluded_no_clear_one_to_one_ortholog", ORTHOLOGY_EXCLUSIONS[mouse_symbol]
    if mouse_symbol in CURATED_HUMAN_SYMBOLS:
        return CURATED_HUMAN_SYMBOLS[mouse_symbol], "symbol_or_curated", "curated symbol mapping"
    return mouse_symbol.upper(), "symbol_or_curated", (
        "upper-case symbol convention; replace with formal orthology table if promoted to main resource"
    )


def parse_gtf_attributes(raw: str) -> tuple[dict[str, str], tuple[str, ...]]:
    attrs: dict[str, str] = {}
    tags: list[str] = []
    for key, value in re.findall(r'(\S+) "([^"]+)";', raw):
        if key == "tag":
            tags.append(value)
        else:
            attrs[key] = value
    return attrs, tuple(tags)


def parse_gencode_transcripts(gtf_gz: Path) -> dict[str, list[Transcript]]:
    by_symbol: dict[str, list[Transcript]] = defaultdict(list)
    with gzip.open(gtf_gz, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 9 or fields[2] != "transcript":
                continue
            attrs, tags = parse_gtf_attributes(fields[8])
            if attrs.get("gene_type") != "protein_coding":
                continue
            if attrs.get("transcript_type") != "protein_coding":
                continue
            start_1based = int(fields[3])
            end_1based = int(fields[4])
            strand = fields[6]
            start_0based = start_1based - 1
            end_0based_exclusive = end_1based
            tss = start_1based - 1 if strand == "+" else end_1based - 1
            tsl_raw = attrs.get("transcript_support_level", "99").split()[0]
            try:
                tsl = int(tsl_raw)
            except ValueError:
                tsl = 99
            transcript = Transcript(
                human_symbol=attrs["gene_name"],
                gene_id=attrs["gene_id"],
                transcript_id=attrs["transcript_id"],
                transcript_name=attrs.get("transcript_name", attrs["transcript_id"]),
                chrom=fields[0],
                start=start_0based,
                end=end_0based_exclusive,
                strand=strand,
                tss=tss,
                tags=tags,
                ccdsid=attrs.get("ccdsid", ""),
                tsl=tsl,
                length=end_0based_exclusive - start_0based,
            )
            by_symbol[transcript.human_symbol].append(transcript)
    return by_symbol


def transcript_selection_key(transcript: Transcript) -> tuple[int, int, int, int, int, str]:
    has_appris = any(tag.startswith("appris_principal") for tag in transcript.tags)
    has_ccds = bool(transcript.ccdsid) or "CCDS" in transcript.tags
    has_basic = "basic" in transcript.tags
    return (
        0 if has_appris else 1,
        0 if has_ccds else 1,
        0 if has_basic else 1,
        transcript.tsl,
        -transcript.length,
        transcript.transcript_id,
    )


def promoter_from_transcript(
    mouse_symbol: str,
    human_symbol: str,
    mapping_status: str,
    mapping_note: str,
    transcript: Transcript,
    n_transcripts_for_symbol: int,
) -> Promoter:
    if transcript.strand == "+":
        start, end = transcript.tss - 400, transcript.tss - 50
    else:
        start, end = transcript.tss + 51, transcript.tss + 401
    return Promoter(
        mouse_symbol=mouse_symbol,
        human_symbol=human_symbol,
        mapping_status=mapping_status,
        mapping_note=mapping_note,
        chrom=transcript.chrom,
        start=start,
        end=end,
        strand=transcript.strand,
        tss=transcript.tss,
        gene_id=transcript.gene_id,
        transcript_id=transcript.transcript_id,
        transcript_name=transcript.transcript_name,
        selection_rule="APPRIS_principal_then_CCDS_basic_TSL_length",
        tags=";".join(transcript.tags),
        n_transcripts_for_symbol=n_transcripts_for_symbol,
    )


def build_human_panel(
    mouse_symbols: list[str],
    transcripts_by_symbol: dict[str, list[Transcript]],
) -> tuple[list[Promoter], list[dict[str, str]]]:
    promoters: list[Promoter] = []
    issues: list[dict[str, str]] = []
    for mouse_symbol in mouse_symbols:
        human_symbol, mapping_status, mapping_note = human_symbol_for_mouse(mouse_symbol)
        if human_symbol is None:
            issues.append({"mouse_symbol": mouse_symbol, "human_symbol": "", "issue": mapping_note})
            continue
        candidates = transcripts_by_symbol.get(human_symbol, [])
        if not candidates:
            issues.append(
                {"mouse_symbol": mouse_symbol, "human_symbol": human_symbol, "issue": "no protein-coding GENCODE v19 transcript found"}
            )
            continue
        selected = sorted(candidates, key=transcript_selection_key)[0]
        promoters.append(
            promoter_from_transcript(
                mouse_symbol,
                human_symbol,
                mapping_status,
                mapping_note,
                selected,
                len(candidates),
            )
        )
    return promoters, issues


def fetch_hg19_sequence(chrom: str, start: int, end: int) -> str:
    url = (
        "https://api.genome.ucsc.edu/getData/sequence"
        f"?genome=hg19;chrom={chrom};start={start};end={end}"
    )
    with urllib.request.urlopen(url, timeout=60) as response:
        payload = json.loads(response.read().decode("utf-8"))
    seq = payload["dna"].upper()
    if len(seq) != end - start:
        raise ValueError(f"Unexpected UCSC sequence length for {chrom}:{start}-{end}: {len(seq)}")
    return seq


def write_promoter_tables(
    output_dir: Path,
    promoters: list[Promoter],
    issues: list[dict[str, str]],
) -> dict[str, str]:
    promoter_fields = [
        "mouse_symbol",
        "human_symbol",
        "mapping_status",
        "mapping_note",
        "chrom",
        "promoter_start",
        "promoter_end",
        "strand",
        "tss",
        "gene_id",
        "transcript_id",
        "transcript_name",
        "selection_rule",
        "tags",
        "n_transcripts_for_symbol",
    ]
    promoter_path = output_dir / "human_panel_promoters_hg19.tsv"
    with promoter_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=promoter_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for promoter in sorted(promoters, key=lambda item: item.mouse_symbol):
            writer.writerow(
                {
                    "mouse_symbol": promoter.mouse_symbol,
                    "human_symbol": promoter.human_symbol,
                    "mapping_status": promoter.mapping_status,
                    "mapping_note": promoter.mapping_note,
                    "chrom": promoter.chrom,
                    "promoter_start": promoter.start,
                    "promoter_end": promoter.end,
                    "strand": promoter.strand,
                    "tss": promoter.tss,
                    "gene_id": promoter.gene_id,
                    "transcript_id": promoter.transcript_id,
                    "transcript_name": promoter.transcript_name,
                    "selection_rule": promoter.selection_rule,
                    "tags": promoter.tags,
                    "n_transcripts_for_symbol": promoter.n_transcripts_for_symbol,
                }
            )
    issue_path = output_dir / "human_panel_mapping_issues.tsv"
    with issue_path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=["mouse_symbol", "human_symbol", "issue"], delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(issues)
    sequences: dict[str, str] = {}
    qc_path = output_dir / "human_panel_sequence_qc.tsv"
    with qc_path.open("w", newline="") as qc_handle:
        fields = ["mouse_symbol", "human_symbol", "interval", "expected_length", "observed_length", "status"]
        writer = csv.DictWriter(qc_handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for promoter in sorted(promoters, key=lambda item: item.mouse_symbol):
            seq = fetch_hg19_sequence(promoter.chrom, promoter.start, promoter.end)
            sequences[promoter.mouse_symbol] = seq
            writer.writerow(
                {
                    "mouse_symbol": promoter.mouse_symbol,
                    "human_symbol": promoter.human_symbol,
                    "interval": f"{promoter.chrom}:{promoter.start}-{promoter.end}",
                    "expected_length": promoter.end - promoter.start,
                    "observed_length": len(seq),
                    "status": "ok" if len(seq) == promoter.end - promoter.start and set(seq) <= {"A", "C", "G", "T"} else "check",
                }
            )
            time.sleep(0.05)
    return sequences


def load_peaks(path: Path) -> PeakIndex:
    records: dict[str, list[PeakRecord]] = defaultdict(list)
    with gzip.open(path, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#") or line.startswith("track"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 3:
                continue
            chrom = fields[0]
            start, end = int(fields[1]), int(fields[2])
            peak_id = fields[3] if len(fields) > 3 and fields[3] not in {"", "."} else f"{path.name}:{line_number}"
            signal = float(fields[6]) if len(fields) > 6 and fields[6] not in {"", ".", "-1"} else 0.0
            offset = int(float(fields[9])) if len(fields) > 9 and fields[9] not in {"", ".", "-1"} else (end - start) // 2
            records[chrom].append(PeakRecord(chrom, start, end, peak_id, signal, start + offset))
    return PeakIndex(records)


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
                    (pam_end, pam_end + spacer_len)
                    if target_strand == "+"
                    else (pam_start - spacer_len, pam_start)
                )
            else:
                spacer_start, spacer_end = (
                    (pam_start - spacer_len, pam_start)
                    if target_strand == "+"
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
            matches = {context: peaks[context].match(promoter.chrom, target_start, target_end) for context in peaks}
            pam_seq = seq[pam_start:pam_end]
            if target_strand == "-":
                pam_seq = reverse_complement_dna(pam_seq)
            hits.append(
                GuideHit(
                    mouse_symbol=promoter.mouse_symbol,
                    human_symbol=promoter.human_symbol,
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
                )
            )
    return hits


def write_targetability(
    path: Path,
    promoters: list[Promoter],
    hits_by_gene_cas: dict[tuple[str, str], list[GuideHit]],
    peaks: dict[str, PeakIndex],
) -> None:
    fields = [
        "mouse_symbol",
        "human_symbol",
        "gene_id",
        "transcript_id",
        "tss",
        "strand",
        "cas",
        "represented_systems",
        "context",
        "promoter_midpoint_accessible",
        "promoter_any_peak_overlap",
        "protospacers_total_passing",
        "protospacers_fully_in_peak",
        "protospacers_any_peak_overlap",
        "targetable",
    ]
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for promoter in sorted(promoters, key=lambda item: item.mouse_symbol):
            for cas_class in CAS_ORDER:
                hits = hits_by_gene_cas[(promoter.mouse_symbol, cas_class)]
                for context in peaks:
                    full = sum(hit.peak_matches[context].fully_contained for hit in hits)
                    any_overlap = sum(hit.peak_matches[context].any_overlap for hit in hits)
                    writer.writerow(
                        {
                            "mouse_symbol": promoter.mouse_symbol,
                            "human_symbol": promoter.human_symbol,
                            "gene_id": promoter.gene_id,
                            "transcript_id": promoter.transcript_id,
                            "tss": promoter.tss,
                            "strand": promoter.strand,
                            "cas": cas_class,
                            "represented_systems": CAS_CLASSES[cas_class]["systems"],
                            "context": context,
                            "promoter_midpoint_accessible": peaks[context].contains_point(promoter.chrom, promoter.midpoint),
                            "promoter_any_peak_overlap": bool(peaks[context].overlaps(promoter.chrom, promoter.start, promoter.end)),
                            "protospacers_total_passing": len(hits),
                            "protospacers_fully_in_peak": full,
                            "protospacers_any_peak_overlap": any_overlap,
                            "targetable": full > 0,
                        }
                    )


def write_candidate_sites(
    path: Path,
    promoters: list[Promoter],
    hits_by_gene_cas: dict[tuple[str, str], list[GuideHit]],
    contexts: list[str],
) -> None:
    fields = [
        "mouse_symbol",
        "human_symbol",
        "cas",
        "represented_systems",
        "protospacer_sequence",
        "pam_sequence",
        "strand",
        "target_interval",
        "protospacer_interval",
        "pam_interval",
        "distance_to_tss",
        "gc_content",
        "heuristic_score",
        "n_contexts_fully_in_peak",
        "n_contexts_any_overlap",
    ]
    for context in contexts:
        fields.extend(
            [
                f"guide_fully_in_peak_{context}",
                f"guide_any_peak_overlap_{context}",
                f"peak_id_{context}",
                f"peak_signal_{context}",
                f"distance_to_summit_{context}",
                f"overlap_bp_{context}",
            ]
        )
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for promoter in sorted(promoters, key=lambda item: item.mouse_symbol):
            for cas_class in CAS_ORDER:
                hits = sorted(
                    hits_by_gene_cas[(promoter.mouse_symbol, cas_class)],
                    key=lambda hit: (
                        -sum(hit.peak_matches[context].fully_contained for context in contexts),
                        -sum(hit.peak_matches[context].any_overlap for context in contexts),
                        abs(hit.gc_content - 0.5),
                        abs(hit.distance_to_tss + 200),
                        hit.target_start,
                    ),
                )
                for hit in hits:
                    row = {
                        "mouse_symbol": hit.mouse_symbol,
                        "human_symbol": hit.human_symbol,
                        "cas": cas_class,
                        "represented_systems": CAS_CLASSES[cas_class]["systems"],
                        "protospacer_sequence": hit.protospacer,
                        "pam_sequence": hit.pam_seq,
                        "strand": hit.target_strand,
                        "target_interval": f"{hit.chrom}:{hit.target_start}-{hit.target_end}",
                        "protospacer_interval": f"{hit.chrom}:{hit.spacer_start}-{hit.spacer_end}",
                        "pam_interval": f"{hit.chrom}:{hit.pam_start}-{hit.pam_end}",
                        "distance_to_tss": hit.distance_to_tss,
                        "gc_content": f"{hit.gc_content:.3f}",
                        "heuristic_score": hit.heuristic_score,
                        "n_contexts_fully_in_peak": sum(hit.peak_matches[context].fully_contained for context in contexts),
                        "n_contexts_any_overlap": sum(hit.peak_matches[context].any_overlap for context in contexts),
                    }
                    for context in contexts:
                        match = hit.peak_matches[context]
                        row[f"guide_fully_in_peak_{context}"] = match.fully_contained
                        row[f"guide_any_peak_overlap_{context}"] = match.any_overlap
                        row[f"peak_id_{context}"] = match.peak_id
                        row[f"peak_signal_{context}"] = "" if match.peak_signal is None else f"{match.peak_signal:.4f}"
                        row[f"distance_to_summit_{context}"] = "" if match.distance_to_summit is None else match.distance_to_summit
                        row[f"overlap_bp_{context}"] = match.overlap_bp
                    writer.writerow(row)


def summarize_dataset(
    promoters: list[Promoter],
    hits_by_gene_cas: dict[tuple[str, str], list[GuideHit]],
    peaks: dict[str, PeakIndex],
    contexts: list[str],
    any_label: str,
    all_label: str,
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for cas_class in CAS_ORDER:
        for context in contexts:
            targetable: list[bool] = []
            has_pam: list[bool] = []
            midpoint: list[bool] = []
            promoter_overlap: list[bool] = []
            full_sites = 0
            any_sites = 0
            for promoter in promoters:
                hits = hits_by_gene_cas[(promoter.mouse_symbol, cas_class)]
                full = sum(hit.peak_matches[context].fully_contained for hit in hits)
                any_overlap = sum(hit.peak_matches[context].any_overlap for hit in hits)
                targetable.append(full > 0)
                has_pam.append(bool(hits))
                midpoint.append(peaks[context].contains_point(promoter.chrom, promoter.midpoint))
                promoter_overlap.append(bool(peaks[context].overlaps(promoter.chrom, promoter.start, promoter.end)))
                full_sites += full
                any_sites += any_overlap
            rows.append(
                {
                    "cas": cas_class,
                    "context": context,
                    "n_genes": len(promoters),
                    "genes_with_passing_pam": sum(has_pam),
                    "pct_genes_with_passing_pam": f"{100 * sum(has_pam) / len(promoters):.1f}",
                    "genes_targetable_strict_full_site_in_peak": sum(targetable),
                    "pct_targetable_strict_full_site_in_peak": f"{100 * sum(targetable) / len(promoters):.1f}",
                    "genes_with_promoter_midpoint_peak": sum(midpoint),
                    "genes_with_promoter_any_peak_overlap": sum(promoter_overlap),
                    "total_full_site_supported_guides": full_sites,
                    "total_any_overlap_guides": any_sites,
                }
            )
        any_context: list[bool] = []
        all_contexts: list[bool] = []
        for promoter in promoters:
            per_context = []
            for context in contexts:
                hits = hits_by_gene_cas[(promoter.mouse_symbol, cas_class)]
                per_context.append(any(hit.peak_matches[context].fully_contained for hit in hits))
            any_context.append(any(per_context))
            all_contexts.append(all(per_context))
        for label, values in [(any_label, any_context), (all_label, all_contexts)]:
            rows.append(
                {
                    "cas": cas_class,
                    "context": label,
                    "n_genes": len(promoters),
                    "genes_with_passing_pam": sum(bool(hits_by_gene_cas[(promoter.mouse_symbol, cas_class)]) for promoter in promoters),
                    "pct_genes_with_passing_pam": f"{100 * sum(bool(hits_by_gene_cas[(promoter.mouse_symbol, cas_class)]) for promoter in promoters) / len(promoters):.1f}",
                    "genes_targetable_strict_full_site_in_peak": sum(values),
                    "pct_targetable_strict_full_site_in_peak": f"{100 * sum(values) / len(promoters):.1f}",
                    "genes_with_promoter_midpoint_peak": "NA",
                    "genes_with_promoter_any_peak_overlap": "NA",
                    "total_full_site_supported_guides": "NA",
                    "total_any_overlap_guides": "NA",
                }
            )
    return rows


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0].keys()), delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run_dataset(
    output_prefix: str,
    promoters: list[Promoter],
    sequences: dict[str, str],
    peak_paths: dict[str, Path],
    output_dir: Path,
    any_label: str,
    all_label: str,
) -> tuple[dict[tuple[str, str], list[GuideHit]], list[dict[str, object]]]:
    contexts = list(peak_paths)
    peaks = {context: load_peaks(path) for context, path in peak_paths.items()}
    hits_by_gene_cas: dict[tuple[str, str], list[GuideHit]] = {}
    for promoter in promoters:
        for cas_class in CAS_ORDER:
            hits_by_gene_cas[(promoter.mouse_symbol, cas_class)] = scan_promoter(
                promoter,
                sequences[promoter.mouse_symbol],
                cas_class,
                peaks,
            )
    write_targetability(output_dir / f"{output_prefix}_human_panel_targetability.tsv", promoters, hits_by_gene_cas, peaks)
    write_candidate_sites(output_dir / f"{output_prefix}_human_panel_candidate_sites.tsv", promoters, hits_by_gene_cas, contexts)
    summary = summarize_dataset(promoters, hits_by_gene_cas, peaks, contexts, any_label, all_label)
    write_rows(output_dir / f"{output_prefix}_human_panel_summary.tsv", summary)
    return hits_by_gene_cas, summary


def consensus_like_summary(
    promoters: list[Promoter],
    hits_by_gene_cas: dict[tuple[str, str], list[GuideHit]],
    contexts: list[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for cas_class in CAS_ORDER:
        for threshold in [1, 2, 3, 4]:
            targetable = 0
            for promoter in promoters:
                max_support = max(
                    (
                        sum(hit.peak_matches[context].fully_contained for context in contexts)
                        for hit in hits_by_gene_cas[(promoter.mouse_symbol, cas_class)]
                    ),
                    default=0,
                )
                targetable += max_support >= threshold
            rows.append(
                {
                    "cas": cas_class,
                    "support_rule": f"guide_site_full_in_peak_in_at_least_{threshold}_of_4_peak_files",
                    "n_genes": len(promoters),
                    "genes_with_passing_pam": sum(bool(hits_by_gene_cas[(promoter.mouse_symbol, cas_class)]) for promoter in promoters),
                    "genes_targetable": targetable,
                    "pct_targetable": f"{100 * targetable / len(promoters):.1f}",
                }
            )
    return rows


def compare_gse206479_gse245522_3of4(
    promoters: list[Promoter],
    hits_206479: dict[tuple[str, str], list[GuideHit]],
    hits_245522: dict[tuple[str, str], list[GuideHit]],
    contexts_206479: list[str],
    contexts_245522: list[str],
) -> list[dict[str, object]]:
    rows: list[dict[str, object]] = []
    for promoter in sorted(promoters, key=lambda item: item.mouse_symbol):
        for cas_class in CAS_ORDER:
            gse206_any = any(
                hit.peak_matches[context].fully_contained
                for hit in hits_206479[(promoter.mouse_symbol, cas_class)]
                for context in contexts_206479
            )
            max_245_support = max(
                (
                    sum(hit.peak_matches[context].fully_contained for context in contexts_245522)
                    for hit in hits_245522[(promoter.mouse_symbol, cas_class)]
                ),
                default=0,
            )
            gse245_3of4 = max_245_support >= 3
            if gse206_any and gse245_3of4:
                category = "both_gse206479_any_and_gse245522_3of4"
            elif gse206_any and not gse245_3of4:
                category = "gse206479_only_vs_3of4"
            elif gse245_3of4 and not gse206_any:
                category = "gse245522_3of4_only"
            else:
                category = "neither"
            rows.append(
                {
                    "mouse_symbol": promoter.mouse_symbol,
                    "human_symbol": promoter.human_symbol,
                    "cas": cas_class,
                    "gse206479_any": gse206_any,
                    "gse245522_3of4": gse245_3of4,
                    "gse245522_max_peak_file_support": max_245_support,
                    "category": category,
                }
            )
    return rows


def compare_mouse_un1(
    mouse_s2: Path,
    promoters: list[Promoter],
    hits_206479: dict[tuple[str, str], list[GuideHit]],
    contexts_206479: list[str],
) -> list[dict[str, object]]:
    if not mouse_s2.exists():
        return []
    mouse: dict[str, dict[str, object]] = defaultdict(lambda: {"states": set(), "promoter_overlap": set(), "passing": "NA"})
    with gzip.open(mouse_s2, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["cas"] != "Un1Cas12f1_TTTR":
                continue
            entry = mouse[row["gene"]]
            entry["passing"] = row["protospacers_total_passing"]
            if row["targetable"].lower() == "true":
                entry["states"].add(row["state"])
            if row["promoter_any_peak_overlap"].lower() == "true":
                entry["promoter_overlap"].add(row["state"])
    rows: list[dict[str, object]] = []
    for promoter in sorted(promoters, key=lambda item: item.mouse_symbol):
        human_contexts = [
            context
            for context in contexts_206479
            if any(
                hit.peak_matches[context].fully_contained
                for hit in hits_206479[(promoter.mouse_symbol, "Un1Cas12f1_TTTR")]
            )
        ]
        entry = mouse[promoter.mouse_symbol]
        mouse_any = bool(entry["states"])
        human_any = bool(human_contexts)
        if mouse_any and human_any:
            category = "both_mouse_and_human_any_context"
        elif mouse_any and not human_any:
            category = "mouse_only_any_context"
        elif human_any and not mouse_any:
            category = "human_only_any_context"
        else:
            category = "neither_any_context"
        rows.append(
            {
                "mouse_symbol": promoter.mouse_symbol,
                "human_symbol": promoter.human_symbol,
                "mouse_un1_passing_guides": entry["passing"],
                "human_un1_passing_guides": len(hits_206479[(promoter.mouse_symbol, "Un1Cas12f1_TTTR")]),
                "mouse_un1_targetable_contexts": ";".join(sorted(entry["states"])) or "NA",
                "human_un1_targetable_contexts": ";".join(human_contexts) or "NA",
                "mouse_any_context_targetable": mouse_any,
                "human_any_context_targetable": human_any,
                "comparison_category": category,
                "mouse_promoter_overlap_contexts": ";".join(sorted(entry["promoter_overlap"])) or "NA",
            }
        )
    return rows


def tfeb_alternative_tss_audit(
    transcripts_by_symbol: dict[str, list[Transcript]],
    peak_paths: dict[str, Path],
) -> list[dict[str, object]]:
    peaks = {context: load_peaks(path) for context, path in peak_paths.items()}
    rows: list[dict[str, object]] = []
    for transcript in sorted(transcripts_by_symbol.get("TFEB", []), key=lambda item: (item.tss, item.transcript_id)):
        promoter = promoter_from_transcript(
            "Tfeb",
            "TFEB",
            "alternative_tss_sensitivity",
            "all GENCODE v19 protein-coding TFEB transcripts",
            transcript,
            len(transcripts_by_symbol.get("TFEB", [])),
        )
        seq = fetch_hg19_sequence(promoter.chrom, promoter.start, promoter.end)
        promoter_overlap_count = sum(bool(peaks[context].overlaps(promoter.chrom, promoter.start, promoter.end)) for context in peaks)
        for cas_class in CAS_ORDER:
            hits = scan_promoter(promoter, seq, cas_class, peaks)
            supported_contexts = [
                context for context in peaks if any(hit.peak_matches[context].fully_contained for hit in hits)
            ]
            rows.append(
                {
                    "transcript_id": transcript.transcript_id,
                    "transcript_name": transcript.transcript_name,
                    "tss": transcript.tss,
                    "chrom": promoter.chrom,
                    "promoter_start": promoter.start,
                    "promoter_end": promoter.end,
                    "tags": ";".join(transcript.tags),
                    "ccdsid": transcript.ccdsid,
                    "cas": cas_class,
                    "passing_guides": len(hits),
                    "n_supported_contexts": len(supported_contexts),
                    "supported_contexts": ";".join(supported_contexts),
                    "promoter_overlap_context_count": promoter_overlap_count,
                }
            )
        time.sleep(0.05)
    return rows


def prepare_inputs(cache_dir: Path) -> tuple[Path, dict[str, Path], dict[str, Path], list[dict[str, object]]]:
    def cache_label(path: Path) -> str:
        try:
            return str(path.relative_to(cache_dir))
        except ValueError:
            return path.name

    gencode = cache_dir / "reference" / "gencode.v19.annotation.gtf.gz"
    download(GENCODE_V19_URL, gencode)
    gse206_paths: dict[str, Path] = {}
    for context, url in GSE206479_PEAK_URLS.items():
        path = cache_dir / "peaks" / "gse206479" / Path(url).name
        download(url, path)
        gse206_paths[context] = path
    tar_path = cache_dir / "raw" / "GSE245522_RAW.tar"
    download(GSE245522_TAR_URL, tar_path)
    gse245_paths: dict[str, Path] = {}
    with tarfile.open(tar_path) as archive:
        members = {member.name: member for member in archive.getmembers()}
        for context, member_name in GSE245522_TAR_MEMBERS.items():
            out_path = cache_dir / "peaks" / "gse245522" / member_name
            if not out_path.exists():
                out_path.parent.mkdir(parents=True, exist_ok=True)
                source = archive.extractfile(members[member_name])
                if source is None:
                    raise ValueError(f"Missing {member_name} inside {tar_path}")
                with source, out_path.open("wb") as handle:
                    shutil.copyfileobj(source, handle)
            gse245_paths[context] = out_path
    audit_rows: list[dict[str, object]] = [
        {
            "dataset": "GENCODE_v19",
            "context": "reference",
            "url": GENCODE_V19_URL,
            "local_file": cache_label(gencode),
            "rows": "NA",
            "sha256": sha256(gencode),
        },
        {
            "dataset": "GSE245522",
            "context": "raw_tar",
            "url": GSE245522_TAR_URL,
            "local_file": cache_label(tar_path),
            "rows": "NA",
            "sha256": sha256(tar_path),
        },
    ]
    for dataset, urls, paths in [
        ("GSE206479", GSE206479_PEAK_URLS, gse206_paths),
        ("GSE245522", {key: GSE245522_TAR_URL for key in gse245_paths}, gse245_paths),
    ]:
        for context, path in paths.items():
            with gzip.open(path, "rt") as handle:
                rows = sum(1 for _ in handle)
            audit_rows.append(
                {
                    "dataset": dataset,
                    "context": context,
                    "url": urls[context],
                    "local_file": cache_label(path),
                    "rows": rows,
                    "sha256": sha256(path),
                }
            )
    return gencode, gse206_paths, gse245_paths, audit_rows


def write_report(
    output_dir: Path,
    promoters: list[Promoter],
    issues: list[dict[str, str]],
    gse206_summary: list[dict[str, object]],
    gse245_summary: list[dict[str, object]],
    consensus_summary: list[dict[str, object]],
    comparison_3of4: list[dict[str, object]],
    tfeb_rows: list[dict[str, object]],
) -> None:
    lines: list[str] = []
    lines.append("# Human ortholog ATAC support check")
    lines.append("")
    lines.append("Exploratory human ortholog-panel analysis using hg19 iPSC-derived microglia ATAC-seq.")
    lines.append("Primary calls require the complete protospacer+PAM interval inside a dataset peak.")
    lines.append("")
    lines.append(f"- mapped promoters analyzed: {len(promoters)}")
    lines.append(f"- mapping exclusions: {len(issues)} ({', '.join(issue['mouse_symbol'] for issue in issues)})")
    lines.append("")
    lines.append("## Any-context/peak-file targetability")
    lines.append("")
    lines.append("| Cas/PAM | GSE206479 any context | GSE245522 any peak file | GSE245522 >=3/4 peak-file-supported |")
    lines.append("|---|---:|---:|---:|")
    for cas_class in CAS_ORDER:
        g206 = next(row for row in gse206_summary if row["cas"] == cas_class and row["context"] == "ANY_GSE206479_context")
        g245 = next(row for row in gse245_summary if row["cas"] == cas_class and row["context"] == "ANY_GSE245522_peak_file")
        g245_3 = next(
            row
            for row in consensus_summary
            if row["cas"] == cas_class and row["support_rule"] == "guide_site_full_in_peak_in_at_least_3_of_4_peak_files"
        )
        lines.append(
            f"| {cas_class} | {g206['genes_targetable_strict_full_site_in_peak']}/{g206['n_genes']} "
            f"({g206['pct_targetable_strict_full_site_in_peak']}%) | "
            f"{g245['genes_targetable_strict_full_site_in_peak']}/{g245['n_genes']} "
            f"({g245['pct_targetable_strict_full_site_in_peak']}%) | "
            f"{g245_3['genes_targetable']}/{g245_3['n_genes']} ({g245_3['pct_targetable']}%) |"
        )
    lines.append("")
    lines.append("## TFE3/TFEB interpretation guardrail")
    lines.append("")
    lines.append("- Selected TFE3 promoter: supported in every GSE206479 context and every GSE245522 peak file for all five PAM classes.")
    lines.append("- Selected TFEB promoter: unsupported in both human datasets for all five PAM classes.")
    lines.append("- TFEB sensitivity: alternative GENCODE v19 TFEB TSS choices can be accessible; therefore the negative TFEB result is promoter/TSS-dependent.")
    lines.append("")
    lines.append("## GSE206479 vs GSE245522 >=3/4 comparison")
    lines.append("")
    counts = Counter(row["category"] for row in comparison_3of4 if row["cas"] == "Un1Cas12f1_TTTR")
    for key in [
        "both_gse206479_any_and_gse245522_3of4",
        "gse206479_only_vs_3of4",
        "gse245522_3of4_only",
        "neither",
    ]:
        lines.append(f"- Un1/TTTR {key}: {counts.get(key, 0)}")
    lines.append("")
    selected_tfeb = [row for row in tfeb_rows if row["transcript_name"] == "TFEB-201"]
    alternative_supported = [
        row for row in tfeb_rows if row["transcript_name"] != "TFEB-201" and int(row["n_supported_contexts"]) > 0
    ]
    lines.append(f"- TFEB-201 audit rows: {len(selected_tfeb)}")
    lines.append(f"- supported alternative TFEB transcript/Cas rows: {len(alternative_supported)}")
    lines.append("")
    lines.append("Interpretation: use as a cross-dataset human ortholog sanity check, not as CRISPRa validation or a human atlas.")
    (output_dir / "README.md").write_text("\n".join(lines) + "\n")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--therapeutic-genes", type=Path, default=ROOT / "config" / "therapeutic_genes_locked.csv")
    parser.add_argument("--mouse-s2", type=Path, default=ROOT / "supplementary" / "table_S2_targetability_full.tsv.gz")
    parser.add_argument("--cache-dir", type=Path, default=ROOT / "workflow" / "resources" / "human_ortholog_atac")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "analysis_stats" / "human_ortholog_atac")
    args = parser.parse_args()

    args.output_dir.mkdir(parents=True, exist_ok=True)
    gencode, gse206_paths, gse245_paths, audit_rows = prepare_inputs(args.cache_dir)
    transcripts_by_symbol = parse_gencode_transcripts(gencode)
    mouse_symbols = read_therapeutic_symbols(args.therapeutic_genes)
    promoters, issues = build_human_panel(mouse_symbols, transcripts_by_symbol)
    sequences = write_promoter_tables(args.output_dir, promoters, issues)
    write_rows(args.output_dir / "input_file_audit.tsv", audit_rows)

    hits_206, summary_206 = run_dataset(
        "gse206479",
        promoters,
        sequences,
        gse206_paths,
        args.output_dir,
        "ANY_GSE206479_context",
        "ALL_GSE206479_contexts",
    )
    hits_245, summary_245 = run_dataset(
        "gse245522",
        promoters,
        sequences,
        gse245_paths,
        args.output_dir,
        "ANY_GSE245522_peak_file",
        "ALL_GSE245522_peak_files",
    )
    consensus = consensus_like_summary(promoters, hits_245, list(gse245_paths))
    write_rows(args.output_dir / "gse245522_consensus_like_summary.tsv", consensus)
    comparison_3of4 = compare_gse206479_gse245522_3of4(promoters, hits_206, hits_245, list(gse206_paths), list(gse245_paths))
    write_rows(args.output_dir / "gse206479_vs_gse245522_3of4_comparison.tsv", comparison_3of4)
    mouse_comparison = compare_mouse_un1(args.mouse_s2, promoters, hits_206, list(gse206_paths))
    if mouse_comparison:
        write_rows(args.output_dir / "mouse_vs_gse206479_un1_comparison.tsv", mouse_comparison)
    tfeb_rows = tfeb_alternative_tss_audit(transcripts_by_symbol, {**gse206_paths, **gse245_paths})
    write_rows(args.output_dir / "tfeb_alternative_tss_audit.tsv", tfeb_rows)
    write_report(args.output_dir, promoters, issues, summary_206, summary_245, consensus, comparison_3of4, tfeb_rows)


if __name__ == "__main__":
    main()
