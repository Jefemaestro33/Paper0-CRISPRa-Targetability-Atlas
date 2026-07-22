#!/usr/bin/env python3
"""Build auditable CRISPRa promoter definitions from a GENCODE GTF.

The primary definition uses the transcript tagged ``Ensembl_canonical``.
Two prespecified sensitivities are emitted: APPRIS principal and the legacy
most-5-prime transcript among GENCODE ``basic`` transcripts.  All choices and
fallbacks are recorded in a transcript-selection table.
"""
from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path


ATTR_RE = re.compile(r'([A-Za-z0-9_]+) "([^"]*)";')
APPRIS_RE = re.compile(r"appris_principal(?:_(\d+))?$")


@dataclass(frozen=True)
class Transcript:
    chrom: str
    start: int
    end: int
    strand: str
    gene_id: str
    gene_name: str
    transcript_id: str
    transcript_name: str
    tags: tuple[str, ...]
    support_level: str

    @property
    def tss(self) -> int:
        # BED-compatible zero-based coordinate of the first transcribed base.
        return self.start if self.strand == "+" else self.end - 1

    @property
    def length(self) -> int:
        return self.end - self.start

    @property
    def appris_rank(self) -> int:
        ranks = []
        for tag in self.tags:
            match = APPRIS_RE.match(tag)
            if match:
                ranks.append(int(match.group(1) or 1))
        return min(ranks) if ranks else 99

    @property
    def tsl_rank(self) -> int:
        value = self.support_level.split()[0]
        return int(value) if value.isdigit() else 99


def parse_attrs(text: str) -> tuple[dict[str, str], tuple[str, ...]]:
    pairs = ATTR_RE.findall(text)
    attrs: dict[str, str] = {}
    tags: list[str] = []
    for key, value in pairs:
        if key == "tag":
            tags.append(value)
        else:
            attrs[key] = value
    return attrs, tuple(tags)


def read_transcripts(gtf: Path, excluded_chromosomes: set[str]) -> dict[str, list[Transcript]]:
    by_gene: dict[str, list[Transcript]] = defaultdict(list)
    with gtf.open() as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "transcript":
                continue
            if fields[0] in excluded_chromosomes:
                continue
            attrs, tags = parse_attrs(fields[8])
            if attrs.get("gene_type") != "protein_coding":
                continue
            transcript = Transcript(
                chrom=fields[0],
                start=int(fields[3]) - 1,
                end=int(fields[4]),
                strand=fields[6],
                gene_id=attrs["gene_id"],
                gene_name=attrs["gene_name"],
                transcript_id=attrs["transcript_id"],
                transcript_name=attrs.get("transcript_name", ""),
                tags=tags,
                support_level=attrs.get("transcript_support_level", "NA"),
            )
            by_gene[transcript.gene_name].append(transcript)
    return by_gene


def canonical_sort_key(tx: Transcript) -> tuple:
    return (
        0 if "Ensembl_canonical" in tx.tags else 1,
        tx.appris_rank,
        tx.tsl_rank,
        -tx.length,
        tx.transcript_id,
    )


def appris_sort_key(tx: Transcript) -> tuple:
    return (
        tx.appris_rank,
        0 if "Ensembl_canonical" in tx.tags else 1,
        tx.tsl_rank,
        -tx.length,
        tx.transcript_id,
    )


def legacy_sort_key(tx: Transcript) -> tuple:
    # Most 5-prime genomic TSS in transcriptional orientation.
    oriented_tss = tx.tss if tx.strand == "+" else -tx.tss
    return (oriented_tss, tx.tsl_rank, -tx.length, tx.transcript_id)


def choose_transcript(transcripts: list[Transcript], definition: str) -> tuple[Transcript, str]:
    if definition == "ensembl_canonical":
        candidates = [tx for tx in transcripts if "Ensembl_canonical" in tx.tags]
        if candidates:
            return sorted(candidates, key=canonical_sort_key)[0], "Ensembl_canonical"
        candidates = [tx for tx in transcripts if tx.appris_rank < 99]
        if candidates:
            return sorted(candidates, key=appris_sort_key)[0], "fallback_APPRIS"
        candidates = [tx for tx in transcripts if "basic" in tx.tags] or transcripts
        return sorted(candidates, key=canonical_sort_key)[0], "fallback_basic"
    if definition == "appris_principal":
        candidates = [tx for tx in transcripts if tx.appris_rank < 99]
        if candidates:
            return sorted(candidates, key=appris_sort_key)[0], "APPRIS_principal"
        candidates = [tx for tx in transcripts if "Ensembl_canonical" in tx.tags]
        if candidates:
            return sorted(candidates, key=canonical_sort_key)[0], "fallback_Ensembl_canonical"
        candidates = [tx for tx in transcripts if "basic" in tx.tags] or transcripts
        return sorted(candidates, key=canonical_sort_key)[0], "fallback_basic"
    if definition == "legacy_most5_basic":
        candidates = [tx for tx in transcripts if "basic" in tx.tags] or transcripts
        return sorted(candidates, key=legacy_sort_key)[0], "legacy_most5_basic"
    raise ValueError(f"Unknown TSS definition: {definition}")


def promoter_interval(tx: Transcript, upstream_far: int, upstream_near: int) -> tuple[int, int]:
    if tx.strand == "+":
        start, end = tx.tss - upstream_far, tx.tss - upstream_near
    else:
        start, end = tx.tss + upstream_near + 1, tx.tss + upstream_far + 1
    return max(0, start), max(0, end)


def write_outputs(by_gene: dict[str, list[Transcript]], outdir: Path) -> None:
    outdir.mkdir(parents=True, exist_ok=True)
    definitions = ("ensembl_canonical", "appris_principal", "legacy_most5_basic")
    selected: dict[tuple[str, str], tuple[Transcript, str]] = {}
    for gene, transcripts in by_gene.items():
        for definition in definitions:
            selected[(gene, definition)] = choose_transcript(transcripts, definition)

    selection_path = outdir / "tss_selection.tsv"
    with selection_path.open("w", newline="") as handle:
        fields = [
            "gene_symbol", "gene_id", "definition", "selection_source", "chrom", "tss",
            "strand", "transcript_id", "transcript_name", "transcript_support_level", "tags",
            "n_annotated_transcripts",
        ]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for gene in sorted(by_gene):
            for definition in definitions:
                tx, source = selected[(gene, definition)]
                writer.writerow({
                    "gene_symbol": gene,
                    "gene_id": tx.gene_id,
                    "definition": definition,
                    "selection_source": source,
                    "chrom": tx.chrom,
                    "tss": tx.tss,
                    "strand": tx.strand,
                    "transcript_id": tx.transcript_id,
                    "transcript_name": tx.transcript_name,
                    "transcript_support_level": tx.support_level,
                    "tags": ";".join(tx.tags),
                    "n_annotated_transcripts": len(by_gene[gene]),
                })

    for definition in definitions:
        bed_path = outdir / f"promoters_{definition}.bed"
        with bed_path.open("w") as handle:
            for gene in sorted(by_gene):
                tx, source = selected[(gene, definition)]
                start, end = promoter_interval(tx, 400, 50)
                fields = [
                    tx.chrom, str(start), str(end), gene, "0", tx.strand, str(tx.tss),
                    tx.gene_id, tx.transcript_id, definition, source,
                ]
                handle.write("\t".join(fields) + "\n")

    summary_path = outdir / "tss_definition_summary.tsv"
    with summary_path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["comparison", "same_tss", "different_tss", "total_genes"])
        primary = {gene: selected[(gene, "ensembl_canonical")][0].tss for gene in by_gene}
        for definition in ("appris_principal", "legacy_most5_basic"):
            same = sum(primary[gene] == selected[(gene, definition)][0].tss for gene in by_gene)
            writer.writerow([f"ensembl_canonical_vs_{definition}", same, len(by_gene) - same, len(by_gene)])


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gtf", type=Path, required=True)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument(
        "--exclude-chrom",
        action="append",
        default=["chrM"],
        help="Chromosome to exclude; may be supplied more than once (default: chrM).",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    by_gene = read_transcripts(args.gtf, set(args.exclude_chrom))
    if not by_gene:
        raise SystemExit("No protein-coding transcripts were parsed from the GTF")
    write_outputs(by_gene, args.outdir)
    print(f"Prepared TSS definitions for {len(by_gene):,} protein-coding genes")


if __name__ == "__main__":
    main()
