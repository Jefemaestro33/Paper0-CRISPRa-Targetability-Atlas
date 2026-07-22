#!/usr/bin/env python3
"""Fail closed when a release-level atlas invariant is violated."""
from __future__ import annotations

import argparse
import csv
import gzip
from collections import Counter
from pathlib import Path


EXPECTED_GENES = 21_599
EXPECTED_CLASSES = {
    "Un1Cas12f1_TTTR": (20, 4),
    "SaCas9_NNGRRT": (21, 6),
    "SpCas9_NGG": (20, 3),
    "CjCas9_NNNVRYM": (22, 7),
    "Nme2Cas9_NNNNCC": (24, 6),
}
EXPECTED_SYSTEMS = {
    "Un1Cas12f1_TTTR": "HEAL;SminiCRa",
    "SaCas9_NNGRRT": "SaCas9_CRISPRa",
    "SpCas9_NGG": "SpCas9_CRISPRa",
    "CjCas9_NNNVRYM": "MiniCAFE",
    "Nme2Cas9_NNNNCC": "proposed_dNme2Cas9_activator",
}
EXPECTED_CONTEXTS = {
    "homeostatic", "PP_control", "PL_acute_LPS", "LL_tolerized",
    "sham_WT", "stroke_WT",
}


def truth(value: str) -> bool:
    return value.strip().lower() == "true"


def open_text(path: Path):
    """Open a plain or gzip-compressed release table in text mode."""
    return gzip.open(path, "rt", newline="") if path.suffix == ".gz" else path.open(newline="")


def parse_interval(value: str) -> tuple[str, int, int]:
    chrom, coordinates = value.split(":", 1)
    start, end = map(int, coordinates.split("-", 1))
    if end <= start:
        raise AssertionError(f"Invalid half-open interval: {value}")
    return chrom, start, end


def validate_atlas(path: Path) -> list[str]:
    genes: set[str] = set()
    classes: set[str] = set()
    contexts: set[str] = set()
    keys: set[tuple[str, str, str]] = set()
    therapeutic: set[str] = set()
    row_count = 0
    with open_text(path) as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            row_count += 1
            key = (row["gene"], row["cas"], row["state"])
            assert key not in keys, f"Duplicate atlas key: {key}"
            keys.add(key)
            genes.add(row["gene"])
            classes.add(row["cas"])
            contexts.add(row["state"])
            assert row["represented_systems"] == EXPECTED_SYSTEMS[row["cas"]], (
                f"Represented-system mismatch: {key}"
            )
            if truth(row["is_therapeutic"]):
                therapeutic.add(row["gene"])
            total = int(row["protospacers_total_passing"])
            complete = int(row["protospacers_fully_in_peak"])
            partial = int(row["protospacers_any_peak_overlap"])
            assert 0 <= complete <= partial <= total, f"Non-nested counts: {key}"
            assert truth(row["targetable"]) == (complete > 0), f"Targetability mismatch: {key}"
            assert not truth(row["targetable"]) or total > 0, f"Positive call without candidate: {key}"
    expected_rows = EXPECTED_GENES * len(EXPECTED_CLASSES) * len(EXPECTED_CONTEXTS)
    assert row_count == expected_rows, f"Atlas rows {row_count:,} != {expected_rows:,}"
    assert len(genes) == EXPECTED_GENES, f"Genes {len(genes):,} != {EXPECTED_GENES:,}"
    assert classes == set(EXPECTED_CLASSES), f"Unexpected classes: {classes}"
    assert contexts == EXPECTED_CONTEXTS, f"Unexpected contexts: {contexts}"
    assert len(therapeutic) == 55, f"Therapeutic genes {len(therapeutic)} != 55"
    return [
        f"atlas_rows={row_count}", f"genes={len(genes)}",
        f"classes={len(classes)}", f"contexts={len(contexts)}",
        f"therapeutic_genes={len(therapeutic)}",
    ]


def validate_candidates(path: Path) -> list[str]:
    count = 0
    ranks: Counter[tuple[str, str]] = Counter()
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle):
            count += 1
            assert row["candidate_class"] in {
                "constitutive_guide_site_candidate",
                "dataset_context_conditional_candidate",
                "limited_context_candidate",
                "no_primary_peak_support",
            }, f"Unexpected candidate class in row {count + 1}"
            cas_class = row["nuclease_pam_class"]
            spacer_length, pam_length = EXPECTED_CLASSES[cas_class]
            target = parse_interval(row["target_interval"])
            spacer = parse_interval(row["protospacer_interval"])
            pam = parse_interval(row["pam_interval"])
            assert target[0] == spacer[0] == pam[0], f"Chromosome mismatch in candidate row {count + 1}"
            assert spacer[2] - spacer[1] == spacer_length, f"Spacer length mismatch in row {count + 1}"
            assert pam[2] - pam[1] == pam_length, f"PAM length mismatch in row {count + 1}"
            assert target[1] == min(spacer[1], pam[1]) and target[2] == max(spacer[2], pam[2]), (
                f"Target interval is not spacer+PAM union in row {count + 1}"
            )
            assert target[2] - target[1] == spacer_length + pam_length, f"Non-adjacent target in row {count + 1}"
            rank = int(row["rank"])
            assert 1 <= rank <= 5
            ranks[(row["gene_symbol"], cas_class)] += 1
            if cas_class == "Un1Cas12f1_TTTR":
                assert row["pam_sequence"] in {"TTTA", "TTTG"}, f"Invalid TTTR: {row['pam_sequence']}"
                assert row["off_target_status"] == "preliminary_complete_alignment_PAM_aware_to_3_mismatches"
            else:
                assert row["off_target_status"] == "not_evaluated_requires_nuclease_specific_model"
    assert count > 0, "Candidate table is empty"
    assert max(ranks.values()) <= 5
    return [f"candidate_rows={count}", f"candidate_gene_class_groups={len(ranks)}"]


def validate_replicate_evidence(path: Path, candidate_path: Path) -> list[str]:
    with candidate_path.open(newline="") as handle:
        candidate_count = sum(1 for _ in csv.DictReader(handle))
    keys = set()
    rows = 0
    with path.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows += 1
            key = (int(row["candidate_row"]), row["run_accession"])
            assert key not in keys, f"Duplicate replicate-evidence key: {key}"
            keys.add(key)
            complete = truth(row["guide_fully_in_replicate_peak"])
            overlap = truth(row["guide_any_replicate_peak_overlap"])
            assert not complete or overlap, f"Complete support without overlap: {key}"
            assert int(row["overlap_bp"]) >= 0
    expected = candidate_count * 13
    assert rows == expected, f"Replicate-evidence rows {rows} != {expected}"
    return [f"replicate_evidence_rows={rows}"]


def validate_qc(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    run_rows = [row for row in rows if row["level"] == "run"]
    condition_rows = [row for row in rows if row["level"] == "condition"]
    assert len(run_rows) == 13, f"Run QC rows {len(run_rows)} != 13"
    assert len(condition_rows) == 6, f"Condition QC rows {len(condition_rows)} != 6"
    required = [
        "bowtie2_overall_alignment_rate_pct", "picard_percent_duplication",
        "usable_fragments_or_reads", "frip", "replicate_peak_count",
        "tss_enrichment_max",
    ]
    for row in run_rows:
        for field in required:
            assert row.get(field, "") not in {"", "NA", "nan"}, f"Missing {field} for {row['run_accession']}"
    assert {row["condition"] for row in condition_rows} == EXPECTED_CONTEXTS
    return [f"qc_run_rows={len(run_rows)}", f"qc_condition_rows={len(condition_rows)}"]


def validate_reference(selection: Path, blacklist: Path) -> list[str]:
    with selection.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    definitions = Counter(row["definition"] for row in rows)
    assert set(definitions) == {"ensembl_canonical", "appris_principal", "legacy_most5_basic"}
    assert all(value == EXPECTED_GENES for value in definitions.values()), f"TSS counts: {definitions}"
    blacklist_rows = []
    with blacklist.open() as handle:
        for line in handle:
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            assert int(fields[1]) < int(fields[2]), f"Invalid blacklist interval: {line.rstrip()}"
            blacklist_rows.append(fields)
    assert len(blacklist_rows) == 3360, f"Blacklist rows {len(blacklist_rows)} != 3360"
    return [f"tss_rows={len(rows)}", f"blacklist_mm39_rows={len(blacklist_rows)}"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
    parser.add_argument("--replicate-evidence", type=Path, required=True)
    parser.add_argument("--qc", type=Path, required=True)
    parser.add_argument("--tss-selection", type=Path, required=True)
    parser.add_argument("--blacklist", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    details = []
    details.extend(validate_atlas(args.atlas))
    details.extend(validate_candidates(args.candidates))
    details.extend(validate_replicate_evidence(args.replicate_evidence, args.candidates))
    details.extend(validate_qc(args.qc))
    details.extend(validate_reference(args.tss_selection, args.blacklist))
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("PASS\n" + "\n".join(details) + "\n")
    print(f"Release validation passed; wrote {args.output}")


if __name__ == "__main__":
    main()
