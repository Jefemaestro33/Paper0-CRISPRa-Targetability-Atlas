#!/usr/bin/env python3
"""Summarize how many targeting classes support each gene.

This addresses the reviewer request to report whether genes are supported by
only one, two, or multiple Cas/PAM targeting classes.  The script uses the
release Table S2 and writes compact count summaries for the genome-wide set and
the locked 55-gene panel.
"""
from __future__ import annotations

import argparse
import csv
import gzip
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_INPUT = ROOT / "supplementary/table_S2_targetability_full.tsv.gz"
DEFAULT_OUTPUT = ROOT / "analysis_stats/cas_multiplicity_summary.tsv"


def as_bool(value: str) -> bool:
    return str(value).strip().lower() in {"true", "1", "yes"}


def write_distribution(
    writer: csv.DictWriter,
    *,
    scope: str,
    metric: str,
    state: str,
    genes: set[str],
    supported_by_gene: dict[str, set[str]],
) -> None:
    total = len(genes)
    counts = {index: 0 for index in range(6)}
    for gene in genes:
        n_classes = len(supported_by_gene.get(gene, set()))
        counts[n_classes] = counts.get(n_classes, 0) + 1
    for n_classes in range(6):
        genes_count = counts.get(n_classes, 0)
        writer.writerow(
            {
                "scope": scope,
                "metric": metric,
                "state": state,
                "n_supported_classes": n_classes,
                "genes": genes_count,
                "total_genes": total,
                "percent_genes": f"{100 * genes_count / total:.6f}",
            }
        )


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", type=Path, default=DEFAULT_INPUT)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()

    genes_by_scope: dict[str, set[str]] = {"genome": set(), "locked_panel": set()}
    sequence_supported: dict[str, dict[str, set[str]]] = {
        "genome": defaultdict(set),
        "locked_panel": defaultdict(set),
    }
    primary_by_state: dict[str, dict[str, dict[str, set[str]]]] = {
        "genome": defaultdict(lambda: defaultdict(set)),
        "locked_panel": defaultdict(lambda: defaultdict(set)),
    }
    any_context_primary: dict[str, dict[str, set[str]]] = {
        "genome": defaultdict(set),
        "locked_panel": defaultdict(set),
    }

    opener = gzip.open if args.input.suffix == ".gz" else open
    with opener(args.input, "rt", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row["gene"]
            cas = row["cas"]
            state = row["state"]
            scopes = ["genome"]
            genes_by_scope["genome"].add(gene)
            if as_bool(row["is_therapeutic"]):
                scopes.append("locked_panel")
                genes_by_scope["locked_panel"].add(gene)

            has_sequence_candidate = int(row["protospacers_total_passing"]) > 0
            has_primary_site = as_bool(row["targetable"])
            for scope in scopes:
                if has_sequence_candidate:
                    sequence_supported[scope][gene].add(cas)
                if has_primary_site:
                    primary_by_state[scope][state][gene].add(cas)
                    any_context_primary[scope][gene].add(cas)

    args.output.parent.mkdir(parents=True, exist_ok=True)
    with args.output.open("w", newline="") as handle:
        fieldnames = [
            "scope",
            "metric",
            "state",
            "n_supported_classes",
            "genes",
            "total_genes",
            "percent_genes",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for scope in ("genome", "locked_panel"):
            write_distribution(
                writer,
                scope=scope,
                metric="sequence_candidate",
                state="sequence_window",
                genes=genes_by_scope[scope],
                supported_by_gene=sequence_supported[scope],
            )
            for state in sorted(primary_by_state[scope]):
                write_distribution(
                    writer,
                    scope=scope,
                    metric="primary_targetable",
                    state=state,
                    genes=genes_by_scope[scope],
                    supported_by_gene=primary_by_state[scope][state],
                )
            write_distribution(
                writer,
                scope=scope,
                metric="primary_targetable_any_context",
                state="any_context",
                genes=genes_by_scope[scope],
                supported_by_gene=any_context_primary[scope],
            )

    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
