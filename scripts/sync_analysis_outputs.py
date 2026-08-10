#!/usr/bin/env python3
"""Synchronize compact supplementary tables from corrected targetability outputs."""
from __future__ import annotations

import csv
from collections import Counter
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "supplementary"
STATS = ROOT / "analysis_stats"
TABLE_S1 = SUPP / "table_S1_therapeutic_genes.csv"
LOCKED_PANEL = ROOT / "config/therapeutic_genes_locked.csv"
TABLE_S2 = SUPP / "table_S2_targetability_full.tsv.gz"
TABLE_S5 = SUPP / "table_S5_accessibility_dynamics.csv"
TABLE_S7 = SUPP / "table_S7_statistical_tests.csv"
PRIMARY_CLASS = "Un1Cas12f1_TTTR"
STATES = [
    "homeostatic", "PP_control", "PL_acute_LPS", "LL_tolerized",
    "sham_WT", "stroke_WT",
]


def truth(value) -> bool:
    return str(value).lower() == "true"


def classify(calls: dict[str, bool]) -> tuple[str, int]:
    study_support = sum([
        calls["homeostatic"],
        any(calls[state] for state in ("PP_control", "PL_acute_LPS", "LL_tolerized")),
        any(calls[state] for state in ("sham_WT", "stroke_WT")),
    ])
    if all(calls.values()):
        return "all_six_contexts", study_support
    if not any(calls.values()):
        return "no_surveyed_context", study_support
    if study_support >= 2:
        return "multi_study_support", study_support
    return "single_study_support", study_support


def note(pattern: str, calls: dict[str, bool]) -> str:
    open_states = [state for state in STATES if calls[state]]
    if pattern == "all_six_contexts":
        return "At least one complete Un1Cas12f1 protospacer+PAM site was contained in a primary ATAC peak in all six surveyed contexts."
    if pattern == "no_surveyed_context":
        return "No complete Un1Cas12f1 protospacer+PAM site was contained in a primary ATAC peak in the surveyed contexts; this is not evidence of functional non-activatability."
    return "Guide-site support was observed in: " + "; ".join(open_states) + ". Cross-study differences are descriptive and may reflect biology, protocol, depth, or cell composition."


def update_s1(patterns: dict[str, str]) -> tuple[list[str], dict[str, dict]]:
    with LOCKED_PANEL.open(newline="") as handle:
        original = list(csv.DictReader(handle))
    fields = ["gene_symbol", "gene_name", "category", "priority", "justification", "guide_site_support_pattern"]
    with TABLE_S1.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in original:
            gene = row["gene_symbol"]
            writer.writerow({
                "gene_symbol": gene, "gene_name": row.get("gene_name", ""),
                "category": row.get("category", ""), "priority": row.get("priority", ""),
                "justification": row.get("justification", ""),
                "guide_site_support_pattern": patterns[gene],
            })
    rows = {row["gene_symbol"]: row for row in original}
    return [row["gene_symbol"] for row in original], rows


def write_s5(atlas: pd.DataFrame, genes: list[str], metadata: dict[str, dict]) -> dict[str, str]:
    subset = atlas[(atlas["cas"] == PRIMARY_CLASS) & atlas["gene"].isin(genes)].copy()
    by_key = {(row.gene, row.state): row for row in subset.itertuples()}
    patterns: dict[str, str] = {}
    rows = []
    for gene in genes:
        calls = {state: bool(by_key[(gene, state)].targetable) for state in STATES}
        any_overlap = {state: bool(by_key[(gene, state)].promoter_any_peak_overlap) for state in STATES}
        midpoint = {state: bool(by_key[(gene, state)].promoter_midpoint_accessible) for state in STATES}
        pattern, studies = classify(calls)
        patterns[gene] = pattern
        row = {
            "gene_symbol": gene, "category": metadata[gene].get("category", ""),
            "guide_site_support_pattern": pattern,
            "n_contexts_guide_targetable": sum(calls.values()), "n_studies_with_support": studies,
            "notes": note(pattern, calls),
        }
        for state in STATES:
            row[f"guide_targetable_{state}"] = calls[state]
            row[f"promoter_any_peak_overlap_{state}"] = any_overlap[state]
            row[f"promoter_midpoint_accessible_{state}"] = midpoint[state]
        rows.append(row)
    pd.DataFrame(rows).to_csv(TABLE_S5, index=False)
    return patterns


def write_s7() -> None:
    rows: list[dict] = []
    sources = [
        ("panel_resampling_stability", STATS / "bootstrap_cis.tsv"),
        ("pam_to_guide_site_loss", STATS / "pam_chromatin_loss.tsv"),
        ("matched_promoter_null", STATS / "matched_promoter_null.tsv"),
        ("analysis_sensitivity", STATS / "sensitivity_summary.tsv"),
        ("cas_class_multiplicity", STATS / "cas_multiplicity_summary.tsv"),
    ]
    for analysis, path in sources:
        if not path.exists():
            continue
        frame = pd.read_csv(path, sep="\t")
        for record in frame.to_dict(orient="records"):
            rows.append({"analysis": analysis, "record": len(rows) + 1, "details": "; ".join(f"{key}={value}" for key, value in record.items())})
    pd.DataFrame(rows, columns=["analysis", "record", "details"]).to_csv(TABLE_S7, index=False)


def main() -> None:
    atlas = pd.read_csv(TABLE_S2, sep="\t", low_memory=False)
    for column in ("targetable", "promoter_any_peak_overlap", "promoter_midpoint_accessible"):
        atlas[column] = atlas[column].map(truth)
    with LOCKED_PANEL.open(newline="") as handle:
        original = list(csv.DictReader(handle))
    genes = [row["gene_symbol"] for row in original]
    metadata = {row["gene_symbol"]: row for row in original}
    patterns = write_s5(atlas, genes, metadata)
    update_s1(patterns)
    write_s7()
    counts = Counter(patterns.values())
    print("Updated Tables S1, S5, and S7 from guide-site-aware primary calls")
    for label, count in sorted(counts.items()):
        print(f"  {label}: {count}")


if __name__ == "__main__":
    main()
