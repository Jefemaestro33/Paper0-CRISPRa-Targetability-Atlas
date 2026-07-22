#!/usr/bin/env python3
"""Statistical and robustness analyses for the guide-site-aware atlas.

The curated 55-gene panel is the complete panel of interest, not a random
sample of all therapeutic genes.  Bootstrap intervals are therefore labelled
as panel-resampling stability intervals.  A separate matched-promoter null
compares the panel with non-panel promoters matched on promoter GC, annotated
transcript count, and Un1Cas12f1 passing-protospacer count.
"""
from __future__ import annotations

import argparse
import csv
import math
import subprocess
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import pysam


ROOT = Path(__file__).resolve().parents[1]
STATES = [
    "homeostatic", "PP_control", "PL_acute_LPS", "LL_tolerized",
    "sham_WT", "stroke_WT",
]
PRIMARY_CLASS = "Un1Cas12f1_TTTR"
VARIANTS = {
    "primary_reproducible_peaks_canonical_tss": ROOT / "supplementary/table_S2_targetability_full.tsv.gz",
    "matched_depth_genrich": ROOT / "workflow/results/sensitivity/matched_depth_genrich_atlas.tsv",
    "matched_depth_macs3": ROOT / "workflow/results/sensitivity/matched_depth_macs3_atlas.tsv",
    "appris_principal_tss": ROOT / "workflow/results/sensitivity/tss/appris_principal_atlas.tsv",
    "legacy_most5_basic_tss": ROOT / "workflow/results/sensitivity/tss/legacy_most5_basic_atlas.tsv",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=VARIANTS["primary_reproducible_peaks_canonical_tss"])
    parser.add_argument("--genes", type=Path, default=ROOT / "config/therapeutic_genes_locked.csv")
    parser.add_argument("--promoters", type=Path, default=ROOT / "reference/promoters_ensembl_canonical.bed")
    parser.add_argument("--tss-selection", type=Path, default=ROOT / "reference/tss_selection.tsv")
    parser.add_argument("--fasta", type=Path, default=ROOT / "workflow/resources/mm39.fa")
    parser.add_argument("--outdir", type=Path, default=ROOT / "analysis_stats")
    parser.add_argument("--iterations", type=int, default=10_000)
    parser.add_argument("--seed", type=int, default=1729)
    return parser.parse_args()


def read_panel(path: Path) -> list[str]:
    with path.open(newline="") as handle:
        return [row["gene_symbol"] for row in csv.DictReader(handle)]


def read_atlas(path: Path) -> pd.DataFrame:
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    for column in ("targetable", "is_therapeutic", "promoter_midpoint_accessible", "promoter_any_peak_overlap"):
        if column in frame:
            frame[column] = frame[column].astype(str).str.lower().eq("true")
    return frame


def percentile_interval(values: np.ndarray) -> tuple[float, float]:
    return tuple(np.quantile(values, [0.025, 0.975]).tolist())


def panel_summaries(frame: pd.DataFrame, panel: list[str], outdir: Path, iterations: int, rng: np.random.Generator) -> None:
    subset = frame[frame["gene"].isin(panel)].copy()
    rows: list[dict] = []
    loss_rows: list[dict] = []
    for (cas, state), group in subset.groupby(["cas", "state"], sort=False):
        group = group.drop_duplicates("gene").set_index("gene").reindex(panel)
        values = group["targetable"].fillna(False).to_numpy(dtype=bool)
        resamples = rng.choice(values.astype(float), size=(iterations, len(panel)), replace=True).mean(axis=1)
        lower, upper = percentile_interval(resamples)
        n_targetable = int(values.sum())
        rows.append({
            "cas": cas, "state": state, "n_targetable": n_targetable,
            "n_panel": len(panel), "proportion": n_targetable / len(panel),
            "stability_interval_lower": lower, "stability_interval_upper": upper,
            "resamples": iterations, "seed": rng.bit_generator._seed_seq.entropy,
            "interpretation": "percentile interval from resampling the fixed curated panel; not population inference",
        })
        has_pam = group["protospacers_total_passing"].fillna(0).to_numpy(dtype=int) > 0
        both = int((has_pam & values).sum())
        pam_only = int((has_pam & ~values).sum())
        impossible = int((~has_pam & values).sum())
        absent = int((~has_pam & ~values).sum())
        loss_rows.append({
            "cas": cas, "state": state, "pam_bearing_and_targetable": both,
            "pam_bearing_not_targetable": pam_only,
            "targetable_without_pam_check": impossible, "pam_absent": absent,
            "loss_fraction_among_pam_bearing": pam_only / (both + pam_only) if both + pam_only else math.nan,
        })
    pd.DataFrame(rows).to_csv(outdir / "bootstrap_cis.tsv", sep="\t", index=False)
    pd.DataFrame(loss_rows).to_csv(outdir / "pam_chromatin_loss.tsv", sep="\t", index=False)


def promoter_gc(promoters: Path, fasta: Path) -> dict[str, float]:
    if not Path(str(fasta) + ".fai").exists():
        subprocess.run(["samtools", "faidx", str(fasta)], check=True)
    genome = pysam.FastaFile(str(fasta))
    values: dict[str, float] = {}
    with promoters.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) < 4:
                continue
            sequence = genome.fetch(fields[0], int(fields[1]), int(fields[2])).upper()
            called = [base for base in sequence if base in "ACGT"]
            values[fields[3]] = (sum(base in "GC" for base in called) / len(called)) if called else math.nan
    genome.close()
    return values


def matched_promoter_null(
    frame: pd.DataFrame,
    panel: list[str],
    promoters: Path,
    selection: Path,
    fasta: Path,
    outdir: Path,
    iterations: int,
    rng: np.random.Generator,
) -> None:
    if not fasta.exists():
        print(f"Matched-promoter null skipped: FASTA not found at {fasta}")
        return
    primary = frame[(frame["cas"] == PRIMARY_CLASS) & (frame["state"] == STATES[0])].drop_duplicates("gene")
    features = primary[["gene", "protospacers_total_passing"]].rename(columns={"protospacers_total_passing": "pam_count"})
    tss = pd.read_csv(selection, sep="\t")
    tss = tss[tss["definition"] == "ensembl_canonical"][["gene_symbol", "n_annotated_transcripts"]]
    tss = tss.rename(columns={"gene_symbol": "gene", "n_annotated_transcripts": "transcript_count"})
    features = features.merge(tss, on="gene", how="left")
    gc = promoter_gc(promoters, fasta)
    features["gc"] = features["gene"].map(gc)
    features = features.dropna().reset_index(drop=True)
    features["gc_bin"] = pd.qcut(features["gc"], 10, labels=False, duplicates="drop")
    features["tx_bin"] = pd.cut(features["transcript_count"], [0, 1, 2, 4, 8, np.inf], labels=False, include_lowest=True)
    features["pam_bin"] = pd.cut(features["pam_count"], [-1, 0, 2, 5, 9, np.inf], labels=False)
    features["stratum"] = list(zip(features["gc_bin"], features["tx_bin"], features["pam_bin"]))

    panel_set = set(panel)
    background = features[~features["gene"].isin(panel_set)].copy().reset_index(drop=True)
    panel_features = features[features["gene"].isin(panel_set)].set_index("gene").reindex(panel)
    if panel_features.isna().any(axis=None):
        missing = panel_features[panel_features["gc"].isna()].index.tolist()
        raise ValueError(f"Missing matching features for panel genes: {missing}")

    exact = defaultdict(list)
    for index, row in background.iterrows():
        exact[row["stratum"]].append(index)
    scaled = background[["gc", "transcript_count", "pam_count"]].to_numpy(float)
    scales = np.nanstd(scaled, axis=0)
    scales[scales == 0] = 1
    candidate_lists: list[np.ndarray] = []
    for _, row in panel_features.iterrows():
        candidates = exact.get(row["stratum"], [])
        if len(candidates) < 25:
            target = row[["gc", "transcript_count", "pam_count"]].to_numpy(float)
            distance = np.sqrt((((scaled - target) / scales) ** 2).sum(axis=1))
            candidates = np.argsort(distance)[:200].tolist()
        candidate_lists.append(np.asarray(candidates, dtype=int))

    state_matrix = (
        frame[frame["cas"] == PRIMARY_CLASS]
        .pivot_table(index="gene", columns="state", values="targetable", aggfunc="first")
        .reindex(columns=STATES).fillna(False)
    )
    background_genes = background["gene"].to_numpy()
    null = np.zeros((iterations, len(STATES)), dtype=float)
    for iteration in range(iterations):
        selected_indices: list[int] = []
        used: set[int] = set()
        for candidates in candidate_lists:
            available = candidates[~np.isin(candidates, list(used))]
            source = available if len(available) else candidates
            choice = int(rng.choice(source))
            selected_indices.append(choice)
            used.add(choice)
        genes = background.loc[selected_indices, "gene"].tolist()
        null[iteration] = state_matrix.reindex(genes).fillna(False).to_numpy(dtype=float).mean(axis=0)

    observed = state_matrix.reindex(panel).fillna(False).to_numpy(dtype=float).mean(axis=0)
    rows = []
    for position, state in enumerate(STATES):
        values = null[:, position]
        lower, upper = percentile_interval(values)
        upper_p = (int((values >= observed[position]).sum()) + 1) / (iterations + 1)
        lower_p = (int((values <= observed[position]).sum()) + 1) / (iterations + 1)
        rows.append({
            "state": state, "observed_panel_proportion": observed[position],
            "matched_null_mean": values.mean(), "matched_null_lower_2.5pct": lower,
            "matched_null_upper_97.5pct": upper, "empirical_p_upper": upper_p,
            "empirical_p_lower": lower_p, "empirical_p_two_sided": min(1.0, 2 * min(upper_p, lower_p)),
            "iterations": iterations, "matching_variables": "promoter_GC_decile;annotated_transcript_count_bin;Un1Cas12f1_passing_protospacer_count_bin",
            "limitations": "not matched on expression or mappability because comparable measurements were unavailable for all promoters",
        })
    pd.DataFrame(rows).to_csv(outdir / "matched_promoter_null.tsv", sep="\t", index=False)


def sensitivity_summaries(primary: pd.DataFrame, panel: list[str], outdir: Path) -> None:
    frames: dict[str, pd.DataFrame] = {"primary_reproducible_peaks_canonical_tss": primary}
    for label, path in VARIANTS.items():
        if label not in frames and path.exists():
            frames[label] = read_atlas(path)
    primary_lookup = primary.set_index(["gene", "cas", "state"])["targetable"]
    summary_rows = []
    stability_rows = []
    for label, frame in frames.items():
        for scope, scoped in (("therapeutic_panel", frame[frame["gene"].isin(panel)]), ("genome_wide", frame)):
            for (cas, state), group in scoped.groupby(["cas", "state"], sort=False):
                group = group.drop_duplicates("gene")
                n = len(group)
                n_targetable = int(group["targetable"].sum())
                if label == "primary_reproducible_peaks_canonical_tss":
                    changed = 0
                    delta = 0.0
                else:
                    indexed = group.set_index(["gene", "cas", "state"])["targetable"]
                    common = indexed.index.intersection(primary_lookup.index)
                    changed = int((indexed.loc[common] != primary_lookup.loc[common]).sum())
                    primary_prop = float(primary_lookup.loc[common].mean()) if len(common) else math.nan
                    delta = n_targetable / n - primary_prop if n else math.nan
                summary_rows.append({
                    "analysis_variant": label, "scope": scope, "cas": cas, "state": state,
                    "n_targetable": n_targetable, "n_total": n,
                    "proportion": n_targetable / n if n else math.nan,
                    "delta_vs_primary": delta, "n_changed_vs_primary": changed,
                })
    for gene in panel:
        for state in STATES:
            calls = {}
            for label, frame in frames.items():
                subset = frame[(frame["gene"] == gene) & (frame["cas"] == PRIMARY_CLASS) & (frame["state"] == state)]
                calls[label] = bool(subset.iloc[0]["targetable"]) if len(subset) else False
            stability_rows.append({
                "gene": gene, "state": state, **calls,
                "n_variants_targetable": sum(calls.values()), "n_variants": len(calls),
                "unanimous": len(set(calls.values())) == 1,
            })
    pd.DataFrame(summary_rows).to_csv(outdir / "sensitivity_summary.tsv", sep="\t", index=False)
    pd.DataFrame(stability_rows).to_csv(outdir / "therapeutic_gene_stability.tsv", sep="\t", index=False)


def main() -> None:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    panel = read_panel(args.genes)
    atlas = read_atlas(args.atlas)
    rng = np.random.default_rng(args.seed)
    panel_summaries(atlas, panel, args.outdir, args.iterations, rng)
    sensitivity_summaries(atlas, panel, args.outdir)
    matched_promoter_null(atlas, panel, args.promoters, args.tss_selection, args.fasta, args.outdir, args.iterations, rng)
    print(f"Wrote corrected statistical outputs to {args.outdir}")


if __name__ == "__main__":
    main()
