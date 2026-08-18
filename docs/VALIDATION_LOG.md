# Validation log for revision package v2.1.4

Date: 2026-08-17

This log records the local validation steps run after the v2.1.4 citation and
integrity audit. The
goal was to verify that the revised manuscript, figure code, source tables, peak
artifacts, and release checks are internally consistent before creating the
final GitHub/Zenodo release.

## Commands run

```bash
conda run -n paper0 python scripts/summarize_results.py
```

Result: PASS. Regenerated `manuscript/results_macros.tex`,
`analysis_stats/key_results.md`, and `analysis_stats/key_counts.json`.

```bash
conda run -n paper0 bash scripts/build_figures.sh
```

Result: PASS. Regenerated all main and supplementary PDF/PNG figures.

```bash
conda run -n paper0 python scripts/audit_reference_metadata.py manuscript/references.bib
```

Result: PASS. Thirty DOI-backed BibTeX entries were checked against Crossref for
title, journal, volume, pages/article number, year, and explicitly listed author
metadata. The two non-DOI entries (`genrich2019`, `krueger2015trimgalore`) were
retained as software/provenance references.

```bash
conda run -n paper0 python scripts/audit_manuscript_integrity.py --repo .
```

Result: PASS. The audit checked first-citation order for Figures 1--6 and
Supplementary Figures S1--S7, presence of Table S1--S7 citations, figure
script/output mapping for all 13 figures, Additional-file paths, and common
revision-residue patterns.

```bash
conda run -n paper0 latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Run directory: `manuscript/`.

Result: PASS. Built `manuscript/main.pdf` with resolved citations. The remaining
LaTeX messages are line-breaking warnings (`Underfull \hbox`), not missing
references or compilation errors.

```bash
conda run -n paper0 python scripts/validate_release.py \
  --atlas supplementary/table_S2_targetability_full.tsv.gz \
  --candidates supplementary/table_S3_candidate_protospacers.csv \
  --replicate-evidence supplementary/table_S3_replicate_evidence.tsv \
  --qc supplementary/table_S4_atac_qc.csv \
  --tss-selection reference/tss_selection.tsv \
  --blacklist supplementary/table_S6_blacklist_mm39_lifted.bed \
  --output analysis_stats/release_validation.txt
```

Result: PASS. The validation output was written to
`analysis_stats/release_validation.txt`.

```bash
conda run -n paper0 python -c "<direct test_core_logic.py runner>"
```

Result: 13/13 PASS. The local `paper0` environment did not expose a `pytest`
binary, so the `test_*` functions in `tests/test_core_logic.py` were imported
and executed directly.

## Direct regression tests passed

- `test_candidate_rank_signal_ignores_unsupported_contexts`
- `test_candidate_resource_does_not_label_predictions_as_recommendations`
- `test_historical_midpoint_is_minus_225_on_both_strands`
- `test_iupac_reverse_complement_for_tttr`
- `test_matched_depth_allocation_preserves_total_and_respects_capacity`
- `test_peak_match_requires_complete_guide_for_primary_call`
- `test_promoter_interval_respects_transcriptional_orientation`
- `test_reciprocal_consensus_does_not_use_transitive_chain`
- `test_release_atlas_is_gzip_compressed_and_stream_readable`
- `test_run_level_peak_evidence_distinguishes_complete_and_partial_support`
- `test_tracked_blacklists_have_verified_source_and_valid_mm39_intervals`
- `test_un1cas12f_scanning_reports_both_orientations_and_full_intervals`
- `test_workflow_temporary_paths_are_expanded_per_wildcard_job`

## Compute notes

A full raw-data rebuild requires approximately 17 GiB of compressed FASTQ
downloads, approximately 300 GB of free disk, and at least 8 vCPU. Peak calling
is the memory-sensitive step; on constrained machines, serialize the peak-heavy
rules as shown in `README.md`.
