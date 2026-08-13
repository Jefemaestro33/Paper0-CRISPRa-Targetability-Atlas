# Archival release checklist for the BMC Genomics revision

Purpose: ensure the revised package satisfies the reviewer/editor request for a
versioned, persistent, reconstructable resource.

## Current status after the 2026-08-13 VM rebuild

Completed:

- The configured Snakemake workflow completed on the Google Cloud VM
  `paper0-cpu-runner` (`e2-standard-8`, 300 GB disk) after rerunning
  memory-heavy rules with serialized thread allocation.
- Final Snakemake dry-run reported that all requested files are present and up
  to date.
- `analysis_stats/release_validation.txt` reports `PASS`.
- The ignored Figure 6 BigWig intermediates were regenerated and restored
  locally under `workflow/results/bigwig/`.
- Primary, replicate, pooled, and matched-depth peak outputs were restored
  locally for archival packaging.
- VM Snakemake logs were restored locally under `workflow/logs/`.
- The archival package and release manifest scripts now include selected
  ignored workflow products when they are present locally.

Still pending before journal resubmission:

- Upload the final package to a persistent archive and update the manuscript
  with the final release URL/DOI.

## Required before resubmission

- Create a clean commit containing the revised manuscript, scripts, source
  tables, figures, documentation, validation outputs, and manifest.
- Build a local archival package from the clean approved commit:

```bash
python3 scripts/prepare_archival_package.py
```

- Confirm ignored workflow intermediates needed to rebuild Figure 6 are present:
  `workflow/results/bigwig/*.bw`.
- Rerun the full figure build:

```bash
PYTHON_BIN=/path/to/pinned/env/bin/python ./scripts/build_figures.sh
```

- Run the full pytest suite inside the pinned environment from
  `environment.yml`:

```bash
micromamba create -f environment.yml
micromamba activate paper0-atlas
pytest -q tests
```

- Run the primary release validator:

```bash
python3 scripts/validate_release.py \
  --atlas supplementary/table_S2_targetability_full.tsv.gz \
  --candidates supplementary/table_S3_candidate_protospacers.csv \
  --replicate-evidence supplementary/table_S3_replicate_evidence.tsv \
  --qc supplementary/table_S4_atac_qc.csv \
  --tss-selection reference/tss_selection.tsv \
  --blacklist supplementary/table_S6_blacklist_mm39_lifted.bed \
  --output analysis_stats/release_validation.txt
```

- Regenerate the release manifest after the final PDF and outputs settle:

```bash
python3 scripts/build_release_manifest.py
```

## Files/classes that should be archived

### Source and workflow

- `workflow/Snakefile`
- `workflow/scripts/`
- `scripts/`
- `environment.yml`
- `config/samples.tsv`
- `config/therapeutic_genes_locked.csv`
- `reference/tss_selection.tsv`
- `reference/tss_definition_summary.tsv`
- `reference/promoters_*.bed` if present in the final source tree.

### Primary outputs

- `supplementary/table_S1_therapeutic_genes.csv`
- `supplementary/table_S2_targetability_full.tsv.gz`
- `supplementary/table_S3_candidate_protospacers.csv`
- `supplementary/table_S3_replicate_evidence.tsv`
- `supplementary/table_S3_offtarget_alignments.tsv`
- `supplementary/table_S4_atac_qc.csv`
- `supplementary/table_S5_accessibility_dynamics.csv`
- `supplementary/table_S6_blacklist_mm10_original.bed`
- `supplementary/table_S6_blacklist_mm39_lifted.bed`
- `supplementary/table_S7_statistical_tests.csv`

### Analysis statistics and revision checks

- `analysis_stats/release_validation.txt`
- `analysis_stats/key_results.md`
- `analysis_stats/cas_multiplicity_summary.tsv`
- `analysis_stats/human_ortholog_atac/`
- `analysis_stats/matched_depth_summary.tsv`
- `analysis_stats/peak_caller_concordance.tsv`
- `analysis_stats/shared_peak_signal.tsv`
- `analysis_stats/sensitivity_summary.tsv`
- `analysis_stats/within_panel_functional_associations.tsv`

### Figures and manuscript

- `figures/output/*.pdf`
- `figures/scripts/`
- `manuscript/main.tex`
- `manuscript/main.pdf`
- `manuscript/references.bib`
- `manuscript/results_macros.tex`

### Provenance and integrity

- `release_manifest.sha256`
- `LICENSE`
- `CITATION.cff` if present.
- `docs/ANALYSIS_CONTRACT.md`
- `docs/OUTPUT_SCHEMA.md`
- `docs/REVISION_AUDIT.md`
- `docs/POST_REVIEW_VALIDATION_LOG.md`
- `docs/HUMAN_ORTHOLOG_ATAC_CHECK.md`
- `docs/REVISION_DECISIONS.md`
- `docs/FIGURE_VISUAL_AUDIT.md`
- `scripts/prepare_archival_package.py`

## Inputs to archive or identify with persistent checksums

- GRCm39/mm39 FASTA: include exact source URL and checksum; include the file
  itself if archive size permits.
- GENCODE vM33 annotation and derived promoter/TSS tables.
- Primary murine ATAC peak files used for Table S2 calls.
- Matched-depth Genrich/MACS3 peak files used for sensitivity analyses.
- BigWig signal tracks used for Figure 6 or a documented full-rebuild path to
  recreate them.
- GSE206479 and GSE245522 processed human peak inputs or exact URLs and
  SHA-256 hashes from `analysis_stats/human_ortholog_atac/input_file_audit.tsv`.

## Data availability update after DOI assignment

Replace the provisional Data availability sentence:

> Source code, frozen configuration, source tables, figure code, and the
> release manifest will be deposited in a versioned GitHub/Zenodo archival
> release for the revised submission.

with the final GitHub release URL and Zenodo DOI.
