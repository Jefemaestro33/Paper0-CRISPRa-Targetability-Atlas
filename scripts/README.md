# Analysis entry points

- `prepare_reference.py`: parses GENCODE vM33 and records the primary
  Ensembl-canonical TSS plus APPRIS and legacy sensitivities.
- `prepare_blacklist.py`: maps the byte-verified official ENCODE mm10
  blacklist to mm39 and records mapped and unmapped interval provenance.
- `rebuild_atlas_strict_iupac.py`: scans both strands and makes the primary
  complete-protospacer-plus-PAM-in-peak calls.
- `analysis_statistics.py`: fixed-panel resampling stability, PAM-to-guide-site
  attrition, TSS/depth/caller sensitivity, and promoter-matched null.
- `sync_analysis_outputs.py`: derives compact supplementary tables from the
  canonical genome-wide matrix.
- `within_panel_functional_associations.py`: descriptive, within-panel
  functional-category associations; it does not use an all-gene GO background.
- `summarize_results.py`: emits a readable result audit and LaTeX macros from
  final tables so the manuscript cannot silently retain obsolete numbers.
- `validate_release.py`: fails the workflow if final dimensions, nested calls,
  guide coordinates, QC coverage, TSS counts, or blacklist validity drift.
- `build_release_manifest.py`: hashes every non-ignored release file after the
  final PDF and outputs have been regenerated.
- `build_figures.sh`: rebuilds all main and supplementary figure files.

Raw-to-result orchestration and the lower-level ATAC-seq utilities are in
`workflow/Snakefile` and `workflow/scripts/`, including the normalized
candidate-by-run peak-evidence generator.
