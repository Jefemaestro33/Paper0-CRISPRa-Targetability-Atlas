# Output schema

## Table S2: genome-wide targetability matrix

Each row is one gene × targeting class × dataset context.

- `gene`, `gene_id`, `transcript_id`, `tss`, `strand`: selected reference locus.
- `tss_definition`, `selection_source`: TSS rule and any fallback.
- `cas`: one of five distinct nuclease/PAM targeting classes.
- `represented_systems`: CRISPRa configurations sharing that targeting rule.
- `state`: legacy machine-readable field name for the six dataset contexts; it
  must not be interpreted as a balanced causal state variable.
- `promoter_midpoint_accessible`: midpoint sensitivity only.
- `promoter_any_peak_overlap`: broad promoter-accessibility sensitivity only.
- `protospacers_total_passing`: sequence-filtered candidates in the window.
- `protospacers_fully_in_peak`: candidates whose complete spacer+PAM interval
  is contained in a primary peak.
- `protospacers_any_peak_overlap`: candidates with any partial peak overlap.
- `targetable`: primary operational call (`protospacers_fully_in_peak > 0`).
- `is_therapeutic`: membership in the frozen panel.

## Table S3: candidate protospacers

Each row is a ranked candidate site. Coordinates are zero-based, half-open,
GRCm39 intervals.

- `target_interval`: complete protospacer+PAM interval.
- `protospacer_interval`, `pam_interval`: component intervals.
- `distance_to_tss`: signed transcription-oriented distance of target centre.
- `guide_fully_in_peak_*`: primary context-specific guide evidence.
- `guide_any_peak_overlap_*`: partial-overlap sensitivity.
- `peak_id_*`, `peak_signal_*`, `distance_to_summit_*`: supporting peak
  provenance.
- `pam_valid_offtargets_*`: preliminary Un1Cas12f1-only exact/mismatch screen;
  blank for classes lacking a class-specific implementation.

`table_S3_replicate_evidence.tsv` is the normalized candidate-by-run companion.
It records complete and partial interval support, overlap length, run-level peak
identifier and signal, and summit distance for all 13 input runs.

## Table S4: ATAC quality and provenance

Run rows report alignment, Picard duplication, usable units, run-level peaks,
FRiP, TSS enrichment, and paired-end fragment summaries. Condition rows report
the primary peak rule and peak count. Technical runs are explicitly identified.

## Analysis statistics

- `bootstrap_cis.tsv`: fixed-panel resampling stability intervals.
- `pam_chromatin_loss.tsv`: nested sequence-to-primary-call attrition counts.
- `matched_promoter_null.tsv`: matched non-panel promoter comparison.
- `cas_multiplicity_summary.tsv`: genes supported by exactly 0--5 targeting
  classes at sequence, per-context primary-call, and any-context primary-call
  layers.
- `sensitivity_summary.tsv`: TSS/depth/caller deltas versus primary.
- `therapeutic_gene_stability.tsv`: gene-level calls under all variants.
- `within_panel_functional_associations.tsv`: Fisher tests within the fixed
  55-gene universe.
- `matched_depth_summary.tsv`: requested and realized depth, seeds, source-unit
  counts, and replicate-aware downsampling rule for every context.
- `peak_caller_concordance.tsv`: promoter-set Jaccard and binary concordance for
  matched-depth Genrich and MACS3 calls.
- `shared_peak_signal.tsv`: fragment/read counts and CPM values on the merged
  six-context primary-peak universe.

## Blacklist provenance

`reference/blacklist_liftover_provenance.tsv` records every official ENCODE
mm10 v2 interval (ENCFF543DDX), its UCSC liftOver result in mm39, and explicit
unmapped status. The tracked mm10 input has SHA-256
`9638bbeb4be8d99ddf56f1b70700f6a9336ce7f54d87032bd262465ecf3bfac7`.
