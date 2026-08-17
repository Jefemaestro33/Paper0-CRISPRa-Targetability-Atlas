# Release v2.1.3 notes

Release v2.1.3 is the BMC Genomics major-revision package after the second
pre-submission audit of v2.1.2.

## Main changes since v2.1.2

- Added version-controlled compressed primary, replicate, and pooled ATAC-seq
  peak files under `workflow/results/peaks/`.
- Replaced hard-coded figure scope labels with generated counts from
  `analysis_stats/key_counts.json`.
- Corrected reference metadata against Crossref, including author names for
  CRISPRa/chromatin-design and microglia citations requested by reviewers.
- Clarified that GSE245522 human iPSC-derived microglia peak files represent
  input-size/preparation files rather than independent biological replicates.
- Expanded the `Tfeb`/`Tfe3` example to state explicitly that canonical `Tfeb`
  is annotation-sensitive: APPRIS and legacy TSS definitions restore support in
  five of six murine contexts.
- Added the corresponding human `TFEB` transcript-level guardrail:
  ENST00000373033.1 shares APPRIS-principal/CCDS evidence with the selected
  transcript and is supported in all eight human peak contexts for four of five
  targeting classes.
- Reported the direction and magnitude of matched-depth MACS3 sensitivity and
  stated that cross-context rank ordering is descriptive.
- Repaired missing figure/table citations and expanded the additional-file
  inventory to enumerate all supplementary tables and supplementary figures.
- Improved figure layout and keys for Figures 2, 3, 4, 6, and S7.
- Added compute-envelope notes for full raw-data rebuilds.

## Persistent identifiers

- GitHub tag: `v2.1.3`
- Zenodo concept DOI: `10.5281/zenodo.21970940`
- The version-specific Zenodo record is associated with this release after
  archival publication.
