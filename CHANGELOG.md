# Changelog

## 2.1.3 — major-revision response hardening

- Added compressed primary, replicate, and pooled ATAC peak files to the
  version-controlled resource so the primary calls can be inspected from a
  clone without unpacking the full archival tarball.
- Replaced hard-coded figure scope labels with generated release counts and
  added a tracked `analysis_stats/key_counts.json` source.
- Corrected reference metadata against Crossref, including author fields.
- Reworked figure layout, legends, and encodings for Figures 2, 3, 4, 6, and
  Supplementary Figure S7.
- Expanded manuscript reporting for batch/caller sensitivity, TSS annotation
  dependence at Tfeb/TFEB, same-study sham--stroke interpretation, human
  GSE245522 peak-file provenance, panel curation, additional-file declarations,
  compute requirements, and data availability.
- Updated release validation documentation and regenerated manuscript, figures,
  validation output, and SHA-256 manifest.

## 2.1.2 — manuscript accountability wording cleanup

- Refined the manuscript reproducibility/accountability wording.
- Rebuilt the manuscript PDF and SHA-256 release manifest.
- Repackaged the archival release from the updated public source tree.

## 2.1.1 — public-tree hygiene release

- Removed internal revision-planning and response-draft documents from the
  public source tree while preserving reproducibility-facing documentation.
- Corrected Zenodo creator ORCID metadata.
- Regenerated the SHA-256 release manifest after the public-tree cleanup.

## 2.0.0 — guide-site-aware major revision

- Replaced promoter-midpoint targetability with complete
  protospacer-plus-PAM containment in a primary ATAC peak.
- Rebuilt the analysis from 13 checksum-verified ENA runs, with retry-on-HTTP-
  error downloads and an MD5 gate before any FASTQ enters preprocessing.
- Added biological-replicate consensus peaks and explicit technical-run labels.
- Added pooled-library duplicate removal across technical lanes without
  deduplicating across independent biological replicates.
- Added Ensembl-canonical, APPRIS, and legacy TSS analyses.
- Added matched-depth Genrich/MACS3 sensitivities and continuous shared-universe
  peak quantification.
- Matched-depth analyses preserve biological replicate structure: runs are
  downsampled separately and consensus peaks are rebuilt. Technical-only
  contexts remain pooled. MACS3 uses native paired-fragment intervals for
  paired-end libraries and shift/extension only for single-end libraries.
- Corrected the comparison to five nuclease/PAM classes and six CRISPRa system
  configurations.
- Replaced genome-wide shuffled peaks and all-gene GO testing with a
  promoter-matched null and within-panel descriptive associations.
- Removed inherited result-derived accessibility labels from the frozen panel
  input so the locked gene metadata no longer encodes the old classification.
- Added guide-specific peak provenance and a preliminary PAM-aware Un1Cas12f1
  off-target screen.
- Added normalized candidate-by-run evidence so biological-replicate and
  technical-lane peak support can be inspected independently.
- Isolated every wildcard-dependent temporary BAM by run or condition and
  added a regression test against concurrent-path collisions.
- Standardized predictive output terminology on candidate protospacers and
  `candidate_class`; no output is labelled as an experimental recommendation.
- Replaced an inherited invalid mm39 blacklist lift with a validated UCSC
  liftOver of the byte-verified ENCODE mm10 v2 source and recorded all 75
  unmapped intervals.
- Made paired-end blacklist filtering fragment-aware and declared the FASTA
  index as an explicit workflow product.
- Disabled Picard optical-duplicate parsing because archived SRR read names do
  not retain flowcell coordinates; coordinate-based duplicate removal remains.
- Corrected a one-base negative-strand offset in the historical midpoint
  sensitivity so it represents transcription-oriented -225 on both strands.
- Rebuilt the manuscript, figures, supplementary tables, tests, environment,
  license, and complete workflow documentation.
- Added a fail-closed release validator for atlas dimensions, nested calls,
  guide coordinates, TTTR identity, QC completeness, TSS counts, and blacklist
  validity.
- Added a SHA-256 manifest for the complete release package.
- Stored the 647,970-row genome-wide atlas as deterministic gzip
  (`table_S2_targetability_full.tsv.gz`) so the complete release remains below
  GitHub's per-file size limit without Git LFS.
