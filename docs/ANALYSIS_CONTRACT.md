# Analysis contract for revision 2.0

This document records the inferential and operational choices used for the
guide-site-aware reanalysis.

## Primary analysis

- Reference: GRCm39/mm39 and GENCODE mouse vM33.
- Gene universe: nuclear protein-coding genes with a selected transcript.
- TSS: `Ensembl_canonical`; deterministic APPRIS/basic fallbacks are recorded.
- Window: transcription-oriented −400 to −50 bp from the selected TSS.
- Primary site call: the complete protospacer plus PAM interval is contained in
  a primary ATAC peak.
- Biological peak support: qualifying peak overlap in at least two independent
  replicates, with 50% reciprocal overlap and pairwise intersection.
- Technical-only inputs: sequencing lanes from the same source preparation are
  pooled with a common library identity, deduplicated across lanes, and
  labelled without a biological-replication claim.
- Frozen focused universe: the 55 genes in
  `config/therapeutic_genes_locked.csv`.

## Prespecified sensitivities

- APPRIS-principal TSS.
- Legacy most-5-prime GENCODE-basic TSS.
- Any guide/peak overlap.
- Promoter midpoint covered by a peak.
- Any promoter-window/peak overlap.
- Replicate-aware deterministic matched depth with Genrich.
- Replicate-aware deterministic matched depth with MACS3.
- Shared-universe continuous peak signal.

## Interpretive constraints

- Six columns are dataset contexts, not a balanced causal state experiment.
- ATAC support is a prioritization feature, not a universal requirement.
- A negative resource call is not evidence that a locus cannot be activated.
- Panel resampling measures stability of the fixed panel, not inference to all
  possible therapeutic genes.
- Candidate rankings and off-target counts require experimental validation.

Random seed for resampling, matched panels, and downsampling: 1729, with
deterministic per-condition offsets where stated in source tables.

## Post-review exploratory human ortholog check

The human ortholog-panel analysis in `analysis_stats/human_ortholog_atac/` is
not part of the primary murine resource definition.  It is a revision-response
robustness check using GENCODE v19/hg19 promoters and public iPSC-/hPSC-derived
microglia ATAC peak files from GSE206479 and GSE245522.  It preserves the same
PAM classes, sequence filters, and complete-site-in-peak primary call, but it
does not create a human atlas and does not validate CRISPRa activity.
