# Closed revision decisions

Date: 2026-08-11

This document records the implementation choices for the discretionary
major-revision items. These are execution decisions for the next revision pass,
not new scientific claims.

## 1. Human ortholog-panel analysis placement

Decision: keep the human ortholog-panel ATAC check as a brief main-text Result
and Method, with detailed documentation in `docs/HUMAN_ORTHOLOG_ATAC_CHECK.md`
and outputs under `analysis_stats/human_ortholog_atac/`.

Rationale: this directly addresses the reviewer concern about murine-only
translation while keeping the scope limited. The manuscript must continue to
call this a limited exploratory sanity check, not a human atlas and not CRISPRa
validation.

## 2. Additional human datasets

Decision: do not add another human dataset in this revision unless a clearly
superior, immediately usable microglial ATAC dataset is identified later.

Rationale: two independent iPSC-/hPSC-derived microglia peak resources already
provide a limited translational check. Expanding to a third dataset or to a
genome-wide human reconstruction would increase scope and invite new review
questions without being required to answer the current comments.

## 3. Archival release contents

Decision: prepare a GitHub/Zenodo-style archival package containing source code,
workflow files, environment files, manuscript source/PDF, figures, reference
tables, supplementary tables, analysis statistics, validation outputs, and
checksums.

Minimum expected archived material:

- tracked source/workflow/scripts/docs;
- `reference/` promoter/TSS provenance tables;
- all `supplementary/` publication tables;
- `analysis_stats/` release and revision outputs;
- `figures/output/*.pdf` and figure scripts;
- `manuscript/main.tex`, `main.pdf`, `references.bib`, and generated macros;
- `release_manifest.sha256`, `LICENSE`, and `CITATION.cff`.

## 4. BigWig/signal tracks

Decision: archive `workflow/results/bigwig/*.bw` if they are present when the
final release package is built. If they are absent, the package must explicitly
record that Figure 6 signal tracks require a full Snakemake rebuild or restored
workflow intermediates.

Current local status: the Paper0 bigWigs are not present in this checkout, and
Snakemake is not installed in the currently available local environments. A
2026-08-11 ENA header check found that the compressed FASTQs alone total
16.90 GiB, while the local filesystem had about 12 GiB available; the full
workflow would require additional reference, index, BAM, peak-calling, and
bigWig space. Figure 6 therefore remains a final-release rebuild dependency
rather than a completed local artifact regeneration step.

## 5. Figures

Decision: do not perform a full figure redesign. Use a complete visual audit
plus targeted corrections.

Already implemented:

- Figure 1: removed visible "atlas" language and fixed clipped text boxes.
- Figure 3: improved spacing and clarified color encoding/Un1Cas12f1-only panel
  scope in the legend.
- Figure 4: clarified filled-circle versus red-outline encoding in the legend.
- Figure 6: legend states that contexts share a y-axis scale within each gene
  column and that the two gene columns are scaled independently; direct PNG/PDF
  inspection confirms numeric y-axis scales and a global CPM-normalized signal
  label are present.

Remaining final-release action: regenerate Figure 6 after bigWigs are restored
or recreated, then rerun the visual audit if rebuilding the release artifacts.
