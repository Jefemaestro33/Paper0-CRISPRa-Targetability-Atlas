# Human ortholog ATAC support check

This document records the post-review exploratory human analysis added after
the BMC Genomics major-revision decision.  It is intended to support the
response to reviewers, not to convert the murine resource into a human atlas.

## Scope

- Question addressed: do selected murine prioritization patterns survive in
  independent human microglia-like ATAC-seq datasets?
- Human universe: ortholog-mapped version of the locked 55-gene focused panel.
- Mapped genes: 53.
- Excluded genes: `Siglech` and `Chil3`, because a direct one-to-one human
  ortholog was not used for this exploratory check.
- Reference: GENCODE v19 / GRCh37-hg19.
- Promoter window: transcription-oriented −400 to −50 bp from the selected TSS.
- Primary call: complete protospacer+PAM interval fully contained in a
  context- or replicate-specific ATAC peak.
- PAM classes and sequence filters: identical to the murine guide-site-aware
  rebuild (`scripts/rebuild_atlas_strict_iupac.py`).

The output is archived under `analysis_stats/human_ortholog_atac/` and can be
regenerated with:

```bash
python scripts/human_ortholog_atac_check.py
```

Raw public inputs are cached under `workflow/resources/human_ortholog_atac/`
when the script is run.  That cache is intentionally ignored by Git.

## Human datasets

### GSE206479

Used because it provides processed hg19 ATAC peak calls for human
hPSC-derived microglia in resting and IFN-beta-stimulated contexts across two
cell lines (`WTC11`, `H1`).  The processed files are IDR optimal narrowPeak
sets from the GEO supplementary files.

Contexts used:

- `WTC11_rest`
- `WTC11_IFNb`
- `H1_rest`
- `H1_IFNb`

### GSE245522

Used as an independent human iPSC-derived microglia ATAC-seq dataset with
processed hg19 MACS2 narrowPeak files.  The GEO TAR contains four peak files
(`10K`, `31K`, `100K`, `100K2`) and no consensus peak file, so the analysis
treats them as replicate/input-size peak contexts and reports both
any-replicate support and support by at least 3 of 4 replicates.

Contexts used:

- `10K_peak_file`
- `31K_peak_file`
- `100K_peak_file`
- `100K2_peak_file`

## Integrated output files

- `analysis_stats/human_ortholog_atac/README.md`: compact result summary.
- `analysis_stats/human_ortholog_atac/input_file_audit.tsv`: input URLs,
  row counts, and SHA-256 hashes.
- `analysis_stats/human_ortholog_atac/human_panel_promoters_hg19.tsv`: selected
  human promoters and transcript-selection metadata.
- `analysis_stats/human_ortholog_atac/human_panel_mapping_issues.tsv`: excluded
  panel genes.
- `analysis_stats/human_ortholog_atac/human_panel_sequence_qc.tsv`: promoter
  sequence-length and base-content checks.
- `analysis_stats/human_ortholog_atac/gse206479_human_panel_targetability.tsv`
  and `gse245522_human_panel_targetability.tsv`: gene x Cas x context calls.
- `analysis_stats/human_ortholog_atac/gse206479_human_panel_candidate_sites.tsv`
  and `gse245522_human_panel_candidate_sites.tsv`: candidate-site-level calls.
- `analysis_stats/human_ortholog_atac/gse245522_consensus_like_summary.tsv`:
  1/4, 2/4, 3/4, and 4/4 replicate-support summaries.
- `analysis_stats/human_ortholog_atac/gse206479_vs_gse245522_3of4_comparison.tsv`:
  cross-human-dataset comparison using GSE206479 any-context support and
  GSE245522 >=3/4 peak-file support.
- `analysis_stats/human_ortholog_atac/mouse_vs_gse206479_un1_comparison.tsv`:
  Un1/TTTR comparison between the murine table S2 and the first human dataset.
- `analysis_stats/human_ortholog_atac/tfeb_alternative_tss_audit.tsv`:
  transcript/TSS sensitivity audit for human `TFEB`.

## Key numerical result

Any-context/replicate strict targetability across the 53 mapped genes:

| Cas/PAM | GSE206479 any context | GSE245522 any peak file | GSE245522 >=3/4 peak-file-supported |
|---|---:|---:|---:|
| Un1Cas12f1_TTTR | 32/53 (60.4%) | 20/53 (37.7%) | 15/53 (28.3%) |
| SaCas9_NNGRRT | 38/53 (71.7%) | 29/53 (54.7%) | 21/53 (39.6%) |
| SpCas9_NGG | 42/53 (79.2%) | 36/53 (67.9%) | 29/53 (54.7%) |
| CjCas9_NNNVRYM | 42/53 (79.2%) | 35/53 (66.0%) | 26/53 (49.1%) |
| Nme2Cas9_NNNNCC | 42/53 (79.2%) | 36/53 (67.9%) | 26/53 (49.1%) |

The exact percentages are dataset-dependent.  The robust qualitative pattern
is that PAM availability remains broader than accessibility-supported
targetability.

## TFEB/TFE3 interpretation

The selected human `TFE3` promoter is supported in every GSE206479 context and
every GSE245522 peak file for all five PAM classes.  The selected human
`TFEB` promoter is unsupported in both human datasets for all five PAM classes.

This is useful, but it must be stated precisely:

- the result supports a cross-dataset, cross-species prioritization hypothesis;
- it does not demonstrate CRISPRa efficacy;
- it does not show that all human `TFEB` promoters are inaccessible;
- it applies to the promoter/TSS selected by the same deterministic rule used
  for this exploratory analysis.

The `TFEB` alternative-TSS audit found that several GENCODE v19 `TFEB`
transcripts have accessible promoter windows.  Therefore the correct language
is not "human TFEB is not targetable"; the correct language is:

> Using the same single-promoter selection rule as the murine resource, the
> selected human TFEB promoter lacked ATAC-supported guide sites in two
> independent human iPSC-derived microglia datasets, whereas TFE3 was
> consistently supported.  A transcript-level sensitivity analysis showed that
> TFEB support is promoter/TSS-dependent.

## Internal audit status

The integrated script was run against the cached public inputs and regenerated
the tracked outputs.  Independent checks confirmed:

- 53 unique mapped promoters;
- all promoter windows are 350 bp;
- no ambiguous bases in fetched promoter sequences;
- selected `TFEB` and `TFE3` sequences match UCSC hg19;
- expected targetability table dimensions: 53 genes x 5 PAM classes x 4
  contexts = 1,060 data rows per dataset;
- expected candidate-site rows: 6,187 data rows per dataset;
- an independent reimplementation of the scan/intersection logic produced 0
  mismatches across 1,060 GSE206479 rows and 0 mismatches across 1,060
  GSE245522 rows.

## Recommended use in the revision

Use this as a small, controlled addition:

1. Add one subsection or supplement-labelled analysis, not a new human atlas.
2. Frame it as a human ortholog-panel robustness/sanity check.
3. Explicitly state that both human datasets are iPSC-/hPSC-derived microglia
   and do not substitute for primary adult human microglia or functional
   CRISPRa validation.
4. Use the `TFE3`/selected-`TFEB` contrast as a worked example of hypothesis
   generation, with the TSS-dependence caveat included in the same paragraph.
