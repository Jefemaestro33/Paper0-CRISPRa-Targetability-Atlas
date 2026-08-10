# Human ortholog ATAC support check

Exploratory human ortholog-panel analysis using hg19 iPSC-derived microglia ATAC-seq.
Primary calls require the complete protospacer+PAM interval inside a dataset peak.

- mapped promoters analyzed: 53
- mapping exclusions: 2 (Siglech, Chil3)

## Any-context/replicate targetability

| Cas/PAM | GSE206479 any context | GSE245522 any replicate | GSE245522 >=3/4 replicate-supported |
|---|---:|---:|---:|
| Un1Cas12f1_TTTR | 32/53 (60.4%) | 20/53 (37.7%) | 15/53 (28.3%) |
| SaCas9_NNGRRT | 38/53 (71.7%) | 29/53 (54.7%) | 21/53 (39.6%) |
| SpCas9_NGG | 42/53 (79.2%) | 36/53 (67.9%) | 29/53 (54.7%) |
| CjCas9_NNNVRYM | 42/53 (79.2%) | 35/53 (66.0%) | 26/53 (49.1%) |
| Nme2Cas9_NNNNCC | 42/53 (79.2%) | 36/53 (67.9%) | 26/53 (49.1%) |

## TFE3/TFEB interpretation guardrail

- Selected TFE3 promoter: supported in every GSE206479 context and every GSE245522 replicate for all five PAM classes.
- Selected TFEB promoter: unsupported in both human datasets for all five PAM classes.
- TFEB sensitivity: alternative GENCODE v19 TFEB TSS choices can be accessible; therefore the negative TFEB result is promoter/TSS-dependent.

## GSE206479 vs GSE245522 >=3/4 comparison

- Un1/TTTR both_gse206479_any_and_gse245522_3of4: 15
- Un1/TTTR gse206479_only_vs_3of4: 17
- Un1/TTTR gse245522_3of4_only: 0
- Un1/TTTR neither: 21

- TFEB-201 audit rows: 5
- supported alternative TFEB transcript/Cas rows: 50

Interpretation: use as a cross-dataset human ortholog sanity check, not as CRISPRa validation or a human atlas.
