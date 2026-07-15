# Paper 0: CRISPRa Targetability Atlas in Murine Microglia

This repository contains the manuscript, analysis scripts, figures, and
supplementary tables for a genome-wide CRISPRa targetability atlas across murine
microglial states.

## Manuscript

**Title:** A genome-wide CRISPRa targetability atlas across murine microglial
states integrating PAM availability with chromatin accessibility

**Author:** Ernest Darell Zermeno

The study is a predictive in silico genomics resource. It integrates promoter
PAM availability for six Cas orthologs with ATAC-seq-derived chromatin
accessibility across six murine microglial states. The atlas is intended to
prioritize experimentally testable CRISPRa hypotheses, not to establish
validated CRISPRa activity or therapeutic efficacy.

## Data Sources

The analysis uses publicly available ATAC-seq datasets:

| Dataset | GEO | States used |
|---|---|---|
| Gosselin et al. | GSE89960 | homeostatic adult microglia |
| Zhang X / Kracht et al. | GSE175578 | naive, acute LPS, LPS-tolerized |
| Zhang L et al. | GSE220041 | post-surgical sham, post-ischemic stroke |

Genome assembly: GRCm39/mm39

Gene annotation: GENCODE vM33

## Repository Structure

```text
Paper0-CRISPRa-Targetability-Atlas/
├── manuscript/
│   ├── main.tex
│   ├── main.pdf
│   └── references.bib
├── figures/
│   ├── scripts/
│   └── output/
├── supplementary/
│   ├── table_S1_therapeutic_genes.csv
│   ├── table_S2_targetability_full.tsv
│   ├── table_S3_sgrna_recommendations.csv
│   ├── table_S4_atac_qc.csv
│   ├── table_S5_accessibility_dynamics.csv
│   ├── table_S6_blacklist_mm10_original.bed
│   ├── table_S6_blacklist_mm39_lifted.bed
│   └── table_S7_statistical_tests.csv
├── analysis_results/
│   ├── browser_tracks_Tfeb.pdf/png
│   ├── browser_tracks_Tfe3.pdf/png
│   ├── frip_results.tsv
│   ├── tss_enrichment.tsv
│   ├── fragsize_*.tsv/png
│   ├── idr_*.txt
│   └── bigwigs/
├── analysis_stats/
│   ├── bootstrap_cis.tsv
│   ├── pam_chromatin_loss.tsv
│   └── permutation_results.tsv
└── scripts/
    ├── README.md
    ├── rebuild_atlas_strict_iupac.py
    ├── analysis_statistics.py
    └── sync_analysis_outputs.py
```

## Supplementary Files

The manuscript lists the following supplementary files:

| File | Description |
|---|---|
| `table_S1_therapeutic_genes.csv` | Curated 55-gene therapeutic panel |
| `table_S2_targetability_full.tsv` | Full gene x Cas x state atlas matrix used by figure scripts and audit checks |
| `table_S3_sgrna_recommendations.csv` | Candidate sgRNAs with atlas state flags and non-validated heuristic scores |
| `table_S4_atac_qc.csv` | ATAC-seq QC metrics |
| `table_S5_accessibility_dynamics.csv` | Six-state accessibility classification for therapeutic genes |
| `table_S6_blacklist_*.bed` | ENCODE blacklist files used in filtering |
| `table_S7_statistical_tests.csv` | Bootstrap intervals, descriptive PAM-to-chromatin loss counts, and peak-shuffle random-placement results |

`table_S2_targetability_full.tsv` is approximately 39 MB and is tracked so that
the main figures and numerical audit checks can run from a fresh clone.

## Main Analyses

The atlas evaluates:

- promoter PAM availability in the CRISPRa-optimal window (-400 to -50 bp from
  TSS);
- six Cas orthologs relevant to CRISPRa design;
- ATAC-seq promoter accessibility across six microglial states;
- predicted PAM + chromatin targetability;
- accessibility dynamics among 55 curated therapeutic genes;
- robustness/transparency analyses including bootstrap confidence intervals,
  descriptive PAM-to-chromatin loss counts, peak-shuffle random-placement
  analysis, and accessibility-criterion sensitivity analysis.

## Reproducibility Notes

The matrix-level canonical rebuild is `scripts/rebuild_atlas_strict_iupac.py`.
It scans both strands with IUPAC-aware reverse complements, so Cas12f PAMs are
reported as functional `TTTA`/`TTTG` PAMs in Table S3. The raw genome FASTA,
promoter BED, and ATAC peak files are supplied as command-line inputs because
the mouse genome and raw/intermediate sequencing files are not stored in this
repository. The GEO accessions and reference resources are listed above and in
the manuscript.

To compile the manuscript locally:

```bash
cd manuscript
pdflatex -interaction=nonstopmode main.tex
bibtex main
pdflatex -interaction=nonstopmode main.tex
pdflatex -interaction=nonstopmode main.tex
```

To regenerate the atlas matrix and sgRNA table from local reference inputs:

```bash
python3 scripts/rebuild_atlas_strict_iupac.py \
  --fasta path/to/mm39.fa \
  --promoters path/to/promoters_crispra_optimal.bed
```

To regenerate statistical tables from Table S2 and sync the smaller
supplementary tables:

```bash
python3 scripts/analysis_statistics.py
python3 scripts/sync_analysis_outputs.py
```

## Important Limitations

This repository supports a predictive computational resource. The atlas does not
experimentally validate CRISPRa activation, sgRNA efficacy, or therapeutic
benefit. Chromatin accessibility should be interpreted as a necessary but not
sufficient condition for CRISPRa activity.

The peak-shuffle analysis uses genome-wide shuffled peaks and is not a
promoter-matched null model. It is best interpreted as a random-placement sanity
check showing that observed promoter overlaps are not compatible with uniformly
shuffled peaks. It does not test whether therapeutic promoters are enriched
relative to matched promoters, nor whether each promoter-level prediction is
experimentally functional.

## License

No license has been assigned yet. Contact the author before reuse beyond
citation.
