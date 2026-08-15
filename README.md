# CRISPRa guide-site targetability dataset in murine microglia

This repository contains the manuscript and the end-to-end computational
workflow for a genome-wide, microglia-focused CRISPRa guide-site targetability
resource.
The resource scans a fixed promoter window, identifies candidate
protospacer--PAM intervals, and asks whether each **complete interval** is
contained in an ATAC-seq peak. It is a prioritization resource: a positive
call does not demonstrate CRISPRa efficacy, safety, or therapeutic benefit.

## What the revision changes

The primary analysis no longer uses coverage of an arbitrary promoter midpoint.
It now requires the full protospacer plus PAM to be contained in the relevant
peak. Promoter-midpoint coverage, any promoter/peak overlap, and any partial
guide/peak overlap are retained only as sensitivity analyses.

Other prespecified safeguards include:

- Ensembl-canonical TSSs as the primary definition, with APPRIS-principal and
  the legacy most-5-prime GENCODE-basic choice as sensitivities;
- biological-replicate consensus peaks for Gosselin, sham, and stroke data;
- explicit treatment of the Zhang X/Kracht runs as technical runs, without a
  biological reproducibility claim, including pooled-library deduplication
  across technical lanes;
- deterministic matched-depth reanalysis with Genrich and MACS3;
- a shared peak universe with continuous signal quantification;
- a promoter-matched null and fixed-panel bootstrap stability intervals;
- PAM-aware, exhaustive Bowtie screening through three mismatches for the
  Un1Cas12f1 candidate subset, labelled as a preliminary computational screen;
- exact source accession, checksum, software, random-seed, and TSS provenance.
- byte-verified ENCODE mm10 blacklist provenance and a validated UCSC
  mm10-to-mm39 liftOver (3,360 mapped and 75 unmapped intervals).

The targeting comparison comprises **five distinct nuclease/PAM classes**
across six modelled configurations. HEAL and SminiCRa are two Un1Cas12f1
activation architectures represented by the same TTTR targeting class, so they
are not counted as independent targeting rules. The dNme2Cas9 activator is a
proposed sequence-level configuration, not a validated CRISPRa platform.

## Public input data

| Study | Accession | Contexts used | Replication used here |
|---|---|---|---|
| Gosselin et al. | GSE89960 | ex vivo adult homeostatic microglia | 2 biological replicates |
| Zhang X/Kracht et al. | GSE175578 | PBS/PBS control, PBS/LPS, LPS/LPS | paired technical runs only |
| Zhang et al. | GSE220041 | sham and post-ischemic stroke | 2 sham and 3 stroke biological replicates |

The exact 13 ENA run URLs and MD5 checksums are frozen in
[`config/samples.tsv`](config/samples.tsv). References are GRCm39/mm39 and
GENCODE mouse vM33.

The six context columns are best understood as **surveyed dataset contexts**,
not as a balanced causal experiment. Across-study differences remain confounded
by laboratory, isolation, sequencing layout, depth, and cell composition. The
within-study contrasts are PBS/PBS--PBS/LPS--LPS/LPS and sham--stroke.

## Reproduce the analysis

Create and activate the pinned environment, then run Snakemake:

```bash
micromamba create -f environment.yml
micromamba activate paper0-atlas
snakemake -s workflow/Snakefile --cores 8 --rerun-incomplete
```

The workflow downloads and checksum-verifies the public FASTQs and reference
files, performs uniform ATAC-seq processing, rebuilds all targetability tables and
sensitivities, and regenerates the figures. Raw data, reference indexes, BAMs,
bigWigs, and other intermediates are intentionally ignored by Git.

Run the core definition regression tests with:

```bash
pytest -q tests
```

## Repository map

```text
.
├── config/
│   ├── samples.tsv                         # ENA provenance and checksums
│   └── therapeutic_genes_locked.csv       # frozen 55-gene panel metadata
├── workflow/
│   ├── Snakefile                          # complete raw-to-figure DAG
│   └── scripts/                           # ATAC QC, consensus, depth/caller tests
├── scripts/
│   ├── prepare_reference.py               # canonical and sensitivity TSSs
│   ├── prepare_blacklist.py               # validated mm10-to-mm39 blacklist lift
│   ├── rebuild_atlas_strict_iupac.py      # guide-site-aware matrix rebuild
│   ├── analysis_statistics.py             # stability and matched-null analyses
│   ├── sync_analysis_outputs.py           # compact supplementary tables
│   ├── summarize_results.py               # manuscript macros and result audit
│   └── validate_release.py                 # fail-closed release invariants
├── reference/                             # tracked TSS provenance tables/BEDs
├── docs/                                  # frozen analysis contract and schemas
├── supplementary/                         # publication source tables
├── analysis_stats/                        # statistical and sensitivity outputs
├── figures/
│   ├── scripts/
│   └── output/
├── manuscript/
│   ├── main.tex
│   ├── main.pdf
│   ├── references.bib
│   └── results_macros.tex                 # workflow-generated values
├── tests/
├── environment.yml
├── release_manifest.sha256                # SHA-256 for tracked and selected archived artifacts
└── LICENSE
```

## Primary operational definition

For a gene, nuclease/PAM class, and dataset context, the primary call is true
when at least one sequence-filtered candidate has its complete protospacer and
PAM interval inside the applicable primary ATAC-seq peak set within the
transcription-oriented window $-400$ to $-50$ bp from the selected TSS.

For biologically replicated datasets, primary peaks contain sequence supported
by at least two independent replicate peak calls after 50% reciprocal overlap.
For the Zhang X/Kracht contexts, technical runs are pooled and the resulting
calls are explicitly labelled as lacking biological replication.

## Main outputs

| File | Contents |
|---|---|
| `table_S1_therapeutic_genes.csv` | Frozen therapeutic panel and provenance |
| `table_S2_targetability_full.tsv.gz` | Gzip-compressed genome-wide gene × class × context targetability matrix (directly readable by pandas and standard command-line tools) |
| `table_S3_candidate_protospacers.csv` | Candidate sites, peak evidence, TSS coordinates, and preliminary off-target summary |
| `table_S3_offtarget_alignments.tsv` | PAM-valid Un1Cas12f1 alignments through three mismatches |
| `table_S3_replicate_evidence.tsv` | Candidate-by-run peak support, signal, peak ID, and summit distance |
| `table_S4_atac_qc.csv` | Per-run QC and primary-peak provenance |
| `table_S5_accessibility_dynamics.csv` | Descriptive cross-context panel patterns |
| `table_S7_statistical_tests.csv` | Fixed-panel stability, Cas-class multiplicity, and sensitivity summaries |
| `analysis_stats/cas_multiplicity_summary.tsv` | Genes supported by exactly 0--5 targeting classes at sequence and primary-call layers |
| `analysis_stats/matched_depth_summary.tsv` | Requested/realized matched depth, seeds, and replicate-handling provenance |
| `analysis_stats/peak_caller_concordance.tsv` | Matched-depth Genrich/MACS3 promoter-call concordance |
| `analysis_stats/shared_peak_signal.tsv` | Continuous counts and CPM on the shared primary-peak universe |
| `analysis_stats/release_validation.txt` | Machine-checked release invariants and expected dimensions |
| `release_manifest.sha256` | Cryptographic hashes for tracked files plus selected ignored archival artifacts present locally |
| `reference/source_input_audit.tsv` | Source URLs, byte counts, and SHA-256 hashes for compressed reference inputs packaged with the release |

## Interpretation limits

ATAC accessibility is a useful guide-site prioritization variable, not a binary
requirement for successful activation. The resource does not model local
nucleosome dynamics, 3D contacts, activation-domain potency, guide RNA folding,
cell delivery, paralog-specific biology, or in vivo benefit. Cross-study
contrasts cannot be assigned to biological state alone, and the technical-run
contexts do not provide biological replication. All candidate protospacers
require nuclease-specific computational review and experimental validation.

## Author and license

Ernest Darell Zermeño (ORCID 0009-0002-1721-4526). Released under the
[MIT License](LICENSE).
