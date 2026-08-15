# Post-review validation log

Date range: 2026-08-10 to 2026-08-13

This log records the local checks run after integrating the BMC Genomics
major-revision planning material and the exploratory human ortholog-panel ATAC
analysis.

## Full configured workflow rebuild on Google Cloud VM

Date: 2026-08-13

Purpose: close the local-environment blockers noted below by rebuilding the
configured murine workflow on a machine with enough CPU, memory, and disk to
regenerate the ignored workflow intermediates needed for final reproducibility
and Figure 6 signal-track support.

Execution environment:

- Google Cloud project: `pli-biomarker`.
- VM: `paper0-cpu-runner`, `us-east1-b`, `e2-standard-8`.
- Disk: 300 GB persistent standard disk.
- Project checkout used for the rebuild: commit `56f0d17`.
- Linux lock/environment used on the VM: `environment-linux-64.lock`.
- Python: 3.11.15.
- Snakemake: 8.25.5.
- `bamCoverage`: 3.5.5.
- `Genrich`: available in the activated workflow environment.

The first full run progressed to 52 of 125 jobs and failed during
`sample_peak` for `SRR14667444` because the Genrich process was killed by the
operating system, consistent with an out-of-memory failure under the initial
parallel scheduling.

The workflow was then resumed with memory-heavy rules effectively serialized:

```bash
snakemake -s workflow/Snakefile \
  --cores 8 \
  --rerun-incomplete \
  --latency-wait 60 \
  --printshellcmds \
  --set-threads \
    sample_peak=8 \
    pooled_peak=8 \
    primary_consensus_peak=8 \
    primary_technical_peak=8 \
    condition_bam_technical=8 \
    condition_bam_biological=8 \
    condition_bigwig=8
```

Result: the resumed run completed successfully with 73 of 73 jobs finished and
0 errors. A final Snakemake dry-run reported:

```text
Nothing to be done (all requested files are present and up to date).
```

VM outputs produced and copied into the local working tree:

- 6 of 6 condition-level BigWig tracks in `workflow/results/bigwig/`.
- 6 of 6 primary peak files in `workflow/results/peaks/primary/`.
- 13 replicate peak files in `workflow/results/peaks/replicate/`.
- pooled and matched-depth peak files used for sensitivity checks.
- 13 PDF figures and 13 PNG figure renderings.
- VM Snakemake logs in `workflow/logs/`.
- regenerated supplementary tables and analysis statistics.
- regenerated `manuscript/results_macros.tex`.

The primary release validator reported:

```text
PASS
atlas_rows=647970
genes=21599
classes=5
contexts=6
therapeutic_genes=55
candidate_rows=1262
candidate_gene_class_groups=268
replicate_evidence_rows=16406
qc_run_rows=13
qc_condition_rows=6
tss_rows=64797
blacklist_mm39_rows=3360
```

Interpretation: the previous local blocker for rebuilding Figure 6 signal
intermediates and completing the full configured murine workflow is resolved.
Large derived workflow products remain git-ignored, but they are now included
by the local archival-package and release-manifest scripts when present.

## Checks that passed

### Python syntax check

Command:

```bash
python3 -m py_compile \
  scripts/human_ortholog_atac_check.py \
  scripts/rebuild_atlas_strict_iupac.py \
  scripts/analysis_statistics.py \
  scripts/summarize_results.py \
  scripts/validate_release.py
```

Result: passed with exit code 0.

### Human ortholog-panel output invariants

Validated:

- `human_panel_promoters_hg19.tsv`: 53 mapped promoters.
- `human_panel_mapping_issues.tsv`: 2 excluded mouse symbols (`Siglech`,
  `Chil3`).
- `gse206479_human_panel_targetability.tsv`: 1,060 data rows.
- `gse245522_human_panel_targetability.tsv`: 1,060 data rows.
- `gse206479_human_panel_candidate_sites.tsv`: 6,187 data rows.
- `gse245522_human_panel_candidate_sites.tsv`: 6,187 data rows.
- selected human `TFEB` promoter is not targetable in either human dataset
  under the single-promoter strict rule.
- selected human `TFE3` promoter is Un1Cas12f1-supported in all four
  GSE206479 contexts and all four GSE245522 replicate peak files.

Result: passed.

### Cas-class multiplicity summary

Command:

```bash
python3 scripts/cas_multiplicity_summary.py
```

Result: passed; `analysis_stats/cas_multiplicity_summary.tsv` was regenerated
and synchronized into Table S7 as `cas_class_multiplicity` records.

### Primary murine release validation

Command:

```bash
python3 scripts/validate_release.py \
  --atlas supplementary/table_S2_targetability_full.tsv.gz \
  --candidates supplementary/table_S3_candidate_protospacers.csv \
  --replicate-evidence supplementary/table_S3_replicate_evidence.tsv \
  --qc supplementary/table_S4_atac_qc.csv \
  --tss-selection reference/tss_selection.tsv \
  --blacklist supplementary/table_S6_blacklist_mm39_lifted.bed \
  --output analysis_stats/release_validation.txt
```

Result: passed; `analysis_stats/release_validation.txt` was regenerated.

### Core regression tests

The local `paper0` environment does not provide the `pytest` command/module, but
the current test file contains plain `test_*` functions with direct assertions.
Those tests were executed through a small import-and-call runner using the
`paper0` environment, which has the project scientific dependencies.

Result: passed; 13 tests executed.

### Release manifest

Command:

```bash
python3 scripts/build_release_manifest.py
```

Result: passed; `release_manifest.sha256` was regenerated after the integrated
source/output changes.

### Manuscript compile

Command:

```bash
cd manuscript
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

Result: passed; `manuscript/main.pdf` was regenerated. Final log check found
no LaTeX errors and no unresolved citations. Remaining warnings were underfull
boxes caused by long file paths/URLs.

### Figure regeneration and visual audit

Figure scripts were made portable through the `PYTHON_BIN` environment
variable. The default system `python3` lacks the plotting stack; the
`/Users/darellplascencia/tesis_env/bin/python` environment has
matplotlib/pandas/numpy but lacks `pyBigWig`; the conda `paper0` environment has
the plotting stack and `pyBigWig`.

Command used for the most complete local attempt:

```bash
PYTHON_BIN=/opt/homebrew/Caskroom/miniconda/base/envs/paper0/bin/python \
  ./scripts/build_figures.sh
```

Result: Figures 1--5 were regenerated. The script stopped at Figure 6 because
the ignored workflow bigWig intermediates are not present locally:

```text
workflow/results/bigwig/homeostatic.bw
```

Interpretation: Figure 6 remains reproducible from the full Snakemake workflow,
but the current lightweight local checkout cannot rebuild it unless the
`workflow/results/bigwig/*.bw` intermediates are recreated or restored. Figure 1
and Figure 3 were rendered to PNG and visually inspected after correction:

- Figure 1 no longer contains the visible "atlas" label and the previously
  clipped text boxes were fixed.
- Figure 3 spacing was improved and the gene-level color encoding is described
  in the manuscript legend.

The full current figure set was inspected as a contact sheet. The audit is
recorded in `docs/FIGURE_VISUAL_AUDIT.md`.

Follow-up check on 2026-08-11: the tracked Figure 6 PNG/PDF were inspected
directly and already include numeric y-axis scales plus the global
`CPM-normalized ATAC signal` label; the manuscript legend states that the six
contexts share a y-axis within each gene column and that the two gene columns
are scaled independently. Thus the reviewer-facing y-axis/scale issue is
addressed in the current tracked figure/legend. The remaining Figure 6 issue is
not visual content but regeneration of the ignored signal intermediates for a
final archival rebuild.

The FASTQ inputs listed in `config/samples.tsv` were checked against ENA HTTP
headers before attempting a raw-data rebuild. The compressed FASTQs alone total
18,145,869,953 bytes (16.90 GiB), while the local filesystem had only about
12 GiB available. Because the workflow would also need mm39 references, bowtie
indices, sorted BAMs, duplicate-marked/filter BAMs, condition BAMs, peak-calling
temporaries, and bigWigs, a full local FASTQ-to-bigWig rebuild is not viable on
the current disk state without freeing or attaching substantially more space.

## Checks blocked by local environment

### Pytest command-line runner

Attempted commands:

```bash
python3 -m pytest -q tests
/Users/darellplascencia/tesis_env/bin/python -m pytest -q tests
```

The standard pytest command failed before collecting tests because `pytest` is
not installed in the checked Python environments:

```text
No module named pytest
```

Interpretation: this is an environment/dependency blocker for the pytest CLI,
not a failing test assertion. The same test functions passed under the manual
runner above. The next full validation pass should still be run inside the
pinned project environment from `environment.yml` or after installing `pytest`
into the active environment.

## Remaining final release work

- Upload the final archive to the selected persistent repository and replace the
  provisional data-availability language with the final DOI/release URL.
- Optionally run the pytest command-line suite inside a pinned environment that
  includes `pytest`; the direct assertion runner already passed the current
  plain test functions.

## v2.1.0 release-packaging audit

Date: 2026-08-15

Actions completed:

- Prepared `.zenodo.json` and updated `CITATION.cff` metadata for v2.1.0.
- Added a companion reference-input packaging step for the compressed public
  reference files used by the workflow.
- Downloaded and archived:
  - `mm39.fa.gz`: 870,543,764 bytes; SHA-256
    `e558d498b49ee50b2ae14262f65c55b0feae733e833d1f4f0da0c0b789fffe59`.
  - `gencode.vM33.annotation.gtf.gz`: 29,297,942 bytes; SHA-256
    `6b0218f0c587591053a099797f049db5f6fcfa1a3cff6dae912db3e045c0fd90`.
- Wrote the tracked reference-source audit table:
  `reference/source_input_audit.tsv`.
- Built the local companion reference archive:
  `release_archives/paper0_reference_inputs_20260815T060647Z.tar.gz`.
- Regenerated `release_manifest.sha256`; it now contains 210 hashes.
- Re-ran the primary release validator; result: `PASS`.
- Confirmed manuscript PDF is up to date with `latexmk`.
- Confirmed Python syntax for release, manifest, validation, human-check,
  targetability, statistics, and summarization scripts.
- The local Python environments still lack the `pytest` module, but the direct
  assertion runner executed the 13 plain `test_*` functions successfully in the
  project scientific environment.

Interpretation: the final release package now covers both the derived peak and
signal artifacts requested by reviewers and the compressed public FASTA/GTF
reference inputs needed to reconstruct the workflow.
