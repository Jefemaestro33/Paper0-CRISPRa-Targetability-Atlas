# Post-review validation log

Date: 2026-08-10

This log records the local checks run after integrating the BMC Genomics
major-revision planning material and the exploratory human ortholog-panel ATAC
analysis.

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

## Remaining validation work

- Recreate or restore `workflow/results/bigwig/*.bw` and rerun the complete
  figure build, including Figure 6, before final archival release.
- Run the full raw-to-result Snakemake workflow before final resubmission or
  archival release.
