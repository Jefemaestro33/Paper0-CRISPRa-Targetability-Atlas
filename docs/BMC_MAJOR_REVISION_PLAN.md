# BMC Genomics major-revision plan

Submission ID: `1a7999a3-44dc-402d-ae20-02bcc4ab2ffe`

Decision: major revision.

This plan maps the editor and reviewer requests to the concrete changes needed
for the next revision.  It is an internal execution document, not the final
point-by-point response.

## Global revision position

The safest revision strategy is to make the manuscript a rigorous, reproducible
murine microglia CRISPRa targetability dataset and prioritization framework,
not a broad "atlas" or a conceptual CRISPRa-biology paper.

Core framing to preserve:

- genome-wide murine promoter/PAM scan;
- guide-site-aware integration with microglial ATAC contexts;
- compact-nuclease relevance for vector design;
- hypothesis generation and prioritization, not functional validation.

Core framing to avoid:

- "atlas" as the primary identity of the resource unless a browsable/curated
  resource is actually delivered;
- causal biological-state claims across datasets from different labs;
- statements implying ATAC accessibility is necessary or sufficient for CRISPRa;
- claims that the human exploratory check validates therapeutic efficacy.

## Work packages

### 1. Batch/dataset confounding and biological-state language

Reviewer/editor concern: current state-dependent claims are difficult to
separate from technical/lab/batch effects because the six microglial contexts
come from different studies.

Required changes:

- Replace causal "state-dependent" claims with "dataset-context-aware" or
  "context-associated" language where cross-study comparisons are involved.
- Keep within-study contrasts separate from cross-study descriptive ranges.
- Add a limitations paragraph explaining lab, sample-preparation, sorting,
  sequencing-depth, and peak-calling effects.
- Add or highlight existing matched-depth/caller/TSS sensitivity analyses and
  explicitly state what they do and do not solve.
- In figures/tables, label contexts with provenance and avoid implying a
  balanced experimental design.

Implementation status:

- Manuscript framing now treats the six columns as dataset contexts rather than
  exchangeable biological states.
- Existing matched-depth, caller, TSS, and shared-universe signal analyses are
  retained as robustness checks that diagnose technical sensitivity but do not
  remove study confounding.
- Remaining work: figure/legend audit to ensure every panel label matches this
  narrowed interpretation.

Response stance:

- Accept the concern.
- State that the revision narrows the claim.
- Emphasize that the resource is for prioritization under public-data
  constraints, not for causal inference about microglial biology.

### 2. "Atlas" language and title

Reviewer/editor concern: "atlas" implies a curated, browsable community
resource; the current deliverable is a compressed table/resource package.

Required changes:

- Retitle away from "atlas".
- Remove adjectival phrases such as "atlas midpoint criterion", "atlas
  candidates", and "atlas language".
- Use "dataset", "resource", "framework", "targetability matrix", or
  "prioritization resource" instead.

Candidate revised title:

> A genome-wide, state-aware CRISPRa targetability dataset and prioritization
> framework for murine microglia

Alternative, more explicit title:

> Genome-wide murine microglial CRISPRa targetability prioritization by PAM
> availability and chromatin accessibility

Response stance:

- Accept.
- State that the title and manuscript terminology were revised to match the
  delivered resource.

Implementation status:

- Manuscript title has been changed to a targetability dataset/framework title.
- Main manuscript and root README have been cleaned of affirmative "atlas"
  framing, except immutable repository/script/environment names and explicit
  caveats such as "not a human atlas."

### 3. Resource availability and reproducibility

Reviewer/editor concern: key inputs are not fully deposited; GitHub alone is
not enough; release/version/license/persistent identifier are needed.

Required changes:

- Prepare a versioned GitHub release for the revised package.
- Create a Zenodo archival record/DOI for the release package.
- Include or archive:
  - reference FASTA or exact source URL plus checksum;
  - GENCODE-derived promoter BEDs;
  - primary ATAC peak files used by the workflow;
  - bigWig/source signal files or exact archived links/checksums;
  - locked panel and supplementary tables;
  - scripts/workflow/environment locks;
  - release manifest.
- Confirm license and citation metadata are present.
- If the web tool is not live and stable, remove it from claims or soften it
  to future work.

Response stance:

- Accept.
- Explain that the revised package includes a release manifest and archival DOI.
- Do not promise an interactive web tool unless it is actually live.

Implementation status:

- `release_manifest.sha256` is regenerated locally from all non-ignored release
  files.
- `docs/ARCHIVAL_RELEASE_CHECKLIST.md` records the files, inputs, checksums,
  figure intermediates, and final Data availability update needed before
  resubmission.
- `scripts/prepare_archival_package.py` creates a local GitHub/Zenodo-style
  tarball from tracked files and records whether bigWigs were included.
- Remaining work: create the actual versioned GitHub/Zenodo release after the
  final source state is approved.

### 4. Human data / translational relevance

Reviewer 1 concern: murine-only data limits clinical relevance.

Completed preparatory work:

- Added exploratory human ortholog-panel check across:
  - GSE206479 hPSC/iPSC-derived microglia ATAC, resting/IFN-beta contexts;
  - GSE245522 iPSC-derived microglia ATAC, four processed peak files.
- Integrated reproducible script:
  - `scripts/human_ortholog_atac_check.py`
- Integrated outputs:
  - `analysis_stats/human_ortholog_atac/`
- Documented caveats:
  - `docs/HUMAN_ORTHOLOG_ATAC_CHECK.md`
- Prepared archival release checklist:
  - `docs/ARCHIVAL_RELEASE_CHECKLIST.md`

Required manuscript use:

- Add as a limited exploratory translational check, probably in supplement or
  a compact Results paragraph.
- Do not call it a human atlas.
- State that both datasets are iPSC-/hPSC-derived, not primary adult human
  microglia.
- State that mouse and human contexts are not equivalent biological states.
- Use the `TFE3` / selected-`TFEB` pattern as an example of concrete
  hypothesis generation, with the `TFEB` TSS-dependence caveat.

Response stance:

- Partially accept by adding available human data.
- Clearly explain why the analysis is limited to an ortholog panel rather than
  genome-wide human reconstruction in this revision.

Implementation status:

- Compact Results and Methods subsections have been drafted in
  `manuscript/main.tex`.
- Data availability now records the GSE206479/GSE245522 processed peak inputs,
  hashes, and regeneration script.
- Human-source citations were added for Yang et al. and Booms et al.

### 5. Novelty and positioning against prior CRISPRa/chromatin-aware work

Reviewer/editor concern: chromatin accessibility effects on CRISPRa are known;
the manuscript should not imply conceptual novelty.

Required changes:

- Add explicit positioning relative to:
  - Horlbeck-style CRISPRa guide-design literature;
  - chromatin-aware guide-design frameworks/resources;
  - crisprVerse or similar guide-design ecosystems;
  - Dräger et al. / human iPSC-derived microglia CRISPRi/a platform.
- State the contribution affirmatively:
  - quantitative, genome-wide murine promoter/PAM x microglial ATAC
    targetability matrix;
  - compact-nuclease comparison;
  - state/dataset-context-aware prioritization for murine microglia;
  - reproducible resource and worked prioritization examples.
- Remove or soften negative self-framing such as "not a new principle" from
  prominent positions; keep it as a limitation if needed.

Response stance:

- Accept the novelty boundary.
- Argue practical/resource value, not conceptual discovery.

Implementation status:

- Background now cites established CRISPRi/a guide-design and chromatin-aware
  resources and positions Dräger et al. as the closest functional microglia
  CRISPRi/a platform.
- Remaining work: final response letter should explicitly say that the revised
  claim is practical/resource value, not discovery of a new CRISPRa principle.

### 6. Therapeutic 55-gene panel rationale

Reviewer concern: the curated therapeutic gene set appears arbitrary.

Required changes:

- Add a clear curation paragraph:
  - categories used;
  - inclusion logic;
  - literature basis;
  - locked before result-derived analyses;
  - role as a focused illustrative panel, not a universal therapeutic list.
- Add citations for key categories/genes where not already present.
- Ensure `config/therapeutic_genes_locked.csv` and Table S1 carry enough
  metadata.

Response stance:

- Accept.
- Reframe as a locked focused panel used for illustration and fixed-panel
  stability analysis.

Implementation status:

- Methods now states the panel categories, cites representative literature, and
  records that no gene was added or removed based on revised targetability
  results.
- Table S1 retains row-level category, priority, justification, and generated
  support pattern.

### 7. Statistical procedures and permutation testing

Reviewer concern: permutation/null procedures are insufficiently described.

Required changes:

- In Methods, explicitly describe:
  - what variable is shuffled or sampled;
  - number of permutations/bootstrap draws;
  - random seed;
  - matching constraints;
  - test statistic;
  - p-value or interval calculation;
  - BH correction family where used.
- Confirm outputs in `analysis_stats/` match the Methods text.
- Avoid overinterpreting fixed-panel resampling as population inference.

Response stance:

- Accept.
- Add reproducibility detail and cite output files.

Implementation status:

- Methods now describe fixed-panel bootstrap, matched-promoter sampling,
  empirical p-value calculation, random seed, and within-panel Fisher/BH family.
- Remaining work: verify the response letter mirrors the exact Methods language.

Reviewer 1 minor extension:

- A Cas-class multiplicity summary was added to report how many genes are
  supported by exactly 0--5 targeting classes at the sequence layer, per-context
  primary-call layer, and any-context primary-call layer.

### 8. Figures, tables, legends, and ordering

Reviewer concern: overlapping labels, small text, missing y-axis/scale, missing
legend explanations, inconsistent numbering/order.

Required changes:

- Regenerate publication-quality vector figures.
- Fix Figure 3A/3C/4C overlapping text.
- Fix Figure 4B font/readability.
- Add y-axis/scale labels to track/locus panels.
- Explain color coding in legends.
- Clarify missing Cas variants in any panel where variants are filtered or not
  displayed.
- Ensure every table/figure appears after first citation and numbering is
  consistent.

Response stance:

- Accept.
- Treat figure clarity as a resource-quality issue.

### 9. Web tool references

Reviewer/editor concern: an under-development interactive tool cannot support a
publication claim.

Required changes:

- If no live stable tool exists, remove from title/abstract/results and move to
  a short future-work sentence or remove entirely.
- If a tool is deployed later, include stable URL, archive/version, and minimal
  usage documentation.

Response stance:

- Accept.

### 10. "AI/hard-coded artifacts" / computational QC

Reviewer 3 concern: terminal-style outputs, hard-coded numerical values,
inconsistent styling, and low-quality exports suggest insufficient verification.

Required changes:

- Do not argue about AI usage unless specifically required by journal policy.
- Answer as a reproducibility/QC concern:
  - remove terminal-style artifacts from manuscript;
  - ensure figure numbers are generated from data or macros;
  - document validation scripts;
  - retain release manifest/checksums;
  - add an internal audit note if useful.
- Confirm no manuscript conclusions depend on hard-coded figure values.

Response stance:

- Accept the quality-control implication.
- State that all values in text/figures are regenerated from version-controlled
  tables/scripts.

Implementation status:

- A reproducible human-check script and input hash table were added.
- Release validation and manifest generation are used as explicit QC artifacts.
- Figure 1 and Figure 3 were regenerated and visually inspected after fixing
  terminology/spacing; Figures 1--5 can be regenerated in the local `paper0`
  environment.
- Remaining work: recreate or restore ignored `workflow/results/bigwig/*.bw`
  intermediates and rerun the complete figure build including Figure 6.

### 11. Reviewer 3 novelty/validation critique

Reviewer 3 is harsh but not fully wrong: the work is computational and lacks
functional validation.

Required changes:

- Add explicit limitation that predicted targetability is not CRISPRa efficacy.
- Remove/soften "actionable" phrasing.
- Present examples as hypotheses requiring guide-level and functional testing.
- Use the human ortholog check to demonstrate that the framework generates a
  concrete cross-dataset hypothesis, not to claim validation.

Response stance:

- Concede the validation boundary.
- Defend the resource value and practical prioritization contribution.
- Do not over-argue novelty.

## Execution order

1. Freeze current repo state and integrated human outputs.
2. Update manuscript title and terminology away from "atlas".
3. Add/adjust Results subsection for dataset-context-aware interpretation.
4. Add human ortholog-panel exploratory subsection or supplement paragraph.
5. Rewrite Discussion limitations around:
   - batch/dataset confounding;
   - murine-only primary resource;
   - human iPSC-derived exploratory check;
   - lack of CRISPRa validation;
   - transcript/TSS dependence.
6. Expand Methods:
   - peak support;
   - permutation/null testing;
   - panel curation;
   - human ortholog exploratory check;
   - reproducibility/release.
7. Regenerate and audit figures/tables.
8. Prepare Zenodo/release package.
9. Draft point-by-point response.

## Current status

- Human exploratory check: integrated and documented.
- Discretionary revision decisions: closed in `docs/REVISION_DECISIONS.md`.
- Manuscript edits: partially integrated; title, human Results/Methods,
  availability, terminology cleanup, panel-rationale expansion, and Cas-class
  multiplicity reporting are in place.
- Figure regeneration: partial; Figures 1--5 regenerated, Figure 6 requires
  ignored bigWig intermediates or a full Snakemake rebuild; current visual
  audit recorded in `docs/FIGURE_VISUAL_AUDIT.md`.
- Release/Zenodo archival package: local packaging script added; actual remote
  DOI/release pending.
- Response letter: pending.
- GitHub visibility/push: not changed locally; requires explicit confirmation.
