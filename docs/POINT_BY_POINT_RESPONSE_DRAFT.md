# Point-by-point response draft

Submission ID: `1a7999a3-44dc-402d-ae20-02bcc4ab2ffe`

Status: working draft for the revised submission. Line numbers must be updated
after the final manuscript export.

## Editorial summary

We thank the editor and reviewers for the detailed assessment. We have revised
the manuscript to present the work as a murine microglial CRISPRa
targetability dataset and prioritization framework rather than as a broad
atlas or a study establishing a new general principle of CRISPRa biology. The
revision narrows cross-study claims, strengthens reproducibility, adds an
exploratory human ortholog-panel ATAC check, expands methodological reporting,
and improves figure legends and quality-control documentation.

## Editor comments

### State-dependent targetability and batch effects

Response draft: We agree that the public ATAC-seq inputs do not constitute a
balanced state-comparison experiment. The revised manuscript now labels the six
columns as dataset contexts and explicitly states that cross-context
differences confound biology with laboratory, isolation, sequencing, depth, and
population-composition effects. We retained matched-depth, peak-caller, TSS,
and shared-signal analyses as robustness checks, while clarifying that these do
not remove study confounding.

Revision locations: Abstract, Results, Discussion, Methods, Figure legends,
README, `docs/ANALYSIS_CONTRACT.md`.

### Resource scope and use of "atlas"

Response draft: We agree that the term "atlas" overstated the current
deliverable. The title and text have been revised to use "dataset,"
"resource," "framework," and "targetability matrix." The manuscript no longer
presents a future web tool as part of the current resource.

Revision locations: title, Abstract, Background, Results, Discussion, Data
availability, Figure 1, README.

### Resource availability and reproducibility

Response draft: We expanded the reproducibility record with version-controlled
source tables, validation outputs, SHA-256 manifests, and explicit input
provenance. The revised Data availability section now states that a
GitHub/Zenodo archival release will accompany the revised submission. Remaining
release work is tracked in the validation log and must be completed before
resubmission.

Revision locations: Data availability, `release_manifest.sha256`,
`docs/POST_REVIEW_VALIDATION_LOG.md`.

### Positioning relative to existing chromatin-aware design work

Response draft: We revised the Background to state that chromatin-aware CRISPR
guide design is established, cited relevant design frameworks, and positioned
the contribution as a practical murine microglial promoter/PAM-by-ATAC
prioritization matrix. We also now discuss Dräger et al. as the closest
functional human iPSC-derived microglia CRISPRi/a comparator.

Revision locations: Background, Discussion, references.

### Methodological reporting, organization, figure quality, and QC

Response draft: We expanded the statistical Methods, clarified the curated
panel rationale, improved figure legends, added Cas-class multiplicity records,
and documented local validation. All manuscript values are linked to
version-controlled tables, scripts, macros, or validation outputs.

Revision locations: Methods, Figure legends, Table S7,
`docs/POST_REVIEW_VALIDATION_LOG.md`.

## Reviewer 1

### Major comment 1: biological state versus batch effects

Response draft: We agree and have revised the claim. The manuscript now treats
the six columns as surveyed dataset contexts, not as a balanced causal
comparison of biological states. We explicitly separate within-study contrasts
from across-study descriptive differences and explain the remaining batch and
sample-preparation confounding.

Status: addressed in text; final response should include exact lines.

### Major comment 2: murine-only resource and human relevance

Response draft: We agree that the primary resource is murine and preclinical.
To provide a limited translational sanity check without overextending the
resource, we added an exploratory 53-gene human ortholog-panel analysis using
processed public iPSC-/hPSC-derived microglia ATAC-seq peaks from GSE206479 and
GSE245522. The analysis is explicitly not a human atlas and not CRISPRa
validation.

Status: addressed with `scripts/human_ortholog_atac_check.py`,
`analysis_stats/human_ortholog_atac/`, and manuscript Results/Methods.

### Major comment 3: interactive web tool

Response draft: We removed claims about an under-development interactive tool
from the manuscript. Any future web tool will be described only if live and
stable.

Status: addressed.

### Minor comment 1: chromatin accessibility citations

Response draft: We added citations to established CRISPRi/a and guide-design
work linking local context/chromatin-aware design to CRISPR activity.

Status: addressed.

### Minor comment 2: coverage by one/two/multiple Cas variants

Response draft: We added a Cas-class multiplicity summary reporting genes
supported by exactly 0--5 targeting classes at the sequence layer, per-context
primary-call layer, and any-context primary-call layer. These records are
included in Table S7 and `analysis_stats/cas_multiplicity_summary.tsv`.

Status: addressed.

### Minor comment 3: potential redundancy in criteria

Response draft: We clarified that sequence-filtered candidate definition,
complete guide-site containment, promoter-midpoint coverage, and promoter/peak
overlap are distinct operational quantities. Promoter-level quantities are
retained as sensitivities and are not relabelled as guide targetability.

Status: addressed.

### Minor comment 4: 55-gene panel rationale

Response draft: We expanded the Methods to describe the curated panel
categories, representative literature, locked status before revised
targetability analysis, and intended use as a focused panel rather than a
random sample of therapeutic genes.

Status: addressed in text; final response should cite Table S1.

### Minor comments 5--9: legends, missing Cas variants, text overlap, and y-axis

Response draft: We revised Figure 3 and Figure 4 legends, made explicit that
Figure 3C is intentionally restricted to Un1Cas12f1 while all targeting classes
are reported in Figure 3B/Table S2, improved Figure 3 spacing, and clarified
Figure 6 scale information. The current tracked Figure 6 includes numeric
y-axis scales and a CPM-normalized signal label; the figure legend further
states how y-axis scales are shared within gene columns. Regenerating Figure 6
from source requires restoration or recreation of the ignored bigWig
intermediates before the final archival rebuild.

Status: addressed for the manuscript/response; final reproducibility rebuild
pending for archival packaging if derived bigWigs are included.

## Reviewer 2

### Comment 1: "atlas" terminology

Response draft: We agree and revised the title and manuscript terminology away
from "atlas." We now use dataset/resource/framework/targetability matrix and no
longer use "atlas" adjectivally for methodological choices.

Status: addressed, except immutable repository/script/environment names.

### Comment 2: reconstruction and archival repository

Response draft: We strengthened reproducibility documentation, checksums, and
release validation. We agree that a persistent archival release is needed and
will deposit the revised source/output package with GitHub/Zenodo before
resubmission. Peak files, promoter BEDs, and signal-file availability remain
explicit release-package tasks.

Status: partially addressed; archival release pending.

### Comment 3: contribution stated negatively and Dräger positioning

Response draft: We revised the contribution to be affirmative: a quantitative,
genome-wide murine microglial promoter/PAM-by-ATAC targetability matrix and
prioritization framework. We softened the gap claim and positioned the work
relative to Dräger et al.

Status: addressed in text; final response should include exact lines.

## Reviewer 3

### Novelty and conceptual claim

Response draft: We agree that chromatin accessibility as a determinant of
CRISPRi/a activity is established and have revised the manuscript so it does
not claim a new general principle. The contribution is now framed as a
practical, reproducible murine microglial resource and candidate-prioritization
framework.

Status: addressed.

### Secondary computational analysis and lack of functional validation

Response draft: We agree that the work is computational and does not validate
CRISPRa efficacy. The revised manuscript states this limitation throughout and
uses "candidate protospacers" and "hypotheses" rather than experimental
recommendations or therapeutic claims.

Status: addressed.

### Hypothesis-generation value

Response draft: We added a limited human ortholog-panel ATAC check and
reframed the TFE3/TFEB example as a concrete follow-up hypothesis. We also
state that selected-promoter TFEB support is TSS-dependent and must not be read
as a gene-wide negative claim.

Status: addressed in text, with limitations.

### Therapeutic panel rationale

Response draft: We expanded the panel-rationale Methods paragraph and retained
row-level panel metadata in Table S1. We emphasize that the panel is fixed and
curated, not a random therapeutic-gene population.

Status: addressed; final response should cite Table S1.

### Permutation/statistical details

Response draft: We added the number of resamples, seed, matching variables,
empirical p-value calculation, and BH correction family, and clarified that
fixed-panel resampling intervals are stability summaries rather than population
confidence intervals.

Status: addressed.

### Descriptive Results and caveats

Response draft: We moved stronger caveats into the Discussion, reduced
overclaiming, and emphasize the interpretive limits of each analysis layer.

Status: partially addressed; final text pass still recommended.

### Presentation, hard-coded artifacts, and QC

Response draft: We regenerated available figures, improved legends, added a
validation log, and retained a release manifest. We also disclose generative-AI
use in a short Methods statement while emphasizing author verification and
responsibility.

Status: partially addressed; full workflow/figure rebuild pending.
