# Figure visual audit

Date: 2026-08-15

Audit mode: visual inspection of all current main and supplementary figure PNG
previews generated from the rebuilt PDF figures after the full VM workflow
rebuild. This audit is intended to catch obvious clipping, missing labels, and
reviewer-flagged readability issues before the final public release.

Contact sheet used locally:

```text
/private/tmp/paper0_all_figures_contact_sheet.png
```

## Main figures

- Figure 1: passed after correction. Text boxes are no longer clipped and the
  visible output label now says "Genome-wide targetability matrix" rather than
  "atlas."
- Figure 2: passed for current revision. Dense but readable; no obvious clipped
  labels. Multi-Cas coverage is now also reported numerically in Table S7.
- Figure 3: passed after correction. Panel C remains dense because it shows all
  55 locked genes, but labels are readable in the vector/PDF version and the
  legend now explains teal/beige cells and the Un1Cas12f1-only scope.
- Figure 4: passed for current revision. Panel D contains small row labels but
  functions as a compact stability matrix; the legend now explains panel B
  encoding.
- Figure 5: passed for current revision. Candidate labels are readable in the
  vector/PDF version.
- Figure 6: passed after the VM rebuild. The restored
  `workflow/results/bigwig/*.bw` intermediates regenerate the track layer; the
  tracked PDF includes y-axis scales and the CPM-normalized signal label.

## Supplementary figures

- Figure S1: passed for vector/PDF inspection; small x-axis labels are expected
  because this is a run-level QC summary.
- Figure S2: passed.
- Figure S3: passed.
- Figure S4: passed.
- Figure S5: passed; dense labels are appropriate for a supplementary
  association summary.
- Figure S6: passed.
- Figure S7: passed.

## Final release status

The reviewer-flagged figure issues are resolved in the current tracked figure
set. If any figure script is edited again before resubmission, rerun the
contact-sheet check and inspect Figure 6 specifically for y-axis scale, signal
track readability, and candidate/peak overlay placement.
