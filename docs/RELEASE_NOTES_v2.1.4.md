# Release v2.1.4 notes

Date: 2026-08-17

This release closes the final pre-submission citation and internal-consistency
audit after the v2.1.3 archival package.

## Main changes since v2.1.3

- Corrected remaining Crossref-audited author metadata in the TFEB reference
  (`settembre2011tfeb`) and added the SAMtools group author.
- Added missing issue-number metadata for DOI-backed references where Crossref
  reports an issue number.
- Reordered main-text figure references so Figures 1--6 are first cited in
  sequence.
- Added an explicit manual-selection procedure for the frozen curated gene
  panel.
- Added reusable reference-metadata and manuscript-integrity audit scripts.

## Validation

- `scripts/audit_reference_metadata.py manuscript/references.bib`: PASS.
- `scripts/audit_manuscript_integrity.py --repo .`: PASS.
- `scripts/validate_release.py`: PASS.
- Direct execution of the core regression tests: 13/13 PASS.
- `shasum -a 256 -c release_manifest.sha256`: PASS.

## Persistent identifiers

- GitHub tag: `v2.1.4`
- Zenodo concept DOI: `10.5281/zenodo.21970940`
- The version-specific Zenodo record is associated with this release after
  archival publication.
