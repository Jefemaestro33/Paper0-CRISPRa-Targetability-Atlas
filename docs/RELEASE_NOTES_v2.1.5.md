# Release v2.1.5 notes

Date: 2026-08-17

This release closes the final archival-consistency checks after v2.1.4.

## Main changes since v2.1.4

- Restored ignored archival artifacts locally before repackaging so
  `release_manifest.sha256` verifies completely in the release workspace.
- Updated `prepare_archival_package.py` to exclude files already tracked by Git
  from the optional ignored-artifact pass, preventing duplicate peak files in the
  archival tarball.
- Clarified the manuscript text so the candidate resource is cited as Tables
  S3a--S3c, matching the Additional files inventory.
- Rebuilt manuscript outputs, figures, validation files, and release archives.

## Validation

- `scripts/audit_reference_metadata.py manuscript/references.bib`: PASS.
- `scripts/audit_manuscript_integrity.py --repo .`: PASS.
- `scripts/validate_release.py`: PASS.
- Direct execution of the core regression tests: 13/13 PASS.
- `shasum -a 256 -c release_manifest.sha256`: PASS.
- `tar -tzf <v2.1.5 archival tarball> | sort | uniq -d`: no duplicate archive
  members.

## Persistent identifiers

- GitHub tag: `v2.1.5`
- Zenodo concept DOI: `10.5281/zenodo.21970940`
- The version-specific Zenodo record is associated with this release after
  archival publication.
