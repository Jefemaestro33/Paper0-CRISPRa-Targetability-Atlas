# Revision 2.0 audit ledger

This ledger maps the methodological audit to concrete, version-controlled
changes.  It records implementation status; numerical conclusions are not
considered final until the raw-to-figure workflow completes and the validation
checks pass.

| Audit issue | Resolution in revision 2.0 | Primary evidence |
|---|---|---|
| Midpoint-based targetability | Replaced by complete protospacer+PAM containment in a primary peak; midpoint and overlap rules are sensitivities only. | `scripts/rebuild_atlas_strict_iupac.py`, Table S2 |
| State versus dataset confounding | Six columns are labelled dataset contexts; only within-study contrasts are identified and cross-study ranges remain descriptive. | manuscript, README |
| Peak-set provenance | Biological datasets use pairwise 50% reciprocal-overlap intersections supported by at least two independent replicates; technical-only preparations are pooled and explicitly labelled. | `workflow/scripts/consensus_peaks.py`, Table S4 |
| Depth and caller dependence | Biological runs are deterministically downsampled at replicate level to the minimum condition-total count and rebuilt with the same consensus rule; technical pools remain explicitly separate. Both Genrich and MACS3 are evaluated, and requested/realized depths, seeds, caller concordance, and shared-universe signal are retained in the release. | `workflow/scripts/matched_depth_analysis.py`, `analysis_stats/` |
| TSS choice | Ensembl canonical is primary; APPRIS principal and the prior most-5-prime GENCODE-basic rule are rebuilt as sensitivities. | `scripts/prepare_reference.py`, `reference/` |
| Six ortholog claim | Reframed as five nuclease/PAM targeting classes represented by six CRISPRa configurations; HEAL and SminiCRa share one class. | Table 1, Table S2 |
| Absolute chromatin language | Accessibility is described as a prioritization feature, not a necessary or sufficient condition for CRISPRa. | manuscript |
| Statistical interpretation | Bootstrap output is explicitly fixed-panel resampling stability; a promoter-matched non-panel null and within-panel Fisher tests are separate. | `scripts/analysis_statistics.py`, `scripts/within_panel_functional_associations.py` |
| Functional enrichment universe | Genome-wide GO claims were removed; functional associations use the frozen 55-gene universe with BH correction across the complete family. | Table S7, Supplementary Figure S5 |
| Result-derived panel metadata | Legacy accessibility-pattern labels were removed from the frozen panel input; revised support patterns are generated only from the rebuilt atlas. | `config/therapeutic_genes_locked.csv`, Table S5 |
| Candidate recommendations | Renamed candidate protospacers and supplied exact target, spacer, PAM, TSS, peak, summit, signal, context, and individual-run provenance. | Table S3 and run-level companion |
| Preliminary off-target screen | Un1Cas12f1 candidates receive exhaustive Bowtie alignments through three mismatches followed by orientation-aware TTTR validation and locus annotation; limitations are explicit. | Table S3 alignments |
| Genome-wide presentation | Genome-wide sequence and guide-site coverage are moved into the main results and main Figure 2. | manuscript, Figure 2 |
| Reproducibility | Added a raw-FASTQ Snakemake DAG, checksums, pinned environment, tests, fail-closed release invariants, analysis contract, output schema, license, citation metadata, deterministic gzip for the full atlas, and release files. | repository root, `workflow/`, `tests/`, `analysis_stats/release_validation.txt` |
| Figure legibility | Figures are regenerated as individual vector PDFs; locus panels use canonical-TSS-centred zooms with peaks and candidate positions. | `figures/output/` |
| Reference metadata | Every DOI is checked against Crossref and the identified HEAL, SminiCRa, crisprVerse, Dräger, and stroke metadata errors are corrected. | `scripts/audit_reference_dois.sh`, `manuscript/references.bib` |

Additional implementation defects were found during execution and fixed
before the final rebuild:

1. paired-end blacklist removal now drops an intact fragment when either mate
   overlaps a blacklist interval, avoiding flagged-but-orphaned reads in QC and
   downsampling;
2. the FASTA index is a declared workflow product, preventing concurrent atlas
   jobs from racing while creating the same `.fai` file.
3. the inherited mm39 blacklist contained 18 inverted intervals; it was
   replaced by a fresh, validated liftOver from the byte-verified official
   ENCODE mm10 source, with every mapped and unmapped interval recorded.
4. technical sequencing lanes from the same biological preparation are now
   assigned a shared library identity and deduplicated again after pooling,
   while independent biological replicates are never deduplicated together.
5. the MACS3 matched-depth sensitivity uses native fragment intervals in
   paired-end mode and applies shift/extension parameters only to single-end
   inputs;
6. wildcard-dependent temporary BAM paths are expanded per run or condition,
   preventing concurrent peak jobs from writing to one literal temporary path;
   a regression test now fails if this isolation is lost; and
7. candidate-table terminology now uses `candidate_class` throughout rather
   than implying that computationally ranked protospacers are recommendations.
