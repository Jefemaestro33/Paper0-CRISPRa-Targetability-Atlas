"""Small regression tests for the atlas definitions most likely to drift."""
from __future__ import annotations

import importlib.util
import gzip
import hashlib
import sys
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module


atlas = load_module("atlas", ROOT / "scripts/rebuild_atlas_strict_iupac.py")
consensus = load_module("consensus", ROOT / "workflow/scripts/consensus_peaks.py")
reference = load_module("reference", ROOT / "scripts/prepare_reference.py")
matched_depth = load_module("matched_depth", ROOT / "workflow/scripts/matched_depth_analysis.py")
replicate_evidence = load_module(
    "replicate_evidence", ROOT / "workflow/scripts/candidate_replicate_evidence.py"
)


def test_promoter_interval_respects_transcriptional_orientation():
    plus = reference.Transcript("chr1", 999, 2000, "+", "g", "G", "t1", "", (), "1")
    minus = reference.Transcript("chr1", 999, 2000, "-", "g", "G", "t2", "", (), "1")
    assert reference.promoter_interval(plus, 400, 50) == (599, 949)
    assert reference.promoter_interval(minus, 400, 50) == (2050, 2400)


def test_historical_midpoint_is_minus_225_on_both_strands():
    plus = atlas.Promoter("chr1", 600, 950, "G", "+", 1000, "g", "t", "d", "s")
    minus = atlas.Promoter("chr1", 1051, 1401, "G", "-", 1000, "g", "t", "d", "s")
    assert plus.midpoint == 775
    assert minus.midpoint == 1225
    assert atlas.oriented_distance(plus, plus.midpoint) == -225
    assert atlas.oriented_distance(minus, minus.midpoint) == -225


def test_peak_match_requires_complete_guide_for_primary_call():
    peak = atlas.PeakRecord("chr1", 100, 124, "p1", 10.0, 112)
    index = atlas.PeakIndex({"chr1": [peak]})
    assert index.match("chr1", 101, 124).fully_contained
    partial = index.match("chr1", 101, 125)
    assert partial.any_overlap
    assert not partial.fully_contained


def test_candidate_rank_signal_ignores_unsupported_contexts():
    supported = atlas.PeakMatch(True, True, 24, "p1", 7.5, 2)
    unsupported = atlas.PeakMatch(False, False, 0, "", None, None)
    hit = atlas.GuideHit(
        "Gene", "Un1Cas12f1_TTTR", "chr1", "+", 100, 104, "TTTA",
        104, 124, "ACGT" * 5, 0.5, "1.000", -200,
        {state: supported if state == atlas.STATES[0] else unsupported for state in atlas.STATES},
    )
    assert atlas.minimum_primary_signal(hit) == 7.5


def test_matched_depth_allocation_preserves_total_and_respects_capacity():
    assert matched_depth.allocate_target([10, 10], 15) == [8, 7]
    assert matched_depth.allocate_target([3, 20], 15) == [3, 12]
    allocation = matched_depth.allocate_target([1, 2, 30], 20)
    assert sum(allocation) == 20
    assert all(value <= capacity for value, capacity in zip(allocation, [1, 2, 30], strict=True))


def test_run_level_peak_evidence_distinguishes_complete_and_partial_support():
    with tempfile.TemporaryDirectory() as directory:
        peaks = Path(directory) / "replicate.narrowPeak"
        peaks.write_text("chr1\t100\t124\tp1\t10\t.\t8.5\t-1\t-1\t12\n")
        index = replicate_evidence.PeakIndex(peaks)
        complete = index.evidence("chr1", 101, 124)
        partial = index.evidence("chr1", 101, 125)
        assert complete["guide_fully_in_replicate_peak"]
        assert complete["replicate_peak_signal"] == "8.500000"
        assert partial["guide_any_replicate_peak_overlap"]
        assert not partial["guide_fully_in_replicate_peak"]


def test_reciprocal_consensus_does_not_use_transitive_chain():
    # A--B and B--C pass 50% reciprocal overlap, but A and C do not overlap.
    # The result must retain two locally supported intervals instead of dropping
    # the locus after intersecting all three peaks.
    peaks = [
        consensus.Peak(1, "chr1", 0, 100, 10.0, 50),
        consensus.Peak(2, "chr1", 40, 140, 9.0, 90),
        consensus.Peak(3, "chr1", 80, 180, 8.0, 130),
    ]
    intersections = consensus.qualifying_pair_intersections(peaks, 0.50)
    merged = consensus.merge_intersections(intersections)
    assert [(row[1], row[2]) for row in merged] == [(40, 140)]
    assert merged[0][-1] == 3


def test_iupac_reverse_complement_for_tttr():
    assert atlas.reverse_complement_iupac("TTTR") == "YAAA"


def test_un1cas12f_scanning_reports_both_orientations_and_full_intervals():
    sequence = list("CG" * 175)
    sequence[10:14] = list("TTTA")
    sequence[14:34] = list("ACGT" * 5)
    sequence[80:100] = list("ACGT" * 5)
    sequence[100:104] = list("CAAA")  # reverse-complement PAM is TTTG
    promoter = atlas.Promoter(
        "chr1", 1000, 1350, "Gene", "+", 1400, "g", "t", "canonical", "test"
    )
    peak = atlas.PeakRecord("chr1", 1000, 1350, "p", 5.0, 1175)
    peaks = {state: atlas.PeakIndex({"chr1": [peak]}) for state in atlas.STATES}
    hits = atlas.scan_promoter(promoter, "".join(sequence), "Un1Cas12f1_TTTR", peaks)
    expected = {(hit.target_strand, hit.pam_seq, hit.target_start, hit.target_end) for hit in hits}
    assert ("+", "TTTA", 1010, 1034) in expected
    assert ("-", "TTTG", 1080, 1104) in expected
    assert all(hit.peak_matches[atlas.STATES[0]].fully_contained for hit in hits)


def test_tracked_blacklists_have_verified_source_and_valid_mm39_intervals():
    source = ROOT / "supplementary/table_S6_blacklist_mm10_original.bed"
    lifted = ROOT / "supplementary/table_S6_blacklist_mm39_lifted.bed"
    assert hashlib.sha256(source.read_bytes()).hexdigest() == (
        "9638bbeb4be8d99ddf56f1b70700f6a9336ce7f54d87032bd262465ecf3bfac7"
    )
    rows = [line.rstrip("\n").split("\t") for line in lifted.open() if line.strip()]
    assert len(rows) == 3360
    assert all(int(row[1]) < int(row[2]) for row in rows)


def test_workflow_temporary_paths_are_expanded_per_wildcard_job():
    snakefile = (ROOT / "workflow/Snakefile").read_text()
    assert 'temp=lambda wc: WORK / f"results/peaks/replicate/{wc.run}.name_sorted.bam"' in snakefile
    assert 'merged=lambda wc: WORK / f"results/bam/condition/{wc.condition}.pre_pool_dedup.bam"' in snakefile
    assert 'regrouped=lambda wc: WORK / f"results/bam/condition/{wc.condition}.regrouped.bam"' in snakefile
    assert 'temp=lambda wc: WORK / f"results/peaks/pooled/{wc.condition}.name_sorted.bam"' in snakefile


def test_candidate_resource_does_not_label_predictions_as_recommendations():
    source = (ROOT / "scripts/rebuild_atlas_strict_iupac.py").read_text()
    assert '"candidate_class"' in source
    assert '"recommendation_class"' not in source


def test_release_atlas_is_gzip_compressed_and_stream_readable():
    path = ROOT / "supplementary/table_S2_targetability_full.tsv.gz"
    assert path.read_bytes()[:2] == b"\x1f\x8b"
    with gzip.open(path, "rt") as handle:
        assert handle.readline().rstrip("\n").split("\t")[:4] == [
            "gene", "gene_id", "transcript_id", "tss",
        ]
