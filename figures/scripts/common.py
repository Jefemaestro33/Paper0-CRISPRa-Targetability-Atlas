"""Shared plotting constants and loaders for the revised targetability figures."""
import json
from pathlib import Path

import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
OUTPUT = ROOT / "figures/output"
ATLAS = ROOT / "supplementary/table_S2_targetability_full.tsv.gz"
PANEL = ROOT / "config/therapeutic_genes_locked.csv"
CANDIDATES = ROOT / "supplementary/table_S3_candidate_protospacers.csv"
KEY_COUNTS = ROOT / "analysis_stats/key_counts.json"
STATES = ["homeostatic", "PP_control", "PL_acute_LPS", "LL_tolerized", "sham_WT", "stroke_WT"]
STATE_LABELS = ["Gosselin\nhomeostatic", "Zhang X\nPBS/PBS", "Zhang X\nPBS/LPS", "Zhang X\nLPS/LPS", "Zhang L\nsham", "Zhang L\nstroke"]
CAS_ORDER = ["Un1Cas12f1_TTTR", "SaCas9_NNGRRT", "SpCas9_NGG", "CjCas9_NNNVRYM", "Nme2Cas9_NNNNCC"]
CAS_LABELS = {
    "Un1Cas12f1_TTTR": "Un1Cas12f1\nTTTR\n(HEAL/SminiCRa)",
    "SaCas9_NNGRRT": "SaCas9\nNNGRRT",
    "SpCas9_NGG": "SpCas9\nNGG",
    "CjCas9_NNNVRYM": "CjCas9\nNNNVRYM",
    "Nme2Cas9_NNNNCC": "Nme2Cas9\nNNNNCC",
}
COLORS = {
    "Un1Cas12f1_TTTR": "#D1495B", "SaCas9_NNGRRT": "#2A9D8F",
    "SpCas9_NGG": "#264653", "CjCas9_NNNVRYM": "#E9C46A",
    "Nme2Cas9_NNNNCC": "#577590",
}


def load_atlas(path=ATLAS):
    frame = pd.read_csv(path, sep="\t", low_memory=False)
    for column in ("targetable", "is_therapeutic", "promoter_midpoint_accessible", "promoter_any_peak_overlap"):
        if column in frame:
            frame[column] = frame[column].astype(str).str.lower().eq("true")
    return frame


def load_panel():
    return pd.read_csv(PANEL)


def load_key_counts():
    with KEY_COUNTS.open() as handle:
        return json.load(handle)


def save(fig, name):
    OUTPUT.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUTPUT / f"{name}.pdf", dpi=300, bbox_inches="tight")
    fig.savefig(OUTPUT / f"{name}.png", dpi=300, bbox_inches="tight")
