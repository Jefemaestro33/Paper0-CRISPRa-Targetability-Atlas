#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from common import load_key_counts, save


def box(ax, x, y, text, color, width=2.55, height=0.92, fontsize=5.4):
    patch = FancyBboxPatch((x-width/2, y-height/2), width, height, boxstyle="round,pad=0.08", facecolor=color, edgecolor="#222", linewidth=1)
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color="white", weight="bold", linespacing=1.05)


def arrow(ax, start, end):
    ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 1.6, "color": "#333"})


def main():
    counts = load_key_counts()
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.35), gridspec_kw={"width_ratios": [2.35, 1]})
    ax = axes[0]
    ax.set_xlim(-0.35, 12.25); ax.set_ylim(0, 7); ax.axis("off")
    box(ax, 1.4, 5.7, f"GENCODE vM33\n{counts['n_genes']:,} genes", "#457B9D")
    box(ax, 4.45, 5.7, "Canonical TSS primary\nAPPRIS + legacy\nsensitivity", "#457B9D")
    box(ax, 7.5, 5.7, "−400 to −50 bp\nboth strands", "#457B9D")
    box(ax, 10.55, 5.7, f"{counts['n_targeting_classes']} nuclease/PAM classes\n{counts['n_crispra_configurations']} CRISPRa\nconfigurations", "#264653")
    for x in (2.72, 5.77, 8.82): arrow(ax, (x, 5.7), (x+0.4, 5.7))
    box(ax, 1.4, 3.7, f"{counts['n_runs']} public\nATAC-seq runs\n{counts['n_studies']} studies; {counts['n_contexts']} contexts", "#E76F51")
    box(ax, 4.45, 3.7, "Uniform trim + mm39\nalignment\ndeduplicate/filter", "#E76F51")
    box(ax, 7.5, 3.7, "Replicate consensus\ntechnical pools\nexplicitly labelled", "#C8553D")
    box(ax, 10.55, 3.7, "Matched depth\nGenrich + MACS3", "#C8553D")
    for x in (2.72, 5.77, 8.82): arrow(ax, (x, 3.7), (x+0.4, 3.7))
    arrow(ax, (8.7, 5.25), (6.6, 2.35)); arrow(ax, (8.7, 3.25), (6.6, 2.35))
    box(ax, 6, 1.9, "Primary call\ncomplete protospacer + PAM inside a primary peak\nbiological consensus or labelled technical pool", "#1D3557", width=6.5, height=1.18, fontsize=5.3)
    arrow(ax, (6, 1.35), (6, 0.75))
    box(ax, 6, 0.48, "Genome-wide targetability matrix\n+ guide-specific candidates\nrobustness/provenance tables", "#2A9D8F", width=6.5, height=0.86, fontsize=5.25)
    ax.set_title("A  Corrected analysis architecture", loc="left", weight="bold", fontsize=9)

    ax = axes[1]; ax.axis("off")
    stats = [
        (f"{counts['n_genes']:,}", "protein-coding promoters"),
        (str(counts["n_targeting_classes"]), "nuclease/PAM classes"),
        (str(counts["n_crispra_configurations"]), "CRISPRa configurations"),
        (str(counts["n_contexts"]), "surveyed dataset contexts"),
        (str(counts["n_runs"]), "raw sequencing runs"),
        (str(counts["n_panel"]), "locked panel genes"),
        (str(counts["n_tss_definitions"]), "TSS definitions"),
        (str(counts["n_matched_depth_peak_callers"]), "matched-depth peak callers"),
    ]
    for i, (number, label) in enumerate(stats):
        y = 0.93 - i*0.115
        ax.text(0.27, y, number, transform=ax.transAxes, ha="right", va="center", fontsize=12, weight="bold", color="#264653")
        ax.text(0.31, y, label, transform=ax.transAxes, ha="left", va="center", fontsize=6.8)
    ax.set_title("B  Scope", loc="left", weight="bold", fontsize=9)
    fig.subplots_adjust(left=.015, right=.995, top=.91, bottom=.035, wspace=.12)
    save(fig, "fig1")


if __name__ == "__main__": main()
