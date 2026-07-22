#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from common import save


def box(ax, x, y, text, color, width=2.55, height=0.92, fontsize=5.8):
    patch = FancyBboxPatch((x-width/2, y-height/2), width, height, boxstyle="round,pad=0.08", facecolor=color, edgecolor="#222", linewidth=1)
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=fontsize, color="white", weight="bold", linespacing=1.05)


def arrow(ax, start, end):
    ax.annotate("", xy=end, xytext=start, arrowprops={"arrowstyle": "->", "lw": 1.6, "color": "#333"})


def main():
    fig, axes = plt.subplots(1, 2, figsize=(7.4, 4.35), gridspec_kw={"width_ratios": [2.35, 1]})
    ax = axes[0]
    ax.set_xlim(-0.1, 12.1); ax.set_ylim(0, 7); ax.axis("off")
    box(ax, 1.4, 5.7, "GENCODE vM33\n21,599 genes", "#457B9D")
    box(ax, 4.45, 5.7, "Canonical TSS primary\nAPPRIS + legacy\nsensitivity", "#457B9D")
    box(ax, 7.5, 5.7, "−400 to −50 bp\nboth strands", "#457B9D")
    box(ax, 10.55, 5.7, "5 nuclease/PAM classes\n6 CRISPRa\nconfigurations", "#264653")
    for x in (2.72, 5.77, 8.82): arrow(ax, (x, 5.7), (x+0.4, 5.7))
    box(ax, 1.4, 3.7, "13 public ATAC-seq runs\n3 studies\n6 dataset contexts", "#E76F51")
    box(ax, 4.45, 3.7, "Uniform trim + mm39\nalignment\ndeduplicate + filter", "#E76F51")
    box(ax, 7.5, 3.7, "Replicate consensus\ntechnical pools\nexplicitly labelled", "#C8553D")
    box(ax, 10.55, 3.7, "Matched depth\nGenrich + MACS3", "#C8553D")
    for x in (2.72, 5.77, 8.82): arrow(ax, (x, 3.7), (x+0.4, 3.7))
    arrow(ax, (8.7, 5.25), (6.6, 2.35)); arrow(ax, (8.7, 3.25), (6.6, 2.35))
    box(ax, 6, 1.9, "Primary call\ncomplete protospacer + PAM inside a primary peak\nbiological consensus or explicitly labelled technical pool", "#1D3557", width=6.5, height=1.18, fontsize=5.8)
    arrow(ax, (6, 1.35), (6, 0.75))
    box(ax, 6, 0.48, "Genome-wide atlas + guide-specific candidates\nrobustness and provenance tables", "#2A9D8F", width=6.5, height=0.78, fontsize=5.8)
    ax.set_title("A  Corrected analysis architecture", loc="left", weight="bold", fontsize=9)

    ax = axes[1]; ax.axis("off")
    stats = [("21,599", "protein-coding promoters"), ("5", "nuclease/PAM classes"), ("6", "CRISPRa configurations"), ("6", "surveyed dataset contexts"), ("13", "raw sequencing runs"), ("55", "locked therapeutic genes"), ("3", "TSS definitions"), ("2", "matched-depth peak callers")]
    for i, (number, label) in enumerate(stats):
        y = 0.93 - i*0.115
        ax.text(0.27, y, number, transform=ax.transAxes, ha="right", va="center", fontsize=12, weight="bold", color="#264653")
        ax.text(0.31, y, label, transform=ax.transAxes, ha="left", va="center", fontsize=6.8)
    ax.set_title("B  Scope", loc="left", weight="bold", fontsize=9)
    fig.subplots_adjust(left=.015, right=.995, top=.91, bottom=.035, wspace=.12)
    save(fig, "fig1")


if __name__ == "__main__": main()
