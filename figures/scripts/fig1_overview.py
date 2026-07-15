#!/usr/bin/env python3
"""
Figure 1: Overview of the CRISPRa targetability atlas and pipeline
Panel A: Pipeline schematic
Panel B: Atlas summary statistics
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from pathlib import Path

OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'fig1.pdf'

BOX_STYLE = "round,pad=0.2"


def draw_box(ax, x, y, w, h, text, color, fontsize=8, textcolor='white'):
    bbox = FancyBboxPatch((x - w/2, y - h/2), w, h,
                           boxstyle=BOX_STYLE, facecolor=color,
                           edgecolor='black', linewidth=1.2)
    ax.add_patch(bbox)
    ax.text(x, y, text, ha='center', va='center', fontsize=fontsize,
            fontweight='bold', color=textcolor, linespacing=1.3)


def draw_arrow(ax, x1, y1, x2, y2):
    ax.annotate('', xy=(x2, y2), xytext=(x1, y1),
                arrowprops=dict(arrowstyle='->', color='#333333', lw=2))


def main():
    fig = plt.figure(figsize=(16, 7))
    gs = fig.add_gridspec(1, 2, width_ratios=[2, 1], wspace=0.15)

    # === Panel A: Pipeline ===
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.axis('off')
    ax_a.set_xlim(-1, 11)
    ax_a.set_ylim(-0.5, 6.5)

    # Phase 1 boxes
    draw_box(ax_a, 1.5, 5.5, 2.8, 0.9, 'GENCODE vM33\n21,599 genes', '#457B9D')
    draw_box(ax_a, 5.0, 5.5, 2.8, 0.9, 'Promoter windows\nTSS -400 to -50 bp', '#457B9D')
    draw_box(ax_a, 8.5, 5.5, 2.8, 0.9, 'PAM scanning\n6 Cas orthologs\nboth strands', '#264653')
    draw_arrow(ax_a, 2.9, 5.5, 3.6, 5.5)
    draw_arrow(ax_a, 6.4, 5.5, 7.1, 5.5)

    # Phase 1 label
    ax_a.text(5.0, 6.2, 'PHASE 1: PAM Availability', fontsize=11, fontweight='bold',
              ha='center', color='#264653')

    # Phase 2 boxes
    draw_box(ax_a, 1.5, 3.5, 2.8, 0.9, 'ATAC-seq\n3 datasets, 3 labs\n6 microglial states', '#E76F51')
    draw_box(ax_a, 5.0, 3.5, 2.8, 0.9, 'Trim → Align → Filter\nBowtie2, MAPQ≥30\nchrM, blacklist', '#E76F51')
    draw_box(ax_a, 8.5, 3.5, 2.8, 0.9, 'Peak calling\nGenrich\nATAC-seq mode', '#E63946')
    draw_arrow(ax_a, 2.9, 3.5, 3.6, 3.5)
    draw_arrow(ax_a, 6.4, 3.5, 7.1, 3.5)

    ax_a.text(5.0, 4.2, 'PHASE 2: Chromatin Accessibility', fontsize=11, fontweight='bold',
              ha='center', color='#E63946')

    # Integration
    draw_arrow(ax_a, 8.5, 5.0, 5.0, 2.3)
    draw_arrow(ax_a, 8.5, 3.0, 5.0, 2.3)
    draw_box(ax_a, 5.0, 1.5, 4.0, 1.2,
             'INTEGRATION\nTargetability = PAM × Chromatin\nper gene × Cas × state', '#1D3557', fontsize=9)

    # Output
    draw_arrow(ax_a, 5.0, 0.9, 5.0, 0.2)
    draw_box(ax_a, 5.0, -0.2, 4.0, 0.6,
             'CRISPRa Targetability Atlas + sgRNA recommendations', '#2A9D8F', fontsize=8)

    ax_a.set_title('A  Computational pipeline', fontsize=13, fontweight='bold', loc='left')

    # === Panel B: Atlas numbers ===
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.axis('off')

    stats = [
        ('21,599', 'protein-coding genes'),
        ('6', 'Cas orthologs'),
        ('6', 'microglial states'),
        ('3', 'independent laboratories'),
        ('55', 'curated therapeutic genes'),
        ('2.72M', 'passing PAM candidates'),
        ('~295K', 'ATAC-seq peaks (total)'),
        ('100%', 'publicly available data'),
    ]

    y_start = 0.92
    for i, (number, label) in enumerate(stats):
        y = y_start - i * 0.115
        ax_b.text(0.15, y, number, transform=ax_b.transAxes, fontsize=18,
                  fontweight='bold', color='#264653', ha='right', va='center')
        ax_b.text(0.2, y, label, transform=ax_b.transAxes, fontsize=10,
                  color='#333333', ha='left', va='center')

    ax_b.set_title('B  Atlas summary', fontsize=13, fontweight='bold', loc='left')

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
