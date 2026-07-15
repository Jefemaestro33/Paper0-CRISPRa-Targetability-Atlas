#!/usr/bin/env python3
"""
Figure 5: Practical framework for state-aware CRISPRa vector design
Panel A: Decision tree — uniform styling
Panel B: TFE3 case study — sgRNAs without misleading score column
"""
import csv
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

ROOT = Path(__file__).resolve().parent.parent.parent
SGRNA_FILE = ROOT / 'supplementary' / 'table_S3_sgrna_recommendations.csv'
OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'fig5.pdf'

BOX_COLOR = '#264653'


def main():
    # Load TFE3 sgRNAs
    tfe3_sgrnas = []
    with open(SGRNA_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            if (
                row['gene_symbol'] == 'Tfe3'
                and row['cas_ortholog'] == 'HEAL_Un1Cas12f1'
                and row['atlas_n_targetable_states'] == '6'
            ):
                tfe3_sgrnas.append(row)

    fig = plt.figure(figsize=(16, 8))
    gs = fig.add_gridspec(1, 2, wspace=0.35, width_ratios=[1, 1.2])

    # === Panel A: Decision tree — uniform box style ===
    ax_a = fig.add_subplot(gs[0, 0])
    ax_a.axis('off')
    ax_a.set_xlim(0, 10)
    ax_a.set_ylim(0, 10)

    boxes = [
        (5, 8.8, 'Step 1\nSelect Cas by packaging\nconstraints (NOT by PAM)'),
        (5, 6.8, 'Step 2\nVerify promoter accessibility\nvia ATAC-seq in target cell state'),
        (5, 4.8, 'Step 3\nDesign sgRNAs within\nATAC-seq peaks'),
        (5, 2.8, 'Step 4\nConsider state-dependence\nof target accessibility'),
        (5, 0.8, 'Output\nState-aware CRISPRa vector'),
    ]

    for x, y, text in boxes:
        bbox = FancyBboxPatch((x - 2.5, y - 0.75), 5.0, 1.5,
                               boxstyle="round,pad=0.15", facecolor=BOX_COLOR, alpha=0.9,
                               edgecolor='white', linewidth=2)
        ax_a.add_patch(bbox)
        ax_a.text(x, y, text, ha='center', va='center', fontsize=9,
                  fontweight='bold', color='white', linespacing=1.3)

    for i in range(len(boxes) - 1):
        y_from = boxes[i][1] - 0.75
        y_to = boxes[i+1][1] + 0.75
        ax_a.annotate('', xy=(5, y_to), xytext=(5, y_from),
                       arrowprops=dict(arrowstyle='->', color='#264653', lw=2.5))

    # Side annotations
    ax_a.text(8.5, 8.8, 'HEAL (1.5 kb)\nSminiCRa (1.6 kb)\nfor single-AAV',
              fontsize=7, color='gray', fontstyle='italic', va='center')
    ax_a.text(8.5, 6.8, 'This atlas or\nequivalent ATAC-seq',
              fontsize=7, color='gray', fontstyle='italic', va='center')
    ax_a.text(8.5, 4.8, 'PAMs in peaks,\nnot just in promoter',
              fontsize=7, color='gray', fontstyle='italic', va='center')
    ax_a.text(8.5, 2.8, 'Constitutive vs\ninflammation-gained',
              fontsize=7, color='gray', fontstyle='italic', va='center')

    ax_a.set_title('A  Decision framework for in vivo CRISPRa design',
                    fontsize=11, fontweight='bold', loc='left')

    # === Panel B: TFE3 sgRNA table — no score column ===
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.axis('off')

    if tfe3_sgrnas:
        col_labels = ['#', 'Protospacer (5\' to 3\')', 'Strand', 'GC%', 'Chromatin']
        table_data = []
        for i, sg in enumerate(tfe3_sgrnas[:5], 1):
            gc_pct = f"{float(sg['gc_content'])*100:.0f}%"
            table_data.append([
                str(i),
                sg['protospacer_sequence'],
                '+' if sg['strand'] == '+' else '\u2212',
                gc_pct,
                'Open (all states)',
            ])

        table = ax_b.table(cellText=table_data, colLabels=col_labels,
                           cellLoc='center', loc='upper center',
                           colColours=['#264653'] * 5)
        table.auto_set_font_size(False)
        table.set_fontsize(9)
        table.scale(1, 1.8)

        # Style
        for j in range(len(col_labels)):
            table[0, j].set_text_props(color='white', fontweight='bold')
        for i in range(1, len(table_data) + 1):
            for j in range(len(col_labels)):
                table[i, j].set_facecolor('#F4F1DE' if i % 2 == 0 else 'white')
            # Monospace for protospacer
            table[i, 1].set_text_props(fontfamily='monospace', fontsize=8)

    ax_b.set_title('B  TFE3 sgRNA candidates: HEAL (TTTR) system\n'
                    '     Constitutively accessible across all microglial states',
                    fontsize=11, fontweight='bold', loc='left')

    ax_b.text(0.5, 0.15,
              'HEAL system: dUn1Cas12f1 + activation domain in single AAV\n'
              'All candidates in ATAC-seq peaks across all 6 microglial states (3 labs)\n'
              'Experimental validation required (heuristic pre-selection, not ML-predicted efficacy)',
              transform=ax_b.transAxes, ha='center', fontsize=8, fontstyle='italic', color='gray',
              linespacing=1.5)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
