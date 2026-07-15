#!/usr/bin/env python3
"""
Supplementary Figure S6: Cross-dataset comparison Gosselin vs Zhang sham
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'figS6.pdf'

# Actual data
JACCARD = 0.457
AGREE = 34
DISAGREE = 21
TOTAL = 55

# Discordant genes (20 closed→open, 1 open→closed)
CLOSED_TO_OPEN = 20
OPEN_TO_CLOSED = 1  # Megf10


def main():
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # Panel A: Concordance pie
    ax_a = axes[0]
    sizes = [AGREE, DISAGREE]
    colors = ['#2A9D8F', '#E63946']
    labels = [f'Concordant\n{AGREE}/{TOTAL} ({100*AGREE/TOTAL:.0f}%)',
              f'Discordant\n{DISAGREE}/{TOTAL} ({100*DISAGREE/TOTAL:.0f}%)']
    ax_a.pie(sizes, labels=labels, colors=colors, autopct='', startangle=90,
             textprops={'fontsize': 10, 'fontweight': 'bold'})
    ax_a.set_title('A  Gene-level concordance\n    (55 therapeutic genes)', fontsize=11, fontweight='bold')

    # Panel B: Direction of discordance
    ax_b = axes[1]
    dirs = [CLOSED_TO_OPEN, OPEN_TO_CLOSED]
    dir_labels = ['Gosselin CLOSED\n→ Zhang OPEN', 'Gosselin OPEN\n→ Zhang CLOSED']
    dir_colors = ['#F4A261', '#457B9D']
    bars = ax_b.bar(dir_labels, dirs, color=dir_colors, edgecolor='black', linewidth=0.8)
    for bar, val in zip(bars, dirs):
        ax_b.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.3,
                  str(val), ha='center', fontweight='bold', fontsize=14)
    ax_b.set_ylabel('Number of genes', fontsize=10)
    ax_b.set_title('B  Direction of discordance\n    (sham opens, does not close)', fontsize=11, fontweight='bold')
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)

    # Panel C: Summary text
    ax_c = axes[2]
    ax_c.axis('off')
    text = (
        f"Cross-dataset summary\n"
        f"{'='*35}\n\n"
        f"Jaccard index (1kb bins):  {JACCARD:.3f}\n"
        f"Gene concordance:          {AGREE}/{TOTAL} ({100*AGREE/TOTAL:.1f}%)\n\n"
        f"Discordant genes:          {DISAGREE}\n"
        f"  Closed→Open (sham opens): {CLOSED_TO_OPEN}\n"
        f"  Open→Closed (Megf10):     {OPEN_TO_CLOSED}\n\n"
        f"Interpretation:\n"
        f"Most discordance is Gosselin-closed\n"
        f"and Zhang-sham-open. This is compatible\n"
        f"with post-surgical activation, but may\n"
        f"also reflect protocol, depth, gating,\n"
        f"or batch differences.\n"
        f"The sole exception (Megf10) has a\n"
        f"unique pattern across all states."
    )
    ax_c.text(0.05, 0.95, text, transform=ax_c.transAxes, fontsize=9,
              verticalalignment='top', fontfamily='monospace',
              bbox=dict(boxstyle='round', facecolor='#F4F1DE', alpha=0.8))
    ax_c.set_title('C  Interpretation', fontsize=11, fontweight='bold', loc='left')

    plt.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
