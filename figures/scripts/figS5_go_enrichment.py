#!/usr/bin/env python3
"""
Supplementary Figure S5: GO enrichment by accessibility category
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'figS5.pdf'

# GO results from Enrichr (actual data)
GO_DATA = {
    'Constitutive\n(21 genes)': [
        ('Pos. reg. transcription RNA Pol II', 0.0001),
        ('Embryonic digit morphogenesis', 0.0002),
        ('Reg. macrophage foam cell diff.', 0.0006),
        ('Pos. reg. nitrogen compound metab.', 0.0006),
    ],
    'Inflammation-\ngained (8)': [
        ('Microglial cell activation', 0.00001),
        ('Reg. microglial cell migration', 0.0004),
        ('Synapse pruning', 0.0004),
        ('Macrophage activation (immune)', 0.0006),
    ],
    'Surgical/stroke\ncontext (13)': [
        ('Pos. reg. cellular catabolic process', 0.00001),
        ('Pos. reg. MAPK cascade', 0.00001),
        ('Pos. reg. STAT tyrosine phosph.', 0.00001),
        ('Reg. STAT tyrosine phosph.', 0.00001),
    ],
    'Never\naccessible (12)': [
        ('Pos. reg. miRNA transcription', 0.0006),
        ('Neg. reg. cytokine signaling', 0.0006),
        ('Neg. reg. heterotypic cell adhesion', 0.0009),
        ('Pos. reg. miRNA metabolic process', 0.0006),
    ],
}

CAT_COLORS = {
    'Constitutive\n(21 genes)': '#2A9D8F',
    'Inflammation-\ngained (8)': '#F4A261',
    'Surgical/stroke\ncontext (13)': '#6A994E',
    'Never\naccessible (12)': '#D3D3D3',
}


def main():
    fig, axes = plt.subplots(2, 2, figsize=(16, 10))
    axes = axes.flatten()

    for i, (cat, terms) in enumerate(GO_DATA.items()):
        ax = axes[i]
        names = [t[0] for t in terms]
        logp = [-np.log10(t[1]) for t in terms]
        color = CAT_COLORS[cat]

        y = range(len(terms))
        ax.barh(y, logp, color=color, edgecolor='black', linewidth=0.5)
        ax.set_yticks(y)
        ax.set_yticklabels(names, fontsize=8)
        ax.set_xlabel('-log10(p-value)', fontsize=9)
        ax.set_title(cat, fontsize=11, fontweight='bold')
        ax.invert_yaxis()
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    plt.suptitle('GO Biological Process enrichment by accessibility category',
                 fontsize=13, fontweight='bold', y=1.02)
    plt.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
