#!/usr/bin/env python3
"""
Figure 2: PAM availability is not the bottleneck
Panel A: Barplot — % therapeutic genes with ≥1 PAM, by Cas
Panel B: Heatmap — PAM counts per gene × Cas
Panel C: TFEB detail
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
FULL_FILE = ROOT / 'supplementary' / 'table_S2_targetability_full.tsv'
GENES_FILE = ROOT / 'supplementary' / 'table_S1_therapeutic_genes.csv'
OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'fig2.pdf'

CAS_ORDER = ['SpCas9', 'Nme2Cas9', 'SaCas9', 'SminiCRa_Un1Cas12f1', 'HEAL_Un1Cas12f1', 'CjCas9_MiniCAFE']
CAS_LABELS = {
    'HEAL_Un1Cas12f1': 'HEAL\n(TTTR)',
    'SminiCRa_Un1Cas12f1': 'SminiCRa\n(TTTR)',
    'SaCas9': 'SaCas9\n(NNGRRT)',
    'SpCas9': 'SpCas9\n(NGG)',
    'CjCas9_MiniCAFE': 'CjCas9\n(NNNVRYM)',
    'Nme2Cas9': 'Nme2Cas9\n(NNNNCC)',
}
COLORS = {
    'HEAL_Un1Cas12f1': '#E63946',
    'SminiCRa_Un1Cas12f1': '#F4A261',
    'SaCas9': '#2A9D8F',
    'SpCas9': '#264653',
    'CjCas9_MiniCAFE': '#E9C46A',
    'Nme2Cas9': '#606C38',
}
SINGLE_AAV = {'HEAL_Un1Cas12f1', 'SminiCRa_Un1Cas12f1', 'CjCas9_MiniCAFE'}


def main():
    # Load therapeutic genes
    therapeutic = []
    with open(GENES_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            therapeutic.append(row['gene_symbol'])

    # Load PAM counts from the canonical targetability matrix. PAM counts are
    # state-invariant, so the homeostatic row is sufficient.
    gene_cas_count = defaultdict(lambda: defaultdict(int))
    with open(FULL_FILE) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['is_therapeutic'] == 'True' and row['state'] == 'homeostatic':
                gene_cas_count[row['gene']][row['cas']] = int(row['pams_total_passing'])

    # Panel A: % therapeutic genes with ≥1 PAM
    cas_pct = {}
    for cas in CAS_ORDER:
        n_with_pam = sum(1 for g in therapeutic if gene_cas_count[g][cas] > 0)
        cas_pct[cas] = 100 * n_with_pam / len(therapeutic)

    # Panel B: Heatmap matrix
    matrix = np.zeros((len(therapeutic), len(CAS_ORDER)))
    for i, gene in enumerate(therapeutic):
        for j, cas in enumerate(CAS_ORDER):
            matrix[i, j] = gene_cas_count[gene][cas]

    # Create figure
    fig = plt.figure(figsize=(14, 10))
    gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.5], width_ratios=[1.2, 1], hspace=0.35, wspace=0.3)

    # Panel A
    ax_a = fig.add_subplot(gs[0, 0])
    x = range(len(CAS_ORDER))
    for i, cas in enumerate(CAS_ORDER):
        edge = 'black' if cas in SINGLE_AAV else '#888888'
        lw = 2.5 if cas in SINGLE_AAV else 1
        ax_a.bar(i, cas_pct[cas], color=COLORS[cas], edgecolor=edge, linewidth=lw)
        ax_a.text(i, cas_pct[cas] + 0.5, f"{cas_pct[cas]:.1f}%", ha='center', va='bottom',
                  fontsize=9, fontweight='bold')
    ax_a.set_xticks(x)
    ax_a.set_xticklabels([CAS_LABELS[c] for c in CAS_ORDER], fontsize=8)
    ax_a.set_ylabel('% therapeutic genes\nwith ≥1 passing PAM', fontsize=10)
    ax_a.set_ylim(85, 102)
    ax_a.set_title('A', fontsize=14, fontweight='bold', loc='left')
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)
    single_patch = mpatches.Patch(edgecolor='black', facecolor='lightgray', lw=2.5, label='Single-AAV')
    dual_patch = mpatches.Patch(edgecolor='#888', facecolor='lightgray', lw=1, label='Dual-AAV')
    ax_a.legend(handles=[single_patch, dual_patch], fontsize=8, loc='lower left')

    # Panel B: Heatmap
    ax_b = fig.add_subplot(gs[:, 1])
    im = ax_b.imshow(np.log2(matrix + 1), cmap='YlOrRd', aspect='auto', interpolation='nearest')
    ax_b.set_xticks(range(len(CAS_ORDER)))
    ax_b.set_xticklabels([CAS_LABELS[c].replace('\n', ' ') for c in CAS_ORDER], fontsize=7, rotation=45, ha='right')
    ax_b.set_yticks(range(len(therapeutic)))
    ax_b.set_yticklabels(therapeutic, fontsize=5)
    ax_b.set_title('B', fontsize=14, fontweight='bold', loc='left')
    cbar = plt.colorbar(im, ax=ax_b, shrink=0.5, label='log2(PAM count + 1)')

    # Panel C: TFEB detail
    ax_c = fig.add_subplot(gs[1, 0])
    tfeb_idx = therapeutic.index('Tfeb') if 'Tfeb' in therapeutic else 0
    tfeb_counts = [gene_cas_count['Tfeb'][cas] for cas in CAS_ORDER]
    bars = ax_c.bar(range(len(CAS_ORDER)), tfeb_counts, color=[COLORS[c] for c in CAS_ORDER],
                    edgecolor='black', linewidth=0.5)
    for i, v in enumerate(tfeb_counts):
        ax_c.text(i, v + 0.3, str(v), ha='center', fontsize=9, fontweight='bold')
    ax_c.set_xticks(range(len(CAS_ORDER)))
    ax_c.set_xticklabels([CAS_LABELS[c] for c in CAS_ORDER], fontsize=8)
    ax_c.set_ylabel('Passing PAMs in\nCRISPRa window', fontsize=10)
    ax_c.set_title('C  TFEB PAM availability', fontsize=12, fontweight='bold', loc='left')
    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
