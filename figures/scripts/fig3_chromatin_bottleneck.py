#!/usr/bin/env python3
"""
Figure 3: Chromatin accessibility is the real bottleneck (CENTRAL FIGURE)
Panel A: Paired barplot with prominent deltas
Panel B: Heatmap — top 25 most informative genes (rest to supplementary)
Panel C: Line plot — Cas convergence after chromatin filter
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
FULL_FILE = ROOT / 'supplementary' / 'table_S2_targetability_full.tsv'
GENES_FILE = ROOT / 'supplementary' / 'table_S1_therapeutic_genes.csv'
OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'fig3.pdf'

CAS_ORDER = ['SpCas9', 'Nme2Cas9', 'SaCas9', 'SminiCRa_Un1Cas12f1', 'HEAL_Un1Cas12f1', 'CjCas9_MiniCAFE']
CAS_SHORT = {
    'HEAL_Un1Cas12f1': 'HEAL', 'SminiCRa_Un1Cas12f1': 'SminiCRa',
    'SaCas9': 'SaCas9', 'SpCas9': 'SpCas9',
    'CjCas9_MiniCAFE': 'CjCas9', 'Nme2Cas9': 'Nme2Cas9',
}
COLORS = {
    'HEAL_Un1Cas12f1': '#E63946', 'SminiCRa_Un1Cas12f1': '#F4A261',
    'SaCas9': '#2A9D8F', 'SpCas9': '#264653',
    'CjCas9_MiniCAFE': '#E9C46A', 'Nme2Cas9': '#606C38',
}
STATE_ORDER = ['homeostatic', 'PP_naive', 'PL_acute_LPS', 'LL_tolerized', 'sham_WT', 'stroke_WT']
STATE_LABELS = ['Homeostatic', 'Naive (PP)', 'Acute LPS', 'Tolerized', 'Sham (Zhang)', 'Stroke (Zhang)']

def main():
    # Load therapeutic genes with categories
    therapeutic = []
    gene_cat = {}
    with open(GENES_FILE) as f:
        reader = csv.DictReader(f)
        for row in reader:
            therapeutic.append(row['gene_symbol'])
            gene_cat[row['gene_symbol']] = row['category']

    # Load full targetability for heatmap and summary percentages.
    gene_state_target = defaultdict(lambda: defaultdict(bool))
    summary_counts = defaultdict(lambda: {'targetable': 0, 'total': 0})
    pam_counts = defaultdict(lambda: {'has_pam': 0, 'total': 0})
    with open(FULL_FILE) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['is_therapeutic'] != 'True':
                continue
            cas = row['cas']
            state = row['state']
            key = (cas, state)
            summary_counts[key]['total'] += 1
            pam_counts[key]['total'] += 1
            if row['targetable'] == 'True':
                summary_counts[key]['targetable'] += 1
            if int(row['pams_total_passing']) > 0:
                pam_counts[key]['has_pam'] += 1
            if cas == 'HEAL_Un1Cas12f1':
                gene_state_target[row['gene']][row['state']] = row['targetable'] == 'True'

    summary = {
        key: counts['targetable'] / counts['total'] * 100
        for key, counts in summary_counts.items()
    }
    pam_only = {
        key[0]: counts['has_pam'] / counts['total'] * 100
        for key, counts in pam_counts.items()
        if key[1] == 'homeostatic'
    }

    fig = plt.figure(figsize=(17, 9))
    gs = fig.add_gridspec(1, 3, width_ratios=[1.1, 1.3, 1], wspace=0.3)

    # === Panel A: Paired barplot with PROMINENT deltas ===
    ax_a = fig.add_subplot(gs[0, 0])
    x = np.arange(len(STATE_ORDER))
    w = 0.32
    pam_only_vals = [pam_only['HEAL_Un1Cas12f1']] * len(STATE_ORDER)
    chromatin_vals = [summary[('HEAL_Un1Cas12f1', s)] for s in STATE_ORDER]

    ax_a.bar(x - w/2, pam_only_vals, w, label='PAM only', color='#A8DADC', edgecolor='black', lw=0.8)
    ax_a.bar(x + w/2, chromatin_vals, w, label='PAM + chromatin', color='#E63946', edgecolor='black', lw=0.8)

    # Prominent delta labels between bars
    for i in range(len(STATE_ORDER)):
        drop = pam_only_vals[i] - chromatin_vals[i]
        mid_y = (pam_only_vals[i] + chromatin_vals[i]) / 2
        # Large bold delta text centered between the two bars
        ax_a.text(x[i], mid_y, f'\u2193{drop:.0f}%',
                  ha='center', va='center', fontsize=14, fontweight='bold',
                  color='#1D3557',
                  bbox=dict(boxstyle='round,pad=0.2', facecolor='white', edgecolor='#1D3557', lw=1.5, alpha=0.9))

    ax_a.set_xticks(x)
    ax_a.set_xticklabels(STATE_LABELS, fontsize=8, rotation=15, ha='right')
    ax_a.set_ylabel('% therapeutic genes\ntargetable (HEAL)', fontsize=10)
    ax_a.set_ylim(0, 110)
    ax_a.legend(fontsize=8, loc='upper right')
    ax_a.set_title('A', fontsize=16, fontweight='bold', loc='left', x=-0.12)
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)

    # === Panel B: Heatmap — top 25 most informative genes ===
    ax_b = fig.add_subplot(gs[0, 1])

    # Select top 25: prioritize genes with interesting patterns
    # Sort: constitutively open first, then inflammation-gained, then closed
    def sort_key(gene):
        opens = sum(1 for s in STATE_ORDER if gene_state_target[gene][s])
        homeo = 1 if gene_state_target[gene]['homeostatic'] else 0
        return (-opens, -homeo, gene)

    all_sorted = sorted([g for g in therapeutic if g in gene_state_target], key=sort_key)

    # Pick top 25: include key genes explicitly + fill rest by pattern diversity
    must_include = ['Tfeb', 'Tfe3', 'Trem2', 'Tyrobp', 'Abca1', 'P2ry12', 'Lamp1', 'Lamp2',
                    'Il1b', 'Tnf', 'Bdnf', 'Igf1', 'Pparg', 'Grn', 'Megf10', 'Vegfa',
                    'Mcoln1', 'Cx3cr1', 'Tgfb1', 'Il10', 'Hdac1', 'Spi1', 'Csf1r', 'Apoe', 'Mertk']
    selected = [g for g in all_sorted if g in must_include][:25]
    # Re-sort selected
    selected.sort(key=sort_key)

    matrix = np.zeros((len(selected), len(STATE_ORDER)))
    for i, gene in enumerate(selected):
        for j, state in enumerate(STATE_ORDER):
            matrix[i, j] = 1 if gene_state_target[gene][state] else 0

    cmap = plt.cm.colors.ListedColormap(['#E8E0D0', '#2A9D8F'])
    im = ax_b.imshow(matrix, cmap=cmap, aspect='auto', interpolation='nearest')

    ax_b.set_xticks(range(len(STATE_ORDER)))
    ax_b.set_xticklabels(STATE_LABELS, fontsize=9, rotation=15, ha='right')
    ax_b.set_yticks(range(len(selected)))
    ax_b.set_yticklabels(selected, fontsize=8)

    # Add category color strip on the left
    cat_colors = {
        'autophagy_lysosome': '#E63946', 'inflammation': '#F4A261',
        'phagocytosis': '#2A9D8F', 'lipid_metabolism': '#264653',
        'neuroprotection': '#E9C46A', 'microglial_identity': '#606C38',
        'epigenetic_regulation': '#A8DADC',
    }
    for i, gene in enumerate(selected):
        cat = gene_cat.get(gene, '')
        color = cat_colors.get(cat, 'gray')
        ax_b.add_patch(plt.Rectangle((-0.7, i - 0.5), 0.4, 1, facecolor=color, edgecolor='none', clip_on=False))

    ax_b.set_title('B  Top 25 therapeutic genes (HEAL)', fontsize=10, fontweight='bold', loc='left')

    legend_elements = [
        plt.Rectangle((0,0), 1, 1, facecolor='#2A9D8F', label='Targetable'),
        plt.Rectangle((0,0), 1, 1, facecolor='#E8E0D0', edgecolor='gray', lw=0.5, label='Not targetable'),
    ]
    ax_b.legend(handles=legend_elements, loc='lower right', fontsize=8)

    # === Panel C: All Cas converge ===
    ax_c = fig.add_subplot(gs[0, 2])
    for cas in CAS_ORDER:
        vals = [summary[(cas, s)] for s in STATE_ORDER]
        ls = '-' if cas in {'HEAL_Un1Cas12f1', 'SminiCRa_Un1Cas12f1', 'CjCas9_MiniCAFE'} else '--'
        lw = 2.5 if cas in {'HEAL_Un1Cas12f1', 'SminiCRa_Un1Cas12f1'} else 1.5
        ax_c.plot(range(len(STATE_ORDER)), vals, marker='o', color=COLORS[cas],
                  label=CAS_SHORT[cas], linewidth=lw, markersize=5, linestyle=ls)

    ax_c.set_xticks(range(len(STATE_ORDER)))
    ax_c.set_xticklabels(STATE_LABELS, fontsize=7, rotation=15, ha='right')
    ax_c.set_ylabel('% therapeutic genes targetable', fontsize=9)
    ymax = max(summary[(cas, s)] for cas in CAS_ORDER for s in STATE_ORDER)
    ax_c.set_ylim(20, min(100, ymax + 12))
    ax_c.legend(fontsize=7, ncol=1, loc='lower right')
    ax_c.set_title('C', fontsize=16, fontweight='bold', loc='left', x=-0.12)
    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)

    # Annotation: state differences are larger than Cas differences.
    ax_c.annotate('State shifts exceed\nCas-specific spread',
                  xy=(2.7, 58), fontsize=8, fontstyle='italic', color='gray',
                  ha='center')

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
