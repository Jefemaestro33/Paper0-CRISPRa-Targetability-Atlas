#!/usr/bin/env python3
"""
Supplementary Figure S2: Genome-wide targetability (all 21,599 genes)
Panel A: Distribution of PAM counts per gene per Cas
Panel B: % of all genes targetable (PAM + chromatin) per state
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
OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'figS2.pdf'

CAS_ORDER = ['SpCas9', 'Nme2Cas9', 'SaCas9', 'SminiCRa_Un1Cas12f1', 'HEAL_Un1Cas12f1', 'CjCas9_MiniCAFE']
CAS_SHORT = {'HEAL_Un1Cas12f1': 'HEAL', 'SminiCRa_Un1Cas12f1': 'SminiCRa',
             'SaCas9': 'SaCas9', 'SpCas9': 'SpCas9',
             'CjCas9_MiniCAFE': 'CjCas9', 'Nme2Cas9': 'Nme2Cas9'}
COLORS = {'HEAL_Un1Cas12f1': '#E63946', 'SminiCRa_Un1Cas12f1': '#F4A261',
          'SaCas9': '#2A9D8F', 'SpCas9': '#264653',
          'CjCas9_MiniCAFE': '#E9C46A', 'Nme2Cas9': '#606C38'}
STATE_ORDER = ['homeostatic', 'PP_naive', 'PL_acute_LPS', 'LL_tolerized', 'sham_WT', 'stroke_WT']
STATE_LABELS = ['Homeostatic', 'Naïve', 'Acute LPS', 'Tolerized', 'Sham', 'Stroke']


def main():
    # Load PAM counts and genome-wide targetability from the canonical matrix.
    gene_cas_count = defaultdict(lambda: defaultdict(int))
    gene_state_target = defaultdict(lambda: defaultdict(lambda: defaultdict(bool)))
    with open(FULL_FILE) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['state'] == 'homeostatic':
                gene_cas_count[row['gene']][row['cas']] = int(row['pams_total_passing'])
            gene_state_target[row['gene']][row['cas']][row['state']] = row['targetable'] == 'True'

    fig, axes = plt.subplots(1, 2, figsize=(14, 6))

    # Panel A: Violin/box of PAM counts per Cas
    ax_a = axes[0]
    data_for_box = []
    labels = []
    colors_list = []
    for cas in CAS_ORDER:
        counts = [gene_cas_count[g][cas] for g in gene_state_target]
        data_for_box.append(counts)
        labels.append(CAS_SHORT[cas])
        colors_list.append(COLORS[cas])

    bp = ax_a.boxplot(data_for_box, tick_labels=labels, patch_artist=True, showfliers=False,
                       medianprops=dict(color='black', lw=2))
    for patch, color in zip(bp['boxes'], colors_list):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)

    ax_a.set_ylabel('Passing PAMs per gene\n(CRISPRa optimal window)', fontsize=10)
    ax_a.set_title('A  PAM count distribution (genome-wide)', fontsize=11, fontweight='bold', loc='left')
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)

    # Panel B: % genome-wide genes targetable per Cas per state
    ax_b = axes[1]
    all_genes = set(gene_state_target.keys())
    n_total = len(all_genes)

    x = np.arange(len(STATE_ORDER))
    width = 0.12
    for i, cas in enumerate(CAS_ORDER):
        vals = []
        for state in STATE_ORDER:
            n_target = sum(1 for g in all_genes if gene_state_target[g][cas][state])
            vals.append(100 * n_target / n_total if n_total > 0 else 0)
        offset = (i - len(CAS_ORDER)/2 + 0.5) * width
        ax_b.bar(x + offset, vals, width, label=CAS_SHORT[cas], color=COLORS[cas],
                 edgecolor='black', linewidth=0.3)

    ax_b.set_xticks(x)
    ax_b.set_xticklabels(STATE_LABELS, fontsize=9)
    ax_b.set_ylabel('% of all 21,599 genes targetable', fontsize=10)
    ax_b.set_title('B  Genome-wide targetability by state', fontsize=11, fontweight='bold', loc='left')
    ax_b.legend(fontsize=7, ncol=2)
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)

    plt.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
