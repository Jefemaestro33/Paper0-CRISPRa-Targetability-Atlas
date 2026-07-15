#!/usr/bin/env python3
"""
Figure 4: State-dependent chromatin dynamics
Panel A: Gene classification (constitutive / gained / lost / never)
Panel B: Genes gaining/losing accessibility
Panel C: Locus plots for TFEB, TFE3, TREM2
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
OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'fig4.pdf'

STATE_ORDER = ['homeostatic', 'PP_naive', 'PL_acute_LPS', 'LL_tolerized', 'sham_WT', 'stroke_WT']
STATE_LABELS = ['Homeostatic', 'Naïve', 'Acute LPS', 'Tolerized', 'Sham', 'Stroke']
STATE_COLORS = ['#457B9D', '#A8DADC', '#F4A261', '#E76F51', '#6A994E', '#BC4749']


def load_gene_accessibility(full_file, genes_file):
    therapeutic = []
    with open(genes_file) as f:
        reader = csv.DictReader(f)
        for row in reader:
            therapeutic.append(row['gene_symbol'])

    gene_state = defaultdict(lambda: defaultdict(bool))
    with open(full_file) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['is_therapeutic'] == 'True' and row['cas'] == 'HEAL_Un1Cas12f1':
                gene_state[row['gene']][row['state']] = row['promoter_accessible'] == 'True'

    return therapeutic, gene_state


def classify_genes(therapeutic, gene_state):
    constitutive, inflammation_gained, state_context, never, other_pattern = [], [], [], [], []
    for gene in therapeutic:
        acc = {s: gene_state[gene].get(s, False) for s in STATE_ORDER}
        all_open = all(acc.values())
        n_open = sum(acc.values())
        homeo = acc.get('homeostatic', False)
        sham = acc.get('sham_WT', False)
        stroke = acc.get('stroke_WT', False)

        if all_open:
            constitutive.append(gene)
        elif n_open == 0:
            never.append(gene)
        elif (sham or stroke) and not homeo and not all(acc[s] for s in ['PP_naive','PL_acute_LPS','LL_tolerized']):
            state_context.append(gene)
        elif any(acc[s] for s in STATE_ORDER[1:]) and not homeo:
            inflammation_gained.append(gene)
        else:
            other_pattern.append(gene)
    return constitutive, inflammation_gained, state_context, never, other_pattern


def main():
    therapeutic, gene_state = load_gene_accessibility(FULL_FILE, GENES_FILE)
    constitutive, gained, injury_cond, never, other = classify_genes(therapeutic, gene_state)

    fig = plt.figure(figsize=(16, 12))
    gs = fig.add_gridspec(2, 2, hspace=0.4, wspace=0.3)

    # Panel A: Bar chart of classification
    ax_a = fig.add_subplot(gs[0, 0])
    categories = ['Constitutively\nopen', 'Inflammation-\ngained', 'Surgical/stroke\ncontext', 'Never\naccessible', 'Other\npattern']
    counts = [len(constitutive), len(gained), len(injury_cond), len(never), len(other)]
    colors_pie = ['#2A9D8F', '#F4A261', '#6A994E', '#D3D3D3', '#457B9D']
    bars = ax_a.bar(categories, counts, color=colors_pie, edgecolor='black', linewidth=0.8)
    for bar, count in zip(bars, counts):
        ax_a.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                  str(count), ha='center', fontweight='bold', fontsize=12)
    ax_a.set_ylabel('Number of therapeutic genes', fontsize=11)
    ax_a.set_title('A  Accessibility classification', fontsize=12, fontweight='bold', loc='left')
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)

    # Panel B: Gene lists
    ax_b = fig.add_subplot(gs[0, 1])
    ax_b.axis('off')
    text = (
        f"Constitutively open ({len(constitutive)}):\n"
        f"  {', '.join(sorted(constitutive))}\n\n"
        f"Inflammation-gained ({len(gained)}):\n"
        f"  {', '.join(sorted(gained))}\n\n"
        f"Surgical/stroke-context ({len(injury_cond)}):\n"
        f"  {', '.join(sorted(injury_cond))}\n\n"
        f"Never accessible ({len(never)}):\n"
        f"  {', '.join(sorted(never))}"
        + (f"\n\nOther ({len(other)}): {', '.join(sorted(other))}" if other else "")
    )
    ax_b.text(0.05, 0.95, text, transform=ax_b.transAxes, fontsize=8,
              verticalalignment='top', fontfamily='monospace',
              bbox=dict(boxstyle='round', facecolor='#F4F1DE', alpha=0.8))
    ax_b.set_title('B  Gene lists by category', fontsize=12, fontweight='bold', loc='left')

    # Panel C: TFEB vs TFE3 vs TREM2 accessibility across states
    ax_c = fig.add_subplot(gs[1, 0])
    highlight_genes = ['Tfeb', 'Tfe3', 'Trem2']
    highlight_colors = ['#E63946', '#2A9D8F', '#F4A261']

    x = np.arange(len(STATE_ORDER))
    width = 0.25
    for i, (gene, color) in enumerate(zip(highlight_genes, highlight_colors)):
        vals = [1 if gene_state[gene].get(s, False) else 0 for s in STATE_ORDER]
        ax_c.bar(x + i * width, vals, width, label=gene, color=color, edgecolor='black', linewidth=0.5)

    ax_c.set_xticks(x + width)
    ax_c.set_xticklabels(STATE_LABELS, fontsize=9)
    ax_c.set_yticks([0, 1])
    ax_c.set_yticklabels(['Closed', 'Open'], fontsize=10)
    ax_c.set_title('C  Key targets: accessibility by state', fontsize=12, fontweight='bold', loc='left')
    ax_c.legend(fontsize=10)
    ax_c.spines['top'].set_visible(False)
    ax_c.spines['right'].set_visible(False)

    # Panel D: PAM availability for TFEB and TFE3 (to show PAMs exist but chromatin blocks TFEB)
    ax_d = fig.add_subplot(gs[1, 1])
    cas_order_short = ['HEAL_Un1Cas12f1', 'SminiCRa_Un1Cas12f1', 'SpCas9']
    cas_labels_short = ['HEAL', 'SminiCRa', 'SpCas9']

    # Load PAM counts for these genes from the state-invariant homeostatic rows.
    gene_cas_pam = defaultdict(lambda: defaultdict(int))
    with open(FULL_FILE) as f:
        reader = csv.DictReader(f, delimiter='\t')
        for row in reader:
            if row['gene'] in ['Tfeb', 'Tfe3'] and row['state'] == 'homeostatic':
                gene_cas_pam[row['gene']][row['cas']] = int(row['pams_total_passing'])

    x = np.arange(len(cas_order_short))
    width = 0.35
    tfeb_vals = [gene_cas_pam['Tfeb'][c] for c in cas_order_short]
    tfe3_vals = [gene_cas_pam['Tfe3'][c] for c in cas_order_short]

    ax_d.bar(x - width/2, tfeb_vals, width, label='Tfeb (midpoint open in sham only)', color='#E63946',
             edgecolor='black', linewidth=0.5, hatch='//')
    ax_d.bar(x + width/2, tfe3_vals, width, label='Tfe3 (6/6 OPEN)', color='#2A9D8F',
             edgecolor='black', linewidth=0.5)

    for i in range(len(cas_order_short)):
        ax_d.text(x[i] - width/2, tfeb_vals[i] + 0.5, str(tfeb_vals[i]), ha='center', fontsize=9)
        ax_d.text(x[i] + width/2, tfe3_vals[i] + 0.5, str(tfe3_vals[i]), ha='center', fontsize=9)

    ax_d.set_xticks(x)
    ax_d.set_xticklabels(cas_labels_short, fontsize=10)
    ax_d.set_ylabel('Passing PAMs in promoter', fontsize=10)
    ax_d.set_title('D  PAMs available; midpoint criterion limits Tfeb', fontsize=11, fontweight='bold', loc='left')
    ax_d.legend(fontsize=9)
    ax_d.spines['top'].set_visible(False)
    ax_d.spines['right'].set_visible(False)

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
