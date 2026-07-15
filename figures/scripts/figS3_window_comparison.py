#!/usr/bin/env python3
"""
Supplementary Figure S3: sensitivity to accessibility criterion.

The main atlas uses a stringent midpoint criterion plus PAM-in-peak filtering.
The relaxed comparison uses pre-computed promoter/peak overlap statistics from
the permutation analysis. It is an accessibility sensitivity analysis, not an
experimental CRISPRa efficacy estimate.
"""
import csv
from collections import defaultdict
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent.parent
TABLE_S2 = ROOT / 'supplementary' / 'table_S2_targetability_full.tsv'
PERM_FILE = ROOT / 'analysis_stats' / 'permutation_results.tsv'
GENES_FILE = ROOT / 'supplementary' / 'table_S1_therapeutic_genes.csv'
OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'figS3.pdf'

STATE_ORDER = ['homeostatic', 'PP_naive', 'PL_acute_LPS', 'LL_tolerized', 'sham_WT', 'stroke_WT']
STATE_LABELS = ['Homeostatic', 'Naive (PP)', 'Acute LPS (PL)', 'Tolerized (LL)', 'Sham', 'Stroke']


def load_therapeutic():
    genes = []
    with open(GENES_FILE, newline='') as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            genes.append(row['gene_symbol'])
    return genes


def load_midpoint_targetability(therapeutic):
    therapeutic = set(therapeutic)
    targetable = defaultdict(lambda: defaultdict(bool))
    accessible = defaultdict(lambda: defaultdict(bool))
    with open(TABLE_S2, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            if row['cas'] != 'HEAL_Un1Cas12f1' or row['gene'] not in therapeutic:
                continue
            targetable[row['gene']][row['state']] = row['targetable'] == 'True'
            accessible[row['gene']][row['state']] = row['promoter_accessible'] == 'True'
    return targetable, accessible


def load_relaxed_overlap():
    relaxed = {}
    with open(PERM_FILE, newline='') as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        for row in reader:
            relaxed[row['state']] = float(row['observed_pct'])
    return relaxed


def main():
    therapeutic = load_therapeutic()
    targetable, accessible = load_midpoint_targetability(therapeutic)
    relaxed = load_relaxed_overlap()

    n = len(therapeutic)
    midpoint_pct = []
    midpoint_access_pct = []
    relaxed_pct = []
    for state in STATE_ORDER:
        midpoint_pct.append(sum(targetable[g][state] for g in therapeutic) / n * 100)
        midpoint_access_pct.append(sum(accessible[g][state] for g in therapeutic) / n * 100)
        relaxed_pct.append(relaxed[state])
        print(
            f"{state}: midpoint PAM+peak={midpoint_pct[-1]:.1f}%, "
            f"midpoint accessibility={midpoint_access_pct[-1]:.1f}%, "
            f"relaxed promoter overlap={relaxed_pct[-1]:.1f}%"
        )

    gene_mid_count = {g: sum(targetable[g][s] for s in STATE_ORDER) for g in therapeutic}
    gene_access_count = {g: sum(accessible[g][s] for s in STATE_ORDER) for g in therapeutic}

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.5), gridspec_kw={'width_ratios': [1.3, 1, 1]})

    ax_a = axes[0]
    x = np.arange(len(STATE_ORDER))
    w = 0.34
    ax_a.bar(x - w/2, midpoint_pct, w, label='Stringent midpoint + PAM', color='#264653',
             edgecolor='black', lw=0.8)
    ax_a.bar(x + w/2, relaxed_pct, w, label='Relaxed promoter/peak overlap', color='#A8DADC',
             edgecolor='black', lw=0.8)
    for i, (mid, rel) in enumerate(zip(midpoint_pct, relaxed_pct)):
        ax_a.text(i, max(mid, rel) + 2, f'+{rel-mid:.1f}%', ha='center', fontsize=9,
                  fontweight='bold', color='#2A9D8F')
    ax_a.set_xticks(x)
    ax_a.set_xticklabels(STATE_LABELS, fontsize=8, rotation=20, ha='right')
    ax_a.set_ylabel('% therapeutic genes', fontsize=10)
    ax_a.set_ylim(0, 95)
    ax_a.legend(fontsize=8, loc='upper left')
    ax_a.set_title('A  Accessibility criterion sensitivity', fontsize=11, fontweight='bold', loc='left')
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)

    ax_b = axes[1]
    rng = np.random.default_rng(42)
    mid_vals = np.array([gene_mid_count[g] for g in therapeutic])
    access_vals = np.array([gene_access_count[g] for g in therapeutic])
    ax_b.scatter(
        mid_vals + rng.uniform(-0.12, 0.12, len(mid_vals)),
        access_vals + rng.uniform(-0.12, 0.12, len(access_vals)),
        s=38, alpha=0.75, c='#264653', edgecolors='black', lw=0.3
    )
    ax_b.plot([0, 6], [0, 6], '--', color='gray', alpha=0.5)
    ax_b.set_xlabel('# states targetable\n(midpoint + PAM)', fontsize=10)
    ax_b.set_ylabel('# states accessible\n(midpoint only)', fontsize=10)
    ax_b.set_xlim(-0.5, 6.5)
    ax_b.set_ylim(-0.5, 6.5)
    ax_b.set_xticks(range(7))
    ax_b.set_yticks(range(7))
    ax_b.set_aspect('equal')
    ax_b.set_title('B  PAM filter versus accessibility', fontsize=11, fontweight='bold', loc='left')
    ax_b.spines['top'].set_visible(False)
    ax_b.spines['right'].set_visible(False)

    ax_c = axes[2]
    ax_c.axis('off')
    text = "Criterion comparison (55 therapeutic genes)\n"
    text += "-" * 44 + "\n"
    text += "State          Strict   Relaxed  Difference\n"
    for label, mid, rel in zip(STATE_LABELS, midpoint_pct, relaxed_pct):
        short = label.replace(' (PP)', '').replace(' (PL)', '').replace(' (LL)', '')
        text += f"{short:<13}{mid:>6.1f}%  {rel:>7.1f}%  {rel-mid:>8.1f} pp\n"
    text += "-" * 44 + "\n\n"
    text += "Strict = HEAL PAM+chromatin targetability\n"
    text += "using the midpoint criterion.\n\n"
    text += "Relaxed = any ATAC-seq peak overlap with\n"
    text += "the CRISPRa promoter window. This is an\n"
    text += "accessibility sensitivity analysis, not a\n"
    text += "validated guide-efficacy estimate."
    ax_c.text(0.03, 0.95, text, transform=ax_c.transAxes, fontsize=8.2,
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
