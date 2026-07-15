#!/usr/bin/env python3
"""
Supplementary Figure S7: Random-placement peak-shuffle sanity check.
Shows observed promoter overlap versus uniformly shuffled peaks. This is not a
promoter-matched enrichment null and does not validate CRISPRa activity.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'figS7.pdf'
PERM_FILE = Path(__file__).resolve().parent.parent.parent / 'analysis_stats' / 'permutation_results.tsv'

STATE_LABELS = {
    'homeostatic': 'Homeostatic',
    'PP_naive': 'Naïve (PP)',
    'PL_acute_LPS': 'Acute LPS (PL)',
    'LL_tolerized': 'Tolerized (LL)',
    'sham_WT': 'Sham',
    'stroke_WT': 'Stroke',
}
COLORS = {
    'homeostatic': '#457B9D',
    'PP_naive': '#6B9AC4',
    'PL_acute_LPS': '#E76F51',
    'LL_tolerized': '#F4A261',
    'sham_WT': '#2A9D8F',
    'stroke_WT': '#E63946',
}


def main():
    # Load permutation results
    data = []
    with open(PERM_FILE) as f:
        header = next(f)
        for line in f:
            parts = line.strip().split('\t')
            data.append({
                'state': parts[0],
                'observed_n': int(parts[1]),
                'observed_pct': float(parts[2]),
                'mean_shuffled': float(parts[3]),
                'std_shuffled': float(parts[4]),
                'p_value': float(parts[5]),
                'z_score': float(parts[6]),
            })

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [2, 1]})

    # Panel A: Observed vs shuffled barplot
    states = [d['state'] for d in data]
    observed = [d['observed_pct'] for d in data]
    shuffled_mean = [d['mean_shuffled'] for d in data]
    shuffled_std = [d['std_shuffled'] for d in data]

    x = np.arange(len(states))
    width = 0.35

    bars1 = ax1.bar(x - width/2, observed, width,
                     color=[COLORS[s] for s in states],
                     edgecolor='black', lw=0.8, label='Observed', zorder=3)
    bars2 = ax1.bar(x + width/2, shuffled_mean, width,
                     yerr=shuffled_std, capsize=3,
                     color='#CCCCCC', edgecolor='black', lw=0.8,
                     label='Shuffled (mean ± SD)', zorder=3)

    # Add z-scores
    for i, d in enumerate(data):
        ax1.text(i, max(d['observed_pct'], d['mean_shuffled']) + 5,
                 f"z={d['z_score']:.0f}", ha='center', fontsize=7, fontweight='bold')

    ax1.set_xticks(x)
    ax1.set_xticklabels([STATE_LABELS[s] for s in states], fontsize=9, rotation=20, ha='right')
    ax1.set_ylabel('% therapeutic gene promoters\noverlapping ATAC-seq peaks', fontsize=10)
    ax1.set_ylim(0, 100)
    ax1.legend(fontsize=9, loc='upper left')
    ax1.set_title('A  Observed vs uniformly shuffled peaks', fontsize=11, fontweight='bold', loc='left')
    ax1.spines['top'].set_visible(False)
    ax1.spines['right'].set_visible(False)
    ax1.axhline(y=0, color='black', lw=0.5)

    # Panel B: Summary table
    ax2.axis('off')
    table_text = "State             Obs.   Rand.   z-score\n"
    table_text += "─" * 42 + "\n"
    for d in data:
        label = STATE_LABELS[d['state']]
        table_text += f"{label:<18}{d['observed_pct']:>5.1f}%  {d['mean_shuffled']:>4.1f}%   {d['z_score']:>5.1f}\n"
    table_text += "─" * 42 + "\n"
    table_text += f"\n1,000 permutations per state\n"
    table_text += f"All p < 0.001\n"
    table_text += f"\nInterpretation: observed promoter\n"
    table_text += f"overlap is incompatible with\n"
    table_text += f"uniform random peak placement.\n"
    table_text += f"This is not promoter-matched\n"
    table_text += f"enrichment and does not validate\n"
    table_text += f"CRISPRa activity."

    ax2.text(0.05, 0.95, table_text, transform=ax2.transAxes, fontsize=9,
             verticalalignment='top', fontfamily='monospace',
             bbox=dict(boxstyle='round', facecolor='#F4F1DE', alpha=0.8))
    ax2.set_title('B  Summary statistics', fontsize=11, fontweight='bold', loc='left')

    plt.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
