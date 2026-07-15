#!/usr/bin/env python3
"""
Supplementary Figure S4: HOMER motif enrichment in surgical/stroke-context promoters
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'figS4.pdf'

# HOMER known motif results (from actual analysis)
MOTIFS = [
    {'name': 'T-bet (T-box)', 'pval': 1e-4, 'target_pct': 60.0, 'bg_pct': 0.0, 'qval': 0.0035},
    {'name': 'Hoxc13 (Homeobox)', 'pval': 1e-4, 'target_pct': 60.0, 'bg_pct': 6.3, 'qval': 0.0035},
    {'name': 'bHLH80', 'pval': 1e-4, 'target_pct': 60.0, 'bg_pct': 8.7, 'qval': 0.0035},
    {'name': 'Foxf1 (Forkhead)', 'pval': 1e-3, 'target_pct': 50.0, 'bg_pct': 8.7, 'qval': 0.032},
    {'name': 'PRDM1/Blimp-1 (Zf)', 'pval': 1e-3, 'target_pct': 50.0, 'bg_pct': 6.6, 'qval': 0.032},
    {'name': 'HIC1 (Zf)', 'pval': 1e-3, 'target_pct': 90.0, 'bg_pct': 35.6, 'qval': 0.127},
    {'name': 'Cdx2 (Homeobox)', 'pval': 1e-2, 'target_pct': 40.0, 'bg_pct': 9.0, 'qval': 0.127},
    {'name': 'Hoxd13 (Homeobox)', 'pval': 1e-2, 'target_pct': 40.0, 'bg_pct': 6.3, 'qval': 0.127},
]


def main():
    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={'width_ratios': [1.5, 1]})

    # Panel A: -log10(p-value) barplot with target% annotation
    ax_a = axes[0]
    names = [m['name'] for m in MOTIFS]
    logp = [-np.log10(m['pval']) for m in MOTIFS]
    colors = ['#E63946' if m['qval'] < 0.05 else '#A8DADC' for m in MOTIFS]

    y = range(len(MOTIFS))
    bars = ax_a.barh(y, logp, color=colors, edgecolor='black', linewidth=0.5)

    for i, m in enumerate(MOTIFS):
        ax_a.text(logp[i] + 0.1, i, f"{m['target_pct']:.0f}% vs {m['bg_pct']:.0f}%",
                  va='center', fontsize=8)

    ax_a.set_yticks(y)
    ax_a.set_yticklabels(names, fontsize=9)
    ax_a.set_xlabel('-log10(p-value)', fontsize=10)
    ax_a.set_title('A  Known motifs enriched in\n    surgical/stroke-context promoters (n=10)',
                    fontsize=11, fontweight='bold', loc='left')
    ax_a.axvline(x=-np.log10(0.05), color='gray', linestyle='--', alpha=0.5, label='p=0.05')
    ax_a.invert_yaxis()
    ax_a.spines['top'].set_visible(False)
    ax_a.spines['right'].set_visible(False)

    from matplotlib.patches import Patch
    legend_elements = [Patch(facecolor='#E63946', label='q < 0.05'),
                       Patch(facecolor='#A8DADC', label='q > 0.05')]
    ax_a.legend(handles=legend_elements, fontsize=8, loc='lower right')

    # Panel B: Biological interpretation
    ax_b = axes[1]
    ax_b.axis('off')
    text = (
        "Biologically plausible hits:\n\n"
        "T-bet (Tbx21)\n"
        "  IFN-γ/innate immunity TF\n"
        "  60% of surgical/stroke-context, 0% constitutive\n"
        "  Consistent with IFN-γ signaling post-surgery\n\n"
        "PRDM1 (Blimp-1)\n"
        "  Terminal differentiation of immune cells\n"
        "  50% of surgical/stroke-context, 7% constitutive\n\n"
        "Hoxc13, bHLH80, Foxf1\n"
        "  Less clear microglial function\n"
        "  Reported for completeness\n\n"
        "Limitation: n=10 promoters (3 of 13\n"
        "surgical/stroke-context genes not mapped\n"
        "by HOMER). Low power for rare motifs."
    )
    ax_b.text(0.05, 0.95, text, transform=ax_b.transAxes, fontsize=9,
              verticalalignment='top', fontfamily='monospace',
              bbox=dict(boxstyle='round', facecolor='#F4F1DE', alpha=0.8))
    ax_b.set_title('B  Interpretation', fontsize=11, fontweight='bold', loc='left')

    plt.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f"Saved: {OUTPUT}")


if __name__ == '__main__':
    main()
