#!/usr/bin/env python3
"""
Supplementary Figure S1: Complete ATAC-seq QC
Panels: A) Alignment rates, B) Read counts, C) Fragment sizes,
        D) FRiP, E) TSS enrichment, F) IDR summary
"""
import csv
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'figS1.pdf'
RESULTS_DIR = Path(__file__).resolve().parent.parent.parent / 'analysis_results'
S4_FILE = Path(__file__).resolve().parent.parent.parent / 'supplementary' / 'table_S4_atac_qc.csv'

SAMPLE_ORDER = ['homeostatic', 'PP', 'PL', 'LL', 'sham_WT', 'stroke_WT']
SAMPLE_LABELS = ['Homeostatic', 'Naïve\n(PP)', 'Acute LPS\n(PL)', 'Tolerized\n(LL)', 'Sham', 'Stroke']
COLORS = ['#457B9D', '#6B9AC4', '#E76F51', '#F4A261', '#2A9D8F', '#E63946']


def main():
    fig, axes = plt.subplots(2, 3, figsize=(16, 10))

    # === Panel A: Alignment rates (sourced from Table S4, single source of truth) ===
    ax = axes[0, 0]
    # Fallback mirrors Table S4 in case the file is unavailable (e.g. moved offline)
    align_map = {'homeostatic': 95.3, 'PP': 97.5, 'PL': 96.6, 'LL': 98.1, 'sham_WT': 93.5, 'stroke_WT': 92.8}
    if S4_FILE.exists():
        with open(S4_FILE) as f:
            for row in csv.DictReader(f):
                cond = row['condition']
                rate = float(row['alignment_rate'].strip().rstrip('%'))
                for s in SAMPLE_ORDER:
                    if cond == s or cond.endswith(s):
                        align_map[s] = rate
    align_rates = [align_map[s] for s in SAMPLE_ORDER]
    bars = ax.bar(range(len(SAMPLE_ORDER)), align_rates, color=COLORS, edgecolor='black', lw=0.5)
    ax.set_ylim(85, 100)
    ax.set_ylabel('Alignment rate (%)', fontsize=9)
    ax.set_xticks(range(len(SAMPLE_ORDER)))
    ax.set_xticklabels(SAMPLE_LABELS, fontsize=7, rotation=30, ha='right')
    ax.axhline(y=92, color='gray', linestyle='--', alpha=0.5)
    ax.set_title('A  Alignment rates', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # === Panel B: Read counts ===
    ax = axes[0, 1]
    # From frip_results.tsv
    frip_data = {}
    frip_file = RESULTS_DIR / 'frip_results.tsv'
    if frip_file.exists():
        with open(frip_file) as f:
            next(f)  # header
            for line in f:
                parts = line.strip().split('\t')
                frip_data[parts[0]] = {'total': int(parts[1]), 'in_peaks': int(parts[2]), 'frip': float(parts[3])}

    read_counts = [frip_data.get(s, {}).get('total', 0) / 1e6 for s in SAMPLE_ORDER]
    bars = ax.bar(range(len(SAMPLE_ORDER)), read_counts, color=COLORS, edgecolor='black', lw=0.5)
    ax.set_ylabel('Usable reads (millions)', fontsize=9)
    ax.set_xticks(range(len(SAMPLE_ORDER)))
    ax.set_xticklabels(SAMPLE_LABELS, fontsize=7, rotation=30, ha='right')
    ax.set_title('B  Usable read counts', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # === Panel C: Fragment size distributions ===
    ax = axes[0, 2]
    pe_samples = ['PP', 'PL', 'LL', 'sham_WT', 'stroke_WT']
    pe_labels = ['Naïve', 'Acute LPS', 'Tolerized', 'Sham', 'Stroke']
    pe_colors = ['#6B9AC4', '#E76F51', '#F4A261', '#2A9D8F', '#E63946']

    for sname, label, color in zip(pe_samples, pe_labels, pe_colors):
        frag_file = RESULTS_DIR / f'fragsize_{sname}.tsv'
        if frag_file.exists():
            sizes = []
            counts = []
            with open(frag_file) as f:
                header = next(f)
                for line in f:
                    if line.startswith('#') or line.startswith('Size'):
                        continue
                    parts = line.strip().split('\t')
                    if len(parts) >= 2:
                        try:
                            s = int(parts[0])
                            c = float(parts[1])
                            if 0 < s <= 800:
                                sizes.append(s)
                                counts.append(c)
                        except (ValueError, IndexError):
                            continue
            if sizes:
                total = sum(counts) if sum(counts) > 0 else 1
                counts_norm = [c / total for c in counts]
                ax.plot(sizes, counts_norm, color=color, alpha=0.7, linewidth=1, label=label)

    ax.set_xlabel('Fragment size (bp)', fontsize=9)
    ax.set_ylabel('Density', fontsize=9)
    ax.set_xlim(0, 800)
    ax.legend(fontsize=6, loc='upper right')
    ax.axvline(x=150, color='gray', linestyle=':', alpha=0.3)
    ax.axvline(x=300, color='gray', linestyle=':', alpha=0.3)
    ax.text(75, ax.get_ylim()[1] * 0.9, 'sub-nuc', fontsize=5, ha='center', color='gray')
    ax.text(225, ax.get_ylim()[1] * 0.9, 'mono-nuc', fontsize=5, ha='center', color='gray')
    ax.set_title('C  Fragment size distributions (PE only)', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # === Panel D: FRiP ===
    ax = axes[1, 0]
    frip_vals = [frip_data.get(s, {}).get('frip', 0) for s in SAMPLE_ORDER]
    bars = ax.bar(range(len(SAMPLE_ORDER)), frip_vals, color=COLORS, edgecolor='black', lw=0.5)
    ax.axhline(y=0.3, color='red', linestyle='--', alpha=0.7, label='ENCODE threshold (0.3)')
    ax.set_ylabel('FRiP', fontsize=9)
    ax.set_xticks(range(len(SAMPLE_ORDER)))
    ax.set_xticklabels(SAMPLE_LABELS, fontsize=7, rotation=30, ha='right')
    ax.set_ylim(0, 0.7)
    ax.legend(fontsize=7)
    ax.set_title('D  Fraction of Reads in Peaks', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # === Panel E: TSS enrichment ===
    ax = axes[1, 1]
    tss_data = {}
    tss_file = RESULTS_DIR / 'tss_enrichment.tsv'
    if tss_file.exists():
        with open(tss_file) as f:
            next(f)
            for line in f:
                parts = line.strip().split('\t')
                if len(parts) >= 2:
                    tss_data[parts[0]] = float(parts[1])

    tss_vals = [tss_data.get(s, 0) for s in SAMPLE_ORDER]
    bars = ax.bar(range(len(SAMPLE_ORDER)), tss_vals, color=COLORS, edgecolor='black', lw=0.5)
    ax.axhline(y=6, color='red', linestyle='--', alpha=0.7, label='ENCODE threshold (6)')
    ax.set_ylabel('TSS enrichment score', fontsize=9)
    ax.set_xticks(range(len(SAMPLE_ORDER)))
    ax.set_xticklabels(SAMPLE_LABELS, fontsize=7, rotation=30, ha='right')
    ax.legend(fontsize=7)
    ax.set_title('E  TSS enrichment', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # === Panel F: IDR summary ===
    ax = axes[1, 2]
    idr_datasets = ['Gosselin', 'Zhang\nsham', 'Zhang\nstroke']
    idr_reprod = [33260, 33040, 31840]
    idr_total = [51593, 51066, 45861]
    idr_pct = [r / t * 100 for r, t in zip(idr_reprod, idr_total)]
    idr_colors = ['#457B9D', '#2A9D8F', '#E63946']

    bars = ax.bar(range(3), idr_pct, color=idr_colors, edgecolor='black', lw=0.5)
    for i, (pct, rep, tot) in enumerate(zip(idr_pct, idr_reprod, idr_total)):
        ax.text(i, pct + 1, f'{pct:.1f}%\n({rep:,}/{tot:,})', ha='center', fontsize=7)
    ax.set_ylabel('Peaks passing IDR < 0.05 (%)', fontsize=9)
    ax.set_xticks(range(3))
    ax.set_xticklabels(idr_datasets, fontsize=8)
    ax.set_ylim(0, 85)
    ax.set_title('F  IDR (biological replicates)', fontsize=11, fontweight='bold', loc='left')
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.text(0.5, -0.15, 'Note: Zhang X (PP, PL, LL) has technical replicates only — IDR not applicable',
            transform=ax.transAxes, fontsize=6, ha='center', style='italic', color='gray')

    plt.tight_layout()
    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {OUTPUT}')


if __name__ == '__main__':
    main()
