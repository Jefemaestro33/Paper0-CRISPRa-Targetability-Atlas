#!/usr/bin/env python3
"""
Figure 6: ATAC-seq browser tracks at Tfeb and Tfe3 loci.
Uses pre-generated bigWig files from analysis_results/bigwigs/.
"""
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

try:
    import matplotlib.image as mpimg
    import pyBigWig
    HAS_PYBIGWIG = True
except ImportError:
    HAS_PYBIGWIG = False
    import matplotlib.image as mpimg
    print("WARNING: pyBigWig not installed. Using pre-generated browser-track PNGs from analysis_results/.")

OUTPUT = Path(__file__).resolve().parent.parent / 'output' / 'fig6.pdf'
BWDIR = Path(__file__).resolve().parent.parent.parent / 'analysis_results' / 'bigwigs'
RESULTS_DIR = Path(__file__).resolve().parent.parent.parent / 'analysis_results'

SAMPLES = [
    ('homeostatic', 'Homeostatic\n(Gosselin)', '#457B9D'),
    ('PP', 'Naïve (PP)\n(Zhang X)', '#6B9AC4'),
    ('PL', 'Acute LPS (PL)\n(Zhang X)', '#E76F51'),
    ('LL', 'Tolerized (LL)\n(Zhang X)', '#F4A261'),
    ('sham_WT', 'Post-surgical sham\n(Zhang L)', '#2A9D8F'),
    ('stroke_WT', 'Post-stroke\n(Zhang L)', '#E63946'),
]

GENES = {
    'Tfeb': ('chr17', 48042955, 48108344, 48047955, 48103344),  # chrom, view_start, view_end, gene_start, gene_end
    'Tfe3': ('chrX', 7623799, 7646441, 7628799, 7641441),
}


def plot_gene_panel(fig, axes, gene, chrom, view_start, view_end, gene_start, gene_end):
    """Plot browser tracks for one gene across all samples."""
    all_signals = []
    for sname, label, color in SAMPLES:
        bw_path = BWDIR / f'{sname}.bw'
        if HAS_PYBIGWIG and bw_path.exists():
            bw = pyBigWig.open(str(bw_path))
            signal = np.array(bw.values(chrom, view_start, view_end))
            signal = np.nan_to_num(signal, nan=0.0)
            bw.close()
        else:
            signal = np.zeros(view_end - view_start)
        all_signals.append(signal)

    global_max = max(np.percentile(s, 99.5) for s in all_signals)
    if global_max == 0:
        global_max = 1
    positions = np.linspace(view_start, view_end, len(all_signals[0]))

    tss = gene_start  # + strand
    promo_start = tss - 400
    promo_end = tss - 50

    for i, ((sname, label, color), signal) in enumerate(zip(SAMPLES, all_signals)):
        ax = axes[i]
        ax.fill_between(positions, 0, signal, color=color, alpha=0.85, linewidth=0)
        ax.set_ylim(0, global_max * 1.15)
        ax.set_ylabel(label, fontsize=6.5, rotation=0, ha='right', va='center', labelpad=75)
        ax.set_yticks([])
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['left'].set_visible(False)
        ax.axvspan(gene_start, gene_end, alpha=0.06, color='black')
        ax.axvspan(promo_start, promo_end, alpha=0.15, color='red', zorder=0)

        if i == 0:
            ax.set_title(f'{gene} ({chrom})', fontsize=10, fontweight='bold', loc='left')
            ax.text(0.98, 0.85, f'{global_max:.0f} RPKM', transform=ax.transAxes,
                    fontsize=6, ha='right', va='top', color='gray')

    axes[-1].set_xlabel(f'Position (bp)', fontsize=8)
    # Format x-axis with Mbp
    from matplotlib.ticker import FuncFormatter
    axes[-1].xaxis.set_major_formatter(FuncFormatter(lambda x, p: f'{x/1e6:.3f} Mb'))


def main():
    if not HAS_PYBIGWIG:
        fig, axes = plt.subplots(1, 2, figsize=(16, 7))
        for ax, gene in zip(axes, ['Tfeb', 'Tfe3']):
            png = RESULTS_DIR / f'browser_tracks_{gene}.png'
            if not png.exists():
                raise FileNotFoundError(f'Missing fallback browser-track image: {png}')
            ax.imshow(mpimg.imread(png))
            ax.axis('off')
            ax.set_title(gene, fontsize=12, fontweight='bold', loc='left')
        fig.suptitle('ATAC-seq browser tracks at Tfeb and Tfe3 loci', fontsize=13, fontweight='bold')
        plt.tight_layout()
        OUTPUT.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
        plt.close()
        print(f'Saved: {OUTPUT}')
        return

    fig, all_axes = plt.subplots(len(SAMPLES), 2, figsize=(16, 8),
                                  sharex='col', gridspec_kw={'wspace': 0.3})

    for col, (gene, (chrom, vs, ve, gs, ge)) in enumerate(GENES.items()):
        axes_col = all_axes[:, col]
        plot_gene_panel(fig, axes_col, gene, chrom, vs, ve, gs, ge)

    plt.tight_layout(rect=[0.1, 0.02, 1, 0.98])
    fig.text(0.01, 0.98, 'A', fontsize=14, fontweight='bold', va='top')
    fig.text(0.51, 0.98, 'B', fontsize=14, fontweight='bold', va='top')

    OUTPUT.parent.mkdir(parents=True, exist_ok=True)
    plt.savefig(OUTPUT, dpi=300, bbox_inches='tight')
    plt.close()
    print(f'Saved: {OUTPUT}')


if __name__ == '__main__':
    main()
