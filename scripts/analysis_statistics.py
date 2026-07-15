#!/usr/bin/env python3
"""
Statistical analyses for Paper 0:
1. Random-placement peak shuffle (BEDTools shuffle)
2. Bootstrap 95% CIs on targetability proportions
3. Descriptive PAM-only to PAM+chromatin loss counts
"""
import csv
import subprocess
import shutil
import tempfile
import os
import numpy as np
from collections import defaultdict
from pathlib import Path

# Paths
ROOT = Path(__file__).resolve().parents[1]
TARGETABILITY = ROOT / "supplementary" / "table_S2_targetability_full.tsv"
GENES_FILE = ROOT / "supplementary" / "table_S1_therapeutic_genes.csv"
PEAKS_DIR = ROOT / "data" / "phase2_results" / "peaks"
ZHANG_PEAKS_DIR = ROOT / "data" / "phase2_results" / "peaks_zhang"
PROMOTERS = ROOT / "data" / "processed" / "promoters_crispra_optimal.bed"
GENOME_FILE = ROOT / "analysis_stats" / "mm39.chrom.sizes"
RUN_PERMUTATION = os.environ.get("PAPER0_RUN_PERMUTATION", "0") == "1"

OUTDIR = ROOT / "analysis_stats"
OUTDIR.mkdir(exist_ok=True)

# ============================================================
# Load data
# ============================================================
print("Loading data...")

# Load therapeutic genes
therapeutic = []
with open(GENES_FILE) as f:
    reader = csv.DictReader(f)
    for row in reader:
        therapeutic.append(row['gene_symbol'])
therapeutic_set = set(therapeutic)
print(f"  Therapeutic genes: {len(therapeutic)}")

# Load targetability matrix
# Structure: gene -> cas -> state -> {targetable, promoter_accessible, pams_total_passing}
data = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
states = set()
cas_types = set()

with open(TARGETABILITY) as f:
    reader = csv.DictReader(f, delimiter='\t')
    for row in reader:
        gene = row['gene']
        cas = row['cas']
        state = row['state']
        states.add(state)
        cas_types.add(cas)
        data[gene][cas][state] = {
            'targetable': row['targetable'] == 'True',
            'promoter_accessible': row['promoter_accessible'] == 'True',
            'pams_total': int(row['pams_total_passing']),
            'pams_in_peak': int(row['pams_in_peak']),
            'is_therapeutic': row['is_therapeutic'] == 'True'
        }

states = sorted(states)
cas_types = sorted(cas_types)
print(f"  States: {states}")
print(f"  Cas types: {cas_types}")

# ============================================================
# 1. Bootstrap 95% CIs
# ============================================================
print("\n" + "="*60)
print("1. BOOTSTRAP 95% CIs ON TARGETABILITY PROPORTIONS")
print("="*60)

np.random.seed(42)
N_BOOT = 10000

with open(OUTDIR / "bootstrap_cis.tsv", 'w') as out:
    out.write("cas\tstate\tn_targetable\tn_total\tproportion\tci_lower\tci_upper\n")

    for cas in cas_types:
        for state in states:
            # Get targetability for therapeutic genes
            values = []
            for gene in therapeutic:
                if gene in data and cas in data[gene] and state in data[gene][cas]:
                    values.append(1 if data[gene][cas][state]['targetable'] else 0)

            if not values:
                continue

            values = np.array(values)
            n_target = values.sum()
            n_total = len(values)
            prop = n_target / n_total

            # Bootstrap
            boot_props = []
            for _ in range(N_BOOT):
                sample = np.random.choice(values, size=n_total, replace=True)
                boot_props.append(sample.mean())

            ci_lower = np.percentile(boot_props, 2.5)
            ci_upper = np.percentile(boot_props, 97.5)

            out.write(f"{cas}\t{state}\t{n_target}\t{n_total}\t{prop:.4f}\t{ci_lower:.4f}\t{ci_upper:.4f}\n")

            if cas == 'HEAL_Un1Cas12f1':
                print(f"  {cas} / {state}: {n_target}/{n_total} = {prop*100:.1f}% [95% CI: {ci_lower*100:.1f}–{ci_upper*100:.1f}%]")

# PAM-only proportions (genes with at least 1 passing PAM)
print("\n  PAM-only proportions (HEAL):")
for state in states:
    values = []
    for gene in therapeutic:
        if gene in data and 'HEAL_Un1Cas12f1' in data[gene] and state in data[gene]['HEAL_Un1Cas12f1']:
            has_pam = data[gene]['HEAL_Un1Cas12f1'][state]['pams_total'] > 0
            values.append(1 if has_pam else 0)
    values = np.array(values)
    prop = values.mean()
    boot_props = [np.random.choice(values, size=len(values), replace=True).mean() for _ in range(N_BOOT)]
    ci_lo, ci_hi = np.percentile(boot_props, [2.5, 97.5])
    print(f"    {state}: {values.sum()}/{len(values)} = {prop*100:.1f}% [95% CI: {ci_lo*100:.1f}–{ci_hi*100:.1f}%]")

# ============================================================
# 2. Descriptive PAM-only to PAM+chromatin loss counts
# ============================================================
print("\n" + "="*60)
print("2. PAM-only to PAM+chromatin loss counts")
print("="*60)

with open(OUTDIR / "pam_chromatin_loss.tsv", 'w') as out:
    out.write("cas\tstate\tboth_targetable\tpam_bearing_not_targetable\tchrom_only_impossible\tpam_absent\tlost_fraction_pam_bearing\n")

    for cas in ['HEAL_Un1Cas12f1']:
        for state in states:
            # For each therapeutic gene:
            # PAM-only targetable = has at least 1 passing PAM
            # PAM+chromatin targetable = targetable (PAM in accessible region)
            a = 0  # both targetable
            b = 0  # PAM-only yes, PAM+chromatin no (lost due to chromatin)
            c = 0  # PAM-only no, PAM+chromatin yes (structurally impossible)
            d = 0  # neither

            for gene in therapeutic:
                if gene not in data or cas not in data[gene] or state not in data[gene][cas]:
                    continue
                pam_only = data[gene][cas][state]['pams_total'] > 0
                pam_chrom = data[gene][cas][state]['targetable']

                if pam_only and pam_chrom:
                    a += 1
                elif pam_only and not pam_chrom:
                    b += 1
                elif not pam_only and pam_chrom:
                    c += 1
                else:
                    d += 1

            lost_fraction = b / (a + b) if (a + b) else 0.0
            out.write(f"{cas}\t{state}\t{a}\t{b}\t{c}\t{d}\t{lost_fraction:.4f}\n")
            print(f"  {state}: {b}/{a+b} PAM-bearing promoters not predicted targetable ({lost_fraction*100:.1f}%)")

# ============================================================
# 3. Random-placement peak shuffle
# ============================================================
print("\n" + "="*60)
print("3. RANDOM-PLACEMENT PEAK SHUFFLE")
print("="*60)

if not RUN_PERMUTATION:
    print("  Skipping permutation rebuild by default.")
    print("  Set PAPER0_RUN_PERMUTATION=1 to recompute analysis_stats/permutation_results.tsv.")
    print("\n" + "="*60)
    print("ALL STATISTICS COMPLETE")
    print("="*60)
    print(f"Results in: {OUTDIR}")
    for f in sorted(OUTDIR.iterdir()):
        print(f"  {f.name}")
    raise SystemExit(0)

if shutil.which("bedtools") is None:
    raise SystemExit("bedtools is required for PAPER0_RUN_PERMUTATION=1")

# Find genome chrom.sizes file
genome_file = GENOME_FILE if GENOME_FILE.exists() else None

if genome_file is None:
    # Create from narrowPeak files
    print("  Creating chrom.sizes from peak files...")
    genome_file = OUTDIR / "mm39.chrom.sizes"
    chroms = {}
    for peak_file in PEAKS_DIR.rglob("*.narrowPeak"):
        with open(peak_file) as f:
            for line in f:
                fields = line.strip().split('\t')
                chrom = fields[0]
                end = int(fields[2])
                if chrom not in chroms or end > chroms[chrom]:
                    chroms[chrom] = end
    with open(genome_file, 'w') as f:
        for chrom in sorted(chroms):
            # Add buffer to ensure shuffle works
            f.write(f"{chrom}\t{chroms[chrom] + 10000}\n")

print(f"  Genome file: {genome_file}")

# Map peak files to states
peak_files = {
    'homeostatic': PEAKS_DIR / 'gosselin_2017' / 'homeostatic_peaks.narrowPeak',
    'PP_naive': PEAKS_DIR / 'holtman_wendeln_2022' / 'PP_peaks.narrowPeak',
    'PL_acute_LPS': PEAKS_DIR / 'holtman_wendeln_2022' / 'PL_peaks.narrowPeak',
    'LL_tolerized': PEAKS_DIR / 'holtman_wendeln_2022' / 'LL_peaks.narrowPeak',
    'sham_WT': ZHANG_PEAKS_DIR / 'sham_WT_peaks.narrowPeak',
    'stroke_WT': ZHANG_PEAKS_DIR / 'stroke_WT_peaks.narrowPeak',
}

# Load promoter regions
promoters = {}  # gene -> (chrom, start, end)
if PROMOTERS.exists():
    with open(PROMOTERS) as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 4:
                chrom, start, end, gene = fields[0], int(fields[1]), int(fields[2]), fields[3]
                promoters[gene] = (chrom, start, end)
    print(f"  Promoters loaded: {len(promoters)}")
else:
    print(f"  WARNING: Promoter file not found at {PROMOTERS}")
    print("  Peak-shuffle analysis requires local pipeline intermediates and will be skipped if they are absent.")

# Function to compute targetability from peaks
def compute_targetability_from_peaks(peak_bed, promoters_bed, therapeutic_genes):
    """Intersect peaks with promoters and count targetable therapeutic genes."""
    # Use bedtools intersect
    result = subprocess.run(
        ['bedtools', 'intersect', '-a', str(promoters_bed), '-b', str(peak_bed), '-u'],
        capture_output=True, text=True
    )
    accessible_genes = set()
    for line in result.stdout.strip().split('\n'):
        if line:
            fields = line.split('\t')
            if len(fields) >= 4:
                accessible_genes.add(fields[3])

    targetable = accessible_genes & therapeutic_genes
    return len(targetable), len(therapeutic_genes)

# Run permutations
N_PERM = 1000
print(f"  Running {N_PERM} permutations per state...")

with open(OUTDIR / "permutation_results.tsv", 'w') as out:
    out.write("state\tobserved_targetable\tobserved_pct\tmean_shuffled_pct\tstd_shuffled_pct\tp_value\tz_score\n")

    for state, peak_file in peak_files.items():
        if not peak_file.exists():
            print(f"  WARNING: {peak_file} not found, skipping {state}")
            continue

        # Observed
        obs_n, total = compute_targetability_from_peaks(peak_file, PROMOTERS, therapeutic_set)
        obs_pct = obs_n / total * 100

        # Shuffle peaks N_PERM times
        shuffled_pcts = []
        for i in range(N_PERM):
            with tempfile.NamedTemporaryFile(suffix='.bed', delete=False, mode='w') as tmp:
                tmp_path = tmp.name

            with open(tmp_path, 'w') as f_out:
                subprocess.run(
                    ['bedtools', 'shuffle', '-i', str(peak_file), '-g', str(genome_file), '-chrom', '-noOverlapping'],
                    stdout=f_out, stderr=subprocess.DEVNULL, check=False
                )

            shuf_n, _ = compute_targetability_from_peaks(tmp_path, PROMOTERS, therapeutic_set)
            shuffled_pcts.append(shuf_n / total * 100)
            os.unlink(tmp_path)

            if (i + 1) % 100 == 0:
                print(f"    {state}: {i+1}/{N_PERM} permutations done")

        shuffled_pcts = np.array(shuffled_pcts)
        mean_shuf = shuffled_pcts.mean()
        std_shuf = shuffled_pcts.std()

        # One-tailed random-placement check.
        # H0: peaks are uniformly distributed with respect to therapeutic promoters.
        # H1: observed peaks overlap therapeutic promoters more often than uniformly shuffled peaks.
        # This is not a promoter-matched enrichment null and should not be
        # interpreted as therapeutic-gene-specific enrichment.
        # p-value = Pr(shuffled targetability >= observed | H0)
        #        = fraction of permutations whose shuffled targetability is at least the observed.
        # If observed >> shuffled (peaks are concentrated at therapeutic promoters), p -> 0.
        p_value = (shuffled_pcts >= obs_pct).sum() / N_PERM
        z_score = (obs_pct - mean_shuf) / max(std_shuf, 1e-10)

        out.write(f"{state}\t{obs_n}\t{obs_pct:.1f}\t{mean_shuf:.1f}\t{std_shuf:.1f}\t{p_value:.4f}\t{z_score:.2f}\n")
        print(f"  {state}: observed={obs_pct:.1f}%, shuffled={mean_shuf:.1f}%±{std_shuf:.1f}%, z={z_score:.2f}, p={p_value:.4f}")

print("\n" + "="*60)
print("ALL STATISTICS COMPLETE")
print("="*60)
print(f"Results in: {OUTDIR}")
for f in sorted(OUTDIR.iterdir()):
    print(f"  {f.name}")
