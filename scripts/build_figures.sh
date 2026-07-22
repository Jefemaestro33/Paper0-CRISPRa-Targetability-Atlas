#!/usr/bin/env bash
set -euo pipefail

for script in \
  figures/scripts/fig1_overview.py \
  figures/scripts/fig2_pam_availability.py \
  figures/scripts/fig3_chromatin_bottleneck.py \
  figures/scripts/fig4_state_dynamics.py \
  figures/scripts/fig5_practical_framework.py \
  figures/scripts/fig6_browser_tracks.py \
  figures/scripts/figS1_qc.py \
  figures/scripts/figS2_genomewide.py \
  figures/scripts/figS3_window_comparison.py \
  figures/scripts/figS4_motifs.py \
  figures/scripts/figS5_go_enrichment.py \
  figures/scripts/figS6_cross_validation.py \
  figures/scripts/figS7_permutation.py
do
  python "$script"
done
