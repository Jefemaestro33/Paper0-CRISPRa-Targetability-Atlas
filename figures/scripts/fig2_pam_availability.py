#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import CAS_LABELS, CAS_ORDER, COLORS, STATE_LABELS, STATES, load_atlas, load_panel, save


def main():
    atlas, panel = load_atlas(), load_panel()
    baseline = atlas[atlas.state == STATES[0]].drop_duplicates(["gene", "cas"])
    genome_pam = baseline.assign(has_candidate=baseline.protospacers_total_passing.gt(0)).groupby("cas").has_candidate.mean().mul(100)
    panel_base = baseline[baseline.gene.isin(panel.gene_symbol)]
    panel_pam = panel_base.assign(has_candidate=panel_base.protospacers_total_passing.gt(0)).groupby("cas").has_candidate.mean().mul(100)
    genome_target = atlas.groupby(["cas", "state"]).targetable.mean().mul(100)

    fig = plt.figure(figsize=(7.4, 6.15))
    grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.25], hspace=.40, wspace=.38)
    for axis, series, title in [(fig.add_subplot(grid[0,0]), genome_pam, "A  Genome-wide PAM/protospacer coverage"), (fig.add_subplot(grid[0,1]), panel_pam, "B  Locked 55-gene panel")]:
        values = [series[cas] for cas in CAS_ORDER]
        axis.bar(range(5), values, color=[COLORS[c] for c in CAS_ORDER], edgecolor="black", lw=.6)
        for i, value in enumerate(values):
            axis.text(i, value + 1.0 + (i % 2) * 4.0, f"{value:.1f}%", ha="center", fontsize=7.2, weight="bold")
        axis.set_xticks(range(5)); axis.set_xticklabels([CAS_LABELS[c].split("\n")[0] for c in CAS_ORDER], fontsize=6.5)
        axis.set_ylabel("Promoters with ≥1 candidate (%)", fontsize=7.5); axis.set_ylim(0, 110)
        axis.set_title(title, loc="left", weight="bold", fontsize=9); axis.spines[["top","right"]].set_visible(False)

    ax = fig.add_subplot(grid[1,0])
    width=.15; x=np.arange(len(STATES))
    for j, cas in enumerate(CAS_ORDER):
        vals=[genome_target[(cas,state)] for state in STATES]
        ax.bar(x+(j-2)*width, vals, width, color=COLORS[cas], label=CAS_LABELS[cas].replace("\n"," "), edgecolor="none")
    short_states = ["Homeo.", "PBS/PBS", "PBS/LPS", "LPS/LPS", "Sham", "Stroke"]
    ax.set_xticks(x); ax.set_xticklabels(short_states, fontsize=6.5, rotation=25, ha="right"); ax.set_ylabel("Genome-wide guide-site targetability (%)", fontsize=7.5)
    ax.set_ylim(0, max(genome_target.max() * 1.18, 5))
    ax.set_title("C  Guide-site support by dataset context", loc="left", weight="bold", fontsize=9)
    ax.legend(fontsize=5.2, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.22), frameon=False)
    ax.spines[["top","right"]].set_visible(False)

    ax = fig.add_subplot(grid[1,1])
    data=[]; labels=[]
    for cas in CAS_ORDER:
        data.append(np.log2(baseline.loc[baseline.cas==cas,"protospacers_total_passing"].to_numpy()+1)); labels.append(CAS_LABELS[cas].split("\n")[0])
    ax.boxplot(data, tick_labels=labels, showfliers=False, patch_artist=True, boxprops={"facecolor":"#A8DADC"}, medianprops={"color":"#D1495B","lw":1.5})
    ax.set_ylabel("log2(passing protospacers + 1)", fontsize=7.5); ax.tick_params(axis="x", labelrotation=25, labelsize=6.5)
    ax.set_title("D  Per-promoter candidate abundance", loc="left", weight="bold", fontsize=9); ax.spines[["top","right"]].set_visible(False)
    save(fig, "fig2")


if __name__ == "__main__": main()
