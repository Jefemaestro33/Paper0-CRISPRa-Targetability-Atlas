#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from common import ROOT, STATE_LABELS, STATES, load_atlas, save


TSS_COLORS = {
    "ensembl_canonical": "#457B9D",
    "appris_principal": "#E9C46A",
    "legacy_most5_basic": "#2A9D8F",
}


def main():
    dynamics=pd.read_csv(ROOT/"supplementary/table_S5_accessibility_dynamics.csv")
    sensitivity=pd.read_csv(ROOT/"analysis_stats/therapeutic_gene_stability.tsv",sep="\t")
    tss=pd.read_csv(ROOT/"reference/tss_selection.tsv",sep="\t")
    atlas=load_atlas(); primary=atlas[atlas.cas=="Un1Cas12f1_TTTR"]
    fig=plt.figure(figsize=(7.4,6.25)); grid=fig.add_gridspec(2,2,wspace=.38,hspace=.42)
    ax=fig.add_subplot(grid[0,0]); order=["all_six_contexts","multi_study_support","single_study_support","no_surveyed_context"]
    counts=dynamics.guide_site_support_pattern.value_counts().reindex(order,fill_value=0)
    labels=["All six","≥2\nstudies","One\nstudy","None"]
    ax.bar(labels,counts,color=["#2A9D8F","#457B9D","#E9C46A","#D8D8D8"],edgecolor="black",lw=.5)
    for i,v in enumerate(counts): ax.text(i,v+.4,str(v),ha="center",weight="bold")
    ax.set_ylabel("Genes in locked panel"); ax.set_title("A  Descriptive support patterns",loc="left",weight="bold",fontsize=9); ax.spines[["top","right"]].set_visible(False)

    ax=fig.add_subplot(grid[0,1]); genes=["Tfeb","Tfe3"]
    values=np.zeros((2,6)); overlaps=np.zeros((2,6))
    for i,gene in enumerate(genes):
        rows=primary[primary.gene==gene].set_index("state")
        values[i]=[rows.loc[s,"targetable"] for s in STATES]; overlaps[i]=[rows.loc[s,"promoter_any_peak_overlap"] for s in STATES]
    for i,gene in enumerate(genes):
        ax.scatter(np.arange(6),np.full(6,1-i),s=180,facecolor=np.where(values[i]>0,"#2A9D8F","#EEE8DE"),edgecolor=np.where(overlaps[i]>0,"#D1495B","#777"),linewidth=2)
    short_states = ["Homeo.", "PBS/PBS", "PBS/LPS", "LPS/LPS", "Sham", "Stroke"]
    ax.set_yticks([1,0]); ax.set_yticklabels(genes,fontstyle="italic"); ax.set_xticks(range(6)); ax.set_xticklabels(short_states,fontsize=5.8,rotation=28,ha="right")
    ax.set_xlim(-.5,5.5); ax.set_ylim(-.7,1.7); ax.set_title("B  Tfeb/Tfe3 guide-site and promoter-peak support",loc="left",weight="bold",fontsize=8.2)

    ax=fig.add_subplot(grid[1,0]); subset=tss[(tss.gene_symbol.isin(genes))]
    for i,gene in enumerate(genes):
        rows=subset[subset.gene_symbol==gene]
        base=int(rows[rows.definition=="ensembl_canonical"].tss.iloc[0])
        for j,definition in enumerate(["ensembl_canonical","appris_principal","legacy_most5_basic"]):
            value=int(rows[rows.definition==definition].tss.iloc[0])-base
            ax.scatter(value,i+(j-1)*.16,s=70,color=TSS_COLORS[definition],edgecolor="#222",linewidth=.4,label=definition if i==0 else None)
        if gene == "Tfe3":
            ax.text(-1800, i+.26, "0 bp (all definitions)", ha="right", fontsize=6.3)
        else:
            ax.text(0, i+.18, "0 bp", ha="center", fontsize=6.3)
            ax.text(-48690, i+.20, "−48,690 bp\n(APPRIS/legacy)", ha="center", fontsize=6.3)
    ax.axvline(0,color="#777",lw=.7); ax.set_ylim(-.35,1.4); ax.set_yticks([0,1]); ax.set_yticklabels(genes,fontstyle="italic"); ax.set_xlabel("TSS offset relative to Ensembl canonical")
    ax.legend(fontsize=6.2, frameon=False); ax.set_title("C  TSS-definition sensitivity",loc="left",weight="bold",fontsize=9); ax.spines[["top","right"]].set_visible(False)

    ax=fig.add_subplot(grid[1,1]); sub=sensitivity[sensitivity.gene.isin(genes)]
    variants=[c for c in sub.columns if c not in {"gene","state","n_variants_targetable","n_variants","unanimous"}]
    matrix=[]; labels=[]
    for gene in genes:
        for state in STATES:
            row=sub[(sub.gene==gene)&(sub.state==state)].iloc[0]
            matrix.append([str(row[v]).lower() == "true" for v in variants]); labels.append(f"{gene} | {state}")
    ax.imshow(matrix,aspect="auto",cmap=plt.matplotlib.colors.ListedColormap(["#EEE8DE","#2A9D8F"]),vmin=0,vmax=1)
    variant_labels = ["Primary\ncanonical", "Matched\nGenrich", "Matched\nMACS3", "APPRIS\nTSS", "Legacy\nTSS"]
    ax.set_xticks(range(len(variants))); ax.set_xticklabels(variant_labels,fontsize=6.2)
    ax.set_yticks(range(len(labels))); ax.set_yticklabels(labels,fontsize=5.8); ax.set_title("D  Calls across TSS/depth/caller variants",loc="left",weight="bold",fontsize=8.2)
    save(fig,"fig4")


if __name__ == "__main__": main()
