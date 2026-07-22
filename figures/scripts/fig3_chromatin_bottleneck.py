#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from common import CAS_LABELS, CAS_ORDER, COLORS, STATE_LABELS, STATES, load_atlas, load_panel, save


PRIMARY="Un1Cas12f1_TTTR"


def main():
    atlas, panel = load_atlas(), load_panel()
    subset=atlas[atlas.gene.isin(panel.gene_symbol)]
    one=subset[subset.cas==PRIMARY]
    pam=100*(one[one.state==STATES[0]].protospacers_total_passing>0).mean()
    target=one.groupby("state").targetable.mean().mul(100)

    fig=plt.figure(figsize=(7.4,8.2)); grid=fig.add_gridspec(2,2,height_ratios=[.72,2.15],wspace=.34,hspace=.28)
    ax=fig.add_subplot(grid[0,0]); x=np.arange(6); width=.36
    ax.bar(x-width/2,[pam]*6,width,label="Passing TTTR protospacer",color="#A8DADC",edgecolor="black",lw=.5)
    ax.bar(x+width/2,[target[s] for s in STATES],width,label="Complete guide site in primary peak",color=COLORS[PRIMARY],edgecolor="black",lw=.5)
    for i,state in enumerate(STATES): ax.text(i+width/2,target[state]+1,f"{target[state]:.1f}",ha="center",fontsize=8)
    short_states = ["Homeo.", "PBS/PBS", "PBS/LPS", "LPS/LPS", "Sham", "Stroke"]
    ax.set_xticks(x); ax.set_xticklabels(short_states,fontsize=6.2,rotation=25,ha="right"); ax.set_ylabel("Therapeutic panel (%)"); ax.set_ylim(0,105)
    ax.legend(fontsize=5.8, loc="upper left"); ax.set_title("A  Sequence versus guide-site support",loc="left",weight="bold",fontsize=9); ax.spines[["top","right"]].set_visible(False)

    ax=fig.add_subplot(grid[0,1])
    for cas in CAS_ORDER:
        values=[100*subset[(subset.cas==cas)&(subset.state==state)].targetable.mean() for state in STATES]
        ax.plot(range(6),values,marker="o",lw=1.8,color=COLORS[cas],label=CAS_LABELS[cas].replace("\n"," "))
    ax.set_xticks(range(6)); ax.set_xticklabels(short_states,fontsize=6.2,rotation=25,ha="right"); ax.set_ylabel("Panel guide-site targetability (%)")
    ax.legend(fontsize=5.2,ncol=1,loc="upper left"); ax.set_title("B  Five nuclease/PAM classes",loc="left",weight="bold",fontsize=9); ax.spines[["top","right"]].set_visible(False)

    ax=fig.add_subplot(grid[1,:])
    matrix=one.pivot(index="gene",columns="state",values="targetable").reindex(index=panel.gene_symbol,columns=STATES).fillna(False)
    order=matrix.sum(axis=1).sort_values(ascending=False).index
    matrix=matrix.loc[order]
    ax.imshow(matrix.to_numpy(),aspect="auto",interpolation="nearest",cmap=plt.matplotlib.colors.ListedColormap(["#EEE8DE","#2A9D8F"]),vmin=0,vmax=1)
    ax.set_xticks(range(6)); ax.set_xticklabels(STATE_LABELS,fontsize=7); ax.set_yticks(range(len(matrix))); ax.set_yticklabels(matrix.index,fontsize=6.2)
    for gene in ("Tfeb","Tfe3"):
        if gene in matrix.index: ax.get_yticklabels()[list(matrix.index).index(gene)].set_color("#D1495B"); ax.get_yticklabels()[list(matrix.index).index(gene)].set_weight("bold")
    ax.set_title("C  Guide-specific Un1Cas12f1 support for all 55 locked genes",loc="left",weight="bold",fontsize=9)
    save(fig,"fig3")


if __name__ == "__main__": main()
