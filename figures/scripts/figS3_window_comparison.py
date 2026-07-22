#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from common import STATE_LABELS,STATES,load_atlas,load_panel,save


def main():
    atlas,panel=load_atlas(),load_panel(); one=atlas[(atlas.cas=="Un1Cas12f1_TTTR")&atlas.gene.isin(panel.gene_symbol)]
    metrics={"Complete guide in peak":one.targetable,"Any guide/peak overlap":one.protospacers_any_peak_overlap>0,"Promoter midpoint in peak":one.promoter_midpoint_accessible,"Any promoter/peak overlap":one.promoter_any_peak_overlap}
    fig,axes=plt.subplots(1,2,figsize=(14,5.5)); ax=axes[0]; x=np.arange(6); width=.19
    colors=["#264653","#457B9D","#E9C46A","#D1495B"]
    for j,(label,values) in enumerate(metrics.items()):
        vals=[100*values[one.state==state].mean() for state in STATES]; ax.bar(x+(j-1.5)*width,vals,width,label=label,color=colors[j])
    ax.set_xticks(x); ax.set_xticklabels(STATE_LABELS,fontsize=7); ax.set_ylabel("Locked panel (%)"); ax.legend(fontsize=7); ax.set_title("A  Operational-definition sensitivity",loc="left",weight="bold"); ax.spines[["top","right"]].set_visible(False)
    ax=axes[1]; primary=one.pivot(index="gene",columns="state",values="targetable").reindex(columns=STATES); any_overlap=one.assign(any=one.protospacers_any_peak_overlap>0).pivot(index="gene",columns="state",values="any").reindex(columns=STATES)
    delta=(any_overlap.astype(int)-primary.astype(int)).sum(axis=1).sort_values(ascending=False).head(20)
    ax.barh(range(len(delta)),delta,color="#E9C46A",edgecolor="black",lw=.3); ax.set_yticks(range(len(delta))); ax.set_yticklabels(delta.index,fontsize=7); ax.invert_yaxis(); ax.set_xlabel("Additional contexts under any-overlap sensitivity"); ax.set_title("B  Genes most sensitive to full containment",loc="left",weight="bold"); ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); save(fig,"figS3")


if __name__=="__main__": main()
