#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from common import ROOT,STATE_LABELS,STATES,save


def main():
    frame=pd.read_csv(ROOT/"analysis_stats/sensitivity_summary.tsv",sep="\t")
    frame=frame[(frame.scope=="therapeutic_panel")&(frame.cas=="Un1Cas12f1_TTTR")]
    variants=frame.analysis_variant.drop_duplicates().tolist(); matrix=np.zeros((len(variants),6)); changed=np.zeros_like(matrix)
    for i,variant in enumerate(variants):
        rows=frame[frame.analysis_variant==variant].set_index("state"); matrix[i]=[100*rows.loc[s,"proportion"] for s in STATES]; changed[i]=[rows.loc[s,"n_changed_vs_primary"] for s in STATES]
    fig,axes=plt.subplots(1,2,figsize=(14,5))
    for ax,data,title,fmt in [(axes[0],matrix,"A  Targetability proportion (%)",".1f"),(axes[1],changed,"B  Gene calls changed versus primary",".0f")]:
        im=ax.imshow(data,aspect="auto",cmap="viridis"); ax.set_xticks(range(6)); ax.set_xticklabels(STATE_LABELS,fontsize=7); ax.set_yticks(range(len(variants))); ax.set_yticklabels([v.replace("_"," ") for v in variants],fontsize=7)
        for i in range(data.shape[0]):
            for j in range(data.shape[1]): ax.text(j,i,format(data[i,j],fmt),ha="center",va="center",fontsize=7,color="white" if data[i,j]>data.max()/2 else "black")
        ax.set_title(title,loc="left",weight="bold"); fig.colorbar(im,ax=ax,shrink=.7)
    fig.tight_layout(); save(fig,"figS4")


if __name__=="__main__": main()
