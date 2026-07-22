#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from common import ROOT,STATE_LABELS,STATES,save


def main():
    concord=pd.read_csv(ROOT/"workflow/results/sensitivity/peak_caller_concordance.tsv",sep="\t").set_index("condition").reindex(STATES)
    depth=pd.read_csv(ROOT/"workflow/results/matched_depth/matched_depth_summary.tsv",sep="\t").set_index("condition").reindex(STATES)
    fig,axes=plt.subplots(1,2,figsize=(13,5)); ax=axes[0]; x=np.arange(6)
    ax.bar(x-.18,concord.jaccard,.36,label="Jaccard of peak-overlapped promoters",color="#457B9D"); ax.bar(x+.18,concord.binary_concordance,.36,label="Binary promoter concordance",color="#E9C46A"); ax.set_xticks(x); ax.set_xticklabels(STATE_LABELS,fontsize=7); ax.set_ylim(0,1); ax.legend(fontsize=7); ax.set_title("A  Genrich versus MACS3 at matched depth",loc="left",weight="bold")
    ax=axes[1]; ax.bar(x,depth.sampling_fraction,color="#2A9D8F",edgecolor="black",lw=.3); ax.set_xticks(x); ax.set_xticklabels(STATE_LABELS,fontsize=7); ax.set_ylim(0,1.05); ax.set_ylabel("Deterministic sampling fraction"); ax.set_title("B  Depth equalization",loc="left",weight="bold")
    for ax in axes: ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); save(fig,"figS6")


if __name__=="__main__": main()
