#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from common import ROOT,STATE_LABELS,STATES,load_key_counts,save


def main():
    counts=load_key_counts()
    frame=pd.read_csv(ROOT/"analysis_stats/matched_promoter_null.tsv",sep="\t").set_index("state").reindex(STATES)
    x=np.arange(6); observed=100*frame.observed_panel_proportion; mean=100*frame.matched_null_mean
    lower=100*frame["matched_null_lower_2.5pct"]; upper=100*frame["matched_null_upper_97.5pct"]
    fig,axes=plt.subplots(1,2,figsize=(13,5)); ax=axes[0]
    ax.bar(x-.18,observed,.36,label="Locked panel",color="#D1495B"); ax.bar(x+.18,mean,.36,yerr=np.vstack([mean-lower,upper-mean]),label="Matched non-panel promoters",color="#A8DADC",capsize=3)
    ax.set_xticks(x); ax.set_xticklabels(STATE_LABELS,fontsize=7); ax.set_ylabel("Guide-site targetability (%)"); ax.legend(fontsize=7); ax.set_title("A  Promoter-matched null",loc="left",weight="bold"); ax.spines[["top","right"]].set_visible(False)
    ax=axes[1]; ax.axis("off"); text=f"Matched on:\n• promoter GC decile\n• annotated transcript-count bin\n• passing TTTR-protospacer-count bin\n\n{counts['permutation_resamples']:,} panels; seed {counts['permutation_seed']}\nEmpirical p uses (b+1)/(N+1).\n\nNot matched on expression or mappability;\nthis is a panel-comparison sensitivity,\nnot CRISPRa validation."
    ax.text(.05,.95,text,va="top",fontsize=10,bbox={"boxstyle":"round","facecolor":"#F4F1DE"}); ax.set_title("B  Null specification",loc="left",weight="bold")
    fig.tight_layout(); save(fig,"figS7")


if __name__=="__main__": main()
