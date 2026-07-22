#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from common import ROOT,save


def main():
    frame=pd.read_csv(ROOT/"analysis_stats/within_panel_functional_associations.tsv",sep="\t")
    # Exact Fisher p values remain unchanged. A 0.5 continuity correction is
    # used only to display finite log odds ratios when a contingency cell is 0.
    frame["plot_odds_ratio"]=(frame.in_pattern_and_category+.5)*(frame.outside_pattern_not_category+.5)/((frame.in_pattern_not_category+.5)*(frame.outside_pattern_in_category+.5))
    frame=frame.sort_values("bh_adjusted_p").head(20).sort_values("plot_odds_ratio")
    labels=[f"{p} | {c}" for p,c in zip(frame.guide_site_support_pattern,frame.curated_functional_category)]
    fig,ax=plt.subplots(figsize=(10,7)); values=np.log2(frame.plot_odds_ratio.clip(lower=.05,upper=20)); colors=np.where(frame.bh_adjusted_p<.05,"#D1495B","#A8DADC")
    ax.barh(range(len(frame)),values,color=colors,edgecolor="black",lw=.3); ax.set_yticks(range(len(frame))); ax.set_yticklabels(labels,fontsize=7); ax.axvline(0,color="#555",lw=.8); ax.set_xlabel("log2 odds ratio within the fixed 55-gene panel")
    ax.set_title("Within-panel functional-category associations (Fisher tests; BH across all pattern-category tests)",loc="left",weight="bold"); ax.spines[["top","right"]].set_visible(False)
    ax.text(0,-.09,"Exact Fisher tests use the observed counts; displayed odds ratios use a 0.5 correction. The analysis is descriptive because panel functions were curated a priori.",transform=ax.transAxes,fontsize=8)
    fig.tight_layout(); save(fig,"figS5")


if __name__=="__main__": main()
