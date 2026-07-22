#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch
import pandas as pd

from common import CANDIDATES, STATES, load_atlas, save


def main():
    candidates=pd.read_csv(CANDIDATES)
    rows=candidates[(candidates.gene_symbol=="Tfe3")&(candidates.nuclease_pam_class=="Un1Cas12f1_TTTR")].sort_values("rank").head(5)
    fig,axes=plt.subplots(1,2,figsize=(7.4,4.55),gridspec_kw={"width_ratios":[.9,1.45]})
    ax=axes[0]; ax.axis("off")
    steps=[("1","Choose delivery architecture","Size and payload constraints"),("2","Choose expressed TSS","Canonical + cell-state sensitivity"),("3","Require guide-site support","Complete spacer+PAM in peak"),("4","Inspect robustness","Replicates, depth, caller, context"),("5","Validate experimentally","Activation, function and safety")]
    for i,(number,title,detail) in enumerate(steps):
        y=.9-i*.18; patch=FancyBboxPatch((.08,y-.055),.84,.11,boxstyle="round,pad=.02",facecolor="#EAF2F1",edgecolor="#264653")
        ax.add_patch(patch); ax.text(.13,y,number,ha="center",va="center",fontsize=10,weight="bold",color="#D1495B"); ax.text(.20,y+.018,title,va="center",weight="bold",fontsize=7.3); ax.text(.20,y-.025,detail,va="center",fontsize=6.2,color="#555")
    ax.set_title("A  Evidence-aware prioritization",loc="left",weight="bold",fontsize=9)

    ax=axes[1]; y=range(len(rows)); signals=[]
    for _,row in rows.iterrows():
        values=[
            float(row[f"peak_signal_{state}"])
            for state in STATES
            if str(row[f"guide_fully_in_peak_{state}"]).lower()=="true"
            and pd.notna(row[f"peak_signal_{state}"])
        ]
        signals.append(min([v for v in values if v>0],default=0))
    colors=["#2A9D8F" if n==6 else "#E9C46A" for n in rows.n_primary_states]
    ax.barh(list(y),signals,color=colors,edgecolor="black",lw=.5)
    ax.set_yticks(list(y)); ax.set_yticklabels([f"Rank {int(value)}" for value in rows["rank"]],fontsize=7); ax.invert_yaxis(); ax.set_xlabel("Minimum supporting peak signal across positive contexts")
    for index, (_, row) in enumerate(rows.iterrows()):
        off = row.get("pam_valid_offtargets_total_le3mm", "")
        off_label = int(float(off)) if str(off) not in {"", "nan"} else "NA"
        label = f"{row['protospacer_sequence']}  |  {int(row['n_primary_states'])}/6  |  off≤3: {off_label}"
        ax.text(.02, index, label, transform=ax.get_yaxis_transform(), ha="left", va="center", fontsize=6.1,
                color="white" if colors[index] == "#2A9D8F" else "#333", weight="bold")
    ax.set_xlabel("Minimum supporting peak signal across positive contexts", fontsize=7.2, labelpad=4)
    ax.set_title("B  Tfe3 Un1Cas12f1 candidate protospacers",loc="left",weight="bold",fontsize=9); ax.spines[["top","right"]].set_visible(False)
    fig.text(.47,.035,"Predictive ranking; PAM-aware Bowtie screening does not replace purpose-built\noff-target tools or experimental safety assays.",fontsize=5.8,ha="left")
    fig.subplots_adjust(left=.055, right=.985, top=.90, bottom=.22, wspace=.32)
    save(fig,"fig5")


if __name__ == "__main__": main()
