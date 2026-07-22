#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from common import CAS_LABELS,CAS_ORDER,COLORS,STATE_LABELS,STATES,load_atlas,load_panel,save


def main():
    atlas,panel=load_atlas(),load_panel(); base=atlas[atlas.state==STATES[0]].drop_duplicates(["gene","cas"])
    fig,axes=plt.subplots(1,3,figsize=(16,5))
    ax=axes[0]
    for cas in CAS_ORDER:
        values=base[base.cas==cas].protospacers_total_passing.to_numpy(); bins=np.arange(0,min(40,np.percentile(values,99))+2)
        ax.hist(values,bins=bins,histtype="step",lw=1.5,label=CAS_LABELS[cas].replace("\n"," "),color=COLORS[cas],density=True)
    ax.set_xlabel("Passing protospacers per promoter (display truncated at 99th percentile)"); ax.set_ylabel("Density"); ax.legend(fontsize=6); ax.set_title("A  Genome-wide candidate counts",loc="left",weight="bold")
    ax=axes[1]; panel_set=set(panel.gene_symbol)
    matrix=[]
    for cas in CAS_ORDER:
        rows=base[base.cas==cas]; matrix.append([100*(rows.protospacers_total_passing>0).mean(),100*(rows[rows.gene.isin(panel_set)].protospacers_total_passing>0).mean()])
    matrix=np.asarray(matrix); x=np.arange(5); ax.bar(x-.18,matrix[:,0],.36,label="Genome",color="#A8DADC"); ax.bar(x+.18,matrix[:,1],.36,label="55-gene panel",color="#D1495B")
    ax.set_xticks(x); ax.set_xticklabels([CAS_LABELS[c].split("\n")[0] for c in CAS_ORDER],rotation=25,ha="right",fontsize=7); ax.set_ylim(0,105); ax.set_ylabel("Coverage (%)"); ax.legend(fontsize=7); ax.set_title("B  Panel versus genome",loc="left",weight="bold")
    ax=axes[2]; one=atlas[atlas.cas=="Un1Cas12f1_TTTR"]
    genome=[100*one[one.state==s].targetable.mean() for s in STATES]; focused=[100*one[(one.state==s)&one.gene.isin(panel_set)].targetable.mean() for s in STATES]
    ax.plot(range(6),genome,"o-",label="Genome",color="#264653"); ax.plot(range(6),focused,"o-",label="55-gene panel",color="#D1495B"); ax.set_xticks(range(6)); ax.set_xticklabels(STATE_LABELS,fontsize=6); ax.set_ylabel("Guide-site targetability (%)"); ax.legend(fontsize=7); ax.set_title("C  Un1Cas12f1 primary calls",loc="left",weight="bold")
    for ax in axes: ax.spines[["top","right"]].set_visible(False)
    fig.tight_layout(); save(fig,"figS2")


if __name__=="__main__": main()
