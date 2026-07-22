#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from common import ROOT, save


def main():
    frame=pd.read_csv(ROOT/"supplementary/table_S4_atac_qc.csv"); runs=frame[frame.level=="run"].copy(); labels=runs.run_accession.tolist(); x=np.arange(len(runs))
    fig,axes=plt.subplots(2,3,figsize=(16,9))
    specs=[
        ("bowtie2_overall_alignment_rate_pct","Alignment rate (%)",None),
        ("picard_percent_duplication","Picard duplication fraction",None),
        ("frip","FRiP",0.3),
        ("tss_enrichment_max","TSS enrichment (maximum)",6),
        ("usable_fragments_or_reads","Usable fragments or reads",None),
        ("median_fragment_size","Median fragment size (PE only)",None),
    ]
    for ax,(column,label,threshold),letter in zip(axes.flat,specs,"ABCDEF"):
        values=pd.to_numeric(runs[column],errors="coerce")
        ax.bar(x,values,color=["#457B9D" if r=="biological" else "#E9C46A" for r in runs.replicate_type],edgecolor="black",lw=.3)
        if threshold is not None: ax.axhline(threshold,color="#D1495B",ls="--",lw=1,label=f"reference {threshold}"); ax.legend(fontsize=7)
        ax.set_xticks(x); ax.set_xticklabels(labels,rotation=70,ha="right",fontsize=6); ax.set_ylabel(label); ax.set_title(f"{letter}  {label}",loc="left",weight="bold"); ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("Per-run QC; yellow bars are technical sequencing runs, not biological replicates",fontsize=11,weight="bold")
    fig.tight_layout(); save(fig,"figS1")


if __name__=="__main__": main()
