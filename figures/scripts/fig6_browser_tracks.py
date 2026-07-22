#!/usr/bin/env python3
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import pyBigWig

from common import CANDIDATES, ROOT, STATE_LABELS, STATES, load_atlas, save


def load_peaks(path,chrom,start,end):
    values=[]
    with open(path) as handle:
        for line in handle:
            fields=line.split("\t");
            if fields[0]==chrom and int(fields[1])<end and int(fields[2])>start: values.append((max(start,int(fields[1])),min(end,int(fields[2]))))
    return values


def main():
    atlas=load_atlas(); candidates=pd.read_csv(CANDIDATES)
    genes=["Tfeb","Tfe3"]; metadata={}
    for gene in genes:
        row=atlas[(atlas.gene==gene)&(atlas.cas=="Un1Cas12f1_TTTR")].iloc[0]
        # Chromosome is recovered from candidate intervals because Table S2 is gene-level.
        candidate=candidates[(candidates.gene_symbol==gene)&(candidates.nuclease_pam_class=="Un1Cas12f1_TTTR")].iloc[0]
        chrom=candidate.target_interval.split(":")[0]; tss=int(row.tss); metadata[gene]=(chrom,tss,row.strand)
    fig,axes=plt.subplots(6,2,figsize=(7.4,8.25),sharex="col",sharey="col")
    for col,gene in enumerate(genes):
        chrom,tss,strand=metadata[gene]; start,end=tss-2500,tss+2500
        pstart,pend=(tss-400,tss-50) if strand=="+" else (tss+51,tss+401)
        gene_candidates=candidates[(candidates.gene_symbol==gene)&(candidates.nuclease_pam_class=="Un1Cas12f1_TTTR")]
        for row,state in enumerate(STATES):
            ax=axes[row,col]; bw_path=ROOT/f"workflow/results/bigwig/{state}.bw"
            if not bw_path.exists(): raise FileNotFoundError(f"Required signal track missing: {bw_path}")
            with pyBigWig.open(str(bw_path)) as bw:
                values=np.nan_to_num(np.asarray(bw.values(chrom,start,end,numpy=True),dtype=float))
            x=np.arange(start,end)-tss; ax.fill_between(x,values,color="#607D8B",alpha=.75,lw=0)
            ax.axvspan(pstart-tss,pend-tss,color="#D1495B",alpha=.16)
            for peak_start,peak_end in load_peaks(ROOT/f"workflow/results/peaks/primary/{state}.narrowPeak",chrom,start,end):
                ax.plot([peak_start-tss,peak_end-tss],[values.max()*1.02]*2,color="#D1495B",lw=2.4,solid_capstyle="butt")
            for _,candidate in gene_candidates.iterrows():
                if str(candidate[f"guide_fully_in_peak_{state}"]).lower()=="true":
                    interval=candidate.target_interval.split(":")[1]; a,b=map(int,interval.split("-")); ax.axvline((a+b)/2-tss,color="#2A9D8F",lw=.7,alpha=.8)
            if col==0: ax.set_ylabel(STATE_LABELS[row].replace("\n"," "),fontsize=6.3)
            if row==0: ax.set_title(f"{gene} | {chrom}:{tss:,} | canonical TSS",fontstyle="italic",weight="bold",fontsize=8)
            ax.set_xlim(-2500,2500); ax.spines[["top","right"]].set_visible(False); ax.tick_params(labelsize=6.3)
        axes[-1,col].set_xticks([-2500,-1250,0,1250,2500], ["−2.5","−1.25","0","1.25","2.5"])
        axes[-1,col].set_xlabel("Position relative to canonical TSS (kb)", fontsize=7.5)
    fig.supylabel("CPM-normalized ATAC signal", x=.01, fontsize=8)
    fig.subplots_adjust(left=.13, right=.985, top=.95, bottom=.07, hspace=.12, wspace=.20)
    save(fig,"fig6")


if __name__ == "__main__": main()
