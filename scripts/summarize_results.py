#!/usr/bin/env python3
"""Generate manuscript macros and a human-readable audit from final tables."""
from __future__ import annotations

from pathlib import Path
import pandas as pd


ROOT=Path(__file__).resolve().parents[1]
STATES=["homeostatic","PP_control","PL_acute_LPS","LL_tolerized","sham_WT","stroke_WT"]
STATE_MACRO={"homeostatic":"Homeostatic","PP_control":"PPControl","PL_acute_LPS":"PLAcuteLPS","LL_tolerized":"LLTolerized","sham_WT":"Sham","stroke_WT":"Stroke"}
CAS_MACRO={"Un1Cas12f1_TTTR":"Un","SaCas9_NNGRRT":"Sa","SpCas9_NGG":"Sp","CjCas9_NNNVRYM":"Cj","Nme2Cas9_NNNNCC":"Nme"}


def pct(value): return f"{100*value:.1f}"


def main():
    atlas=pd.read_csv(ROOT/"supplementary/table_S2_targetability_full.tsv.gz",sep="\t",low_memory=False)
    for col in ("targetable","is_therapeutic"): atlas[col]=atlas[col].astype(str).str.lower().eq("true")
    panel=atlas[atlas.is_therapeutic]; base=atlas[atlas.state==STATES[0]].drop_duplicates(["gene","cas"])
    dynamics=pd.read_csv(ROOT/"supplementary/table_S5_accessibility_dynamics.csv")
    qc=pd.read_csv(ROOT/"supplementary/table_S4_atac_qc.csv"); runs=qc[qc.level=="run"]
    sensitivity=pd.read_csv(ROOT/"analysis_stats/sensitivity_summary.tsv",sep="\t")
    multiplicity=pd.read_csv(ROOT/"analysis_stats/cas_multiplicity_summary.tsv",sep="\t")
    tss=pd.read_csv(ROOT/"reference/tss_selection.tsv",sep="\t")
    macros={"NGene":atlas.gene.nunique(),"NPanel":panel.gene.nunique(),"NRun":len(runs)}
    lines=["# Final result audit","",f"Genes: {macros['NGene']:,}; locked panel: {macros['NPanel']}; raw runs: {macros['NRun']}",""]
    for cas,label in CAS_MACRO.items():
        rows=base[base.cas==cas]; panel_rows=rows[rows.is_therapeutic]
        genome_n=int((rows.protospacers_total_passing>0).sum()); panel_n=int((panel_rows.protospacers_total_passing>0).sum())
        macros[f"{label}GenomePamN"]=genome_n; macros[f"{label}GenomePamPct"]=pct(genome_n/len(rows))
        macros[f"{label}PanelPamN"]=panel_n; macros[f"{label}PanelPamPct"]=pct(panel_n/len(panel_rows))
        lines.append(f"- {cas}: genome PAM/protospacer {genome_n}/{len(rows)} ({macros[f'{label}GenomePamPct']}%); panel {panel_n}/{len(panel_rows)} ({macros[f'{label}PanelPamPct']}%).")
        for state in STATES:
            state_rows=panel[(panel.cas==cas)&(panel.state==state)]
            n=int(state_rows.targetable.sum()); macros[f"{label}{STATE_MACRO[state]}TargetN"]=n; macros[f"{label}{STATE_MACRO[state]}TargetPct"]=pct(n/len(state_rows))
    def multiplicity_n(scope: str, metric: str, state: str, n_supported_classes: int) -> int:
        rows=multiplicity[
            (multiplicity.scope==scope)
            & (multiplicity.metric==metric)
            & (multiplicity.state==state)
            & (multiplicity.n_supported_classes==n_supported_classes)
        ]
        if len(rows) != 1:
            raise ValueError(f"Expected one multiplicity row for {scope}/{metric}/{state}/{n_supported_classes}; found {len(rows)}")
        return int(rows.iloc[0].genes)
    macros["MultiGenomeSeqAllFiveN"]=multiplicity_n("genome","sequence_candidate","sequence_window",5)
    macros["MultiPanelSeqAllFiveN"]=multiplicity_n("locked_panel","sequence_candidate","sequence_window",5)
    macros["MultiGenomeSeqOneN"]=multiplicity_n("genome","sequence_candidate","sequence_window",1)
    macros["MultiGenomeSeqNoneN"]=multiplicity_n("genome","sequence_candidate","sequence_window",0)
    macros["MultiGenomeAnyAllFiveN"]=multiplicity_n("genome","primary_targetable_any_context","any_context",5)
    macros["MultiPanelAnyAllFiveN"]=multiplicity_n("locked_panel","primary_targetable_any_context","any_context",5)
    macros["MultiGenomeAnyNoneN"]=multiplicity_n("genome","primary_targetable_any_context","any_context",0)
    macros["MultiPanelAnyNoneN"]=multiplicity_n("locked_panel","primary_targetable_any_context","any_context",0)
    lines.extend([
        "",
        "Cas-class multiplicity:",
        f"- Sequence layer: {macros['MultiGenomeSeqAllFiveN']}/{macros['NGene']} genome genes and {macros['MultiPanelSeqAllFiveN']}/{macros['NPanel']} panel genes have passing candidates for all five targeting classes.",
        f"- Sequence layer: {macros['MultiGenomeSeqOneN']} genome genes have candidates for exactly one class and {macros['MultiGenomeSeqNoneN']} have none.",
        f"- Any-context primary layer: {macros['MultiGenomeAnyAllFiveN']}/{macros['NGene']} genome genes and {macros['MultiPanelAnyAllFiveN']}/{macros['NPanel']} panel genes have support for all five classes.",
        f"- Any-context primary layer: {macros['MultiGenomeAnyNoneN']} genome genes and {macros['MultiPanelAnyNoneN']} panel genes have no support for any class.",
    ])
    un_context_values=[]
    class_spreads=[]
    genome_context_values=[]
    for state in STATES:
        state_panel=panel[panel.state==state].groupby("cas").targetable.mean()
        state_genome=atlas[atlas.state==state].groupby("cas").targetable.mean()
        un_context_values.append(float(state_panel.loc["Un1Cas12f1_TTTR"]))
        class_spreads.append(float(state_panel.max()-state_panel.min()))
        genome_context_values.append(float(state_genome.loc["Un1Cas12f1_TTTR"]))
    macros["UnPanelTargetMinPct"]=pct(min(un_context_values)); macros["UnPanelTargetMaxPct"]=pct(max(un_context_values))
    macros["UnPanelAcrossContextRangePP"]=f"{100*(max(un_context_values)-min(un_context_values)):.1f}"
    macros["PanelClassSpreadMinPP"]=f"{100*min(class_spreads):.1f}"; macros["PanelClassSpreadMaxPP"]=f"{100*max(class_spreads):.1f}"
    macros["UnGenomeTargetMinPct"]=pct(min(genome_context_values)); macros["UnGenomeTargetMaxPct"]=pct(max(genome_context_values))
    lines.extend(["","Un1Cas12f1 panel guide-site calls:"])
    for state in STATES:
        lines.append(f"- {state}: {macros[f'Un{STATE_MACRO[state]}TargetN']}/{macros['NPanel']} ({macros[f'Un{STATE_MACRO[state]}TargetPct']}%).")
    for gene in ("Tfeb","Tfe3"):
        macro_gene = "TfeThree" if gene == "Tfe3" else gene
        rows=panel[(panel.gene==gene)&(panel.cas=="Un1Cas12f1_TTTR")].set_index("state")
        macros[f"{macro_gene}UnPamN"]=int(rows.protospacers_total_passing.iloc[0])
        macros[f"{macro_gene}UnContextsN"]=int(rows.targetable.sum())
        lines.append(f"- {gene}: {macros[f'{macro_gene}UnPamN']} passing TTTR candidates; guide-site support in {macros[f'{macro_gene}UnContextsN']}/6 contexts.")
        gene_tss=tss[tss.gene_symbol==gene].set_index("definition").tss
        macros[f"{macro_gene}CanonicalTSS"]=int(gene_tss["ensembl_canonical"])
        macros[f"{macro_gene}LegacyTSS"]=int(gene_tss["legacy_most5_basic"])
        macros[f"{macro_gene}LegacyOffset"]=abs(int(gene_tss["legacy_most5_basic"]-gene_tss["ensembl_canonical"]))
    pattern_names={"all_six_contexts":"AllSixPatternN","multi_study_support":"MultiStudyPatternN","single_study_support":"SingleStudyPatternN","no_surveyed_context":"NoContextPatternN"}
    for pattern,macro in pattern_names.items(): macros[macro]=int((dynamics.guide_site_support_pattern==pattern).sum())
    numeric_frip=pd.to_numeric(runs.frip,errors="coerce"); numeric_tss=pd.to_numeric(runs.tss_enrichment_max,errors="coerce")
    macros["FripMin"]=f"{numeric_frip.min():.2f}"; macros["FripMax"]=f"{numeric_frip.max():.2f}"; macros["TSSEnrichMin"]=f"{numeric_tss.min():.1f}"; macros["TSSEnrichMax"]=f"{numeric_tss.max():.1f}"
    sensitive=sensitivity[(sensitivity.scope=="therapeutic_panel")&(sensitivity.analysis_variant!="primary_reproducible_peaks_canonical_tss")]
    macros["SensitivityMaxChanged"]=int(sensitive.n_changed_vs_primary.max())
    macros["SensitivityMaxDeltaPP"]=f"{100*sensitive.delta_vs_primary.abs().max():.1f}"
    tss_summary=pd.read_csv(ROOT/"reference/tss_definition_summary.tsv",sep="\t")
    for _,row in tss_summary.iterrows():
        label="Appris" if "appris" in row["comparison"] else "Legacy"
        macros[f"Canonical{label}DifferentN"]=int(row["different_tss"])
    macro_path=ROOT/"manuscript/results_macros.tex"; macro_path.parent.mkdir(exist_ok=True)
    with macro_path.open("w") as handle:
        handle.write("% Auto-generated by scripts/summarize_results.py; do not edit manually.\n")
        for key,value in sorted(macros.items()): handle.write(f"\\newcommand{{\\{key}}}{{{value}}}\n")
    audit=ROOT/"analysis_stats/key_results.md"; audit.write_text("\n".join(lines)+"\n")
    print(f"Wrote {macro_path} and {audit}")


if __name__=="__main__": main()
