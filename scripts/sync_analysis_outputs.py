#!/usr/bin/env python3
"""
Synchronize supplementary tables from the six-state atlas.

The complete targetability matrix (Table S2) is treated as the source of truth
for therapeutic-gene accessibility patterns. This script updates the smaller
supplementary tables so they do not retain obsolete four-state classifications.
"""
import csv
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUPP = ROOT / "supplementary"
STATS = ROOT / "analysis_stats"

TABLE_S2 = SUPP / "table_S2_targetability_full.tsv"
TABLE_S1 = SUPP / "table_S1_therapeutic_genes.csv"
TABLE_S3 = SUPP / "table_S3_sgrna_recommendations.csv"
TABLE_S5 = SUPP / "table_S5_accessibility_dynamics.csv"
TABLE_S7 = SUPP / "table_S7_statistical_tests.csv"
THERAPEUTIC_FILE = TABLE_S1

STATES = [
    "homeostatic",
    "PP_naive",
    "PL_acute_LPS",
    "LL_tolerized",
    "sham_WT",
    "stroke_WT",
]


def bool_text(value):
    return "True" if value else "False"


def load_therapeutic_genes():
    genes = []
    metadata = {}
    with open(THERAPEUTIC_FILE, newline="") as handle:
        delimiter = "\t" if THERAPEUTIC_FILE.suffix == ".tsv" else ","
        reader = csv.DictReader(handle, delimiter=delimiter)
        for row in reader:
            gene = row["gene_symbol"]
            genes.append(gene)
            metadata[gene] = row
    return genes, metadata


def load_table_s2(therapeutic):
    therapeutic = set(therapeutic)
    atlas = defaultdict(lambda: defaultdict(dict))
    with open(TABLE_S2, newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            gene = row["gene"]
            if gene not in therapeutic:
                continue
            cas = row["cas"]
            state = row["state"]
            atlas[gene][cas][state] = {
                "promoter_accessible": row["promoter_accessible"] == "True",
                "targetable": row["targetable"] == "True",
                "pams_total_passing": int(row["pams_total_passing"]),
                "pams_in_peak": int(row["pams_in_peak"]),
            }
    return atlas


def state_accessibility(atlas, gene):
    """Promoter accessibility is Cas-independent; use HEAL rows for consistency."""
    cas_data = atlas[gene]["HEAL_Un1Cas12f1"]
    return {state: cas_data[state]["promoter_accessible"] for state in STATES}


def classify_pattern(accessible):
    open_states = [state for state in STATES if accessible[state]]
    if len(open_states) == 6:
        return "constitutively_open"
    if len(open_states) == 0:
        return "never_accessible"
    if (
        not accessible["homeostatic"]
        and accessible["PP_naive"]
        and accessible["PL_acute_LPS"]
        and accessible["LL_tolerized"]
        and accessible["sham_WT"]
        and accessible["stroke_WT"]
    ):
        return "inflammation_gained"
    if accessible["homeostatic"]:
        return "other_pattern"
    if accessible["sham_WT"] or accessible["stroke_WT"]:
        return "surgical_stroke_context"
    return "other_pattern"


def pattern_note(gene, pattern, accessible):
    labels = {
        "homeostatic": "homeostatic",
        "PP_naive": "naive",
        "PL_acute_LPS": "acute LPS",
        "LL_tolerized": "tolerized",
        "sham_WT": "sham",
        "stroke_WT": "stroke",
    }
    open_states = [labels[state] for state in STATES if accessible[state]]
    if pattern == "constitutively_open":
        return "Promoter accessible in all six surveyed states."
    if pattern == "never_accessible":
        return "No promoter accessibility detected under the atlas midpoint criterion."
    if pattern == "inflammation_gained":
        return "Closed in homeostatic microglia; accessible in activated/injury-associated states."
    if pattern == "surgical_stroke_context":
        return "Accessibility detected in at least one sham or stroke dataset state; this is a predictive state-context label, not evidence of injury-driven accessibility: " + "; ".join(open_states) + "."
    return "Non-canonical pattern: accessible in " + "; ".join(open_states) + "."


def write_table_s1(genes, metadata, classifications):
    fields = [
        "gene_symbol",
        "gene_name",
        "category",
        "priority",
        "justification",
        "accessibility_pattern",
    ]
    with open(TABLE_S1, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for gene in genes:
            row = {field: metadata[gene].get(field, "") for field in fields}
            row["accessibility_pattern"] = classifications[gene]
            writer.writerow(row)


def write_table_s5(genes, metadata, accessibility, classifications):
    fields = [
        "gene_symbol",
        "category",
        "accessibility_pattern",
        "accessible_homeostatic",
        "accessible_PP_naive",
        "accessible_PL_acute_LPS",
        "accessible_LL_tolerized",
        "accessible_sham_WT",
        "accessible_stroke_WT",
        "n_states_accessible",
        "notes",
    ]
    with open(TABLE_S5, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for gene in genes:
            acc = accessibility[gene]
            pattern = classifications[gene]
            writer.writerow({
                "gene_symbol": gene,
                "category": metadata[gene]["category"],
                "accessibility_pattern": pattern,
                "accessible_homeostatic": bool_text(acc["homeostatic"]),
                "accessible_PP_naive": bool_text(acc["PP_naive"]),
                "accessible_PL_acute_LPS": bool_text(acc["PL_acute_LPS"]),
                "accessible_LL_tolerized": bool_text(acc["LL_tolerized"]),
                "accessible_sham_WT": bool_text(acc["sham_WT"]),
                "accessible_stroke_WT": bool_text(acc["stroke_WT"]),
                "n_states_accessible": sum(acc.values()),
                "notes": pattern_note(gene, pattern, acc),
            })


def recommendation_class(targetable_states):
    n_states = sum(targetable_states.values())
    if n_states == 6:
        return "constitutive_atlas_candidate"
    if n_states == 0:
        return "not_targetable_under_atlas"
    if targetable_states["sham_WT"] or targetable_states["stroke_WT"]:
        return "state_conditional_atlas_candidate"
    return "limited_state_atlas_candidate"


def write_table_s3(atlas):
    original_rows = []
    with open(TABLE_S3, newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            original_rows.append(row)

    expanded = []
    seen = set()
    for row in original_rows:
        base = {
            "gene_symbol": row["gene_symbol"],
            "cas_ortholog": row["cas_ortholog"],
            "rank": row["rank"],
            "protospacer_sequence": row["protospacer_sequence"],
            "pam_sequence": row["pam_sequence"],
            "strand": row["strand"],
            "genomic_position": row["genomic_position"],
            "gc_content": row["gc_content"],
            "heuristic_score": row["heuristic_score"],
        }
        key = tuple(base.values())
        if key not in seen:
            seen.add(key)
            expanded.append(base)

        # HEAL and SminiCRa share Un1Cas12f1/TTTR targeting logic. Include HEAL
        # rows explicitly where the original table only listed SminiCRa candidates.
        if row["cas_ortholog"] == "SminiCRa_Un1Cas12f1":
            heal = dict(base)
            heal["cas_ortholog"] = "HEAL_Un1Cas12f1"
            key = tuple(heal.values())
            if key not in seen:
                seen.add(key)
                expanded.append(heal)

    fields = [
        "gene_symbol",
        "cas_ortholog",
        "rank",
        "protospacer_sequence",
        "pam_sequence",
        "strand",
        "genomic_position",
        "gc_content",
        "heuristic_score",
        "atlas_targetable_homeostatic",
        "atlas_targetable_PP_naive",
        "atlas_targetable_PL_acute_LPS",
        "atlas_targetable_LL_tolerized",
        "atlas_targetable_sham_WT",
        "atlas_targetable_stroke_WT",
        "atlas_n_targetable_states",
        "recommendation_class",
        "note",
    ]

    with open(TABLE_S3, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        for row in expanded:
            gene = row["gene_symbol"]
            cas = row["cas_ortholog"]
            if gene not in atlas or cas not in atlas[gene]:
                continue
            targetable = {state: atlas[gene][cas][state]["targetable"] for state in STATES}
            out = dict(row)
            out.update({
                "atlas_targetable_homeostatic": bool_text(targetable["homeostatic"]),
                "atlas_targetable_PP_naive": bool_text(targetable["PP_naive"]),
                "atlas_targetable_PL_acute_LPS": bool_text(targetable["PL_acute_LPS"]),
                "atlas_targetable_LL_tolerized": bool_text(targetable["LL_tolerized"]),
                "atlas_targetable_sham_WT": bool_text(targetable["sham_WT"]),
                "atlas_targetable_stroke_WT": bool_text(targetable["stroke_WT"]),
                "atlas_n_targetable_states": sum(targetable.values()),
                "recommendation_class": recommendation_class(targetable),
                "note": (
                    "Predictive atlas candidate. State flags indicate gene/Cas-level "
                    "PAM+chromatin targetability; guide efficacy and off-target risk "
                    "require experimental/orthogonal validation."
                ),
            })
            writer.writerow(out)


def write_table_s7():
    fields = ["analysis", "cas", "state", "metric", "value", "n", "ci_lower", "ci_upper", "p_value", "z_score"]
    rows = []

    with open(STATS / "bootstrap_cis.tsv", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append({
                "analysis": "bootstrap_ci",
                "cas": row["cas"],
                "state": row["state"],
                "metric": "targetability_proportion",
                "value": row["proportion"],
                "n": f'{row["n_targetable"]}/{row["n_total"]}',
                "ci_lower": row["ci_lower"],
                "ci_upper": row["ci_upper"],
                "p_value": "",
                "z_score": "",
            })

    with open(STATS / "pam_chromatin_loss.tsv", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append({
                "analysis": "pam_chromatin_loss",
                "cas": row["cas"],
                "state": row["state"],
                "metric": "pam_bearing_promoters_not_predicted_targetable",
                "value": (
                    f'{row["pam_bearing_not_targetable"]}/'
                    f'{int(row["both_targetable"]) + int(row["pam_bearing_not_targetable"])}'
                ),
                "n": (
                    f'{row["both_targetable"]} targetable; '
                    f'{row["pam_bearing_not_targetable"]} PAM-bearing not targetable; '
                    f'{row["pam_absent"]} PAM-absent'
                ),
                "ci_lower": "",
                "ci_upper": "",
                "p_value": "",
                "z_score": "",
            })

    with open(STATS / "permutation_results.tsv", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        for row in reader:
            rows.append({
                "analysis": "peak_shuffle_permutation",
                "cas": "",
                "state": row["state"],
                "metric": "therapeutic_promoter_peak_overlap",
                "value": row["observed_pct"],
                "n": row["observed_targetable"],
                "ci_lower": "",
                "ci_upper": "",
                "p_value": row["p_value"],
                "z_score": row["z_score"],
            })

    with open(TABLE_S7, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def main():
    genes, metadata = load_therapeutic_genes()
    atlas = load_table_s2(genes)
    accessibility = {gene: state_accessibility(atlas, gene) for gene in genes}
    classifications = {gene: classify_pattern(accessibility[gene]) for gene in genes}

    write_table_s1(genes, metadata, classifications)
    write_table_s5(genes, metadata, accessibility, classifications)
    write_table_s3(atlas)
    write_table_s7()

    counts = defaultdict(int)
    for pattern in classifications.values():
        counts[pattern] += 1
    print("Updated Tables S1, S3, S5, and S7")
    for pattern in sorted(counts):
        print(f"  {pattern}: {counts[pattern]}")


if __name__ == "__main__":
    main()
