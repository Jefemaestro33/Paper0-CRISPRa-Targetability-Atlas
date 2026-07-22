#!/usr/bin/env python3
"""Test descriptive functional-category associations within the curated panel.

This replaces the previous Enrichr analysis against a genome-wide background,
which was inappropriate because the 55 genes were preselected for therapeutic
microglial functions.  The valid universe here is the fixed 55-gene panel.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


ROOT = Path(__file__).resolve().parents[1]


def bh_adjust(p_values: list[float]) -> np.ndarray:
    values = np.asarray(p_values, dtype=float)
    order = np.argsort(values)
    ranked = values[order]
    adjusted = np.minimum.accumulate((ranked * len(values) / np.arange(1, len(values) + 1))[::-1])[::-1]
    result = np.empty_like(adjusted)
    result[order] = np.minimum(adjusted, 1.0)
    return result


def main() -> None:
    frame = pd.read_csv(ROOT / "supplementary/table_S5_accessibility_dynamics.csv")
    rows = []
    patterns = sorted(frame["guide_site_support_pattern"].dropna().unique())
    categories = sorted(frame["category"].dropna().unique())
    for pattern in patterns:
        for category in categories:
            a = int(((frame.guide_site_support_pattern == pattern) & (frame.category == category)).sum())
            b = int(((frame.guide_site_support_pattern == pattern) & (frame.category != category)).sum())
            c = int(((frame.guide_site_support_pattern != pattern) & (frame.category == category)).sum())
            d = int(((frame.guide_site_support_pattern != pattern) & (frame.category != category)).sum())
            odds, p_value = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            rows.append({
                "guide_site_support_pattern": pattern, "curated_functional_category": category,
                "in_pattern_and_category": a, "in_pattern_not_category": b,
                "outside_pattern_in_category": c, "outside_pattern_not_category": d,
                "odds_ratio": odds, "p_value": p_value,
                "universe": "fixed_55_gene_therapeutic_panel",
            })
    result = pd.DataFrame(rows)
    result["bh_adjusted_p"] = bh_adjust(result.p_value.tolist())
    output = ROOT / "analysis_stats/within_panel_functional_associations.tsv"
    result.to_csv(output, sep="\t", index=False)
    print(f"Wrote {output}")


if __name__ == "__main__":
    main()
