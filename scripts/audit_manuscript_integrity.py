#!/usr/bin/env python3
"""Audit manuscript cross-references, figure assets, and common QC residues."""
from __future__ import annotations

import argparse
import re
from pathlib import Path


MAIN_FIGURES = list(range(1, 7))
SUPP_FIGURES = list(range(1, 8))
SUPP_TABLES = list(range(1, 8))

FIGURE_SCRIPT_MAP = {
    "fig1": "fig1_overview.py",
    "fig2": "fig2_pam_availability.py",
    "fig3": "fig3_chromatin_bottleneck.py",
    "fig4": "fig4_state_dynamics.py",
    "fig5": "fig5_practical_framework.py",
    "fig6": "fig6_browser_tracks.py",
    "figS1": "figS1_qc.py",
    "figS2": "figS2_genomewide.py",
    "figS3": "figS3_window_comparison.py",
    "figS4": "figS4_sensitivity.py",
    "figS5": "figS5_panel_associations.py",
    "figS6": "figS6_cross_validation.py",
    "figS7": "figS7_permutation.py",
}

RISK_PATTERNS = [
    "AI-assisted",
    "artificial intelligence",
    "PENDING",
    "XXXXXXXX",
    "TO_BE_FILLED",
    "local software failure",
    "useful but unsurprising",
    "literature-linked",
    "actionable insight",
    "atlas midpoint",
    "atlas candidates",
    "atlas language",
    "10.5281/zenodo.21971670",
]

SCRIPT_RISK_PATTERNS = [
    "21,599 genes",
    "13 public ATAC-seq runs",
    "10,000 panels; seed",
    "hard-coded",
]


def line_number(text: str, index: int) -> int:
    return text[:index].count("\n") + 1


def first_refs(tex: str, regex: str) -> dict[int, int]:
    refs: dict[int, int] = {}
    for match in re.finditer(regex, tex):
        number = int(match.group(1))
        refs.setdefault(number, line_number(tex, match.start()))
    return refs


def assert_sequence(name: str, refs: dict[int, int], expected: list[int]) -> list[str]:
    errors: list[str] = []
    missing = [number for number in expected if number not in refs]
    if missing:
        errors.append(f"{name}: missing first citation(s): {missing}")
    present = [number for number in expected if number in refs]
    ordered = sorted(present, key=lambda number: refs[number])
    if ordered != present:
        errors.append(
            f"{name}: first-citation order {ordered} by line "
            f"{[refs[number] for number in ordered]} != expected {present}"
        )
    return errors


def audit_additional_files(repo: Path, tex: str) -> list[str]:
    errors: list[str] = []
    additional = re.findall(r"\\item\[Additional file (\d+)\.\].*?\\path\{([^}]+)\}", tex)
    numbers = [int(number) for number, _ in additional]
    if numbers != list(range(1, 18)):
        errors.append(f"Additional files numbered {numbers}, expected 1..17")
    for _, rel in additional:
        candidates = [
            repo / "supplementary" / rel,
            repo / "figures" / "output" / rel,
        ]
        if not any(path.exists() for path in candidates):
            errors.append(f"Additional file path not found: {rel}")
    return errors


def audit_figure_assets(repo: Path) -> list[str]:
    errors: list[str] = []
    scripts_dir = repo / "figures" / "scripts"
    output_dir = repo / "figures" / "output"
    build_script = (repo / "scripts" / "build_figures.sh").read_text()
    for stem, script in FIGURE_SCRIPT_MAP.items():
        script_path = scripts_dir / script
        if not script_path.exists():
            errors.append(f"Missing figure script: {script}")
        if script not in build_script:
            errors.append(f"Figure script not called by build_figures.sh: {script}")
        for suffix in [".pdf", ".png"]:
            output = output_dir / f"{stem}{suffix}"
            if not output.exists():
                errors.append(f"Missing figure output: {output.relative_to(repo)}")
    actual_scripts = sorted(path.name for path in scripts_dir.glob("fig*.py"))
    expected_scripts = sorted(FIGURE_SCRIPT_MAP.values())
    if actual_scripts != expected_scripts:
        errors.append(f"Figure scripts {actual_scripts} != expected {expected_scripts}")
    return errors


def audit_risk_patterns(repo: Path) -> list[str]:
    errors: list[str] = []
    active_paths = [
        repo / "manuscript" / "main.tex",
        repo / "README.md",
        repo / "CITATION.cff",
        repo / ".zenodo.json",
        repo / "docs" / "RELEASE_NOTES_v2.1.5.md",
    ]
    for path in active_paths:
        if not path.exists():
            continue
        text = path.read_text(errors="replace")
        for pattern in RISK_PATTERNS:
            if pattern in text:
                errors.append(f"Risk pattern '{pattern}' found in {path.relative_to(repo)}")
    for path in (repo / "figures" / "scripts").glob("fig*.py"):
        text = path.read_text(errors="replace")
        for pattern in SCRIPT_RISK_PATTERNS:
            if pattern in text:
                errors.append(f"Figure-script risk pattern '{pattern}' found in {path.relative_to(repo)}")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--repo", type=Path, default=Path("."))
    args = parser.parse_args()
    repo = args.repo.resolve()
    tex = (repo / "manuscript" / "main.tex").read_text()

    errors: list[str] = []
    main_refs = first_refs(tex, r"Figure~(\d+)")
    supp_refs = first_refs(tex, r"Supplementary Figure~S(\d+)")
    table_refs = first_refs(tex, r"Table~S(\d+)")

    errors += assert_sequence("Main figures", main_refs, MAIN_FIGURES)
    errors += assert_sequence("Supplementary figures", supp_refs, SUPP_FIGURES)
    missing_tables = [number for number in SUPP_TABLES if number not in table_refs]
    if missing_tables:
        errors.append(f"Supplementary tables missing text citation(s): {missing_tables}")

    for number in MAIN_FIGURES:
        if f"\\textbf{{Figure {number}." not in tex:
            errors.append(f"Missing legend for Figure {number}")
    for number in SUPP_FIGURES:
        if f"\\textbf{{Supplementary Figure S{number}." not in tex:
            errors.append(f"Missing legend for Supplementary Figure S{number}")

    errors += audit_additional_files(repo, tex)
    errors += audit_figure_assets(repo)
    errors += audit_risk_patterns(repo)

    print(f"Main figure first citations: {main_refs}")
    print(f"Supplementary figure first citations: {supp_refs}")
    print(f"Supplementary table first citations: {table_refs}")
    print(f"Figure scripts and outputs checked: {len(FIGURE_SCRIPT_MAP)}")
    if errors:
        print("Manuscript integrity audit FAILED:")
        for error in errors:
            print(f"ERROR: {error}")
        raise SystemExit(1)
    print("Manuscript integrity audit PASS")


if __name__ == "__main__":
    main()
