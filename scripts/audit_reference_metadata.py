#!/usr/bin/env python3
"""Audit BibTeX DOI entries against Crossref metadata.

The audit is intentionally strict about explicitly listed authors, because
incorrect author metadata was a pre-submission QC failure mode. Middle-name
initials are treated as compatible with full middle names, but first names and
family names must match after TeX/accent/punctuation normalization.
"""
from __future__ import annotations

import argparse
import html
import json
import re
import time
import unicodedata
import urllib.parse
import urllib.request
from dataclasses import dataclass
from pathlib import Path


YEAR_EXCEPTIONS = {
    # NAR database issue is cited as 2021; Crossref reports online publication
    # in 2020.
    "frankish2021gencode": ("2021", "2020"),
}


@dataclass
class BibAuthor:
    given: str
    family: str
    group: bool = False


def tex_to_text(value: str) -> str:
    """Best-effort conversion of the simple TeX markup used in references.bib."""
    value = html.unescape(value)
    value = re.sub(r"<[^>]+>", "", value)
    value = value.replace("\\&", "&")
    value = re.sub(r"\\[\"'`^~=.uvHtcdbk]\{?([A-Za-z])\}?", r"\1", value)
    value = re.sub(r"\{\\[\"'`^~=.uvHtcdbk]\s*([A-Za-z])\}", r"\1", value)
    value = re.sub(r"\\[A-Za-z]+", "", value)
    return value.replace("{", "").replace("}", "")


def norm(value: str) -> str:
    value = tex_to_text(value)
    value = unicodedata.normalize("NFKD", value).encode("ascii", "ignore").decode()
    value = value.replace(".", "").replace("'", "").replace("’", "")
    value = re.sub(r"[^A-Za-z0-9]+", " ", value)
    return re.sub(r"\s+", " ", value).strip().lower()


def name_tokens(value: str) -> list[str]:
    tokens: list[str] = []
    for token in re.findall(r"[^\W\d_]+", tex_to_text(value), re.UNICODE):
        if token.isupper() and 1 < len(token) <= 4:
            tokens.extend(token)
        else:
            tokens.append(token)
    return [norm(token) for token in tokens if norm(token)]


def initials_compatible(bib_given: str, crossref_given: str) -> bool:
    left = name_tokens(bib_given)
    right = name_tokens(crossref_given)
    if not left or not right:
        return left == right
    for left_token, right_token in zip(left, right):
        if left_token == right_token:
            continue
        if len(left_token) == 1 and left_token == right_token[:1]:
            continue
        if len(right_token) == 1 and right_token == left_token[:1]:
            continue
        return False
    longer = left if len(left) > len(right) else right
    shorter = right if len(left) > len(right) else left
    if len(longer) > len(shorter):
        return all(len(token) == 1 for token in longer[len(shorter):])
    return True


def family_compatible(bib_family: str, crossref_family: str) -> bool:
    left = name_tokens(bib_family)
    right = name_tokens(crossref_family)
    if left == right:
        return True
    # Crossref occasionally places the first Spanish surname in `given` and the
    # final surname in `family`; treat a suffix family match as compatible.
    return bool(left and right and (left[-len(right):] == right or right[-len(left):] == left))


def text_compatible(left: str, right: str) -> bool:
    norm_left = norm(left)
    norm_right = norm(right)
    return norm_left == norm_right or norm_left.replace(" ", "") == norm_right.replace(" ", "")


def full_name_compatible(bib_author: BibAuthor, cr_author: BibAuthor) -> bool:
    return (
        name_tokens(bib_author.given) + name_tokens(bib_author.family)
        == name_tokens(cr_author.given) + name_tokens(cr_author.family)
    )


def parse_bib_entries(path: Path) -> dict[str, str]:
    text = path.read_text()
    entries: dict[str, str] = {}
    for match in re.finditer(r"@(\w+)\{([^,]+),", text):
        start = match.start()
        key = match.group(2)
        end = text.find("\n}\n", start)
        if end == -1:
            end = len(text)
        entries[key] = text[start : end + 3]
    return entries


def field(entry: str, name: str) -> str | None:
    match = re.search(rf"\b{name}=\{{(.+?)\}}\s*,?\n", entry, re.S)
    if not match:
        return None
    return " ".join(match.group(1).split())


def parse_bib_authors(author_field: str | None) -> tuple[list[BibAuthor], bool]:
    authors: list[BibAuthor] = []
    if not author_field:
        return authors, False
    has_others = False
    for raw in [part.strip() for part in author_field.split(" and ")]:
        if raw.lower() == "others":
            has_others = True
            break
        if raw.startswith("{") and raw.endswith("}"):
            authors.append(BibAuthor("", raw.strip("{}"), group=True))
        elif "," in raw:
            family, given = [part.strip() for part in raw.split(",", 1)]
            authors.append(BibAuthor(given, family))
        else:
            parts = raw.split()
            authors.append(BibAuthor(" ".join(parts[:-1]), parts[-1] if parts else raw))
    return authors, has_others


def crossref_author(author: dict) -> BibAuthor:
    if "name" in author:
        return BibAuthor("", author["name"], group=True)
    return BibAuthor(author.get("given", ""), author.get("family", ""))


def fetch_crossref(doi: str, mailto: str) -> dict:
    params = urllib.parse.urlencode({"mailto": mailto})
    url = f"https://api.crossref.org/works/{urllib.parse.quote(doi)}?{params}"
    request = urllib.request.Request(
        url,
        headers={"User-Agent": f"Paper0-reference-audit/2.1.4 (mailto:{mailto})"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode())["message"]


def compare_entry(key: str, entry: str, crossref: dict) -> list[str]:
    failures: list[str] = []
    warnings: list[str] = []

    bib_authors, has_others = parse_bib_authors(field(entry, "author"))
    crossref_authors = [crossref_author(author) for author in crossref.get("author", [])]
    if has_others:
        compare_authors = zip(bib_authors, crossref_authors[: len(bib_authors)])
        if len(crossref_authors) < len(bib_authors):
            failures.append(f"{key}: Crossref has fewer authors than explicit BibTeX prefix")
    else:
        compare_authors = zip(bib_authors, crossref_authors)
        if len(bib_authors) != len(crossref_authors):
            failures.append(
                f"{key}: author count {len(bib_authors)} != Crossref {len(crossref_authors)}"
            )

    for index, (bib_author, cr_author) in enumerate(compare_authors, start=1):
        if bib_author.group or cr_author.group:
            if norm(bib_author.family) != norm(cr_author.family):
                failures.append(
                    f"{key}: author {index} group '{bib_author.family}' != '{cr_author.family}'"
                )
            continue
        if not family_compatible(bib_author.family, cr_author.family):
            failures.append(
                f"{key}: author {index} family '{bib_author.family}' != '{cr_author.family}'"
            )
        if not initials_compatible(bib_author.given, cr_author.given) and not full_name_compatible(bib_author, cr_author):
            failures.append(
                f"{key}: author {index} given '{bib_author.given}' != '{cr_author.given}'"
            )

    title = field(entry, "title") or ""
    cr_title = (crossref.get("title") or [""])[0]
    if not text_compatible(title, cr_title):
        failures.append(f"{key}: title mismatch")

    journal = field(entry, "journal") or ""
    cr_journal = (crossref.get("container-title") or [""])[0]
    if journal and cr_journal and norm(journal) != norm(cr_journal):
        failures.append(f"{key}: journal mismatch")

    for bib_name, cr_name in [("volume", "volume"), ("pages", "page")]:
        bib_value = field(entry, bib_name) or ""
        cr_value = crossref.get(cr_name) or crossref.get("article-number") or ""
        if bib_value and cr_value and norm(bib_value) != norm(cr_value):
            failures.append(f"{key}: {bib_name} '{bib_value}' != '{cr_value}'")

    bib_year = field(entry, "year") or ""
    cr_year = str((crossref.get("published", {}).get("date-parts") or [[""]])[0][0])
    if (bib_year, cr_year) != YEAR_EXCEPTIONS.get(key, (bib_year, cr_year)):
        if bib_year and cr_year and bib_year != cr_year:
            failures.append(f"{key}: year '{bib_year}' != Crossref '{cr_year}'")

    bib_issue = field(entry, "number")
    cr_issue = crossref.get("issue")
    if cr_issue and not bib_issue:
        warnings.append(f"{key}: BibTeX omits Crossref issue {cr_issue}")

    return failures + [f"WARNING: {warning}" for warning in warnings]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("bib", type=Path, nargs="?", default=Path("manuscript/references.bib"))
    parser.add_argument("--mailto", default="0244552@up.edu.mx")
    args = parser.parse_args()

    entries = parse_bib_entries(args.bib)
    failures: list[str] = []
    warnings: list[str] = []
    doi_count = 0
    no_doi: list[str] = []
    for key, entry in entries.items():
        doi = field(entry, "doi")
        if not doi:
            no_doi.append(key)
            continue
        doi_count += 1
        crossref = fetch_crossref(doi, args.mailto)
        messages = compare_entry(key, entry, crossref)
        for message in messages:
            if message.startswith("WARNING:"):
                warnings.append(message)
            else:
                failures.append(message)
        time.sleep(0.25)

    print(f"BibTeX entries: {len(entries)}")
    print(f"DOI entries audited against Crossref: {doi_count}")
    print(f"Entries without DOI requiring manual provenance check: {', '.join(no_doi) or 'none'}")
    for warning in warnings:
        print(warning)
    if failures:
        print("Reference metadata audit FAILED:")
        for failure in failures:
            print(f"ERROR: {failure}")
        raise SystemExit(1)
    print("Reference metadata audit PASS")


if __name__ == "__main__":
    main()
