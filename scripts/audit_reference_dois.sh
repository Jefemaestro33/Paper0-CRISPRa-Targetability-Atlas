#!/usr/bin/env bash
set -euo pipefail

BIB=${1:-manuscript/references.bib}

while read -r doi; do
  curl -sS --fail --retry 6 --retry-all-errors --retry-delay 2 \
    --user-agent 'Paper0-reference-audit/2.0 (mailto:0244552@up.edu.mx)' \
    "https://api.crossref.org/works/$doi?mailto=0244552@up.edu.mx" \
    | jq -r --arg doi "$doi" '[
        $doi,
        (.message.title[0] // ""),
        (.message["container-title"][0] // ""),
        (.message.volume // ""),
        (.message.issue // ""),
        (.message.page // (.message["article-number"] // "")),
        ((.message.published["date-parts"][0][0] // "") | tostring),
        ([.message.author[]? | ((.given // "") + " " + (.family // ""))] | join("; "))
      ] | @tsv'
  sleep 0.25
done < <(rg -o 'doi=\{[^}]+' "$BIB" | sed 's/doi={//')
