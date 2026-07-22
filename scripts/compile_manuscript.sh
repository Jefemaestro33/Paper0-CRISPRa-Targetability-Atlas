#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd "$(dirname "$0")/.." && pwd)
cd "$ROOT/manuscript"
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
