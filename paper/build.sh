#!/usr/bin/env bash
# Build the paper and its supplement from their TeX sources in a clean directory.
# Usage: bash build.sh OUT_DIR   (run next to the two .tex files; OUT_DIR must be empty or absent)
# The two documents link to each other, so each pass builds both. Three rounds settle all references.
set -euo pipefail
src="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
out="${1:?usage: bash build.sh OUT_DIR}"
mkdir -p "$out"
[ -z "$(ls -A "$out")" ] || { echo "build.sh: $out is not empty" >&2; exit 1; }
cp "$src/friend-tree-exponent.tex" "$src/friend-tree-exponent-supplement.tex" "$out/"
cd "$out"
export SOURCE_DATE_EPOCH=1795564800 FORCE_SOURCE_DATE=1 TZ=UTC
for round in 1 2 3; do
  for doc in friend-tree-exponent-supplement friend-tree-exponent; do
    pdflatex -halt-on-error -interaction=nonstopmode "$doc.tex" > "$doc.pass$round.out" || { tail -40 "$doc.log" >&2; exit 1; }
  done
done
for doc in friend-tree-exponent friend-tree-exponent-supplement; do
  if grep -E 'undefined|Rerun to get' "$doc.log" >/dev/null; then grep -E 'undefined|Rerun' "$doc.log" >&2; echo "build.sh: unresolved references in $doc" >&2; exit 1; fi
done
echo "built: $out/friend-tree-exponent.pdf $out/friend-tree-exponent-supplement.pdf"
