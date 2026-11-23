#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")"
tmp=$(mktemp -d)
trap 'rm -rf -- "$tmp"' EXIT
gcc -O3 -fopenmp ../motif-bracket/cert_build.c -lm -o "$tmp/build"
"$tmp/build" 1 1133 2000 "$tmp/upper.bin"
"$tmp/build" 2 113 200 "$tmp/lower.bin"
python3 ../motif-bracket/verify_python_cross.py --max-tree-size "${1:-14}" --upper "$tmp/upper.bin" --lower "$tmp/lower.bin"
gzip -dc ../banded-upper/recipe.txt.gz > "$tmp/recipe"
g++ -std=c++17 -O2 -fsanitize=undefined -fno-sanitize-recover=all read_weights.cpp -o "$tmp/lookup"
python3 banded_histories.py "$tmp/lookup" "$tmp/recipe" --max-size "${2:-9}"
