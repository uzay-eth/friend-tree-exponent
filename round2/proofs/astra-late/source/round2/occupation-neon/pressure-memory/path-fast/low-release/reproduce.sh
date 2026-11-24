#!/usr/bin/env bash
# From the archive root: bash reproduce.sh PYTHON LOWER UPPER EMPTY_SCRATCH
set -euo pipefail
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
[[ $# == 4 ]] || { printf 'usage: %s PYTHON LOWER UPPER EMPTY_SCRATCH\n' "$0" >&2; exit 2; }
P=$1 L=$(realpath -e "$2") U=$(realpath -e "$3") OUT=$(realpath -e "$4")
[[ -d $OUT && -z $(ls -A -- "$OUT") ]] || { printf 'scratch must be an empty directory\n' >&2; exit 2; }
M=round2/occupation-neon/pressure-memory
R=$M/path-fast
H=$M/hv-lifetime/release
sha256sum --strict -c "$H/release.sha256"
"$P" -EBu "$H/verify.py" --payloads "$L" "$U"
g++ -std=c++20 -O3 -Wall -Wextra -Werror -shared -fPIC "$R/solver.cc" -o "$OUT/path-fast.so"
sha256sum "$OUT/path-fast.so"
run() {
  local tag=$1; shift
  /usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M exit=%x' "$P" -EBu "$R/replay.py" \
    --library "$OUT/path-fast.so" --lower "$L" --upper "$U" "$@" > "$OUT/$tag.log" 2>&1
  printf 'PASS %s\n' "$tag"
}
run parity --python-parity
run overflow --python-parity --include-old-band --side-cap 0 --mother-cap 0
run medium --degrees 9 7 6 10 --arities 3 3 3 2 --entries "$OUT/medium-entries.json"
wide=(--degrees 16 9 9 16 --arities 4 4 4 3 --side-cap 3 --budget 40000000)
run head "${wide[@]}" --parts minimum mother 0 1 --mother-cap 3 --entries "$OUT/wide-head-entries.json"
run tail "${wide[@]}" --parts 2 3 --mother-cap 3 --entries "$OUT/wide-tail-entries.json"
run mother12 "${wide[@]}" --parts mother --mother-cap 12 --entries "$OUT/wide-mother12-entries.json"
for file in medium-entries.json wide-head-entries.json wide-tail-entries.json wide-mother12-entries.json; do
  cmp "$R/$file" "$OUT/$file"
done
printf 'PASS four byte-identical conditional receipt files\n'
/usr/bin/time -f 'wall=%e cpu=%U+%S rss=%M exit=%x' "$P" -EBu "$R/pool.py" --limit 18 \
  "$OUT/wide-head-entries.json" "$OUT/wide-tail-entries.json" \
  --refine-mother "$OUT/wide-mother12-entries.json" > "$OUT/pool.log" 2>&1
printf 'PASS pooled fresh receipts\n'
"$P" -EBu "$R/low-release/verify.py" --compare-replays "$OUT"
