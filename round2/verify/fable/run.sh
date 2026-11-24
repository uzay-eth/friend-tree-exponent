#!/usr/bin/env bash
set -euo pipefail
root=$(cd -- "$(dirname -- "$0")/../../.." && pwd -P)
work=$(mktemp -d "${TMPDIR:-/tmp}/fable.XXXXXXXX")
trap 'rm -rf -- "$work"' EXIT
proofs="$root/round2/proofs/fable"
verify="$root/round2/verify/fable"
hub="$proofs/newground/q3-hub-probability"
order="$proofs/newground/stochastic-domination"
paths="$proofs/newground/q5-low-degree-forest/longest-path"
cc=${CXX:-g++}
py=${PYTHON:-python3}
mode=${1:-quick}

ir0() {
    "$cc" -O2 -UNDEBUG -std=c++17 -frounding-math -Wno-unknown-pragmas "$proofs/src/ir0_rootcount_check.cpp" -o "$work/ir0"
    "$work/ir0"
    printf '%s\n' 'RESULT fable-ir0: certified root counts and planted-extra-zero control for the IR-0 cousin'
}

case "$mode" in
quick)
    "$py" -I -B -c 'import sys; sys.exit(0 if __debug__ else "Python assertions unavailable")'
    printf '%s\n' '== F1: finite-window statistics from authenticated selected moments =='
    (cd "$proofs"; sha256sum --check DATA.sha256)
    "$py" -I -B "$verify/f1.py" "$proofs/data/f1-moments.tsv"
    "$cc" -O2 -UNDEBUG -std=c++17 "$proofs/src/test_feat2.cpp" -o "$work/test_feat2"
    for model in 0 1 2; do
        "$work/test_feat2" 7 6 8 Rsum=3 Rpow=3 Dp=7 Dk=8 L1cap=8 "model=$model"
    done
    printf '%s\n' 'RESULT fable-window: final finite-window estimates, conditional bounds, controls, and small generator checks'

    printf '%s\n' '== F3: exact history certificate and small rooted-DP check =='
    "$cc" -O3 -UNDEBUG -std=c++17 "$hub/exact.cpp" -o "$work/exact"
    for m in 3 4 5 9 10 11; do
        "$work/exact" "$m" < "$hub/exact-$m.in"
    done > "$work/exact.log"
    (cd "$work"; "$py" -I -B "$hub/tails.py") | cmp - "$hub/exact-summary.log"
    printf '%s\n' 'Exact history summary matches through size 11 (3,628,800 histories).'
    "$cc" -O2 -UNDEBUG -std=c++17 "$hub/lead-signed_orbit_dp.cpp" -lgmpxx -lgmp -o "$work/hub-dp"
    "$work/hub-dp" 2 7 1 9 > "$work/hub.out"
    "$py" -I -B - "$work/hub.out" <<'PY'
from fractions import Fraction
import sys
with open(sys.argv[1]) as f:
    rows = {(int(u), int(n)): tuple(map(Fraction, (lo, hi)))
            for kind, u, n, _, _, lo, hi in (s.split() for s in f if s.startswith('H '))}
assert set(rows) == {(u, u+k) for u in range(2,8) for k in (1,2)}
assert all(lo <= hi for lo, hi in rows.values())
assert all(max(rows[u,u+1][0], rows[u,u+2][0]) > 0 for u in range(2,8))
PY
    printf '%s\n' 'RESULT fable-hubs: exact ordering through label 8; size-11 receipt and separate size-9 rooted DP'

    printf '%s\n' '== F5: exact sparse transport and local criterion (modest sizes) =='
    "$cc" -O2 -UNDEBUG -std=c++17 "$order/lead-transport_sparse.cpp" -lgmpxx -lgmp -o "$work/transport"
    "$work/transport" 2 11 > "$work/transport.out"
    cmp "$work/transport.out" "$order/lead-transport_sparse_2_11.log"
    "$cc" -O2 -UNDEBUG -std=c++17 "$order/lead-local_scalar.cpp" -lgmpxx -lgmp -o "$work/local"
    "$work/local" 4 11 > "$work/local.out"
    awk '($0 !~ ("^pair \\(" (NR+3) "," (NR+4) "\\):") || $0 !~ /PASS \(stochastic domination certified by the local criterion\)$/) {bad=1} END {if (bad || NR!=8) exit 1}' "$work/local.out"
    "$py" -I -B "$order/lead-domination_transport.py" 8 | tee "$work/dense.out"
    awk 'END {if (ok != 7) exit 1} /FULL STOCHASTIC DOMINATION CERTIFIED/ {ok++}' "$work/dense.out"
    printf '%s\n' 'RESULT fable-transport: sparse pairs (2,3) to (11,12), dense size 8, and local pairs (4,5) to (11,12)'

    printf '%s\n' '== F6/F7: bare-triple repair fixture and exact degree-2 path census =='
    "$py" -I -B "$proofs/newground/law-of-X/density/referee-check.py"
    printf '%s\n' 'RESULT fable-density-fixtures: exact star-exit and moment arithmetic; unmodified bare-triple counterexample'
    "$cc" -O2 -UNDEBUG -std=c++17 "$paths/sim.cpp" -o "$work/paths"
    "$work/paths" check > "$work/path-check.out"
    "$py" -I -B "$paths/referee-enumerate_law.py" > "$work/path-enum.out"
    printf '%s\n' 'RESULT fable-path-fixtures: degree-2 path exact finite-size enumerations'
    ir0
    ;;
ir0)
    ir0
    ;;
*)
    printf 'Usage: %s [quick|ir0]\nFull-size commands: round2/proofs/fable/README.md\n' "$0" >&2
    exit 2
    ;;
esac
