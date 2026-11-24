#!/usr/bin/env bash
set -euo pipefail

root=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../../.." && pwd)
src=$root/round2/proofs/neon2
q=$src/quintet-cert
m=$src/stationary-marginals
py=${PYTHON:-python3}
cxx=${CXX:-g++}
threads=${NEON2_THREADS:-2}
export OMP_THREAD_LIMIT=${OMP_THREAD_LIMIT:-$threads}
export OMP_WAIT_POLICY=${OMP_WAIT_POLICY:-PASSIVE}
unset PYTHONOPTIMIZE

usage() {
    printf '%s\n' 'Usage: run.sh [quick | q3-full | gold1-build DIR | gold2-build DIR | gold1-check DIR | gold2-check DIR | marginals-build DIR | marginals-check DIR]' >&2
    exit 2
}

mode=${1:-quick}
case "$mode" in
    quick|q3-full) [[ $# -le 1 ]] || usage ;;
    gold[12]-build|gold[12]-check|marginals-build|marginals-check) [[ $# -eq 2 ]] || usage ;;
    *) usage ;;
esac
work=$(mktemp -d "${TMPDIR:-/tmp}/neon2.XXXXXXXX")
trap 'rm -r -- "$work"' EXIT

compile() {
    "$cxx" -std=c++17 -O3 -fopenmp -Wall -Wextra -Werror -pedantic "$1" -o "$2"
}
build_marginal() {
    "$work/build_marginal" "$1" "$2" "$3" 565713486333 1000000000000 565736986679 1000000000000 "$4" >&2
}

q3_full() {
    "$py" -B "$src/q3/cert_q3.py" --size 16 > "$work/cert_n16.txt"
    cmp "$src/q3/cert_n16.txt" "$work/cert_n16.txt"
    result q3 'size-16 exact receipt matched for order through label 8'
}

result() { printf 'RESULT neon2-%s: %s\n' "$1" "$2"; }

case "$mode" in
    quick)
        "$py" -B "$src/limit-density/check_density.py"
        result density 'exact split fixtures and fixed-seed Beta diagnostic'
        "$py" -B "$src/q3/symmetries.py"
        q3_full
        "$py" -B "$src/claim14-proofs/verify_track3.py"
        result fringe 'finite signed-generator, rerooting and stopped-prefix fixtures'
        "$py" -Bu "$src/obstruction/verify_transcendence.py"
        result transcendence 'exact finite response coefficients and valuation fixtures'
        "$py" -B "$src/q5/verify_q5_exact.py"
        result q5 'exact cutoff-score, holding-time and boundary fixtures'
        "$py" -B "$src/density-lp/verify_density_lp.py"
        result density-lp 'packing, rate and exponent fixtures; floating-point diagnostics'
        "$py" -OB "$src/mean-order/cert_mean_order.py"
        result mean-order 'exact finite drift tables and controls'
        "$py" -OB "$src/all-moments-order/cert_all_moments.py"
        result moments-order 'exact boxes for real p in [1,12] and analytic-tail checks'
        "$py" -OB "$src/height-diameter/cert_height.py"
        result height-diameter 'exact depth-potential and branching-process rows'
        "$py" -OB "$src/fractional-moments-order/cert_fractional_moments.py"
        result fractional-moments 'exact boxes for real p in [1/2,1] and controls'

        compile "$q/build_quintet.cpp" "$work/build_quintet"
        compile "$q/verify_quintet.cpp" "$work/verify_quintet"
        for side in lower upper; do
            "$work/build_quintet" tiny "$side" "$work/quintet_$side.bin"
        done
        "$py" -B "$q/cert_check.py" --verify-bin "$work/verify_quintet" --lower "$work/quintet_lower.bin" --upper "$work/quintet_upper.bin" --threads "$threads" --run-mutations > "$work/quintet.json"
        "$py" -B "$q/audit_independent.py" "$work/quintet.json" "$work/quintet_lower.bin" "$work/quintet_upper.bin"
        result quintet-tiny 'tiny builder, exact checker and controls; neither gold receipt replayed'

        compile "$m/build_marginal.cpp" "$work/build_marginal"
        compile "$m/verify_marginal.cpp" "$work/verify_marginal"
        for side in lower upper; do
            build_marginal tiny 1 "$side" "$work/obs1_$side.bin"
        done
        "$py" -OB "$m/cert_marginal.py" --verify-bin "$work/verify_marginal" --lower "$work/obs1_lower.bin" --upper "$work/obs1_upper.bin" --threads "$threads" --run-mutations > "$work/marginal.json"
        result marginals-tiny 'tiny observable-1 builder, exact checker and controls; production receipt not replayed'
        ;;
    q3-full)
        q3_full
        ;;
    gold[12]-build)
        mkdir -p -- "$2"
        compile "$q/build_quintet.cpp" "$work/build_quintet"
        for side in lower upper; do
            [[ ! -e "$2/quintet_$side.bin" ]] || { printf 'Output already exists: %s\n' "$2/quintet_$side.bin" >&2; exit 2; }
            "$work/build_quintet" "${mode%-build}" "$side" "$2/quintet_$side.bin"
        done
        ;;
    gold[12]-check)
        compile "$q/verify_quintet.cpp" "$work/verify_quintet"
        "$py" -B "$q/cert_check.py" --verify-bin "$work/verify_quintet" --lower "$2/quintet_lower.bin" --upper "$2/quintet_upper.bin" --threads "$threads" > "$work/replay.json"
        "$py" -B "$q/audit_independent.py" "$work/replay.json" "$2/quintet_lower.bin" "$2/quintet_upper.bin" >&2
        cat "$work/replay.json"
        if cmp -s "$q/certificate_${mode%-check}.json" "$work/replay.json"; then
            printf '%s\n' 'Historical receipt matched byte for byte.' >&2
        else
            printf '%s\n' 'Historical receipt did not match; only the newly verified endpoints apply.' >&2
        fi
        ;;
    marginals-build)
        mkdir -p -- "$2"
        compile "$m/build_marginal.cpp" "$work/build_marginal"
        for obs in {1..8}; do
            for side in lower upper; do
                [[ ! -e "$2/obs${obs}_$side.bin" ]] || { printf 'Output already exists: %s\n' "$2/obs${obs}_$side.bin" >&2; exit 2; }
                build_marginal high "$obs" "$side" "$2/obs${obs}_$side.bin"
            done
        done
        ;;
    marginals-check)
        compile "$m/verify_marginal.cpp" "$work/verify_marginal"
        "$py" -OB "$m/verify_all_marginals.py" --verify-bin "$work/verify_marginal" --dir "$2" --threads "$threads" > "$work/replay.json"
        cat "$work/replay.json"
        if cmp -s "$m/certificates_all.json" "$work/replay.json"; then
            printf '%s\n' 'Historical receipt matched byte for byte.' >&2
        else
            printf '%s\n' 'Historical receipt did not match; only the newly verified endpoints apply, conditional on the original gold2 bracket.' >&2
        fi
        ;;
esac
