#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1
ulimit -c 0
run() {
    local case=$1; shift
    for mode in 1 2; do
        local log="${case}_m${mode}.log"
        if [[ -f "$log" ]] && grep -q '^TIME .* exit=0$' "$log"; then continue; fi
        /usr/bin/time -f 'TIME wall=%e user=%U sys=%S peak_KiB=%M exit=%x' \
            "${SOLVER:-${TMPDIR:-.}/neon-mixed-diag}" --threads 1 --D0 12 --KT 11 \
            --LF 6000 --KF 100 --DP4 10 --KP4 8 --DYT 9 --KYT 7 --DYB 10 --KYB 8 --DS3 8 --KS3 6 \
            --mode "$mode" --nit 18 --l5 6 "$@" > "$log" 2>&1
    done
}
abc() {
    run zeroA --slot-mask 1
    run zeroB --slot-mask 2
    run zeroC --slot-mask 4
}
rest() {
    run clipA --slot-mask 1 --slot-denom 6
    run clipB --slot-mask 2 --slot-denom 6
    run clipAB --slot-mask 3 --slot-denom 6
    run clipC --slot-mask 4 --slot-denom 6
    run reference
    run zeroABC --slot-mask 7
}
case "${1:-grid}" in
    abc) abc ;;
    rest) rest ;;
    grid) abc; rest ;;
    coupled) run coupled --coarse 6 --coarse-mask "${COARSE_MASK:-3}" ;;
    *) exit 2 ;;
esac
