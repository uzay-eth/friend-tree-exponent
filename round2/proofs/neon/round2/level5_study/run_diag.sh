#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1
run() {
    local d=$1 case=$2 nit=$3; shift 3
    for mode in 1 2; do
        local tag="d${d}_${case}_m${mode}"
        if [[ -f "$tag.log" ]] && grep -q '^TIME .* exit=0$' "$tag.log"; then continue; fi
        /usr/bin/time -f 'TIME wall=%e user=%U sys=%S peak_KiB=%M exit=%x' \
            "${SOLVER:-./diag_slack}" --threads 1 --D0 "$d" --KT "$((d-1))" --mode "$mode" \
            --nit "$nit" "$@" > "$tag.log" 2>&1
    done
}
run_grid() {
for d in 12; do
    for scenario in base center d5-6 d5-7 d5-8 h-center; do
        run "$d" "$scenario" 18 --scenario "$scenario"
    done
done
}

run_extensions() {
    # Exact verification throughput; the certificate is local to this study.
    /usr/bin/time -f 'TIME wall=%e user=%U sys=%S peak_KiB=%M exit=%x' \
        ./diag_exact --threads 1 --mode 2 --D0 12 --KT 11 --nit 18 --out sample.bin > sample_build.log 2>&1
    /usr/bin/time -f 'TIME wall=%e user=%U sys=%S peak_KiB=%M exit=%x' \
        ./verify_gmp sample.bin --threads 1 --tree-check 25 > sample_verify.log 2>&1
}

case "${1:-grid}" in
    grid) run_grid ;;
    extensions) run_extensions ;;
    ladder)
        run 12 l4-large 18 --DP4 12 --KP4 10 --DYT 11 --KYT 9 --DYB 12 --KYB 10 --DS3 10 --KS3 8
        run 14 l4-large 18 --DP4 12 --KP4 10 --DYT 11 --KYT 9 --DYB 12 --KYB 10 --DS3 10 --KS3 8
        ;;
    l5)
        export SOLVER=./solve_l5
        run 12 l5-6 16 --l5 6
        run 12 l5-6-no-inner-slack 16 --l5 6 --l5-split 1
        run 12 l5-6-no-mixed-slack 16 --l5 6 --l5-split 2
        ;;
    *) exit 2 ;;
esac
