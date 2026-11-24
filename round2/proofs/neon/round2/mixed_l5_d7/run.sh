#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1
ulimit -c 0
d=$1 b=$2
shift 2
args=(--D0 12 --KT 11 --LF 6000 --KF 100 --DP4 10 --KP4 8 --DYT 9 --KYT 7
      --DYB 10 --KYB 8 --DS3 8 --KS3 6 --nit 18 --threads 1
      --l5 "$d" --coarse "$d" --coarse-buckets "$b" --coarse-mask 3)
if (($#)); then
    exec "${SOLVER:?}" "${args[@]}" "$@"
fi
for mode in 2 1; do
    /usr/bin/time -f 'TIME wall=%e user=%U sys=%S peak_KiB=%M exit=%x' \
        stdbuf -oL "${SOLVER:?}" "${args[@]}" --mode "$mode" >"d${d}b${b}m${mode}.log" 2>&1
done
