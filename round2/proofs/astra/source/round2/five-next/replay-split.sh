#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../.."
if [[ ${FIVE_NEXT_SPLIT_PINNED:-} != 1 ]]; then
    core=$(python3 -c 'import os; print(min(os.sched_getaffinity(0)))')
    exec env FIVE_NEXT_SPLIT_PINNED=1 taskset -c "$core" nice -n 10 bash round2/five-next/replay-split.sh "$@"
fi
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1
ulimit -v 5242880
dir=round2/five-next
sha256sum -c "$dir/split-inputs.sha256"
case ${1:-cold} in
    frozen)
        python3 -B "$dir/split_integers.py" round2/outputs/five-next-split.log
        ;;
    cold)
        mkdir -p "$dir/build"
        scratch=$(mktemp -d "$dir/build/split-audit.XXXXXXXX")
        root=$PWD
        mkdir -p "$scratch/round2/five-next/build"
        flags=(-std=c++20 -O3 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=undefined)
        for program in reserve gradient_probe census; do
            c++ "${flags[@]}" "$dir/$program.cpp" -o "$scratch/$program"
        done
        python3 -B "$dir/check_slopes.py" "$scratch"
        python3 -B "$dir/check.py" "$scratch" small-split
        (
            cd "$scratch"
            for side in lower upper; do
                /usr/bin/time -f 'wall %e s, max RSS %M KiB' ./reserve "$root/round2/outputs/mc-$side.bin" split 5
            done
        ) 2>&1 | tee "$scratch/output.log"
        python3 -B "$dir/split_integers.py" "$scratch/output.log"
        cp "$scratch/output.log" round2/outputs/five-next-split-cold.log
        ;;
    *) echo 'usage: bash round2/five-next/replay-split.sh [cold|frozen]' >&2; exit 2 ;;
esac
