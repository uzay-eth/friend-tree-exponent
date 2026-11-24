#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")"
ulimit -v 4194304
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
run() {
    printf '\n$ timeout 5s taskset -c 0 python3 -B'
    printf ' %q' "$@"
    printf '\n'
    /usr/bin/time -f 'elapsed_s=%e cpu_user_s=%U max_rss_KiB=%M' \
        timeout 5s taskset -c 0 python3 -B "$@" 2>&1
}
printf 'Single core 0. Address limit %s KiB. Each diagnostic has a five-second wall limit.\n' "$(ulimit -v)"
run excursion_check.py
run excursion_rows.py --size 9 --s .5665 --test A
run excursion_rows.py --size 9 --s .5665 --test B
run excursion_rows.py --size 9 --s .54001 --test B --virtual
run excursion_rows.py --size 9 --s .565615 --test half
run excursion_rows.py --size 9 --s .565809 --test half
run excursion_rows.py --size 7 --s 1 --test scalar
run excursion_constants.py
