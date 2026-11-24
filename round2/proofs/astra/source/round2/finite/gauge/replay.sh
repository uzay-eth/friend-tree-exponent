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
run binary_weight.py
run virtual_weight.py
run constant_obstruction.py
run binary_rows.py --virtual --s .54001
run binary_rows.py --virtual --s .5665
run binary_rows.py --virtual --s .565615
run return_constants.py
run dormant_check.py
run transfer_check.py
run depth_error.py
printf '\nFloating diagnostics for sufficient inequalities only.\n'
run virtual_probe.py
run current_parent_probe.py
