#!/usr/bin/env bash
set -euo pipefail
cd "$(dirname "$0")/../../../.."
name=$1
seconds=$2
shift 2
[[ $name =~ ^[a-z0-9-]+$ && $seconds =~ ^[0-9]+$ && $seconds -le 600 ]]
ulimit -v 4194304
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
cmd=(timeout "${seconds}s" taskset -c 0 python3 -u -B "round2/finite/gauge/low-sector/$1")
shift
cmd+=("$@")
set -o noclobber
{
    date -u
    printf '$'
    printf ' %q' "${cmd[@]}"
    printf '\n'
    /usr/bin/time -f 'wall_seconds %e user_seconds %U sys_seconds %S max_rss_kib %M exit %x' "${cmd[@]}"
    date -u
} > "round2/outputs/finite-low-$name.log" 2>&1
