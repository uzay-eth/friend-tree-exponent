#!/usr/bin/env bash
# Reconstruct the archived run lengths; compare decompressed text.
set -euo pipefail
cd "$(dirname "$0")"
case "${1:-}" in
    first|full) mode=$1 ;;
    *) echo "usage: bash replay.sh first|full" >&2; exit 2 ;;
esac
out="build/replay-$mode"
mkdir -p "$out"
g++ -O2 -o "$out/friend_ct" sim/friend_ct.cpp

replay() {
    local batch=$1 seed=$2 runs=$3 tmax=$4
    local stem="ct_t${tmax}_seed${batch}"
    gzip -dc "outputs/$stem.txt.gz" > "$out/$stem.archive.txt"
    if [[ $mode == first ]]; then
        runs=1
        head -n "$((tmax + 1))" "$out/$stem.archive.txt" > "$out/$stem.expected.txt"
    else
        cp "$out/$stem.archive.txt" "$out/$stem.expected.txt"
    fi
    echo "friend_ct $seed $runs $tmax 1"
    time "$out/friend_ct" "$seed" "$runs" "$tmax" 1 > "$out/$stem.txt"
    cmp "$out/$stem.expected.txt" "$out/$stem.txt"
    echo "$stem: byte-identical ($runs runs)"
}

replay A 20261122001 527 28
replay B 20261122002 294 28
replay C 20261122003 91 31
echo "friend_ct 4242 1 27 1 1"
time "$out/friend_ct" 4242 1 27 1 1 > "$out/hist_t27.txt"
cmp outputs/hist_t27.txt "$out/hist_t27.txt"
echo 'hist_t27: byte-identical'
