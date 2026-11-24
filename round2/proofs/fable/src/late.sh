#!/bin/bash
# late.sh -- late cohort (sensitivity check at a later time scale).
# Usage: late.sh OUTDIR THETAS DEADLINE_EPOCH "LAYOUT KEYS" [NWORKERS] [FIRST_CHUNK]
# Each worker (one core) writes chunks of 20 runs with checkpoints 20 24 28, stream id 9500; worker k uses run
# indices k*10^7 + chunk*20.  Uses the frozen production binary.  A worker stops before starting a chunk once
# the deadline has passed or the file OUTDIR/STOP exists.
out=$1; th=$2; deadline=$3; keys=$4; nw=${5:-2}; first=${6:-0}
BIN=$(dirname "$0")/../bin/ftprod.frozen-0446
mkdir -p $out
worker() {
  k=$1; chunk=$first
  while [ $(date +%s) -lt $deadline ] && [ ! -e $out/STOP ]; do
    $BIN 9500 $((k*10000000 + chunk*20)) 20 $out/late_w${k}_c$(printf %05d $chunk) $th 20 24 28 $keys
    chunk=$((chunk+1))
  done
}
for k in $(seq 0 $((nw-1))); do worker $k & done
wait
