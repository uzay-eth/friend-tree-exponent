#!/bin/bash
# prod.sh -- production driver.  Usage: prod.sh OUTDIR THETAS DEADLINE_EPOCH "LAYOUT KEYS" [FIRST_CHUNK]
# FIRST_CHUNK (default 0) lets a later invocation continue with fresh run indices (chunk numbers must not repeat).
# Starts 4 worker loops (one core each).  Workers 0-2: main window, checkpoints 12 16 20 24, stream id 9100,
# chunks of 100 runs, worker k uses run indices k*10^7 + ...  Worker 3: bias-monitor rungs, alternating chunks of
# 5000 runs with checkpoints 12 16 (stream 9200) and 500 runs with checkpoints 16 20 (stream 9300).
# Every chunk is a separate complete output file, so finished chunks are usable at any time.  A worker stops
# before starting a chunk once the deadline has passed or the file OUTDIR/STOP exists.
out=$1; th=$2; deadline=$3; keys=$4; first=${5:-0}
BIN=$(dirname "$0")/../bin/ftprod
mkdir -p $out
main_worker() {
  k=$1; chunk=$first
  while [ $(date +%s) -lt $deadline ] && [ ! -e $out/STOP ]; do
    r0=$((k*10000000 + chunk*100))
    $BIN 9100 $r0 100 $out/main_w${k}_c$(printf %05d $chunk) $th 12 16 20 24 $keys
    chunk=$((chunk+1))
  done
}
rung_worker() {
  chunk=$first
  while [ $(date +%s) -lt $deadline ] && [ ! -e $out/STOP ]; do
    $BIN 9200 $((chunk*5000)) 5000 $out/r1216_c$(printf %05d $chunk) $th 12 16 $keys
    $BIN 9300 $((chunk*500)) 500 $out/r1620_c$(printf %05d $chunk) $th 16 20 $keys
    chunk=$((chunk+1))
  done
}
for k in 0 1 2; do main_worker $k & done
rung_worker &
wait
