#!/bin/bash
# late_worker.sh -- one additional worker of the late cohort (same stream 9500, same binary, same checkpoints as
# late.sh); worker index K gives run indices K*10^7 + chunk*20.  Usage: late_worker.sh OUTDIR THETAS DEADLINE_EPOCH "KEYS" K [FIRST_CHUNK]
out=$1; th=$2; deadline=$3; keys=$4; k=$5; chunk=${6:-0}
BIN=$(dirname "$0")/../bin/ftprod.frozen-0446
while [ $(date +%s) -lt $deadline ] && [ ! -e $out/STOP ]; do
  $BIN 9500 $((k*10000000 + chunk*20)) 20 $out/late_w${k}_c$(printf %05d $chunk) $th 20 24 28 $keys
  chunk=$((chunk+1))
done
