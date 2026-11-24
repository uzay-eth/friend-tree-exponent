#!/bin/bash
# final.sh -- reduce all raw outputs to data/moments.npz and print the final numbers.  Run from round2/.
# Raw files (outputs/, not archived) -> data/moments.npz (archived) -> the numbers of RESULT.md sections 2.5 and 3.1.
set -e
cd "$(dirname "$0")/.."
PY=.venv/bin/python; MU=0.5657236; MU0=0.4822669183233687
$PY src/moments_reduce.py data/moments.npz \
  "main:outputs/prod/main_w*_c*.txt:400:$MU::500" "mainA:outputs/prod/main_w*_c*.txt:400:$MU::550" "main1:outputs/prod/main_w*_c*.txt:400:$MU::0:499" \
  "r1216:outputs/prod/r1216_c*.txt:10000:$MU::500" "r1216p1:outputs/prod/r1216_c*.txt:10000:$MU::0:499" \
  "r1620:outputs/prod/r1620_c*.txt:1000:$MU::500" "r1620p1:outputs/prod/r1620_c*.txt:1000:$MU::0:499" \
  "late:outputs/late/late_w*_c*.txt:60:$MU" "b1216:outputs/btest/b_c*.txt:10000:$MU" \
  "e0408:outputs/rungs_early/r04_08.txt:0:$MU:10000" "e0812:outputs/rungs_early/r08_12.txt:0:$MU:10000" \
  "e0408five:outputs/rungs_early/r04_08_five.txt:0:$MU:10000" "e0812five:outputs/rungs_early/r08_12_five.txt:0:$MU:10000" \
  "ir0early:outputs/valid/ir0_early.txt:0:$MU0:10000" "ir0ladder:outputs/valid/ir0_ladder.txt:0:$MU0:1000"
$PY src/final_analysis.py data/moments.npz --value-file data/final_value.txt
