#!/usr/bin/env bash
# Run from any directory. All generated files live under simulations/build/demo.
set -euo pipefail
cd "$(dirname "$0")"
out=build/demo
mkdir -p "$out"
g++ -O2 -o "$out/friend_ct" sim/friend_ct.cpp
"$out/friend_ct" 20261123 64 17 1 > "$out/demo.txt"
echo 'Fresh fixed-time sample, seed 20261123 (64 runs through t=17)'
python3 local_slopes.py --sample "$out/demo.txt" 64 17 --table "$out/demo-slopes.tsv"
awk -F '\t' 'NR == 1 { print " t       mean size   instantaneous slope     SE"; next }
    { printf "%2d  %14.1f             %.6f   %.6f\n", $1, $3, $5, $6 }' "$out/demo-slopes.tsv"
python3 analyze_ct.py 12 "$out/demo.txt" > "$out/demo-summary.txt"
tail -2 "$out/demo-summary.txt"
echo 'Full-size archived data'
for name in ct_t28_seedA ct_t28_seedB ct_t31_seedC; do
    gzip -dc "outputs/$name.txt.gz" > "$out/$name.txt"
done
python3 analyze_ct.py 21 "$out"/ct_*.txt > "$out/summary.txt"
diff -u outputs/sim_summary.txt "$out/summary.txt"
tail -2 "$out/summary.txt"
python3 local_slopes.py --table "$out/local_slopes.tsv"
diff -u local_slopes.tsv "$out/local_slopes.tsv"
echo 'Reproduced archived summary and full checkpoint table.'
python3 fit_window.py --table "$out/fixed-size.tsv" > "$out/fits.txt"
diff -u fixed-size/summary.tsv "$out/fixed-size.tsv"
diff -u finite_window.txt "$out/fits.txt"
cat "$out/fits.txt"
