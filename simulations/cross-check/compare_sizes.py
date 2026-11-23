#!/usr/bin/env python3
"""Compare the packaged fixed-size sample with the independent direct chain at the same sizes.

Reads the packaged per-trajectory core counts (fixed-size/core_counts.tsv.gz, one row per run,
one column per size) and one or more outputs of `naive_chain sizes`, pools the latter, and prints a
two-sample z-score for the mean core size at every common size.  Sizes share trajectories on each
side, so the z-scores are correlated across rows.
Usage: compare_sizes.py CORE_COUNTS.tsv.gz sizes_1.txt [sizes_2.txt ...]
"""
import gzip, math, sys

with gzip.open(sys.argv[1], 'rt') as fh:
    head = fh.readline().split()[1:]
    rows = [list(map(int, line.split()[1:])) for line in fh if line.strip()]
packaged = {}
for j, size in enumerate(map(int, head)):
    xs = [r[j] for r in rows]
    m = sum(xs) / len(xs)
    v = sum((x - m) ** 2 for x in xs) / (len(xs) - 1)
    packaged[size] = (len(xs), m, math.sqrt(v / len(xs)))

pool = {}
for path in sys.argv[2:]:
    R = None
    for line in open(path):
        if line.startswith('#'):
            R = int(line.split('runs=')[1].split()[0]); continue
        f = line.split(); n = int(f[1]); d = dict(zip(f[2::2], map(float, f[3::2])))
        p = pool.setdefault(n, [0, 0.0, 0.0])
        p[0] += R; p[1] += R * d['C_mean']; p[2] += (R - 1) * d['C_sd'] ** 2 + R * d['C_mean'] ** 2

print(f"{'N':>8} {'pkg runs':>8} {'pkg mean C':>11} {'+-':>8} {'direct runs':>11} {'direct mean C':>13} {'+-':>8} {'z':>6}")
worst = 0.0
for n in sorted(pool):
    if n not in packaged:
        continue
    R, s, ss = pool[n]
    m = s / R; se = math.sqrt((ss - R * m * m) / (R - 1) / R)
    Rp, mp, sep = packaged[n]
    z = (mp - m) / math.hypot(sep, se)
    worst = max(worst, abs(z))
    print(f"{n:>8} {Rp:>8} {mp:>11.4f} {sep:>8.4f} {R:>11} {m:>13.4f} {se:>8.4f} {z:>6.2f}")
print(f"max |z|: {worst:.2f} (correlated across sizes)")
