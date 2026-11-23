#!/usr/bin/env python3
"""Local exponent from friend_ct output files.

mu_eff(t) = sum_runs S(t) / sum_runs C(t)  (ratio of sample means; dE[C]/dt = E[S] exactly), jackknife error
over runs.  Also the window estimator sum A / sum F over the core firings in (t-1, t] (A = sum of L/D at the
firing vertex just before it fires), and the pooled window estimator over all t >= T0.
usage: analyze_ct.py T0 file1 [file2 ...]
"""
import sys, math
from collections import defaultdict

def jack(num, den):
    n = len(num); sn = sum(num); sd = sum(den)
    est = sn / sd
    if n < 2:
        return est, float('nan')
    th = [(sn - a) / (sd - b) for a, b in zip(num, den)]
    m = sum(th) / n
    return est, math.sqrt((n - 1) / n * sum((x - m) ** 2 for x in th))

def main():
    T0 = float(sys.argv[1])
    runs = defaultdict(dict)   # (file, run) -> t -> row
    for fn in sys.argv[2:]:
        for line in open(fn):
            if line[0] in '#H':
                continue
            p = line.split()
            if len(p) < 9:
                continue
            runs[(fn, int(p[0]))][float(p[1])] = (int(p[2]), int(p[3]), float(p[4]), int(p[5]), int(p[6]), float(p[7]), int(p[8]))
    # keep only complete runs: a run killed before the last checkpoint of its file is dropped
    fmax = defaultdict(float)
    for (fn, r), d in runs.items():
        fmax[fn] = max(fmax[fn], max(d))
    dropped = [k for k, d in runs.items() if max(d) < fmax[k[0]]]
    for k in dropped:
        del runs[k]
    ts = sorted({t for r in runs.values() for t in r})
    tmax_all = max(ts)
    print(f"{len(runs)} complete runs ({len(dropped)} incomplete runs dropped); checkpoints up to t = {tmax_all}")
    print(" t    runs   mean n      mean C      mu_eff=S/C   +-        window A/F   +-        frac L=0")
    for t in ts:
        rs = [r[t] for r in runs.values() if t in r]
        if len(rs) < 10 or t < 8:
            continue
        e1, s1 = jack([x[2] for x in rs], [x[1] for x in rs])
        e2, s2 = jack([x[5] for x in rs], [x[3] for x in rs])
        z = sum(x[6] for x in rs) / sum(x[1] for x in rs)
        print(f"{t:4.0f} {len(rs):6d} {sum(x[0] for x in rs)/len(rs):11.4g} {sum(x[1] for x in rs)/len(rs):11.1f}"
              f"  {e1:.6f} {s1:.6f}   {e2:.6f} {s2:.6f}   {z:.5f}")
    # pooled window estimator for t > T0, per run, jackknife over runs
    num = []; den = []
    for r in runs.values():
        a = sum(v[5] for t, v in r.items() if t > T0); f = sum(v[3] for t, v in r.items() if t > T0)
        if f > 0:
            num.append(a); den.append(f)
    e, s = jack(num, den)
    print(f"pooled window estimator, all core firings at t > {T0}: mu = {e:.7f} +- {s:.7f}  ({len(num)} runs, {sum(den):.3g} firings)")
    num = []; den = []
    for r in runs.values():
        a = sum(v[4] for t, v in r.items() if t > T0); f = sum(v[3] for t, v in r.items() if t > T0)
        if f > 0:
            num.append(a); den.append(f)
    e, s = jack(num, den)
    print(f"same with raw conversion counts Cv/F:                       mu = {e:.7f} +- {s:.7f}")

if __name__ == '__main__':
    main()
