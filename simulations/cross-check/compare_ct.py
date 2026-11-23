#!/usr/bin/env python3
"""Compare the packaged fast CT simulator (sim/friend_ct.cpp, run fresh) with the independent direct
chain observed at an exact Yule(3) size draw (naive_chain yule).  Cross-check, Claude Opus 5.5.

Both estimate E[C(t)] and E[S(t)] for the process started at T_3 at time 0.  Two-sample z-scores use
independent runs on each side.  Also checks the CT mean size against 3 e^t.
Usage: compare_ct.py "out_ct_*.txt" out_yule_t*.txt [--label NAME]
"""
import sys, math, glob

ct = {}
for path in sorted(glob.glob(sys.argv[1])):
    for line in open(path):
        f = line.split()
        if not f or f[0] != 't':
            continue
        t = float(f[1]); d = dict(zip(f[2::2], map(float, f[3::2])))
        a = ct.setdefault(t, dict.fromkeys(['runs', 'sum_n', 'sum_C', 'sum_C2', 'sum_S', 'sum_S2'], 0.0))
        for k in a:
            a[k] += d[k]

def ms(s, s2, R):
    m = s / R; v = (s2 - R * m * m) / (R - 1)
    return m, math.sqrt(max(v, 0.0) / R)

print(f"{'t':>4} {'CT runs':>9} {'CT mean n':>11} {'3e^t':>11} {'CT E[C]':>11} {'+-':>8} {'direct E[C]':>12} {'+-':>8} {'zC':>6} {'CT E[S]':>11} {'direct E[S]':>12} {'zS':>6}  direct-file")
worst = None
for path in sys.argv[2:]:
    for line in open(path):
        if line.startswith('#'):
            hdr = line; R = int(line.split('runs=')[1].split()[0]); mode = line.split()[1]
            continue
        f = line.split(); t = float(f[1]); d = dict(zip(f[2::2], map(float, f[3::2])))
        a = ct.get(t)
        if a is None:
            print(f"{t:>4g}  no CT data"); continue
        Rc = a['runs']
        mc, ec = ms(a['sum_C'], a['sum_C2'], Rc); msv, es = ms(a['sum_S'], a['sum_S2'], Rc)
        zc = (mc - d['C_mean']) / math.hypot(ec, d['C_se'])
        zs = (msv - d['S_mean']) / math.hypot(es, d['S_se'])
        if mode == 'yule' and 'model=0' in hdr:
            worst = max(worst or 0.0, abs(zc), abs(zs))
        print(f"{t:>4g} {int(Rc):>9} {a['sum_n']/Rc:>11.5g} {3*math.exp(t):>11.5g} {mc:>11.5f} {ec:>8.5f} {d['C_mean']:>12.5f} {d['C_se']:>8.5f} {zc:>6.2f} {msv:>11.5f} {d['S_mean']:>12.5f} {zs:>6.2f}  {path} [{mode}, {R} runs]")
print("max |z| over friend-tree (model 0) Yule-matched comparisons: " + ("none in this set" if worst is None else f"{worst:.2f}"))
