"""Sample the annealed fixed-time population (a diagnostic, not a certificate)."""
import argparse
from collections import Counter
from pathlib import Path
import struct
import subprocess

import numpy as np


def trees(binary, reps, horizon, seed):
    with subprocess.Popen([str(Path(binary).resolve()), str(reps), str(horizon), str(seed)],
                          stdout=subprocess.PIPE) as proc:
        for _ in range(reps):
            n, = struct.unpack('=i', proc.stdout.read(4))
            yield np.frombuffer(proc.stdout.read(12*n), dtype=np.int32).reshape(n, 3)
        if proc.wait():
            raise RuntimeError('simulator failed')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary', default='./sim_core')
    ap.add_argument('--reps', type=int, default=2000)
    ap.add_argument('--time', type=float, default=12)
    ap.add_argument('--seed', type=int, default=38129)
    ap.add_argument('--output')
    args = ap.parse_args()
    p, e, chain, wedge = (Counter() for _ in range(4))
    checks = []
    saved = []
    for arr in trees(args.binary, args.reps, args.time, args.seed):
        if args.output:
            saved.append(arr)
        n = len(arr)
        c, d, par = arr.T
        promotion = 1 - c / d
        rho = np.zeros(n)
        np.add.at(rho, par[1:], 1/d[1:])
        rho[1:] += 1/d[par[1:]]
        vals = np.r_[n, sum(promotion), sum(c), sum(1/d), sum(1./d**2),
                     [sum((d-c+rho)[d == k]) - .56572366*sum(d > k) for k in range(2, 11)]]
        checks.append(vals)
        s = list(zip(d.tolist(), c.tolist()))
        p.update(s)
        e.update((s[j], s[i]) for i, j in enumerate(par[1:], 1))
        chain.update((s[par[j]], s[j], s[i]) for i, j in enumerate(par[1:], 1) if j)
        hist = [Counter() for _ in s]
        for i, j in enumerate(par[1:], 1):
            hist[j][s[i]] += 1
        for i, co in enumerate(hist):
            for x, nx in co.items():
                for y, ny in co.items():
                    if ny > (x == y):
                        wedge[s[i], x, y] += nx*(ny-(x == y))
    a = np.array(checks)
    means = a.sum(axis=0)[1:] / a[:, 0].sum()
    residual = a[:, 1:] - a[:, :1]*means
    stderr = residual.std(axis=0, ddof=1) / (len(a)**.5*a[:, 0].mean())
    print(f'cores={int(a[:,0].sum())} trees={len(a)} singleton/edge/chain/wedge keys='
          f'{len(p)}/{len(e)}/{len(chain)}/{len(wedge)}')
    for label, x, se in zip(['b', 'c', '1/D', '1/D^2']+[f'flux[{k}]' for k in range(2, 11)], means, stderr):
        print(f'{label:9} {x:+.9f} (tree-cluster SE {se:.3g})')
    for label, table in [('p', p), ('edge', e), ('chain', chain), ('ordered wedge', wedge)]:
        print(label, 'top 3 per core:', [(s, n/a[:,0].sum()) for s, n in table.most_common(3)])
    if args.output:
        np.savez(args.output, cores=np.concatenate(saved), offsets=np.cumsum([0]+list(map(len, saved))))


if __name__ == '__main__':
    main()
