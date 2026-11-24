#!/usr/bin/env python3
"""Block under-dispersion test (registered 06:31 UTC, see RESULT.md section 3.3 (b)).

Input: chunk files of 5000 consecutive runs with checkpoints 12 and 16 (ftprod output, weight vectors in the order
sane, full, ...).  Residuals e = (num(16)-num(12)) - R (den(16)-den(12)) with R the ratio estimate of the whole
sample, for sane and for full; each standardised; z = z_sane + z_full.  Statistic
    T = sum over blocks (block sum of z)^2 / sum of z^2,
whose value is about 1 - 1/m for m blocks of independent runs.  The lower-tail p-value is the fraction of 1999
random reorderings of the runs with T <= the observed T (same reordering for all quantities).
Usage: block_test.py 'GLOB' [BLOCK]
"""
import sys, glob, re
import numpy as np
files = sorted(glob.glob(sys.argv[1]), key=lambda f: int(re.search(r"_c(\d+)", f).group(1))); bs = int(sys.argv[2]) if len(sys.argv) > 2 else 5000
rows = []
for f in files:
    a = np.loadtxt(f, ndmin=2)
    if a.shape[0] != 2 * bs: continue
    rows.append(a[:, :8])
A = np.concatenate(rows); a1 = A[A[:, 1] == 12]; a2 = A[A[:, 1] == 16]; o = np.argsort(a1[:, 0], kind="stable"); a1, a2 = a1[o], a2[o]
assert np.all(a1[:, 0] == a2[:, 0]) and np.all(np.diff(a1[:, 0]) == 1), "runs must be consecutive"
n = len(a1); m = n // bs
def resid(cn, cd):
    num = a2[:, cn] - a1[:, cn]; den = a2[:, cd] - a1[:, cd]; R = num.sum() / den.sum(); e = num - R * den; return R, e / e.std()
Rs, zs = resid(4, 5); Rf, zf = resid(6, 7); Rp, zp = resid(3, 2)
def T(x): return float((x[:m * bs].reshape(m, bs).sum(axis=1)**2).sum() / (x[:m * bs]**2).sum())
rng = np.random.default_rng(20261124); perms = [rng.permutation(n) for _ in range(1999)]
print(f"{len(rows)} blocks of {bs} runs ({n} runs, run indices {int(a1[0,0])}..{int(a1[-1,0])}); R_sane = {Rs:.10f}, R_full = {Rf:.10f}; correlation of the two residuals {np.corrcoef(zs, zf)[0,1]:.3f}")
for name, x in (("z_sane + z_full (registered statistic)", zs + zf), ("z_sane", zs), ("z_full", zf), ("z_sane - z_full", zs - zf), ("plain", zp)):
    t = T(x); null = np.array([T(x[p]) for p in perms]); print(f"   {name:40s} T = {t:.3f};  reorderings: mean {null.mean():.3f}, sd {null.std():.3f};  lower-tail p = {(np.sum(null <= t) + 1) / (len(null) + 1):.4f}")
