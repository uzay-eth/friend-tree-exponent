#!/usr/bin/env python3
"""Eigenvalues of the Monte Carlo renewal kernel K_s on a grid of complex s (output of residual_spec).

For each s the truncated kernel is completed by a geometric tail (ratio of the last two unit time intervals), the
lumped high-degree columns are distributed onto their endpoint rows by linear interpolation in the mark, and the
eigenvalues of the resulting square matrix are computed.  A zero of det(I - K_s) is an eigenvalue equal to 1.
Printed: for every sigma = Re s the eigenvalue closest to 1 along the omega grid (omega = 0 included, except at s = mu), with a delete-one-batch
jackknife error of the leading eigenvalue.  Usage: residual_spec_eval.py FILE [--table]
"""
import sys, numpy as np
fn = sys.argv[1]; raw = np.fromfile(fn, dtype=np.float64)
T, nb, nf, nr, nc, nsig, nom = raw[:7]; nb, nf, nr, nc, nsig, nom = map(int, (nb, nf, nr, nc, nsig, nom)); p = 7
sig = raw[p:p + nsig]; p += nsig; om = raw[p:p + nom]; p += nom; ra = raw[p:p + nr]; p += nr
ns = nsig * nom; blk = ns * nr * nc * 2; per = 3 * blk + ns * nc * 2
done = (len(raw) - p) // per
def cplx(x, shape): return (x[0::2] + 1j * x[1::2]).reshape(shape)
K = np.zeros((done, ns, nr, nc), complex); KL = np.zeros_like(K); KP = np.zeros_like(K); Km = np.zeros((done, ns, nc), complex)
for b in range(done):
    o = p + b * per
    K[b] = cplx(raw[o:o + blk], (ns, nr, nc)); KL[b] = cplx(raw[o + blk:o + 2 * blk], (ns, nr, nc)); KP[b] = cplx(raw[o + 2 * blk:o + 3 * blk], (ns, nr, nc)); Km[b] = cplx(raw[o + 3 * blk:o + per], (ns, nc))
JX = 40; nex = JX - 1; ends = [(41, 60), (61, 100), (101, 200), (201, 500), (501, None)]
def square(Kmat, Kmean):
    n = len(ra); M = np.zeros((n, n), complex); M[:, :nex] = Kmat[:, :nex]
    for li, (lo, hi) in enumerate(ends):
        col = nex + li; rlo = nex + 2 * li; rhi = rlo + 1; alo = 1.0 / lo; ahi = 0.0 if hi is None else 1.0 / hi
        tot = Kmat[-1, col]; abar = Kmean[col] / tot if abs(tot) > 0 else 0.5 * (alo + ahi); w = (abar - ahi) / (alo - ahi)
        M[:, rlo] += w * Kmat[:, col]; M[:, rhi] += (1 - w) * Kmat[:, col]
    return M
def eig(sel, s):
    k = K[sel, s].mean(axis=0); kl = KL[sel, s].mean(axis=0); kp = KP[sel, s].mean(axis=0); q = kl.sum() / kp.sum()
    ev = np.linalg.eigvals(square(k + kl * q / (1 - q), Km[sel, s].mean(axis=0))); return ev[np.argsort(-np.abs(ev))], q
print(f"{fn}: T={T:g}, batches={done} x {nf} families, sigma grid {list(np.round(sig, 4))}, omega grid {om[0]:g}..{om[-1]:g} ({nom} values)")
allb = np.ones(done, bool); res = {}
for j in range(nsig):
    best = None
    for k in range(nom):
        s = j * nom + k; ev, q = eig(allb, s); d = np.abs(1 - ev); i = int(np.argmin(d)); res[(j, k)] = (ev[:3], q)
        if not (abs(om[k]) < 1e-12 and abs(sig[j] - 0.5657237) < 1e-6) and (best is None or d[i] < best[0]): best = (d[i], om[k], ev[i], ev[0], abs(q))
    lam0 = eig(allb, j * nom)[0][0]; loo = np.array([eig(np.arange(done) != b, j * nom)[0][0].real for b in range(done)])
    se = np.sqrt((done - 1) / done * ((loo - loo.mean())**2).sum()) if done > 1 else float('nan')
    print(f"sigma={sig[j]:.4f}: leading eigenvalue at omega=0: {lam0.real:.4f} +- {se:.4f}; second, third at omega=0: {abs(res[(j,0)][0][1]):.4f}, {abs(res[(j,0)][0][2]):.4f}; tail ratio |q| = {abs(res[(j,0)][1]):.3f};"
          f"  min over the omega grid (the point s = mu itself excluded) of |1 - eigenvalue| = {best[0]:.3f} at omega={best[1]:g} (eigenvalue {best[2].real:+.3f}{best[2].imag:+.3f}i, leading there {abs(best[3]):.3f})")
if "--table" in sys.argv:
    print("leading eigenvalue lambda_1(sigma + i omega) (re, im):")
    for k in range(nom):
        print(f"  omega={om[k]:6.2f}: " + "  ".join(f"{res[(j,k)][0][0].real:+.3f}{res[(j,k)][0][0].imag:+.3f}i" for j in range(nsig)))
