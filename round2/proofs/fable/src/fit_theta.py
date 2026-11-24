#!/usr/bin/env python3
"""Fit the GMM weight vector theta on a TRAINING sample and save it for the production runs (chunked, low memory).

Input: ftvec records written with diffout=1 (one differenced record per run: [dC, dS, da(p), db(p)]), or with two
checkpoints (then the difference last-minus-first is formed here).
theta = S^-1 bbar for the moment vectors d_r = da_r - mu0 db_r, with the plain equation (dS, dC) appended as the
last feature.  Features that are non-zero in fewer than --minfrac of the records get weight 0.
Also reports a K-fold cross-fitted estimate (theta from the other folds) with its delete-one-run jackknife error.
Output: OUT float64[p+1] in the feature layout of ftvec (last entry: weight of the plain equation, normalised to 1).
Usage: fit_theta.py PREFIX[,PREFIX2..] OUT [--minfrac F --ridge R --mu0 M --folds K --cols a:b,c:d --chunk N]
"""
import argparse, numpy as np
ap = argparse.ArgumentParser(); ap.add_argument("pre"); ap.add_argument("out")
ap.add_argument("--minfrac", type=float, default=0.01); ap.add_argument("--ridge", type=float, default=1e-6)
ap.add_argument("--mu0", type=float, default=0.565723); ap.add_argument("--folds", type=int, default=5)
ap.add_argument("--cols", default=""); ap.add_argument("--chunk", type=int, default=4000)
ap.add_argument("--drop", default="", help="comma-separated a:b or a:b:step index ranges of features forced to weight 0")
args = ap.parse_args()
srcs = []
for pre in args.pre.split(","):
    with open(pre + ".meta.txt") as f:
        p, ncp, nrec, runs, flag = map(int, f.readline().split())
    srcs.append((np.memmap(pre + ".rec.bin", dtype=np.float64, mode="r").reshape(nrec, ncp, 2 + 2 * p), ncp))
def chunks():
    off = 0
    for X, ncp in srcs:
        for i in range(0, len(X), args.chunk):
            Z = np.asarray(X[i:i + args.chunk, 0, :]) if ncp == 1 else np.asarray(X[i:i + args.chunk, -1, :]) - np.asarray(X[i:i + args.chunk, 0, :])
            A1 = np.concatenate([Z[:, 2:2 + p], Z[:, 1:2]], axis=1); B1 = np.concatenate([Z[:, 2 + p:], Z[:, 0:1]], axis=1)
            yield off + i, A1, B1
        off += len(X)
n = sum(len(X) for X, _ in srcs); K = args.folds
nz = np.zeros(p + 1); sabs = np.zeros(p + 1)
for off, A1, B1 in chunks(): nz += ((A1 != 0) | (B1 != 0)).sum(axis=0); sabs += np.abs(B1).sum(axis=0)
keep = nz >= args.minfrac * n; keep[-1] = True
if args.cols:
    m = np.zeros(p + 1, bool)
    for x in args.cols.split(","):
        lo, hi = map(int, x.split(":")); m[lo:hi] = True
    m[-1] = True; keep &= m
if args.drop:
    for x in args.drop.split(","):
        f = list(map(int, x.split(":"))); keep[f[0]:f[1]:(f[2] if len(f) > 2 else 1)] = False
    keep[-1] = True
q = int(keep.sum()); sc = sabs[keep] / n; sc[sc == 0] = 1
Mk = np.zeros((K, q, q)); bk = np.zeros((K, q)); ak = np.zeros((K, q))
for off, A1, B1 in chunks():
    A = A1[:, keep] / sc; B = B1[:, keep] / sc; D = A - args.mu0 * B; f = (off + np.arange(len(A))) % K
    for k in range(K):
        s = f == k
        if s.any(): Mk[k] += D[s].T @ D[s]; bk[k] += B[s].sum(axis=0); ak[k] += A[s].sum(axis=0)
def solve(M, b):
    M = M + args.ridge * np.diag(np.diag(M)) + 1e-300 * np.eye(q)
    return np.linalg.solve(M, b)
ths = [solve(Mk.sum(axis=0) - Mk[k], bk.sum(axis=0) - bk[k]) for k in range(K)]
th_all = solve(Mk.sum(axis=0), bk.sum(axis=0))
num = np.zeros(n); den = np.zeros(n); pn = np.zeros(n); pd = np.zeros(n)
for off, A1, B1 in chunks():
    A = A1[:, keep] / sc; B = B1[:, keep] / sc; idx = off + np.arange(len(A)); f = idx % K
    for k in range(K):
        s = f == k
        if s.any(): num[idx[s]] = A[s] @ ths[k]; den[idx[s]] = B[s] @ ths[k]
    pn[idx] = A1[:, -1]; pd[idx] = B1[:, -1]
def jack(a, b):
    est = a.sum() / b.sum(); loo = (a.sum() - a) / (b.sum() - b)
    return est, np.sqrt((len(a) - 1) / len(a) * ((loo - loo.mean())**2).sum())
ce, cse = jack(num, den); pe, pse = jack(pn, pd)
theta = np.zeros(p + 1); theta[keep] = th_all / sc; theta /= theta[-1]
theta.tofile(args.out)
print(f"training runs={n} layout p={p} equations kept q={q} (minfrac={args.minfrac})")
print(f"plain (differenced)      : {pe:.9f} +- {pse:.2e}")
print(f"cross-fitted GMM ({K} folds): {ce:.9f} +- {cse:.2e}   se ratio vs plain = {pse/cse:.2f}")
print(f"theta (all runs) written to {args.out}")
