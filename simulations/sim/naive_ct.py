#!/usr/bin/env python3
"""Naive reference simulator (full tree, discrete steps, no lazy updates) used to validate friend_ct.cpp.

The continuous-time process observed at time t is the discrete tree T_N with N = N(t) independent of the
jump chain, N(t) - 3 ~ NegBin(3, exp(-t)) (Yule process started from 3). So

    E[C(t)] = sum_n P(N(t) = n) E[C_n],     E[S(t)] likewise.

This script simulates the discrete model directly (adjacency lists, uniform V, uniform neighbour W),
records C_n and S_n for all n <= NMAX, and prints the mixtures for the requested times, with standard
errors.  Usage: naive_ct.py RUNS NMAX SEED t1 t2 ...
"""
import sys, random, math

def run(nmax, rng):
    # T_3: path 1 - 0 - 2
    adj = [[1, 2], [0], [0]]
    Cs = [0.0] * (nmax + 1); Ss = [0.0] * (nmax + 1)
    # maintain C and S incrementally
    deg = [2, 1, 1]
    leafcnt = [2, 0, 0]          # number of leaf neighbours
    C = 1; S = 1.0               # S = sum over core of L/D
    n = 3
    Cs[3] = C; Ss[3] = S
    while n < nmax:
        v = rng.randrange(n)
        w = adj[v][rng.randrange(deg[v])]
        # attach new vertex x to w
        x = n
        if deg[w] == 1:
            # w is a leaf with unique neighbour u: conversion
            u = adj[w][0]
            S -= leafcnt[u] / deg[u]
            leafcnt[u] -= 1
            S += leafcnt[u] / deg[u]
            adj[w].append(x); deg[w] = 2; leafcnt[w] = 1
            C += 1; S += 0.5
        else:
            S -= leafcnt[w] / deg[w]
            adj[w].append(x); deg[w] += 1; leafcnt[w] += 1
            S += leafcnt[w] / deg[w]
        adj.append([w]); deg.append(1); leafcnt.append(0)
        n += 1
        Cs[n] = C; Ss[n] = S
    return Cs, Ss

def main():
    runs = int(sys.argv[1]); nmax = int(sys.argv[2]); seed = int(sys.argv[3])
    ts = [float(x) for x in sys.argv[4:]]
    rng = random.Random(seed)
    # weights P(N(t) = n) = C(n-1, 2) p^3 (1-p)^(n-3), p = exp(-t)
    W = []
    for t in ts:
        p = math.exp(-t); w = [0.0] * (nmax + 1)
        for n in range(3, nmax + 1):
            w[n] = math.exp(math.log((n - 1) * (n - 2) / 2) + 3 * math.log(p) + (n - 3) * math.log1p(-p))
        W.append(w)
        print(f"# t={t}: truncated mass beyond NMAX = {1 - sum(w):.3e}")
    acc = [[0.0, 0.0, 0.0, 0.0] for _ in ts]
    for r in range(runs):
        Cs, Ss = run(nmax, rng)
        for i, w in enumerate(W):
            c = sum(w[n] * Cs[n] for n in range(3, nmax + 1))
            s = sum(w[n] * Ss[n] for n in range(3, nmax + 1))
            acc[i][0] += c; acc[i][1] += c * c; acc[i][2] += s; acc[i][3] += s * s
    for i, t in enumerate(ts):
        mc = acc[i][0] / runs; ms = acc[i][2] / runs
        ec = math.sqrt(max(acc[i][1] / runs - mc * mc, 0) / runs)
        es = math.sqrt(max(acc[i][3] / runs - ms * ms, 0) / runs)
        print(f"t={t}  E[C(t)] = {mc:.5f} +- {ec:.5f}   E[S(t)] = {ms:.5f} +- {es:.5f}   (mixture over n <= {nmax})")

if __name__ == "__main__":
    main()
