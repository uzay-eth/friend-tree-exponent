#!/usr/bin/env python3
"""Reduce raw per-run output files of ftprod to per-chunk first and second moments (small enough to archive).

For every complete chunk and every pair of consecutive checkpoints (t_i -> t_{i+1}) each run contributes, for the
plain equation (k = 0: num = S, den = C) and for each weight vector k >= 1,
    e_k = [num_k(t_{i+1}) - num_k(t_i)] - MU_REF * [den_k(t_{i+1}) - den_k(t_i)],     d_k = den_k(t_{i+1}) - den_k(t_i).
A run is the stacked vector u = (e_0, d_0, e_1, d_1, ...), observables inside consecutive windows.  Stored per
chunk: the number of runs n, sum of u, sum of u u^T.  From these follow (final_analysis.py): every ratio estimate
R = MU_REF + sum e / sum d over any window that is a union of consecutive windows, its delete-one-run (linearised)
and delete-one-chunk jackknife standard errors, and every paired difference of two such estimates inside one
sample.  MU_REF is only a centring constant; it keeps the second moments free of cancellation.

Usage: moments_reduce.py OUT.npz SPEC [SPEC ...]
   SPEC = kind:glob:rows_per_chunk:mu_ref[:block[:minchunk[:maxchunk]]]
   rows_per_chunk = runs * checkpoints of a complete chunk file (incomplete files are skipped; 0 = any size);
   block > 0: every file is one long sample, cut into blocks of this many runs (each block is a "chunk");
   minchunk/maxchunk: use only files whose chunk number (the digits after "_c") lies in [minchunk, maxchunk].
"""
import sys, glob, os, re, json
import numpy as np

def reduce_array(a, mu_ref):
    """a: rows (run, t, C, S, num_1, den_1, ...) of one chunk -> (n, checkpoints, sum u, sum u u^T, run ids)"""
    ts = np.unique(a[:, 1]); ncp = len(ts); n = a.shape[0] // ncp
    assert n * ncp == a.shape[0]
    o = np.lexsort((a[:, 1], a[:, 0])); a = a[o].reshape(n, ncp, a.shape[1])
    assert np.all(a[:, :, 0] == a[:, :1, 0]) and np.all(a[:, :, 1] == ts[None, :])
    assert len(np.unique(a[:, 0, 0])) == n, "duplicate run ids inside a chunk"
    nth = (a.shape[2] - 4) // 2
    parts = []
    for i in range(ncp - 1):
        d = a[:, i + 1, :] - a[:, i, :]
        for k in range(nth + 1):
            cn, cd = (3, 2) if k == 0 else (4 + 2 * (k - 1), 5 + 2 * (k - 1))
            parts.append(d[:, cn] - mu_ref * d[:, cd]); parts.append(d[:, cd])
    u = np.stack(parts, axis=1)
    return n, ts, u.sum(axis=0), u.T @ u, a[:, 0, 0]

def main():
    out = sys.argv[1]; store = {}; meta = {}
    for spec in sys.argv[2:]:
        f = spec.split(":"); kind, pat, nrows, mu_ref = f[0], f[1], int(f[2]), float(f[3])
        block = int(f[4]) if len(f) > 4 and f[4] else 0
        cmin = int(f[5]) if len(f) > 5 and f[5] else -1; cmax = int(f[6]) if len(f) > 6 and f[6] else 10**9
        ns, sums, sqs, names, ts0, seen = [], [], [], [], None, set()
        for fn in sorted(glob.glob(pat)):
            m = re.search(r"_c(\d+)", os.path.basename(fn))
            if m and not (cmin <= int(m.group(1)) <= cmax): continue
            try: a = np.loadtxt(fn, ndmin=2)
            except Exception: continue
            if nrows and a.shape[0] != nrows: continue
            if a.shape[0] == 0: continue
            pieces = [a]
            if block:
                ncp = len(np.unique(a[:, 1])); nrun = a.shape[0] // ncp; o = np.lexsort((a[:, 1], a[:, 0])); a = a[o]
                pieces = [a[i * block * ncp:(i + 1) * block * ncp] for i in range(nrun // block)]
            for j, piece in enumerate(pieces):
                n, ts, s, q, ids = reduce_array(piece, mu_ref)
                if ts0 is None: ts0 = ts
                assert np.array_equal(ts, ts0), f"checkpoints differ in {fn}"
                key = set(ids.tolist()); assert not (key & seen), f"run ids repeat across chunks of {kind}: {fn}"
                seen |= key
                ns.append(n); sums.append(s); sqs.append(q); names.append(os.path.basename(fn) + (f"#{j}" if block else ""))
        if not ns: print(f"{kind}: no complete chunks for {pat}"); continue
        store[kind + "_n"] = np.array(ns); store[kind + "_sum"] = np.array(sums); store[kind + "_sq"] = np.array(sqs)
        nth = (len(sums[0]) // (2 * (len(ts0) - 1))) - 1
        meta[kind] = {"checkpoints": [float(x) for x in ts0], "nth": int(nth), "mu_ref": mu_ref, "files": names}
        print(f"{kind}: {len(ns)} chunks, {int(np.sum(ns))} runs, checkpoints {list(ts0)}, weight vectors {nth}")
    store["meta"] = np.array(json.dumps(meta))
    np.savez_compressed(out, **store)
    print(f"wrote {out} ({os.path.getsize(out)/1e6:.2f} MB)")

if __name__ == "__main__":
    main()
