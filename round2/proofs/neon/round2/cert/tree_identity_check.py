#!/usr/bin/env python3
"""
Independent Python Exact-Rational (fractions.Fraction) Tree-Drift Checker for B3-FullH Certificates.
Reads the binary certificate from scratch (no shared C code!) and verifies:
  1. Header & exact byte length.
  2. Exact rational Phi(T) and brute-force generator drift L Phi(T) <= (P/Q) Phi(T)
     (mode=1) or >= (mode=2) by summing over all single-leaf attachments w in T on
     hand-crafted boundary trees (T_3, star D0-1, D0, D0+1, hub k=KT+2, deep P4/Ytop/Ybot/S3)
     plus seeded random friend trees.
"""
import sys, struct, array, random
from fractions import Fraction

MAGIC = 0x4E454F4E42334648

def load_cert(path):
    with open(path, "rb") as fp:
        hdr_raw = fp.read(96)
        (magic, mode, version, P, Q, S, S_r,
         LF, KF, D0, KT, DP4, KP4, DYT, KYT, DYB, KYB, DS3, KS3) = struct.unpack("<Qiiqqqqiiiiiiiiiiii", hdr_raw)
        assert magic == MAGIC, f"Bad magic: {hex(magic)}"
        assert mode in (1, 2) and version == 1
        assert P > 0 and Q > 0 and S > 0 and S_r > 0 and KT >= D0 - 1

        NLF, NKF = LF + 1, KF + 1
        VSTAR1 = 2 * NKF * NLF
        VSTAR0 = VSTAR1 + 1
        NV = VSTAR0 + 1
        NT = 2 * (KT + 2)

        lut = {}
        cnt = 0
        for h in (1, 0):
            for D in range(2, D0):
                for k in range(0, D - h + 1):
                    lut[(h, k, D - h - k)] = cnt
                    cnt += 1
        INNER_N = cnt
        INF_ID = cnt
        T1_BASE = cnt + 1
        T2_BASE = T1_BASE + NT
        T3_BASE = T2_BASE + NT
        NG = T3_BASE + NT

        def tidx(kap, h): return (0 if h else (KT + 2)) + kap
        def t1(k, h): return T1_BASE + tidx(KT + 1 if k > KT else k, h)
        def t2(k, h): return T2_BASE + tidx(KT + 1 if k > KT else k, h)
        def t3(k, h): return T3_BASE + tidx(KT + 1 if k > KT else k, h)

        kind = [0] * NG; K = [0] * NG; H = [0] * NG
        for (h, k, L), idx in lut.items():
            kind[idx] = 0; K[idx] = k; H[idx] = h
        kind[INF_ID] = 1; K[INF_ID] = -1; H[INF_ID] = -1
        for h in (1, 0):
            for kap in range(KT + 2):
                for b, kk in ((t1(kap, h), 2), (t2(kap, h), 3), (t3(kap, h), 4)):
                    kind[b] = kk; K[b] = kap; H[b] = h

        p_to_pi = [-1] * NG; w_to_wi = [-1] * NG; hu_of_g = [-1] * NG; cu_of_g = [-1] * NG
        npars = nchs = NHU = NCU = 0
        for u in range(NG):
            uk, k, h = kind[u], K[u], H[u]
            if uk == 1 or k >= 1: p_to_pi[u] = npars; npars += 1
            if uk == 1 or h == 1: w_to_wi[u] = nchs;  nchs += 1
            if (uk == 0 or uk >= 2) and h == 1 and 1 <= k <= KT: hu_of_g[u] = NHU; NHU += 1
            if (uk == 0 or uk >= 2) and 2 <= k <= KT:            cu_of_g[u] = NCU; NCU += 1

        p4_lut = {}; NP4 = 0
        for k in range(1, KP4 + 1):
            for L in range(0, DP4 - k - 1): p4_lut[(k, L)] = NP4; NP4 += 1
        yt_lut = {}; NYT = 0
        for k in range(2, KYT + 1):
            for L in range(0, DYT - k - 1): yt_lut[(k, L)] = NYT; NYT += 1
        ybu_lut = {}; NYBU = 0
        for h in (1, 0):
            for k in range(2, KYB + 1):
                for L in range(0, DYB - k - h): ybu_lut[(h, k, L)] = NYBU; NYBU += 1
        ybw_lut = {}; NYBW = 0
        for k in range(1, KYB + 1):
            for L in range(0, DYB - k - 1): ybw_lut[(k, L)] = NYBW; NYBW += 1
        s3_lut = {}; NS3 = 0
        for h in (1, 0):
            for k in range(3, KS3 + 1):
                for L in range(0, DS3 - k - h): s3_lut[(h, k, L)] = NS3; NS3 += 1

        yp_len = NT * nchs
        yw_len = npars * (KT + 2)
        g_len  = NG * NG
        hv_len = npars * NHU * nchs
        hc_len = NCU * nchs * (nchs + 1) // 2
        p4_len = npars * NP4 * NP4 * nchs
        yt_len = npars * NYT * nchs * (nchs + 1) // 2
        yb_len = NYBU * NYBW * nchs * nchs
        s3_len = NS3 * nchs * (nchs + 1) * (nchs + 2) // 6

        def rd(n):
            a = array.array("q")
            a.fromfile(fp, n)
            return a

        rpar = rd(INNER_N); rch = rd(INNER_N); trpar = rd(NT); trch = rd(NT)
        c2_f = rd(NT); c3_f = rd(NT); Yp_ref = rd(yp_len); Yw_ref = rd(yw_len)
        f = rd(NV); g = rd(g_len); hV = rd(hv_len); hC = rd(hc_len)
        qP4 = rd(p4_len); qYt = rd(yt_len); qYb = rd(yb_len); qS3 = rd(s3_len)
        assert fp.read(1) == b"", "Trailing bytes in certificate!"

    return locals()

def eval_Sphi(C, par):
    n = len(par)
    deg = [0] * n
    for v in range(1, n):
        deg[v] += 1; deg[par[v]] += 1
    isC = [deg[u] >= 2 for u in range(n)]
    hv = [1 if par[u] >= 0 else 0 for u in range(n)]
    Lv = [0] * n
    ck = [[] for _ in range(n)]
    for v in range(1, n):
        p = par[v]
        if deg[v] == 1: Lv[p] += 1
        else: ck[p].append(v)
    kv = [len(ck[u]) for u in range(n)]

    D0, KT, LF, KF, NLF, NKF, NG = C["D0"], C["KT"], C["LF"], C["KF"], C["NLF"], C["NKF"], C["NG"]
    lut, t1, t2, t3, INF_ID = C["lut"], C["t1"], C["t2"], C["t3"], C["INF_ID"]
    p_to_pi, w_to_wi, hu_of_g, cu_of_g = C["p_to_pi"], C["w_to_wi"], C["hu_of_g"], C["cu_of_g"]
    nchs, NHU, NP4, NYT, NYBW = C["nchs"], C["NHU"], C["NP4"], C["NYT"], C["NYBW"]
    f, g, hV, hC, qP4, qYt, qYb, qS3 = C["f"], C["g"], C["hV"], C["hC"], C["qP4"], C["qYt"], C["qYb"], C["qS3"]

    def vidx(L, k, h):
        if k > KF: return C["VSTAR1"] if h else C["VSTAR0"]
        if L > LF: L = LF
        return (0 if h else (KF + 1) * NLF) + k * NLF + L

    def embed(L, k, h):
        D = L + k + h
        if D < D0: return [(lut[(h, k, L)], Fraction(1))]
        t = Fraction(D0, D); s = 1 - t
        return [(INF_ID, s**3), (t1(k, h), 3 * t * s**2), (t2(k, h), 3 * t**2 * s), (t3(k, h), t**3)]

    def hv_off(p, ua, w): return (p_to_pi[p] * NHU + ua) * nchs + w_to_wi[w]
    def hc_off(ca, w1, w2):
        i1, i2 = w_to_wi[w1], w_to_wi[w2]
        if i1 > i2: i1, i2 = i2, i1
        pn = nchs * (nchs + 1) // 2
        return ca * pn + i1 * nchs - i1 * (i1 - 1) // 2 + (i2 - i1)
    def p4_off(a, b4, c4, d): return ((p_to_pi[a] * NP4 + b4) * NP4 + c4) * nchs + w_to_wi[d]
    def yt_off(p, u4, w1, w2):
        i1, i2 = w_to_wi[w1], w_to_wi[w2]
        if i1 > i2: i1, i2 = i2, i1
        pn = nchs * (nchs + 1) // 2
        return (p_to_pi[p] * NYT + u4) * pn + i1 * nchs - i1 * (i1 - 1) // 2 + (i2 - i1)
    def yb_off(u4, w2_4, w1, c): return ((u4 * NYBW + w2_4) * nchs + w_to_wi[w1]) * nchs + w_to_wi[c]
    def s3_off(u4, w1, w2, w3):
        a, b, c = sorted((w_to_wi[w1], w_to_wi[w2], w_to_wi[w3]))
        tn = nchs * (nchs + 1) * (nchs + 2) // 6
        na, nb = nchs - a, nchs - b
        return u4 * tn + (tn - na * (na + 1) * (na + 2) // 6) + (na * (na + 1) // 2 - nb * (nb + 1) // 2) + (c - b)

    em = [embed(Lv[u], kv[u], hv[u]) if isC[u] else [] for u in range(n)]
    tot = Fraction(0)
    for u in range(n):
        if isC[u]: tot += f[vidx(Lv[u], kv[u], hv[u])]
    for w in range(1, n):
        if isC[w]:
            p = par[w]
            for ip, wp in em[p]:
                for iw, ww in em[w]: tot += wp * ww * g[ip * NG + iw]
    for u in range(1, n):
        if isC[u] and 1 <= kv[u] <= KT:
            p = par[u]
            for w in ck[u]:
                for iu, wu in em[u]:
                    ua = hu_of_g[iu]
                    if ua < 0: continue
                    for ip, wp in em[p]:
                        for iw, ww in em[w]: tot += wu * wp * ww * hV[hv_off(ip, ua, iw)]
    for u in range(n):
        if isC[u] and 2 <= kv[u] <= KT:
            for j1 in range(len(ck[u])):
                for j2 in range(j1 + 1, len(ck[u])):
                    w1, w2 = ck[u][j1], ck[u][j2]
                    for iu, wu in em[u]:
                        ca = cu_of_g[iu]
                        if ca < 0: continue
                        for i1, wt1 in em[w1]:
                            for i2, wt2 in em[w2]: tot += wu * wt1 * wt2 * hC[hc_off(ca, i1, i2)]
    for c in range(1, n):
        if isC[c] and 1 <= kv[c] <= C["KP4"] and deg[c] < C["DP4"]:
            b = par[c]
            if hv[b] == 1 and 1 <= kv[b] <= C["KP4"] and deg[b] < C["DP4"]:
                b4, c4, a = C["p4_lut"][(kv[b], Lv[b])], C["p4_lut"][(kv[c], Lv[c])], par[b]
                for d in ck[c]:
                    for ia, wa in em[a]:
                        for id_, wd in em[d]: tot += wa * wd * qP4[p4_off(ia, b4, c4, id_)]
    for u in range(1, n):
        if isC[u] and 2 <= kv[u] <= C["KYT"] and deg[u] < C["DYT"]:
            u4, p = C["yt_lut"][(kv[u], Lv[u])], par[u]
            for j1 in range(len(ck[u])):
                for j2 in range(j1 + 1, len(ck[u])):
                    w1, w2 = ck[u][j1], ck[u][j2]
                    for ip, wp in em[p]:
                        for i1, wt1 in em[w1]:
                            for i2, wt2 in em[w2]: tot += wp * wt1 * wt2 * qYt[yt_off(ip, u4, i1, i2)]
    for u in range(n):
        if isC[u] and 2 <= kv[u] <= C["KYB"] and deg[u] < C["DYB"]:
            u4 = C["ybu_lut"][(hv[u], kv[u], Lv[u])]
            for w2 in ck[u]:
                if not (1 <= kv[w2] <= C["KYB"] and deg[w2] < C["DYB"]): continue
                w2_4 = C["ybw_lut"][(kv[w2], Lv[w2])]
                for w1 in ck[u]:
                    if w1 == w2: continue
                    for c in ck[w2]:
                        for i1, wt1 in em[w1]:
                            for ic, wtc in em[c]: tot += wt1 * wtc * qYb[yb_off(u4, w2_4, i1, ic)]
    for u in range(n):
        if isC[u] and 3 <= kv[u] <= C["KS3"] and deg[u] < C["DS3"]:
            u4 = C["s3_lut"][(hv[u], kv[u], Lv[u])]
            m = len(ck[u])
            for j1 in range(m):
                for j2 in range(j1 + 1, m):
                    for j3 in range(j2 + 1, m):
                        w1, w2, w3 = ck[u][j1], ck[u][j2], ck[u][j3]
                        for i1, wt1 in em[w1]:
                            for i2, wt2 in em[w2]:
                                for i3, wt3 in em[w3]: tot += wt1 * wt2 * wt3 * qS3[s3_off(u4, i1, i2, i3)]
    return tot

def check_tree(C, par):
    n = len(par)
    deg = [0] * n
    nbrs = [[] for _ in range(n)]
    for v in range(1, n):
        p = par[v]; deg[v] += 1; deg[p] += 1; nbrs[v].append(p); nbrs[p].append(v)
    phi0 = eval_Sphi(C, par)
    lphi = Fraction(0)
    for w in range(n):
        rate = sum(Fraction(1, deg[v]) for v in nbrs[w])
        phi1 = eval_Sphi(C, par + [w])
        lphi += rate * (phi1 - phi0)
    lam = Fraction(C["P"], C["Q"])
    diff = lphi - lam * phi0
    ok = (diff <= 0) if C["mode"] == 1 else (diff >= 0)
    return ok, float(diff / phi0) if phi0 != 0 else 0.0

def main():
    if len(sys.argv) < 2:
        print("Usage: tree_identity_check.py <cert.bin> [ntrees]"); sys.exit(2)
    C = load_cert(sys.argv[1])
    ntrees = int(sys.argv[2]) if len(sys.argv) > 2 else 15
    D0, KT = C["D0"], C["KT"]
    trees = [
        ("T_3", [-1, 0, 0]),
        ("star_D0-1", [-1] + [0] * (D0 - 1)),
        ("star_D0",   [-1] + [0] * D0),
        ("star_D0+2", [-1] + [0] * (D0 + 2)),
        ("hub_KT+2",  [-1] + [0 if i % 2 == 0 else i for i in range(2 * (KT + 2))]),
        ("deep_mixed", [-1, 0, 0, 1, 1, 3, 3, 5, 5, 2, 2, 1, 1, 0, 0, 3, 3, 1, 1, 0, 0]),
    ]
    rng = random.Random(20261124)
    for i in range(ntrees):
        n = 5 + (i * 4) % 36
        par = [-1, 0, 0]; adj = [[1, 2], [0], [0]]
        for j in range(3, n):
            V = rng.randrange(j); W = rng.choice(adj[V])
            par.append(W); adj[W].append(j); adj.append([W])
        trees.append((f"rand_{i}_n{n}", par))

    fails = 0
    for name, par in trees:
        ok, rel = check_tree(C, par)
        if not ok:
            print(f"  FAIL on {name}: rel_diff = {rel:+.6e}")
            fails += 1
    print(f"Python exact-Fraction tree check ({sys.argv[1]}): {len(trees)} trees checked, fails = {fails}")
    sys.exit(0 if fails == 0 else 1)

if __name__ == "__main__":
    main()
