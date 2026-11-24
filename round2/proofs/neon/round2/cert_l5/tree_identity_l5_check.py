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
import argparse, sys, struct, mmap, random
from fractions import Fraction
from functools import cache
from collections import Counter
from math import comb

MAGIC = 0x4E454F4E42334648

@cache
def downward(shapes, shape, k, Dc):
    result = 1
    for c,m in sorted(Counter(shapes[shape]).items()):
        if k < m: return 0
        result *= comb(k,m)*downward(shapes,c,Dc-2,Dc)**m
        k -= m
    return result

def charge_shape(C, i, k, h):
    shapes = C['l5_shapes']
    s = shapes[i]
    return downward(shapes,i,k,C['Dc']) if len(s) > 1 else h*downward(shapes,s[0],k,C['Dc'])

def l5_layout(lut, Dc, buckets):
    shapes, sizes, states, lengths, offsets = [], [], {}, {}, {}

    def partitions(first, left):
        if not left:
            yield ()
        else:
            for i, size in enumerate(sizes[:], 0):
                if i >= first and size <= left:
                    for suffix in partitions(i, left-size):
                        yield (i, *suffix)

    rows = 0
    for n in range(1, 6):
        for children in list(partitions(0, n-1)):
            i = len(shapes)
            shapes.append(children)
            sizes.append(n)
            for h in (0, 1):
                domain = [g for (hh,k,L),g in lut.items() if hh == h and k >= len(children) and h+k+L < Dc]
                if len(children) == 1-h:
                    domain += [-1-h-2*b for b in range(buckets)]
                states[i,h] = {g:j for j,g in enumerate(domain)}
                size = len(domain)
                for c in sorted(set(children)):
                    m = children.count(c)
                    size *= comb(lengths[c,1]+m-1,m)
                lengths[i,h] = size
                if n == 5:
                    offsets[i,h] = rows
                    rows += size
    index = {children:i for i,children in enumerate(shapes)}

    @cache
    def rank(g, h, children):
        shape = tuple(i for i,r in children)
        i = index[shape]
        r = states[i,h][g]
        for c in sorted(set(shape)):
            keys = [v for j,v in children if j == c]
            r = r*comb(lengths[c,1]+len(keys)-1,len(keys)) + sum(comb(v+j,j+1) for j,v in enumerate(keys))
        return i,r

    return rows, offsets, rank, tuple(shapes), states, lengths

def load_cert(path):
    if sys.byteorder != "little":
        raise ValueError("LP64 little-endian host required")
    with open(path, "rb") as fp:
        hdr_raw = fp.read(96)
        (magic, mode, version, P, Q, S, S_r,
         LF, KF, D0, KT, DP4, KP4, DYT, KYT, DYB, KYB, DS3, KS3) = struct.unpack("<Qiiqqqqiiiiiiiiiiii", hdr_raw)
        assert magic == MAGIC, f"Bad magic: {hex(magic)}"
        assert mode in (1, 2) and version == 2
        assert 0 < P <= 10**12 and 0 < Q <= 10**12 and 0 < S <= 10**14 and 0 < S_r <= 10**12
        assert 4 <= D0 <= 40 and D0-1 <= KT <= 64 and KT < KF <= 500 and D0 <= LF <= 30000
        Dc, buckets, mask, reserved, nL5 = struct.unpack("<IIIIQ",fp.read(24))
        assert 4 <= Dc <= D0 and buckets in (1,2) and (buckets == 1 or Dc < D0)
        assert mask in (1,2,3) and reserved == 0
        assert all(Dc <= d <= D0 and max(kmin,Dc-1) <= k <= KT
                   for d,k,kmin in ((DP4,KP4,1),(DYT,KYT,2),(DYB,KYB,2),(DS3,KS3,3)))

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

        rows, l5_offset, l5_rank, l5_shapes, l5_states, l5_lengths = l5_layout(lut,Dc,buckets)
        assert nL5 == rows and rows > 0
        lengths = (INNER_N,INNER_N,NT,NT,NT,NT,yp_len,yw_len,NV,g_len,hv_len,hc_len,
                   p4_len,yt_len,yb_len,s3_len,rows)
        fp.seek(0,2)
        assert fp.tell() == 120+8*sum(lengths), "invalid file length"
        data = memoryview(mmap.mmap(fp.fileno(),0,access=mmap.ACCESS_READ))[120:].cast("q")
        pos = 0

        def rd(n):
            nonlocal pos
            a = data[pos:pos+n]
            pos += n
            return a

        rpar = rd(INNER_N); rch = rd(INNER_N); trpar = rd(NT); trch = rd(NT)
        c2_f = rd(NT); c3_f = rd(NT); Yp_ref = rd(yp_len); Yw_ref = rd(yw_len)
        f = rd(NV); g = rd(g_len); hV = rd(hv_len); hC = rd(hc_len)
        qP4 = rd(p4_len); qYt = rd(yt_len); qYb = rd(yb_len); qS3 = rd(s3_len)
        q5 = rd(rows)

    return locals()

def eval_Sphi(C, par, motifs5=None):
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
    # Enumerate by extending the frontier, then deduplicate vertex subsets.
    # This does not depend on C's rank decoder or L5 row/source equations.
    def subsets(u):
        frontier = {frozenset((u,))}
        for _ in range(4):
            frontier = {s | {v} for s in frontier for p in s for v in ck[p] if v not in s}
        return frontier

    def rank_motif(u, subset, root):
        children = tuple(v for v in ck[u] if v in subset)
        h = hv[u]
        d = deg[u]
        if d >= C["Dc"]:
            if len(children) + (u != root) != 1:
                raise UnsupportedMotif
            b = int(C["buckets"] == 2 and d >= C["D0"])
            h = int(u != root)
            g = -1-h-2*b
        else:
            g = C["lut"][h,kv[u],Lv[u]]
        return C["l5_rank"](g,h,tuple(sorted(rank_motif(v,subset,root) for v in children)))

    counts = Counter()
    for u in range(n):
        if not isC[u]:
            continue
        for s in subsets(u):
            try:
                i,r = rank_motif(u,s,u)
            except UnsupportedMotif:
                continue
            h = hv[u] if deg[u] < C["Dc"] else 0
            value = C["q5"][C["l5_offset"][i,h]+r]
            ap_children = [v for v in ck[u] if v in s]
            charge = u if len(ap_children) > 1 else ap_children[0]
            counts[charge,i] += 1
            if motifs5 is not None:
                motifs5[s] = value
            tot += value
    for (u,i),count in counts.items():
        cap = charge_shape(C,i,kv[u],hv[u])
        assert count <= cap, 'per-shape charge count'
    for u in {u for u,i in counts}:
        ncap = sum(charge_shape(C,i,kv[u],hv[u]) for i in range(8,17))
        assert ncap <= deg[u]*comb(deg[u]-1+4*(C['Dc']-2),3)//4, 'total charge count'
    return tot


class UnsupportedMotif(Exception):
    pass

def comparison(C):
    """Exact coefficient reduction with the v1 charging and L5 envelope."""
    def pair(vals, include_zero=False):
        lo = hi = 0 if include_zero else next(vals)
        for v in vals:
            if v < lo: lo = v
            if v > hi: hi = v
        return lo,hi

    def join(a,b): return min(a[0],b[0]),max(a[1],b[1])
    def add(*pairs): return tuple(map(sum,zip(*pairs)))
    def mul(p,n): return p[0]*n,p[1]*n

    n, N, D0, KT = (C[k] for k in ('nchs','NG','D0','KT'))
    P = [u for u,j in enumerate(C['p_to_pi']) if j >= 0]
    W = [u for u,j in enumerate(C['w_to_wi']) if j >= 0]
    g,hv,hc,p4,yt,yb,s3 = (C[k] for k in ('g','hV','hC','qP4','qYt','qYb','qS3'))
    child = []
    for i,u in enumerate(W):
        chg = pair((g[p*N+u] for p in P))
        phv = pair(iter(hv[i::n]), True)
        pp4 = pair(iter(p4[i::n]), True)
        width = C['NYBW']*n*n
        pyb = (0,0)
        for j,(_,k,_) in enumerate(C['ybu_lut']):
            pyb = join(pyb,mul(pair(iter(yb[j*width+i:(j+1)*width:n]),True),k-1))
        child.append(add(chg,phv,pp4,pyb))

    # For tail children the actual Bernstein weights multiply all four
    # charges together, so taking joint extrema of their summed envelopes
    # is valid and can be tighter than the verifier's separate extrema.
    mtail = []
    for k in range(KT+2):
        ids = (C['INF_ID'], C['t1'](k,1), C['t2'](k,1), C['t3'](k,1))
        vals = [child[C['w_to_wi'][u]] for u in ids]
        mtail.append((min(lo for lo,hi in vals),max(hi for lo,hi in vals)))

    def blocks(arr, nrows, step, repeats=1):
        from itertools import chain
        return [pair(chain.from_iterable(arr[(r*nrows+i)*step:(r*nrows+i+1)*step] for r in range(repeats)))
                for i in range(nrows)]

    ahc = blocks(hc,C['NCU'],comb(n+1,2))
    ayt = blocks(yt,C['NYT'],comb(n+1,2),C['npars'])
    as3 = blocks(s3,C['NS3'],comb(n+2,3))
    l5 = {}
    states, lens, arr = C['l5_states'], C['l5_lengths'], C['q5']
    for (i,h),off in C['l5_offset'].items():
        shape = C['l5_shapes'][i]
        domain = states[i,h]
        if not domain: continue
        root_span = lens[i,h]//len(domain)
        if len(shape) > 1:
            for g,r in domain.items():
                l5[g,i] = pair(iter(arr[off+r*root_span:off+(r+1)*root_span]),True)
        else:
            child_domain = states[shape[0],1]
            if not child_domain: continue
            branch_span = root_span//len(child_domain)
            for g,r in child_domain.items():
                for j in range(len(domain)):
                    start = off+j*root_span+r*branch_span
                    l5[g,i] = join(l5.get((g,i),(0,0)),pair(iter(arr[start:start+branch_span]),True))
    f = C['f']
    def fidx(h,k,L): return (0 if h else (C['KF']+1)*C['NLF'])+k*C['NLF']+L
    lo4, hi = None,None
    def record(h,k,lo,up):
        nonlocal lo4,hi
        lo = 4*lo+(k-h)*C['S']
        lo4 = lo if lo4 is None else min(lo4,lo)
        hi = up if hi is None else max(hi,up)

    envelopes = {g:add(*(mul(l5.get((g,i),(0,0)),charge_shape(C,i,k,h)) for i in range(8,17)))
                 for (h,k,L),g in C['lut'].items() if 2 <= h+k+L < C['Dc']}
    for h in (0,1):
        for k in range(C['KF']+1):
            for L in range(C['LF']+1):
                d = L+k+h
                if d < 2: continue
                v = f[fidx(h,k,L)]
                charges = [(v,v)]
                if h:
                    charges.append(child[C['w_to_wi'][C['lut'][h,k,L]]] if d < D0 else mtail[min(k,KT+1)])
                if 2 <= k <= KT:
                    if d < D0: ch = ahc[C['cu_of_g'][C['lut'][h,k,L]]]
                    else:
                        ch = (0,0)
                        for t in ('t1','t2','t3'):
                            ch = join(ch,ahc[C['cu_of_g'][C[t](k,h)]])
                    charges.append(mul(ch,comb(k,2)))
                if h and (k,L) in C['yt_lut']:
                    charges.append(mul(ayt[C['yt_lut'][k,L]],comb(k,2)))
                if (h,k,L) in C['s3_lut']:
                    charges.append(mul(as3[C['s3_lut'][h,k,L]],comb(k,3)))
                if d < C['Dc']: charges.append(envelopes[C['lut'][h,k,L]])
                record(h,k,*add(*charges))
        v = f[C['VSTAR1'] if h else C['VSTAR0']]
        record(h,C['KF']+1,*add((v,v),mtail[KT+1] if h else (0,0)))
    return Fraction(lo4,4),Fraction(hi)

def check_tree(C, par, bounds, plant=None):
    n = len(par)
    deg = [0] * n
    nbrs = [[] for _ in range(n)]
    for v in range(1, n):
        p = par[v]; deg[v] += 1; deg[p] += 1; nbrs[v].append(p); nbrs[p].append(v)
    old = {}
    phi0 = eval_Sphi(C, par, old)
    lphi = Fraction(0)
    for w in range(n):
        rate = sum(Fraction(1, deg[v]) for v in nbrs[w])
        new = {}
        phi1 = eval_Sphi(C, par + [w], new)
        lphi += rate * (phi1 - phi0)
        if deg[w] == 1 and par[w] >= 0:
            # Promotion creates motifs with the former leaf w, not the
            # just-attached physical leaf. Enumerate from old 4-core sets
            # independently of the forward enumeration on the new tree.
            u = par[w]
            core = {x for x in range(n) if deg[x] >= 2}
            sets = {frozenset((u,))}
            for _ in range(3):
                sets = {s | {v} for s in sets for x in s for v in nbrs[x] if v in core-s}
            predicted = {}
            if deg[u] < C['Dc'] and plant != 'newborn':
                for s in sets:
                    t = s | {w}
                    if all(deg[x] < C['Dc'] or len(t.intersection(nbrs[x])) == 1 for x in s):
                        predicted[t] = new[t]
            assert predicted == {s:v for s,v in new.items() if w in s}, 'newborn ownership'
        if C['buckets'] == 2 and deg[w] == C['D0']-1:
            # Compare the actual shift at the D0-1 crossing with the
            # directional corner, using gain rate separated from promotions.
            for s,v in old.items():
                if w not in s:
                    continue
                delta = new[s]-v
                cap = (0 if plant == 'crossing' else C['D0']-1)*delta
                assert min(0,cap) <= rate*delta <= max(0,cap), 'boundary crossing'
    lam = Fraction(C["P"], C["Q"])
    diff = lphi - lam * phi0
    ok = (diff <= 0) if C["mode"] == 1 else (diff >= 0)
    core = sum(d >= 2 for d in deg)
    ok &= bounds[0]*core <= phi0 <= bounds[1]*core
    return ok, diff

def main():
    ap = argparse.ArgumentParser(description="Independent exact concrete-tree replay of a v2 certificate")
    ap.add_argument("certificate")
    ap.add_argument("ntrees",type=int,nargs="?",default=30)
    ap.add_argument("--scaled-bounds",type=int,nargs=2,metavar=("4Sc0","SC0"))
    ap.add_argument("--plant",choices=("crossing","newborn"),help=argparse.SUPPRESS)
    args = ap.parse_args()
    C = load_cert(args.certificate)
    ntrees = args.ntrees
    if args.scaled_bounds:
        bounds = (Fraction(args.scaled_bounds[0],4), Fraction(args.scaled_bounds[1]))
    else:
        bounds = comparison(C)
        print(f"Python comparison: 4Sc0={4*bounds[0]} SC0={bounds[1]}")
    assert bounds[1] > 0 and (C['mode'] == 2 or bounds[0] > 0)
    assert eval_Sphi(C,[-1,0,0]) > 0
    D0, KT = C["D0"], C["KT"]
    trees = [
        ("T_3", [-1, 0, 0]),
        ("star_D0-1", [-1] + [0] * (D0 - 1)),
        ("star_D0",   [-1] + [0] * D0),
        ("star_D0+2", [-1] + [0] * (D0 + 2)),
        ("hub_KT+2",  [-1] + [0 if i % 2 == 0 else i for i in range(2 * (KT + 2))]),
        ("deep_mixed", [-1, 0, 0, 1, 1, 3, 3, 5, 5, 2, 2, 1, 1, 0, 0, 3, 3, 1, 1, 0, 0]),
        ("split_path", [-1,0,1,2,3,4] + [0]*(D0-2) + [5]*(D0-1) + [2]),
        ("split_descendant", [-1,0,1,2,3,4] + [0,0,2] + [5]*(D0-2)),
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
        ok, rel = check_tree(C, par, bounds, args.plant)
        if not ok:
            print(f"  FAIL on {name}: scaled drift residual = {rel}")
            fails += 1
    print(f"Python exact-Fraction tree check ({args.certificate}): {len(trees)} trees checked, fails = {fails}")
    sys.exit(0 if fails == 0 else 1)

if __name__ == "__main__":
    main()
