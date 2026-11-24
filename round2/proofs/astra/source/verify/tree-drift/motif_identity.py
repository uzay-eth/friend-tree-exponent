#!/usr/bin/env python3
"""Independent exact check of the motif-bracket drift identity.

For each test tree T it computes, in exact Fractions:
  * LPhi(T) by brute force: sum over attachment targets W of wt(W) * (Phi(T+W) - Phi(T));
  * the local expressions V(u), E(p->w), P(p->u->w), Ch(u;{w1,w2}) of the
    written-out decomposition (Identity (I) in proofs/bracket-motifs.md);
and checks
  (a) LPhi(T) == sum of all local expressions          (exact identity),
  (b) every local expression <= lam * term (UPPER) or >= lam * term (LOWER),
  (c) Phi(T) >= c0*C(T) (UPPER) or <= C0*C(T) (LOWER).
Its own loader re-derives the layout from cert_format.h; it does not import repository code.
"""
import argparse
import random
import struct
import sys
from fractions import Fraction as Fr
from functools import lru_cache

from array import array


class Cert:
    def __init__(self, path):
        raw = open(path, "rb").read()
        v = struct.unpack("<Qiiiiqqiiiiiiiiii9ii", raw[:120])
        (self.magic, self.mode, self.P, self.Q, _, self.S, self.Sr, self.LF, self.KF, self.D0, self.KT,
         self.DHU, self.KHU, self.DCU, self.KCU, self.QK, nk) = v[:17]
        self.knots = [Fr(x, self.QK) for x in v[17:17 + nk]]
        self.lam = Fr(self.P, self.Q)
        D0, KT = self.D0, self.KT
        self.inner = {}
        for h in (1, 0):
            for D in range(2, D0):
                for k in range(0, D - h + 1):
                    self.inner[(D - h - k, k, h)] = len(self.inner)
        self.NI = len(self.inner)
        self.INF = self.NI
        self.NG = self.NI + 1 + 2 * (KT + 2)
        kinds = {}
        for s, i in self.inner.items():
            kinds[i] = ("I", s)
        kinds[self.INF] = ("INF", None)
        for h in (1, 0):
            for kap in range(KT + 2):
                kinds[self.t1(kap, h)] = ("T1", (kap, h))
        self.kinds = kinds
        isP = lambda i: kinds[i][0] == "INF" or (kinds[i][0] == "I" and kinds[i][1][1] >= 1) or (kinds[i][0] == "T1" and kinds[i][1][0] >= 1)
        isW = lambda i: kinds[i][0] == "INF" or (kinds[i][0] == "I" and kinds[i][1][2] == 1) or (kinds[i][0] == "T1" and kinds[i][1][1] == 1)
        self.pi = {}
        self.wi = {}
        for i in range(self.NG):
            if isP(i):
                self.pi[i] = len(self.pi)
            if isW(i):
                self.wi[i] = len(self.wi)
        self.hu = {}
        for k in range(1, self.KHU + 1):
            for L in range(self.DHU):
                if L + k + 1 < self.DHU:
                    self.hu[(L, k, 1)] = len(self.hu)
        self.cu = {}
        for h in (1, 0):
            for k in range(2, self.KCU + 1):
                for L in range(self.DCU):
                    if L + k + h < self.DCU:
                        self.cu[(L, k, h)] = len(self.cu)
        NV = 2 * (self.KF + 1) * (self.LF + 1) + 2
        nP, nW = len(self.pi), len(self.wi)
        off = 120
        def take(n):
            nonlocal off
            a = array("q")
            a.frombytes(raw[off:off + 8 * n])
            off += 8 * n
            return a
        NI = self.NI
        self.nk = nk
        self.rpar, self.rch, self.xpar, self.xch, self.xcc = (take(NI) for _ in range(5))
        self.kpc, self.kcp, self.kcc = (take(NI * nk) for _ in range(3))
        self.f = take(NV)
        self.g = take(self.NG * self.NG)
        self.nP, self.nHU = nP, len(self.hu)
        self.hV = take(nP * len(self.hu) * nW)
        self.pairn = nW * (nW + 1) // 2
        self.hC = take(len(self.cu) * self.pairn)
        assert off == len(raw), "trailing bytes"
        self.nW = nW

    def t1(self, kap, h):
        return self.NI + 1 + (0 if h else self.KT + 2) + kap

    # ----- state helpers: a state is (L, k, h) of a core vertex -----
    def mu(self, s):
        L, k, h = s
        D = L + k + h
        if D < self.D0:
            return ((self.inner[s], Fr(1)),)
        t = Fr(self.D0, D)
        return ((self.INF, 1 - t), (self.t1(min(k, self.KT + 1), h), t))

    def F(self, s):
        L, k, h = s
        if k > self.KF:
            i = 2 * (self.KF + 1) * (self.LF + 1) + (0 if h else 1)
        else:
            i = (0 if h else (self.KF + 1) * (self.LF + 1)) + k * (self.LF + 1) + min(L, self.LF)
        return Fr(self.f[i], self.S)

    def tail(self, s):
        return sum(s) >= self.D0

    def G(self, sp, sw):
        if self.tail(sp) and self.tail(sw):
            return Fr(0)
        return sum((a * b * Fr(self.g[i * self.NG + j], self.S) for i, a in self.mu(sp) for j, b in self.mu(sw)), Fr(0))

    def HV(self, sp, su, sw):
        if su not in self.hu:
            return Fr(0)
        ua = self.hu[su]
        return sum((a * b * Fr(self.hV[(self.pi[i] * self.nHU + ua) * self.nW + self.wi[j]], self.S)
                    for i, a in self.mu(sp) for j, b in self.mu(sw)), Fr(0))

    def HC(self, su, s1, s2):
        if su not in self.cu:
            return Fr(0)
        ca = self.cu[su]
        tot = Fr(0)
        for i, a in self.mu(s1):
            for j, b in self.mu(s2):
                x, y = sorted((self.wi[i], self.wi[j]))
                tot += a * b * Fr(self.hC[ca * self.pairn + x * self.nW - x * (x - 1) // 2 + (y - x)], self.S)
        return tot

    def A(self, table, s, x):
        i0 = self.inner[s] * self.nk
        a = [Fr(table[i0 + j], self.S) for j in range(self.nk)]
        xs = self.knots
        for i in range(len(xs) - 1):
            if x <= xs[i + 1]:
                al = (x - xs[i]) / (xs[i + 1] - xs[i])
                return (1 - al) * a[i] + al * a[i + 1]
        raise ValueError(x)

    def isinner(self, s):
        return sum(s) < self.D0


def tree_states(par):
    n = len(par)
    deg = [0] * n
    for v in range(1, n):
        deg[v] += 1
        deg[par[v]] += 1
    core = [v for v in range(n) if deg[v] >= 2]
    L = {v: 0 for v in core}
    ch = {v: [] for v in core}
    for v in range(1, n):
        if deg[v] == 1:
            L[par[v]] += 1
        else:
            ch[par[v]].append(v)
    st = {v: (L[v], len(ch[v]), 0 if par[v] == -1 else 1) for v in core}
    return deg, core, ch, st


def phi(C, par):
    deg, core, ch, st = tree_states(par)
    tot = sum((C.F(st[v]) for v in core), Fr(0))
    for w in core:
        if par[w] == -1:
            continue
        p = par[w]
        tot += C.G(st[p], st[w])
        if par[p] != -1:
            tot += C.HV(st[par[p]], st[p], st[w])
    for u in core:
        c = ch[u]
        for i in range(len(c)):
            for j in range(i + 1, len(c)):
                tot += C.HC(st[u], st[c[i]], st[c[j]])
    return tot, len(core)


def brute_drift(C, par):
    n = len(par)
    deg = [0] * n
    adj = [[] for _ in range(n)]
    for v in range(1, n):
        deg[v] += 1
        deg[par[v]] += 1
        adj[v].append(par[v])
        adj[par[v]].append(v)
    p0, _ = phi(C, par)
    d = Fr(0)
    for W in range(n):
        wt = sum((Fr(1, deg[V]) for V in adj[W]), Fr(0))
        p1, _ = phi(C, par + [W])
        d += wt * (p1 - p0)
    return d, p0


NB = (1, 0, 1)


def decomposition(C, par):
    """Return list of (kind, local_expression, term_value)."""
    deg, core, ch, st = tree_states(par)
    inv = {v: Fr(1, deg[v]) for v in core}
    up = lambda s: (s[0] + 1, s[1], s[2])
    bud = lambda s: (s[0] - 1, s[1] + 1, s[2])
    rho = {}
    for v in core:
        r = sum((inv[c] for c in ch[v]), Fr(0))
        if par[v] != -1:
            r += inv[par[v]]
        rho[v] = r
    rate = {v: st[v][0] + rho[v] for v in core}
    b = {v: Fr(st[v][0], deg[v]) for v in core}
    H = lambda s: s in C.hu
    Cu = lambda s: s in C.cu
    out = []

    def jumps(fn, verts, x):
        """Leaf-gain and budding jumps of fn(states) at vertex x (index into verts)."""
        s = [st[v] for v in verts]
        base = fn(*s)
        sL = list(s); sL[x] = up(s[x])
        dL = fn(*sL) - base
        if s[x][0] >= 1:
            sB = list(s); sB[x] = bud(s[x])
            dB = fn(*sB) - base
        else:
            dB = Fr(0)
        return dL, dB

    # vertex terms
    for u in core:
        s = st[u]
        dLf = C.F(up(s)) - C.F(s)
        if C.isinner(s):
            kept = s[2] * Fr(int(C.rpar[C.inner[s]]), C.Sr) + s[1] * Fr(int(C.rch[C.inner[s]]), C.Sr)
        else:
            kept = rho[u]
        V = (s[0] + kept) * dLf
        if s[0] >= 1:
            V += b[u] * (C.F(bud(s)) - C.F(s) + C.F(NB) + C.G(bud(s), NB))
        out.append(("V", V, C.F(s)))

    # edge terms
    for w in core:
        if par[w] == -1:
            continue
        p = par[w]
        sp, sw = st[p], st[w]
        E = Fr(0)
        gfun = lambda a, c: C.G(a, c)
        dLp, dBp = jumps(gfun, [p, w], 0)
        dLw, dBw = jumps(gfun, [p, w], 1)
        # vertex-to-edge transfers
        if C.isinner(sp):
            E += (inv[w] - Fr(int(C.rch[C.inner[sp]]), C.Sr)) * (C.F(up(sp)) - C.F(sp))
        if C.isinner(sw):
            E += (inv[p] - Fr(int(C.rpar[C.inner[sw]]), C.Sr)) * (C.F(up(sw)) - C.F(sw))
        # p side
        if C.isinner(sp):
            own = sp[0] + inv[w]
            if H(sp):
                own += Fr(int(C.xpar[C.inner[sp]]), C.Sr)
            if Cu(sp):
                own += (sp[1] - 1) * Fr(int(C.xcc[C.inner[sp]]), C.Sr)
            E += own * dLp
            if sp[2] == 1 and not H(sp):
                x = inv[par[p]]
                E += x * dLp - C.A(C.kpc, sp, x)
            if sp[1] >= 2 and not Cu(sp):
                for c in ch[p]:
                    if c != w:
                        E += inv[c] * dLp - C.A(C.kcc, sp, inv[c])
        else:
            E += rate[p] * dLp
        if sp[0] >= 1:
            E += b[p] * dBp
            if Cu(bud(sp)):
                E += b[p] * C.HC(bud(sp), sw, NB)
        # w side
        if C.isinner(sw):
            own = sw[0] + inv[p]
            if H(sw):
                own += sw[1] * Fr(int(C.xch[C.inner[sw]]), C.Sr)
            E += own * dLw
            if sw[1] >= 1 and not H(sw):
                for c in ch[w]:
                    E += inv[c] * dLw - C.A(C.kcp, sw, inv[c])
        else:
            E += rate[w] * dLw
        if sw[0] >= 1:
            E += b[w] * dBw
            if H(bud(sw)):
                E += b[w] * C.HV(sp, bud(sw), NB)
        # A-pieces received from neighbouring motifs
        if C.isinner(sw) and sw[2] == 1 and sw[1] >= 1 and not H(sw):
            E += sw[1] * C.A(C.kpc, sw, inv[p])
        if C.isinner(sp) and sp[2] == 1 and not H(sp):
            E += C.A(C.kcp, sp, inv[w])
        if C.isinner(sp) and sp[1] >= 2 and not Cu(sp):
            E += (sp[1] - 1) * C.A(C.kcc, sp, inv[w])
        out.append(("E", E, C.G(sp, sw)))

    # vertical 2-paths
    for w in core:
        if par[w] == -1 or par[par[w]] == -1:
            continue
        u = par[w]; p = par[u]
        if not H(st[u]):
            continue
        su = st[u]
        Pv = (inv[p] - Fr(int(C.xpar[C.inner[su]]), C.Sr)) * (C.G(up(su), st[w]) - C.G(su, st[w])) \
           + (inv[w] - Fr(int(C.xch[C.inner[su]]), C.Sr)) * (C.G(st[p], up(su)) - C.G(st[p], su))
        verts = [p, u, w]
        for x in range(3):
            dL, dB = jumps(C.HV, verts, x)
            Pv += rate[verts[x]] * dL + b[verts[x]] * dB
        out.append(("P", Pv, C.HV(st[p], su, st[w])))

    # cherries
    for u in core:
        su = st[u]
        if not Cu(su):
            continue
        c = ch[u]
        xc = Fr(int(C.xcc[C.inner[su]]), C.Sr)
        for i in range(len(c)):
            for j in range(i + 1, len(c)):
                w1, w2 = c[i], c[j]
                Chv = (inv[w1] - xc) * (C.G(up(su), st[w2]) - C.G(su, st[w2])) \
                    + (inv[w2] - xc) * (C.G(up(su), st[w1]) - C.G(su, st[w1]))
                verts = [u, w1, w2]
                for x in range(3):
                    dL, dB = jumps(C.HC, verts, x)
                    Chv += rate[verts[x]] * dL + b[verts[x]] * dB
                out.append(("Ch", Chv, C.HC(su, st[w1], st[w2])))
    return out


# ---------------- test trees ----------------
def friend_tree(n, rng):
    par = [-1, 0, 0]
    adj = [[1, 2], [0], [0]]
    while len(par) < n:
        V = rng.randrange(len(par))
        W = rng.choice(adj[V])
        par.append(W)
        adj.append([W])
        adj[W].append(len(par) - 1)
    return par


def build(spec):
    """spec: nested list; a node is (num_leaves, [children nodes]). Root first."""
    par = []
    def rec(node, p):
        idx = len(par)
        par.append(p)
        nl, kids = node
        for _ in range(nl):
            par.append(idx)
        for kd in kids:
            rec(kd, idx)
    rec(spec, -1)
    return par


def random_spec(rng, depth, C):
    """Random core shape stressing H_u, C_u, D0 boundary, K_T clamp, tails."""
    choices_L = [0, 1, 2, 3, 5, 8, 10, 13, 14, 15, 20, 21, 22, 23, 24, 25, 30, 40]
    nkid = 0 if depth == 0 else rng.choice([0, 0, 1, 1, 2, 2, 3, 4, 6, 7, 13, 14])
    kids = [random_spec(rng, depth - 1, C) for _ in range(nkid)]
    L = rng.choice(choices_L)
    return (L, kids)


def fix_spec(spec, is_root=True):
    L, kids = spec
    kids = [fix_spec(k, False) for k in kids]
    D = L + len(kids) + (0 if is_root else 1)
    if D < 2:
        L += 2 - D
    return (L, kids)


def sample_trees(up, seed):
    rng = random.Random(seed)
    trees = [("T3", [-1, 0, 0]), ("K13", [-1, 0, 0, 0]), ("P4", [-1, 0, 0, 1])]
    for i in range(20):
        trees.append((f"friend{i}", friend_tree(rng.randrange(5, 160), rng)))
    made = 0
    while made < 20:
        par = build(fix_spec(random_spec(rng, rng.randrange(1, 4), up)))
        if len(par) <= 260:
            trees.append((f"adv{made}", par))
            made += 1
    return trees


def boundary_trees(up):
    leaves, children = up.LF, up.KF
    return [
        ("root leaf overflow", build((leaves + 1, []))),
        ("child leaf overflow", build((1, [(leaves + 2, [(2, [])])]))),
        ("bud at both cutoffs", build((2, [(leaves, [(1, [])] * children)]))),
        ("root child overflow", build((2, [(1, [])] * (children + 2)))),
        ("child capacity exceeded", build((1, [(3, [(1, [])] * (children + 1))]))),
    ]


def check_group(up, lo, trees, label):
    c0 = Fr(2989767809, 10**10)
    C0 = Fr(17527873930, 10**10)
    stats = {"trees": 0, "terms": 0, "id_fail": 0, "loc_fail": 0, "cmp_fail": 0}
    for name, par in trees:
        for C, cb in ((up, c0), (lo, C0)):
            d, p0 = brute_drift(C, par)
            dec = decomposition(C, par)
            s = sum((e for _, e, _ in dec), Fr(0))
            tv = sum((t for _, _, t in dec), Fr(0))
            assert tv == p0, (name, "term inventory != Phi")
            if s != d:
                stats["id_fail"] += 1
                print("IDENTITY FAIL", name, C.mode, float(s - d))
            for kind, e, t in dec:
                ok = e <= C.lam * t if C.mode == 1 else e >= C.lam * t
                if not ok:
                    stats["loc_fail"] += 1
                    print("LOCAL FAIL", name, C.mode, kind, float(e), float(C.lam * t))
            _, ncore = phi(C, par)
            okc = p0 >= cb * ncore if C.mode == 1 else p0 <= cb * ncore
            if not okc:
                stats["cmp_fail"] += 1
                print("CMP FAIL", name, C.mode)
            okd = d <= C.lam * p0 if C.mode == 1 else d >= C.lam * p0
            if not okd:
                print("DRIFT FAIL", name, C.mode)
                stats["loc_fail"] += 1
            stats["terms"] += len(dec)
        stats["trees"] += 1
    maxn = max(len(p) for _, p in trees)
    print(f"{label} trees={stats['trees']} (max n={maxn}) local terms checked={stats['terms']} "
          f"identity_fail={stats['id_fail']} local_fail={stats['loc_fail']} compare_fail={stats['cmp_fail']}")
    return 1 if (stats["id_fail"] or stats["loc_fail"] or stats["cmp_fail"]) else 0


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('upper', help='shipped upper motif certificate')
    cli.add_argument('lower', help='shipped lower motif certificate')
    args = cli.parse_args()
    up, lo = Cert(args.upper), Cert(args.lower)
    if up.mode != 1 or lo.mode != 2:
        cli.error('supply upper first and lower second')
    failures = 0
    for seed in (2, 3, 4, 5):
        failures += check_group(up, lo, sample_trees(up, seed), f"seed={seed}")
    failures += check_group(up, lo, boundary_trees(up), "boundary")
    return int(failures != 0)


if __name__ == "__main__":
    sys.exit(main())
