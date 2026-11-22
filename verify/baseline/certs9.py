"""LOWER certificate with exactly tracked degrees of 'other' neighbours (REPORT.md, Section 3.7).

State of a core vertex u (a function of the current tree only):
  x = (L, c, m, b, nn):  m = (m_1..m_K0) pendant classes, b = class B, nn = (n_3, ..., n_{Dcap-1}) where
  n_d = number of 'other' core neighbours of u whose CURRENT degree is exactly d; the remaining
  top = c - sum(m) - b - sum(nn) 'other' neighbours have degree >= Dcap.
Events (rate = n * probability per step), in addition to E1, E2 and the pendant/B events of certs6.py:
  E3_K0 : pendant K0 gains a leaf -> state (K0+1, 1), degree K0+2   -> n_{K0+2} + 1 (or top if K0+2 >= Dcap)
  E4_l (l >= 2): pendant l converts -> (l-1, 2), degree l+1          -> n_{l+1} + 1
  E5    : class B gains a leaf -> (1, 2), degree 3                   -> n_3 + 1
  E7_d  : an 'other' neighbour of degree d gains a leaf, rate n_d * s, s = L_w + rho_w in [0, d]
                                                                     -> n_d - 1, n_{d+1} + 1 (or top)
  credit rho_u in [B0 + sum_d n_d/d , B0 + sum_d n_d/d + top/Dcap]
  newborn: u seen as pendant class L-1 (c = 0, 1 <= L-1 <= K0), class B ((L,c) = (1,1)), or as an 'other'
  neighbour of degree D = L + c (n_D = 1 if D < Dcap, else top).
LOWER inequality: min over rho, r (class B), s (E7) of the drift >= delta f.

usage: python3 certs9.py DELTA_LO DELTA_HI STEPS DELTA_CERT Lmax cmax K0 Dcap   (DELTA_CERT '-' = no certificate)
"""
import sys
from fractions import Fraction as Fr
from itertools import product


class Model:
    def __init__(self, Lmax, cmax, K, Dcap):
        assert K >= 1 and Dcap >= 4
        self.Lmax, self.cmax, self.K, self.Dcap = Lmax, cmax, K, Dcap
        self.nd = Dcap - 3                         # tracked degrees 3..Dcap-1
        cl = []
        for c in range(cmax + 1):
            for m in product(range(c + 1), repeat=K):
                if sum(m) > c:
                    continue
                for b in range(c + 1 - sum(m)):
                    for nn in product(range(c + 1), repeat=self.nd):
                        if sum(m) + b + sum(nn) <= c:
                            cl.append((c, m, b, nn))
        def rk(t):
            c, m, b, nn = t
            top = c - sum(m) - b - sum(nn)
            return (sum((l + 1) * m[l] for l in range(K)) + (K + 1) * b
                    + sum((K + 2 + i) * nn[i] for i in range(self.nd)) + (K + 2 + self.nd) * top)
        cl.sort(key=lambda t: (-t[0], -rk(t)))
        self.classes = cl
        self.NB = []
        for L in range(1, Lmax + 2):
            for c in range(0, cmax + 2):
                x = self.nb(L, c)
                if x not in self.NB:
                    self.NB.append(x)

    def add_other(self, nn, d):
        nn = list(nn)
        if d < self.Dcap:
            nn[d - 3] += 1
        return tuple(nn)

    def nb(self, L, c):
        K = self.K
        m = [0] * K; b = 0; nn = tuple([0] * self.nd)
        if c == 0 and 1 <= L - 1 <= K:
            m[L - 2] = 1
        elif (L, c) == (1, 1):
            b = 1
        elif L + c >= 3:
            nn = self.add_other(nn, L + c)
        return (1, 1, tuple(m), b, nn)

    def rho_bounds(self, c, m, b, nn, F=float):
        base = sum(F(m[l]) / (l + 2) for l in range(self.K)) + F(b) / 2 + sum(F(nn[i]) / (i + 3) for i in range(self.nd))
        top = c - sum(m) - b - sum(nn)
        return base, base + F(top) / self.Dcap

    def moves(self, c, m, b, nn):
        """list of (kind, l, target (m, b, nn)); kinds: 'gain' (pendant), 'conv' (pendant), 'B', 'grow' (other)."""
        K = self.K
        out = []
        for l in range(K):
            if m[l] == 0:
                continue
            m2 = list(m); m2[l] -= 1
            if l + 1 < K:
                m2[l + 1] += 1
                out.append(('gain', l, (tuple(m2), b, nn)))
            else:
                out.append(('gain', l, (tuple(m2), b, self.add_other(nn, K + 2))))
            m3 = list(m); m3[l] -= 1
            if l == 0:
                out.append(('conv', l, (tuple(m3), b + 1, nn)))
            else:
                out.append(('conv', l, (tuple(m3), b, self.add_other(nn, l + 2))))
        if b:
            out.append(('B', None, (m, b - 1, self.add_other(nn, 3))))
        for i in range(self.nd):
            if nn[i]:
                d = i + 3
                n2 = list(nn); n2[i] -= 1
                n2 = self.add_other(tuple(n2), d + 1)
                out.append(('grow', i, (m, b, n2)))
        return out


def rate_num(kind, l, m, b, nn, iD):
    """pessimistic (LOWER) numerical rates: B at max 1/D+1/2, grow at max d."""
    if kind == 'gain':
        return m[l] * ((l + 1) + iD)
    if kind == 'conv':
        return m[l] * (l + 1) / (l + 2)
    if kind == 'B':
        return b * (iD + 0.5)
    return nn[l] * (l + 3)


def solve(M, delta, sv):
    Lmax, cmax = M.Lmax, M.cmax
    f = {}
    def g(L, c, m, b, nn):
        if c > cmax:
            return 0.0
        return f[(min(L, Lmax), c, m, b, nn)]
    for (c, m, b, nn) in M.classes:
        mv = M.moves(c, m, b, nn)
        e1m = list(m); e1m[0] += 1; e1m = tuple(e1m)
        lo_r, _ = M.rho_bounds(c, m, b, nn)
        for L in range(Lmax, -1, -1):
            D = L + c
            if D < 2:
                f[(L, c, m, b, nn)] = 0.0
                continue
            iD = 1.0 / D; bL = L / D
            num = bL * ((g(L - 1, c + 1, e1m, b, nn) if L >= 1 else 0.0) + sv[M.nb(L, c)])
            den = bL + delta
            for (kind, l, (m2, b2, n2)) in mv:
                r = rate_num(kind, l, m, b, nn, iD)
                num += r * f[(L, c, m2, b2, n2)]; den += r
            if L < Lmax:
                up = L + lo_r
                num += up * f[(L + 1, c, m, b, nn)]; den += up
            f[(L, c, m, b, nn)] = num / den if den > 0 else 0.0
    return f


def lower_table(M, delta):
    import numpy as np
    gs = [solve(M, delta, {y: (1.0 if y == x else 0.0) for y in M.NB}) for x in M.NB]
    A = np.array([[gs[j][x] for j in range(len(M.NB))] for x in M.NB])
    w, V = np.linalg.eig(A)
    i = int(np.argmax(w.real)); rho = w[i].real
    v = np.abs(V[:, i].real)
    f = {k: sum(v[j] * gs[j][k] for j in range(len(gs))) for k in gs[0]}
    return f, rho


def verify(M, ftab, delta):
    Lmax, cmax = M.Lmax, M.cmax
    F = {k: Fr(v) for k, v in ftab.items()}
    def f(L, c, m, b, nn):
        if c > cmax:
            return Fr(0)
        return F[(min(L, Lmax), c, m, b, nn)]
    snb = {x: f(*x) for x in M.NB}
    n = 0
    for (c, m, b, nn) in M.classes:
        mv = M.moves(c, m, b, nn)
        e1m = list(m); e1m[0] += 1; e1m = tuple(e1m)
        lo_r, hi_r = M.rho_bounds(c, m, b, nn, F=Fr)
        def drift(L, bLs, iDs, sb, dn):
            v = f(L, c, m, b, nn)
            dL = f(L + 1, c, m, b, nn) - v
            tot = min((L + lo_r) * dL, (L + hi_r) * dL)
            if L >= 1:
                tot += min(bL * (dn - v + sb) for bL in bLs)
            for (kind, l, (m2, b2, n2)) in mv:
                d = f(L, c, m2, b2, n2) - v
                if kind == 'gain':
                    tot += min(m[l] * ((l + 1) + iD) * d for iD in iDs)
                elif kind == 'conv':
                    tot += m[l] * Fr(l + 1, l + 2) * d
                elif kind == 'B':
                    tot += min(b * r * d for iD in iDs for r in (iD, iD + Fr(1, 2)))
                else:                      # grow: rate in [0, d]
                    tot += min(Fr(0), nn[l] * (l + 3) * d)
            return tot, v
        for L in range(Lmax + 1):
            D = L + c
            if D < 2:
                continue
            tot, v = drift(L, [Fr(L, D)], [Fr(1, D)], snb[M.nb(L, c)],
                           f(L - 1, c + 1, e1m, b, nn) if L >= 1 else Fr(0))
            assert v >= 0 and tot >= delta * v, ("ineq", L, c, m, b, nn, float(tot - delta * v))
            n += 1
        Dm = Lmax + 1 + c
        assert Lmax >= M.K + 1 and Lmax + 1 >= M.Dcap
        tot, v = drift(Lmax + 1, [Fr(Lmax + 1, Dm), Fr(1)], [Fr(0), Fr(1, Dm)],
                       snb[M.nb(Lmax + 1, c)], f(Lmax, c + 1, e1m, b, nn))
        assert tot >= delta * v, ("tail", c, m, b, nn)
        n += 1
    return snb, max(F.values()), n


if __name__ == '__main__':
    lo, hi, steps = float(sys.argv[1]), float(sys.argv[2]), int(sys.argv[3])
    cert = sys.argv[4]
    Lmax, cmax, K, Dcap = map(int, sys.argv[5:9])
    M = Model(Lmax, cmax, K, Dcap)
    print(f"model: {len(M.classes)} classes, {len(M.NB)} newborn types", flush=True)
    for _ in range(steps):
        mid = (lo + hi) / 2
        ok = lower_table(M, mid)[1] > 1
        lo, hi = (mid, hi) if ok else (lo, mid)
        print(f"  delta={mid:.5f} {'ok' if ok else 'fail'}", flush=True)
    print('lower9', Lmax, cmax, K, Dcap, 'threshold in', (lo, hi), flush=True)
    if cert != '-':
        delta = Fr(cert)
        ftab, rho = lower_table(M, float(lo))
        assert rho >= 1
        snb, sup, n = verify(M, ftab, delta)
        x0 = (2, 0, tuple([0] * K), 0, tuple([0] * M.nd))
        from mpmath import mp, gamma, mpf
        mp.dps = 30
        kap = Fr(ftab[x0]) / sup
        const = mpf(kap.numerator) / kap.denominator * 2 / gamma(3 + mpf(delta.numerator) / delta.denominator)
        print(f"LOWER9(K={K},Dcap={Dcap}) certificate OK: delta={delta} (={float(delta):.4f}) Lmax={Lmax} cmax={cmax} "
              f"rho={rho:.6f} sup f={float(sup):.5f} f(T3 centre)={ftab[x0]:.5f}; {n} checks", flush=True)
        print(f"  => E[C_n] >= {float(kap):.8f} * prod_(m=3..n-1)(1+delta/m) = {float(const):.8f} * Gamma(n+delta)/Gamma(n)  (n>=3)")
