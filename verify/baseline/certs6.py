"""General certificate: state x = (L, c, m, b, P) of a core vertex u (REPORT.md, Section 4).

  L, c, D = L + c as before;
  m = (m_1, ..., m_K): m_l = # core neighbours w in state (l, 1) ("pendant of class l": w's only
      core neighbour is u, w has l leaves, degree l+1);
  b = # core neighbours in state (0, 2) (class B);
  o = c - sum(m) - b "other" core neighbours, each of degree >= 3;
  P in {3..Pmax}: static label (UPPER only; LOWER uses Pmax = 3).  P = min(D_parent at birth, Pmax)
      if the parent was an "other" neighbour of the newborn, else 3.  Requires K <= 2 when Pmax > 3
      (then a parent of degree >= 4 can be neither pendant of class <= K nor class B).

Exact event list for u (rate = n * probability per step), with exact effect on x_u:
  E1   u gains a leaf                rate L + rho                     L -> L+1
  E2   u converts a leaf            rate L/D                          L-1, c+1, m_1+1; newborn nb(x)
  E3_l class-l pendant gains a leaf  rate m_l (l + 1/D)              m_l-1, m_{l+1}+1 (or -> other if l = K)
  E4_1 class-1 pendant converts     rate m_1 / 2                     m_1-1, b+1
  E4_l class-l pendant converts (l>=2) rate m_l * l/(l+1)            m_l-1 (-> other)
  E5   class-B nbr gains a leaf      rate b * r, r in [1/D, 1/D+1/2]  b-1 (-> other)
  credit rho in [ sum_l m_l/(l+1) + b/2 ,  sum_l m_l/(l+1) + b/2 + o/3 ]  (P = 3), and with
         o/3 replaced by (o-1)/3 + 1/P when P >= 4.
  Newborn nb(x): the new vertex w is (1,1); u (after the conversion) is seen by w as:
     pendant class L-1 if c = 0 and 1 <= L-1 <= K;  class B if (L,c) = (1,1);  otherwise "other",
     in which case w's label is P_w = max(3, min(D, Pmax)).

usage: python3 certs6.py lower DELTA_NUM DELTA_CERT Lmax cmax K
       python3 certs6.py upper LAM_NUM  LAM_CERT  Lmax cmax K Pmax
"""
import sys
from fractions import Fraction as Fr
from itertools import product


class Model:
    def __init__(self, mode, Lmax, cmax, K, Pmax):
        assert Pmax == 3 or K <= 2
        self.mode, self.Lmax, self.cmax, self.K, self.Pmax = mode, Lmax, cmax, K, Pmax
        # class tuples (c, m, b), ordered so that every transition target precedes its source
        cl = []
        for c in range(cmax + 1):
            for m in product(range(c + 1), repeat=K):
                if sum(m) > c:
                    continue
                for b in range(c + 1 - sum(m)):
                    cl.append((c, m, b))
        def key(t):
            c, m, b = t
            tot = sum(m) + b
            psi = sum((l + 1) * m[l] for l in range(K)) + (K + 1) * b
            return (-c, tot, -psi)
        cl.sort(key=key)
        self.classes = cl
        self.NB = []
        for x in [self.nb(L, c) for L in range(1, Lmax + 2) for c in range(0, cmax + 2)]:
            if x not in self.NB:
                self.NB.append(x)

    def nb(self, L, c):
        K, Pmax = self.K, self.Pmax
        m = [0] * K
        b = 0
        P = 3
        if c == 0 and 1 <= L - 1 <= K:
            m[L - 2] = 1
        elif (L, c) == (1, 1):
            b = 1
        else:
            P = max(3, min(L + c, Pmax))
        return (1, 1, tuple(m), b, P)

    def rho_bounds(self, c, m, b, P, F=float):
        base = sum(F(m[l]) / (l + 2) for l in range(self.K)) + F(b) / 2
        o = c - sum(m) - b
        if P >= 4 and o >= 1:
            hi = base + F(o - 1) / 3 + F(1) / P
        else:
            hi = base + F(o) / 3
        return base, hi

    def moves(self, c, m, b):
        """class-changing moves (excluding E1, E2): list of (kind, l, new (m,b)); rates given later."""
        K = self.K
        out = []
        for l in range(K):
            if m[l] == 0:
                continue
            m2 = list(m); m2[l] -= 1
            if l + 1 < K:
                m2[l + 1] += 1
            out.append(('gain', l, (tuple(m2), b)))
            m3 = list(m); m3[l] -= 1
            out.append(('conv', l, (tuple(m3), b + 1 if l == 0 else b)))
        if b:
            out.append(('B', None, (m, b - 1)))
        return out


def rate(kind, l, m, b, invD, rB):
    """rate of a class move; for 'B' the per-neighbour rate rB is supplied (interval endpoint)."""
    if kind == 'gain':
        return m[l] * ((l + 1) + invD)
    if kind == 'conv':
        return m[l] * (l + 1) / (l + 2)
    return b * rB


def solve(M, par, sv, G=None):
    """lower: linear value function with newborn rewards sv; upper: floor-1 super-solution."""
    Lmax, cmax, low = M.Lmax, M.cmax, M.mode == 'lower'
    f = {}
    def g(L, c, m, b, P):
        if c > cmax:
            return 0.0 if low else G
        return f[(min(L, Lmax), c, m, b, P)]
    Ps = range(3, M.Pmax + 1)
    for P in Ps:
        for (c, m, b) in M.classes:
            mv = M.moves(c, m, b)
            e1m = list(m); e1m[0] += 1; e1m = tuple(e1m) if M.K else m
            if not low:
                Dm = Lmax + 1 + c
                combos = [(Lmax / (Lmax + c), 1.0 / (Lmax + c), g(Lmax - 1, c + 1, e1m, b, P), sv[M.nb(Lmax, c)])]
                for bb in ((Lmax + 1) / Dm, 1.0):
                    for iD in (0.0, 1.0 / Dm):
                        combos.append((bb, iD, g(Lmax, c + 1, e1m, b, P), sv[M.nb(Lmax + 1, c)]))
                cands = [1.0]
                for (bb, iD, dn, sb) in combos:
                    num = bb * (dn + sb); den = bb + par
                    for (kind, l, (m2, b2)) in mv:
                        r = rate(kind, l, m, b, iD, iD)
                        num += r * g(Lmax, c, m2, b2, P); den += r
                    cands.append(num / den)
                f[(Lmax, c, m, b, P)] = max(cands) * (1 + 1e-12)
                Ls = range(Lmax - 1, -1, -1)
            else:
                Ls = range(Lmax, -1, -1)
            lo_r, hi_r = M.rho_bounds(c, m, b, P)
            for L in Ls:
                D = L + c
                if D < 2:
                    f[(L, c, m, b, P)] = 0.0 if low else 1.0
                    continue
                iD = 1.0 / D
                bL = L / D
                num = bL * ((g(L - 1, c + 1, e1m, b, P) if L >= 1 else 0.0) + sv[M.nb(L, c)])
                den = bL + par
                for (kind, l, (m2, b2)) in mv:
                    r = rate(kind, l, m, b, iD, (iD + 0.5) if low else iD)
                    num += r * f[(L, c, m2, b2, P)]; den += r
                if L < Lmax or not low:
                    up = L + (lo_r if low else hi_r)
                    num += up * g(L + 1, c, m, b, P); den += up
                v = num / den if den > 0 else 0.0
                f[(L, c, m, b, P)] = v if low else max(1.0, v * (1 + 1e-12))
    return f


def lower_table(M, delta):
    import numpy as np
    gs = []
    for x in M.NB:
        sv = {y: (1.0 if y == x else 0.0) for y in M.NB}
        gs.append(solve(M, delta, sv))
    A = np.array([[gs[j][x] for j in range(len(M.NB))] for x in M.NB])
    w, V = np.linalg.eig(A)
    i = int(np.argmax(w.real)); rho = w[i].real
    v = np.abs(V[:, i].real)
    f = {k: sum(v[j] * gs[j][k] for j in range(len(gs))) for k in gs[0]}
    return f, rho


def upper_table(M, lam):
    sv = {x: 1.0 for x in M.NB}
    for _ in range(40):
        G = max(1.0, max(sv.values()) / lam) * (1 + 1e-9)
        f = solve(M, lam, sv, G)
        if all(f[x] <= sv[x] for x in M.NB):
            return f, G
        sv = {x: max(sv[x], f[x]) for x in M.NB}
    base = dict(sv); t = 1.0
    while t < 1e4:
        sv = {x: base[x] * t for x in M.NB}
        G = max(1.0, max(sv.values()) / lam) * (1 + 1e-9)
        f = solve(M, lam, sv, G)
        if all(f[x] <= sv[x] for x in M.NB):
            return f, G
        t *= 1.03
    return None


def verify(M, ftab, par, G=None):
    low = M.mode == 'lower'
    Lmax, cmax = M.Lmax, M.cmax
    F = {k: Fr(v) for k, v in ftab.items()}
    if G is not None:
        G = Fr(G)
    def f(L, c, m, b, P):
        if c > cmax:
            return Fr(0) if low else G
        return F[(min(L, Lmax), c, m, b, P)]
    snb = {x: f(*x) for x in M.NB}
    ext = min if low else max
    if not low:
        assert G >= 1 and max(snb.values()) <= par * G
    n = 0
    for P in range(3, M.Pmax + 1):
        for (c, m, b) in M.classes:
            mv = M.moves(c, m, b)
            e1m = list(m); e1m[0] += 1; e1m = tuple(e1m)
            lo_r, hi_r = M.rho_bounds(c, m, b, P, F=Fr)
            def drift(L, bLs, iDs, sb, dn):
                v = f(L, c, m, b, P)
                dL = f(L + 1, c, m, b, P) - v
                tot = ext((L + lo_r) * dL, (L + hi_r) * dL)
                if L >= 1:
                    tot += ext(bL * (dn - v + sb) for bL in bLs)
                for (kind, l, (m2, b2)) in mv:
                    d = f(L, c, m2, b2, P) - v
                    if kind == 'B':
                        tot += ext(b * r * d for iD in iDs for r in (iD, iD + Fr(1, 2)))
                    elif kind == 'gain':
                        tot += ext(m[l] * ((l + 1) + iD) * d for iD in iDs)
                    else:
                        tot += m[l] * Fr(l + 1, l + 2) * d
                return tot, v
            for L in range(Lmax + 1):
                D = L + c
                if D < 2:
                    continue
                tot, v = drift(L, [Fr(L, D)], [Fr(1, D)], snb[M.nb(L, c)],
                               f(L - 1, c + 1, e1m, b, P) if L >= 1 else Fr(0))
                if low:
                    assert v >= 0 and tot >= par * v, ("ineq", L, c, m, b, P, float(tot - par * v))
                else:
                    assert v >= 1 and tot <= par * v, ("ineq", L, c, m, b, P, float(tot - par * v))
                n += 1
            # tail L > Lmax: newborn type nb(L,c) is the same for all L > Lmax (L-1 > K, D > Pmax)
            assert Lmax >= M.K + 1 and Lmax + 1 >= M.Pmax
            Dm = Lmax + 1 + c
            tot, v = drift(Lmax + 1, [Fr(Lmax + 1, Dm), Fr(1)], [Fr(0), Fr(1, Dm)],
                           snb[M.nb(Lmax + 1, c)], f(Lmax, c + 1, e1m, b, P))
            assert (tot >= par * v) if low else (tot <= par * v), ("tail", c, m, b, P)
            n += 1
    return snb, max(F.values()), n


if __name__ == '__main__':
    kind = sys.argv[1]
    num = float(sys.argv[2]); cert = Fr(sys.argv[3])
    Lmax, cmax, K = int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6])
    Pmax = int(sys.argv[7]) if kind == 'upper' and len(sys.argv) > 7 else 3
    M = Model(kind, Lmax, cmax, K, Pmax)
    if kind == 'lower':
        ftab, rho = lower_table(M, num)
        assert rho >= 1, f"Perron root {rho} < 1"
        snb, sup, n = verify(M, ftab, cert)
        x0 = (2, 0, tuple([0] * K), 0, 3)
        print(f"LOWER(K={K}) certificate OK: delta={cert} (={float(cert):.4f}) Lmax={Lmax} cmax={cmax} rho={rho:.6f} "
              f"sup f={float(sup):.5f} f(T3 centre)={ftab[x0]:.5f}; {n} checks")
    else:
        out = upper_table(M, num)
        assert out is not None, "no numerical super-solution"
        ftab, G = out
        snb, sup, n = verify(M, ftab, cert, G)
        x0 = (2, 0, tuple([0] * K), 0, 3)
        print(f"UPPER(K={K},Pmax={Pmax}) certificate OK: lambda={cert} (={float(cert):.4f}) Lmax={Lmax} cmax={cmax} "
              f"sup f={float(sup):.4f} f(T3 centre)={ftab[x0]:.4f}; {n} checks")


# ------------------------------------------------------------------ linear (Perron) construction of
# UPPER certificates: solve  lam f - A_worst f = 1  (reward 1 gives strictly positive f and slack 1),
# with the newborn values and the value G of the region c > cmax as unknowns of a small linear system.
def solve_upper_lin(M, lam, sv, G, rew):
    Lmax, cmax = M.Lmax, M.cmax
    f = {}
    def g(L, c, m, b, P):
        if c > cmax:
            return G
        return f[(min(L, Lmax), c, m, b, P)]
    for P in range(3, M.Pmax + 1):
        for (c, m, b) in M.classes:
            mv = M.moves(c, m, b)
            e1m = list(m); e1m[0] += 1; e1m = tuple(e1m)
            lo_r, hi_r = M.rho_bounds(c, m, b, P)
            for L in range(Lmax, -1, -1):
                D = L + c
                if D < 2:
                    f[(L, c, m, b, P)] = 0.0
                    continue
                iD = 1.0 / D; bL = L / D
                dn = g(L - 1, c + 1, e1m, b, P) if L >= 1 else 0.0
                if L == Lmax:      # pessimistic tail equation (covers all L >= Lmax): L/D -> 1, 1/D -> 0,
                    iD, bL = 0.0, 1.0                    # E2 target of the tail states L > Lmax
                    dn = g(Lmax, c + 1, e1m, b, P)
                num = rew + bL * (dn + sv[M.nb(L, c)])
                den = bL + lam
                for (kind, l, (m2, b2)) in mv:
                    r = rate(kind, l, m, b, iD, iD)
                    num += r * f[(L, c, m2, b2, P)]; den += r
                if L < Lmax:
                    up = L + hi_r
                    num += up * f[(L + 1, c, m, b, P)]; den += up
                f[(L, c, m, b, P)] = num / den
    return f


def upper_lin_table(M, lam):
    import numpy as np
    zero = {x: 0.0 for x in M.NB}
    h = solve_upper_lin(M, lam, zero, 0.0, 1.0)
    gO = solve_upper_lin(M, lam, zero, 1.0, 0.0)
    gs = [solve_upper_lin(M, lam, {y: (1.0 if y == x else 0.0) for y in M.NB}, 0.0, 0.0) for x in M.NB]
    k = len(M.NB)
    A = np.array([[gs[j][M.NB[i]] for j in range(k)] for i in range(k)])
    go = np.array([gO[x] for x in M.NB])
    jstar = 0
    for _ in range(k + 1):
        w = np.zeros(k); w[jstar] = 1.0          # G = (1 + s_jstar)/lam, valid if s_jstar = max s
        B = A + np.outer(go, w) / lam
        rho = max(abs(np.linalg.eigvals(B)))
        if rho >= 1:
            return None, rho
        s = np.linalg.solve(np.eye(k) - B, np.array([h[x] for x in M.NB]) + go / lam)
        if s[jstar] >= s.max() - 1e-12:
            break
        jstar = int(np.argmax(s))
    G = (1.0 + s[jstar]) / lam
    f = {key: h[key] + sum(s[j] * gs[j][key] for j in range(k)) + G * gO[key] for key in h}
    mn = min(v for key, v in f.items() if key[0] + key[1] >= 2)
    scale = (1.0 + 1e-9) / mn
    f = {key: v * scale for key, v in f.items()}
    return (f, G * scale * (1 + 1e-9)), rho
