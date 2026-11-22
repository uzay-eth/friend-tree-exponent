#!/usr/bin/env python3
"""
Independent referee checks for attempts/strogatz-carbon/sol_as_lower/AS_LOWER_PROOF.md.

Written cold from AS_LOWER_PROOF.md and REPORT.md sections 3.1-3.5 only.
It does NOT import or read verify_as_lower.py or any certificate script.

What is checked (exact rational arithmetic unless stated):
  A. one-step facts on random trees: receiver law, identity (2), REPORT Lemma 2 event list
     for the certs9 state and for the upper (L,c,m1,b,P) state (per-vertex drift (1) by brute
     force vs. event-list formula at true rates, plus the rate intervals), jump counts,
     L>=L_max / degree>=H claims, second-moment bound (3);
  B. branch restriction across a saturated parent (section 2): exterior invariance, exact
     restricted drift bookkeeping (4), identity/bound (5), bound (6), replay rates (section 4);
  C. replay invariance of states at pattern depth >= 2 (section 4);
  D. exact enumeration of all histories up to n = 9 (depth<=1 count, sum rule);
  E. analytic inequalities (7), (10), (12), (13);
  F. Monte Carlo (heuristic only).
Every block has negative controls with planted errors that must be DETECTED.
"""
import math
import random
import sys
from fractions import Fraction as Fr

FAILS = []
NEG_MISSED = []


def report(ok, msg):
    print(("[PASS] " if ok else "[FAIL] ") + msg)
    if not ok:
        FAILS.append(msg)


def neg(detected, msg):
    print(("[NEG-CONTROL DETECTED] " if detected else "[NEG-CONTROL MISSED] ") + msg)
    if not detected:
        NEG_MISSED.append(msg)


# ---------------------------------------------------------------------------
# Tree with history (creator = parent in the tree rooted at the centre r of T_3)
# ---------------------------------------------------------------------------
class Tree:
    def __init__(self):
        self.adj = [[1, 2], [0], [0]]   # T_3: centre 0 with leaves 1, 2
        self.par = [-1, 0, 0]
        self.dep = [0, 1, 1]
        self.birth = {}                 # v -> (L_p, c_p, D_p) of creator right after promotion

    @property
    def n(self):
        return len(self.adj)

    def deg(self, v):
        return len(self.adj[v])

    def Lc(self, v):
        L = 0
        for w in self.adj[v]:
            if len(self.adj[w]) == 1:
                L += 1
        return L, len(self.adj[v]) - L

    def add_leaf(self, w):
        promoted = len(self.adj[w]) == 1
        new = len(self.adj)
        self.adj[w].append(new)
        self.adj.append([w])
        self.par.append(w)
        self.dep.append(self.dep[w] + 1)
        if promoted:
            p = self.adj[w][0]
            Lp, cp = self.Lc(p)
            self.birth[w] = (Lp, cp, len(self.adj[p]))
        return new, promoted

    def undo(self, w, promoted):
        self.adj.pop()
        self.adj[w].pop()
        self.par.pop()
        self.dep.pop()
        if promoted:
            del self.birth[w]

    def copy(self):
        t = Tree.__new__(Tree)
        t.adj = [list(a) for a in self.adj]
        t.par = list(self.par)
        t.dep = list(self.dep)
        t.birth = dict(self.birth)
        return t


def all_Lc(T):
    return [T.Lc(v) for v in range(T.n)]


def recv_probs(T):
    """Backward formula of section 1: P(W=w) = n^{-1} sum_{z~w} 1/D_z."""
    n = T.n
    return [Fr(1, n) * sum(Fr(1, T.deg(z)) for z in T.adj[w]) for w in range(n)]


def recv_probs_forward(T):
    """Direct enumeration of (V, W): V uniform, W uniform neighbour of V."""
    n = T.n
    pi = [Fr(0)] * n
    for V in range(n):
        d = T.deg(V)
        for W in T.adj[V]:
            pi[W] += Fr(1, n * d)
    return pi


def credits(T, LC):
    rho = [Fr(0)] * T.n
    for u in range(T.n):
        if T.deg(u) >= 2:
            rho[u] = sum((Fr(1, T.deg(v)) for v in T.adj[u] if T.deg(v) >= 2), Fr(0))
    return rho


def rft_receiver(T, rng):
    V = rng.randrange(T.n)
    return rng.choice(T.adj[V])


def grow(T, steps, rng, bias=0.0):
    for _ in range(steps):
        if bias and rng.random() < bias:
            core = [u for u in range(T.n) if T.deg(u) >= 2]
            w = rng.choice(core)
        else:
            w = rft_receiver(T, rng)
        T.add_leaf(w)
    return T


# ---------------------------------------------------------------------------
# State functions
# ---------------------------------------------------------------------------
def st9(T, LC, u, K0, Dcap):
    """certs9 state (L, c, (m_1..m_K0), b, (n_3..n_{Dcap-1})) of core vertex u; top implicit."""
    L, c = LC[u]
    m = [0] * (K0 + 1)
    b = 0
    nd = [0] * max(Dcap, 3)
    for w in T.adj[u]:
        Lw, cw = LC[w]
        dw = Lw + cw
        if dw == 1:
            continue
        if cw == 1 and 1 <= Lw <= K0:
            m[Lw] += 1
        elif Lw == 0 and cw == 2:
            b += 1
        else:
            assert dw >= 3, "an 'other' neighbour of degree < 3"
            if dw < Dcap:
                nd[dw] += 1
    return (L, c, tuple(m[1:]), b, tuple(nd[3:Dcap]))


def label_P(T, u, Pmax, K0u=1):
    if u not in T.birth:
        return 3                         # centre of T_3: no claim
    Lp, cp, Dp = T.birth[u]
    pendant = (cp == 1 and 1 <= Lp <= K0u)
    classB = (Lp == 0 and cp == 2)
    if pendant or classB:
        return 3
    return min(Dp, Pmax)


def stU(T, LC, u, Pmax):
    """upper state (L, c, m1, b, P) with K0 = 1 and static creator label P."""
    L, c = LC[u]
    m1 = b = 0
    for w in T.adj[u]:
        Lw, cw = LC[w]
        if Lw + cw == 1:
            continue
        if cw == 1 and Lw == 1:
            m1 += 1
        elif Lw == 0 and cw == 2:
            b += 1
    return (L, c, m1, b, label_P(T, u, Pmax))


def st_r2(T, LC, u, K0, Dcap):
    """PLANTED non-radius-1 table: own certs9 state plus the parent's full certs9 state."""
    p = T.par[u]
    ps = st9(T, LC, p, K0, Dcap) if (p >= 0 and T.deg(p) >= 2) else None
    return (st9(T, LC, u, K0, Dcap), ps)


def trunc9(s, Lmax, cmax):
    L, c = s[0], s[1]
    if c > cmax:
        return 'Z'
    return (min(L, Lmax),) + s[1:]


def truncU(s, Lmax, cmax):
    L, c = s[0], s[1]
    if c > cmax:
        return 'G'
    return (min(L, Lmax),) + s[1:]


class RandTable:
    """Random table, lazily filled; values in [lo, hi] with denominator 10^6."""

    def __init__(self, rng, lo, hi, special=None):
        self.rng, self.lo, self.hi = rng, Fr(lo), Fr(hi)
        self.d = {}
        if special:
            self.d.update(special)

    def __call__(self, s):
        v = self.d.get(s)
        if v is None:
            v = self.lo + (self.hi - self.lo) * Fr(self.rng.randint(0, 10**6), 10**6)
            self.d[s] = v
        return v


# ---------------------------------------------------------------------------
# Event lists (REPORT Lemma 2) - transitions computed from the state tuple only
# ---------------------------------------------------------------------------
def _enter(nd, d, Dcap):
    if d < Dcap:
        nd[d - 3] += 1


def nb9(s, K0, Dcap, planted=None):
    L, c = s[0], s[1]
    D = L + c
    Lp, cp = L - 1, c + 1
    m = [0] * K0
    b = 0
    nd = [0] * max(Dcap - 3, 0)
    if planted == 'nb':                        # planted: newborn always sees creator as top
        return (1, 1, tuple(m), 0, tuple(nd))
    if cp == 1 and 1 <= Lp <= K0:
        m[Lp - 1] = 1
    elif Lp == 0 and cp == 2:
        b = 1
    else:
        _enter(nd, D, Dcap)
    return (1, 1, tuple(m), b, tuple(nd))


def predicted9(T, LC, rho, u, s, f, K0, Dcap, planted=None, stats=None):
    """sum_events rate*(f(new)-f(x)) + (L/D) f(nb) at TRUE rates; also checks rate intervals."""
    L, c, m, b, nd = s
    D = L + c
    fs = f(s)
    tot = (L + rho[u]) * (f((L + 1, c, m, b, nd)) - fs)                     # E1
    if L > 0:                                                              # E2
        m2 = list(m)
        if K0 >= 1:
            m2[0] += 1
        tot += Fr(L, D) * (f((L - 1, c + 1, tuple(m2), b, nd)) - fs)
        tot += Fr(L, D) * f(nb9(s, K0, Dcap, planted))
    # credit interval of REPORT 3.3 (tracked-degree variant)
    B0 = sum(Fr(m[l - 1], l + 1) for l in range(1, K0 + 1)) + Fr(b, 2)
    lo = B0 + sum(Fr(nd[d - 3], d) for d in range(3, Dcap))
    top = c - sum(m) - b - sum(nd)
    hi = lo + Fr(top, Dcap)
    if stats is not None:
        stats['rho'] += 1
        if not (lo <= rho[u] <= hi):
            stats['rho_bad'] += 1
    for w in T.adj[u]:
        Lw, cw = LC[w]
        dw = Lw + cw
        if dw == 1:
            continue
        g = Lw + rho[w]                       # true leaf-gain rate of w (times n)
        if cw == 1 and 1 <= Lw <= K0:         # pendant of class l
            l = Lw
            if stats is not None:
                stats['rate'] += 2
                if g != l + Fr(1, D):
                    stats['rate_bad'] += 1
                if Fr(Lw, dw) != Fr(l, l + 1):
                    stats['rate_bad'] += 1
            m2 = list(m); m2[l - 1] -= 1; nd2 = list(nd)                    # E3_l
            if l < K0:
                m2[l] += 1
            else:
                _enter(nd2, K0 + 2, Dcap)
            tot += g * (f((L, c, tuple(m2), b, tuple(nd2))) - fs)
            r4 = Fr(Lw, dw)                                                  # E4_l
            if planted == 'E4rate' and l == 1:
                r4 = Fr(1, 3)
            m2 = list(m); m2[l - 1] -= 1; nd2 = list(nd); b2 = b
            if l == 1:
                b2 += 1
            else:
                _enter(nd2, l + 1, Dcap)
            tot += r4 * (f((L, c, tuple(m2), b2, tuple(nd2))) - fs)
        elif Lw == 0 and cw == 2:             # class B                       E5
            if stats is not None:
                stats['rate'] += 1
                if not (Fr(1, D) <= g <= Fr(1, D) + Fr(1, 2)):
                    stats['rate_bad'] += 1
            nd2 = list(nd)
            _enter(nd2, 3, Dcap)
            tot += g * (f((L, c, m, b - 1, tuple(nd2))) - fs)
        else:                                 # other
            if dw < Dcap and planted != 'noE7':                              # E7_d
                if stats is not None:
                    stats['rate'] += 1
                    if not (0 <= g <= dw):
                        stats['rate_bad'] += 1
                nd2 = list(nd); nd2[dw - 3] -= 1
                _enter(nd2, dw + 1, Dcap)
                tot += g * (f((L, c, m, b, tuple(nd2))) - fs)
    return tot


def nbU(s, Pmax):
    L, c = s[0], s[1]
    D = L + c
    Lp, cp = L - 1, c + 1
    if cp == 1 and Lp == 1:
        return (1, 1, 1, 0, 3)
    if Lp == 0 and cp == 2:
        return (1, 1, 0, 1, 3)
    return (1, 1, 0, 0, min(D, Pmax))


def predictedU(T, LC, rho, u, s, f, Pmax, planted=None, stats=None):
    L, c, m1, b, P = s
    D = L + c
    fs = f(s)
    tot = (L + rho[u]) * (f((L + 1, c, m1, b, P)) - fs)                    # E1
    if L > 0:                                                              # E2
        tot += Fr(L, D) * (f((L - 1, c + 1, m1 + 1, b, P)) - fs)
        tot += Fr(L, D) * f(nbU(s, Pmax))
    o = c - m1 - b
    B0 = Fr(m1, 2) + Fr(b, 2)
    if P >= 4:
        lo, hi = B0, B0 + Fr(o - 1, 3) + Fr(1, P)
    else:
        lo, hi = B0, B0 + Fr(o, 3)
    if stats is not None:
        stats['rho'] += 1
        if not (lo <= rho[u] <= hi) or (P >= 4 and o < 1):
            stats['rho_bad'] += 1
    for w in T.adj[u]:
        Lw, cw = LC[w]
        dw = Lw + cw
        if dw == 1:
            continue
        g = Lw + rho[w]
        if cw == 1 and Lw == 1:                                             # E3_1, E4_1
            if stats is not None:
                stats['rate'] += 1
                if g != 1 + Fr(1, D):
                    stats['rate_bad'] += 1
            tot += g * (f((L, c, m1 - 1, b, P)) - fs)
            r4 = Fr(1, 2) if planted != 'E4rate' else Fr(1, 3)
            tot += r4 * (f((L, c, m1 - 1, b + 1, P)) - fs)
        elif Lw == 0 and cw == 2:                                           # E5
            if stats is not None:
                stats['rate'] += 1
                if not (Fr(1, D) <= g <= Fr(1, D) + Fr(1, 2)):
                    stats['rate_bad'] += 1
            tot += g * (f((L, c, m1, b - 1, P)) - fs)
    return tot


# ---------------------------------------------------------------------------
# One-step enumeration
# ---------------------------------------------------------------------------
def enumerate_step(T, statefuns):
    """For every receiver w: dict name -> (old states of old core, new states, newborn, promoted, creator).
    statefuns: name -> function(T, LC, u)."""
    n = T.n
    LC0 = all_Lc(T)
    core0 = [u for u in range(n) if T.deg(u) >= 2]
    old = {nm: {u: fn(T, LC0, u) for u in core0} for nm, fn in statefuns.items()}
    outs = []
    for w in range(n):
        creator = T.adj[w][0] if T.deg(w) == 1 else None
        new, promoted = T.add_leaf(w)
        LC1 = all_Lc(T)
        newst = {}
        for nm, fn in statefuns.items():
            d = {u: fn(T, LC1, u) for u in core0}
            if promoted:
                d[w] = fn(T, LC1, w)
            newst[nm] = d
        T.undo(w, promoted)
        outs.append((w, promoted, creator, newst))
    return LC0, core0, old, outs


def brute_drift(T, core0, old, outs, nm, f, pi):
    n = T.n
    acc = {u: Fr(0) for u in core0}
    for (w, promoted, creator, newst) in outs:
        d = newst[nm]
        for u in core0:
            val = f(d[u]) - f(old[nm][u])
            if promoted and creator == u:
                val += f(d[w])
            acc[u] += pi[w] * val
    return {u: n * acc[u] for u in core0}


# ---------------------------------------------------------------------------
# PART A
# ---------------------------------------------------------------------------
def part_A(rng):
    print("\n=== PART A: one-step facts on random trees (exact) ===")
    params9 = [(1, 4), (2, 6), (3, 5), (2, 4), (2, 5)]
    paramsU = [3, 6, 20]
    ntrees = 0
    recv_bad = recv_neg = 0
    id2_bad = 0
    drift_cmp = drift_bad = 0
    neg_hits = {'E4rate': 0, 'nb': 0, 'noE7': 0, 'U-E4rate': 0, 'PA-law': 0}
    stats = {'rho': 0, 'rho_bad': 0, 'rate': 0, 'rate_bad': 0}
    statsU = {'rho': 0, 'rho_bad': 0, 'rate': 0, 'rate_bad': 0}
    jump_bad = {'J1': 0, 'J2': 0, 'J3': 0, 'J4a': 0, 'J4b': 0, 'J4c': 0, 'J4d': 0, 'J4U': 0}
    jump_tests = {k: 0 for k in jump_bad}
    neg_J3 = 0
    sm_bad = sm_tests = 0
    gainfrac_bad = 0
    worst_ratio_3 = 0.0
    for t in range(170):
        T = Tree()
        size = rng.randint(1, 55)
        grow(T, size, rng, bias=rng.choice([0.0, 0.0, 0.25, 0.5]))
        ntrees += 1
        n = T.n
        pi = recv_probs(T)
        pif = recv_probs_forward(T)
        if pi != pif or sum(pi) != 1:
            recv_bad += 1
        pa = [Fr(T.deg(w), 2 * (n - 1)) for w in range(n)]
        if pa != pif:
            recv_neg += 1
        LC = all_Lc(T)
        rho = credits(T, LC)
        core = [u for u in range(n) if T.deg(u) >= 2]
        C = len(core)
        # identity (2) and S_n <= C_n
        lhs = sum(rho[u] for u in core)
        rhs = sum(Fr(LC[v][1], T.deg(v)) for v in core)
        S = sum(Fr(LC[u][0], T.deg(u)) for u in core)
        if lhs != rhs or lhs > C or S > C:
            id2_bad += 1
        # receiver rates: core u gets L_u + rho_u, leaves of u get L_u/D_u in total
        for u in core:
            if n * pi[u] != LC[u][0] + rho[u]:
                recv_bad += 1
            leaves = [x for x in T.adj[u] if T.deg(x) == 1]
            if n * sum((pi[x] for x in leaves), Fr(0)) != Fr(LC[u][0], T.deg(u)):
                recv_bad += 1
            # stage-1 claim: gain rate >= promotion rate, so gain fraction >= 1/2
            if LC[u][0] >= 1 and (LC[u][0] + rho[u]) < Fr(LC[u][0], T.deg(u)):
                gainfrac_bad += 1
        # --- per-vertex drift (1) brute force vs event list ---
        sfs = {}
        for (K0, Dc) in params9:
            sfs[('9', K0, Dc)] = (lambda K0, Dc: (lambda T_, LC_, u: st9(T_, LC_, u, K0, Dc)))(K0, Dc)
        for Pm in paramsU:
            sfs[('U', Pm)] = (lambda Pm: (lambda T_, LC_, u: stU(T_, LC_, u, Pm)))(Pm)
        LC0, core0, old, outs = enumerate_step(T, sfs)
        for (K0, Dc) in params9:
            nm = ('9', K0, Dc)
            f = RandTable(rng, 0, 1)
            bd = brute_drift(T, core0, old, outs, nm, f, pi)
            for u in core0:
                pr = predicted9(T, LC0, rho, u, old[nm][u], f, K0, Dc, stats=stats)
                drift_cmp += 1
                if pr != bd[u]:
                    drift_bad += 1
                for pl in ('E4rate', 'nb', 'noE7'):
                    if predicted9(T, LC0, rho, u, old[nm][u], f, K0, Dc, planted=pl) != bd[u]:
                        neg_hits[pl] += 1
        for Pm in paramsU:
            nm = ('U', Pm)
            f = RandTable(rng, 1, 5)
            bd = brute_drift(T, core0, old, outs, nm, f, pi)
            for u in core0:
                pr = predictedU(T, LC0, rho, u, old[nm][u], f, Pm, stats=statsU)
                drift_cmp += 1
                if pr != bd[u]:
                    drift_bad += 1
                if predictedU(T, LC0, rho, u, old[nm][u], f, Pm, planted='E4rate') != bd[u]:
                    neg_hits['U-E4rate'] += 1
        # --- jump counts (raw and truncated) ---
        for (K0, Dc) in params9:
            nm = ('9', K0, Dc)
            Lmax = max(K0 + 1, Dc - 1)          # smallest L_max with L_max>K0, L_max+1>=D_cap
            cmax = rng.randint(2, 6)
            H = max(Lmax, Dc, K0 + 2, 3)
            for (w, promoted, creator, newst) in outs:
                d = newst[nm]
                changed = [u for u in core0 if d[u] != old[nm][u]]
                tchanged = [u for u in core0 if trunc9(d[u], Lmax, cmax) != trunc9(old[nm][u], Lmax, cmax)]
                if not promoted:
                    jump_tests['J1'] += 1
                    if not set(changed) <= set([w] + T.adj[w]) or len(changed) > max(3, Dc):
                        jump_bad['J1'] += 1
                    Lw = LC0[w][0]
                    if Lw >= Lmax:
                        jump_tests['J3'] += 1
                        if tchanged:
                            jump_bad['J3'] += 1
                    if T.deg(w) >= H:
                        jump_tests['J4a'] += 1
                        if [u for u in changed if u != w]:
                            jump_bad['J4a'] += 1
                    if Lw >= H:
                        jump_tests['J4b'] += 1
                        if tchanged:
                            jump_bad['J4b'] += 1
                    jump_tests['J4d'] += 1
                    if len(tchanged) > H + 1:
                        jump_bad['J4d'] += 1
                else:
                    jump_tests['J2'] += 1
                    cnt = len(changed) + 1     # newborn summand
                    others = [u for u in changed if u != creator]
                    pend = (LC0[creator][1] == 1 and LC0[creator][0] <= K0)
                    if cnt > 3 or len(others) > 1 or (others and not pend):
                        jump_bad['J2'] += 1
                    if T.deg(creator) >= H:
                        jump_tests['J4c'] += 1
                        if others:
                            jump_bad['J4c'] += 1
                    jump_tests['J4d'] += 1
                    if len(tchanged) + 1 > H + 1:
                        jump_bad['J4d'] += 1
            # second-moment bound (3) with a random truncated table in [0, 1]
            f = RandTable(rng, 0, 1, special={'Z': Fr(0)})
            Bj = max(3, Dc)
            m2 = Fr(0); pchg = Fr(0); plow = Fr(0)
            for (w, promoted, creator, newst) in outs:
                d = newst[nm]
                dphi = sum((f(trunc9(d[u], Lmax, cmax)) - f(trunc9(old[nm][u], Lmax, cmax)) for u in core0), Fr(0))
                if promoted:
                    dphi += f(trunc9(d[w], Lmax, cmax))
                m2 += pi[w] * dphi * dphi
                if dphi != 0:
                    pchg += pi[w]
                if not promoted and LC0[w][0] < Lmax:
                    plow += pi[w]
            sm_tests += 1
            if not (m2 <= Fr(Bj * Bj * (Lmax + 2) * C, n) and plow <= Fr((Lmax + 1) * C, n)
                    and pchg <= Fr((Lmax + 2) * C, n)):
                sm_bad += 1
            worst_ratio_3 = max(worst_ratio_3, float(m2 / Fr(Bj * Bj * (Lmax + 2) * C, n)))
        # upper table: degree >= H gains / promotions change no neighbour entry
        for Pm in paramsU:
            nm = ('U', Pm)
            H = max(Pm, 3)
            for (w, promoted, creator, newst) in outs:
                d = newst[nm]
                changed = [u for u in core0 if d[u] != old[nm][u]]
                x = w if not promoted else creator
                if T.deg(x) >= H:
                    jump_tests['J4U'] += 1
                    if [u for u in changed if u != x]:
                        jump_bad['J4U'] += 1
        # negative control for J3: violate L_max + 1 >= D_cap  (K0=1, D_cap=6, L_max=2)
        sfneg = {'n': lambda T_, LC_, u: st9(T_, LC_, u, 1, 6)}
        _, c0n, oldn, outsn = enumerate_step(T, sfneg)
        for (w, promoted, creator, newst) in outsn:
            if not promoted and LC0[w][0] >= 2:
                if any(trunc9(newst['n'][u], 2, 50) != trunc9(oldn['n'][u], 2, 50) for u in c0n):
                    neg_J3 += 1
    report(recv_bad == 0, f"receiver law P(W=w)=n^-1 sum_(z~w) 1/D_z (backward = forward enumeration), "
                          f"core rate L+rho, leaf rate L/D: {ntrees} trees, bad={recv_bad}")
    neg(recv_neg == ntrees, f"planted preferential-attachment law D_w/(2(n-1)) differs from true law on "
                            f"{recv_neg}/{ntrees} trees")
    report(id2_bad == 0, f"identity (2) sum_u rho_u = sum_v c_v/D_v <= C_n and S_n <= C_n: bad={id2_bad}")
    report(gainfrac_bad == 0, f"stage-1: gain rate L+rho >= promotion rate L/D for every core vertex: bad={gainfrac_bad}")
    report(drift_bad == 0, f"per-vertex drift (1) by brute force == REPORT event-list formula at true rates "
                           f"(certs9 (K0,Dcap) in {params9}; upper P_max in {paramsU}): "
                           f"{drift_cmp} vertex comparisons, mismatches={drift_bad}")
    report(stats['rho_bad'] == 0 and stats['rate_bad'] == 0,
           f"certs9 true rates inside REPORT intervals: credit {stats['rho']} checks bad={stats['rho_bad']}; "
           f"E3/E4/E5/E7 {stats['rate']} checks bad={stats['rate_bad']}")
    report(statsU['rho_bad'] == 0 and statsU['rate_bad'] == 0,
           f"upper (L,c,m1,b,P) true rates inside REPORT intervals: credit {statsU['rho']} checks "
           f"bad={statsU['rho_bad']}; E3/E5 {statsU['rate']} checks bad={statsU['rate_bad']}")
    for k, v in neg_hits.items():
        if k != 'PA-law':
            neg(v > 0, f"planted event-list error '{k}' changes the formula: mismatching vertices={v}")
    names = {'J1': 'gain: changed entries within receiver+neighbours, count <= max(3,Dcap)',
             'J2': 'promotion: <=3 changed summands, old neighbour changes only if creator was pendant',
             'J3': 'gain at core w with L_w >= L_max (L_max>K0, L_max+1>=Dcap): no truncated entry changes',
             'J4a': 'gain at core w with D_w >= H: no neighbour entry changes',
             'J4b': 'gain at core w with L_w >= H: no truncated entry changes',
             'J4c': 'promotion at creator with D >= H: no old neighbour entry changes',
             'J4d': 'every step: <= H+1 changed truncated summands',
             'J4U': 'upper table: receiver/creator of degree >= H changes no neighbour entry'}
    for k in jump_bad:
        report(jump_bad[k] == 0, f"{k} {names[k]}: tests={jump_tests[k]} bad={jump_bad[k]}")
    neg(neg_J3 > 0, f"planted L_max+1 < D_cap (K0=1,Dcap=6,L_max=2): gains at L>=L_max that change a "
                    f"truncated entry = {neg_J3}")
    report(sm_bad == 0, f"second-moment bound (3) E[(dPhi)^2] <= B^2(L_max+2)C/n, B=F max(3,Dcap), plus "
                        f"P(W core, L_W<L_max) <= (L_max+1)C/n: tests={sm_tests} bad={sm_bad}; "
                        f"max observed ratio LHS/RHS = {worst_ratio_3:.4f}")


# ---------------------------------------------------------------------------
# PART B: branch restriction
# ---------------------------------------------------------------------------
def side_of(T, p, v):
    """vertices on the v-side of edge {p, v}"""
    seen = {v}
    stack = [v]
    while stack:
        x = stack.pop()
        for y in T.adj[x]:
            if y != p and y not in seen:
                seen.add(y)
                stack.append(y)
    return seen


def build_branch(rng, Dp_target):
    T = Tree()
    grow(T, rng.randint(2, 30), rng, bias=0.2)
    core = [u for u in range(T.n) if T.deg(u) >= 2]
    p = rng.choice(core)
    while T.deg(p) < Dp_target or not [x for x in T.adj[p] if T.deg(x) == 1]:
        T.add_leaf(p)
    v = rng.choice([x for x in T.adj[p] if T.deg(x) == 1])
    T.add_leaf(v)                      # the qualifying promotion of v at p
    return T, p, v


def part_B(rng):
    print("\n=== PART B: restriction across a saturated parent edge (section 2), exact ===")
    K0, Dc, Lmax, cmax = 2, 5, 6, 4
    Pm, LmaxU, cmaxU = 5, 6, 4
    H = max(Lmax, Dc, K0 + 2, 3, Pm, LmaxU)       # = 6
    F, U = 1, 5
    tests = dict(ext=0, book=0, b5=0, b6=0, rep=0)
    bad = dict(ext=0, book=0, b5=0, b6=0, rep=0)
    worst6 = 0.0
    neg_ext = neg_book = 0
    for trial in range(60):
        saturated = trial % 4 != 3
        T, p, v = build_branch(rng, H if saturated else 3)
        if not saturated and T.deg(p) >= Dc:
            continue
        fX = RandTable(rng, 0, F, special={'Z': Fr(0)})
        fY = RandTable(rng, 1, U, special={'G': Fr(rng.randint(1, 4))})
        for step in range(rng.randint(10, 40)):
            R = side_of(T, p, v)
            # choose next step: inside, at p / leaf of p, or RFT
            x = rng.random()
            if x < 0.45:
                w = rng.choice(sorted(R))
            elif x < 0.65:
                w = rng.choice([p] + [y for y in T.adj[p] if T.deg(y) == 1 and y not in R] or [p])
            else:
                w = rft_receiver(T, rng)
            if not saturated and w == p and T.deg(p) + 1 >= Dc:
                w = rng.choice(sorted(R))      # keep p unsaturated for the negative control
            # ----- checks at the current tree -----
            n = T.n
            pi = recv_probs(T)
            LC = all_Lc(T)
            rho = credits(T, LC)
            I = sorted(u for u in R if T.deg(u) >= 2)
            N = len(I)
            k = len(R)
            sfs = {'X': lambda T_, LC_, u: trunc9(st9(T_, LC_, u, K0, Dc), Lmax, cmax),
                   'Y': lambda T_, LC_, u: truncU(stU(T_, LC_, u, Pm), LmaxU, cmaxU),
                   'raw9': lambda T_, LC_, u: st9(T_, LC_, u, K0, Dc),
                   'rawU': lambda T_, LC_, u: stU(T_, LC_, u, Pm)}
            LC0, core0, old, outs = enumerate_step(T, sfs)
            Iset = set(I)
            m2X = m2Y = Fr(0)
            pchg = Fr(0)
            EdX_book = Fr(0)
            for (wr, promoted, creator, newst) in outs:
                inside = wr in R
                if not inside:
                    tests['ext'] += 1
                    ch = [u for u in I if newst['raw9'][u] != old['raw9'][u] or
                          newst['rawU'][u] != old['rawU'][u]]
                    if ch:
                        if saturated:
                            bad['ext'] += 1
                        else:
                            neg_ext += 1
                # exact bookkeeping: new inside core = I plus newborn iff creator in I
                T.add_leaf(wr)
                promoted2 = promoted
                R1 = side_of(T, p, v)
                I1 = set(u for u in R1 if T.deg(u) >= 2)
                T.undo(wr, promoted2)
                expect = set(Iset)
                if promoted and creator in Iset:
                    expect.add(wr)
                tests['book'] += 1
                if I1 != expect:
                    bad['book'] += 1
                for nm, f, acc in (('X', fX, 'x'), ('Y', fY, 'y')):
                    Xold = sum((f(old[nm][u]) for u in I), Fr(0))
                    Xnew = sum((f(newst[nm][u]) for u in I1 if u in newst[nm]), Fr(0))
                    per = sum((f(newst[nm][u]) - f(old[nm][u]) for u in I), Fr(0))
                    if promoted and creator in Iset:
                        per += f(newst[nm][wr])
                    if Xnew - Xold != per:
                        bad['book'] += 1
                    if nm == 'X':
                        m2X += pi[wr] * (Xnew - Xold) ** 2
                        if Xnew != Xold:
                            pchg += pi[wr]
                        EdX_book += pi[wr] * (Xnew - Xold)
                        # planted bookkeeping error: forget the newborn
                        if promoted and creator in Iset and fX(newst['X'][wr]) != 0:
                            neg_book += 1
                    else:
                        m2Y += pi[wr] * (Xnew - Xold) ** 2
            # E[dX] equals (1/n) * sum over I of the per-vertex quantity (1)
            bdX = brute_drift(T, core0, old, outs, 'X', fX, pi)
            tests['book'] += 1
            if EdX_book != sum((bdX[u] for u in I), Fr(0)) / n:
                bad['book'] += 1
            # (5)
            lhs = sum((rho[u] for u in I), Fr(0))
            mid = sum((Fr(len([u for u in T.adj[z] if u in Iset]), T.deg(z)) for z in I), Fr(0)) + Fr(1, T.deg(p))
            tests['b5'] += 1
            if lhs != mid or lhs > N + 1 or N + 1 > 2 * N:
                bad['b5'] += 1
            # (6) and the (H+2)N/n change probability (only claimed for saturated p)
            if saturated:
                BF = (H + 1) * F
                BU = (H + 1) * U
                tests['b6'] += 1
                if not (pchg <= Fr((H + 2) * N, n) and m2X <= Fr((H + 3) * BF * BF * N, n)
                        and m2Y <= Fr((H + 3) * BU * BU * N, n)):
                    bad['b6'] += 1
                worst6 = max(worst6, float(m2X / Fr((H + 3) * BF * BF * N, n)))
            # replay rates (section 4)
            tests['rep'] += 1
            totR = sum((pi[x] for x in R), Fr(0))
            if any(pi[x] < Fr(1, n * k) for x in R) or totR > Fr(k + 1, n):
                bad['rep'] += 1
            # ---- advance ----
            T.add_leaf(w)
    report(bad['ext'] == 0, f"saturated parent (D_p >= H={H} >= D_cap={Dc}, P_max={Pm}): every exterior receiver "
                            f"(incl. p and leaves of p) leaves all raw inside certs9 and upper states unchanged: "
                            f"exterior outcomes={tests['ext']} bad={bad['ext']}")
    neg(neg_ext > 0, f"planted unsaturated parent (D_p < D_cap): exterior outcomes that change an inside state "
                     f"= {neg_ext}")
    report(bad['book'] == 0, f"restricted bookkeeping (4): I_(n+1) = I_n + newborn iff creator in I_n; "
                             f"X_(n+1)-X_n = sum over I_n of summand changes + newborn, for X and Y, every outcome; "
                             f"E[dX] = n^-1 sum_(u in I) (1)_u: checks={tests['book']} bad={bad['book']}")
    neg(neg_book > 0, f"planted bookkeeping error (newborn omitted) changes dX on {neg_book} outcomes")
    report(bad['b5'] == 0, f"(5): sum_(u in I) rho_u = sum_(z in I) |N(z) cap I|/D_z + 1/D_p <= N+1 <= 2N: "
                           f"checks={tests['b5']} bad={bad['b5']}")
    report(bad['b6'] == 0, f"(6): P(X changes) <= (H+2)N/n, E[dX^2] <= (H+3)B_F^2 N/n, E[dY^2] <= (H+3)B_U^2 N/n: "
                           f"checks={tests['b6']} bad={bad['b6']}; max ratio for X = {worst6:.5f}")
    report(bad['rep'] == 0, f"replay rates: each inside vertex has P(W=x) >= 1/(nk), P(W in branch) <= (k+1)/n: "
                            f"checks={tests['rep']} bad={bad['rep']}")


# ---------------------------------------------------------------------------
# PART C: replay invariance
# ---------------------------------------------------------------------------
def part_C(rng):
    print("\n=== PART C: replay of a pattern inside a saturated branch (section 4) ===")
    params9 = [(1, 4), (2, 6), (2, 5), (3, 7)]
    Pm = 8
    H = 9
    cmp_deep = bad_deep = 0
    root_bad = 0
    diff_d1 = 0
    diffU_d1 = 0
    neg_r2 = 0
    size_bad = 0
    for trial in range(120):
        kstar = rng.randint(6, 28)
        P = Tree()
        T, p, v = build_branch(rng, H)
        a = T.adj[v][1]                       # the leaf created at the promotion
        b, _ = T.add_leaf(v)                  # "first attach another leaf to v"
        phi = {0: v, 1: a, 2: b}
        inside = {v, a, b}
        for i in range(kstar - 3):
            # exterior noise, including p and leaves of p
            for _ in range(rng.randint(0, 3)):
                x = rng.random()
                if x < 0.4:
                    w = p
                else:
                    w = rft_receiver(T, rng)
                    if w in inside:
                        w = p
                T.add_leaf(w)
            # pattern step
            if rng.random() < 0.5:
                wp = rft_receiver(P, rng)
            else:
                wp = rng.randrange(P.n)
            newp, _ = P.add_leaf(wp)
            newh, _ = T.add_leaf(phi[wp])
            phi[newp] = newh
            inside.add(newh)
            # compare after every step
            LP = all_Lc(P)
            LT = all_Lc(T)
            for xp in range(P.n):
                if P.deg(xp) < 2:
                    if T.deg(phi[xp]) >= 2:
                        bad_deep += 1
                    continue
                xh = phi[xp]
                if xp == 0:
                    if LP[0][0] != LT[xh][0] or LP[0][1] + 1 != LT[xh][1]:
                        root_bad += 1
                    continue
                for (K0, Dc) in params9:
                    s1 = st9(P, LP, xp, K0, Dc)
                    s2 = st9(T, LT, xh, K0, Dc)
                    if P.dep[xp] >= 2:
                        cmp_deep += 1
                        if s1 != s2:
                            bad_deep += 1
                        if st_r2(P, LP, xp, K0, Dc) != st_r2(T, LT, xh, K0, Dc):
                            neg_r2 += 1
                    elif s1 != s2:
                        diff_d1 += 1
                u1 = stU(P, LP, xp, Pm)
                u2 = stU(T, LT, xh, Pm)
                if P.dep[xp] >= 2:
                    cmp_deep += 1
                    if u1 != u2:
                        bad_deep += 1
                elif u1 != u2:
                    diffU_d1 += 1
        if len(side_of(T, p, v)) != kstar or P.n != kstar:
            size_bad += 1
    report(bad_deep == 0, f"depth>=2 pattern vertices: raw certs9 states ((K0,Dcap) in {params9}) and upper "
                          f"(L,c,m1,b,P<= {Pm}) states identical in pattern and replayed branch, after every replay "
                          f"step with exterior noise at p: comparisons={cmp_deep} mismatches={bad_deep}")
    report(root_bad == 0, f"copy of r: same L, one extra core neighbour: bad={root_bad}")
    report(size_bad == 0, f"branch has exactly k_* vertices at the end: bad={size_bad}")
    print(f"[INFO] depth-1 vertices whose certs9 state differs (expected, the proof only uses depth>=2): {diff_d1}; "
          f"upper-state differences at depth 1: {diffU_d1}")
    neg(diff_d1 > 0, "depth-1 states do change in the copy, so the restriction to depth>=2 is needed")
    neg(neg_r2 > 0, f"planted radius-2 table (own state + parent's full state): depth>=2 mismatches = {neg_r2}")


# ---------------------------------------------------------------------------
# PART D: exact enumeration of all histories
# ---------------------------------------------------------------------------
def part_D():
    print("\n=== PART D: exact enumeration of all receiver histories from T_3 up to n = 9 ===")
    NMAX = 9
    EC = {}
    ED1 = {}
    ES = {}
    ELr = {}
    nhist = [0]

    def rec(T, prob):
        n = T.n
        LC = all_Lc(T)
        core = [u for u in range(n) if T.deg(u) >= 2]
        d1 = sum(1 for u in core if T.dep[u] <= 1)
        EC[n] = EC.get(n, Fr(0)) + prob * len(core)
        ED1[n] = ED1.get(n, Fr(0)) + prob * d1
        if n == NMAX:
            nhist[0] += 1
            return
        S = sum(Fr(LC[u][0], T.deg(u)) for u in core)
        ES[n] = ES.get(n, Fr(0)) + prob * S
        ELr[n] = ELr.get(n, Fr(0)) + prob * Fr(LC[0][0], T.deg(0))
        pi = recv_probs(T)
        for w in range(n):
            new, promoted = T.add_leaf(w)
            rec(T, prob * pi[w])
            T.undo(w, promoted)

    rec(Tree(), Fr(1))
    ok_sum = all(EC[n + 1] - EC[n] == ES[n] / n for n in range(3, NMAX))
    ok_d1 = all(ED1[n] <= 1 + sum(Fr(1, j) for j in range(3, n)) for n in range(3, NMAX + 1))
    ok_d1inc = all(ED1[n + 1] - ED1[n] == ELr[n] / n and ELr[n] <= 1 for n in range(3, NMAX))
    print(f"[INFO] histories enumerated to n={NMAX}: {nhist[0]}")
    for n in range(3, NMAX + 1):
        print(f"[INFO] n={n}: E[C_n]={float(EC[n]):.6f}  E[#core depth<=1]={float(ED1[n]):.6f}  "
              f"bound 1+H={float(1 + sum(Fr(1, j) for j in range(3, n))):.6f}")
    report(ok_sum, "sum rule E[C_(n+1)] - E[C_n] = E[S_n]/n exactly, n=3..8")
    report(ok_d1 and ok_d1inc, "E[#core of depth<=1] increments by E[L_r/D_r]/n <= 1/n and is <= 1+sum_(j=3)^(n-1) 1/j")
    report(ED1[NMAX] < EC[NMAX], "some reachable history has a core vertex of depth >= 2 (E[C] > E[depth<=1 core])")


# ---------------------------------------------------------------------------
# PART E: analytic inequalities
# ---------------------------------------------------------------------------
def logP(a, m, n):
    return sum(math.log1p(a / j) for j in range(m, n))


def part_E():
    print("\n=== PART E: analytic inequalities (7), (10), (12), (13) ===")
    bad7 = 0
    n7 = 0
    negU = negL = 0
    for a in (0.1, 0.5, 0.543, 0.586, 1.0, 1.086, 1.7, 3.0):
        for m in (3, 4, 5, 7, 10, 31, 100, 1000):
            lp = 0.0
            j = m
            for n in sorted(set([m, m + 1, m + 2, m + 5] + [int(m * 1.3 ** e) for e in range(1, 45)])):
                while j < n:
                    lp += math.log1p(a / j)
                    j += 1
                base = a * math.log(n / m)
                n7 += 1
                if not (base - a * a - 1e-12 <= lp <= base + a + 1e-12):
                    bad7 += 1
                if lp > base + 1e-12:
                    negU += 1
                if lp < base - 1e-12:
                    negL += 1
    report(bad7 == 0, f"(7) e^(-a^2)(n/m)^a <= P_a(m,n) <= e^a (n/m)^a, a in [0.1,3], m>=3: {n7} cases, bad={bad7}")
    neg(negU > 0, f"planted (7) without the e^a factor (P <= (n/m)^a) fails in {negU} cases")
    neg(negL > 0, f"planted (7) without the e^(-a^2) factor (P >= (n/m)^a) fails in {negL} cases")
    # harmonic comparison and discrepancy bound
    badh = 0
    for m in (3, 4, 10, 100):
        s = 0.0
        for n in range(m + 1, 5000):
            s += 1.0 / (n - 1)
            dlt = s - math.log(n / m)
            if not (-1e-12 <= dlt <= 1.0 / m + 1e-12):
                badh += 1
    badd = 0
    for a in (0.543, 0.586, 2.0):
        for j in range(3, 20000):
            if abs(math.log1p(a / j) - a * math.log1p(1.0 / j)) > (a * a + a) / (2 * j * j) + 1e-15:
                badd += 1
    report(badh == 0 and badd == 0, f"harmonic 0 <= sum 1/j - log(n/m) <= 1/m and |log(1+a/j)-a log(1+1/j)| <= "
                                    f"(a^2+a)/(2j^2): bad={badh + badd}")
    # Gamma identity / finite positive limit of P_a(m,n)/n^a
    badg = 0
    for a in (0.543, 0.586):
        for n in (4, 10, 1000, 100000):
            lhs = logP(a, 3, n)
            rhs = math.log(2) + math.lgamma(n + a) - math.lgamma(n) - math.lgamma(3 + a)
            if abs(lhs - rhs) > 1e-9:
                badg += 1
    report(badg == 0, "P_a(3,n) = 2 Gamma(n+a)/(Gamma(n) Gamma(3+a)) (so P_a(3,n)/n^a has a positive finite limit)")
    # (10)
    d, lam = 0.543, 0.586
    q = 2 * d - lam
    const = math.exp(lam + 2 * d * d)
    bad10 = 0
    for m in (3, 10, 100, 1000):
        J = 300000
        lpl = 0.0
        lpd = 0.0
        S = 0.0
        termbad = 0
        for j in range(m, J):
            lpd1 = lpd + math.log1p(d / j)                      # log P_delta(m, j+1)
            term = math.exp(lpl - 2 * lpd1) / j
            if term > const * (j / m) ** (-q) / j * (1 + 1e-12):
                termbad += 1
            S += term
            lpl += math.log1p(lam / j)
            lpd = lpd1
        tail = const * ((J - 1) / m) ** (-q) / q              # sum_{j>=J} g(j) <= int_{J-1}^inf g
        bound = const * (1.0 / m + 1.0 / q)
        ok = termbad == 0 and S + tail <= bound
        if not ok:
            bad10 += 1
        print(f"[INFO] (10) m={m}: partial sum to J={J} = {S:.5f}, + tail bound {tail:.5f} = {S + tail:.5f} "
              f"<= e^(lam+2del^2)(1/m+1/q) = {bound:.5f}")
    report(bad10 == 0, "(10) per-term bound and total bound with delta=0.543, lambda=0.586, q=0.5")
    # negative control q = 0 : lambda = 2 delta  -> divergence
    lam2 = 2 * d
    lpl = lpd = 0.0
    S = 0.0
    checkpoints = {}
    for j in range(3, 2000000):
        lpd1 = lpd + math.log1p(d / j)
        S += math.exp(lpl - 2 * lpd1) / j
        lpl += math.log1p(lam2 / j)
        lpd = lpd1
        if j + 1 in (10**3, 10**4, 10**5, 10**6, 2 * 10**6 - 1):
            checkpoints[j + 1] = S
    grows = all(checkpoints[b] - checkpoints[a] > 0.5 for a, b in ((10**3, 10**4), (10**4, 10**5), (10**5, 10**6)))
    neg(grows, "planted lambda = 2 delta (q=0): the sum in (10) grows ~log J: " +
        ", ".join(f"J={k}: {v:.3f}" for k, v in sorted(checkpoints.items())))
    # (12) Cauchy-Schwarz lemma
    rng = random.Random(12)
    bad12 = 0
    neg12 = 0
    for _ in range(20000):
        k = rng.randint(2, 6)
        vals = [Fr(rng.randint(-40, 200), 10) for _ in range(k)]
        ws = [Fr(rng.randint(1, 20)) for _ in range(k)]
        tot = sum(ws)
        ps = [w / tot for w in ws]
        a = sum(p * x for p, x in zip(ps, vals))
        if a <= 0:
            continue
        EJ2 = sum(p * x * x for p, x in zip(ps, vals))
        pr = sum(p for p, x in zip(ps, vals) if x > a / 2)
        if pr < a * a / (4 * EJ2):
            bad12 += 1
        if pr < Fr(1, 4):
            neg12 += 1
    report(bad12 == 0, "(12) for any J with E J = a > 0: P(J > a/2) >= a^2/(4 E J^2) on 20000 random laws: "
                       f"bad={bad12}")
    neg(neg12 > 0, f"planted variance-free version P(J > a/2) >= 1/4 fails on {neg12} laws")
    # (13)
    for ks in (4, 6, 10):
        pp = Fr(1)
        for k in range(2, ks):
            pp *= Fr(1, k * (k + 1))
        print(f"[INFO] (13) k_*={ks}: p_pattern = {pp} ~ {float(pp):.3e} > 0")


# ---------------------------------------------------------------------------
# PART F: Monte Carlo (heuristic only)
# ---------------------------------------------------------------------------
def part_F(rng):
    print("\n=== PART F: Monte Carlo (heuristic only, proves nothing) ===")
    runs = 12
    NM = 200000
    cps = [10**3, 10**4, 10**5, 2 * 10**5]
    table = {c: [] for c in cps}
    for r in range(runs):
        adj = [[1], [0]]
        C = 0
        n = 2
        while n < NM:
            V = rng.randrange(n)
            W = rng.choice(adj[V])
            if len(adj[W]) == 1:
                C += 1
            adj[W].append(n)
            adj.append([W])
            n += 1
            if n in table:
                table[n].append(C / n ** 0.543)
    for c in cps:
        xs = sorted(table[c])
        print(f"[INFO] n={c}: C_n/n^0.543 over {runs} runs: min={xs[0]:.3f} median={xs[len(xs)//2]:.3f} "
              f"max={xs[-1]:.3f}")
    # stage 1 of section 4 with a small H
    Hs = 8
    trials = 150
    succ = 0
    bound_sum = 0.0
    for t in range(trials):
        T = grow(Tree(), 200, rng)
        leaf = rng.choice([x for x in range(T.n) if T.deg(x) == 1])
        u = T.adj[leaf][0]
        D0 = T.deg(u)
        bound_sum += 2.0 ** (-max(Hs - D0, 0))
        ok = None
        while ok is None:
            if T.deg(u) >= Hs:
                ok = True
                break
            w = rft_receiver(T, rng)
            if T.deg(w) == 1 and T.adj[w][0] == u:
                ok = False
            T.add_leaf(w)
        succ += ok
    print(f"[INFO] stage 1 (H={Hs}): reached D_u >= H before any promotion at u in {succ}/{trials} trials; "
          f"mean of the proof's lower bound 2^-(H-D_u) = {bound_sum / trials:.4f}")
    report(succ / trials >= bound_sum / trials - 0.1, "stage-1 empirical success frequency is above the proof's bound")


def main():
    seed = int(sys.argv[1]) if len(sys.argv) > 1 else 20261122
    rng = random.Random(seed)
    print(f"seed = {seed}")
    part_A(rng)
    part_B(rng)
    part_C(rng)
    part_D()
    part_E()
    part_F(rng)
    print("\n=== SUMMARY ===")
    print(f"failures: {len(FAILS)}")
    for m in FAILS:
        print("  FAIL:", m)
    print(f"negative controls missed: {len(NEG_MISSED)}")
    for m in NEG_MISSED:
        print("  MISSED:", m)
    print("ALL CHECKS PASSED" if not FAILS and not NEG_MISSED else "SOME CHECKS FAILED")


if __name__ == '__main__':
    main()
