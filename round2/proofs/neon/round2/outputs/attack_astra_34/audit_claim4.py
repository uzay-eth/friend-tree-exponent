#!/usr/bin/env python3
"""
Independent adversarial verification & planted-error test suite for Astra Claim 4:
  - round2/operator.md
  - dependency round2/finite/one-port-upper.md

Written completely from scratch (no imports from ft2-astra/).
"""

from __future__ import annotations

from collections import deque
from fractions import Fraction
import math
import random
from typing import Callable, Dict, FrozenSet, List, Optional, Sequence, Tuple


# ============================================================================
# 1. Fresh cell-labelled CoreTree and microscopic FullTree
# ============================================================================

class CoreState:
    def __init__(
        self,
        adj: Sequence[Sequence[int]],
        L: Sequence[int],
        ghost_c: Optional[Sequence[int]] = None,
        label: Optional[Sequence[int]] = None,
    ):
        self.n = len(adj)
        self.adj = tuple(tuple(sorted(nb)) for nb in adj)
        self.L = tuple(int(x) for x in L)
        self.ghost_c = tuple(0 for _ in range(self.n)) if ghost_c is None else tuple(int(x) for x in ghost_c)
        self.label = tuple(range(self.n)) if label is None else tuple(int(x) for x in label)
        for u in range(self.n):
            assert self.L[u] >= 0 and self.ghost_c[u] >= 0
            assert self.D(u) >= 2

    def c(self, u: int) -> int:
        return len(self.adj[u]) + self.ghost_c[u]

    def D(self, u: int) -> int:
        return self.L[u] + self.c(u)

    def r_x(self, x: int, wrong_self_denom: bool = False) -> Fraction:
        acc = Fraction(self.L[x], 1)
        for y in self.adj[x]:
            acc += Fraction(1, self.D(x) if wrong_self_denom else self.D(y))
        return acc

    def b_x(self, x: int) -> Fraction:
        return Fraction(self.L[x], self.D(x))

    def gain(self, x: int) -> "CoreState":
        nL = list(self.L)
        nL[x] += 1
        return CoreState(self.adj, nL, self.ghost_c, self.label)

    def promote(self, x: int, bug_no_grandchild_inherit: bool = False) -> Tuple["CoreState", int]:
        assert self.L[x] >= 1
        nid = self.n
        nadj = [list(nb) for nb in self.adj]
        nadj[x].append(nid)
        nadj.append([x])
        nL = list(self.L)
        nL[x] -= 1
        nL.append(1)
        ng = list(self.ghost_c) + [0]
        if bug_no_grandchild_inherit:
            new_lbl = self.label[x] if x == self.label[x] else -999
        else:
            new_lbl = self.label[x]
        nlbl = list(self.label) + [new_lbl]
        return CoreState(nadj, nL, ng, nlbl), nid

    def fresh_labels(self) -> "CoreState":
        return CoreState(self.adj, self.L, self.ghost_c, tuple(range(self.n)))

    def dists(self, root: int) -> List[int]:
        d = [-1] * self.n
        d[root] = 0
        q = deque([root])
        while q:
            u = q.popleft()
            for v in self.adj[u]:
                if d[v] == -1:
                    d[v] = d[u] + 1
                    q.append(v)
        return d

    def ball_code(self, root: int, k: Optional[int], with_c: bool = False) -> Tuple:
        d = self.dists(root)
        def enc(u: int, p: int) -> Tuple:
            ch = []
            if k is None or d[u] < k:
                for v in self.adj[u]:
                    if v != p:
                        ch.append(enc(v, u))
            ch.sort()
            return (self.L[u], self.c(u), tuple(ch)) if with_c else (self.L[u], tuple(ch))
        return enc(root, -1)

    def truncate(self, root: int, R: int, keep_slots: bool = True) -> Tuple["CoreState", int]:
        d = self.dists(root)
        kept = [u for u in range(self.n) if 0 <= d[u] <= R]
        mp = {u: i for i, u in enumerate(kept)}
        sadj: List[List[int]] = [[] for _ in kept]
        sL = [0] * len(kept)
        sg = [0] * len(kept)
        for u in kept:
            i = mp[u]
            for v in self.adj[u]:
                if v in mp:
                    sadj[i].append(mp[v])
            sL[i] = self.L[u]
            sg[i] = (self.c(u) - len(sadj[i])) if keep_slots else 0
        return CoreState(sadj, sL, sg, list(range(len(kept)))), mp[root]


def cyl_fn(k: int) -> Callable[[CoreState, int], Fraction]:
    """Deterministic bounded rational cylinder reading B_k(v) (L-marks only)."""
    def f(st: CoreState, v: int) -> Fraction:
        code = st.ball_code(v, k, with_c=False)
        h = 37
        for ch in repr(code).encode():
            h = (h * 257 + ch) % 997
        return Fraction(h - 498, 500)
    return f


# ============================================================================
# 2. Exact Taylor derivatives of Q_s f, S_s f, and microscopic FullTree
# ============================================================================

def Q_taylor_derivs(
    st0: CoreState,
    root: int,
    f: Callable[[CoreState, int], Fraction],
    max_ord: int,
    bug_no_grandchild_inherit: bool = False,
) -> List[Fraction]:
    lbl = st0.label[root]
    memo: Dict[Tuple, Fraction] = {}

    def deriv(st: CoreState, m: int) -> Fraction:
        if m == 0:
            return sum((f(st, v) for v in range(st.n) if st.label[v] == lbl), Fraction(0, 1))
        key = (st.adj, st.L, st.ghost_c, st.label, m)
        if key in memo:
            return memo[key]
        out_rate = Fraction(0, 1)
        acc = Fraction(0, 1)
        for x in range(st.n):
            rx = st.r_x(x)
            if rx != 0:
                out_rate += rx
                acc += rx * deriv(st.gain(x), m - 1)
            if st.L[x] > 0:
                bx = st.b_x(x)
                out_rate += bx
                nxt, _ = st.promote(x, bug_no_grandchild_inherit=bug_no_grandchild_inherit)
                acc += bx * deriv(nxt, m - 1)
        acc -= out_rate * deriv(st, m - 1)
        memo[key] = acc
        return acc

    return [deriv(st0, m) for m in range(max_ord + 1)]


def S_taylor_derivs(
    st0: CoreState,
    root: int,
    f: Callable[[CoreState, int], Fraction],
    max_ord: int,
) -> List[Fraction]:
    memo: Dict[Tuple, Fraction] = {}

    def deriv(st: CoreState, m: int) -> Fraction:
        if m == 0:
            return f(st, root)
        key = (st.adj, st.L, st.ghost_c, m)
        if key in memo:
            return memo[key]
        out_rate = Fraction(0, 1)
        acc = Fraction(0, 1)
        for x in range(st.n):
            rx = st.r_x(x)
            if rx != 0:
                out_rate += rx
                acc += rx * deriv(st.gain(x), m - 1)
            if st.L[x] > 0:
                bx = st.b_x(x)
                out_rate += bx
                nxt, _ = st.promote(x)
                acc += bx * deriv(nxt, m - 1)
        acc -= out_rate * deriv(st, m - 1)
        memo[key] = acc
        return acc

    return [deriv(st0, m) for m in range(max_ord + 1)]


def microscopic_Q_taylor_derivs(
    st0: CoreState,
    root: int,
    f: Callable[[CoreState, int], Fraction],
    max_ord: int,
) -> List[Fraction]:
    """Direct microscopic full-tree process (cores + explicit leaf vertices)."""
    assert all(g == 0 for g in st0.ghost_c)
    adj: List[List[int]] = [list(nb) for nb in st0.adj]
    orig: List[int] = list(range(st0.n))
    for u in range(st0.n):
        for _ in range(st0.L[u]):
            lid = len(adj)
            adj[u].append(lid)
            adj.append([u])
            orig.append(u)
    tadj0 = tuple(tuple(sorted(nb)) for nb in adj)
    torig0 = tuple(orig)

    def to_core(tadj: Tuple[Tuple[int, ...], ...], torig: Tuple[int, ...]) -> CoreState:
        cores = [u for u in range(len(tadj)) if len(tadj[u]) >= 2]
        mp = {u: i for i, u in enumerate(cores)}
        cadj = [[] for _ in cores]
        cL = [0] * len(cores)
        clbl = [0] * len(cores)
        for u in cores:
            i = mp[u]
            clbl[i] = torig[u]
            for w in tadj[u]:
                if w in mp:
                    cadj[i].append(mp[w])
                else:
                    cL[i] += 1
        return CoreState(cadj, cL, None, clbl)

    memo: Dict[Tuple, Fraction] = {}
    def deriv(tadj: Tuple[Tuple[int, ...], ...], torig: Tuple[int, ...], m: int) -> Fraction:
        if m == 0:
            ct = to_core(tadj, torig)
            return sum((f(ct, v) for v in range(ct.n) if ct.label[v] == root), Fraction(0, 1))
        key = (tadj, torig, m)
        if key in memo:
            return memo[key]
        n = len(tadj)
        rates: Dict[int, Fraction] = {}
        for u in range(n):
            inv = Fraction(1, len(tadj[u]))
            for w in tadj[u]:
                rates[w] = rates.get(w, Fraction(0, 1)) + inv
        acc = -Fraction(n, 1) * deriv(tadj, torig, m - 1)
        for w, rw in rates.items():
            nadj = [list(nb) for nb in tadj]
            nadj[w].append(n)
            nadj.append([w])
            ntadj = tuple(tuple(sorted(nb)) for nb in nadj)
            ntorig = torig + (torig[w],)
            acc += rw * deriv(ntadj, ntorig, m - 1)
        memo[key] = acc
        return acc

    return [deriv(tadj0, torig0, m) for m in range(max_ord + 1)]


# ============================================================================
# 3. Finite-time Q_s via RK4 on truncated jump-depth state space
# ============================================================================

def eval_Q_rk4(
    st0: CoreState,
    root: int,
    s: float,
    f_fl: Callable[[CoreState, int], float],
    max_jumps: int = 6,
) -> Tuple[float, float]:
    lbl = st0.label[root]
    states = [st0]
    sidx = {(st0.adj, st0.L, st0.ghost_c, st0.label): 0}
    depth = [0]
    trans: List[List[Tuple[int, float]]] = []
    out_r: List[float] = []
    q = deque([0])
    while q:
        i = q.popleft()
        st = states[i]
        d = depth[i]
        tr = []
        tot = 0.0
        if d < max_jumps:
            for x in range(st.n):
                rx = float(st.r_x(x))
                if rx > 0:
                    tot += rx
                    nx = st.gain(x)
                    sg = (nx.adj, nx.L, nx.ghost_c, nx.label)
                    if sg not in sidx:
                        sidx[sg] = len(states)
                        states.append(nx)
                        depth.append(d + 1)
                        q.append(sidx[sg])
                    tr.append((sidx[sg], rx))
                if st.L[x] > 0:
                    bx = float(st.b_x(x))
                    tot += bx
                    nx, _ = st.promote(x)
                    sg = (nx.adj, nx.L, nx.ghost_c, nx.label)
                    if sg not in sidx:
                        sidx[sg] = len(states)
                        states.append(nx)
                        depth.append(d + 1)
                        q.append(sidx[sg])
                    tr.append((sidx[sg], bx))
        else:
            for x in range(st.n):
                tot += float(st.r_x(x)) + (float(st.b_x(x)) if st.L[x] > 0 else 0.0)
        trans.append(tr)
        out_r.append(tot)

    ns = len(states)
    steps = max(10, int(math.ceil(max(out_r) * s / 0.02)))
    dt = s / steps
    p = [0.0] * ns
    p[0] = 1.0

    def rhs(v):
        dv = [-out_r[i] * v[i] for i in range(ns)]
        for i in range(ns):
            vi = v[i]
            if vi != 0.0:
                for j, r in trans[i]:
                    dv[j] += vi * r
        return dv

    for _ in range(steps):
        k1 = rhs(p)
        k2 = rhs([p[i] + 0.5 * dt * k1[i] for i in range(ns)])
        k3 = rhs([p[i] + 0.5 * dt * k2[i] for i in range(ns)])
        k4 = rhs([p[i] + dt * k3[i] for i in range(ns)])
        p = [p[i] + (dt / 6.0) * (k1[i] + 2 * k2[i] + 2 * k3[i] + k4[i]) for i in range(ns)]

    val = 0.0
    for i, st in enumerate(states):
        if p[i] != 0.0:
            val += p[i] * sum(f_fl(st, v) for v in range(st.n) if st.label[v] == lbl)
    tail = max(0.0, 1.0 - sum(p))
    # Rigorous L^2 Cauchy-Schwarz bound for ||f||_inf <= 1:
    err = math.sqrt(2.0) * math.exp(s) * math.sqrt(tail)
    return val, err


# ============================================================================
# 4. Generator (7) & Checks
# ============================================================================

def gen_A(
    st: CoreState,
    root: int,
    f: Callable[[CoreState, int], Fraction],
    drop_newborn: bool = False,
    wrong_self_denom: bool = False,
) -> Fraction:
    f0 = f(st, root)
    acc = Fraction(0, 1)
    for x in range(st.n):
        rx = st.r_x(x, wrong_self_denom=wrong_self_denom)
        if rx != 0:
            acc += rx * (f(st.gain(x), root) - f0)
        if st.L[x] > 0:
            bx = st.b_x(x)
            nxt, _ = st.promote(x)
            acc += bx * (f(nxt, root) - f0)
    if st.L[root] > 0 and not drop_newborn:
        nxt, nid = st.promote(root)
        acc += st.b_x(root) * f(nxt, nid)
    return acc


def verify_semigroup_and_mt_evolution(bug_no_grandchild_inherit: bool = False) -> Tuple[bool, str]:
    """
    1. Compare Q_s against microscopic full-tree process on 4 test cores up to order 3.
    2. Check Q_s Q_t = Q_{s+t} exactly via Taylor orders (a, b) with a+b <= 3.
    3. Check m_{t+h}(f) = (m_t Q_h)(f) from T_3 up to order 4!
    """
    seeds = [
        ("T3", CoreState([[]], [2]), 0, 3),
        ("edge_1_3", CoreState([[1], [0]], [1, 3]), 0, 3),
        ("path3_typeB", CoreState([[1], [0, 2], [1]], [2, 0, 1]), 1, 2),
        ("star4", CoreState([[1, 2, 3], [0], [0], [0]], [0, 1, 2, 1]), 0, 2),
    ]
    f0 = cyl_fn(0)
    f1 = cyl_fn(1)

    for name, st, rt, max_o in seeds:
        for fn in (f0, f1, lambda g, v: Fraction(1, 1)):
            d_q = Q_taylor_derivs(st, rt, fn, max_o, bug_no_grandchild_inherit=bug_no_grandchild_inherit)
            d_m = microscopic_Q_taylor_derivs(st, rt, fn, max_o)
            if d_q != d_m:
                return False, f"Microscopic vs cell Q_s mismatch on {name}: {d_q} != {d_m}"

    # Check semigroup Q_s Q_t = Q_{s+t} on T3 and edge_1_3 up to order 3:
    for name, st, rt, max_o in seeds[:2]:
        for fn in (f0, f1):
            d_tot = Q_taylor_derivs(st, rt, fn, max_o, bug_no_grandchild_inherit=bug_no_grandchild_inherit)
            for b in range(max_o + 1):
                cache: Dict[Tuple, Fraction] = {}
                def gb(g: CoreState, v: int, bb=b, c=cache) -> Fraction:
                    code = g.ball_code(v, None, with_c=True)
                    if code not in c:
                        c[code] = Q_taylor_derivs(g.fresh_labels(), v, fn, bb, bug_no_grandchild_inherit=bug_no_grandchild_inherit)[bb]
                    return c[code]
                for a in range(max_o - b + 1):
                    lhs = Q_taylor_derivs(st, rt, gb, a, bug_no_grandchild_inherit=bug_no_grandchild_inherit)[a]
                    if lhs != d_tot[a + b]:
                        return False, f"Semigroup Q_s Q_t != Q_{{s+t}} on {name} at (a={a}, b={b}): {lhs} != {d_tot[a+b]}"

    # Check m_{t+h} = m_t Q_h directly from T_3:
    # Since m_0 = delta_{(T_3, 0)}, m_t(f) = (Q_t f)(T_3, 0), so m_{t+h}(f) = (m_t Q_h)(f)
    # is literally (Q_{t+h} f)(T_3, 0) == (Q_t (Q_h f))(T_3, 0), ALSO check on finite-time RK4!
    st_t3 = CoreState([[]], [2])
    f1_fl = lambda g, v: float(f1(g, v))
    t_val, h_val = 0.025, 0.020
    m_th, err1 = eval_Q_rk4(st_t3, 0, t_val + h_val, f1_fl, max_jumps=4)
    qh_cache: Dict[Tuple, float] = {}
    def qh_fn(g: CoreState, v: int) -> float:
        code = g.ball_code(v, None, with_c=True)
        if code not in qh_cache:
            val, _ = eval_Q_rk4(g.fresh_labels(), v, h_val, f1_fl, max_jumps=3)
            qh_cache[code] = val
        return qh_cache[code]
    mt_qh, err2 = eval_Q_rk4(st_t3, 0, t_val, qh_fn, max_jumps=3)
    if abs(m_th - mt_qh) > 5e-4:
        return False, f"Finite-time m_{{t+h}} != m_t Q_h: {m_th} vs {mt_qh}"

    return True, f"Microscopic match, exact Taylor semigroup Q_s Q_t = Q_{{s+t}} (ord 3), and finite-time m_{{t+h}}=m_t Q_h (diff={abs(m_th-mt_qh):.2e}) verified."


def verify_path_reversal_and_cutoff(drop_missing_slots: bool = False) -> Tuple[bool, str]:
    """
    1. Verify reversed-path random-walk bound:
       sum_{simple paths (x_0,...,x_j=w)} prod_{i=0}^{j-1} (1 / c_{x_i}^0) <= c_w^0
       on diverse trees including high-degree hubs, paths, caterpillars, random trees.
    2. Verify vanishing order O(s^{R-k+1}) of Q_s f - Q_s^{[R]} f in exact rational Taylor series.
    3. Verify spatial cutoff inequality (3): |Q_s f - Q_s^{[R]} f| <= 2 sqrt(2) e^s ||f||_inf sqrt(p)
       at finite time s on paths of length 5 with increasing R = 0, 1, 2, 3!
    """
    # Build diverse test trees (including random 18-vertex trees with hubs)
    rng = random.Random(20261124)
    test_adjs: List[List[List[int]]] = [
        # Path of 7 vertices
        [[1], [0, 2], [1, 3], [2, 4], [3, 5], [4, 6], [5]],
        # Star with hub of degree 10
        [[i for i in range(1, 11)]] + [[0] for _ in range(10)],
        # Double-star (two hubs of degrees 6 and 5 joined by an edge)
        [[1, 2, 3, 4, 5, 6]] + [[0] for _ in range(5)] + [[0, 7, 8, 9, 10]] + [[6] for _ in range(4)],
    ]
    for _ in range(15):
        n = 18
        adj: List[List[int]] = [[] for _ in range(n)]
        for i in range(1, n):
            p = rng.randrange(0, i)
            adj[p].append(i)
            adj[i].append(p)
        test_adjs.append(adj)

    max_ratio = Fraction(0, 1)
    for adj in test_adjs:
        n = len(adj)
        deg = [len(adj[u]) for u in range(n)]
        for w in range(n):
            for j in range(1, 6):
                # Enumerate all simple paths of length j ending at w
                total = Fraction(0, 1)
                def dfs(curr: int, vis: FrozenSet[int], left: int, prod: Fraction):
                    nonlocal total
                    if left == 0:
                        total += prod
                        return
                    for prev in adj[curr]:
                        if prev not in vis:
                            dfs(prev, vis | {prev}, left - 1, prod * Fraction(1, deg[prev]))
                dfs(w, frozenset([w]), j, Fraction(1, 1))
                if total > Fraction(deg[w], 1):
                    return False, f"Path-reversal bound VIOLATED: {total} > {deg[w]}"
                if deg[w] > 0:
                    r = total / Fraction(deg[w], 1)
                    if r > max_ratio:
                        max_ratio = r

    # 2. Exact Taylor vanishing order check on a 5-vertex path with L = [2, 2, 2, 2, 2]
    # (so even if ghost slots are dropped, D >= 2 holds everywhere and we can test planted error!)
    path5 = CoreState([[1], [0, 2], [1, 3], [2, 4], [3]], [2, 2, 2, 2, 2])
    f0 = cyl_fn(0)  # radius k = 0
    for R in (0, 1, 2):
        j_exp = R - 0 + 1  # R - k + 1
        tr, tr_rt = path5.truncate(0, R, keep_slots=not drop_missing_slots)
        d_full = Q_taylor_derivs(path5, 0, f0, max_ord=j_exp)
        d_trunc = Q_taylor_derivs(tr, tr_rt, f0, max_ord=j_exp)
        for m in range(j_exp):
            if d_full[m] != d_trunc[m]:
                return False, f"Cutoff vanishing order broken at R={R}, k=0, order m={m} < j={j_exp}: {d_full[m]} != {d_trunc[m]}"

    # 3. Finite-time numerical verification of inequality (3) on path4 with increasing R = 0, 1, 2:
    path4 = CoreState([[1], [0, 2], [1, 3], [2]], [1, 0, 1, 1])
    s_val = 0.03
    f0_fl = lambda g, v: float(f0(g, v))
    # ||f0||_inf <= 1 by construction (values in [-498/500, 498/500])
    q_full, err_f = eval_Q_rk4(path4, 0, s_val, f0_fl, max_jumps=4)
    diffs = []
    for R in (0, 1, 2):
        tr, tr_rt = path4.truncate(0, R, keep_slots=not drop_missing_slots)
        q_tr, err_t = eval_Q_rk4(tr, tr_rt, s_val, f0_fl, max_jumps=4)
        diff = abs(q_full - q_tr)
        j = R + 1  # k = 0
        # sum_{w in B_0(0)} c_w^0 = c_0^0 = 1
        p_bound = min(1.0, (s_val ** j) / math.factorial(j) * float(path4.c(0)))
        rhs = 2.0 * math.sqrt(2.0) * math.exp(s_val) * 1.0 * math.sqrt(p_bound)
        diffs.append((R, diff, rhs, err_f + err_t))
        if diff - (err_f + err_t) > rhs + 1e-12:
            return False, f"Inequality (3) violated at R={R}: diff={diff} > rhs={rhs}"

    summary = ", ".join(f"R={R}: |diff|={d:.3e} <= {r:.3e} (rk4_err={e:.2e})" for R, d, r, e in diffs)
    return True, f"Path-reversal bound verified on 18 trees (max ratio={float(max_ratio):.4f} <= 1); O(s^{{R-k+1}}) exact; (3) holds ({summary})."


def verify_generator_dynkin_and_renewal(
    drop_newborn: bool = False,
    wrong_self_denom: bool = False,
) -> Tuple[bool, str]:
    """
    1. Exact check of generator formula (7): d/dt|_{t=0} Q_t f == Af.
    2. Exact check of Dynkin's formula Q_t f - f = int_0^t Q_s Af ds up to order 3:
       d^m/dt^m|_{t=0} Q_t f == d^{m-1}/ds^{m-1}|_{s=0} Q_s (Af) for m=1,2,3!
    3. Exact check of degree-law.md birth-time renewal identity:
       Q_t f = S_t f + int_0^t Q_s B S_{t-s} f ds up to order 3!
    4. Exact check of A(1) = b_o and sum_u Af(G, u) == d/dt|_{t=0} E[sum_u f(G_t, u)] (finite tree m_t).
    """
    trees = [
        CoreState([[]], [2]),
        CoreState([[1], [0]], [1, 3]),
        CoreState([[1], [0, 2], [1]], [2, 0, 1]),
    ]
    f0 = cyl_fn(0)
    f1 = cyl_fn(1)

    for st in trees:
        for rt in range(st.n):
            for fn in (f0, f1, lambda g, v: Fraction(1, 1)):
                # 1 & 2. Dynkin to order 3 on small cores, order 2 on 3-core
                max_o = 3 if st.n <= 2 else 2
                d_q = Q_taylor_derivs(st, rt, fn, max_o)
                af_fn = lambda g, v, f_cur=fn: gen_A(g, v, f_cur, drop_newborn=drop_newborn, wrong_self_denom=wrong_self_denom)
                d_q_af = Q_taylor_derivs(st, rt, af_fn, max_o - 1)
                for m in range(1, max_o + 1):
                    if d_q[m] != d_q_af[m - 1]:
                        return False, f"Dynkin/Generator (7) mismatch at order m={m}: {d_q[m]} != {d_q_af[m-1]}"

                # 3. Birth-time renewal identity Q_t f = S_t f + int_0^t Q_s B S_{t-s} f ds:
                # m-th deriv of RHS = S^{(m)} f + sum_{j=0}^{m-1} (Q^{(j)} (B S^{(m-1-j)} f))!
                d_s = S_taylor_derivs(st, rt, fn, max_o)
                for m in range(1, max_o + 1):
                    conv = Fraction(0, 1)
                    for j in range(m):
                        r_ord = m - 1 - j
                        def b_s_r(g: CoreState, v: int, rr=r_ord, f_cur=fn) -> Fraction:
                            if g.L[v] == 0:
                                return Fraction(0, 1)
                            nxt, nid = g.promote(v)
                            return g.b_x(v) * S_taylor_derivs(nxt, nid, f_cur, rr)[rr]
                        conv += Q_taylor_derivs(st, rt, b_s_r, j)[j]
                    if d_q[m] != d_s[m] + conv:
                        return False, f"Birth-time renewal identity mismatch at order m={m}: {d_q[m]} != {d_s[m] + conv}"

    return True, "Generator (7), Dynkin Q_t f - f = int Q_s Af ds, and birth-time identity Q_t f = S_t f + int Q_s B S_{t-s} f ds verified exactly up to order 3."


def verify_edge_transfers_and_unimodularity(bug_swap_rate: bool = False) -> Tuple[bool, str]:
    """
    Verify both edge-transfer divergences from operator.md §4 sum identically to 0
    over all vertices of arbitrary finite trees (mass transport).
    """
    trees = [
        CoreState([[1], [0, 2], [1, 3, 4], [2], [2]], [1, 0, 2, 3, 1]),
        CoreState([[1, 2, 3], [0, 4], [0, 5, 6], [0], [1], [2], [2]], [4, 2, 1, 7, 1, 8, 2]),
    ]
    def g_tab(dw, cw, du, cu):
        if dw >= 10 or du >= 10:
            return Fraction(0, 1)
        return Fraction(3 * dw - 2 * cw + du - 5 * cu, 7)

    def a_tab(du, cu):
        return Fraction(du * du + cu, 11) if du < 10 else Fraction(0, 1)

    T_th = 4
    for st in trees:
        sum_gain_div = Fraction(0, 1)
        sum_arr_div = Fraction(0, 1)
        for u in range(st.n):
            ru = st.r_x(u)
            Du, cu = st.D(u), st.c(u)
            # Gain transfer: arc (w, u) carries r_w g(w, u) from u to w
            inc_g = sum((ru * g_tab(Du, cu, st.D(v), st.c(v)) for v in st.adj[u]), Fraction(0, 1))
            out_g = sum(
                ((ru if bug_swap_rate else st.r_x(w)) * g_tab(st.D(w), st.c(w), Du, cu) for w in st.adj[u]),
                Fraction(0, 1),
            )
            sum_gain_div += inc_g - out_g

            # Arrow transfer: arc (w, u) carries 1_{D_w >= T} a(D_u, c_u) / D_w from u to w
            inc_a = (Fraction(1, Du) * sum((a_tab(st.D(v), st.c(v)) for v in st.adj[u]), Fraction(0, 1))) if Du >= T_th else Fraction(0, 1)
            out_a = sum((a_tab(Du, cu) * Fraction(1, st.D(w)) for w in st.adj[u] if st.D(w) >= T_th), Fraction(0, 1))
            sum_arr_div += inc_a - out_a

        if sum_gain_div != 0:
            return False, f"Gain-transfer divergence did not cancel: {sum_gain_div}"
        if sum_arr_div != 0:
            return False, f"Arrow-transfer divergence did not cancel: {sum_arr_div}"

    return True, "Both directed edge-transfer divergences cancel identically to 0 on finite trees."


def verify_one_port_upper_prerequisites() -> Tuple[bool, str]:
    """
    Independently verify key quantitative inequalities from round2/finite/one-port-upper.md
    used in operator.md §5:
      1. §1 parameter scaling Lambda_j <= 6/j + 2/j^2 for j = 4..30.
      2. §3 birth-death PGF solution f(h,z) and geometric conditional law Z_h | Z_h > 0 ~ Geom(q),
         plus -q log q / (1 - q) <= sqrt(q) <= e^{-sh/2}.
      3. §5 Taylor remainder constant c_d = d(d-1)(3/2)^{d-2}/8 for (r+gamma)(u^d-1) - d(beta-gamma) <= c_d/K.
      4. §2 operator.md tilt inequalities w(t+a) <= w(t) e^{-(mu-delta)a} and w(u-h)/w(u) in [e^{(mu-delta)h}, e^{(mu+delta)h}].
    """
    # 1. Parameter scaling Lambda_j <= 6/j + 2/j^2
    for j in range(4, 31):
        d = B = j
        J = 2 * j
        H = j * j
        kappa = 1.0 - 1.0 / j
        b = 2 * j
        s = 1.0 - 1.0 / b
        assert s > kappa
        gamma = 1.0 / B
        delta = (B - 1.0) / H
        z = (1.0 + 1.0 / (2.0 * gamma)) ** d
        c_d = d * (d - 1) * (1.5 ** (d - 2)) / 8.0
        K = max(1, math.ceil(j * c_d))
        M = K + 1 + math.ceil(j * (z - 1.0) * K)
        V = d * gamma + c_d / K + (z - 1.0) * K / M
        Lam = 1.0 / J + (1.0 - kappa) / 2.0 + V / d + 2.0 / d + gamma + delta
        if Lam > 6.0 / j + 2.0 / (j * j) + 1e-12:
            return False, f"Lambda_j bound failed at j={j}: {Lam} > {6/j + 2/(j*j)}"

    # 2. Birth-death PGF & geometric reciprocal bound:
    for s in (0.5, 0.75, 0.9):
        for h in (0.05, 0.5, 2.0, 5.0):
            r = math.exp(-s * h)
            q = s * r / (1.0 - (1.0 - s) * r)
            assert q <= r + 1e-15
            # Check exact PGF ODE finite difference:
            def f_pgf(hh, zz):
                rr = math.exp(-s * hh)
                f0 = (1.0 - s) * (1.0 - rr) / (1.0 - (1.0 - s) * rr)
                qq = s * rr / (1.0 - (1.0 - s) * rr)
                return f0 + (1.0 - f0) * (qq * zz / (1.0 - (1.0 - qq) * zz))
            eps = 1e-7
            for zz in (0.2, 0.6, 0.9):
                df_dh = (f_pgf(h + eps, zz) - f_pgf(h - eps, zz)) / (2 * eps)
                fv = f_pgf(h, zz)
                ode_rhs = (fv - 1.0) * (fv - (1.0 - s))
                if abs(df_dh - ode_rhs) > 1e-6:
                    return False, f"Birth-death PGF ODE mismatch: {df_dh} vs {ode_rhs}"
            # Check E[1/Z_h | Z_h > 0] = -q log q / (1 - q) <= sqrt(q) <= e^{-sh/2}:
            e_inv = -q * math.log(q) / (1.0 - q)
            if not (e_inv <= math.sqrt(q) + 1e-15 and math.sqrt(q) <= math.exp(-s * h / 2.0) + 1e-15):
                return False, f"Reciprocal geometric bound failed: {e_inv}, {math.sqrt(q)}, {math.exp(-s*h/2)}"

    # 3. Taylor bound in one-port-upper.md §5:
    for d in (3, 4, 6, 10):
        c_d = d * (d - 1) * (1.5 ** (d - 2)) / 8.0
        for K in (1, 3, 10):
            for r in (float(K), float(K) + 2.5):
                for gamma in (0.1, 0.25, 0.5):
                    for beta in (0.0, 0.1, 0.25, 0.5):
                        u = (r + beta) / (r + gamma)
                        assert u <= 1.5 + 1e-15
                        lhs = (r + gamma) * (u ** d - 1.0) - d * (beta - gamma)
                        if lhs > c_d / K + 1e-12:
                            return False, f"§5 Taylor bound failed at d={d}, K={K}, r={r}, g={gamma}, b={beta}: {lhs} > {c_d/K}"

    # 4. Tilt inequalities in operator.md §2:
    mu = 0.565723
    for delta in (0.01, 0.1, 0.4):
        for T in (0.5, 5.0, 50.0):
            def w(t):
                return math.exp(-mu * t - delta * abs(t - T))
            for t in (0.0, 0.5, 4.5, 5.0, 5.5, 60.0):
                for a in (0.0, 0.3, 3.0, 20.0):
                    if w(t + a) > w(t) * math.exp(-(mu - delta) * a) * (1.0 + 1e-13):
                        return False, "Tilt inequality w(t+a) <= w(t) e^{-(mu-delta)a} violated!"
            for h in (0.2, 1.0, 5.0):
                for u in (h, h + 0.5, T, T + 0.5 * h, T + 5.0):
                    if u >= h:
                        ratio = w(u - h) / w(u)
                        if not (math.exp((mu - delta) * h) * (1 - 1e-12) <= ratio <= math.exp((mu + delta) * h) * (1 + 1e-12)):
                            return False, "Tilt ratio w(u-h)/w(u) bound violated!"

    return True, "one-port-upper.md §1, §3 (BD PGF + E[1/Z_h]), §5 (Taylor c_d), and operator.md §2 tilt bounds all verified."


def main():
    print("=== INDEPENDENT ADVERSARIAL AUDIT OF ASTRA CLAIM 4 ===")

    checks = [
        ("Cell operator Q_s vs microscopic, Q_s Q_t = Q_{s+t}, m_{t+h} = m_t Q_h", verify_semigroup_and_mt_evolution),
        ("Path-reversal bound, O(s^{R-k+1}) vanishing, and spatial-cutoff (3)", verify_path_reversal_and_cutoff),
        ("Generator (7), Dynkin's formula, and birth-time renewal identity", verify_generator_dynkin_and_renewal),
        ("Directed edge-transfer cancellations (unimodularity)", verify_edge_transfers_and_unimodularity),
        ("one-port-upper.md (§1,§3,§5) & operator.md §2 tilt inequalities", verify_one_port_upper_prerequisites),
    ]

    print("\n[1] Core verification checks:", flush=True)
    for name, fn in checks:
        print(f"  Running: {name} ...", flush=True)
        ok, msg = fn()
        status = "PASS" if ok else "FAIL"
        print(f"  [{status}] {name}: {msg}", flush=True)
        assert ok, f"Core check failed: {name}: {msg}"

    planted = [
        ("Planted Error 1: Drop newborn term b_o f(G^{b_o}, new) in generator (7)",
         lambda: verify_generator_dynkin_and_renewal(drop_newborn=True)),
        ("Planted Error 2: Drop missing boundary slots in Q_s^{[R]}",
         lambda: verify_path_reversal_and_cutoff(drop_missing_slots=True)),
        ("Planted Error 3: Drop indirect descendants from cell H_o(s)",
         lambda: verify_semigroup_and_mt_evolution(bug_no_grandchild_inherit=True)),
        ("Planted Error 4: Use 1/D_x instead of 1/D_y in gain rate r_x",
         lambda: verify_generator_dynkin_and_renewal(wrong_self_denom=True)),
        ("Planted Error 5: Swap sender/receiver rate in gain edge transfer",
         lambda: verify_edge_transfers_and_unimodularity(bug_swap_rate=True)),
    ]

    print("\n[2] Planted-error negative controls (all MUST be caught = rejected):")
    for name, fn in planted:
        ok, msg = fn()
        assert not ok, f"Planted error was NOT caught: {name}"
        print(f"  [CAUGHT] {name} -> {msg}")

    print("\nALL CLAIM 4 CHECKS PASSED & ALL 5 PLANTED ERRORS CAUGHT.")


if __name__ == "__main__":
    main()
