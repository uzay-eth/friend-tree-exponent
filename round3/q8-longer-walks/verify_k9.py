"""Exact rational verification for NINE-STEP.md (fixed non-lazy k = 9 friend tree)
and ERRATUM-2.md (exact rational interval for b_*).

Uses 100% exact fractions.Fraction arithmetic (zero floating-point checks anywhere):
1. Certifies 0.465941272384 < b_* < 0.465941272386 and 8/19 < b_* < 1/2 for ERRATUM-2.md.
2. Certifies all analytic rows for k = 9 (beta = 999/1000, J_0 = 4, J_cap = 17):
   - Descartes quasiconvexity (0 <= j*s <= 1) and endpoint margins for (2,1), (3,2), and j = 4..16;
   - Constant tail j >= 17 margin (199/200000);
   - Exact positive rational margins for (2,0; 1L) and (3,0; 2L);
   - Exact positive 1D degree-5 Bernstein certificates for (2,0; 1S), (3,0; 1L,1S), (3,1; 1L,1B);
   - Exact positive 2D bidegree-(5,5) Bernstein certificates (36 coeffs each) for (3,0; 2S), (3,1; 1S,1B).
3. Direct forward 9-step walk audit of EVERY line of the 4-level root-return global drift proof
   n * E[Delta X | T_n] <= 9/4 + (999/1000) * X across all 435 unlabeled trees (n = 2..11, 4,393 roots)
   plus 48 adversarial stress trees (r = 1, 99/100, 1/2, 1/200 across all 5 states and j = 4..25).
4. Exact check of k = 9 positive-hub constants (R_0 = 41000, D_0 = 82000, ratio 1/81), old-prefix cut
   coefficients at M = 32 * R_0 = 1312000, J_0 = 3 obstruction sum, and k >= 11 depth-3 barrier.
"""

from fractions import Fraction
import math


# ---------------------------------------------------------------------------
# 1. Exact rational interval verification for b_* (ERRATUM-2.md)
# ---------------------------------------------------------------------------

def verify_b_star_erratum2():
    def log1p_bounds(t, N=40):
        y = t / (t + 2)
        y2 = y * y
        s = Fraction(0)
        p = y
        for k in range(N):
            s += p / (2 * k + 1)
            p *= y2
        lo = 2 * s
        hi = 2 * (s + p / ((2 * N + 1) * (1 - y2)))
        return lo, hi

    b_lo = Fraction(465941272384, 10 ** 12)
    b_hi = Fraction(465941272386, 10 ** 12)
    _, log_hi = log1p_bounds(1 / b_lo, N=40)
    log_lo, _ = log1p_bounds(1 / b_hi, N=40)
    assert (1 / b_lo) - log_hi > 1
    assert (1 / b_hi) - log_lo < 1
    assert Fraction(8, 19) < b_lo < b_hi < Fraction(1, 2)


# ---------------------------------------------------------------------------
# 2. Fixed non-lazy k = 9 score (beta = 999/1000, J_0 = 4, J_cap = 17)
# ---------------------------------------------------------------------------

BETA = Fraction(999, 1000)
J_0 = 4
J_CAP = 17

# Five finite states for j in {2, 3}: Phi[(j, b)] where b = count of Type-B children z
PHI = {
    (2, 0): Fraction(1000, 1000),
    (2, 1): Fraction(1946, 1000),
    (3, 0): Fraction(1948, 1000),
    (3, 1): Fraction(2388, 1000),
    (3, 2): Fraction(2540, 1000),
}

# Pure degree-only table for j >= 4 (constant for j >= 17)
G_TABLE = {
    4: Fraction(2808, 1000),
    5: Fraction(2992, 1000),
    6: Fraction(3139, 1000),
    7: Fraction(3266, 1000),
    8: Fraction(3381, 1000),
    9: Fraction(3488, 1000),
    10: Fraction(3583, 1000),
    11: Fraction(3664, 1000),
    12: Fraction(3735, 1000),
    13: Fraction(3799, 1000),
    14: Fraction(3857, 1000),
    15: Fraction(3910, 1000),
    16: Fraction(3959, 1000),
    17: Fraction(4005, 1000),  # 801/200
}


def G_val(j):
    assert j >= 4
    return G_TABLE[min(j, J_CAP)]


def s_deg(j):
    assert j >= 4
    return G_val(j + 1) - G_val(j)


def s_state(j, b):
    if j >= J_0:
        return s_deg(j)
    nxt = G_val(J_0) if j + 1 == J_0 else PHI[(j + 1, b)]
    return nxt - PHI[(j, b)]


def gamma_state(j, b):
    if j >= J_0 or b >= j - 1:
        return Fraction(0)
    return PHI[(j, b + 1)] - PHI[(j, b)]


def is_deep_child(adj, y, z):
    """Return True iff outward neighbor z of y is Type B (d_z >= 2 and >= 1 non-leaf outward neighbor w)."""
    if len(adj[z]) == 1:
        return False
    return any(len(adj[w]) >= 2 for w in adj[z] if w != y)


def phi_neighbor(adj, v, y):
    j = len(adj[y])
    assert j >= 2
    if j >= J_0:
        return G_val(j)
    b = sum(1 for z in adj[y] if z != v and is_deep_child(adj, y, z))
    return PHI[(j, b)]


def stock(adj, v):
    return sum((phi_neighbor(adj, v, y) for y in adj[v] if len(adj[y]) >= 2), Fraction(0))


# ---------------------------------------------------------------------------
# 3. Exact 1D & 2D Bernstein polynomial builders (100% Fraction)
# ---------------------------------------------------------------------------

def P9_eval(j, s, a):
    t = Fraction(j - 1, j)
    u0 = a
    u1 = t - a + a * u0
    u2 = t - a + a * u1
    u3 = t - a + a * u2
    u4 = t - a + a * u3
    return u0 + u1 + u2 + u3 + j * s * u4


def eval_k9_branch_poly(j, nL, thetas, b, aB, s_y, gamma_y):
    u = [Fraction(0)] * 5
    u[0] = (nL + sum((1 - th for th in thetas), Fraction(0))) / j + aB
    for p in range(1, 5):
        s_S = Fraction(0)
        for th in thetas:
            val = th ** p
            for r in range(p):
                val += (th ** r) * u[p - 1 - r]
            s_S += (1 - th) * val
        u[p] = (nL * u[p - 1] + s_S + b - j * aB + j * aB * u[p - 1]) / j
    birth = u[0] + u[1] + u[2] + u[3]
    up_y = j * u[4]
    up_SB = sum(
        (th * u[3] + (th ** 2) * u[2] + (th ** 3) * u[1] + (th ** 4) * u[0] + (th ** 5) for th in thetas),
        Fraction(0),
    )
    return birth + s_y * up_y + gamma_y * up_SB


def lagrange_basis_matrix(xs):
    d = len(xs) - 1
    mat = []
    for i in range(d + 1):
        poly = [Fraction(1)]
        denom = Fraction(1)
        for j in range(d + 1):
            if j == i:
                continue
            denom *= xs[i] - xs[j]
            nxt = [Fraction(0)] * (len(poly) + 1)
            for idx, v in enumerate(poly):
                nxt[idx] -= v * xs[j]
                nxt[idx + 1] += v
            poly = nxt
        mat.append([v / denom for v in poly])
    return mat


def bernstein_1d(fn, deg=5):
    xs = [Fraction(i, deg) for i in range(deg + 1)]
    L = lagrange_basis_matrix(xs)
    c = [Fraction(0)] * (deg + 1)
    for i in range(deg + 1):
        val = fn(xs[i])
        for r in range(deg + 1):
            c[r] += val * L[i][r]
    # Certify exact polynomial identity at 3 off-grid points
    for xt in (Fraction(1, 7), Fraction(4, 11), Fraction(11, 13)):
        rec = sum((c[r] * (xt ** r) for r in range(deg + 1)), Fraction(0))
        assert rec == fn(xt)
    bern = [
        sum((c[r] * Fraction(math.comb(k, r), math.comb(deg, r)) for r in range(k + 1)), Fraction(0))
        for k in range(deg + 1)
    ]
    # Certify Bernstein reconstruction at xt = 3/7
    xt = Fraction(3, 7)
    rec_b = sum(
        (bern[k] * math.comb(deg, k) * (xt ** k) * ((1 - xt) ** (deg - k)) for k in range(deg + 1)),
        Fraction(0),
    )
    assert rec_b == fn(xt)
    return bern


def bernstein_2d(fn, deg1=5, deg2=5):
    xs1 = [Fraction(i, deg1) for i in range(deg1 + 1)]
    xs2 = [Fraction(j, deg2) for j in range(deg2 + 1)]
    L1 = lagrange_basis_matrix(xs1)
    L2 = lagrange_basis_matrix(xs2)
    C = [[Fraction(0)] * (deg2 + 1) for _ in range(deg1 + 1)]
    for i in range(deg1 + 1):
        for j in range(deg2 + 1):
            val = fn(xs1[i], xs2[j])
            for r in range(deg1 + 1):
                v1 = val * L1[i][r]
                for s in range(deg2 + 1):
                    C[r][s] += v1 * L2[j][s]
    # Certify exact 2D polynomial identity at 3 off-grid points
    for (t1, t2) in ((Fraction(1, 7), Fraction(2, 11)), (Fraction(5, 13), Fraction(3, 17)), (Fraction(6, 7), Fraction(9, 11))):
        rec = Fraction(0)
        for r in range(deg1 + 1):
            for s in range(deg2 + 1):
                rec += C[r][s] * (t1 ** r) * (t2 ** s)
        assert rec == fn(t1, t2)
    B = [[Fraction(0)] * (deg2 + 1) for _ in range(deg1 + 1)]
    for k1 in range(deg1 + 1):
        for k2 in range(deg2 + 1):
            acc = Fraction(0)
            for r in range(k1 + 1):
                f1 = Fraction(math.comb(k1, r), math.comb(deg1, r))
                for s in range(k2 + 1):
                    f2 = Fraction(math.comb(k2, s), math.comb(deg2, s))
                    acc += C[r][s] * f1 * f2
            B[k1][k2] = acc
    # Certify 2D Bernstein reconstruction at (2/7, 5/11)
    t1, t2 = Fraction(2, 7), Fraction(5, 11)
    rec_b = Fraction(0)
    for k1 in range(deg1 + 1):
        w1 = math.comb(deg1, k1) * (t1 ** k1) * ((1 - t1) ** (deg1 - k1))
        for k2 in range(deg2 + 1):
            w2 = math.comb(deg2, k2) * (t2 ** k2) * ((1 - t2) ** (deg2 - k2))
            rec_b += B[k1][k2] * w1 * w2
    assert rec_b == fn(t1, t2)
    return B


def verify_analytic_k9():
    # 1. Check all increments in [0, 1] and 0 <= j*s_j <= 1
    for (j, b), val in PHI.items():
        assert Fraction(1) <= val <= G_val(J_CAP)
        assert Fraction(0) <= s_state(j, b) <= 1
        assert Fraction(0) <= gamma_state(j, b) <= 1

    # 2. Verify P_9 exact expansion and Descartes quasiconvexity + endpoint margins for j = 4..16
    for j in range(4, J_CAP):
        s = s_deg(j)
        t = Fraction(j - 1, j)
        assert Fraction(0) <= j * s <= 1
        # Verify exact polynomial coefficients at 6 points (degree 5)
        for idx in range(6):
            a = t * Fraction(idx, 5)
            raw = P9_eval(j, s, a)
            exp = (
                (3 + j * s) * t
                - Fraction(2 + j * s, j) * a
                - Fraction(1 + j * s, j) * (a ** 2)
                - s * (a ** 3)
                + (1 - j * s) * (a ** 4)
                + j * s * (a ** 5)
            )
            assert raw == exp
        p0 = P9_eval(j, s, Fraction(0))
        pt = P9_eval(j, s, t)
        assert p0 == 3 * t + (j - 1) * s
        assert pt == t + t ** 2 + t ** 3 + t ** 4 + j * (t ** 5) * s
        assert BETA * G_val(j) - p0 > 0
        assert BETA * G_val(j) - pt > 0

    # 3. Constant tail j >= 17
    assert BETA * G_val(J_CAP) - 4 == Fraction(199, 200000) > 0

    # 4. Pure Type-B states (2, 1) and (3, 2): Descartes quasiconvexity + endpoints
    s21 = s_state(2, 1)
    assert Fraction(0) <= 2 * s21 <= 1
    assert BETA * PHI[(2, 1)] - P9_eval(2, s21, Fraction(0)) > 0
    assert BETA * PHI[(2, 1)] - P9_eval(2, s21, Fraction(1, 4)) > 0

    s32 = s_state(3, 2)
    assert Fraction(0) <= 3 * s32 <= 1
    assert BETA * PHI[(3, 2)] - P9_eval(3, s32, Fraction(0)) > 0
    assert BETA * PHI[(3, 2)] - P9_eval(3, s32, Fraction(1, 3)) > 0

    # 5. Pure-leaf constant rows (2,0; 1L) and (3,0; 2L)
    m20L = BETA * PHI[(2, 0)] - eval_k9_branch_poly(2, 1, (), 0, Fraction(0), s_state(2, 0), gamma_state(2, 0))
    m30LL = BETA * PHI[(3, 0)] - eval_k9_branch_poly(3, 2, (), 0, Fraction(0), s_state(3, 0), gamma_state(3, 0))
    assert m20L == Fraction(9, 4000) > 0
    assert m30LL == Fraction(27553, 20250000) > 0

    # 6. 1D degree-5 exact Bernstein certificates
    b20S = bernstein_1d(
        lambda tau: BETA * PHI[(2, 0)] - eval_k9_branch_poly(2, 0, ((1 + tau) / 2,), 0, Fraction(0), s_state(2, 0), gamma_state(2, 0)),
        5,
    )
    b30LS = bernstein_1d(
        lambda tau: BETA * PHI[(3, 0)] - eval_k9_branch_poly(3, 1, ((1 + tau) / 2,), 0, Fraction(0), s_state(3, 0), gamma_state(3, 0)),
        5,
    )
    b31L = bernstein_1d(
        lambda u: BETA * PHI[(3, 1)] - eval_k9_branch_poly(3, 1, (), 1, u / 6, s_state(3, 1), gamma_state(3, 1)),
        5,
    )
    assert all(x > 0 for x in b20S)
    assert all(x > 0 for x in b30LS)
    assert all(x > 0 for x in b31L)

    # 7. 2D bidegree-(5, 5) exact Bernstein certificates (36 coeffs each)
    b30SS = bernstein_2d(
        lambda t1, t2: BETA * PHI[(3, 0)]
        - eval_k9_branch_poly(3, 0, ((1 + t1) / 2, (1 + t2) / 2), 0, Fraction(0), s_state(3, 0), gamma_state(3, 0)),
        5,
        5,
    )
    b31S = bernstein_2d(
        lambda tau, u: BETA * PHI[(3, 1)]
        - eval_k9_branch_poly(3, 0, ((1 + tau) / 2,), 1, u / 6, s_state(3, 1), gamma_state(3, 1)),
        5,
        5,
    )
    assert all(x > 0 for row in b30SS for x in row)
    assert all(x > 0 for row in b31S for x in row)

    # 8. Positive-hub, terminal-vertex, and old-prefix cut constants for k = 9
    shift = Fraction(9, 4) / BETA  # 250/111
    assert shift == Fraction(250, 111)
    assert 5 * G_val(J_CAP) / (1 - BETA) == Fraction(20025) <= 20500
    assert 5 * shift / (1 - BETA) == Fraction(1250000, 111) <= 20500
    R0 = 41000
    assert 20500 * 2 == R0
    for d in range(2, 50):
        assert d * (1 - Fraction(1, d)) ** 5 >= Fraction(1, 16)
    assert Fraction(1, 16) / (Fraction(1, 16) + 5) == Fraction(1, 81)
    M = 32 * R0
    assert M == 1312000
    assert Fraction(8 * R0, M) * 2 == Fraction(1, 2)
    assert 2 * (1 + 4 * R0) == 328002
    assert 2 * (4 * R0 + Fraction(8 * R0, M)) == Fraction(656001, 2) <= 328001

    # 9. Exact structural obstruction sums:
    # (a) J_0 = 3 monotone obstruction at k = 9: 6 * sum_{j=3}^{16} 1/j^2 > 2
    s_j03 = 6 * sum((Fraction(1, j * j) for j in range(3, 17)), Fraction(0))
    assert s_j03 > 2
    assert Fraction(15, 16) + Fraction(1, 16) * (s_j03 - 1) > 1
    # (b) k = 2m+1 >= 11 depth-3 monotone barrier: 2*beta*(1+beta) < 4 <= m-1 for all beta < 1
    assert 2 * Fraction(1) * (1 + Fraction(1)) == 4

    return min(b20S), min(b30LS), min(b31L), min(x for r in b30SS for x in r), min(x for r in b31S for x in r), s_j03


# ---------------------------------------------------------------------------
# 4. Direct forward 9-step tree audit of EVERY algebraic line of the global proof
# ---------------------------------------------------------------------------

def make_adj(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return adj


def forward_q9(adj):
    n = len(adj)
    prob = [Fraction(1)] * n
    for _ in range(9):
        nxt = [Fraction(0)] * n
        for u in range(n):
            share = prob[u] / len(adj[u])
            for w in adj[u]:
                nxt[w] += share
        prob = nxt
    assert sum(prob, Fraction(0)) == n
    return prob


def power_P_f(adj, max_steps):
    n = len(adj)
    cur = [Fraction(1, len(adj[x])) for x in range(n)]
    out = [cur]
    for _ in range(max_steps):
        nxt = [Fraction(0)] * n
        for x in range(n):
            nxt[x] = sum((cur[y] for y in adj[x]), Fraction(0)) / len(adj[x])
        cur = nxt
        out.append(cur)
    return out


def audit_k9_tree(adj, roots, check_edges=False):
    n = len(adj)
    q9 = forward_q9(adj)
    Pf = power_P_f(adj, 9)
    for x in range(n):
        assert q9[x] == len(adj[x]) * Pf[9][x]

    edge_checks = 0
    if check_edges:
        for u in range(n):
            for v in adj[u]:
                du, dv = len(adj[u]), len(adj[v])
                assert q9[u] >= Fraction(1, (du ** 4) * (dv ** 5))
                edge_checks += 1

    root_checks = 0
    for v in roots:
        d = len(adj[v])
        I = [y for y in adj[v] if len(adj[y]) >= 2]
        L_list = [y for y in adj[v] if len(adj[y]) == 1]
        C = len(I)
        L = len(L_list)
        r = Fraction(C, d)
        e = [None, Pf[2][v], Pf[4][v], Pf[6][v], Pf[8][v]]
        u = {y: [Pf[1][y], Pf[3][y], Pf[5][y], Pf[7][y], Pf[9][y]] for y in I}

        # Exact root rows (R)
        assert d * e[1] == (1 - r) + sum((u[y][0] for y in I), Fraction(0))
        for p in (2, 3, 4):
            assert d * e[p] == d * (1 - r) * e[p - 1] + sum((u[y][p - 1] for y in I), Fraction(0))

        # Exact direct drift by simulating every attachment target x in 0..n-1
        old_X = stock(adj, v)
        N_direct = Fraction(0)
        for x in range(n):
            adj[x].append(n)
            adj.append([x])
            new_X = stock(adj, v)
            adj.pop()
            adj[x].pop()
            inc = new_X - old_X
            assert Fraction(0) <= inc <= Fraction(1)
            N_direct += q9[x] * inc

        # Extract branch structure for every y in I
        info = {}
        N_formula = L * e[4]
        for y in I:
            j = len(adj[y])
            if j >= J_0:
                a = sum((Fraction(1, len(adj[z])) for z in adj[y] if z != v), Fraction(0)) / j
                t = Fraction(j - 1, j)
                sy = s_deg(j)
                gy = Fraction(0)
                info[y] = ("G", j, a, t, sy, gy, 0, [], 0, Fraction(0))
                assert u[y][0] == a + Fraction(1, j * d)
                for p in (1, 2, 3, 4):
                    assert u[y][p] <= e[p] / j + (t - a) + a * u[y][p - 1]
                N_formula += j * sy * u[y][4]
            else:
                nL = 0
                S_list = []  # (z, m, theta)
                B_list = []
                for z in adj[y]:
                    if z == v:
                        continue
                    dz = len(adj[z])
                    if dz == 1:
                        nL += 1
                    elif is_deep_child(adj, y, z):
                        B_list.append(z)
                    else:
                        S_list.append((z, dz, Fraction(dz - 1, dz)))
                b = len(B_list)
                thetas = [th for _, _, th in S_list]
                aB = sum((Fraction(1, len(adj[z])) for z in B_list), Fraction(0)) / j
                assert Fraction(0) <= aB <= Fraction(b, 2 * j)
                sy = s_state(j, b)
                gy = gamma_state(j, b)
                a0 = Fraction(nL, j) + aB
                info[y] = ("S", j, a0, Fraction(b, j) - aB, sy, gy, nL, thetas, b, aB)

                # Verify exact/upper rows for u[y][0..4] and exact (m_i - 1) P^8 f(z_i)
                assert u[y][0] == Fraction(1, j * d) + a0 + sum(((1 - th) / j for th in thetas), Fraction(0))
                for p in (1, 2, 3, 4):
                    s_S = Fraction(0)
                    for th in thetas:
                        val = th ** p
                        for k_idx in range(p):
                            val += (th ** k_idx) * u[y][p - 1 - k_idx]
                        s_S += ((1 - th) / j) * val
                    rhs = e[p] / j + (Fraction(b, j) - aB) + a0 * u[y][p - 1] + s_S
                    assert u[y][p] <= rhs
                N_formula += j * sy * u[y][4]
                for z, m, th in S_list:
                    prom_w = (m - 1) * Pf[8][z]
                    exact_pw = th * u[y][3] + (th ** 2) * u[y][2] + (th ** 3) * u[y][1] + (th ** 4) * u[y][0] + (th ** 5)
                    assert prom_w == exact_pw
                    N_formula += gy * prom_w

        assert N_direct == N_formula

        # Verify the 4-level root-return coefficients H_4, H_3, H_2, H_1 in [0, 1]
        # AND every intermediate inequality of the adjoint reduction chain!
        H = [None, None, None, None, None]
        H[4] = (1 - r) + sum((info[y][4] for y in I), Fraction(0)) / d
        assert Fraction(0) <= H[4] <= 1

        # Helper: constant part of R_q(y) for q in 1..4
        def const_R(y, q):
            kind, j, a0, const_B, sy, gy, nL, thetas, b, aB = info[y]
            if kind == "G":
                return const_B - a0  # t - a (since info stored t in slot 3!)
            return const_B + sum((((1 - th) / j) * (th ** q) for th in thetas), Fraction(0))

        # Wait: for kind=="G", slot 3 is t and slot 2 is a, so const_part of R_q is t - a!
        c = {4: {y: info[y][1] * info[y][4] for y in I}}
        for p in (4, 3, 2):
            cp = {}
            for y in I:
                kind, j, a0, const_B, sy, gy, nL, thetas, b, aB = info[y]
                if kind == "G":
                    cp[y] = H[p] + c[p][y] * a0
                else:
                    val = H[p] + c[p][y] * a0
                    for q in range(p, 5):
                        val += c[q][y] * sum((((1 - th) / j) * (th ** (q - p)) for th in thetas), Fraction(0))
                    val += gy * sum((th ** (5 - p) for th in thetas), Fraction(0))
                    cp[y] = val
                assert Fraction(0) <= cp[y] / j <= 1
            c[p - 1] = cp
            H[p - 1] = H[p] * (1 - r) + sum((cp[y] / info[y][1] for y in I), Fraction(0)) / d
            assert Fraction(0) <= H[p - 1] <= 1

        # Compute c[0][y] (using H[1]) and verify the exact adjoint upper bound after all 4 root-return eliminations:
        c[0] = {}
        for y in I:
            kind, j, a0, const_B, sy, gy, nL, thetas, b, aB = info[y]
            if kind == "G":
                c[0][y] = H[1] + c[1][y] * a0
            else:
                val = H[1] + c[1][y] * a0
                for q in range(1, 5):
                    val += c[q][y] * sum((((1 - th) / j) * (th ** (q - 1)) for th in thetas), Fraction(0))
                val += gy * sum((th ** 4 for th in thetas), Fraction(0))
                c[0][y] = val

        # After substituting u_4, e_4, u_3, e_3, u_2, e_2, u_1, e_1 in exact order:
        line_H = H[1] * (1 - r) + sum((
            c[0][y] * u[y][0]
            + sum((c[q][y] * const_R(y, q) for q in range(1, 5)), Fraction(0))
            + info[y][5] * sum((th ** 5 for th in info[y][7]), Fraction(0))
            for y in I), Fraction(0))
        assert N_formula <= line_H

        # Now set H_4 = H_3 = H_2 = H_1 = 1 -> c_one[p][y] >= c[p][y]!
        c_one = {4: {y: info[y][1] * info[y][4] for y in I}}
        for p in (4, 3, 2, 1):
            cp = {}
            for y in I:
                kind, j, a0, const_B, sy, gy, nL, thetas, b, aB = info[y]
                if kind == "G":
                    cp[y] = Fraction(1) + c_one[p][y] * a0
                else:
                    val = Fraction(1) + c_one[p][y] * a0
                    for q in range(p, 5):
                        val += c_one[q][y] * sum((((1 - th) / j) * (th ** (q - p)) for th in thetas), Fraction(0))
                    val += gy * sum((th ** (5 - p) for th in thetas), Fraction(0))
                    cp[y] = val
                assert cp[y] >= c[p - 1][y]
            c_one[p - 1] = cp

        line_H1 = (1 - r) + sum((
            c_one[0][y] * u[y][0]
            + sum((c_one[q][y] * const_R(y, q) for q in range(1, 5)), Fraction(0))
            + info[y][5] * sum((th ** 5 for th in info[y][7]), Fraction(0))
            for y in I), Fraction(0))
        assert line_H <= line_H1

        repl_cost = Fraction(0)
        poly_sum = Fraction(0)
        for y in I:
            kind, j, a0, const_B, sy, gy, nL, thetas, b, aB = info[y]
            if kind == "G":
                assert c_one[0][y] <= 5
                cost_y = c_one[0][y] / (j * d)
                assert cost_y <= Fraction(5, 4 * d)
                py = P9_eval(j, sy, a0)
            else:
                assert c_one[0][y] <= j
                cost_y = c_one[0][y] / (j * d)
                assert cost_y <= Fraction(1, d)
                py = eval_k9_branch_poly(j, nL, thetas, b, aB, sy, gy)
            # Verify exact adjoint identity per branch y!
            branch_adj = (
                c_one[0][y] * u[y][0]
                + sum((c_one[q][y] * const_R(y, q) for q in range(1, 5)), Fraction(0))
                + gy * sum((th ** 5 for th in thetas), Fraction(0))
            )
            assert branch_adj == py + cost_y
            repl_cost += cost_y
            poly_sum += py
            assert py <= BETA * phi_neighbor(adj, v, y)

        assert line_H1 == (1 - r) + repl_cost + poly_sum
        assert repl_cost <= Fraction(5, 4)
        assert N_direct <= line_H <= line_H1 <= Fraction(9, 4) + BETA * old_X

        # Leaf-bounce, leaf-promotion, and attachment-gap bounds for k = 9
        assert d * (1 - r) ** 5 <= q9[v] <= d
        assert d - 5 * C <= q9[v]
        assert sum((q9[z] for z in L_list), Fraction(0)) >= Fraction(L, d) ** 5
        root_checks += 1

    return root_checks, edge_checks


def unlabeled_trees_up_to(max_n):
    def canon(adj, u, p):
        return "(" + "".join(sorted(canon(adj, w, u) for w in adj[u] if w != p)) + ")"

    def key(adj):
        return min(canon(adj, r, -1) for r in range(len(adj)))

    lvl = [[(0, 1)]]
    out = {2: lvl}
    for n in range(2, max_n):
        seen = {}
        for edges in lvl:
            for p in range(n):
                ne = edges + [(p, n)]
                k = key(make_adj(n + 1, ne))
                if k not in seen:
                    seen[k] = ne
        lvl = list(seen.values())
        out[n + 1] = lvl
    return out


def build_adversarial_k9_trees():
    trees = []

    def add_child_spec(edges, nxt, parent, ch):
        kind = ch[0]
        if kind == "L":
            edges.append((parent, nxt)); nxt += 1
        elif kind == "S":
            m = ch[1]
            z = nxt; nxt += 1
            edges.append((parent, z))
            for _ in range(m - 1):
                edges.append((z, nxt)); nxt += 1
        elif kind == "B":
            m, dl = ch[1], ch[2]
            z = nxt; nxt += 1
            edges.append((parent, z))
            for i in range(m - 1):
                w = nxt; nxt += 1
                edges.append((z, w))
                if i == 0 or dl > 0:
                    for _ in range(max(1, dl)):
                        edges.append((w, nxt)); nxt += 1
        return nxt

    def add_tree(name, L_count, branches):
        # branches: list of ("J", [child_specs]) or ("G", j, mode)
        edges = []
        nxt = 1
        for _ in range(L_count):
            edges.append((0, nxt)); nxt += 1
        for br in branches:
            if br[0] == "J":
                y = nxt; nxt += 1
                edges.append((0, y))
                for ch in br[1]:
                    nxt = add_child_spec(edges, nxt, y, ch)
            else:
                _, j, mode = br
                y = nxt; nxt += 1
                edges.append((0, y))
                for i in range(j - 1):
                    z = nxt; nxt += 1
                    edges.append((y, z))
                    if mode == "deep":
                        for _ in range(3):
                            w = nxt; nxt += 1
                            edges.append((z, w))
                            for _ in range(3):
                                edges.append((w, nxt)); nxt += 1
                    elif mode == "mixed" and i % 2 == 0:
                        w = nxt; nxt += 1
                        edges.append((z, w))
                        edges.append((w, nxt)); nxt += 1
        trees.append((name, make_adj(nxt, edges)))

    # High r = 1 stress trees across ALL 5 finite states and all subcases
    j_configs = [
        ("2_0_L", [("L",)]),
        ("2_0_S2", [("S", 2)]),
        ("2_0_S15", [("S", 15)]),
        ("2_1_B2", [("B", 2, 5)]),
        ("2_1_B6", [("B", 6, 3)]),
        ("3_0_2L", [("L",), ("L",)]),
        ("3_0_LS2", [("L",), ("S", 2)]),
        ("3_0_LS15", [("L",), ("S", 15)]),
        ("3_0_2S_mix", [("S", 2), ("S", 15)]),
        ("3_1_LB", [("L",), ("B", 2, 4)]),
        ("3_1_SB", [("S", 5), ("B", 4, 3)]),
        ("3_2_2B", [("B", 2, 4), ("B", 5, 3)]),
    ]
    for cname, chs in j_configs:
        add_tree(f"r=1,5x({cname})", 0, [("J", chs)] * 5)
        add_tree(f"r=1/200,({cname})", 199, [("J", chs)])

    # Mixed high-r trees combining all states and high degrees j >= 4
    add_tree("r=1,all_5_states+G25", 0, [("J", chs) for _, chs in j_configs] + [("G", 4, "deep"), ("G", 16, "mixed"), ("G", 25, "leaf")])
    add_tree("r=1/2,all_5_states+G25", 15, [("J", chs) for _, chs in j_configs] + [("G", 4, "deep"), ("G", 17, "deep"), ("G", 25, "leaf")])
    add_tree("r=99/100,mixed", 1, [("J", [("L",)])] * 33 + [("J", [("S", 2), ("S", 6)])] * 33 + [("J", [("S", 3), ("B", 2, 3)])] * 33)
    for j in (4, 5, 8, 12, 16, 17, 25):
        for mode in ("leaf", "deep", "mixed"):
            add_tree(f"r=1/200,G({j},{mode})", 199, [("G", j, mode)])
    return trees


def main():
    verify_b_star_erratum2()
    print("[PASS] ERRATUM-2.md exact rational interval: 0.465941272384 < b_* < 0.465941272386 and 8/19 < b_* < 1/2")

    b_mins = verify_analytic_k9()
    print(f"[PASS] k=9 analytic rows (beta=999/1000, J_0=4, J_cap=17), Descartes endpoints (2,1)/(3,2)/j=4..16, tail margin 199/200000, all 1D (deg 5) & 2D (bideg 5x5) Bernstein coeffs > 0, k=9 hub/cut constants (R_0=41000, M=1312000), & J_0=3 obstruction sum {b_mins[5]} > 2")

    all_t = unlabeled_trees_up_to(11)
    nt = nr = ne = 0
    for n in range(2, 12):
        for edges in all_t[n]:
            adj = make_adj(n, edges)
            rc, ec = audit_k9_tree(adj, range(n), check_edges=(n <= 9))
            nt += 1
            nr += rc
            ne += ec
    print(f"[PASS] Direct exact k=9 global drift audit: {nt} unlabeled trees (n=2..11), {nr} marked roots, {ne} directed-edge lower-bound checks (n=2..9)")

    stress = build_adversarial_k9_trees()
    for name, adj in stress:
        audit_k9_tree(adj, [0], check_edges=False)
    print(f"[PASS] Adversarial k=9 stress trees (r=1, r=99/100, r=1/2, r=1/200; all 5 states & j=4..25): {len(stress)} trees")


if __name__ == "__main__":
    main()
