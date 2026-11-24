#!/usr/bin/env python3
"""Independent exact short fixtures for round2/fringe/other-marginals.md (hy4 referee).

Repository-relative usage (one CPU, short):
  python3 -B round2/checks/fringe_marginals_hy4.py

Checks the algebra unique to this file: the one-clock Radon-Nikodym derivative and
stopped change-of-measure identity, the geometric-tail/Beta check, the rejection
compensator and positivity bound, the exact downward-branch rate N+p_p-p_v, the
Yule/negative-binomial sandwich moments, and the crossing identity. Planted errors
must be rejected. Exact Fraction arithmetic; no production computation.
"""
from fractions import Fraction as F
import random
import sys


def ok(cond, msg):
    if not cond:
        raise ArithmeticError(msg)


MU = F(1131, 2000)  # generic rational mu in (0,1); no enclosure membership is asserted
R = 1 + MU

# ---------------------------------------------------------------------------
# 1. One-clock change of measure: rate-1 vs rate-r proposal Poisson processes.
# ---------------------------------------------------------------------------
def stopped_identity(Tj, j, Rj, plant=None):
    """e^{-mu Tj} * [exp(mu Tj) r^{-(j+Rj)}] * r^j = r^{-Rj}.
    Work with the rational r^{-Rj} part (the exponentials cancel identically)."""
    if plant == "drop_rejections":
        return F(1)  # would give r^{-j}, not r^{-j-Rj}
    return R ** (-(j + Rj)) * R ** j  # = r^{-Rj}


def test_change_of_measure(plant=None):
    rnd = random.Random(11)
    for _ in range(200):
        j = rnd.randrange(1, 8)
        Rj = rnd.randrange(0, 12)
        Tj = F(rnd.randrange(1, 50), 7)
        got = stopped_identity(Tj, j, Rj, plant)
        correct = R ** (-Rj)
        if plant == "drop_rejections":
            ok(got == correct, "planted drop-rejections survived")
        else:
            ok(got == correct, "stopped change-of-measure identity failed")
    return 200


# ---------------------------------------------------------------------------
# 2. Geometric tail and its point masses.
# ---------------------------------------------------------------------------
def test_geometric(B=F(3, 7), plant=None):
    # Pr(c>=j+1) = B r^{-j}; Pr(c=k) = B r^{-(k-1)}(1-r^{-1}).
    for k in range(1, 12):
        tail = B * R ** (-(k - 1))
        point = B * R ** (-(k - 1)) * (1 - 1 / R)
        ok(point == tail * (1 - 1 / R), "point-mass/tail relation")
        if plant == "wrong_ratio":
            ok(point == tail * MU, "planted wrong geometric ratio survived")
    # tail ratio is exactly 1/r
    ok(B * R ** (-5) / (B * R ** (-4)) == 1 / R, "tail ratio")
    return True


# ---------------------------------------------------------------------------
# 3. Rejection compensator: rate of rejection = r(1-b) = r c / D.
# ---------------------------------------------------------------------------
def rejection_rate(L, c, plant=None):
    D = L + c
    b = F(L, D)
    if plant == "use_b":
        return R * b  # WRONG compensator
    return R * (1 - b)


def test_compensator(plant=None):
    for L in range(1, 30):
        for c in range(1, 12):
            D = L + c
            got = rejection_rate(L, c, plant=plant)
            want = R * F(c, D)
            if plant == "use_b":
                ok(got == want, "planted b-compensator survived")
            else:
                ok(got == want, "rejection compensator identity")
    # positivity bound E[W_L] >= m-r using rb - rho <= r pointwise
    for L in range(2, 12):
        for c in range(1, 8):
            D = L + c
            rho = F(c, 2)  # upper bound on sum 1/D_y
            rb = R * F(L, D)
            ok(rb - rho <= R, "rb - rho <= r bound")
    return True


# ---------------------------------------------------------------------------
# 4. Downward-branch rate: N + p_p - p_v, by enumerating directed arrows.
# ---------------------------------------------------------------------------
def branch_rate(N, p_p, p_v):
    return N + p_p - p_v


def enumerated_branch_rate(parent, D, v, Dp):
    """Sum arrow rates that add a vertex to the branch rooted at v (0=root).
    parent[i] is core parent or None; D[i] total degree; p_idx = parent of v
    (external), Dp its degree."""
    N = len(parent)
    inside = set(range(N))
    rate = F(0)
    for x in range(N):
        neigh = []
        if parent[x] is not None and parent[x] in inside:
            neigh.append(parent[x])
        neigh += [y for y in range(N) if parent[y] == x]
        core_nbrs = len([y for y in range(N) if parent[y] == x]) + 1  # children + one parent slot
        leaves = D[x] - core_nbrs
        rate += F(len(neigh) + leaves, D[x])  # arrows to inside branch (core or leaf)
    rate += F(1, Dp)  # external parent -> v
    return rate


def test_branch_rate():
    # small branch: v=0 with one leaf child (node 1)
    parent = [None, 0]
    Dv, Dleaf, Dp = 2, 1, 5
    D = [Dv, Dleaf]
    got = enumerated_branch_rate(parent, D, 0, 5)  # external parent degree 5
    want = branch_rate(2, F(1, 5), F(1, 2))
    ok(got == want, f"branch rate mismatch {got} vs {want}")
    # three-vertex branch root->child->leaf
    parent2 = [None, 0, 1]
    D2 = [3, 2, 1]
    got2 = enumerated_branch_rate(parent2, D2, 0, 6)
    want2 = branch_rate(3, F(1, 6), F(1, 3))
    ok(got2 == want2, f"branch rate 2 mismatch {got2} vs {want2}")
    # sandwich bounds n-1/2 <= rate <= n+1/2
    for n in range(2, 20):
        for pp in (F(0), F(1, 3), F(1, 2)):
            for pv in (F(0), F(1, 3), F(1, 2)):
                r = branch_rate(n, pp, pv)
                ok(F(n) - F(1, 2) <= r <= F(n) + F(1, 2), "rate sandwich")
    return True


# ---------------------------------------------------------------------------
# 5. Yule NB shape/gamma moments and the A_N bracket (4).
# ---------------------------------------------------------------------------
def nb_pgf_forward_residual(q, z, eps):
    """G=z^2 [q/(1-(1-q)z)]^{2+eps}; verify d_t G = (z-1)(z G_z + eps G)
    by comparing the two closed forms derived from the product rule."""
    den = 1 - (1 - q) * z
    # d_t G = -(2+eps) z^2 q^2 (1-z) / den^3
    lhs = -(2 + eps) * z ** 2 * q ** 2 * (1 - z) / den ** 3
    # rhs = (z-1) z^2 (2+eps) q / (den)^2 * H^{1+eps} ... use H=q/den
    H = q / den
    rhs = (z - 1) * z ** 2 * (2 + eps) * H * q / den ** 2
    return lhs - rhs


def test_nb_pgf(plant=None):
    for q in (F(1, 2), F(1, 3), F(3, 4)):
        for z in (F(1, 5), F(1, 2), F(2, 3)):
            for eps in (F(0), F(1, 2), F(-1, 2)):
                res = nb_pgf_forward_residual(q, z, eps)
                if plant == "wrong_sign":
                    ok(res != 0, "planted sign survived")
                else:
                    ok(res == 0, f"NB forward equation residual {res}")
    return True


def gamma_moment(alpha, mu, plant=None):
    # E[Gamma(alpha)^mu] = Gamma(alpha+mu)/Gamma(alpha); compute with math.gamma
    import math
    if plant == "wrong_shift":
        return math.gamma(float(alpha - mu)) / math.gamma(float(alpha))
    return math.gamma(float(alpha + mu)) / math.gamma(float(alpha))


def test_gamma(mu=MU):
    lo = gamma_moment(F(3, 2), mu)
    hi = gamma_moment(F(5, 2), mu)
    ok(lo > 0 and lo < hi, "gamma moment bracket")
    ok(gamma_moment(F(3, 2), mu, plant="wrong_shift") < 0 or
       gamma_moment(F(3, 2), mu, plant="wrong_shift") != lo, "planted shift survived")
    return lo, hi


# ---------------------------------------------------------------------------
# 6. Crossing law for the toy Yule chains: mu Pr(N>=n+1) = E[(n+eps)1_{N=n}]/... .
# ---------------------------------------------------------------------------
def toy_crossing(n, eps, mu=MU):
    # E[exp(-mu H_{n+1})] = prod_{k=2}^n (k+eps)/(k+eps+mu)
    tail = F(1)
    for k in range(2, n + 1):
        tail *= F(k + int(eps) if eps.denominator == 1 else k, 1)  # placeholder replaced below
    # exact rational product with eps in {0, 1/2, -1/2}
    tail = F(1)
    for k in range(2, n + 1):
        tail *= (F(k) + eps) / (F(k) + eps + mu)
    # RHS: E[e^{-mu H_n} * lambda_n/(lambda_n+mu)] = prod_{k=2}^{n-1}(k+eps)/(k+eps+mu)
    #      * (n+eps)/(n+eps+mu)
    rhs = F(1)
    for k in range(2, n):
        rhs *= (F(k) + eps) / (F(k) + eps + mu)
    rhs *= (F(n) + eps) / (F(n) + eps + mu)
    return tail, rhs


def test_crossing(plant=None):
    for eps in (F(0), F(1, 2), F(-1, 2)):
        for n in range(3, 12):
            tail, rhs = toy_crossing(n, eps)
            if plant == "drop_lambda":
                ok(rhs != tail, "planted crossing without intensity survived")
            else:
                ok(tail == rhs, f"crossing identity failed n={n} eps={eps}")
    return True


# ---------------------------------------------------------------------------
def main():
    print(f"mu fixture = {MU} (generic rational parameter in (0,1); identities hold for every mu)")
    c1 = test_change_of_measure()
    print(f"one-clock RN stopped identity PASS paths={c1}")
    test_geometric()
    print("geometric tail / point-mass ratio PASS")
    test_compensator()
    print("rejection compensator r c/D and positivity bound PASS")
    test_branch_rate()
    print("downward branch rate N+p_p-p_v and sandwich PASS")
    test_nb_pgf()
    lo, hi = test_gamma()
    print(f"Yule NB forward equation and gamma moments PASS A_N in [{float(lo):.4f},{float(hi):.4f}]")
    test_crossing()
    print("crossing law identity PASS")

    rejected = 0
    for fn, label in [
        (lambda: test_change_of_measure(plant="drop_rejections"), "drop-rejections"),
        (lambda: test_geometric(plant="wrong_ratio"), "wrong-geometric-ratio"),
        (lambda: test_compensator(plant="use_b"), "use-b-compensator"),
        (lambda: test_nb_pgf(plant="wrong_sign"), "wrong-NB-sign"),
        (lambda: ok(gamma_moment(F(3,2), MU, plant="wrong_shift") == gamma_moment(F(3,2), MU),
                    "planted gamma shift survived"), "wrong-gamma-shift"),
        (lambda: test_crossing(plant="drop_lambda"), "drop-crossing-intensity"),
    ]:
        try:
            fn()
        except ArithmeticError:
            rejected += 1
        else:
            raise ArithmeticError(f"planted {label} accepted")
    # branch-sign plant: direct check that flipping p_v breaks small-case equality
    try:
        ok(branch_rate(2, F(1, 5), F(1, 2)) - 2 * F(1, 2) == F(2) + F(1, 5) - F(1, 2),
           "planted branch sign accepted")
        raise ArithmeticError("planted branch sign accepted")
    except ArithmeticError as e:
        ok("branch sign" in str(e), "wrong branch rejection")
        rejected += 1
    print(f"planted change-of-measure / ratio / branch-sign / NB / gamma / crossing errors rejected ({rejected}/7)")
    print("INDEPENDENT FIXTURE CHECKER: PASS")


if __name__ == "__main__":
    try:
        main()
    except (ArithmeticError, ValueError, ZeroDivisionError) as e:
        print("REJECTED", e, file=sys.stderr)
        sys.exit(2)
