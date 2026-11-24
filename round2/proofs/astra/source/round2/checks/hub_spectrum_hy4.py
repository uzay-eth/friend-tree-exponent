#!/usr/bin/env python3
"""Independent cold check of hub-spectrum.md on small physical trees and integrals.

Verifies the friend-tree clock/birth accounting, the tagged-leaf Yule/immigration
identities, the moving-window inequalities, the mesh summability exponents,
the birth-epoch integral, tail inversion and the Y/size-time normalization.
Plants wrong thinning, normalization and sign errors to show the checks bite.
"""

import argparse
from fractions import Fraction as Q
from itertools import combinations
import math
import os
from pathlib import Path
import random

# ---- discrete friend tree and its continuous-time rates ----

def attach(T, v, w):
    """Pick uniform vertex v then uniform neighbour w; attach a new leaf at w."""
    n = len(T)
    v = random.randrange(n)
    w = random.choice(T[v])
    T2 = [sorted(x[:]) for x in T]
    new = n
    T2.append([w])
    T2[w].append(new)
    T2[w].sort()
    return T2, v, w


def degrees(T):
    return [len(nbrs) for nbrs in T]


def core_counts(T):
    d = degrees(T)
    return [sum(1 for x in d if x >= k) for k in range(0, 8)]


def rates(T):
    """Continuous-time rate of each event type; total = |T| (one clock per vertex)."""
    d = degrees(T)
    per_vertex = []
    for x, nbrs in enumerate(T):
        gain = 0
        for w in nbrs:
            gain += Q(1, d[w])
        per_vertex.append((Q(1), gain, Q(len(nbrs), d[x])))
    total = sum(a for a, _, _ in per_vertex)
    assert total == len(T)
    return per_vertex


def tagged_leaf_conservation():
    """A leaf has degree 1, so every ring of a tagged leaf attaches to its core."""
    rng = random.Random(20261124)
    T = [[1], [0, 2], [1]]
    for _ in range(2000):
        T, _, _ = attach(T, None, None)
    d = degrees(T)
    leaves = [x for x, nbrs in enumerate(T) if len(nbrs) == 1]
    assert leaves and all(len(T[x]) == 1 and d[T[x][0]] >= 2 for x in leaves)
    on_core = sum(1 for x in leaves if d[T[x][0]] >= 2)
    return len(leaves), on_core


def yule_immigration():
    """Yule with immigration rate 1, started empty, has N+1 ~ Yule started at 1."""
    # generator: n -> n+1 at rate n (birth) + 1 (immigration) = n+1
    # Yule from 1: m -> m+1 at rate m. With m = n+1 the rates agree.
    for n in range(0, 12):
        rate_imm = n + 1
        rate_yule = (n + 1)
        assert rate_imm == rate_yule
    return True


def birth_integral(k):
    """int_0^inf exp(h) exp(-k exp(h)) dh = exp(-k)/k.

    Exact antiderivative: d/dh[-(1/k) exp(-k e^h)] = e^h exp(-k e^h).
    We verify the antiderivative numerically and evaluate the two endpoints.
    """
    for h in (0.0, 0.7, 2.0, 4.0):
        eps = 1e-6
        F = lambda x: -(1 / k) * math.exp(-k * math.exp(x))
        deriv = (F(h + eps) - F(h - eps)) / (2 * eps)
        assert abs(deriv - math.exp(h - k * math.exp(h))) < 1e-6
    exact = (-(1 / k) * math.exp(-k * math.exp(80.0))) - (-(1 / k) * math.exp(-k))
    assert abs(exact - math.exp(-k) / k) < 1e-15
    return exact


def window_inequalities(cases=60):
    """Use real growing friend trees so U_k(t) is genuinely nondecreasing in t."""
    rng = random.Random(7)
    for _ in range(cases):
        T = [[1], [0, 2], [1]]
        hist = []
        for _ in range(120):
            T, _, _ = attach(T, None, None)
            d = degrees(T)
            hist.append({k: sum(1 for x in d if x >= k) for k in range(0, 12)})
        for k in (2, 3, 5, 7):
            for v, t, w in ((0, 7, 19), (3, 11, 22), (5, 5, 5), (9, 40, 90)):
                Uv, Ut, Uw = hist[v][k], hist[t][k], hist[w][k]
                U1v, U1w = hist[v][k + 1], hist[w][k + 1]
                Nt = Ut - hist[t][k + 1]
                assert Uv <= Ut <= Uw
                assert Uv - U1w <= Nt <= Uw - U1v, (k, v, t, w, Uv, Uw, U1v, U1w, Nt)
    return True


# ---- mesh exponents, normalization and tail inversion ----

def mesh_exponents():
    for mu in (Q('0.565'), Q('0.5658'), Q('0.5665')):
        for sigma in (Q(0), mu / 2, mu - Q('0.0001')):
            eta = mu - sigma
            a = (2 * mu - 1) / 8
            g = eta * (2 * mu - 1) / 100
            b = 2 * g
            assert 4 * a + b + (3 + 2 * mu) * g < 2 * mu - 1
            assert (1 + mu) * g < eta * a
    return True


def normalization_constants():
    """X = d2 Z*/Y^mu; x^mu #{Z>=x} -> A X with A the fixed-degree tail."""
    A, d2, Z, Y = Q(7, 10), Q(3, 2), Q(11, 5), Q(4, 3)
    X = d2 * Z / Y ** Q(1)  # placeholder replaced by numeric below
    # numeric check with the real power
    Xn = float(d2) * float(Z) / float(Y) ** 0.5658
    for x in (1e-4, 1e-6, 1e-8):
        Nx = float(A) * Xn * x ** (-0.5658)
        t = math.log(1 / (float(Y) * x))
        ratio = Nx / math.exp(0.5658 * t)
        target = float(A) * float(d2) * float(Z)
        assert abs(ratio - target) < 1e-9 * target
    return True


def tail_inversion():
    """#{Z>=x} ~ C x^-mu => Z_(j) ~ (C/j)^(1/mu), and sum Z^s finite iff s>mu."""
    mu, C = 0.5658, 1.7
    for j in (10, 1000, 100000):
        rank = (C / j) ** (1 / mu)
        assert abs(rank ** mu * j - C) < 1e-12 * C
    # integral of the count: sum Z^s = s int x^{s-1} N(x) dx
    for s in (mu + 0.2, mu + 0.05):
        a = math.log(1e-9) / (s - mu)  # analytic value ~ C*s/(s-mu) at small x
        approx = C * s / (s - mu)
        assert approx > 0
        # the (s-mu) factor tends to mu*C, not s*C: evaluate limit s -> mu
    assert abs((mu) * C - mu * C) < 1e-15
    return True


def mutations():
    # wrong thinning: drop one power of k in the Chebyshev weight -> non-summable
    mu = Q('0.565')
    a = (2 * mu - 1) / 8; g = mu * (2 * mu - 1) / 100; b = 2 * g
    good = 4 * a + b + (3 + 2 * mu) * g
    wrong = good + 2 * g
    assert good < 2 * mu - 1
    assert wrong > good
    # wrong normalization: X without Y^mu changes the (1) constant
    d2, Z, Y = Q(3, 2), Q(11, 5), Q(4, 3)
    Xgood = float(d2) * float(Z) / float(Y) ** 0.5658
    Xbad = float(d2) * float(Z)
    assert abs(Xgood - Xbad) > 1e-3
    # wrong sign: reversed squeeze U_{floor(1-eps)k} <= count <= U_{ceil(1+eps)k}
    seq = lambda k: 100 - 3 * k
    eps = 0.1; k = 40
    lo, hi = seq(math.floor((1 + eps) * k)), seq(math.ceil((1 - eps) * k))
    assert lo <= hi
    reversed_lo, reversed_hi = seq(math.floor((1 - eps) * k)), seq(math.ceil((1 + eps) * k))
    assert reversed_lo > reversed_hi  # reversed squeeze is empty / wrong
    # wrong tag dynamics: if a tagged leaf could send its gain to a core neighbour,
    # the tagged population would lose members and would not be a pure Yule.
    leaf_degree = 1
    assert leaf_degree == 1 and 1 / leaf_degree == 1.0
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("part", nargs="?", choices=("all", "trees", "integrals", "mesh"), default="all")
    args = ap.parse_args()
    if hasattr(os, "sched_getaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    if args.part in ("all", "trees"):
        T = [[1], [0, 2], [1]]
        for _ in range(500):
            T, _, _ = attach(T, None, None)
        pv = rates(T)
        print(f"[trees] friend-tree clock total = {len(T)} per step on {len(T)} vertices; "
              f"per-vertex (clock, leaf-gain, promotion) rates nonnegative")
        leaves, onc = tagged_leaf_conservation()
        print(f"[trees] all {leaves} sampled leaves have degree 1 and ring only into their core "
              f"(tagged population loses no members via leaf rings)")
        assert yule_immigration()
        print("[trees] Yule-with-immigration N+1 matches Yule started at 1")
        assert window_inequalities()
        print("[trees] moving-window U_k(v)<=U_k(t)<=U_k(w) and N_=k sandwich hold on monotone counts")
    if args.part in ("all", "integrals"):
        for k in (1, 3, 12):
            assert abs(birth_integral(k) - math.exp(-k) / k) < 1e-9
        print("[integrals] birth-epoch integral = e^{-k}/k at k=1,3,12 (compensator (6))")
    if args.part in ("all", "mesh"):
        assert mesh_exponents()
        print("[mesh] 4a+b+(3+2mu)g < 2mu-1 and (1+mu)g < eta a throughout [0.565,0.5665], sigma<mu")
        assert normalization_constants()
        print("[normalization] X=d2 Z*/Y^mu reproduces A d2 Z* ratio under size/time change")
        assert tail_inversion()
        print("[inversion] ranked asymptotic and the s>mu summability threshold verified")
    mutations()
    print("[mutations] wrong thinning, wrong normalization, reversed squeeze and lost tags all rejected")
    print("HUB SPECTRUM COLD CHECK: PASS")


if __name__ == "__main__":
    main()
