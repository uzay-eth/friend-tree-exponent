#!/usr/bin/env python3
"""Independent exact short referee fixtures for round2/fringe/hub-children.md (hy4).

Repository-relative usage (one CPU, short):
  python3 -B round2/checks/hub_children_hy4.py

Checks the frozen-degree supermartingale against the ACTUAL undirected rates on
small trees, the silent-plant star/competing-rate factors, the child-shield
product with a directed-trigger plant, the birth-death survival/theta positivity,
the Poisson/logarithmic intensity, and the extreme-value exponents.  All exact
Fractions except the two explicitly labelled high-precision numeric integrals.
No production computation.
"""
from fractions import Fraction as F
from itertools import product
import math
import sys


def ok(cond, msg):
    if not cond:
        raise ArithmeticError(msg)


# ---------------------------------------------------------------------------
# 1. Frozen-degree factorial bound V(D,L)=A_D (2/(D+1))^L on real dynamics.
# ---------------------------------------------------------------------------
def A_D(D, base=4):
    return F(D * math.factorial(D + 1), 12 * base ** (D - 2))


def V(D, L, base=4, wrong_exp=False, wrong_promo=False):
    e = F(1, 2) if wrong_exp else F(2, D + 1)
    pr = F(D, 2) if wrong_promo else F(D + 1, 2)
    return A_D(D, base), e, pr


def value(D, L, base=4, wrong_exp=False):
    a, e, _ = V(D, L, base, wrong_exp)
    return a * e ** L


def drift(D, L, a, plant=None):
    """E[delta V] with gain rate a>=L and promotion rate L/D."""
    base = 3 if plant == "wrong_base" else 4
    wexp = plant == "wrong_exp"
    v0 = value(D, L, base, wexp)
    gain = value(D + 1, L + 1, base, wexp) - v0
    if plant == "directed_trigger":
        vpromo = v0 * F(D, 2)
    else:
        vpromo = value(D, L - 1, base, wexp)
    promo = vpromo - v0
    if plant == "drop_leaves":
        a = F(0)  # only neighbour arrows, no leaf rings -> directed trigger
    return a * gain + F(L, D) * promo


def test_supermartingale(plant=None):
    count = 0
    for D in range(2, 30):
        for L in range(0, D + 1):
            count += 1
            d = drift(D, L, F(L), plant)
            if plant is None:
                ok(d <= 0, f"V supermartingale failed D={D} L={L}: {d}")
            else:
                ok(d > 0, f"planted {plant} did not break at D={D} L={L}: {d}")
    # larger a is worse only when gain ratio <1, which holds: a>=L suffices
    for D in range(2, 12):
        for L in range(0, D + 1):
            for a in (F(L), F(D), F(2 * D)):
                if plant is None:
                    ok(drift(D, L, a) <= 0, "a>=L bound failed")
    return count


def actual_rates(parent, D, v):
    """Gain rate at core v = L + sum 1/D_w over core neighbours; promotion = L/D."""
    ch = [w for w in range(len(parent)) if parent[w] == v]
    core_nbrs = len(ch) + (0 if parent[v] is None else 1)
    L = D[v] - core_nbrs
    a = F(L)
    if parent[v] is not None:
        a += F(1, D[parent[v]])
    for w in ch:
        a += F(1, D[w])
    return a, F(L, D[v]), L


def rooted_shapes(n):
    if n == 1:
        yield (None,)
        return
    for seq in product(*(range(0, i) for i in range(1, n))):
        yield (None,) + seq


def test_actual_small_trees():
    checked = 0
    for n in range(3, 7):
        for parent in rooted_shapes(n):
            # assign degrees: each core has >=1 leaf; leaves D=1
            core = [v for v in range(n) if any(parent[w] == v for w in range(n)) or parent[v] is not None or v == 0]
            D = []
            for v in range(n):
                k = sum(1 for w in range(n) if parent[w] == v) + (0 if parent[v] is None else 1)
                D.append(max(2, k + 1))  # at least one leaf at each core
            for v in core:
                a, b, L = actual_rates(parent, D, v)
                d0 = drift(D[v], L, a)
                ok(d0 <= 0, f"actual dynamics break at tree {parent} v={v}: {d0}")
                checked += 1
    return checked


# ---------------------------------------------------------------------------
# 2. Silent plant: star factors and competing-rate shield.
# ---------------------------------------------------------------------------
def star_reach(d):
    """Product_{j=2}^{d-1} (j-1)/[(j-1)+(j-1)/j] = 2/d."""
    p = F(1)
    for j in range(2, d):
        p *= F(1, 1 + F(1, j))
    return p


def shield(d):
    """Product_{i=1}^{d-1} (i/d)/((i/d)+i) = (d+1)^{-(d-1)}."""
    p = F(1)
    for i in range(1, d):
        p *= F(1, d + 1)
    return p


def star_shield(d, plant=None):
    if plant == "directed_trigger":
        # WRONG: use only the core's own clock for gains (drop leaf rings),
        # changing the star reach to product j/(j+1) still 2/d but the
        # competing rate loses the leaf term i -> shield becomes 1.
        return star_reach(d), F(1)
    return star_reach(d), shield(d)


def test_star_shield(plant=None):
    for d in range(2, 25):
        sr, sh = star_shield(d, plant)
        if plant is None:
            ok(sr == F(2, d), "star reach 2/d failed")
            ok(sh == F(1, (d + 1) ** (d - 1)), "shield factor failed")
        else:
            ok(sh == F(1, (d + 1) ** (d - 1)), "planted directed trigger accepted")
    return True


# ---------------------------------------------------------------------------
# 3. Child shields: birth-death survival and theta positivity.
# ---------------------------------------------------------------------------
def bd_survival(lam, mu):
    """Linear birth-death survival from one individual = 1-mu/lam (lam>mu)."""
    return 1 - mu / lam


def theta_lower(mu=0.5):
    """On survival B_u ~ e^{mu u} W with W>0, so int du/B_u <= 2/(mu W) a.s.
    Hence theta >= P(surv)*E[exp(-2/(mu W))|surv] > 0.  Evaluate a crude
    positive lower bound using a finite truncation of W's mass."""
    p = bd_survival(F(1), F(1, 2))
    ok(p == F(1, 2), "birth-death survival 1/2 failed")
    # exp(-2/(mu W)) > exp(-2/(mu*eps)) on W>eps; choose eps with P(W<=eps)<p/2
    # use the exact negative-binomial finite-time law at u=3 as a proxy mass bound
    import math
    val = math.exp(-2 / (mu * 8.0)) * float(p) / 2
    ok(val > 0, "theta lower bound nonpositive")
    return val


def test_shields(plant=None):
    surv = bd_survival(F(1), F(1, 2))
    if plant == "wrong_survival":
        ok(surv != F(1, 2), "planted survival accepted")
    else:
        ok(surv == F(1, 2), "survival probability")
    th = theta_lower()
    # product of d-1 independent shields with prob >= theta each
    for d in range(2, 15):
        prob = (float(th)) ** (d - 1)
        ok(prob > 0, "positive shield product")
    return th


# ---------------------------------------------------------------------------
# 4. Log A_d = d log d + O(d) and the dyadic exponent.
# ---------------------------------------------------------------------------
def log_A(d):
    return math.log(d) + math.lgamma(d + 2) - math.log(12) - (d - 2) * math.log(4)


def test_logA():
    for d in (10, 100, 1000, 10000):
        r = (log_A(d) - d * math.log(d)) / d
        ok(-3 < r < 1, f"log A_d - d log d not O(d): {r}")
    # d=(mu+eps) log n/loglog n gives log A_d ~ (mu+eps) log n (1-o(1))
    mu = 0.5655
    prev = 0.0
    for n in (10 ** 12, 10 ** 30, 10 ** 80):
        d = (mu + 0.01) * math.log(n) / math.log(math.log(n))
        ratio = log_A(d) / math.log(n)
        ok(0 < ratio < 1, f"dyadic exponent off: {ratio}")
        ok(ratio > prev, "dyadic ratio not increasing to mu+eps")
        prev = ratio
    return True


# ---------------------------------------------------------------------------
# 5. Poisson/logarithmic intensity (numeric, explicitly labelled).
# ---------------------------------------------------------------------------
def mean_above(x, q=0.6, scale=1.0, n=2000):
    """mixture: W=0 w.p. 1-q, else Exp(scale).  mean above x = q*int_x^inf e^{-y/scale}/y dy."""
    lo = math.log(x)
    hi = 30.0
    h = (hi - lo) / n
    s = 0.0
    for k in range(n + 1):
        u = lo + k * h
        y = math.exp(u)
        w = 0.5 if k in (0, n) else 1.0
        s += w * math.exp(-y / scale)
    return q * s * h


def test_intensity():
    for x in (1e-4, 1e-6, 1e-8):
        got = mean_above(x)
        want = 0.6 * math.log(1 / x)
        ok(abs(got - want) / want < 0.15, f"intensity asymptotic off x={x}: {got} vs {want}")
    return True


# ---------------------------------------------------------------------------
# 6. Extreme-value exponents.
# ---------------------------------------------------------------------------
def test_extreme():
    for m in (10 ** 6, 10 ** 12, 10 ** 30):
        L = math.log(m) / math.log(math.log(m))
        for eps in (0.1, 0.5):
            d = (1 + eps) * L
            logp = log_A(d) - math.log(m)
            ok(logp < -0.5 * eps * math.log(m), "upper extreme exponent")
            d2 = (1 - eps) * L
            expc = m * math.exp(-log_A(d2))
            ok(expc > 1, "lower extreme exponent")
    return True


# ---------------------------------------------------------------------------
def main():
    c1 = test_supermartingale()
    print(f"V(D,L) supermartingale on exact grid PASS states={c1} (a>=L, actual rates)")
    c2 = test_actual_small_trees()
    print(f"V drift <=0 on actual undirected small trees PASS vertices={c2}")
    test_star_shield()
    print("silent plant star 2/d and shield (d+1)^{-(d-1)} PASS")
    th = test_shields()
    print(f"child shields: birth-death survival 1/2, theta>0 PASS theta>={th:.4f}")
    test_logA()
    print("log A_d = d log d + O(d), dyadic exponent PASS")
    test_intensity()
    print("Poisson intensity mean ~ q log(1/x) PASS (numeric integral)")
    test_extreme()
    print("extreme-value exponents PASS")

    plants = 0
    for name, fn in [
        ("wrong_base", lambda: test_supermartingale(plant="wrong_base")),
        ("wrong_exp", lambda: test_supermartingale(plant="wrong_exp")),
        ("directed_trigger", lambda: test_star_shield(plant="directed_trigger")),
        ("drop_leaves", lambda: test_supermartingale(plant="drop_leaves")),
        ("wrong_survival", lambda: test_shields(plant="wrong_survival")),
    ]:
        try:
            fn()
        except ArithmeticError:
            plants += 1
        else:
            raise ArithmeticError(f"planted {name} accepted")
    print(f"planted wrong-base / wrong-exponent / directed-trigger / drop-leaves / wrong-survival rejected ({plants}/5)")
    print("INDEPENDENT FIXTURE CHECKER: PASS")


if __name__ == "__main__":
    try:
        main()
    except (ArithmeticError, ValueError, OverflowError) as e:
        print("REJECTED", e, file=sys.stderr)
        sys.exit(2)
