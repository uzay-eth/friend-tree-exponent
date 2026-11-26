#!/usr/bin/env python3
"""Exact finite checks for Theorem 15 of the paper (Appendix E; standard library only).

Model: the original friend tree from one edge. At size n pick a vertex v uniformly,
then a uniform neighbour u of v, and attach a new leaf to u.

Part A  exhaustive exact law of T_n over unlabelled trees, 2 <= n <= NMAX (default 15):
        total probability 1; a_n = E[C_n]; prime support of den(a_n) (fact F2);
        the near-star recurrences of Lemma E.3 against the exhaustive law for 7 <= n <= NMAX;
        the window law of Lemma E.2 (v_p(a_n) = -(n-p-2) and its leading residue)
        wherever p+4 <= n <= min(2p-1, NMAX).
Part B  the near-star recurrences run forward from the exhaustive size-7 values to size p+3
        for every prime 17 <= p <= PMAX (default 1200): condition (C1), the residues of
        Proposition E.4, and Lambda_p = -4 for p = 2 (mod 3).
Part C  the closed-form identities used in the proof of Proposition E.4 (exact, at enough points
        to prove each rational-function or polynomial identity).
This is evidence for the finite steps; the proofs are in Appendix E of the paper.
Usage: python3 -I check_window.py [NMAX] [PMAX]
"""
import sys
from fractions import Fraction as F

NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 15
PMAX = int(sys.argv[2]) if len(sys.argv) > 2 else 1200
assert NMAX >= 9


def canon(adj):
    """Canonical string of an unlabelled free tree (AHU encoding rooted at the centre(s))."""
    n = len(adj)
    deg = [len(a) for a in adj]
    layer = [v for v in range(n) if deg[v] <= 1]
    left, d = n, deg[:]
    while left > 2:
        left -= len(layer)
        nxt = []
        for v in layer:
            for u in adj[v]:
                d[u] -= 1
                if d[u] == 1:
                    nxt.append(u)
        layer = nxt
    def enc(v, parent):
        return "(" + "".join(sorted(enc(u, v) for u in adj[v] if u != parent)) + ")"
    return min(enc(c, -1) for c in layer)


def from_canon(s):
    adj, stack = [], []
    for ch in s:
        if ch == "(":
            adj.append([])
            v = len(adj) - 1
            if stack:
                adj[stack[-1]].append(v)
                adj[v].append(stack[-1])
            stack.append(v)
        else:
            stack.pop()
    return adj


def cores(adj):
    return sum(1 for a in adj if len(a) >= 2)


def shape(adj):
    """Near-star label for trees with K = n-1-maxdeg <= 2 (n >= 7), else None."""
    n = len(adj)
    deg = [len(a) for a in adj]
    md = max(deg)
    K = n - 1 - md
    if K > 2:
        return None
    if K == 0:
        return "S"
    if K == 1:
        return "B"
    h = deg.index(md)
    nb = [u for u in adj[h] if deg[u] >= 2]
    if len(nb) == 2:
        return "b"
    return "a" if deg[nb[0]] == 3 else "c"


def vp(x, p):
    if x == 0:
        return None
    k, a, b = 0, x.numerator, x.denominator
    while a % p == 0:
        a //= p; k += 1
    while b % p == 0:
        b //= p; k -= 1
    return k


def res(x, p):
    assert x.denominator % p != 0
    return x.numerator * pow(x.denominator, -1, p) % p


def is_prime(q):
    return q >= 2 and all(q % r for r in range(2, int(q ** 0.5) + 1))


def step(n, v):
    """Lemma E.3 recurrences: (S,B,a,b,c) at size n >= 7 -> size n+1."""
    S, B, a, b, c = v
    return [S * F(n - 1, n),
            B * (n - F(5, 2)) / n + S / n,
            a * (n - F(11, 3)) / n + B * (F(1, n - 2) + 1) / n,
            b * F(n - 4, n) + B * F(n - 3, n - 2) / n,
            c * (n - F(7, 2)) / n + B / (2 * n)]


fails = 0
def check(ok, msg):
    global fails
    if not ok:
        fails += 1
        print("FAIL", msg)

# ---------------- Part A: exhaustive law -----------------
law = {canon([[1], [0]]): F(1)}
laws = {2: law}
a = {2: F(0)}
near = {}
for n in range(2, NMAX):
    new = {}
    for s, pr in law.items():
        adj = from_canon(s)
        deg = [len(x) for x in adj]
        for u in range(n):
            w = sum(F(1, deg[v]) for v in adj[u]) / n
            nadj = [list(x) for x in adj] + [[u]]
            nadj[u].append(n)
            t = canon(nadj)
            new[t] = new.get(t, F(0)) + pr * w
    law = new
    m = n + 1
    laws[m] = law
    check(sum(law.values()) == 1, f"total probability at n={m}")
    a[m] = sum(pr * cores(from_canon(s)) for s, pr in law.items())
    if m >= 7:
        tab = dict.fromkeys("SBabc", F(0))
        for s, pr in law.items():
            k = shape(from_canon(s))
            if k:
                tab[k] += pr
        near[m] = [tab[k] for k in "SBabc"]
print("Part A: exhaustive law, n <= %d (classes at NMAX: %d)" % (NMAX, len(law)))
print("  a_3 = %s, a_4 = %s, a_%d = %s" % (a[3], a[4], NMAX, a[NMAX]))
# F2: strip all primes < m from den(a_m); the rest must be 1
for m in range(3, NMAX + 1):
    den = a[m].denominator
    for q in range(2, m):
        if is_prime(q):
            while den % q == 0:
                den //= q
    check(den == 1, f"F2: den(a_{m}) has a prime >= {m}")
print("  F2 prime support: every prime dividing den(a_n) is < n, 3 <= n <= %d" % NMAX)
v = near[7]
print("  size-7 near-star law S,B,a,b,c =", ", ".join(map(str, v)))
check(v == [F(1, 3), F(71, 288), F(781, 5184), F(541, 8640), F(23, 360)], "size-7 values")
for n in range(7, NMAX):
    v = step(n, v)
    check(v == near[n + 1], f"Lemma E.3 recurrences vs exhaustive law at n={n + 1}")
print("  Lemma E.3 recurrences reproduce the exhaustive law of S,B,a,b,c for 7 <= n <= %d" % NMAX)
for p in range(5, NMAX - 3):
    if not is_prime(p):
        continue
    # Lemma E.2 directly from the exhaustive law: (C1) on every tree with a vertex of degree >= p
    # at size p+3, Lambda_p from its definition, then the predicted valuation and residue.
    lam = 0
    for s, pr in laws[p + 3].items():
        adj = from_canon(s)
        deg = [len(x) for x in adj]
        if max(deg) >= p:
            check(vp(pr, p) >= -1, f"(C1) exhaustive p={p}")
        if p in deg:
            h = deg.index(p)
            k = sum(1 for u in adj[h] if deg[u] >= 2)
            lam += res(p * pr, p) * k
    lam %= p
    for n in range(p + 4, min(2 * p - 1, NMAX) + 1):
        pred = lam * (-1) ** (n - p - 3)
        for j in range(p + 3, n):
            pred = pred * pow(j, -1, p)
        pred %= p
        x = a[n] * F(p) ** (n - p - 2)
        if lam:
            check(vp(a[n], p) == -(n - p - 2) and res(x, p) == pred, f"window p={p} n={n}")
        else:
            check(vp(a[n], p) > -(n - p - 2), f"window p={p} n={n} (Lambda_p = 0)")
        print("  window p=%d n=%d: Lambda_p = %d mod p, v_p(a_n) = %d, residue of p^(n-p-2) a_n %s" % (
            p, n, lam, vp(a[n], p), "= predicted %d" % pred if lam else "n/a"))

# ---------------- Part B: forward pass to p+3 -----------------
primes = [q for q in range(17, PMAX + 1) if is_prime(q)]
want = set()
for q in primes:
    want |= {q, q + 2, q + 3}
snap, v, n = {}, near[7], 7
while n < PMAX + 3:
    v = step(n, v); n += 1
    if n in want:
        snap[n] = v
cnt2 = 0
for p in primes:
    S, B, xa, xb, xc = snap[p + 3]
    check(all(vp(x, p) >= -1 for x in (S, B, xa, xb, xc)), f"(C1) p={p}")
    check(vp(S, p) == 0 and vp(B, p) >= 0, f"S,B integral at p+3, p={p}")
    _, Bp, ap, bp, cp = snap[p]
    check(res(Bp, p) == res(F(-4, 3), p), f"B_p p={p}")
    check(res(cp, p) == res(F(-4, 15), p), f"x_c(p) p={p}")
    check(res(bp, p) == res(F(-5, 9), p), f"x_b(p) p={p}")
    check(res(snap[p + 2][1], p) == res(F(-4, 3), p), f"B_(p+2) p={p}")
    check(res(p * xc, p) == 0 and res(p * xb, p) == res(F(-4, 3), p), f"res x_b, x_c p={p}")
    lam = (res(p * xa, p) + 2 * res(p * xb, p) + res(p * xc, p)) % p
    if p % 3 == 2:
        cnt2 += 1
        check(res(ap, p) == res(F(-1, 10), p), f"x_a(p) p={p}")
        check(res(p * xa, p) == res(F(-4, 3), p) and lam == res(F(-4), p), f"Lambda_p p={p}")
print("Part B: %d primes 17 <= p <= %d: (C1), B_p, B_(p+2), x_b(p), x_c(p), res x_b, res x_c checked;" % (len(primes), PMAX))
print("  %d primes p = 2 (mod 3): x_a(p) = -1/10, res(p x_a(p+3)) = -4/3, Lambda_p = -4 checked" % cnt2)

# ---------------- Part C: identities -----------------
K = F(71, 288) - F(2, 9)
check(K == F(7, 288), "K")
pts = range(7, 40)
for n in pts:
    Bp = F(4, 3 * (n - 1)); Bq = F(4, 3 * n)
    check(n * Bq == Bp * (n - F(5, 2)) + F(2, n - 1), "broom particular solution")
    Ra = lambda m: F(4, 5 * (m - 2)) - F(3, 10 * (m - 1))
    check(n * Ra(n + 1) == F(4, 3 * (n - 2)) + Ra(n) * (n - F(11, 3)), "telescoper R_a")
    Rc = lambda m: F(4, 15 * (m - 1))
    check(n * Rc(n + 1) == F(2, 3 * (n - 1)) + Rc(n) * (n - F(7, 2)), "telescoper R_c")
G = lambda j: 24 - F(148, 3) * j + F(100, 3) * j**2 - F(26, 3) * j**3 + F(2, 3) * j**4
for j in range(-10, 10):   # degree-5 polynomial identity: 20 points prove it
    check((j - F(5, 2)) * G(j + 1) - j * G(j) == j * (j - 1) * (j - 3) ** 2, "Gosper certificate")
check(G(7) == -60, "G(7)")
check(15 * (F(541, 8640) - F(14, 135) - K * G(7) / 360) == F(-5, 9), "x_b(p) constant")
check(F(4, 3) + K * F(-768, 7) == F(-4, 3), "beta constant")
for p in [q for q in range(17, 400) if is_prime(q)]:
    h = F(1)
    for m in range(7, p):
        h *= (m - F(5, 2)) / m
    check(vp(h, p) == 1 and res(h / p, p) == res(F(-1024, 35), p), f"Wilson step h_B(p)/p, p={p}")
print("Part C: broom particular solution, R_a, R_c, Gosper certificate, constants, Wilson step (17 <= p < 400) checked")
print("RESULT:", "PASS" if fails == 0 else "FAIL (%d)" % fails)
sys.exit(1 if fails else 0)
