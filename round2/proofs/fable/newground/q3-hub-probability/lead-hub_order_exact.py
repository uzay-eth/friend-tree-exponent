#!/usr/bin/env python3
"""Second, independent implementation (the lead's) of the certificate for board claim 3(a):
     p_1 = p_2 > p_3 > ... > p_8,   p_u = P(Z_u > 0),  Z_u = lim D_u(n)/n   in the random friend tree.

Method.  Enumerate ALL histories up to NMAX vertices with exact rational probabilities.  For a vertex v of the tree
T_m the conditional probability P(Z_v > 0 | T_m) depends only on the tree rooted at v up to isomorphism; call that
rooted shape s and the probability p(s).  Then  p_u - p_w = sum_s (q_u(s) - q_w(s)) p(s),  q_u(s) = P((T_m, u) has
shape s), and with rigorous bounds lb(s) <= p(s) <= ub(s):
     p_u - p_w >= sum_s [ (q_u - q_w)^+ lb(s) - (q_u - q_w)^- ub(s) ].
Bounds used (each re-derived by the lead; D = degree, L = number of leaf neighbours, c = D - L):
  lb, core root with L >= 1:  B(L, c) = (sum_{k<L} a_k) / (sum_{k>=0} a_k),  a_k = 1/((c+1)...(c+k)):
       the leaf count dominates the birth-death chain with rates k (up) and k/(k+c) (down), and that chain escapes to
       infinity with this probability; from many leaves the vertex is a hub with probability close to 1.
  lb, any root v with a neighbour w:  1 - ub(w)  (Theorem 3.2 of Addario-Berry et al.: Z_v + Z_w > 0 a.s. on every
       edge).  Switch off with NOEDGE=1 to see what the certificate gives without that theorem.
  ub, leaf root with parent w (k = L_w, root included; c = c_w):  P(root is ever promoted) <= F_c(k), the solution of
       (k(k+c) + k) F(k) = k(k+c) F(k+1) + (k-1) F(k-1) + 1 for k <= K with F(K+1) = 1/K  (valid stopped bound for
       every c >= 0, because 1/(k-1) is a supersolution for k >= 3); it is checked to be non-increasing in k, which
       makes it a supersolution of the true dynamics (more births and more core neighbours only help).
  ub, core root with L = 0:  it must first receive a leaf from a neighbour w:  sum over w of 1/max(1, L_w - 1), <= 1.
  ub, core root with L >= 1:  1.
Usage: hub_order_exact.py [NMAX]   (default 9; 10 takes a few minutes)
"""
import sys, os
from fractions import Fraction as Q
from functools import lru_cache
NMAX = int(sys.argv[1]) if len(sys.argv) > 1 else 9; NOEDGE = os.environ.get("NOEDGE") == "1"; U = 8; K = 40

# ---------- classifier tables (exact rationals) ----------
@lru_cache(None)
def B(l, c):
    if l == 0: return Q(0)
    a = [Q(1)]
    for k in range(1, 40): a.append(a[-1] / (c + k))
    tail = a[-1] / (c + 40) * Q(c + 41, c + 40)              # sum_{k>=40} a_k <= a_40 (1 + 1/(c+41) + ...) <= a_40 (c+41)/(c+40)
    return sum(a[:l]) / (sum(a) + tail)                      # denominator rounded up: a lower bound
@lru_cache(None)
def Ftab(c):
    # unknowns F(1..K), boundary F(K+1) = 1/K; tridiagonal system solved by forward elimination: F(k) = al_k + be_k F(k+1)
    al = {}; be = {}
    for k in range(1, K + 1):
        d = Q(k * (k + c) + k); up = Q(k * (k + c)); lo = Q(k - 1)
        if k == 1: al[k] = Q(1) / d; be[k] = up / d
        else:
            den = d - lo * be[k - 1]; al[k] = (1 + lo * al[k - 1]) / den; be[k] = up / den
    F = {K + 1: Q(1, K)}
    for k in range(K, 0, -1): F[k] = al[k] + be[k] * F[k + 1]
    assert all(F[k] >= F[k + 1] for k in range(1, K + 1)) and all(0 <= F[k] <= 1 for k in F), c
    return F
def Fup(k, c):
    return Ftab(c)[k] if k <= K else Q(1, k - 1)

# ---------- rooted shapes ----------
def canon(adj, v, parent=0):
    return tuple(sorted(canon(adj, w, v) for w in adj[v - 1] if w != parent))
def ub_shape(s):
    """upper bound for the hub probability of the root of the rooted shape s (nested sorted tuples)"""
    if len(s) == 1:                                           # leaf root; s[0] = its parent seen from the root
        w = s[0]; k = 1 + sum(1 for x in w if x == ()); c = sum(1 for x in w if x != ())
        return Fup(k, c)
    L = sum(1 for x in s if x == ())
    if L >= 1: return Q(1)
    return min(Q(1), sum(Q(1, max(1, sum(1 for y in w if y == ()) - 1)) for w in s))
def reroot(s, i):
    """the same tree rooted at the i-th child of the root"""
    rest = tuple(sorted(s[:i] + s[i + 1:])); return tuple(sorted(s[i] + (rest,)))
def lb_shape(s):
    L = sum(1 for x in s if x == ()); c = len(s) - L
    best = B(L, c) if len(s) >= 2 else Q(0)
    if not NOEDGE:
        for i in range(len(s)): best = max(best, 1 - ub_shape(reroot(s, i)))
    return best

# ---------- enumeration ----------
level = [(Q(1), ((2,), (1,)))]
for n in range(2, NMAX):
    nxt = []
    for p, adj in level:
        deg = [len(x) for x in adj]
        for w in range(1, n + 1):
            pw = sum(Q(1, deg[u - 1]) for u in adj[w - 1]) / n
            new = [list(x) for x in adj]; new[w - 1].append(n + 1); new.append([w])
            nxt.append((p * pw, tuple(tuple(x) for x in new)))
    level = nxt
assert sum(p for p, _ in level) == 1
q = {}
for p, adj in level:
    for u in range(1, min(U, NMAX) + 1):
        s = canon(adj, u); q.setdefault(s, [Q(0)] * (U + 1))[u] += p
print(f"size {NMAX}: {len(level)} histories, {len(q)} rooted shapes" + ("  (edge theorem not used)" if NOEDGE else ""))
lb = {s: lb_shape(s) for s in q}; ub = {s: ub_shape(s) for s in q}
assert all(0 <= lb[s] <= ub[s] <= 1 for s in q)
for u in range(1, min(U, NMAX) + 1):
    lo = sum(q[s][u] * lb[s] for s in q); hi = sum(q[s][u] * ub[s] for s in q)
    print(f"   p_{u} in [{float(lo):.6f}, {float(hi):.6f}]")
ok = True
for u in range(1, min(U, NMAX)):
    d = sum((q[s][u] - q[s][u + 1]) * (lb[s] if q[s][u] >= q[s][u + 1] else ub[s]) for s in q)
    e = sum((q[s][u] - q[s][u + 1]) * (ub[s] if q[s][u] >= q[s][u + 1] else lb[s]) for s in q)
    print(f"   p_{u} - p_{u+1} in [{float(d):+.6f}, {float(e):+.6f}]" + ("   (exactly 0 by symmetry)" if d == 0 == e else "   CERTIFIED > 0" if d > 0 else "   not decided"))
    ok &= (d > 0) or (u == 1 and d == 0 == e)
print("ordering p_1 = p_2 > p_3 > ... > p_%d certified: %s" % (min(U, NMAX), "YES" if ok else "NO"))
