> Produced by Claude Fable 5.2.

# Cross-charged drift certificates for the random friend tree

This note defines a family of computer-checkable certificates that bound the growth exponent of
E[C_n] from both sides, and proves that a certificate implies the bound. It refines the vertex-local
certificates of the earlier attempt (`attempts/x-strogatz/REPORT.md` §3) in one respect: the uncertain
"credit" terms are no longer bounded vertex by vertex in the worst case only; a part of each term is moved,
through identities that hold on every tree, to the vertex that knows it exactly ("cross-charging").
The amount moved is a free parameter (a multiplier), found by linear programming; the proof does not depend
on how it was found.

## 1. Model and notation

T_2 is one edge; given T_n, pick V uniformly, W uniformly among the neighbours of V, attach a new leaf to W.
For n >= 3 every leaf is adjacent to a core vertex (degree >= 2). For a core vertex u: L_u leaf neighbours,
c_u core neighbours, D_u = L_u + c_u. One step does one of two things (Lemma 1 of the earlier report):

* **gain at w** (w core): probability (1/n)(L_w + sum over core neighbours z of w of 1/D_z); L_w -> L_w + 1.
  We split it by the source: "own" (a leaf of w was picked as V), n*prob = L_w, and "from z" (the core
  neighbour z was picked as V and chose w), n*prob = 1/D_z.
* **conversion at w** (w core, L_w >= 1): n*prob = L_w/D_w; w goes from (L_w, c_w) to (L_w - 1, c_w + 1) and
  one of its leaves becomes a core vertex ("newborn") with (L, c) = (1, 1).

Degrees and c never decrease.

## 2. States

Fix integers Dcap >= 3, c1 >= 1, Lmax >= Dcap, cmax >= c1.

**Class.** For a core vertex w with c_w >= 1: cls(w) = (L_w, c_w) if L_w + c_w < Dcap and c_w <= c1, and
cls(w) = TOP otherwise. The class is a function of (L_w, c_w). A non-TOP class (l, k) determines the degree
l + k. A TOP vertex stays TOP for ever, and its degree is at least d_T := min(Dcap, c1 + 1).
Let A2 be the set of non-TOP classes (l, k) with k >= 2.

**State of a core vertex u.**
* detailed, if c_u <= c1: x_u = (L_u, c_u, S_u), where S_u is the multiset of the classes of the core
  neighbours of u (c_u elements);
* coarse, if c1 < c_u <= cmax: x_u = (L_u, c_u);
* beyond, if c_u > cmax: x_u = G (one state).

A coarse or beyond vertex is TOP (c_u > c1), and so is every vertex with c_u >= 1 and L_u >= Lmax
(a vertex with c_u = 0 has no core neighbour and needs no class).

**Potential.** f is a real function on states with f(L, .) = f(Lmax, .) for all L > Lmax (so f is given by a
finite table), and Phi(T) = sum over core vertices u of f(x_u). In "upper mode" f(G) is a number G; in "lower
mode" f(G) = 0.

**Multipliers.** Real numbers kappa(b) for every class b (including TOP), and nu(a, b) for a in A2 and b a
non-TOP class. Put nu(a, TOP) = 0 and nu(a, b) = 0 for a not in A2.

## 3. The drift identity

For a core vertex u write N(u) for its set of core neighbours, and

    rho_top(u) = sum of 1/D_w over w in N(u) with cls(w) = TOP,
    rho'(w \ u) = sum of 1/D_z over z in N(w), z != u          (credit that w gets from its neighbours other than u).

**Row of a detailed vertex** u = (L, c, S), with D = L + c, b = cls(u) (if c >= 1), and the conventions
f^up = f(L+1, c, S), Df := f^up - f(x_u) (zero if L >= Lmax). Define

    Row(u) =   (L + sum_{a in S, a != TOP} 1/D_a) * Df                                   [own gains, exact credit]
             + sum_{a in S, a != TOP} (1/D_a) * N(b; S - a)                              [nu received, exact part]
             + rho_top(u) * ( Df + N(b; S - TOP) - kappa(b) )                            [credit from TOP neighbours]
             + (L/D) * ( f(L-1, c+1, S + {(1,1)}) - f(x_u) + f(1, 1, {cls(L-1, c+1)}) )  [conversion, newborn]
             + sum over neighbours w in N(u) with cls(w) = a = (l, k) != TOP of
                   (l + 1/D) * ( f(x_u[a -> a+]) - f(x_u) )                              [w gains: own leaves, from u]
                 + rho'(w \ u) * ( f(x_u[a -> a+]) - f(x_u) - nu(a, b) )                 [w gains from others]
                 + (l/(l+k)) * ( f(x_u[a -> a-]) - f(x_u) )                              [w converts]
             + [b = TOP] * (1/D) * sum_{b' in S} kappa(b')                               [payments of a TOP giver]

where a+ = cls(l+1, k), a- = cls(l-1, k+1) (the term is absent if l = 0), x_u[a -> a'] replaces one copy of a
in S by a', and N(b; S') := sum over b' in S' of nu(b, b') (zero unless b is in A2).
Conventions for terms with rate zero: the conversion line is absent when L = 0; the TOP-credit line is absent
when S contains no TOP; when c = 0 (the centre of a star) there is no class b and all kappa- and nu-terms are
absent. In the conversion line f(L-1, c+1, .) means the coarse value f(L-1, c+1) when c + 1 > c1 and f(G) when
c + 1 > cmax.

**Row of a coarse or beyond vertex** u = (L, c), c > c1 (it is TOP; its neighbours' classes are not part of
the state):

    Row(u) =   L * Df + (rho(u) - rho_top(u)) * Df + rho_top(u) * (Df - kappa(TOP))
             + (L/D) * ( f(L-1, c+1) - f(x_u) + f(1, 1, {TOP}) )
             + (1/D) * sum over v in N(u) of kappa(cls v),

with rho(u) the total credit of u, and for a beyond vertex f(x_u) = f(L-1, c+1) = f(G), Df = 0. The conversion
line is absent when L = 0, and f(L-1, c+1) means f(G) when c + 1 > cmax.

**Lemma 1 (drift identity).** For every tree T_n, n >= 3, and all values of f, kappa, nu:

    n * E[Phi(T_{n+1}) - Phi(T_n) | T_n]  =  sum over core vertices u of Row(u).

*Proof.* First take kappa = nu = 0. A gain at w (rate L_w + rho(w)) changes x_w by L_w -> L_w + 1, and changes
x_u for the core neighbours u of w whose state records cls(w), when cls(w) changes: these are the detailed
neighbours u, and only if cls(w) != TOP (TOP is absorbing); the class changes from a = (l, k) to a+. A
conversion at w (rate L_w/D_w) changes x_w to (L_w - 1, c_w + 1, S_w + {(1,1)}) (the newborn has class (1,1)),
creates the newborn in state (1, 1, {cls(L_w - 1, c_w + 1)}), and changes cls(w) from (l, k) to a- as seen by
the detailed neighbours of w, if cls(w) != TOP. Nothing else changes: all other pairs (L, c) are untouched,
hence all other classes and states. Summing rate times change over all events and regrouping by the vertex whose
f-value changes gives the sum of the rows with kappa = nu = 0; in the row of u the rate of "w gains" is
l + 1/D_u + rho'(w \ u) (own leaves of w, u itself, the other neighbours of w), and rho(u) is split into the
exact part (non-TOP neighbours, whose degree is known from their class) and rho_top(u).

The kappa terms add up to zero: sum_u kappa(cls u) rho_top(u) = sum over ordered pairs (u, w) of adjacent core
vertices with cls(w) = TOP of kappa(cls u)/D_w = sum over TOP vertices w of (1/D_w) sum_{v in N(w)}
kappa(cls v); the left side is subtracted in the rows of the receivers u (third line of the detailed row, with
b = cls u; coarse and beyond rows with cls u = TOP), the right side is added in the rows of the TOP vertices w
(last line of each row).

The nu terms add up to zero: fix a core vertex w with cls(w) = a in A2 (then w is detailed and non-TOP). The
neighbours u of w subtract sum_{u in N(w)} nu(a, cls u) rho'(w \ u) in total (in a coarse or beyond neighbour u
nothing is subtracted, consistently with nu(a, TOP) = 0; a detailed TOP neighbour subtracts nu(a, TOP) = 0).
This equals sum_{z in N(w)} (1/D_z) sum_{u in N(w), u != z} nu(a, cls u) = sum_{z in N(w)} (1/D_z) N(a; S_w - cls z),
which is what the row of w adds (second line for non-TOP z, the N-part of the third line for TOP z). ∎

The identity, including both cancellations, is checked in exact rational arithmetic on random trees with
random values of f, kappa, nu by `xcharge.py test` (the script enumerates every move (V, W) of the tree).

## 4. Intervals

In the row of a detailed vertex u everything is a function of x_u except
* rho_top(u), which lies in [0, n_top / d_T] (n_top = number of TOP neighbours; each has degree >= d_T), and
* for each non-TOP class a = (l, k) with k >= 2, the sum of rho'(w \ u) over the m_a neighbours w of class a,
  which lies in [0, m_a (k-1)/2] (w has k - 1 other core neighbours, each of degree >= 2).
In the row of a coarse or beyond vertex with c core neighbours, n_top of them TOP: rho_top in [0, n_top/d_T],
rho - rho_top in [0, (c - n_top)/2], and sum_{v} kappa(cls v) lies between c * min_b kappa(b) and
c * max_b kappa(b).

## 5. Certificates

**Upper certificate at lambda.** Numbers f >= 1 on all states, G >= 1, multipliers with kappa(TOP) >= 0, and
Kp := max(0, max_b kappa(b)), such that for every state x in the box (L <= Lmax, c <= cmax):

    main(x) + [coarse] (c/D) Kp + sum_j max(0, hi_j * e_j(x)) [+ coarse: max(0, (c/2) Df, (c/d_T)(Df - kappa(TOP)))]
        <=  lambda * f(x),                                                                             (U)

where main(x) collects the exact terms of the row, and (hi_j, e_j) are the interval bounds and coefficients of
the uncertain terms of §4 [hi = n_top/d_T with e = Df + N - kappa(b); hi = m_a (k-1)/2 with
e = f(x[a -> a+]) - f(x) - nu(a, b)]; moreover

* **(tail)** for every (c, S), respectively every coarse c, the same inequality for all L > Lmax; for these
  states Df = 0, f does not depend on L, the newborn sees a TOP creator, and the left side is an affine function
  of y = 1/D on (0, 1/(Lmax + 1 + c)] (the rates are L/D = 1 - c y, 1/D = y, and constants; the coefficients
  hi_j and e_j do not depend on y). So it is enough to check y = 0 and y = 1/(Lmax + 1 + c).
* **(beyond)** f(1, 1, {TOP}) <= lambda G and Kp <= lambda G.

**Lemma 2.** An upper certificate at lambda implies n E[Phi_{n+1} - Phi_n | T_n] <= lambda Phi_n for all trees,
hence E[C_n] <= E[Phi_n] <= f(x_0) prod_{m=3}^{n-1} (1 + lambda/m) for n >= 3, where x_0 = (2, 0, {}) is the
state of the centre of T_3.

*Proof.* By Lemma 1 it suffices that Row(u) <= lambda f(x_u) for every core vertex of every tree. For a box
vertex this is (U), because each uncertain term is (a number in [0, hi_j]) times e_j, which is at most
max(0, hi_j e_j); for a coarse vertex (rho - rho_top) Df + rho_top (Df - kappa(TOP)) <= (c - n_top)(1/2) Df^+
+ n_top (1/d_T)(Df - kappa(TOP))^+ <= max(0, (c/2) Df, (c/d_T)(Df - kappa(TOP))), and the payments are at most
(c/D) Kp. For L > Lmax use (tail). For a beyond vertex Row(u) = (L/D) f(1,1,{TOP}) - kappa(TOP) rho_top(u)
+ (1/D) sum_v kappa(cls v) <= (L/D) f(1,1,{TOP}) + (c/D) Kp <= max(f(1,1,{TOP}), Kp) <= lambda G, using
kappa(TOP) >= 0 and L/D + c/D = 1. Finally Phi_n >= C_n because f >= 1, and Phi_3 = f(x_0). ∎

**Lower certificate at delta, 0 < delta <= 1.** Numbers f >= 0 on the box states, not all zero, f(G) = 0,
multipliers with kappa(b) >= 0 for all b and kappa(TOP) = 0, such that for every box state and both endpoints
of every tail family

    main(x) + sum_j min(0, hi_j * e_j(x)) [+ coarse: min(0, (c/2) Df)]  >=  delta * f(x).                (L)

(The payment term of coarse vertices is omitted: it is >= 0.)

**Lemma 3.** A lower certificate at delta implies n E[Phi_{n+1} - Phi_n | T_n] >= delta Phi_n for all trees,
hence E[C_n] >= E[Phi_n]/sup f >= (f(x_0)/sup f) prod_{m=3}^{n-1} (1 + delta/m).

*Proof.* As for Lemma 2 with the inequalities reversed. For a coarse vertex the uncertain credit term is
rho Df >= min(0, (c/2) Df) (kappa(TOP) = 0) and the payments are >= 0. For a beyond vertex
Row(u) = (L/D) f(1,1,{TOP}) + (1/D) sum_v kappa(cls v) >= 0 = delta f(G). ∎

## 6. Checking a certificate

`xverify.py CERT.pkl VALUE` converts the stored floating-point values to exact rationals, clips them to the
sign constraints above, and checks every inequality (U) or (L), the tail endpoints and the beyond conditions
in exact rational arithmetic. The certificate tables are produced by `xcharge.py` (a linear program in f,
kappa, nu and auxiliary variables for the max/min terms, solved by HiGHS at a value slightly on the safe side
of the certified one); how they are produced is irrelevant for the proof.

## 7. Almost-sure version

A lower certificate satisfies the hypotheses (H1)-(H5) of `as_lower.md` (Theorem B) provided every newborn
state (1, 1, {b}) has f > 0. (H3) with b0 = Dcap: a gain at w changes x_w and, if w is non-TOP, the states of
its detailed core neighbours (a non-TOP vertex has at most Dcap - 1 core neighbours); a conversion at w changes
x_w, creates one vertex, and changes the states of the core neighbours of w if w was non-TOP. (H5): for a vertex
v with f(x_v) > 0 (so c_v <= cmax) the events that can change its f-value are its own gain when L_v < Lmax
(rate <= Lmax + cmax/2), its own conversion (rate <= 1), and, for each of its at most cmax core neighbours w,
a gain at a non-TOP w (rate L_w + rho_w <= D_w < Dcap) or a conversion at w (rate <= 1); so
R = Lmax + cmax/2 + 1 + cmax (Dcap + 1) works. `xverify.py` prints the smallest newborn value.
