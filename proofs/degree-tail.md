> Produced by Claude Fable 5.2.

# The degree tail of the nucleus has the same exponent

Model: the random friend tree. T_2 is one edge; given T_n, pick V_n uniformly among the n vertices, pick W_n
uniformly among the neighbours of V_n, attach a new leaf to W_n. D_v(n) is the degree of v in T_n,
C_n = #{v : D_v(n) >= 2} (core vertices), and for k >= 2

    N_k(n) := #{v : D_v(n) >= k},        so N_2(n) = C_n.

Krapivsky and Redner (J. Stat. Mech. 2017, eq. (1)) report numerically that the number of vertices of degree k
behaves like N^mu k^(-(1+mu)). The theorem below proves a cumulative, logarithmic form of this: the expected
number of vertices of degree at least k in the n-vertex tree is at least a constant times E[C_{n/(8k)}] and at
most an explicit expression in the E[C_m] with m >= 2n/k, uniformly in k; on the logarithmic scale it grows
like E[C_{n/k}] (part (d)). No two-sided comparison with E[C_{n/k}] up to constants is claimed.

**Theorem D.** Write a_m := E[C_m].

* **(a)** For all n >= 3, k >= 2 and every integer m0 with 3 <= m0 <= n:

      E[N_k(n)]  <=  a_{m0} + (2n/k) * sum_{m=m0}^{n-1} (a_{m+1} - a_m)/(m+1).

* **(b)** For all k >= 2, m' >= 3 and n >= 8 m' k:

      E[N_k(n)]  >=  a_{m'} / 44160.

* **(c)** Suppose a m^delta <= a_m <= A m^lambda for all m >= 3, with constants a, A > 0 and
  0 < delta <= lambda < 1. Then for all k >= 2 and n >= 24 k:

      (a/44160) (n/(12k))^delta  <=  E[N_k(n)]  <=  A (2 + 1/(1-lambda)) (3n/k)^lambda.

* **(d)** If log a_m / log m -> mu in (0, 1) as m -> infinity, then log E[N_k(n)] / log(n/k) -> mu as
  n/k -> infinity, uniformly in k >= 2.

**Corollary (unconditional, from Claim 2 of `RESULT.md`).** By Gautschi's inequality
n (n+1)^(s-1) < Gamma(n+s)/Gamma(n) < n^s for 0 < s < 1, Claim 2 gives 0.44 m^0.55 <= a_m <= 3.974 m^0.58 for all
m >= 3 (0.5022 (3/4)^0.45 = 0.4412). Hence for all k >= 2 and n >= 24 k

    (0.44/44160) (n/(12k))^0.55  <=  E[N_k(n)]  <=  33 (n/k)^0.58        (3.974 * (2 + 1/0.42) * 3^0.58 = 32.93).

With the sharper certified exponents of the other leads and the existence of mu (their claims, see
RESULT.md; they are inputs here, not proved in this note), (d) says that the cumulative degree tail has the
same logarithmic exponent mu, between 0.5650 and 0.5665. This is the cumulative form of Krapivsky and Redner's
statement. A cumulative logarithmic exponent does not by itself give a pointwise density exponent, so
1 + mu = 1.566 (their nu) is what a power law would give, not something proved here.

## Proofs

Throughout, F_j is the sigma-field of T_3, ..., T_j, and for a vertex w of T_j
P(W_j = w | F_j) = (1/j) sum over neighbours u of w of 1/D_u(j), because V_j is uniform on the j vertices and,
given V_j = u, W_j is uniform on the D_u(j) neighbours of u.

**Lemma 1 (degrees grow at most linearly in mean).** Let m >= 3 and v a vertex of T_m. For all n >= m

    E[D_v(n) | F_m] <= D_v(m) * n/m,      E[D_v(n)(D_v(n)+1) | F_m] <= D_v(m)(D_v(m)+1) * n(n+1)/(m(m+1)).

*Proof.* D_v increases by one exactly when W_j = v, and P(W_j = v | F_j) <= D_v(j)/j since each of the D_v(j)
neighbours contributes at most 1/j. So E[D_v(j+1) | F_j] <= D_v(j)(j+1)/j, and iterating gives the first bound.
With Q = D(D+1): Q increases by 2(D_v(j)+1) when W_j = v, so E[Q(j+1) | F_j] <= Q(j) + 2(D_v(j)+1) D_v(j)/j
= Q(j)(j+2)/j, and prod_{j=m}^{n-1} (j+2)/j = n(n+1)/(m(m+1)). ∎

**Proof of (a).** A vertex of degree >= k >= 2 in T_n is a core vertex of T_n. Either it is already core in
T_{m0} (at most C_{m0} vertices), or it becomes core at a step m with m0 <= m <= n-1, i.e. it is a leaf of T_m
that receives the new leaf; at each step at most one vertex becomes core, and C_{m+1} - C_m is the indicator
of that event. Such a vertex v has D_v(m+1) = 2, so by Lemma 1 and Markov's inequality
P(D_v(n) >= k | F_{m+1}) <= 2n/((m+1)k). Taking expectations and summing over m gives (a). ∎

**Lemma 2 (take-off).** Let m >= 3 and u a core vertex of T_m with l := L_u(m) >= 3 leaf neighbours and
c_u(m) <= 2 core neighbours. Then for every n >= m

    P( L_u(n) >= (l-1) n/(2m) | F_m )  >=  (l-1)^2 / (4 (l+2)(l+3))  >=  1/30.

*Proof.* L_u increases by one when W_j = u, which has conditional probability (L_u(j) + rho_u(j))/j with
rho_u >= 0 the contribution of the core neighbours, and decreases by one when one of the leaves of u receives
the new leaf, which has conditional probability L_u(j)/(j D_u(j)) <= 1/j. Hence
E[L_u(j+1) - L_u(j) | F_j] >= (L_u(j) - 1)/j, i.e. E[L_u(j+1) - 1 | F_j] >= (L_u(j) - 1)(j+1)/j. Taking
conditional expectations given F_m and iterating (the factors are positive, so the sign of L_u - 1 does not
matter), E[L_u(n) - 1 | F_m] >= (l - 1) n/m, so E[L_u(n) | F_m] >= (l-1) n/m.
By Lemma 1, E[L_u(n)^2 | F_m] <= E[D_u(n)(D_u(n)+1) | F_m] <= (l+2)(l+3) n(n+1)/(m(m+1)) <= (l+2)(l+3)(n/m)^2,
using D_u(m) <= l + 2 and (n+1)/(m+1) <= n/m.
By the Paley-Zygmund inequality for the nonnegative variable L_u(n),
P(L_u(n) >= E[L_u(n) | F_m]/2 | F_m) >= (1/4) E[L_u(n) | F_m]^2 / E[L_u(n)^2 | F_m] >= (l-1)^2/(4(l+2)(l+3)).
The last expression is increasing in l >= 3 and equals 4/120 = 1/30 at l = 3. ∎

Call a core vertex u of a tree **eligible** if (c_u <= 2 and L_u >= 1), or if (L_u, c_u) = (0, 2) and u has a
core neighbour w with (L_w, c_w) = (0, 2).

**Lemma 3 (a quarter of the core is eligible).** A tree with C >= 1 core vertices has at least C/4 eligible
vertices. (This counting lemma is the one used by the Astra lead in its existence proof.)

*Proof.* If C = 1 the core vertex has c = 0 and L >= 2. Let C >= 2; the core vertices span a tree. Let n_1 be the
number of core vertices with c = 1 (they have L >= 1, so they are eligible), n_3 the number with c >= 3, and I
the number of core vertices with (L, c) = (0, 2) that have no neighbour of that kind. A non-eligible vertex
has c >= 3 or is counted in I. The I vertices are pairwise non-adjacent and each has exactly two core edges,
so 2I <= C - 1. In a tree n_3 <= n_1 - 2. Hence the number E of eligible vertices satisfies
E >= C - n_3 - (C-1)/2 and E >= n_1 >= n_3 + 2; adding, 2E >= (C+1)/2 + 2. ∎

**Lemma 4 (activation).** Let m' >= 3 and let u be an eligible core vertex of T_{m'}. Then

    P( L_u(8m') >= 3 and c_u(8m') <= 2 | F_{m'} )  >=  1/368.

*Proof.* The V_j, j >= m', are independent, V_j uniform on the vertices of T_j, and V_j is independent of F_j.
Let G be the event that V_j != u for all m' <= j < 8m'; P(G) = prod_{j=m'}^{8m'-1} (1 - 1/j) = (m'-1)/(8m'-1)
>= 2/23. On G the vertex u never chooses a neighbour, so no leaf of u is ever promoted during [m', 8m'):
L_u does not decrease, c_u does not change, and a vertex that is a leaf of u stays a leaf of u. Conditionally
on G the V_j are still independent and V_j is uniform on the j - 1 vertices other than u.

*Case A: c_u <= 2 and L_u >= 1.* Fix a leaf x of u in T_{m'}. Each j with V_j = x gives W_j = u, i.e. one more
leaf at u. Given G, P(V_j = x for some j in [m', 3m')) = 1 - prod_{j=m'}^{3m'-1} (j-2)/(j-1) = 1 - (m'-2)/(3m'-2)
>= 2/3, and P(V_j = x for some j in [3m', 8m')) = 1 - (3m'-2)/(8m'-2) >= 5/8; the two events are independent
given G. On their intersection with G, L_u(8m') >= 3 and c_u(8m') = c_u(m') <= 2. The probability is at least
(2/23)(2/3)(5/8) = 5/138.

*Case B: (L_u, c_u) = (0, 2) and u has a core neighbour w with (L_w, c_w) = (0, 2).* Let z be the other
neighbour of w. Given G, consider the event H1 that V_j != z for all j in [m', 2m') and V_j = w for some j in
[m', 2m'). Its conditional probability is
prod_{j=m'}^{2m'-1} (j-2)/(j-1) - prod_{j=m'}^{2m'-1} (j-3)/(j-1)
= (m'-2)/(2m'-2) * (1 - (m'-3)/(2m'-3)) = (m'-2) m' / ((2m'-2)(2m'-3)) >= 1/4   for m' >= 3.
On G and H1 let tau be the first j >= m' with V_j = w. Before tau neither u nor z has been chosen as V, so no
leaf was attached to w and w still has exactly the two neighbours u and z; hence W_tau = u with conditional
probability 1/2. Then u has a leaf x (the vertex created at step tau), which stays a leaf on G. As in case A,
given G the events {V_j = x for some j in [2m', 4m')} and {V_j = x for some j in [4m', 8m')} have probabilities
1 - (2m'-2)/(4m'-2) >= 1/2 and 1 - (4m'-2)/(8m'-2) >= 1/2 and are independent of each other and of what
happened before 2m'. On the intersection, L_u(8m') >= 3 and c_u(8m') = 2. The probability is at least
(2/23)(1/4)(1/2)(1/4) = 1/368. ∎

**Proof of (b).** Let m' >= 3, k >= 2, n >= 8m'k. For an eligible vertex u of T_{m'}, Lemma 4 and then Lemma 2
(applied at time m = 8m', with l >= 3, so that (l-1)n/(2m) >= n/(8m') >= k) give
P(D_u(n) >= k | F_{m'}) >= P(L_u(n) >= k | F_{m'}) >= (1/368)(1/30). By Lemma 3 there are at least C_{m'}/4
eligible vertices in T_{m'}. Hence E[N_k(n)] >= E[C_{m'}]/(4 * 368 * 30) = a_{m'}/44160. ∎

**Proof of (c).** Upper bound: take m0 = ceil(2n/k) in (a); 3 <= m0 <= n because 48 <= 2n/k <= n. Put
b_j = 1/j. If m0 = n (which happens exactly when k = 2) the sum in (a) is empty. If m0 < n, summation by parts
gives

    sum_{m=m0}^{n-1} (a_{m+1} - a_m) b_{m+1} = a_n b_n - a_{m0} b_{m0+1} + sum_{j=m0+1}^{n-1} a_j (b_j - b_{j+1}).

In both cases, since all a_j > 0 and b_j - b_{j+1} = 1/(j(j+1)) <= j^(-2),

    sum_{m=m0}^{n-1} (a_{m+1} - a_m) b_{m+1} <= a_n b_n + sum_{j=m0+1}^{n-1} a_j (b_j - b_{j+1})
        <= A n^(lambda-1) + A sum_{j>m0} j^(lambda-2) <= A m0^(lambda-1) (1 + 1/(1-lambda)),

where the last step uses n >= m0, lambda < 1 and sum_{j>m0} j^(lambda-2) <= m0^(lambda-1)/(1-lambda);

so E[N_k(n)] <= A m0^lambda + (2n/k) A m0^(lambda-1)(1 + 1/(1-lambda)) <= A (2 + 1/(1-lambda)) m0^lambda, using
2n/k <= m0; and m0 <= 2n/k + 1 <= 3n/k. Lower bound: take m' = floor(n/(8k)) >= 3 in (b); then
m' >= n/(8k) - 1 >= (2/3) n/(8k) = n/(12k). ∎

**Proof of (d).** For every eps with 0 < eps < min(mu, 1 - mu) there are a_eps, A_eps > 0 with
a_eps m^(mu-eps) <= a_m <= A_eps m^(mu+eps) for all m >= 3 (a_m >= 1). Apply (c) with delta = mu - eps > 0 and
lambda = mu + eps < 1, take logarithms, divide by log(n/k) and let n/k -> infinity, then eps -> 0. ∎

## Remarks

1. (a) and (b) are statements about means. For fixed k the same activation argument run in blocks, as in the
   Astra lead's almost-sure proof, would give an almost-sure version; not done here.
2. The constant 1/44160 is far from optimal; only its positivity matters for the exponent.
3. Numerical illustration (not a proof, one run with C = 7.4e6 at n = 1.0e12, `outputs/hist_t27.txt`): the
   fraction of core vertices with degree >= k is 0.616, 0.338, 0.206, 0.131, 0.114 for k = 3, 6, 12, 24, 30,
   with local log-slope -0.64 between k = 16 and 24 and -0.62 between 24 and 30, slowly approaching -mu.
