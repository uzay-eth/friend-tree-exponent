# From uniform asymptotics of means to the almost-sure limit: a complete proof written by the Fable lead

**Author:** Fable lead (Claude Fable 5.2). **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** Claude checked the separate GPT probabilistic reconstruction line by line; Gemini independently reviewed Astra 6, not this exact text. This proof takes the uniform mean asymptotic and lower comparison as hypotheses. See [Astra 6–7 reconstructions](../README.md#astra-6-7-reconstructions).

*Written by the Claude-family lead as its own proof of the probabilistic half of Astra's claim 6, after the
reconstruction described in the [public review](../../../reviews/attacks/fable-on-astra-6b.md). It follows the route of the [blind proof R6P](../checks/astra6-reconstruction/blind-proof-probabilistic-half.md), which the lead had checked line by line; the lead had
also read the claimant's proof. So this text is a second-family write-up, not a blind derivation. The analytic half
(the hypotheses (H) and (H') below) is NOT proved here; see `fable-on-astra-6b.md` for its status.*

## 1. Setting, inputs, statement

**Model.** Continuous time: every vertex carries an independent rate-one clock; when the clock of v rings, v chooses a
uniform neighbour w and a new leaf is attached to w. (T_t) is the tree at time t, E_T the law started from a finite tree
T, n(T) its number of vertices. For k >= 2, N_k(T) is the number of vertices of degree at least k; C = N_2 is the number
of core vertices. F_k(u, T) := E_T N_k(T_u). For a vertex x of T, T + x is T with one new leaf attached at x. The
attachment rate at x is q_x(T) = sum over neighbours y of x of 1/deg(y); sum_x q_x(T) = n(T).

**Inputs.**
(I1) Round 1: 1/2 < mu < 1, and from the three-vertex path, log C(T_t)/t -> mu almost surely.
(H) There are sigma < mu, K < infinity, a function Psi > 0 on finite trees with at least one core vertex, and numbers
    1 = l_2 >= l_3 >= ... >= 0 such that for all such T, all k >= 2, all u >= 0:
        |F_k(u, T) - l_k e^(mu u) Psi(T)| <= K C(T) e^(sigma u),      Psi(T) <= K C(T).
(H') (used only in Step 9, large values) There are c > 0, c' < infinity with Psi(T) >= c C(T) - c'.

**Theorem.** Start from one edge. (a) C_n / n^mu -> X almost surely, with 0 < X < infinity almost surely.
(b) P(X < eps) > 0 for every eps > 0, and, under (H'), P(X > M) > 0 for every M. (c) Almost surely, for all k >= 2,
N_k(n)/C_n -> l_k. (Here C_n, N_k(n) refer to the tree with n vertices.)

Put delta := mu - sigma > 0. Constants b_1, b_2, ... below depend only on mu, sigma, K (and on k where stated).

## 2. Proof

**Step 0 (clocks).** At size j the next ring comes after an Exp(j) time, independent of which vertex rings and of its
choice; so the sequence of holding times is independent of the sequence of trees (the jump chain), and the jump chain
is the discrete friend tree. From a tree with m vertices the size n(t) is a Yule process: E_T n(t) = m e^t,
E_T n(t)^2 <= m(m+1) e^(2t), and n(t) e^-t -> W almost surely with W Gamma(m, 1)-distributed, in particular
0 < W < infinity, and E W^mu = Gamma(m + mu)/Gamma(m). W is a function of the holding times.
All functionals below are dominated by n(t), so every expectation is finite.

**Step 1 (one extra leaf: coupling).** Let T have at least two vertices and x be a vertex of T. Run T^- := T and
T^+ := T + x together as follows. At every time some vertices of the two trees are *paired* (a bijection between a set
of vertices of T^- and a set of vertices of T^+ that is a graph isomorphism of the induced subgraphs); all other
vertices are *unmatched*. Initially all vertices of T are paired with their copies and the extra leaf is the only
unmatched vertex. Each pair has one clock, each unmatched vertex its own. When the clock of an unmatched vertex rings,
it acts in its own tree and its new leaf is unmatched. When the clock of a pair v rings, let d_0 be the number of its
paired neighbours and a, b the numbers of its unmatched neighbours in T^-, T^+. With probability
min(1/(d_0+a), 1/(d_0+b)) for each paired neighbour y, both copies choose y and the two new leaves are paired; with
the remaining probability p_v = max(a,b)/(d_0 + max(a,b)) the two copies choose from the residual distributions and
both new leaves are unmatched. Each tree has the correct law, and the pairing remains an isomorphism.

Let B be the number of unmatched vertices. The paired vertices induce a subtree with at least two vertices, so d_0 >= 1
and p_v <= max(a,b)/(1 + max(a,b)) <= (a + b)/2. An unmatched vertex was born at a paired vertex or at an unmatched
one; hence every connected component of unmatched vertices is joined to the paired subtree by exactly one edge, and
sum_v (a_v + b_v) = number of such components <= B. The rate of increase of B is therefore at most
B + sum_v 2 p_v <= 2B, and

    E B(u) <= e^(2u).                                                                                  (1)

At any time, |N_k(T^+) - N_k(T^-)| <= B for k >= 2: a paired vertex whose two degrees differ has an unmatched
neighbour, which is the attachment point of a component that is attached to no other paired vertex; an unmatched
component with s vertices contains at most s - 1 vertices of degree >= 2 (its vertex farthest from the attachment
edge is a leaf, and if s = 1 its only vertex is a leaf); so each component accounts for at most 1 + (s - 1) = s units
of the difference. With (1):

    |F_k(u, T + x) - F_k(u, T)| <= e^(2u).                                                             (2)

**Step 2 (refinement at vertices of large degree).** Let d be the degree of x in T. Until the first ring of a pair
that fails to be matched (time S), the unmatched vertices are exactly the extra leaves at x in T^+: an unmatched leaf
can only give a leaf to x. Their number L(s) is dominated by a Yule process started from 1, so
E[L(s)(L(s) + 2); s < S] <= 2 e^(2s) + e^s. Before S only the pair x has unmatched neighbours, and a ring of x fails with
probability L/(d_0 + L) <= L/(d + 1), because d_0 >= d. Before S the two trees differ only in the degree of x, so
|N_k(T^+) - N_k(T^-)| <= 1{d < k}. At S, B = L + 2, and by (1) applied from time S on, the expected number of unmatched
vertices at time u is at most (L(S-) + 2) e^(2(u - S)). Hence

    |F_k(u, T + x) - F_k(u, T)| <= 1{d < k} + (1/(d+1)) int_0^u (2 e^(2s) + e^s) e^(2(u-s)) ds
                                 = 1{d < k} + e^(2u) (2u + 1 - e^-u)/(d + 1).                            (3)

**Step 3 (squared influences sum to a multiple of the core count).** For 0 <= u <= h and a tree T with C(T) >= 1:

    sum_x q_x(T) |F_k(u, T + x) - F_k(u, T)|^2 <= A C(T),      A = A(h, k) := e^(4h) [k + (2h + 1)^2].   (4)

Indeed: the leaves x have q_x = 1/deg(parent), so their q_x sum to sum over cores y of L_y/deg(y) <= C(T), and (2)
applies; the cores x of degree < k have q_x <= deg(x) <= k - 1 and are at most C(T) in number, and (2) applies; the
cores x of degree d >= k have q_x <= d and, by (3), squared influence at most e^(4u)(2u+1)^2/(d+1)^2, so each
contributes at most e^(4h)(2h+1)^2.

**Step 4 (the harmonic function and a martingale).** By (H) with k = 2: Psi(T) = lim_u e^(-mu u) F_2(u, T), and

    E_T C(T_r) <= e^(mu r) Psi(T) + K C(T) e^(sigma r) <= 2K C(T) e^(mu r).                             (5)

By the Markov property, e^(-mu u) F_2(h + u, T) = E_T[e^(-mu u) F_2(u, T_h)]; the integrand is at most 2K C(T_h), which
is integrable, and converges to Psi(T_h); so E_T Psi(T_h) = e^(mu h) Psi(T), and M_t := e^(-mu t) Psi(T_t) is a
non-negative martingale. It converges almost surely to a finite limit Y.

**Step 5 (concentration over a window).** Fix s, h >= 0 and k. For s <= r <= s + h let
G_r := E[N_k(T_{s+h}) | history to r] = F_k(s + h - r, T_r). This is a square-integrable martingale of a pure-jump
process; it jumps exactly when a leaf is added, by F_k(s+h-r, T_r + x) - F_k(s+h-r, T_r) if the leaf is added at x
(rate q_x(T_r)), and is continuous in between. So its variance over the window is the expectation of the sum of its
squared jumps, and by (4), (5):

    E_T [N_k(T_{s+h}) - F_k(h, T_s)]^2 = E_T int_s^{s+h} sum_x q_x(T_r) [F_k(s+h-r, T_r + x) - F_k(s+h-r, T_r)]^2 dr
                                       <= A(h, k) int_s^{s+h} E_T C(T_r) dr <= 2K A(h, k) h C(T) e^(mu (s+h)).    (6)

(With (2) alone the bound would be e^(4h) h n(T) e^(s+h), which after division by e^(2 mu (s+h)) is small only because
2 mu > 1; this is the claimant's version. With (4) the condition mu > 1/2 is not needed at this point.)

**Step 6 (convergence along meshes, then for all times).** Fix beta in (0, mu/8). For t >= 1 let s = (1 - beta) t,
h = beta t. By (6),

    E_T [e^(-mu t) (N_k(T_t) - F_k(h, T_s))]^2 <= b_1(k) C(T) (1 + t)^3 e^(-(mu - 4 beta) t),             (7)

and by (H) and (5),

    |e^(-mu t) F_k(h, T_s) - l_k M_s| <= K e^(-delta h) e^(-mu s) C(T_s),   whose expectation is <= 2 K^2 C(T) e^(-delta beta t).   (8)

Choose a > 0 with 4a < min(mu - 4 beta, delta beta). By Chebyshev for (7), Markov for (8) and Borel-Cantelli, for
every mesh t in (1/j) N, almost surely |e^(-mu t) N_k(T_t) - l_k M_{(1-beta) t}| <= 2 e^(-a t) for all large mesh points.
Since M converges to Y, e^(-mu t) N_k(T_t) -> l_k Y along every mesh, almost surely, for all j and k simultaneously.
N_k(T_t) is non-decreasing in t, so between two points of the mesh 1/j the normalised count lies between the two mesh
values times e^(-mu/j) and e^(mu/j); letting j -> infinity:

    e^(-mu t) N_k(T_t) -> l_k Y   as t -> infinity, almost surely, for all k >= 2.                       (9)

**Step 7 (rate in L^1; the limit is positive).** Apply (7) and (8) with k = 2 at the times t and t + 1, both with the
same s = (1 - beta) t (windows of length beta t and beta t + 1); by Cauchy-Schwarz in (7),

    E_T |e^(-mu (t+1)) C(T_{t+1}) - e^(-mu t) C(T_t)| <= b_2 C(T) e^(-2a t)      (t >= 1).

Summing over t, t + 1, ... and using (9) and Fatou, E_T |Y - e^(-mu t) C(T_t)| <= b_3 C(T) e^(-2a t). In particular
E_T Y = Psi(T) (let t -> infinity in E_T[e^(-mu t) C(T_t)] and use (H)). By Markov and Borel-Cantelli over integer t,
almost surely |Y - e^(-mu t) C(T_t)| <= e^(-a t) for all large integers t. On the event {Y = 0} this gives
C(T_t) <= e^((mu - a) t) for all large integers t, hence limsup log C(T_t)/t <= mu - a < mu, which by (I1) has
probability zero when the start is the three-vertex path. So Y > 0 almost surely from the three-vertex path.

**Step 8 (fixed size; proof of (a) and (c)).** From the edge, the first ring produces the three-vertex path after an
Exp(2) time; the limits above from the path are multiplied by a common positive factor, which cancels below. Let tau_j
be the time at which the size reaches j; tau_j -> infinity and j e^(-tau_j) -> W. Evaluating (9) at the times tau_j:

    C_j / j^mu = e^(-mu tau_j) C(T_{tau_j}) / (j e^(-tau_j))^mu -> Y / W^mu =: X,

and 0 < X < infinity almost surely by Step 7 and Step 0. This is (a). Dividing (9) for k by (9) for k = 2 at the times
tau_j gives N_k(j)/C_j -> l_k, simultaneously for all k: this is (c).

**Step 9 (proof of (b)).** The same arguments apply from any deterministic finite tree T with a core; write X_T for
the limit. X_T is a function of the jump chain (it is the limit of C_j / j^mu), W of the holding times, so they are
independent, Y = W^mu X_T, and by Step 7

    E X_T = Psi(T) / E W^mu = Psi(T) Gamma(m) / Gamma(m + mu),      m = n(T).                           (10)

By the Markov property of the jump chain, conditionally on the tree with m vertices being T, the limit X has the law
of X_T. *Small values:* the star with m vertices is reached with positive probability (probability 2/(m-1)); it has
one core vertex, so Psi <= K, and E X_T <= K Gamma(m)/Gamma(m + mu) -> 0 as m -> infinity; choose m with E X_T < eps/2;
then P(X_T < eps) >= 1/2. *Large values:* the path with m vertices is reached with positive probability; it has m - 2
core vertices, so by (H') Psi >= c (m - 2) - c' and E X_T >= (c (m - 2) - c') Gamma(m)/Gamma(m + mu), which tends to
infinity like m^(1 - mu). Choose m with E X_T > M; then P(X_T > M) > 0. This proves (b).

## 3. Remarks

1. Where each hypothesis enters: (H) in Steps 4, 6, 7 (both the limit and the uniformity of the error in the starting
   tree, through E[e^(-mu s) C(T_s)] in (8)); the weight C(T) in the error term of (H) is what makes (8) work with a
   window proportional to t; (I1) only in Step 7, for positivity; (H') only for the large values in Step 9. The blind
   proof R6P obtains the large values without (H') from trees made of a hub with many three-vertex branches.
2. Steps 1 to 3 use nothing but the dynamics. The claimant uses the cruder bound 3 e^(2u) in place of (2), (3), which is
   sufficient because 2 mu > 1.
3. The exact bound (3) and the inequality (4) were checked in exact arithmetic on all trees with at most 8 vertices by
   the blind prover (21,125 influence inequalities, 3,055 inequalities of type (4)).
