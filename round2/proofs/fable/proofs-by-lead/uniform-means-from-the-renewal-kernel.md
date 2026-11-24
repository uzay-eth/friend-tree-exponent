# Uniform asymptotics of means from the marked-fringe renewal structure: a proof written by the Fable lead

**Author:** Fable lead (Claude Fable 5.2). **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** GPT-6.2 Sol Ultrafast cold referee checked the weight inequalities and generator and supplied the strict-positivity repair. Uses the marked-fringe renewal structure and a contour shift. The Gemini audit addresses Astra 6, not this exact write-up. See [Astra 6–7 reconstructions](../README.md#astra-6-7-reconstructions).

*Written by the Claude-family lead as its own proof of the analytic half of Astra's claim 6 (hypotheses (H) and (H') of
the [probabilistic companion](almost-sure-limit-from-uniform-means.md)), after the reconstruction described in the [public review](../../../reviews/attacks/fable-on-astra-6b.md).
The construction (marks, copies, the weight W and its constants) is the claimant's. The lead had read the claimant's
proof and checked the blind proofs R2 and R6A (`checks/astra6-reconstruction/`); this text is a second-family write-up
in the lead's own words, not a blind derivation. Standard theorems are cited by name with the hypotheses checked.*

## 0. Setting

**Physical process.** Continuous-time friend tree: every vertex rings at rate 1 and gives a new leaf to a uniformly
chosen neighbour. A core is a vertex of degree >= 2. A *promotion* is a ring of a core x that chooses one of its leaf
neighbours, which thereby becomes a core with one leaf. Round-1 input: (I0) for the tree started from the
three-vertex path T_3, (1/t) log E C(T_t) -> mu, and 113/200 <= mu <= 1133/2000. (This is the continuous-time form of
E C_n = n^(mu+o(1)): the size is an independent Yule process, with n(t) e^-t bounded in every L^p.)

**Marked trees.** Let m(D) = min(D, 20). The tree is rooted at a core o (the *founder*; for T_3 its centre). Every other
core v carries a mark a_v in {1/2, ..., 1/20}: when v is promoted by its parent u, a_v := 1/m(D_u), with D_u the
degree of u at that moment (a promotion does not change D_u). Degrees never decrease, so at all later times

    1/D_u <= a_v      (u the parent of v).                                                             (0.1)

For an arbitrary starting tree T we allow any marks with (0.1) ("compatible marks").

**Types.** A *type* S is a finite rooted tree of cores, each with a number of leaves, with a mark a_S in
{0, 1/2, ..., 1/20} at the root, marks satisfying (0.1) at the other cores, and a *parent slot* at the root: the slot
counts in the degree of the root and in its number c of core neighbours, absorbs what is sent to it and never sends.
For a non-founder core v of a marked tree, the *fringe* S_v is the subtree of v and its descendants with their leaves
and marks, with the slot in place of the parent. For a type S write D, L, c for the degree, the number of leaves and
the number of core neighbours (slot included) of its root, p_S = 1/D, b_S = L/D, n(S) for its number of vertices.
Y denotes the fringe of a core child of the root, Y^+ the same with one more leaf at its root, and J_i (i = 2, ..., 20)
the *newborn*: one core with mark 1/i, one leaf and the slot.

**Intrinsic dynamics L^0.** Inside a type the physical rules act: a leaf rings and gives its parent a leaf; a core x
rings and chooses uniformly among its D_x neighbours: a leaf (promotion; the new core gets the mark 1/m(D_x)), a core
child (which gains a leaf), or its parent (which gains a leaf; if x is the root the parent is the slot and nothing
happens).

**Particle systems.** *Residual system:* each particle carries a type; it evolves by L^0, dies at rate a_S, and for
every core child Y of its root it emits, keeping itself, a particle Y^+ at rate p_S and a particle Y at rate
a_Y - p_S (non-negative by (0.1)). *Full system:* in addition it emits a newborn J_{m(D)} at rate b_S (D the degree of
its root). R_t f(S) and P_t f(S) denote the expected sum of f over the particles alive at time t in the residual and
the full system started from one particle S. Their generators are

    R f(S) = L^0 f(S) - a_S f(S) + sum_Y [(a_Y - p_S) f(Y) + p_S f(Y^+)],        (R + B) f,   B f(S) = b_S f(J_{m(D)}).

**Weight.** d = 3/20, q = 9/25, rho = q + d = 51/100, z(a) = a/(a + d); h(1) = 1, h(2) = 17/20, h(3) = 39/50,
h(c) = 7/10 for c >= 4, and h(0) = 2 (used for a founder only). For a core x of a type, w_x is the product of z(a_y)
over the cores y on the path from the root (excluded) to x (included). Put

    W(S) = (a_S + d)^-1 sum over cores x of S of w_x h(c_x),        Z(T) = sum over cores x of T of w_x h(c_x)  (rooted tree T).

Note rho < 113/200 <= mu. We write ||F||_W = sup_S |F(S)|/W(S).

## 1. The weight (exact inequalities)

**Lemma 1.** (a) 14/13 <= W(S) <= (number of cores of S)/d. (b) (a_S + d) W(S) = h(c_root) + sum_Y a_Y W(Y), and
W(Y^+) = W(Y). (c) L^0 W <= q W, and for the physical generator on rooted marked trees, L Z <= q Z.
(d) R W <= rho W - h(c_root) <= rho W, and B W <= 1/d <= W/d.

*Proof.* (a) w_x <= 1 and h <= 1 on cores with c >= 1, and a_S + d >= d; the root alone contributes at least
(7/10)/(1/2 + 3/20) = 14/13. (b) For x in the fringe Y the weight computed in S is z(a_Y) times the weight computed
in Y, and z(a)(a + d) = a, so the cores of Y contribute z(a_Y) (a_Y + d) W(Y) = a_Y W(Y) to (a_S + d) W(S); the root
contributes h(c_root). A leaf at the root of Y changes neither marks nor numbers of core neighbours. (c) A leaf gain changes no mark and no c_x: W is unchanged. A promotion at a core x of degree k with
c core neighbours changes the sum by w_x [h(c+1) - h(c) + z(1/m(k))] (x gets one more core neighbour; the new core
has c = 1, h = 1, and path factor w_x z(1/m(k))), and has rate (k - c)/k. So it suffices that

    ((k - c)/k) [h(c+1) - h(c) + z(1/m(k))] <= q h(c)      for all integers c >= 0, k >= max(2, c).            (1.1)

For c = 0 the bracket is -1 + z < 0. For c = 1, 2, 3 and c >= 4 the left side, as a function of k <= 20, is maximal at
k = 3, 6, 8, 11 with values 313/870, 289/950, 103/440, 140/583, against q h(c) = 9/25, 153/500, 351/1250, 63/250; for
k >= 20, z = 1/4 and the left side is at most 1/4 + h(c+1) - h(c) <= 1/4 < 63/250 <= q h(c) when c >= 4, and at most
1/10, 9/50, 17/100 for c = 1, 2, 3. (Recomputed exactly by the lead; the tightest case is 313/870 < 9/25, margin
1/4350.) Summing (1.1) over the cores with the factors w_x gives both statements of (c). (d) Since W(Y^+) = W(Y),
R W = L^0 W - a_S W + sum_Y a_Y W(Y), and by (b) sum_Y a_Y W(Y) = (a_S + d) W - h(c_root); hence
R W = L^0 W + d W - h(c_root) <= rho W - h(c_root). And B W = b_S W(J) = b_S/(1/m(D) + d) <= 1/d. []

**Lemma 2 (consequences).** Both particle systems are non-explosive, and for t >= 0

    R_t W <= e^(rho t) W,      P_t W <= e^((rho + 1/d) t) W,      E_T Z(T_t) <= e^(q t) Z(T) <= 2 e^(q t) C(T).        (1.2)

*Proof.* Standard: for the residual system the sum of W over the living particles has drift at most rho times
itself (Lemma 1(d)); stopped at the first time the number of particles or the size of a particle exceeds N,
e^(-rho t) times that sum is a non-negative supermartingale; since W >= 1 per particle, the number of particles
exceeds N before time t with probability at most e^(rho t) W(S)/N; the size of each particle is dominated by a Yule
process; letting N -> infinity (Fatou) gives non-explosion and the first bound. The same with rho + 1/d for the full
system, and with q for Z, where h <= 2 gives Z <= 2 C. []

## 2. The first-moment identity

For a rooted marked tree T and a function F on types put G_F(T) = sum over non-founder cores v of F(S_v), and

    I(T) F = b_o F(J_{m(D_o)}) + sum over core children Y of o of [(a_Y - p_o) F(Y) + p_o F(Y^+)],

with p_o = 1/D_o, b_o = L_o/D_o for the founder (no slot). Put I_{T,u} F = E_T I(T_u) F.

**Lemma 3.** L G_F = G_{(R+B)F} + I(.) F for the physical generator L, and for |f| <= W and every rooted marked tree T,

    E_T G_f(T_t) = sum over non-founder cores v of T of P_t f(S_v) + int_0^t I_{T,u} P_{t-u} f du.              (2.1)

Moreover I(T) W <= 1/d + Z(T), so I_{T,u} W <= 1/d + e^(q u) Z(T) <= (1/d + 2) e^(q u) C(T).

*Proof.* Consider what a physical event does to the sum G_F. The fringe S_v of a non-founder core v changes by
every event inside it exactly as under L^0, and by one more event: a leaf sent to v by its parent u, at rate 1/D_u,
which turns S_v into S_v^+. (A ring of v that chooses its parent changes nothing in S_v, like the slot.) A promotion at
a core x (founder or not), rate L_x/D_x, in addition creates the new fringe J_{m(D_x)}. Therefore

    L G_F(T) = sum_{v != o} L^0 F(S_v) + sum_{v != o} (1/D_{u(v)}) [F(S_v^+) - F(S_v)] + sum over all cores x of b_x F(J_{m(D_x)}).

On the other side, in sum_{v != o} R F(S_v) + I(T) F the term -a_v F(S_v) of the particle v cancels against the terms
(a_v - p) F(S_v) + p F(S_v^+) = a_v F(S_v) + p [F(S_v^+) - F(S_v)] contributed by its parent u (p = 1/D_u), whether u is
a non-founder core (through R F(S_u)) or the founder (through I(T) F); what remains is exactly the middle sum above, and
the terms B F(S_x) and b_o F(J) give the last sum. This is the generator identity. For (2.1): u -> P_u f(S) solves the
backward equation d/du P_u f = (R + B) P_u f (first-jump decomposition; the rates out of a fixed type are finite and
|P_u f| <= e^((rho + 1/d) u) W by (1.2)). Apply Dynkin's formula to u -> G_{P_{t-u} f}(T_u) on [0, t]: the time derivative
cancels G_{(R+B) P_{t-u} f}, leaving I(T_u) P_{t-u} f. The integrability needed is polynomial in the size of the tree
(each of at most n(T_u) terms is at most a constant times n(T_u), and so are the rates), and the size is a Yule
process with all moments. Finally I(T) W = b_o W(J) + sum_Y a_Y W(Y) <= 1/d + (Z(T) - h(c_o)) by Lemma 1(b) applied at
the founder, and (1.2). []

## 3. The newborn kernel

For a vector g_0 on {2, ..., 20} put B_0 g_0(S) = b_S g_0(m(D_root(S))), so that |B_0 g_0| <= max|g_0| <= max|g_0| W.
Define the 19 x 19 matrices

    k_ij(t) = R_t [B_0 e_j](J_i),      K_ij(z) = int_0^infinity e^(-z t) k_ij(t) dt      (Re z > rho):

k_ij(t) dt is the expected number of newborns of label j emitted in dt by the residual family of one newborn of
label i.

**Lemma 4.** (a) 0 <= k_ij(t) <= e^(rho t)/d. (b) k is locally Lipschitz with |k'_ij(t)| <= 42 e^(rho t)/d almost
everywhere; hence for Re z >= rho + eps, |K_ij(z)| <= (1 + 42/eps)/(d |z|). (c) k_ij(t) > 0 for all t > 0, and
k_i2(t) >= (1/2) e^(-(1/i + 3/2) t). (d) Branching identity: for |f| <= W,

    P_t f = R_t f + int_0^t R_u [B_0 x(t - u)] du,      x_i(t) := P_t f(J_i),

so x = h_f + k * x with h_f,i(t) = R_t f(J_i), and x = h_f + v * h_f with v = sum_{j >= 1} k^{*j} (convolution powers;
the series converges since |k^{*j}(t)| <= e^(rho t) t^(j-1)/(d^j (j-1)!) entrywise up to the factor 19^(j-1)).

*Proof.* (a) Lemma 2 and W(J_i) = 1/(1/i + d) <= 1/d. (b) Let F = B_0 g_0 with |g_0| <= 1. We bound |R F(S)|. Under L^0,
F(S) changes only when the root gains a leaf (rate at most D; the change of b = L/D is at most 1/(D+1), and the label
m(D) changes only if D <= 19, a change of at most 2 at rate at most 19) or promotes a leaf (rate b <= 1, change at most
1): at most 1 + 38 + 1. The death term is at most 1/2. The copies contribute at most
sum_Y a_Y <= (1/2 + d) sum_Y z(a_Y) <= (1/2 + d)(a_S + d) W(S)/(7/10) <= (13/20)^2 (10/7) W(S), because
(a_S + d) W(S) >= sum_Y z(a_Y) h(c_Y) >= (7/10) sum_Y z(a_Y). Altogether |R F| <= 42 W. Then R_t F = F + int_0^t R_s R F ds
gives |d/dt k(t) g_0| <= 42 R_t W(J_i) <= 42 e^(rho t)/d, and an integration by parts gives the bound on K(z)
(|k_ij(0)| <= 1/2). (c) From J_i the original particle can, without dying, have its leaf ring j - 2 times in a row
(reaching the star of degree j) with positive probability density at every time, and in that state it emits newborns
of label m(j) = j at rate (j-1)/j; for j = 20 use degree 20. At J_i the first change of state or death has hazard
1 + 1/2 + 1/i: the leaf rings at rate 1, the root rings and picks the leaf at rate 1/2 (a ring that picks the slot
changes nothing), and death has rate 1/i. As long as none of these has happened the particle is J_i, of degree 2, and
emits newborns of label 2 at rate b = 1/2; this gives the contribution (1/2) e^(-(1/i + 3/2) t) to k_i2(t), in every
row i. (d) Split the full system
according to whether a particle descends from a direct birth: those that do not form the residual system; each direct
birth of label j, emitted at intensity R_u[B_0 e_j], starts an independent full system. Iterate the equation
x = h_f + k * x; the remainder k^{*j} * x tends to 0 locally uniformly. []

**Lemma 5 (criticality).** For real s > rho, K(s) has positive entries, is continuous and entrywise strictly
decreasing in s, and its spectral radius r(K(s)) satisfies r(K(s)) < 1 for s > mu, r(K(mu)) = 1.

*Proof.* Positivity, continuity and strict decrease follow from Lemma 4 (a), (c). Let s > rho and
Lambda(s) = int e^(-s t) E_{T_3} G_1(T_t) dt in [0, infinity], the transform of the mean number of non-founder cores.
By (2.1) with f = 1 (the initial sum is empty) and Tonelli, Lambda(s) = int_0^infinity e^(-s u) I_u [hat P_s 1] du, with
hat P_s 1(S) = int e^(-s t) P_t 1(S) dt. On newborns, Lemma 4(d) with f = 1 gives x = sum_{j >= 0} k^{*j} * h in the
time domain (the term j = 0 is h); all entries are non-negative, so by Tonelli hat x = sum_{j >= 0} K(s)^j hat h in
[0, infinity]^19, where 1/(s + 1/2) <= hat h_i <= 1/(d (s - rho)). If r(K(s)) < 1 the series converges,
hat P_s 1 <= const W by Lemma 4(d) and (1.2), and Lambda(s) < infinity because I_u W <= 1/d + 2 e^(q u) and s > q. If
r(K(s)) >= 1, let g > 0 be the Perron vector of the positive matrix K(s); hat h >= c g, so every coordinate of hat x
is infinite; and I_u charges J_2 with positive weight for every u (with positive probability nothing has happened by
time u, and then the founder of T_3, of degree 2 with two leaves, emits J_2 at rate 1); so Lambda(s) = infinity. By
(I0), Lambda(s) < infinity for s > mu and = infinity for s < mu. Hence r(K(s)) < 1 for s > mu and r(K(s)) >= 1 for
rho < s < mu; by continuity r(K(mu)) = 1. []

Let g and l be the right and left Perron vectors of K(mu), g, l > 0, max g = 1, l g = 1, and
alpha := - l K'(mu) g = l [int t e^(-mu t) k(t) dt] g > 0.

**Lemma 6 (renewal remainder).** There are sigma in (rho, mu) and A < infinity such that, with
D := g l / alpha (a 19 x 19 matrix),  |v(t) - e^(mu t) D| <= A e^(sigma t) entrywise for all t >= 0.

*Proof.* (i) *The line Re z = mu.* Let z = mu + i tau, tau != 0, and suppose K(z) y = y, y != 0. Then
|y| <= K(mu) |y| entrywise; applying l > 0 gives equality, so |y| is a multiple of g, in particular y_2 != 0, and
for every i there is equality in |sum_j K_ij(z) y_j| <= sum_j K_ij(mu) |y_j|, which forces |K_i2(z)| = K_i2(mu). But
k_i2 has the density component (1/2) e^(-(1/i + 3/2) t) by Lemma 4(c), so
|K_i2(mu + i tau)| <= K_i2(mu) - (1/2)[1/A_i - 1/|A_i + i tau|] < K_i2(mu), A_i = mu + 1/i + 3/2. Contradiction: I - K(z)
is invertible on the line except at z = mu. (ii) *Re z > mu:* |K(z)| <= K(Re z) entrywise and r(K(Re z)) < 1, so
I - K(z) is invertible. (iii) *Near mu:* 1 is a simple eigenvalue of K(mu) (Perron-Frobenius, positive matrix), K is
analytic in Re z > rho, and the derivative of the eigenvalue at mu is l K'(mu) g = -alpha != 0; hence (I - K(z))^-1 has a
simple pole at mu with residue g l/alpha = D. (iv) *High frequencies:* by Lemma 4(b), for Re z >= sigma_0 := (rho + mu)/2
and |Im z| >= Y_0 (Y_0 large) the norm of K(z) is at most 1/2. (v) *A strip:* det(I - K(z)) is analytic and not
identically zero; in the compact rectangle sigma_0 <= Re z <= mu, |Im z| <= Y_0 it has finitely many zeros, and by (i)
the only one on the right edge is mu; so there is sigma in (sigma_0, mu) with no zero in sigma <= Re z <= mu other than
mu. (vi) *Inversion:* for real c > mu the transform of v is K(z)(I - K(z))^-1 (sum of the transforms of k^{*j}, convergent
because r(K(c)) < 1), and that of v - k is F(z) := K(z)^2 (I - K(z))^-1, which by (iv) and Lemma 4(b) is O(|z|^-2) on
sigma <= Re z <= c away from mu. v - k is continuous, so the Laplace inversion formula holds with an absolutely
convergent integral on Re z = c; the horizontal segments vanish by the O(|z|^-2) bound, and moving the line to
Re z = sigma crosses only the simple pole at mu, where e^(z t) F(z) has residue e^(mu t) K(mu)^2 g l/alpha = e^(mu t) D
(because K(mu) g = g). The integral over Re z = sigma is O(e^(sigma t)). So v(t) - k(t) = e^(mu t) D + O(e^(sigma t)),
and k(t) = O(e^(rho t)). []

## 4. Uniform asymptotics

For `|f| <= W` put `d(f) := l [int_0^infinity e^(-mu u) h_f(u) du] / alpha` and `H(S) := int_0^infinity e^(-mu u) R_u [B_0 g](S) du`.
By (1.2): |h_f(u)| <= e^(rho u)/d, so |d(f)| <= M_d := (sum_i l_i)/(alpha d (mu - rho)), and 0 <= H <= W/(mu - rho).

**Lemma 7 (planted states).** There is A_P < infinity such that for all types S, all |f| <= W and all t >= 0,

    |P_t f(S) - e^(mu t) H(S) d(f)| <= A_P e^(sigma t) W(S).                                             (4.1)

*Proof.* Vectors carry the maximum norm and matrices the maximum row sum, written |D|; the constant A of Lemma 6 is
enlarged by the factor 19 so that it bounds the matrix norm of the remainder. On newborns, x = h_f + v * h_f and
Lemma 6 give
x(t) = e^(mu t) D int_0^t e^(-mu u) h_f(u) du + O(e^(sigma t)) = e^(mu t) g d(f) + O(e^(sigma t)), with an error bounded by
A_J e^(sigma t), A_J := 1/d + A/(d (sigma - rho)) + |D|/(d (mu - rho)) (the three terms: h_f itself, the remainder of v
convolved with h_f, and the tail of the integral), uniformly in f. For a general type, by Lemma 4(d),

    P_t f(S) = R_t f(S) + int_0^t R_u [B_0 x(t - u)](S) du
             = e^(mu t) d(f) int_0^t e^(-mu u) R_u [B_0 g](S) du + (terms bounded by e^(rho t) W(S) + A_J int_0^t e^(rho u) e^(sigma (t-u)) du W(S)),

using |R_u [B_0 y]| <= max|y| R_u W <= max|y| e^(rho u) W. Completing the integral to infinity costs
M_d e^(mu t) int_t^infinity e^(-(mu - rho) u) du W(S) = M_d e^(rho t) W(S)/(mu - rho). So (4.1) holds with
A_P = 1 + A_J/(sigma - rho) + M_d/(mu - rho). The uniformity in S comes from R_u W <= e^(rho u) W alone. []

**Lemma 8 (sums of weights).** For every rooted marked tree T: sum over non-founder cores v of W(S_v) <= C_0 C(T),
C_0 := 1/(d (1 - 10/13)) = 260/9.

*Proof.* A core x contributes to W(S_v) for v = x and for each non-founder ancestor v of x, with the factor
(a_v + d)^-1 times the product of z over the path from v to x, at most (1/d)(10/13)^j at distance j, because
z(a) <= z(1/2) = 10/13; and h <= 1. Sum the geometric series. []

**Theorem 9 (= hypothesis (H)).** Define, for a rooted marked tree T with at least one core,

    Psi_0(T) := sum over non-founder cores v of H(S_v) + int_0^infinity e^(-mu u) I_{T,u} H du.

There is C_1 < infinity such that for all such T, all |f| <= W, all t >= 0:

    |E_T G_f(T_t) - e^(mu t) Psi_0(T) d(f)| <= C_1 C(T) e^(sigma t),     and     0 < Psi_0(T) <= C_2 C(T).         (4.2)

Consequently, with f_k(S) = 1{D_root(S) >= k}, d_k := d(f_k), Psi := d_2 Psi_0 and l_k := d_k/d_2:

    |E_T N_k(T_t) - l_k e^(mu t) Psi(T)| <= (C_1 + 1) C(T) e^(sigma t),      1 = l_2 >= l_3 >= ... > 0,       (4.3)

and Psi(T) = lim_t e^(-mu t) E_T C(T_t) does not depend on the choice of the root and of the marks.

*Proof.* In (2.1) insert (4.1): the sum over v is e^(mu t) d(f) sum_v H(S_v) up to A_P e^(sigma t) sum_v W(S_v)
<= A_P C_0 C(T) e^(sigma t) (Lemma 8); the integral is e^(mu t) d(f) int_0^t e^(-mu u) I_{T,u} H du up to
A_P int_0^t e^(sigma (t-u)) I_{T,u} W du <= A_P (1/d + 2) C(T) e^(sigma t)/(sigma - q) (Lemma 3; q < rho < sigma);
completing the integral to infinity costs M_d (1/d + 2) C(T) e^(q t)/((mu - rho)(mu - q)). This is (4.2), and the bound
on Psi_0 follows from H <= W/(mu - rho), Lemma 8 and Lemma 3. For (4.3): the root of S_v has the degree of v in the tree
(the slot stands for the parent), so N_k(T_t) = G_{f_k}(T_t) + 1{D_o(t) >= k}, and the founder costs at most 1 <= C(T).
0 <= f_k <= 1 <= W and f_k decreases in k, so d_k decreases; d_k > 0 because h_{f_k}(u) > 0 for u > 0 (the original
particle reaches degree k with positive probability) and l > 0. Strict positivity of Psi_0: H(J_i) =
(K(mu) g)_i = g_i > 0 for every label i, and for every u > 0 the founder has a leaf at time u with positive
probability (it has one at time 0, or else all its neighbours are cores and with positive probability the only event
up to time u is a ring of one of them that selects the founder); then b_o > 0, so the direct term of I_{T,u} H is
positive for every u > 0 and so is its integral in Psi_0(T). The last statement: E_T C(T_t) involves neither root nor
marks. []

**Theorem 10 (= hypothesis (H'): lower comparison).** There are c > 0, c' < infinity with Psi(T) >= c C(T) - c'.

*Proof.* Since I_{T,u} H >= 0, Psi_0(T) >= sum_{v != o} H(S_v). (i) If the root of S has c <= 4 and at least one leaf,
then H(S) >= g_min/(5 (mu + 3/2)), g_min := min g > 0: with probability at least e^(-3u/2) the original particle has
neither died (rate at most 1/2) nor had its root ring (rate 1) by time u; until then its root keeps c and has at least
one leaf, so b >= 1/5 and B_0 g >= g_min/5. (ii) If the root is bare (no leaf), has c <= 4, and has a core child Y that
is also bare with c_Y <= 4, then H(S) >= (1/(4 mu + 19)) times the bound of (i). Ignore the emitted particles and
follow the original one. Until the first change at one of the two bare vertices, or the death of the particle, the
child sends a leaf to the root at the constant rate lambda = 1/c_Y >= 1/4. The competing events (death, rate at most
1/2; a leaf sent to the root or to Y by another of their at most 8 core neighbours, each at rate at most 1/2) have
rates that may vary with what happens elsewhere, but their sum is at most 9/2 at all times. So the probability S(u)
that neither the send nor a competing event has occurred by time u satisfies S(u) >= e^(-(lambda + 9/2) u), and the
discounted probability that the send comes first is at least
int_0^infinity e^(-mu u) lambda S(u) du >= lambda/(mu + lambda + 9/2) >= (1/4)/(mu + 19/4) = 1/(4 mu + 19). At that
moment the root has a leaf and still c <= 4, and (i) applies by the strong Markov property. (iii) Counting.
The cores form a tree, so sum over cores of c_x = 2 (C - 1), and at most 2(C-1)/5 cores have c >= 5 ("high"). Among the
non-founder, non-high cores (at least C - 1 - 2(C-1)/5 = 3(C-1)/5 of them), those with a leaf satisfy (i). Each of the
others is bare, hence has c = D >= 2 and a core child; choose one child for each; the children chosen are distinct. If
the chosen child is high, discard the parent (at most 2(C-1)/5 discarded). If the child is not high and has a leaf,
credit the parent with the child's fringe (i); if it is bare, the parent satisfies (ii). At least (C-1)/5 cores
remain, each credited with a fringe S_v with H(S_v) >= h_B > 0 (the bound of (ii)), and each fringe is used at most
twice. Hence sum_{v != o} H(S_v) >= h_B (C - 1)/10. []

## 5. Remarks

0. Scope. Theorems 9 and 10 are for starting trees with at least one core (three or more vertices). For the single
   edge C = 0 and the bounds cannot hold with the factor C(T); condition on its first event (an exponential time of
   rate 2, after which the tree is T_3), which gives Psi(edge) = 2 Psi(T_3)/(2 + mu). The statement and the argument for
   strict positivity of Psi, and the wording of four steps (Lemma 4(c), the series in Lemma 5, the norms in Lemma 7,
   case (ii) of Theorem 10), were corrected after a hostile reading by a GPT-family referee
   (`checks/lead-proof-uniform-means-referee-gpt-6.2-sol-ultrafast.md`), who found no other gap and confirmed the
   inequalities (1.1) in all 228 finite cases, the generator identity of Lemma 3 on 2,310 exact comparisons, and the
   drift of the weight on 5,820 marked types.

1. Constants. Everything except sigma and A is explicit; sigma comes from a compactness argument in Lemma 6(v). The
   blind proof R6A makes sigma explicit (with an astronomically small mu - sigma); nothing in the probabilistic half
   needs a value.
2. Where the round-1 input enters: only in Lemma 5, through the growth rate of E C(T_t), and through mu > 51/100.
3. The one place where the blind proof R6A erred is Lemma 6(i): the density component is in column 2 of every row
   (a newborn of any label has degree 2 until its first event), not on the diagonal.
4. Exact checks of the inequalities of sections 1 and 2 exist in three independent programs (blind prover R2, the
   comparison child, the lead's own recomputation of (1.1)); the first-moment identity (2.1) was checked to order 6 in
   t in exact arithmetic by R2 and on 1,039 marked trees by a DeepSeek child earlier.
