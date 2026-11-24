# Stochastic domination between adjacent labels: a coupling lemma and exact certificates

**Authors:** Fable lead (Claude Fable 5.2 for the original transport; the lead continued as Claude Opus 5.5), with Fable's GPT-6.2 Sol Ultrafast proof author for the coupling, strictness and endpoint tail. **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** Astra (GPT) and Neon (Gemini) audited the base argument at distinct ranges. GPT referee and a separate implementation checked parts of the local extension; the full computation shares shape-law inputs. See [F5](../../README.md#f5).

*Found by a GPT-6.2 Sol Ultrafast child of the Fable lead (24 Nov 2026, about 12:21 UTC) while exploring the second
half of question 3 of Addario-Berry, Briend, Devroye, Donderwinkel, Kerriou and Lugosi. This note is the lead's own
write-up (Claude family): the lemma was re-derived and its proof written here independently of the child's text, and
the certificates of section 4 were recomputed by a program of the lead that shares nothing with the child's.*

Model: the random friend tree in discrete time. Given the tree T with n vertices, the next leaf is attached at the
vertex w with probability r_T(w)/n, where r_T(w) is the sum of 1/D_x over the neighbours x of w (pick a uniform vertex,
then a uniform neighbour); the sum of r_T over all vertices is n. D_v(n) is the degree of v at size n, Z_v = lim D_v(n)/n.

## 1. The order

For rooted trees write ch(v) for the set of children of a vertex v (all neighbours, for the root). Define two
relations between a vertex alpha of a rooted tree A and a vertex beta of a rooted tree B, by mutual recursion:

    E(alpha, beta):  there is an injection phi of ch(beta) into ch(alpha) with O(phi(b), b) for every b in ch(beta);
    O(alpha, beta):  there is an injection psi of ch(alpha) into ch(beta) with E(a, psi(a)) for every a in ch(alpha).

(E for a pair at even depth, O for a pair at odd depth; both hold when the set to be injected is empty.) Write
(A, a) >= (B, b) if E(a, b) holds for the two roots. Unfolding the recursion: there is a rooted subtree M common to
both trees, containing the roots (the *matched* vertices), such that at a matched vertex of even depth every child in B
is matched (A may have further children, with anything below them), and at a matched vertex of odd depth every child
in A is matched (B may have further children). The relation is reflexive, and it is transitive (compose the two
matchings: at even depth the children of the third tree are matched into the second and those into the first; at odd
depth the other way round). It does not require the trees to have the same size.

Examples. (i) If u is a vertex of a tree T with a leaf neighbour v, then (T, u) >= (T, v): match u with v (roots) and
the child v of u with the child u of v; at the odd pair the vertex of the first tree is a leaf, so nothing is to be
injected. (ii) More generally, a rooted tree whose root has a leaf child is >= every rooted tree whose root is a leaf:
match the leaf child with the unique child of the other root. (iii) Deleting a whole child branch at a vertex of even
depth makes a rooted tree smaller in this order; deleting one at a vertex of odd depth makes it larger.

## 2. The coupling lemma

**Lemma.** Let (A, a) >= (B, b) with the same number n_0 of vertices. Then the two friend-tree chains (A_n), (B_n)
started from A and B can be defined on one probability space so that (A_n, a) >= (B_n, b) for every n >= n_0. In
particular D_a(A_n) >= D_b(B_n) for all n, and Z_a >= Z_b.

*Proof.* Fix a matching M as in section 1 and call a matched pair even or odd by its depth.

*Degrees.* At an odd pair (alpha, beta) all children of alpha are matched to distinct children of beta, and both have a
parent, so D_A(alpha) <= D_B(beta). At an even pair all children of beta are matched to distinct children of alpha (and
both or neither have a parent), so D_A(alpha) >= D_B(beta).

*Rates.* At an even pair (alpha, beta): every neighbour of beta in B is matched (its parent, if any, and all its
children), each such neighbour is the B-member of an odd pair, and its partner in A is a neighbour of alpha with a
degree that is not larger. Hence r_B(beta), the sum of the reciprocal degrees of the neighbours of beta, is at most
the sum of the reciprocal degrees of their partners, which is at most r_A(alpha). At an odd pair (alpha, beta), in the
same way with the roles exchanged (every neighbour of alpha in A is matched, to a neighbour of beta of degree not
larger): r_A(alpha) <= r_B(beta).

*One step.* Both trees have n vertices. For every matched pair (alpha, beta) attach, with probability
min(r_A(alpha), r_B(beta))/n, a new leaf at alpha in A and a new leaf at beta in B, and match the two new leaves. With
the remaining probability rho attach in A a leaf at a vertex chosen with probabilities proportional to
(r_A(alpha) - r_B(beta)) at the even pairs and to r_A at the unmatched vertices of A, and in B a leaf at a vertex
chosen with probabilities proportional to (r_B(beta) - r_A(alpha)) at the odd pairs and to r_B at the unmatched
vertices of B, in any joint way; these new leaves stay unmatched. Because the rates sum to n in each tree, the two
residual masses are both equal to rho, and each tree receives its leaf at w with probability r(w)/n: both marginals
are friend-tree steps.

*Invariance.* After a common step at an even pair, the two new leaves form an odd pair whose A-member has no children;
after a common step at an odd pair they form an even pair whose B-member has no children; the injection at the pair
where they were attached is extended by one element. After a residual step, A has gained a child at an even matched
vertex or somewhere in an unmatched part, B at an odd matched vertex or in an unmatched part: no condition is
affected, because A never gains an unmatched child at an odd pair (there all its attachment probability is in the
common part) and B never gains one at an even pair. So the matching, enlarged by the pairs of new leaves, satisfies the
conditions of section 1 after the step, and the argument repeats. The root pair is even, hence the degree inequality. []

**Corollary 1 (inside one tree).** In every finite tree, a vertex with at least one leaf neighbour stochastically
dominates every leaf of the tree: for such u and any leaf v, the processes can be coupled with D_u(n) >= D_v(n) for all
n; in particular P(Z_u > x) >= P(Z_v > x) for all x. (Examples (i), (ii): the two rooted trees are the same tree.)

## 3. The certificate for two labels

**Proposition.** Let u < w be labels and m >= w. Let q_u, q_w be the laws of the rooted shapes (T_m, u), (T_m, w). If
there is a coupling of q_u and q_w carried by pairs (s, t) with s >= t, then the degree processes (D_u(n))_{n >= m}
and (D_w(n))_{n >= m} can be coupled with D_u(n) >= D_w(n) for all n >= m. In particular Z_u stochastically dominates
Z_w, and D_u(n) stochastically dominates D_w(n) for every n >= m.

*Proof.* Given the rooted shape of (T_m, u), the future of the tree seen from u is a friend-tree chain from that shape
(Markov property; the transition probabilities depend on the shape only). Draw (s, t) from the coupling and run the
coupled chains of the lemma from s and t. []

Such a coupling exists if and only if a transportation problem on the finite set of rooted shapes of size m is
feasible, i.e. the maximal flow from "supply q_u(s)" to "demand q_w(t)" along the pairs s >= t equals 1. All data are
rational, so this is decided exactly. Two remarks used by the programs: common mass may be cancelled first (the
relation is reflexive and transitive, so a feasible plan can be rerouted to keep the diagonal), and a flow may pass
through shapes of smaller sizes along the moves of example (iii) (by transitivity every such path joins comparable
shapes).

## 4. Results

**Z_1 = Z_2 in law, and Z_u stochastically dominates Z_{u+1} for every u = 2, ..., 19** (each pair certified
by two independent programs with identical exact integers). Equivalently P(Z_u > x) >= P(Z_{u+1} > x) for all x in
[0, 1) and all these u, with the stronger process-level statement of the proposition for all n >= u + 1; by
transitivity Z_u dominates Z_w for all u < w <= 20.

The first two pairs need no computation: at size 3 the label 2 is the centre or an end of the path with probability
1/2 each and the label 3 is an end; at size 4 the labels 3 and 4 are either two leaves of the same vertex, or two ends
of the path, or an interior vertex and the end next to it; examples (i) and reflexivity give the couplings.

*Reduction used by the sparse programs.* At the birth size N = u + 1 the later label is a leaf, so every shape with
q_u < q_{u+1} is leaf-rooted. The older label is a leaf or has degree 2 with the later label as a leaf child; in the
second case its shape is >= every leaf-rooted shape (example (ii)). After cancelling common mass it is therefore
enough to transport ALL the positive mass that sits on leaf-rooted shapes into the negative mass; the mass on the
degree-2 shapes ("universal") fills what remains. For u <= 5 there is no positive leaf-rooted mass at all, so those
pairs follow from example (ii) alone.

Three programs:

* `lead-domination_transport.py` (lead; dense): all histories enumerated with exact probabilities, comparability of
  every pair of shapes by recursive bipartite matching, own exact max-flow between the two full laws. Sizes 4 to 10,
  every pair (u, u+1) with u + 1 <= size, 42 instances: max-flow = 1 exactly in each (logs `lead-transport_*.log`).
* `lead-transport_sparse.cpp` (lead; sparse): exact signed measure q_u - q_{u+1} at size u + 1 from the lead's
  shape-propagation code, flow network on leaf-rooted shapes of sizes at most u + 1 whose arcs are the branch moves of
  example (iii), exact integer max-flow (logs `lead-transport_sparse_*.log`).
* `child-transport.cpp` (GPT-6.2 Sol Ultrafast child; sparse, written first): its own network on leaf-rooted shapes
  and its own max-flow (`child-transport-all.txt`). For the exact laws it reuses the lead's shape-propagation code, so
  the two sparse programs are independent in the network and in the flow, not in the laws. The laws are checked
  separately: against brute force over all histories (dense program, sizes up to 10) and, as integrated against
  classifier bounds, against the other child's independent program on 150 intervals of the hub-order certificate
  (`../q3-hub-probability/EXTENSION.md`, sizes up to 23).

| pair | tree size | leaf-rooted shapes with positive / negative mass | network (lead): shapes, arcs | mass to transport / universal mass (same integers in both sparse programs) | result |
|---|---:|---:|---:|---|---|
| (2,3) | 3 | 0 / 1 | 2, 1 | no positive leaf-rooted mass | flow complete: certified |
| (3,4) | 4 | 0 / 1 | 2, 1 | no positive leaf-rooted mass | flow complete: certified |
| (4,5) | 5 | 0 / 3 | 6, 7 | no positive leaf-rooted mass | flow complete: certified |
| (5,6) | 6 | 0 / 5 | 12, 17 | no positive leaf-rooted mass | flow complete: certified |
| (6,7) | 7 | 2 / 13 | 30, 56 | ratio 0.02095 | flow complete: certified |
| (7,8) | 8 | 5 / 29 | 67, 156 | ratio 0.02832 | flow complete: certified |
| (8,9) | 9 | 15 / 71 | 162, 450 | ratio 0.03757 | flow complete: certified |
| (9,10) | 10 | 40 / 170 | 389, 1,271 | ratio 0.04315 | flow complete: certified |
| (10,11) | 11 | 111 / 424 | 971, 3,627 | ratio 0.04778 | flow complete: certified |
| (11,12) | 12 | 297 / 1,063 | 2,441, 10,328 | ratio 0.05090 | flow complete: certified |
| (12,13) | 13 | 796 / 2,732 | 6,264, 29,587 | ratio 0.05340 | flow complete: certified |
| (13,14) | 14 | 2,117 / 7,092 | 16,228, 84,738 | ratio 0.05531 | flow complete: certified |
| (14,15) | 15 | 5,672 / 18,649 | 42,564, 243,123 | ratio 0.05693 | flow complete: certified |
| (15,16) | 16 | 15,270 / 49,387 | 112,547, 698,281 | ratio 0.05824 | flow complete: certified |
| (16,17) | 17 | 41,347 / 131,866 | 300,084, 2,008,280 | ratio 0.05947 | flow complete: certified |
| (17,18) | 18 | 112,397 / 354,240 | 805,246, 5,782,789 | ratio 0.06052 | flow complete: certified |
| (18,19) | 19 | 306,848 / 957,339 | 2,173,996, 16,671,457 | ratio 0.06146 | flow complete: certified |
| (19,20) | 20 | 840,822 / 2,599,950 | 5,899,410, 48,112,650 | ratio 0.06227 | flow complete: certified |

One program only (the lead's sparse program; 4 minutes, 10.8 GB): the pair (20,21) at size 21 also passes (2,312,218
positive and 7,093,568 negative leaf-rooted shapes; network of 16,084,239 shapes and 138,992,643 arcs; ratio of the
mass to transport to the universal mass 0.06302). So Z_20 dominates Z_21 as well, with one certificate.

All three programs agree wherever they overlap (dense and sparse of the lead on the pairs up to (9,10); the two sparse
programs on every pair of the table, integer for integer).

### 4.1 A scalar certificate (lead)

The full order is not needed. Two special cases of it suffice for a certificate at the birth size u + 1:

* a root with a leaf child is >= every leaf-rooted shape (example (ii));
* a leaf whose parent has k other neighbours is >= a leaf whose parent has at least k other LEAF neighbours (match
  every other neighbour of the first parent with a leaf of the second; at those pairs nothing is required).

Let d = q_u - q_{u+1} on rooted shapes of size u + 1, P_k the total of d^+ over leaf-rooted shapes whose root's parent
has exactly k other neighbours, and N_j the total of d^- over (leaf-rooted) shapes whose root's parent has exactly j
other leaf neighbours. If for every K the sum of P_k over k >= K is at most the sum of N_j over j >= K, the positive
leaf-rooted mass can be transported by the monotone coupling of the two indices, the universal mass fills the rest,
and Z_u stochastically dominates Z_{u+1}. This needs the signed table only (no network), so it runs in the sharded
program `../q3-hub-probability/lead-signed_orbit_dp.cpp` (lines "T u size ..." of its output).

Result: the condition holds for every pair computed, and the largest ratio (sum of P_k, k >= K)/(sum of N_j, j >= K)
is attained at K = 2 each time:

    u        6      7      8      9      10     11     12     13     14     15     16     17     18     19
    ratio    0.208  0.250  0.328  0.367  0.398  0.419  0.434  0.446  0.456  0.463  0.470  0.476  0.481  0.485

    and 0.488, 0.492, 0.494 at u = 20, 21, 22 (sharded run of `../q3-hub-probability/lead-signed_orbit_dp.cpp`, 41 minutes,
    8 GB; logs `lead-dp3_scalar_*.log`).

The increments decrease (0.0068, 0.0056, 0.0049, 0.0041 for the last four steps), which suggests a limit near 0.5 and
therefore that the scalar condition, and with it the domination, holds for all labels. This is an extrapolation, not
a proof: a proof for all large u would need the corresponding inequalities between stationary frequencies of local
configurations (they appear to hold with a margin of a factor 2) and a way to join them to the finite range.

### 4.2 A local criterion (lead): the coupling written out, and a certificate that needs only the law of T_(u-1)

Fix u >= 4, let T = T_(u-1) with n = u - 1 vertices, and let a = u, b = u + 1 be attached at x and then at y (y a vertex
of T, or y = a). Write r for the rates in T, eps_x = 1/(D_x (D_x + 1)) and rho(w) = r(w) eps_w = r(w)/(D_w (D_w + 1));
for a leaf of a core p, rho = 1/(2 D_p). The probability of (x, y) given T is w(x, y)/(n (n + 1)) with

    w(x, x) = r(x)(r(x) + 1),   w(x, a) = r(x)/(D_x + 1),   w(x, y) = r(x)(r(y) - eps_x) if y ~ x,   w(x, y) = r(x) r(y) otherwise

(an attachment at x lowers the rate at each neighbour of x by eps_x, raises the rate at x by 1, and the new leaf has
rate 1/(D_x + 1)). Let S_T(x; y) be the tree T with a leaf at x and a leaf at y, rooted at the leaf on x.

*Pieces of the coupling, given T (masses in units of 1/(n (n + 1))).*

1. y = x: the older and the younger label have the same rooted shape. Mass w(x, x), coupled identically.
2. x, y distinct and not adjacent: the older label in the outcome (x, y) and the younger label in the outcome (y, x) both
   have the shape S_T(x; y), and w(x, y) = w(y, x). Coupled identically.
3. x ~ y: the older label in (x, y) and the younger label in (y, x) both have the shape S_T(x; y), with masses
   w(x, y) and w(y, x), and w(x, y) - w(y, x) = rho(y) - rho(x). The common part is coupled identically. There remains,
   for every edge {x, z} with rho(z) > rho(x), the mass delta = rho(z) - rho(x) of S_T(x; z) on the side of the OLDER
   label (its leaf sits at the endpoint with the smaller rho) and the same mass of S_T(z; x) on the side of the
   YOUNGER label.
4. y = a: mass U_x = r(x)/(D_x + 1) = D_x rho(x); the older label has degree 2 and a leaf child (universal shape), the
   younger one is the leaf of a vertex of degree 2 hanging at x (shape K_T(x)).

*The three kinds of edges with rho(z) > rho(x).*

* x is a leaf of T and z = p its parent (this happens when rho(p) > 1/(2 D_p)). The older-side shape S_T(x; p) is
  isomorphic to K_T(p) (exchange the new leaf at p with the leaf x, and the two vertices of degree 2). These masses,
  L_p (rho(p) - 1/(2 D_p))^+ in total at p, are coupled identically with part of the mass U_p = D_p rho(p) of K_T(p),
  which is larger. The younger-side shape S_T(p; x) is a leaf at p, and p has L_p - 1 other leaf neighbours there:
  SUPPLY at level L_p - 1.
* z is a leaf of T and x its parent (when 1/(2 D_x) > rho(x)). Older side: a leaf at x, which has D_x other
  neighbours: DEMAND at level D_x, total L_x (1/(2 D_x) - rho(x))^+ at x. Younger side: a leaf at the former leaf z
  (no other leaf neighbour); left to the universal mass.
* x and z both cores: DEMAND delta at level D_x (older side, leaf at x), SUPPLY delta at level L_z (younger side,
  leaf at z, which has L_z other leaf neighbours).

A leaf whose parent has k other neighbours is >= a leaf whose parent has at least k other leaf neighbours. So if, with

    dem_K(T) = sum over cores x with D_x >= K of [ sum over core neighbours z of (rho(z) - rho(x))^+ + L_x (1/(2 D_x) - rho(x))^+ ],
    sup_K(T) = sum over cores p with L_p - 1 >= K of L_p (rho(p) - 1/(2 D_p))^+ + sum over core-core edges {x, z} with rho(z) > rho(x), L_z >= K, of (rho(z) - rho(x)),

**E dem_K(T_(u-1)) <= E sup_K(T_(u-1)) for every K >= 2,** then all the demand can be assigned to supply of at least
the same level (the levels are nested, so this cumulative condition is Hall's condition), mixing over the trees T;
what remains on the younger side is leaf-rooted and has the same total mass as the universal mass, which is >= all of
it. This is a coupling of the two laws carried by comparable pairs, and Z_u stochastically dominates Z_(u+1) by the
proposition of section 3. (For K <= 1 the condition follows from K = 2, since every demand has level at least 2.)

The two expectations are those of sums of local terms over the vertices and edges of T_(u-1). Program
`lead-local_scalar.cpp` computes them exactly from the law of the unrooted tree T_(u-1) (all free trees of that size).

Result (exact; logs `lead-local_scalar_*.log`): the criterion holds for every pair from (4,5) to (25,26), and the
largest ratio E dem_K/E sup_K is attained at K = 2 each time:

    u        6      8      10     12     14     16     18     20     22     24     25
    ratio    0.208  0.367  0.436  0.476  0.505  0.526  0.542  0.554  0.563  0.571  0.574

A GPT-6.2 Sol Ultrafast child implemented the same criterion separately (its own code for the local terms; the laws
from the lead's propagation code): it also passes for every pair through (25,26). A second GPT child refereed this
section cold (`../../checks/stochastic-domination-local-criterion-referee-gpt-6.2-sol-ultrafast.md`): CORRECT; its own
exact enumeration matches the lead's numbers for n = 3 to 11, it built the coupling explicitly for u <= 8 and checked
the marginals and the comparability of every coupled pair, and dropping the reserve for K_T(p) makes that explicit
construction fail (plan mass 71/72 at u = 4) although the scalar test alone would still pass. Its one remark on the
extrapolation below: uniform convergence of local averages is not enough; control relative to the tail masses is needed.

On the lead's large simulated trees (30 to 300 trees of 10^9 to 6 x 10^9 vertices in total per file; program
`lead-stationary_hall.py`) the same ratios, computed tree by tree without any cancellation, are 0.640 at K = 2, 0.526
at K = 3 and 4, 0.60 at K = 10, 0.65 at K = 20, 0.71 at K = 50, 0.74 at K = 100, 0.80 at K = 1000, and 0.83 to 0.85
for K between 10^4 and 3 x 10^5: all below 1. So in the stationary regime the criterion appears to hold at every
level, with at least 15% to spare. The value at large K has a simple explanation (a remark, not a proof). For a
vertex x of very large degree the supply is L_x (rho(x) - 1/(2 D_x)), which tends to 1/2, and the demand is the sum
of rho(z) over its core children z. Children are promoted at rate about 1 and then grow like inert newborns, so the
mean demand is the integral over the age t of E rho(z_t). The degree of z increases at rate r(z), hence the integral
of rho(z_t) = r(z_t)/(D (D + 1)) over the time z spends at degree D has mean q_D/(D (D + 1)), q_D being the
probability that an inert newborn ever reaches degree D + 1. So the ratio tends to 2 sum_{D >= 2} q_D/(D (D + 1)),
which would be exactly 1 if every child were certain to keep growing, and is below 1 because some children stop
(with q_2 = 0.93, q_3 = 0.85, ... decreasing to the hub probability 0.7 of an inert newborn it is about 0.83, the
value seen in the simulation). A proof of the domination for all sufficiently large labels would follow from
(i) these strict inequalities for the stationary frequencies of local configurations and (ii) the convergence of the
means of local sums divided by E C_n, uniformly in the level K; neither is proved here.

## 5. Strictness, and the tail at 1 for every label

**Strictness.** For every certified pair, P(Z_u > x) > P(Z_{u+1} > x) for every x in [0, 1). Proof (child's; checked
by the lead). The coupling gives positive mass to pairs (s, t) in which the root of s has degree at least 2 and the
root of t is a leaf (the universal mass). In such a pair the root of s has an unmatched child, so its residual mass
is positive; the parent of the root of t is matched to a child of the root of s, and it receives the reciprocal degree
1 from the leaf root of t against at most 1/2 on the other side, so its residual mass is positive too. Choosing the
two residual sites independently, with positive probability the next k steps attach to the root of s and to the parent
of the root of t. Then the root of s has at least k leaves, so E[Z | tree] >= (k - 1)/m with m the size, which tends to 1,
while the root of t is still a leaf whose parent has k other leaves: E[Z | tree] <= 1/m, and the probability that it is
ever promoted is at most 1/(k - 1) (the supermartingale 1/max(k - 1, 1) of the first certificate). For large k the
event {Z of the first > x >= Z of the second} has positive probability, for x > 0 by Markov's inequality on both
sides and for x = 0 by the promotion bound; since the first is never smaller than the second, the inequality between
the tail probabilities is strict.

**A stochastic upper bound and the exact tail exponent (all labels v >= 3; child's, checked by the lead).** The parent
of v has degree at least 2 from the birth of v on, so the attachment rate at v is at most D_v - 1/2. Hence D_v(n) - 1/2
is dominated by the red weight of a Polya urn started at size v with weights (1/2, v - 1/2), and

    Z_v is stochastically smaller than a Beta(1/2, v - 1/2) variable;        P(Z_v > x) <= C_v (1 - x)^(v - 1/2).

Conversely, attach v to a leaf (positive probability) and then every new vertex to v up to size m: at size j this has
probability (j - v + 1/2)/j, the product is of order m^-(v - 1/2), and at size m the label v has m - v leaves, so
E[1 - Z_v | tree] <= (v + 1)/m; with m = 2 (v + 1)/(1 - x) this gives P(Z_v > x) >= c_v (1 - x)^(v - 1/2). So

    P(Z_v > 1 - eps) is of exact order eps^(v - 1/2),

and for EVERY u >= 3 there is x_u < 1 with P(Z_u > x) > P(Z_{u+1} > x) for all x in (x_u, 1). (No value of x_u is
computed.)

### 5.1 Towards all labels: what is known (second exploration child, report `child2-REPORT.md`)

* The child implemented the criterion of section 4.2 in its own way (after the identity S_T(x; p) = K_T(p), sources
  binned by the number k of other neighbours of the root's parent, targets by the number of other leaves; nested Hall
  inequalities on the law of the unrooted tree T_(u-1)) and evaluated it exactly through size 26: every pair from (2,3)
  to (25,26) passes; the smallest Hall slack at (25,26) is 0.00667 times the universal mass, at the cut 21; two of its
  integer implementations with different accumulator widths agree byte for byte. This agrees with the lead's program.
* *No proof by bounded lookback.* For every k >= 2, the coupling cannot be built conditionally on the tree k steps
  before the birth of the older label: from the 5-vertex broom (a star with one arm subdivided) a positive singleton
  shape arises that has no comparable negative shape among the outcomes of that seed, for every k (exact
  computation, and an argument: only the centre can reach the largest parent degree, and the negative shapes at that
  degree have two non-leaf branches where the positive one has one). Exhaustively over seeds of size 10, 60, 61 and 62
  of the 106 seeds fail at lookbacks 2, 3, 4. So any proof must use the mixture over different histories, as the
  unconditional criterion does.
* *The extreme cut holds for every n.* At the cut j = n - 2 (n = u - 1) the only supply is the broom's centre facing
  its arm and the only demand the star's centre after a promotion; with the exact probabilities of the star and the
  broom (p_star = 2/(n - 1), and a recursion for the broom) the inequality holds for all n, with a gap of order 7/(9n)
  in the natural units. A margin that tends to 0 at the extreme cuts therefore does not by itself indicate failure.
* Per tree, without averaging, the inequality fails at moving large cut-offs on every simulated path (256 paths of 10^4
  to 10^6 vertices); averaged over trees (section 4.2, the lead's simulation) it holds at all levels with room. The
  question is entirely about the averaging.

### 5.2 Two structurally different routes, both closed (14:02 UTC; GPT children, reports `route-sdh-REPORT.md`, `route-sdi-REPORT.md`)

* *Coupling whole histories.* (a) For every u >= 4 no co-adapted coupling (one that gives each chain its correct next
  step given the joint past) can keep the order at size u + 1: the younger chain's star states would need more
  broom-tip states than the older chain has. (b) For every u >= 4 no deterministic measure-preserving permutation of
  the finite attachment histories works: one class of equally likely histories has 2(u - 1) younger broom-tip outcomes
  and only two eligible older ones. (c) Comparing the two labels in the same history fails, and swapping their birth
  sites does not preserve probabilities. (d) For every u >= 6 a double star rules out couplings that keep the final
  unrooted tree. (e) An unrestricted coupling of whole histories is equivalent to the endpoint transport of section 3
  (condition each marginal on its endpoint), so this route gives nothing beyond the certificate. Proved on the way: the
  inequality for the principal up-set of the broom tip, for every u >= 4.
* *Induction through the first promotion.* The natural joint order on (time of the first promotion, rooted shape at
  that time) fails for every pair (u, u + 1) with u >= 3: at the deadline u + 2 an explicit up-set has more mass for the
  younger label, by 1/(2 (u - 2)(u - 1) u (u + 1)). Proved on the way: the coupling lemma extends to unequal starting
  sizes (a larger shape at an earlier size suffices), and from any common finite seed the older newborn is strictly
  more likely to be promoted at its first step. The tails of the parent's leaf stock at birth cross for every u >= 3.
* Both children checked their statements with exact enumeration (u up to 8 or 10) and wrote them themselves (GPT
  family); the lead has not re-derived them. They are dead ends for these particular routes, not counterexamples to
  the domination.

### 5.3 The route through the stationary local limit, closed (14:27 UTC; report `route-sdj-REPORT.md`)

Proved by the child: the expected Hall scores E dem_K(T_n)/E C_n and E sup_K(T_n)/E C_n converge to their stationary
values uniformly in the cut K in ABSOLUTE error, and the stationary supply is asymptotic to (A/2) K^(-mu). But at the
extreme cut K = n - 2 the finite-size supply divided by the stationary one tends to zero, so no uniform RELATIVE
transfer from the stationary law exists, and a small explicit threshold needs a signed finite-size estimate at
fixed and moderate cuts plus a separate treatment of near-extreme configurations on their own scale. The proposed
inert occupation coefficient was bounded below 92/93, but neither its identification with the physical limit nor the
stationary gap at K = 2 was certified. No threshold below 10^6 resulted. Astra's claim 24 gives domination for all
labels above 2^(2^(2^123)); the labels 26 to that threshold remain open.

## 6. What is not proved

* The statement for all labels on the whole range of x. The certificate is one finite computation per pair. No
  general construction of the coupling is known; conditionally on the tree two steps before it can fail (child's
  example: from the 5-vertex seed whose root has three children, one of them with a child, the maximal flow at size 7
  is 1210/1220), so a general proof must mix different histories. In the computed range the mass that has to be
  transported is a stable fraction, about 6%, of the universal mass (0.062 at the pair (19,20)).
* Domination between non-adjacent labels follows by transitivity of stochastic domination, for the certified range only.
