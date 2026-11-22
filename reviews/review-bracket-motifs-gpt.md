VERDICT: CORRECT

# Cold referee report

The two pinned certificates certify the finite-n inequalities in Theorem 1 for every n ≥ 3. The expectation exponent has liminf at least 113/200 and limsup at most 1133/2000. The latter is below 17/30 by exactly 1/6000. The same upper exponent holds almost surely for the limsup of log C_n / log n, as in Corollary 2.

## Inputs and reproducibility

- Proof commit 07919c1db87eda4d746febd4b4eaef13bd1eaaf8, retrieved with `git -C . show 07919c1:attempts/strogatz-exact/PROOF.md`.
- `nice -n 10 sha256sum build/upper_5665.bin build/lower_5650.bin` returned `3a4aef888a5b567a8ada05f9546f98c979d361da9f9b33ba6e794eb91fb30abe` for upper and `9e966fba9aee6eae3da659b693d3bb9c481c1cea66635b5218d0887f182a60d9` for lower. Each file is 138446792 bytes.
- The first read of the permitted header emitted its indexing helpers along with the header struct. The verifier derives state order and table indexing from Section 5 and uses the struct for the binary header layout. The analytic and drift code comes from the pinned proof. The supplied exact expectation table was used for the n = 3 through 21 cross-check.
- Scratch sources and exact logs are under [ref-gpt](../build/ref-gpt). The primary verifier is [replay_gmp.cpp](../build/ref-gpt/replay_gmp.cpp). The checked fixed-width version [replay.cpp](../build/ref-gpt/replay.cpp) has the same row formulas and corroborates the GMP arithmetic. The independently parsed and implemented tree evaluator is [trees.py](../build/ref-gpt/trees.py).

Run from build/ref-gpt. The builds and both certificate runs each exited 0. The fixed-width build had undefined-behavior instrumentation. Both `cmp` calls exited 0 with empty output. The direct tree run exited 0. All verifier and tree arithmetic used integer and rational operations.

```sh
nice -n 10 taskset -c 0 g++ -std=c++20 -O3 -Wall -Wextra -Werror replay_gmp.cpp -lgmpxx -lgmp -o replay_gmp
nice -n 10 taskset -c 0 g++ -std=c++20 -O3 -Wall -Wextra -Werror -fsanitize=undefined -fno-sanitize-recover=all replay.cpp -o replay
nice -n 10 taskset -c 0 ./replay ../upper_5665.bin upper > upper.out
nice -n 10 taskset -c 1 ./replay ../lower_5650.bin lower > lower.out
nice -n 10 taskset -c 0 ./replay_gmp ../upper_5665.bin upper > gmp-upper.out
nice -n 10 taskset -c 1 ./replay_gmp ../lower_5650.bin lower > gmp-lower.out
cmp upper.out gmp-upper.out
cmp lower.out gmp-lower.out
nice -n 10 taskset -c 2 python3 -u trees.py > trees.out
```

## Analytic audit

### Rooting, generator and the two middle domains

T₃ has exactly one degree-two vertex. Its degree persists. The path between two core vertices has degree at least two at every internal vertex, so the induced core is a connected subtree containing the root. For n ≥ 3 a leaf's only neighbor is core. Every non-root vertex of the core has a core parent, every leaf is a child of its core neighbor, and the three disjoint contributions to a core degree are L, k and h. The parent-credit summand is omitted at the root.

Enumerating ordered choices (v,w) gives the unscaled attachment probability at w as the sum of 1/(n deg(v)) over its neighbors. A core target u has scaled rate L_u + ρ_u because each leaf neighbor contributes 1. Each leaf child of u has scaled rate 1/D_u, hence the total budding rate L_u/D_u. The attachment outcomes of budding at different leaf children of the same u are isomorphic as rooted trees and have equal jumps of this potential. Leaf gain increments L and D, while budding decrements L, increments k and leaves D fixed. The newborn state is (1,0,1). Summing the weights of all ordered choices gives exactly n, so the conditional-expectation multiplier in Section 2 is correct.

D and k are nondecreasing. Existing H membership can end on leaf gain crossing D_hu or budding crossing K_hu. Existing C membership can end on leaf gain crossing D_cu or budding crossing K_cu. After either exit the exceeded threshold is permanent. The sole entrance to H is k = 0 to 1 for a non-root vertex of degree below D_hu. There is then no existing child endpoint. The sole entrance to C is k = 1 to 2 at degree below D_cu, with either root indicator. There is then no existing unordered pair. A successful entrance creates exactly the newborn vertical path or cherry counted in the proof. Budding can also create new terms for a middle that stays in its domain. The newborn vertex starts outside both middle domains. Leaf gain creates no core edge, vertical path or cherry and can deactivate middle terms. Both domains preserve the corresponding zero extension once an upper threshold is crossed.

### Lemma 1

For x between adjacent knots, xy − A(x) is the same convex combination of the two knot residuals as x is of the two knot abscissae. This proves both inequalities for arbitrary knot values. The first interval is exactly affine because its right endpoint is exactly 1/D₀. Beyond the first interval the verifier evaluates the appropriate segment at each inner reciprocal degree. The maximum of the affine residuals in y is convex and the minimum is concave. Thus an upper conjugate evaluated at an affine tail jump is at most the chord of its endpoint upper conjugates, with the lower version reversed. Every multiplicity on those terms is nonnegative.

### Lemma 2

For a tail gain, t changes by −t/(D+1), so multiplication by L + ρ yields t γ times the INF-minus-T₁ difference, with γ = (L+ρ)/(D+1). Since each core neighbor has degree at least two, 0 ≤ ρ ≤ (h+k)/2 and γ ≤ (D−(h+k)/2)/(D+1) < 1. For fixed k, (D−k−h)/(D+1) = 1−(k+h+1)/(D+1) is increasing with D. Comparing with its value at D₀ and with zero gives the stated lower bound whenever κ determines k. For the saturated κ class the universal lower bound is zero. Similarly, b = 1−(k+h)/D is increasing in D for fixed k and lies between max(0, D₀−k−h)/D₀ and 1 for the unsaturated class. At k = K_T the new interpolation index is exactly K_T+1. At k > K_T a bud preserves both κ and D, so the change of each existing term through that endpoint is zero. Newborn terms from that bud still enter the vertex row.

### Lemma 3

For each fixed second endpoint, the first-coordinate jump is a convex combination of two numbers. Multiplication by t₁/D₀ preserves the order. Recombining the two second-coordinate bounds gives exactly the two-sided bi-affine estimate in the proof. In a vertical path or cherry with two tail endpoints, the two credits each require this argument, once along each axis. Every occurrence of the other endpoint retains its interpolation weights. The validity of the sum follows term by term, even if the extremizing endpoints differ. An interior point can beat the literal endpoints for a quadratic, for example t(1−t)/D₀ with t = 1/2. The proof uses t times the extremum of the affine jump before multiplication, which bounds that example too.

### Lemma 4 and Section 4.5

- **Vertex and incoming/outgoing credits.** At an inner u, splitting each credit into a center and a residual leaves exactly h r_par + k r_ch times the f gain jump in the vertex sum. The parent residual occurs once on the incoming edge when h = 1. Each child residual occurs once on its outgoing edge. At a tail vertex the f gain jump retains the full interval [L, L+(h+k)/2]. Taking the endpoint extremum of rate times the signed jump works in either direction.
- **Outgoing g gain at inner u.** The rate consists of L, its selected child's credit, the parent credit when h = 1, and each of the k−1 other child credits. The selected child's credit is direct for an inner child and uses Lemma 3 for a tail child. For the parent credit, active H contributes one retained x_par on this edge and one residual to the selected vertical path. Outside H, the outgoing edge takes A*_pc and the incoming edge takes A_pc for each child, hence the k factor on the latter. The other child credits give k−1 copies of x_cc on this edge when C is active, with the two crossed residuals on each unordered pair. Outside C they give k−1 copies each of A_cc and A*_cc on each edge after reindexing the ordered pairs. Root and k = 1 cases remove the parent and other-child credits, respectively.
- **Incoming g gain at inner u.** The selected parent credit is direct for an inner parent and quadratic-bounded for a tail parent. There are k child credits. Active H leaves k x_ch on the incoming edge and one residual on each vertical path. Outside H the incoming edge takes k A*_cp and each outgoing edge takes one A_cp. The factors and subscripts match Section 4.5.
- **Tail and zero edges.** Leaf gain at a tail endpoint follows Lemma 2 and moves toward INF. A tail bud either follows the adjacent κ endpoint or has zero jump on old edges. A tail neighbor's credit to an inner gain produces the quadratic term with the difference at INF and at T₁, including when the inner vertex crosses D₀−1 to D₀. The latter successor is T₁(min(k,K_T+1),h) by Section 3. For a tail–tail edge all old and new g values are zero, no f credit is transferred, and neither endpoint can be a middle with an active newborn term. Such an edge contributes zero to the edge rows.
- **Budding and transfers.** The vertex row gets the change of its old f plus the newborn f and its new edge g, evaluated at the post-bud state. Existing edges get the g differences from each of their endpoints. When u' belongs to H, the incoming edge gets exactly one new path p→u'→nb. When u' belongs to C each previously existing child edge gets exactly one pair with nb. There is no double assignment even when both domains are active. Changes of all existing vertical paths and cherries are accounted for by their endpoint or middle jumps. When a middle exits, its post-transition value is zero. On entering there are zero existing paths or pairs at the old middle, respectively, so the edge newborn terms complete the accounting.
- **h_V rates.** An inner top endpoint sees 1/D_u and at most (h_p+k_p−1)/2 from other core neighbors. An inner bottom endpoint sees the same fixed credit and at most k_w/2 extra. The middle has two specified neighbors with at most (k_u−1)/2 external credit. The source is the sum of the two residuals transferred from g, with the same signs. The middle's leaf gain and budding use the existing path only. Tail endpoint gains and buds use Lemma 2. Tail credit to the middle uses Lemma 3 in the respective endpoint coordinate.
- **h_C rates.** Each inner child sees 1/D_u with at most k_i/2 extra. The apex has two specified children with at most (h_u+k_u−2)/2 external credit. The two ordered source terms occur once each inside the single unordered-pair row. Equal endpoint state IDs can represent distinct children. The row then still contains both endpoint transitions, each with its own rate. The apex gain and bud operate on the existing pair, and the separately assigned newborn pairs complete the bud. Both tail child coordinates use the same conditional interpolation argument as in h_V.
- **Summing finite rows.** The transferred sources in h_V and h_C are affine in each endpoint because the inverse-degree factor uses the opposite endpoint from the g gain jump. Fixed-rate jumps have the same property. Independently taking extrema of signed rate products and the conjugate chord gives a bound in the desired direction. The only extra tail products involve an endpoint's own t in its credit to a neighboring gain, for which Lemma 3 gives the row bound. Consequently each local drift term is bounded by the same nonnegative corner weights applied to its row, and each right-hand side averages to λ times the corresponding term of Φ. Recombining sums yields the stated direction for the entire drift.
- **f saturation.** For L < L_F, the successor after leaf gain is still represented exactly in the L coordinate and the f projection handles a bud at k = K_F. For exact L = L_F, the gain jump is zero, but the bud uses L_F−1. For L ≥ L_F+1 the post-bud projected L is L_F. In this region the budding rate is bracketed by the rate at L_F+1 and 1. The newborn g varies affinely with t over the interval ending at its value for post-bud L_F, so the product is bounded by its four corners. When k > K_F neither f jump changes its overflow value. The child-count clamp is saturated and the largest possible t is achieved at (0,K_F+1,h). Here the rate belongs to [0,1], and the tested rate and newborn-g rectangle contains every case.

### Lemma 5

An edge is assigned to its child, a vertical path to its bottom vertex, and each unordered cherry to its apex. The first two assignments give at most one of each term to a non-root core vertex. An active apex receives exactly k(k−1)/2 terms. Extrema taken with zero allow absent terms and give the indicated bound regardless of the signs of present terms. Multilinear corner weights are nonnegative and sum to one. For the bottom vertex's own tail coordinate, separately comparing INF and the applicable T₁ against each of the two extrema remains valid even when their optimizing corners differ.

For L > L_F the f representative has degree at least D₀ and the same κ and h. For k > K_F the representative uses the same overflow f and saturated κ, and every such physical vertex is tail. Neither region permits an active cherry apex. At D ≥ D₀ the comparison bounds the full t interval, including t = 1 and the formal limit t = 0. The reconstructed extrema cover all projected classes and the strict signs of the constants were checked.

### Theorem 1 and Corollary 2

The conditional drift relation gives the one-step inequality with factor 1+λ/n. Applying total expectation and inducting from the single-core T₃ gives the claimed products. The pointwise inequalities divide by positive constants. In particular the positive lower initial value propagates through the expectation inequality. The exact product has order n^λ by the gamma-ratio asymptotic, or by summing log(1+λ/j) = λ/j + O(j⁻²). This proves the two liminf and limsup inequalities and the strict refutation of 17/30 for the expectation exponent.

For n_k = 2^k, Markov bounds the probability of C_{n_k} > n_k^(1133/2000+ε) by O(2^(−kε)). The first Borel–Cantelli lemma applies to the summable series. Core counts are nondecreasing, so the same eventual bound between dyadic times holds with the constant 2^(1133/2000+ε). A countable intersection over rational ε decreasing to zero gives Corollary 2. The theorem permits distinct liminf and limsup. The external convergence results mentioned in the proof's remark were outside the audit.

## Non-fatal imprecisions

- Section 4.5 abbreviates the degree-crossing successor as T₁(k_v,h_v). Its κ argument must be min(k_v,K_T+1). Section 3 prescribes this projection, and the direct tree case with inner k = 15 and D = 23 exercises it.
- The generic nesting list should include D₀ > 2 for nb to be inner and D₀ ≤ K_F+1 for the VSTAR and inner-projection claims. The pinned values D₀ = 24 and K_F = 150 satisfy both conditions.
- Lemma 1's real-domain interpolant should have real codomain. Its values at rational arguments are rational. The root parent term uses the convention that it is absent when h = 0. The box formulas use f after projection, and the direction symbols denote non-strict inequalities, as the lower certificate has equality rows.
- Active H and C membership can end. The invariant asserted and used in the proof is permanence after exit, with the stated initial entrances by budding. The deterministic T₃ description concerns its rooted shape, for which the potential is deterministic.

## Exact replay and constants

Each certificate passed 484138 vertex, 160545 edge, 13519800 vertical-path and 3164625 cherry inequalities, for 17329108 per certificate and 34658216 combined. Header fields, knot order and endpoints, state counts, domain nesting, successor membership, byte length, and the 841 dense tail–tail zeros per certificate passed. The zeros include all 405 admissible tail–tail edge pairs. The cherry payload is indexed by unordered pairs. The same packed entry serves both orders.

| Mode | λ | Comparison constant | Φ(T₃) | Exact prefactor |
| --- | --- | --- | --- | --- |
| Upper | 1133/2000 | c₀ = 2989767809/10000000000 | 14606753031/10000000000 | 14606753031/2989767809 |
| Lower | 113/200 | C₀ = 1752787393/1000000000 | 1826088077/1250000000 | 7304352308/8763936965 |

The comparison scan tested 241900 projected classes per certificate. The upper minimum occurs at (L,k,h) = (0,1,1). The lower maximum occurs at (668,0,1). Aggregating table extrema with zero as specified in Lemma 5 gives precisely the constants in the table.

## Direct tree tests and finite-n checks

The tree evaluator enumerated the ordered vertex–neighbor choices with weight 1/deg(v), accumulated these exact weights by destination, and rebuilt Φ from the full augmented tree at every destination. It also checked the resulting per-target rates against Section 2. It used neither the local drift rows nor the piecewise-linear centering tables.

The 90 cases comprised 60 named trees and 30 random friend trees using seed 20939. The named tests cover the two T₄ shapes, strict two-tail endpoints for h_V and for h_C at root and non-root apices, tail–tail edges, the degree crossing with and without clamping, L = 0 chains, both root indicators at the K_T, K_F and L_F boundaries, and entries, exits and persistence for both middle domains. The suite asserted observation of 37 named coverage conditions. Both certificates passed the 180 directional drift inequalities and 180 pointwise Φ-to-core-count checks, over 12764 attachment destinations and 25348 ordered choices per certificate.

The exact product brackets held for all 19 supplied expectation rows, or 38 inequalities. The evaluator independently enumerated the rooted labeled distribution up through n = 8, checked unit total probability and recovered the supplied expectations for those six sizes. For n = 9 through 21 it checked the bracket against the supplied exact rational means.

## False controls

Each control changes a loaded copy in memory after header validation. Apply the command form used above with the control name as a final argument. Both verifier builds returned exit 2 on each control and their logs agreed byte-for-byte.

| Mode and control | Change | Rejected rows or structural reason |
| --- | --- | --- |
| Upper `stricter` | λ = 1132/2000 | 1236682 |
| Lower `stricter` | λ = 114/200 | 2519701 |
| Upper `f` | Subtract 1 from f(2,0,0) | 1 vertex row |
| Upper `g` | Subtract 1 from g[(1,1,0),nb] | 1 edge row |
| Upper `hV` | Subtract 10 from h_V[(1,1,0),(0,1,1),nb] | 1 vertical row |
| Upper `hC` | Subtract 10 from h_C[(0,2,0),nb,nb] | 1 cherry row |
| Upper `tail-tail` | Set g[INF,INF] to 1/S | Tail–tail zero check at IDs (572,572) |

## Exact command outputs

Slack is RHS − LHS for upper and LHS − RHS for lower. The row for f gives (L,k,h,case), and the other rows give certificate endpoint or middle IDs in table order, with unused trailing zeros. The tree outputs give Φ(T), direct scaled drift D and the directional slack in actual rational units.

### GMP certificate replay

```text
mode=upper control=none
structure OK states=601 inner=572 parents=555 children=290 middles=84 apices=75 tail_tail_zero=405 dense_tail_tail_zero=841 coefficient_denominator=200783583000000
f checks=484138 failures=0 min_directional_slack=24776800469/2500000000000000 at=0,1,1,0
g checks=160545 failures=0 min_directional_slack=2507610629/13125000000000000 at=276,210,0,0
hV checks=13519800 failures=0 min_directional_slack=2621009681/57500000000000000 at=257,56,252,0
hC checks=3164625 failures=0 min_directional_slack=5253400669/115000000000000000 at=11,252,253,0
total=17329108 failures=0
comparison states=241900 extremizer=0,1,1 constant=2989767809/10000000000 Phi(T3)=14606753031/10000000000 prefactor=14606753031/2989767809
mode=lower control=none
structure OK states=601 inner=572 parents=555 children=290 middles=84 apices=75 tail_tail_zero=405 dense_tail_tail_zero=841 coefficient_denominator=200783583000000
f checks=484138 failures=0 min_directional_slack=0/1 at=0,151,0,0
g checks=160545 failures=0 min_directional_slack=32318417927/165000000000000000 at=279,232,0,0
hV checks=13519800 failures=0 min_directional_slack=18244173547/402500000000000000 at=505,94,252,0
hC checks=3164625 failures=0 min_directional_slack=2515030267/55000000000000000 at=47,192,230,0
total=17329108 failures=0
comparison states=241900 extremizer=668,0,1 constant=1752787393/1000000000 Phi(T3)=1826088077/1250000000 prefactor=7304352308/8763936965
```

### GMP false controls

```text

stricter exit=2
mode=upper control=stricter
structure OK states=601 inner=572 parents=555 children=290 middles=84 apices=75 tail_tail_zero=405 dense_tail_tail_zero=841 coefficient_denominator=200783583000000
f checks=484138 failures=270238 min_directional_slack=-350144045539/475000000000000 at=799,150,1,1
g checks=160545 failures=37346 min_directional_slack=-938950940297/10000000000000000 at=1,1,0,0
hV checks=13519800 failures=831761 min_directional_slack=-539543568339/10000000000000000 at=1,1,572,0
hC checks=3164625 failures=97337 min_directional_slack=-6009145187/375000000000000 at=4,1,572,0
total=17329108 failures=1236682

f exit=2
mode=upper control=f
structure OK states=601 inner=572 parents=555 children=290 middles=84 apices=75 tail_tail_zero=405 dense_tail_tail_zero=841 coefficient_denominator=200783583000000
f checks=484138 failures=1 min_directional_slack=-71327707045877/20000000000000 at=2,0,0,0
g checks=160545 failures=0 min_directional_slack=2507610629/13125000000000000 at=276,210,0,0
hV checks=13519800 failures=0 min_directional_slack=2621009681/57500000000000000 at=257,56,252,0
hC checks=3164625 failures=0 min_directional_slack=5253400669/115000000000000000 at=11,252,253,0
total=17329108 failures=1

g exit=2
mode=upper control=g
structure OK states=601 inner=572 parents=555 children=290 middles=84 apices=75 tail_tail_zero=405 dense_tail_tail_zero=841 coefficient_denominator=200783583000000
f checks=484138 failures=0 min_directional_slack=24776800469/2500000000000000 at=0,1,1,0
g checks=160545 failures=1 min_directional_slack=-22832498996869457/5000000000000000 at=276,0,0,0
hV checks=13519800 failures=0 min_directional_slack=2621009681/57500000000000000 at=257,56,252,0
hC checks=3164625 failures=0 min_directional_slack=5253400669/115000000000000000 at=11,252,253,0
total=17329108 failures=1

hV exit=2
mode=upper control=hV
structure OK states=601 inner=572 parents=555 children=290 middles=84 apices=75 tail_tail_zero=405 dense_tail_tail_zero=841 coefficient_denominator=200783583000000
f checks=484138 failures=0 min_directional_slack=24776800469/2500000000000000 at=0,1,1,0
g checks=160545 failures=0 min_directional_slack=2507610629/13125000000000000 at=276,210,0,0
hV checks=13519800 failures=1 min_directional_slack=-278324999748724947/5000000000000000 at=276,1,0,0
hC checks=3164625 failures=0 min_directional_slack=5253400669/115000000000000000 at=11,252,253,0
total=17329108 failures=1

hC exit=2
mode=upper control=hC
structure OK states=601 inner=572 parents=555 children=290 middles=84 apices=75 tail_tail_zero=405 dense_tail_tail_zero=841 coefficient_denominator=200783583000000
f checks=484138 failures=0 min_directional_slack=24776800469/2500000000000000 at=0,1,1,0
g checks=160545 failures=0 min_directional_slack=2507610629/13125000000000000 at=276,210,0,0
hV checks=13519800 failures=0 min_directional_slack=2621009681/57500000000000000 at=257,56,252,0
hC checks=3164625 failures=1 min_directional_slack=-5566499994959219/100000000000000 at=277,0,0,0
total=17329108 failures=1

tail-tail exit=2
mode=upper control=tail-tail
rejected: nonzero tail-tail g at 572,572

lower-stricter exit=2
mode=lower control=stricter
structure OK states=601 inner=572 parents=555 children=290 middles=84 apices=75 tail_tail_zero=405 dense_tail_tail_zero=841 coefficient_denominator=200783583000000
f checks=484138 failures=385453 min_directional_slack=-95447973453/11500000000000 at=368,0,0,0
g checks=160545 failures=11603 min_directional_slack=-9311383127637/10000000000000000 at=1,1,0,0
hV checks=13519800 failures=1989664 min_directional_slack=-2062650252617/5000000000000000 at=1,1,572,0
hC checks=3164625 failures=132981 min_directional_slack=-15730152467/375000000000000 at=4,1,572,0
total=17329108 failures=2519701
```

### Direct tree drift and finite-n arithmetic

```text
T3 n=3 C=1 ordered=4 targets=3 up: Phi=14606753031/10000000000 D=1654715823/2000000000 slack=2292954123/20000000000000 | lo: Phi=1826088077/1250000000 D=1650907733/2000000000 slack=3878481/62500000000 | 
T4-star n=4 C=1 ordered=6 targets=4 up: Phi=3882153933/2500000000 D=4397906807/5000000000 slack=573599089/5000000000000 | lo: Phi=7767119667/5000000000 D=1755494601/2000000000 slack=62778129/1000000000000 | 
T4-path n=4 C=2 ordered=6 targets=4 up: Phi=2629575843/1250000000 D=23831857331/20000000000 slack=639187/4882812500 | lo: Phi=4202434769/2000000000 D=5936306411/5000000000 slack=29383983/400000000000 | C-entry-h0 H-entry
root-star-L23 n=24 C=1 ordered=46 targets=24 up: Phi=16621183459/10000000000 D=9414774201/10000000000 slack=2252457047/20000000000000 | lo: Phi=16624593951/10000000000 D=9393520923/10000000000 slack=125068137/2000000000000 | D0-1->D0
child-star-L23 n=26 C=2 ordered=50 targets=26 up: Phi=12969123381/5000000000 D=2203255340251/1500000000000 slack=16943566999/30000000000000 | lo: Phi=25905079973/10000000000 D=878982496939/600000000000 slack=8002858543/6000000000000 | C-entry-h0 D0->D0+1 H-entry-blocked-by-D
root-star-L24 n=25 C=1 ordered=48 targets=25 up: Phi=8313186751/5000000000 D=2354428619/2500000000 slack=1126112883/10000000000000 | lo: Phi=103936663/62500000 D=4698249831/5000000000 slack=1563317/25000000000 | D0->D0+1
child-star-L24 n=27 C=2 ordered=52 targets=27 up: Phi=648758185311/250000000000 D=238798266402221/162500000000000 slack=3628655356879/6500000000000000 | lo: Phi=32397294093/12500000000 D=238160175514879/162500000000000 slack=101025200897/81250000000000 | C-entry-h0 H-entry-blocked-by-D strict-tail
root-star-L25 n=26 C=1 ordered=50 targets=26 up: Phi=3326226033/2000000000 D=58877557887/62500000000 slack=450390621/4000000000000 | lo: Phi=1663470729/1000000000 D=234980873707/250000000000 slack=62532943/1000000000000 | strict-tail
child-star-L25 n=28 C=2 ordered=54 targets=28 up: Phi=337498902063/130000000000 D=77643076197/52812500000 slack=1864451877927/3380000000000000 | lo: Phi=168542039589/65000000000 D=495569611263/338000000000 slack=196549475259/169000000000000 | C-entry-h0 H-entry-blocked-by-D strict-tail
root-star-L800 n=801 C=1 ordered=1600 targets=801 up: Phi=16739294053/10000000000 D=59260469241/62500000000 slack=2270004929/20000000000000 | lo: Phi=1045204199/625000000 D=945151403329/1000000000000 slack=286807433/1000000000000 | LF-bud-h0 strict-tail
child-star-L800 n=803 C=2 ordered=1604 targets=803 up: Phi=466884337481/178000000000 D=509251934331269791/343042668000000000 slack=474736635644889719/343042668000000000000 | lo: Phi=1164548435649/445000000000 D=2536650085559315611/1715213340000000000 slack=28156893220562339/85760667000000000000 | C-entry-h0 H-entry-blocked-by-D LF-bud-h1 strict-tail
root-star-L801 n=802 C=1 ordered=1602 targets=802 up: Phi=16739294053/10000000000 D=101264416507/106800000000 slack=599757917083/5340000000000000 | lo: Phi=1045204199/625000000 D=420592126991/445000000000 slack=398068179/1390625000000 | LF-overflow-bud-h0 strict-tail
child-star-L801 n=804 C=2 ordered=1606 targets=804 up: Phi=10518011123177/4010000000000 D=76674330609543659/51649281200000000 slack=3563146266796930623/2582464060000000000000 | lo: Phi=5247009807451/2005000000000 D=238702486413740351/161404003750000000 slack=42511540755468911/129123203000000000000 | C-entry-h0 H-entry-blocked-by-D LF-overflow-bud-h1 strict-tail
root-star-L802 n=803 C=1 ordered=1604 targets=803 up: Phi=16739294053/10000000000 D=3802155261401/4010000000000 slack=903162179649/8020000000000000 | lo: Phi=1045204199/625000000 D=1895026889061/2005000000000 slack=7167178619/25062500000000 | LF-overflow-bud-h0 strict-tail
child-star-L802 n=805 C=2 ordered=1608 targets=805 up: Phi=21062249306063/8030000000000 D=1282690560805720733/864044060000000000 slack=1191134073535626579/864044060000000000000 | lo: Phi=10507103308963/4015000000000 D=1277846915318442581/864044060000000000 slack=14274606738554031/43202203000000000000 | C-entry-h0 H-entry-blocked-by-D LF-overflow-bud-h1 strict-tail
tail-clamp-k12-h0 n=37 C=13 ordered=72 targets=37 up: Phi=25735794259/2000000000 D=1805233874313/250000000000 slack=274912906439/4000000000000 | lo: Phi=127172607507/10000000000 D=1835147968543/250000000000 slack=310679100053/2000000000000 | D0->D0+1 H-entry KT->KT+1-h0
tail-clamp-k12-h1 n=38 C=14 ordered=74 targets=38 up: Phi=138015056447/10000000000 D=3095704614037/400000000000 slack=1585828252601/20000000000000 | lo: Phi=5450337301/400000000 D=126234941393/16000000000 slack=955411997/5000000000 | C-entry-h0 D0->D0+1 H-entry KT->KT+1-h1
tail-clamp-k13-h0 n=38 C=14 ordered=74 targets=38 up: Phi=138470791231/10000000000 D=3070146671741/400000000000 slack=3380072877673/20000000000000 | lo: Phi=136262387681/10000000000 D=9467773863413/1200000000000 slack=572959946603/3000000000000 | D0->D0+1 H-entry KT-overflow-bud-h0
tail-clamp-k13-h1 n=39 C=15 ordered=76 targets=39 up: Phi=36951240461/2500000000 D=49176808537879/6000000000000 slack=5310489964483/30000000000000 | lo: Phi=3632926477/250000000 D=50663277564109/6000000000000 slack=1400794535989/6000000000000 | C-entry-h0 D0->D0+1 H-entry KT-overflow-bud-h1
tail-clamp-k14-h0 n=39 C=15 ordered=76 targets=39 up: Phi=36949189427/2500000000 D=24586949685319/3000000000000 slack=1327773217889/7500000000000 | lo: Phi=1453221021/100000000 D=791640011737/93750000000 slack=350192034817/1500000000000 | D0->D0+1 H-entry KT-overflow-bud-h0
tail-clamp-k14-h1 n=40 C=16 ordered=78 targets=40 up: Phi=31424707759/2000000000 D=6384307006/732421875 slack=2212095686537/12000000000000 | lo: Phi=38585080003/2500000000 D=4502252168471/500000000000 slack=35534532033/125000000000 | C-entry-h0 D0->D0+1 H-entry KT-overflow-bud-h1
box-k150-h0-L0 n=301 C=151 ordered=600 targets=301 up: Phi=1436406330237/10000000000 D=3063394042701409/37750000000000 slack=673180779823951/3020000000000000 | lo: Phi=355572374091/2500000000 D=3256000158143087/37750000000000 slack=444868897171441/75500000000000 | H-entry strict-tail
box-k150-h0-L1 n=302 C=151 ordered=602 targets=302 up: Phi=216904832784499/1510000000000 D=4625886405948889/57002500000000 slack=101638259679330417/456020000000000000 | lo: Phi=53720139968539/377500000000 D=4750647787046459/57002500000000 slack=33502809126110843/11400500000000000 | H-entry KF->KF+1-h0 KT-overflow-bud-h0 strict-tail
box-k150-h1-L0 n=303 C=152 ordered=604 targets=303 up: Phi=109165044178691/755000000000 D=353883445833078741/4332190000000000 slack=965935978158113707/4332190000000000000 | lo: Phi=108097456850049/755000000000 D=188029210828529427/2166095000000000 slack=2560970947290549747/433219000000000000 | C-entry-h0 H-entry strict-tail
box-k150-h1-L1 n=304 C=152 ordered=606 targets=304 up: Phi=27472932527159/190000000000 D=2359200735206467/28880000000000 slack=1609634710535043/7220000000000000 | lo: Phi=27217736802701/190000000000 D=1211319505603303/14440000000000 slack=2129494364766103/722000000000000 | C-entry-h0 H-entry KF->KF+1-h1 KT-overflow-bud-h1 strict-tail
box-k151-h0-L0 n=303 C=152 ordered=604 targets=303 up: Phi=1452387532881/10000000000 D=244933984198763/3020000000000 slack=3544832089117123/3020000000000000 | lo: Phi=715879763957/5000000000 D=4647426499481991/57380000000000 slack=285253138532021/2869000000000000 | H-entry KF-overflow-h0 strict-tail
box-k151-h0-L1 n=304 C=152 ordered=606 targets=304 up: Phi=90772473159/625000000 D=2342418262498849/28880000000000 slack=4214689701100411/3610000000000000 | lo: Phi=27203207934883/190000000000 D=2339228527924837/28880000000000 slack=18856440481781/180500000000000 | H-entry KF-overflow-h0 KT-overflow-bud-h0 strict-tail
box-k151-h1-L0 n=305 C=153 ordered=608 targets=305 up: Phi=27775621454123/190000000000 D=120247958092753793/1472880000000000 slack=432226431999748621/368220000000000000 | lo: Phi=6845816027317/47500000000 D=120081911453889197/1472880000000000 slack=3662516174711729/36822000000000000 | C-entry-h0 H-entry KF-overflow-h1 strict-tail
box-k151-h1-L1 n=306 C=153 ordered=610 targets=306 up: Phi=18638549684773/127500000000 D=849445482528181/10404000000000 slack=15184463275012009/13005000000000000 | lo: Phi=18375462710023/127500000000 D=4241346411064871/52020000000000 slack=67934339379613/650250000000000 | C-entry-h0 H-entry KF-overflow-h1 KT-overflow-bud-h1 strict-tail
box-k152-h0-L0 n=305 C=153 ordered=608 targets=305 up: Phi=1461866690079/10000000000 D=2472188562263/30281250000 slack=22749136255542283/19380000000000000 | lo: Phi=180153694933/1250000000 D=13166922792529/161500000000 slack=8031684904933/80750000000000 | H-entry KF-overflow-h0 strict-tail
box-k152-h0-L1 n=306 C=153 ordered=610 targets=306 up: Phi=74553784435621/510000000000 D=637080369283781/7803000000000 slack=182216592454844729/156060000000000000 | lo: Phi=1531294012737/10625000000 D=3181021051237909/39015000000000 slack=50949861158873/487687500000000 | H-entry KF-overflow-h0 KT-overflow-bud-h0 strict-tail
box-k152-h1-L0 n=307 C=154 ordered=612 targets=307 up: Phi=12506508752433/85000000000 D=141072226561764493/1716660000000000 slack=1007582648059521361/858330000000000000 | lo: Phi=7398541737043/51000000000 D=986128778285198431/12016620000000000 slack=11956899788279431/120166200000000000 | C-entry-h0 H-entry KF-overflow-h1 strict-tail
box-k152-h1-L1 n=308 C=154 ordered=614 targets=308 up: Phi=8092295045903/55000000000 D=38067399701151/463203125000 slack=69232222950037361/59290000000000000 | lo: Phi=55851290447243/385000000000 D=4865815349677449/59290000000000 slack=619456786283557/5929000000000000 | C-entry-h0 H-entry KF-overflow-h1 KT-overflow-bud-h1 strict-tail
kF+LF n=1103 C=152 ordered=2204 targets=1103 up: Phi=459271273070047/3170000000000 D=3454594005918057409/42205380000000000 slack=9404917931276752907/42205380000000000000 | lo: Phi=455186097688633/3170000000000 D=1712230832871130421/21102690000000000 slack=36521262831107647/4220538000000000000 | C-entry-h0 H-entry KF->KF+1-h1 KT-overflow-bud-h1 LF-bud-h1 strict-tail
both-overflows n=1106 C=153 ordered=2210 targets=1106 up: Phi=347440565912719/2382500000000 D=47588389271272266371/577620924000000000 slack=163141383896763013679/722026155000000000000 | lo: Phi=343014317356469/2382500000000 D=237255259502124379201/2888104620000000000 slack=29050672443391923803/36101307750000000000 | C-entry-h0 H-entry KF-overflow-h1 KT-overflow-bud-h1 LF-overflow-bud-h1 strict-tail
tail-tail-edge n=51 C=2 ordered=100 targets=51 up: Phi=32785820167/10000000000 D=784096456748559/422500000000000 slack=4958834128187/3380000000000000 | lo: Phi=16388429551/5000000000 D=978152546464179/528125000000000 slack=974393327257/4225000000000000 | C-entry-blocked-by-D-h0 H-entry-blocked-by-D strict-tail tail-tail-g
tail-inner-tail-V-D15 n=65 C=3 ordered=128 targets=65 up: Phi=15754222405079/3250000000000 D=156544317721321027/57037500000000000 slack=3413719866190877/2281500000000000000 | lo: Phi=7872493614073/1625000000000 D=156145132000216661/57037500000000000 slack=43549785455923/114075000000000000 | C-entry-blocked-by-D-h0 C-entry-blocked-by-D-h1 H-bud-persist H-entry-blocked-by-D H-exit-D strict-tail two-strict-tail-hV
tail-inner-tail-V-D2 n=52 C=3 ordered=102 targets=52 up: Phi=10800262150541/3250000000000 D=91652454968279/48750000000000 slack=49109062377259/19500000000000000 | lo: Phi=5392064301927/1625000000000 D=2381172423352373/1267500000000000 slack=48896854931441/12675000000000000 | C-entry-blocked-by-D-h0 H-entry-blocked-by-D L0-chain strict-tail two-strict-tail-hV
nonroot-tail-inner-tail-V n=55 C=4 ordered=108 targets=55 up: Phi=534799824187/105625000000 D=283363978219193/98865000000000 slack=52604939254657/24716250000000000 | lo: Phi=8536882522811/1690000000000 D=847693496578247/296595000000000 slack=479027129046107/118638000000000000 | C-entry-blocked-by-D-h1 C-entry-h0 C-entry-h1 H-bud-persist H-entry-blocked-by-D strict-tail two-strict-tail-hV
tail-tail-C-root-D11 n=61 C=3 ordered=120 targets=61 up: Phi=390317175017/81250000000 D=1896167131665503/697125000000000 slack=99681970687669/69712500000000000 | lo: Phi=600115459151/125000000000 D=14548397783638663/5362500000000000 slack=4998454014299/10725000000000000 | C-bud-persist C-exit-D-h0 H-entry-blocked-by-D strict-tail two-strict-tail-hC-h0
tail-tail-C-nonroot n=63 C=5 ordered=124 targets=63 up: Phi=21469315368859/3250000000000 D=86918255664688093/23237500000000000 slack=1706780159642601/929500000000000000 | lo: Phi=21428713534463/3250000000000 D=259755821428677383/69712500000000000 slack=223539704546681/278850000000000000 | C-bud-persist C-entry-h0 C-exit-D-h1 H-bud-persist H-entry H-entry-blocked-by-D strict-tail two-strict-tail-hC-h1
mixed-C-with-crossing n=51 C=3 ordered=100 targets=51 up: Phi=536706179201/130000000000 D=188750441306743/80730000000000 slack=9371394544861/12420000000000000 | lo: Phi=535133488917/130000000000 D=2441919745992319/1049490000000000 slack=208959067419467/209898000000000000 | C-bud-persist D0-1->D0 H-entry-blocked-by-D strict-tail
L0-chain n=8 C=6 ordered=14 targets=8 up: Phi=21691431939/5000000000 D=12286847523/5000000000 slack=2697340887/10000000000000 | lo: Phi=21360938791/5000000000 D=50181283633/20000000000 slack=95278098267/1000000000000 | C-entry-h0 H-entry L0-chain
H-entry-degree-2 n=4 C=2 ordered=6 targets=4 up: Phi=2629575843/1250000000 D=23831857331/20000000000 slack=639187/4882812500 | lo: Phi=4202434769/2000000000 D=5936306411/5000000000 slack=29383983/400000000000 | C-entry-h0 H-entry
H-entry-degree-15 n=17 C=2 ordered=32 targets=17 up: Phi=1287533457/500000000 D=218790322087/150000000000 slack=519778603/3000000000000 | lo: Phi=1028738087/400000000 D=435956848487/300000000000 slack=116336483/1200000000000 | C-entry-h0 H-entry
H-entry-degree-16 n=18 C=2 ordered=34 targets=18 up: Phi=25782451877/10000000000 D=116832184027/80000000000 slack=3471969891/20000000000000 | lo: Phi=643750241/250000000 D=29099453437/20000000000 slack=9712719/100000000000 | C-entry-h0 H-entry-blocked-by-D
H-entry-degree-23 n=25 C=2 ordered=48 targets=25 up: Phi=5184979579/2000000000 D=67549635233/46000000000 slack=16112383161/92000000000000 | lo: Phi=12945880683/5000000000 D=336485966031/230000000000 slack=2252707983/23000000000000 | C-entry-h0 D0-1->D0 H-entry-blocked-by-D
C-entry-h0-degree-2 n=4 C=2 ordered=6 targets=4 up: Phi=2629575843/1250000000 D=23831857331/20000000000 slack=639187/4882812500 | lo: Phi=4202434769/2000000000 D=5936306411/5000000000 slack=29383983/400000000000 | C-entry-h0 H-entry
C-entry-h0-degree-11 n=13 C=2 ordered=24 targets=13 up: Phi=6388468487/2500000000 D=318440549891/220000000000 slack=9345280731/55000000000000 | lo: Phi=12767945257/5000000000 D=158715999899/110000000000 slack=1044035449/11000000000000 | C-entry-h0 H-entry
C-entry-h0-degree-12 n=14 C=2 ordered=26 targets=14 up: Phi=25613608051/10000000000 D=43525211047/30000000000 slack=10231671349/60000000000000 | lo: Phi=2559528091/1000000000 D=10846714603/7500000000 slack=57145391/600000000000 | C-entry-blocked-by-D-h0 H-entry
C-entry-h1-degree-3 n=6 C=3 ordered=10 targets=6 up: Phi=6096125763/2000000000 D=34530927639/20000000000 slack=724961679/4000000000000 | lo: Phi=1900938453/625000000 D=51556538647/30000000000 slack=77195041/750000000000 | C-entry-h0 C-entry-h1 H-bud-persist H-entry
C-entry-h1-degree-11 n=14 C=3 ordered=26 targets=14 up: Phi=34620471157/10000000000 D=431424738949/220000000000 slack=50193080691/220000000000000 | lo: Phi=17281764237/5000000000 D=85930570937/44000000000 slack=1409787659/11000000000000 | C-entry-h0 C-entry-h1 H-bud-persist H-entry
C-entry-h1-degree-12 n=15 C=3 ordered=28 targets=15 up: Phi=2170007227/625000000 D=14749988137/7500000000 slack=860496073/3750000000000 | lo: Phi=3466168559/1000000000 D=39170281707/20000000000 slack=25769903/200000000000 | C-entry-blocked-by-D-h1 C-entry-h0 H-bud-persist H-entry
H-exit-k n=20 C=10 ordered=38 targets=20 up: Phi=9756321473/1000000000 D=260868622711/50000000000 slack=619167320469/2000000000000 | lo: Phi=48075460689/5000000000 D=2184279311/400000000 slack=28171219643/1000000000000 | C-entry-h0 H-entry H-exit-k PL-CC
H-continues n=6 C=3 ordered=10 targets=6 up: Phi=6096125763/2000000000 D=34530927639/20000000000 slack=724961679/4000000000000 | lo: Phi=1900938453/625000000 D=51556538647/30000000000 slack=77195041/750000000000 | C-entry-h0 C-entry-h1 H-bud-persist H-entry
C-exit-k-root n=14 C=7 ordered=26 targets=14 up: Phi=68336228993/10000000000 D=9677064047/2500000000 slack=8435073069/20000000000000 | lo: Phi=33975717607/5000000000 D=33802341999/8750000000 slack=167080972063/7000000000000 | C-exit-k-h0 H-entry
C-exit-k-child n=16 C=8 ordered=30 targets=16 up: Phi=77812519751/10000000000 D=176303740287/40000000000 slack=9714734383/20000000000000 | lo: Phi=7719277371/1000000000 D=71633253601/16000000000 slack=46274654179/400000000000 | C-entry-h0 C-exit-k-h1 H-bud-persist H-entry
C-continues n=8 C=4 ordered=14 targets=8 up: Phi=39903520497/10000000000 D=90411840177/40000000000 slack=4768634601/20000000000000 | lo: Phi=7958388941/2000000000 D=11253443829/5000000000 slack=977555987/400000000000 | C-bud-persist C-entry-h0 H-bud-persist H-entry
inner-k-clamped-crossing n=40 C=17 ordered=78 targets=40 up: Phi=83136667127/5000000000 D=4261923108467/460000000000 slack=35496854428993/230000000000000 | lo: Phi=40796821751/2500000000 D=2120945589577/230000000000 slack=16339748001/11500000000000 | C-entry-h0 D0-1->D0 H-entry PL-CC PL-PC-CP inner-crossing-clamped
both-PL-external n=30 C=12 ordered=58 targets=30 up: Phi=118253643521/10000000000 D=238413168763/36000000000 slack=13766559168637/180000000000000 | lo: Phi=117225702789/10000000000 D=1192269522161/180000000000 slack=8412479687/18000000000000 | C-entry-h0 H-entry PL-CC PL-PC-CP
random-00 n=31 C=5 ordered=60 targets=31 up: Phi=61989891877/10000000000 D=1264042877261/360000000000 slack=89488839269/180000000000000 | lo: Phi=30959366207/5000000000 D=279890683563/80000000000 slack=450326293/2000000000000 | C-bud-persist C-entry-blocked-by-D-h1 C-entry-h1 H-bud-persist H-entry PL-PC-CP
random-01 n=37 C=1 ordered=72 targets=37 up: Phi=16665721061/10000000000 D=14160008807/15000000000 slack=6750658339/60000000000000 | lo: Phi=3333830283/2000000000 D=9418695749/10000000000 slack=25007981/400000000000 | strict-tail
random-02 n=13 C=1 ordered=24 targets=13 up: Phi=8251180219/5000000000 D=1869491581/2000000000 slack=1129283127/10000000000000 | lo: Phi=4126319651/2500000000 D=4663054231/5000000000 slack=31302537/500000000000 | 
random-03 n=38 C=9 ordered=74 targets=38 up: Phi=11279023051/1250000000 D=6795581142551/1330000000000 slack=729418894389/332500000000000 | lo: Phi=281232337/31250000 D=6764486401509/1330000000000 slack=9305665361/6650000000000 | C-bud-persist H-bud-persist H-entry L0-chain PL-CC
random-04 n=30 C=5 ordered=58 targets=30 up: Phi=51513971799/10000000000 D=23281066330969/7980000000000 slack=6700358289533/7980000000000000 | lo: Phi=51405765021/10000000000 D=23183689346357/7980000000000 slack=637207133873/798000000000000 | C-entry-blocked-by-D-h0 C-entry-h1 H-bud-persist H-entry L0-chain
random-05 n=28 C=5 ordered=54 targets=28 up: Phi=411936477/78125000 D=29841576481/10000000000 slack=3595167403/1250000000000 | lo: Phi=10507777893/2000000000 D=1785703327181/600000000000 slack=1853989727/240000000000 | C-entry-h1 H-bud-persist H-entry L0-chain PL-CC
random-06 n=22 C=3 ordered=42 targets=22 up: Phi=43672287723/10000000000 D=371062700681/150000000000 slack=17025698077/60000000000000 | lo: Phi=21827204737/5000000000 D=37000644613/15000000000 slack=706516757/3000000000000 | H-entry PL-CC
random-07 n=63 C=5 ordered=124 targets=63 up: Phi=1593031769219/250000000000 D=1610372887037041/446250000000000 slack=201928230653939/178500000000000000 | lo: Phi=1590617662229/250000000000 D=4585348449642833/1275000000000000 slack=3967311859939/2550000000000000 | C-bud-persist H-bud-persist H-entry strict-tail
random-08 n=66 C=14 ordered=130 targets=66 up: Phi=175957388571/12500000000 D=740353929795941/93000000000000 slack=15802915719587/1162500000000000 | lo: Phi=350238710393/25000000000 D=98359132856723/12400000000000 slack=5205916404717/310000000000000 | C-bud-persist C-entry-blocked-by-D-h1 C-entry-h1 H-bud-persist H-entry L0-chain strict-tail
random-09 n=20 C=1 ordered=38 targets=20 up: Phi=518583189/312500000 D=9399749153/10000000000 slack=140862149/1250000000000 | lo: Phi=8298890953/5000000000 D=586148269/625000000 slack=62552711/1000000000000 | 
random-10 n=31 C=1 ordered=60 targets=31 up: Phi=16650056231/10000000000 D=23577828439/25000000000 slack=2250958523/20000000000000 | lo: Phi=16653740837/10000000000 D=2940621507/3125000000 slack=125049899/2000000000000 | strict-tail
random-11 n=70 C=8 ordered=138 targets=70 up: Phi=388848615083/47500000000 D=1940299515369277/418760000000000 slack=3345946073451/817890625000000 | lo: Phi=1549290566823/190000000000 D=199924238015221/43320000000000 slack=17231359854107/2166000000000000 | H-entry L0-chain strict-tail
random-12 n=11 C=3 ordered=20 targets=11 up: Phi=9978991821/2500000000 D=18087901509/8000000000 slack=314822517/1250000000000 | lo: Phi=39903665157/10000000000 D=1352819352883/600000000000 slack=851040607/6000000000000 | C-bud-persist H-entry
random-13 n=54 C=9 ordered=106 targets=54 up: Phi=1190885760399/130000000000 D=99875505224251/19266000000000 slack=1056660557751647/192660000000000000 | lo: Phi=1185964065313/130000000000 D=994403574661531/192660000000000 slack=136028385299671/19266000000000000 | C-entry-blocked-by-D-h1 H-entry L0-chain PL-CC strict-tail
random-14 n=31 C=3 ordered=60 targets=31 up: Phi=499082604913/130000000000 D=19841431617887/9126000000000 slack=62351390746579/91260000000000000 | lo: Phi=19168706861/5000000000 D=845120646493/390000000000 slack=35573512873/39000000000000 | C-entry-blocked-by-D-h0 C-entry-h1 H-bud-persist H-entry strict-tail
random-15 n=50 C=5 ordered=98 targets=50 up: Phi=63023922319/11000000000 D=667220223424571/205700000000000 slack=8536977157429/4114000000000000 | lo: Phi=125868468629/22000000000 D=83171859188971/25712500000000 slack=1772883447673/822800000000000 | C-entry-h1 H-bud-persist H-entry L0-chain strict-tail
random-16 n=37 C=2 ordered=72 targets=37 up: Phi=910842801293/350000000000 D=2256710936443/1531250000000 slack=2419260437183/4900000000000000 | lo: Phi=455146528381/175000000000 D=1125613258391/765625000000 slack=175338735749/245000000000000 | C-entry-blocked-by-D-h0 H-entry strict-tail
random-17 n=12 C=1 ordered=22 targets=12 up: Phi=16478810219/10000000000 D=4667058069/5000000000 slack=2259702127/20000000000000 | lo: Phi=16481761161/10000000000 D=4656410613/5000000000 slack=125234007/2000000000000 | 
random-18 n=51 C=10 ordered=100 targets=51 up: Phi=752546673/62500000 D=73135233957011/10725000000000 slack=52204226653/26812500000000 | lo: Phi=6009502591/500000000 D=388539332012891/57200000000000 slack=21824508203/11440000000000 | C-bud-persist C-entry-h1 C-exit-D-h0 H-bud-persist H-entry L0-chain PL-CC
random-19 n=29 C=3 ordered=56 targets=29 up: Phi=22089855577/5000000000 D=4404387582459/1760000000000 slack=63292304927/220000000000000 | lo: Phi=5519196641/1250000000 D=6586372173671/2640000000000 slack=10630147463/66000000000000 | C-entry-blocked-by-D-h0 C-entry-h1 H-bud-persist H-entry
random-20 n=27 C=8 ordered=52 targets=27 up: Phi=92350677043/10000000000 D=5752641204181/1100000000000 slack=438247150709/220000000000000 | lo: Phi=92170884529/10000000000 D=11458648606323/2200000000000 slack=18076593683/22000000000000 | C-bud-persist C-exit-D-h1 H-bud-persist H-entry
random-21 n=47 C=8 ordered=92 targets=47 up: Phi=794595322183/90000000000 D=1517649283289/303750000000 slack=25076968276153/4860000000000000 | lo: Phi=396439066061/45000000000 D=16148631128041/3240000000000 slack=134312004247/20250000000000 | C-bud-persist C-entry-h1 H-bud-persist H-entry L0-chain strict-tail
random-22 n=19 C=1 ordered=36 targets=19 up: Phi=8293044763/5000000000 D=9394892581/10000000000 slack=1127135479/10000000000000 | lo: Phi=16589150409/10000000000 D=4686747769/5000000000 slack=125111383/2000000000000 | 
random-23 n=40 C=6 ordered=78 targets=40 up: Phi=37677140431/5000000000 D=36268385716877/8505000000000 slack=75856950503423/17010000000000000 | lo: Phi=56437073759/7500000000 D=48305028806659/11340000000000 slack=287391986939/35437500000000 | C-entry-h1 H-bud-persist H-entry strict-tail
random-24 n=9 C=1 ordered=16 targets=9 up: Phi=818441077/500000000 D=9271804949/10000000000 slack=113245341/1000000000000 | lo: Phi=16372176293/10000000000 D=9250906311/10000000000 slack=125341091/2000000000000 | 
random-25 n=32 C=4 ordered=62 targets=32 up: Phi=11458786151/2500000000 D=10852228842947/4180000000000 slack=348973461597/1045000000000000 | lo: Phi=11435365607/2500000000 D=270103335439/104500000000 slack=34305898481/104500000000000 | C-entry-blocked-by-D-h1 C-entry-h0 H-entry L0-chain PL-PC-CP
random-26 n=15 C=3 ordered=28 targets=15 up: Phi=30096913439/10000000000 D=204577133953/120000000000 slack=10841802661/60000000000000 | lo: Phi=15014391389/5000000000 D=13573830929/8000000000 slack=1603737/15625000000 | C-entry-blocked-by-D-h0 H-entry L0-chain
random-27 n=27 C=1 ordered=52 targets=27 up: Phi=3327102331/2000000000 D=122497589167/130000000000 slack=5854566499/52000000000000 | lo: Phi=3327829523/2000000000 D=61111333931/65000000000 slack=325145193/5200000000000 | strict-tail
random-28 n=65 C=12 ordered=128 targets=65 up: Phi=130638352977/10000000000 D=375304295218837/50820000000000 slack=797382999356081/50820000000000000 | lo: Phi=715714140693/55000000000 D=62446279210279/8470000000000 slack=17199182858107/847000000000000 | C-bud-persist C-entry-h1 H-bud-persist H-entry L0-chain PL-CC strict-tail
random-29 n=10 C=2 ordered=18 targets=10 up: Phi=202205493/80000000 D=57268021251/40000000000 slack=5347713/32000000000 | lo: Phi=1262990421/500000000 D=114181806709/80000000000 slack=37363253/400000000000 | C-entry-h0 H-entry
coverage 37 required checks exercised: C-bud-persist C-entry-blocked-by-D-h0 C-entry-blocked-by-D-h1 C-entry-h0 C-entry-h1 C-exit-D-h0 C-exit-D-h1 C-exit-k-h0 C-exit-k-h1 D0-1->D0 D0->D0+1 H-bud-persist H-entry H-entry-blocked-by-D H-exit-D H-exit-k KF->KF+1-h0 KF->KF+1-h1 KF-overflow-h0 KF-overflow-h1 KT->KT+1-h0 KT->KT+1-h1 KT-overflow-bud-h0 KT-overflow-bud-h1 L0-chain LF-bud-h0 LF-bud-h1 LF-overflow-bud-h0 LF-overflow-bud-h1 PL-CC PL-PC-CP inner-crossing-clamped strict-tail tail-tail-g two-strict-tail-hC-h0 two-strict-tail-hC-h1 two-strict-tail-hV
trees=90 directional_drifts=180 targets_per_certificate=12764 ordered_choices_per_certificate=25348 pointwise_comparisons=180 failures=0
exact-enumeration n=3 states=1 E=1
exact-enumeration n=4 states=3 E=4/3
exact-enumeration n=5 states=12 E=19/12
exact-enumeration n=6 states=60 E=647/360
exact-enumeration n=7 states=360 E=17189/8640
exact-enumeration n=8 states=2520 E=4717423/2177280
finite-n n=3 lower=7304352308/8763936965 E=1 upper=14606753031/2989767809 PASS
finite-n n=4 lower=1302000798901/1314590544750 E=4/3 upper=4961427112863/854219374000 PASS
finite-n n=5 lower=1188726729396613/1051672435800000 E=19/12 upper=45312713821777779/6833754992000000 PASS
finite-n n=6 lower=63002516658020489/50079639800000000 E=647/360 upper=504466442977852013607/68337549920000000000 PASS
finite-n n=7 lower=82722304371980902057/60095567760000000000 E=17189/8640 upper=2208385931876043498233577/273350199680000000000000 PASS
finite-n n=8 lower=125158846514807104812241/84133794864000000000000 E=4717423/2177280 upper=33419504307080166258768720741/3826902795520000000000000000 PASS
finite-n n=9 lower=71465701359954856847789611/44871357260800000000000000 E=101568829/43545600 upper=572576367293204488511484492455553/61230444728320000000000000000000 PASS
finite-n n=10 lower=136713886701593641149821525843/80768443069440000000000000000 E=200684675729/80621568000 upper=1217233737268986830965581421572455061/122460889456640000000000000000000000 PASS
finite-n n=11 lower=288876442600467363749572884106259/161536886138880000000000000000000 E=446709364753753/169305292800000 upper=3674828652815071242685090311727241829159/349888255590400000000000000000000000000 PASS
finite-n n=12 lower=6749204158938192043967293746846233/3589708580864000000000000000000000 E=2175595493439775699/782190452736000000 upper=7728164656870094823366744925562389566721377/699776511180800000000000000000000000000000 PASS
finite-n n=13 lower=2422964293058810943784258455117797647/1230757227724800000000000000000000000 E=863016923161025341999/295667991134208000000 upper=64743987440372031065225466738053178993469456047/5598212089446400000000000000000000000000000000 PASS
finite-n n=14 lower=505654009774504160806668706825737308947/246151445544960000000000000000000000000 E=1182295484917290649261473233/387443335582266163200000000 upper=135130662401508793760981737615661300433062057763327/11196424178892800000000000000000000000000000000000 PASS
finite-n n=15 lower=490990043491043540143275314327790926987537/229741349175296000000000000000000000000000 E=418365039667608049635216883810549/131564133463670121037824000000000 upper=3936761587743155688638680961957060665516396928819005491/313499877008998400000000000000000000000000000000000000 PASS
finite-n n=16 lower=1528452005387618540466016053502413155712202681/689224047525888000000000000000000000000000000 E=2259642509981131170779006369702211737/683804583677425454094090240000000000 upper=40854399503735888684796018129536389899840661861640699317101/3134998770089984000000000000000000000000000000000000000000 PASS
finite-n n=17 lower=5063761493849180224563911185253494784874527482153/2205516952082841600000000000000000000000000000000 E=17287848151143794747769099286795627457612027287/5046595988971459310332853429993472000000000000 upper=1353628818757281199793346468685929206551420649461741290473507433/100319960642879488000000000000000000000000000000000000000000000 PASS
finite-n n=18 lower=348803806429258237821431764584226023122827745976539/147034463472189440000000000000000000000000000000000 E=2191102148433572559496221254607828432463286999047319/618321052399156726204326001091032173772800000000000 upper=399639002431929078927223878019686981628328249391087031581560812131/28662845897965568000000000000000000000000000000000000000000000000 PASS
finite-n n=19 lower=1295108533271835837030976141901231223855059420810889307/529324068499881984000000000000000000000000000000000000 E=5732108499635476087835196297054462943247004863583847762453/1566690289049250828667904890451936865052262400000000000000 upper=1648866119700535831978289362500559632089412542737692749302010848540047/114651383591862272000000000000000000000000000000000000000000000000000 PASS
finite-n n=20 lower=723965670098956232900315663322788254134978216233287122613/287347351471364505600000000000000000000000000000000000000 E=2588997551427461141383863803537365961038824299890823463539916684413/686520175275134242800419706889084121927323666612224000000000000000 upper=3396056729591635195410863032775494741187104264997585808338715291364087329/229302767183724544000000000000000000000000000000000000000000000000000000 PASS
finite-n n=21 lower=330852311235222998435444258138514232139685044818612215034141/127709933987273113600000000000000000000000000000000000000000 E=8060798970277915903320915972897680147282202467385602270877224958802337893/2076891346250126117286176382601163467171069659509371699200000000000000000 upper=139690001458292730492835029127154425189249159732145697054396376079679004103757/9172110687348981760000000000000000000000000000000000000000000000000000000000 PASS
finite-n rows=19 exact_inequalities=38 failures=0 independently_enumerated_E_rows=6
```
