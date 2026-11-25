> **Public checkpoint export note:** Edited export of `order-astra/REVIEW-Q3.md` (frozen original SHA-256 `07798e101efd98e19a47d2ffdaf49cb6f184e1e9f0217a65b17aa438bd21fa40`). Changes are limited to local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

HOLDS — Neon's strict pre-placement order for every old tree and its marginal founder-to-later order survive the audit. Two local statements need wording corrections. The finite coverage count is 11,307 accepted box–size pairs.

# Q3 cross-family review by Order-Astra

**Verdict and scope.** Write b=n+1 and m=n+2. For each finite old tree of size n≥2, and each 0≤x<1, the conditional tail of Z_b strictly exceeds that of Z_m. Averaging gives the marginal order for b≥3. The additional founder argument gives the strict marginal comparison for b=2. Swapping the two initial labels gives equality in law for 1 and 2. These statements concern the one-step random-friend rule started from an edge. The comparison conditions on the tree before either of the distinguished placements. At x=1 both tails vanish. For x<0 both equal one.

**Additional founder estimate.** Section 5 gives a separate 14-row proof that the marginal gap between labels 2 and 3 lies strictly above (15/272)(1−x) and at most 1−x. It handles all thresholds including zero. This author-checked corollary awaits independent review.

**Material reviewed.** Neon attributes the following board sources to local commit `4edb12d84426bb5837ea679475ddecd6f5883989`. The board proof exactly matches the owner's SHA-256. These hashes identify the inputs to the review.

| File | SHA-256 |
|---|---|
| PROOF.md | `e85e9cebfbde5439fad9246da2fec2c6eb7f014b09248e960c20d016411ce191` |
| cert_all_m.cpp | `90689e44d08fdd94bb9edc76cae6d14669a5ac39e84d4da35d19a82d8a517c80` |
| cert_small_m.cpp | `a10dce484242acb727973baf48c1bdac448dce9e761c506f55326a66b58c3210` |
| audit_primitives.py | `f397012211e2e14ffaa64193148285b3c90d3be797e8d7109e4e883969c1773a` |
| planted_defects.py | `6679dc7296fbe5fe287074ef8dfdf76fa77dcc7fe446c664aa87e33667114e92` |

The four source files were read from the authorized Neon board export, `audit-4edb12d/research_order/`. **UNSUPPORTED by local Git evidence** is the narrow assertion that Astra independently authenticated their membership in that commit. The named Git object is absent from Astra's checkout, and the board export has no Git tree/commit object or bundle. The producer explicitly attests to the commit. Mathematical and execution verdicts below apply to the exact bytes above.

**Order of work.** My own local commit `e432fa7d679fc016039e313122aff82f7795d63e` and its fresh combined replay both preceded re-opening Neon in this turn. The replay matches the frozen receipt byte for byte. My earlier mathematical route and its component checks preceded the previous cross-review. This renewed review benefits from that prior exposure. Both exact versions of the Astra proof and the witness are pinned. See `Q3-PIN.md` (internal review pin record, not exported).

## 1. Pre-placement identity and the two physical bounds

**HOLDS — §2 and Lemmas 3.1–3.2.** Given the old tree, put w(u)=Σ_{v∼u}1/d(v), s(u)=w(u)/[d(u)(d(u)+1)]. After multiplying by n(n+1), a serial event at u has coefficient w(u)/(d(u)+1)=d(u)s(u). On an old edge uv the two ordered parallel events have masses w(u)w(v)−s(u) and w(u)w(v)−s(v). Swapping the two new labels swaps their rooted future laws. When s(u)>s(v), the unmatched later label is on u and the unmatched earlier one on v. It enters the earlier-minus-later difference with a minus sign times the later-minus-earlier payoff. Twins and nonadjacent parallel hosts cancel. These are identities on each fixed tree.

For a tag of degree d with j tracked leaves, the true gain weight is at least j and the tracked loss weight is j/d≤1. The Poisson mixture of integer beta tails satisfies the birth–death recursion with weights j and 1 by the urn identity and k/k!=1/(k−1)!. Its increase with j gives the lower submartingale. For u>x the binomial tail tends to zero whenever Z<u, since tracked stock divided by size is at most normalized degree. Taking the bounded reverse-Fatou inequality and then u↓x gives P(Z>x), including x=0 even if Z has an atom there. For the upper law, the attachment weight at size s can increase by at most the number of subsequent gains of the tag. Common uniforms therefore dominate that count by a real-initial-stock Pólya urn. The initial degree offset vanishes on division by size. Both comparisons are conditional on the full present tree. The limits exist by bounded nonnegative degree-share supermartingales.

**WRONG — the ancillary strict parameter monotonicities immediately before Lemma 3.1 at x=0.** For example B_{4,1/2}(0)=B_{5,1/2}(0)=B_{4,1}(0)=1. Use non-strict monotonicities on [0,1) and strict parameter monotonicities on (0,1). This substitution suffices everywhere in the proof. Strictness of the final zero-threshold inequality comes from its own ledger.

**HOLDS — Lemma 3.3 and §4.3 as bounds on a physical pending tag, with a repair to the first-jump sentence.** If the tag is a leaf on a host of degree k with j other leaf neighbours, the only state changes before its first child have weights 1/k (promotion), j/k (sibling loss), and g (host growth). Divide by the actual tree size N at each step. The last weight is j+1 plus reciprocal degrees of the other k−j−1 neighbours, hence lies between j+1 and (k+j+1)/2. An unrecorded exterior event may alter g. For a proposed price u the conditional drift numerator is
```
(pay−u)/k + (j/k)(loss−u) + g(grow−u).
```
It is affine in g for each current exterior. Clearing the positive sum of rates gives the same two endpoint tests for the first-relevant-event bound. The displayed mixture in the proof is fractional-linear in g. **WRONG** is its description as affine or any use of a constant g between relevant jumps. The pointwise drift proof works at every physical step. Taking a minimum with an independently valid current stopping price is legitimate by stopping when that price is selected. One may also prove the static recursion by descending k and ascending j, bounding each next relevant payoff at its actual jump time; conditional on still waiting, the promotion weight divided by N is at least 1/(kN), so waiting forever for a relevant event at fixed k has probability zero. A tag may nevertheless remain pending through infinitely many *different* states.

On the branch with a permanent nonleaf neighbour one can use c₁(j)=1/(j+1). At j≥1 its maximum drift is bounded by 1/k−1/(j+2)≤0. The adverse endpoint of every signed edge has such a neighbour, even when the original opposite host was a leaf, since it just received the other mark. In the founder star j=k−1 this invariant fails. At (k,j)=(3,2), g=3, c₁'s exact drift is 1/12. The broader cap c₂(j)=min(1,2/(j+1)) works there. At j≥2 its drift at g=j+1 is at most 1/k−2/(j+2)≤0 (at the threshold j=2 the same formula applies), and at j=0,1 the current cap is one. A transition from the star branch into the invariant branch can only lower this upper price. The finite table is defined through k=Kmax+1. For larger k read the analytic extension c₁ on the invariant branch and c₂ otherwise. This fills the omitted domain convention in the unrestricted wording of Lemma 3.3.

No child ever implies share zero. A probability of ever being promoted is thus an upper bound on the positive-share probability, also allowing extinction after promotion. At y=0 the certificate actually uses this defective upper. The beta cap alone equals one there. At positive y the upper table uses the minimum of the hub cap, the locked beta cap and the continuation bound. At promotion the tag has weight 1+1/k, after its next gain 2+1/k, after its only leaf is lost 1/2+1/k. Its host-growth weight ranges from 1/2 to k−1/2. All are physical transitions with unchanged tag at exterior events. For k above the table, the j=0 top-row cap is valid for larger hosts because this boundary is stopped at a beta cap of shape 1/(Kmax+1) or at one; these dominate the larger hosts' beta price. This deduction requires the actual top-row stopping construction.

**HOLDS — favourable states and stopping boundaries.** The positive-stock lower row V(d,j) tracks a chosen subset of real leaves, so it applies even on a star with more leaves than tracked. Extra real gains are favourable by the checked gain monotonicity. At the last tracked loss the former leaf has degree two and stays adjacent, the input to R(d,2). In R(d,k) the selected neighbour sends with weight 1/k. Its own birth weight is at most k−1+1/d. Its growth is unfavourable and any extra tag gain is favourable by the checked value inequalities. Extra untracked stock can be ignored. At the return cutoff zero is a valid lower. H(p) uses the exact serial middle and both endpoints of its host-growth range. At its cutoff one can forget the host and retain V(2,1). At the degree cutoff the Poisson stock floor is valid at the actual stopping size. Thus none of these rows assumes that a projected state specifies the exterior law.

## 2. The discrete-time bridge and finite-time induction

**HOLDS — Lemma 4.3, including a crossing into a stopping slab.** Fix N and a physical vector of rates. Write Qf(s)=Σr(s,t)(f(t)−f(s)). For an interior step of clock length δ, the change in expected certificate value multiplied by N is
```
−Nδ σ_i(s) + Q f(t+δ)(s),
```
where σ_i is the nonnegative backwards slope. Lower clock log(N/m) has Nδ<1. Upper clock log((N−1)/(m−1)) has Nδ>1. In both cases Qf is affine within a slab, so the two nodal generator comparisons suffice. All m≥200 have δ<1/199<1/128.

For a crossing let θ=t+δ−t_k. At fixed N and δ the exact scaled drift is
```
−N(δ−θ)σ_in(s) − Nθσ_out(s) + Q(f_k−θσ_out)(s).
```
Evaluate this affine expression at θ=0 and θ=δ. At the latter endpoint an outgoing continuing row supplies its nodal comparisons. An outgoing frozen lower row instead needs Q(f_k−θσ_out)(s)≥σ_out(s) throughout the possible overshoot. The producer tests its two endpoints. For a lower crossing N+1≥m exp(t_k), so δ≤1/N≤1/(m exp(t_k)−1). This explains both the minus one and the downward-rounded exponential used in the denominator. A lower current frozen row is already below the right-knot stock floor throughout the slab. The terminal lower floor uses one extra knot and is valid for the whole last overshoot.

For the upper last crossing let F be the stopped values. Its expected payoff is F_s+QF_s/N. The incoming terminal value dominates F_s+u max(0,QF_s), with u≥1/N from the preceding knot. If two adjacent slabs give distinct upper prices for an absorbing payoff, their minimum bounds it at a crossing: the primitive at the earlier knot is valid at all later sizes, including throughout the outgoing slab. Substitute that common price at both endpoints. This can only improve each upper comparison. Finite-time induction then bounds the true continuation at either clock's first stopping size. Every path hits the finite horizon if it has not stopped earlier.

I tested the completed integer arrays against the generator inequalities directly. The new assertions compute Q from stored **left and right values** after all rounding, then compare with their finite-difference slope. They also test both extremal overshoots when the outgoing lower row is frozen. The source under review supplies the arrays. These assertions do not use its candidate-slope calculation. In particular they catch any illicit truncation of a negative return slope to zero. An analytic justification for that truncation is that the right return value starts at zero and stays below the rate-weighted mean of its two successor values. Moving the target values backwards increases them, and the left-end inequality preserves this invariant.

**HOLDS — §5 for 4≤m≤199.** Adjacent degrees d,k have d+k≤N; a promoted serial middle with host k has k≤N−2; a pending leaf can have k=N−1 on a star. The finite sweep uses respectively these three bounds. From k=N−1 the next host size k+1=N is allowed at size N+1. For lower rows, a current value above the stock price must obey every physical lower expectation. For upper rows, a current value below both independent stopping caps must dominate every upper expectation. Floors and ceilings are applied after multiplying out the denominators Nd, Ndk or 2Nk. The finite checker below verifies these implications on the completed rows. At the starting N, the exact-size solver either inherits the full N=200 table or uses raw pointwise stock/beta stops, with zero return lower. The latter construction is valid at any feasible starting size; its cutoff heuristic only affects tightness.

## 3. Arithmetic, geometry and endpoints

**HOLDS — Lemmas 4.1–4.2 and the integer enclosure directions on the certified domain.** For the lower clock, let a=(200 exp(t)−1)/199. Evaluate the decreasing beta tail at the *integer* ceil(1+(m−1)a). Its exponential factor uses the upper bound a+1/199 and its falling-factorial factors use lower bounds (a−s/199)₊. These two uses concern the same selected integer and produce a uniform lower on the original sum.

For the upper, evaluate at the smaller real size 1+(m−1)exp(t) and smaller threshold 1−exp(−τ_lo/(m−1)). Substitution 1−u=(1−x)exp(−v/B), B=N_eff−r, yields e^{-L} times the gamma ratio and the integral of (w+B(1−exp(−v/B)))^{r−1} against exp(−v), with w=B(exp(y)−1) and L=(N_eff−1)y. For r<1 the gamma ratio is ≤1 by concavity of the rth power of a Gamma(B,1) variable. Up to cutoff V, replace B(1−exp(−v/B)) below by v(1−V/(2B₀)); beyond V use its value at V. Extending the first integral to infinity overcounts, safely. For 1≤r<3 the direction switches; use B(1−exp(−v/B))≤v and the gamma recurrence plus the same fractional Jensen bound. These establish the directions of the incomplete-gamma calls.

For the odd continued-fraction truncation, take A_i=∫e^{-vt}t^i/(1+t)^{i+1−a}dt and B_i with one further denominator power, integrated over t≥0. Then A_i=B_i+A_{i+1}, and integration by parts yields v A_i=iB_{i−1}−(i+1−a)B_i, with the first term 1 at i=0. The exact odd denominators iB_{i−1}/A_i exceed v, and the even ones are A_i/B_i. Two successive reciprocal updates preserve order. Starting with v at odd depth gives a lower first denominator and hence an upper scaled upper-gamma value. The small-v lower incomplete-gamma Taylor sum is truncated on its lower side before subtraction. At 4≤N<60 the alternative monotone quadrature uses an analytic upper over [0,1/64] for the negative power, left or right values according to the sign of r−1 thereafter, and a one-sided tail beyond the final node. At N=60 the larger-size method is already valid.

I checked the scale constants, signed division rules, bounds on the exponential remainder and the log atanh remainder. Both code paths double the log remainder as well as the sum. There is a small rounding corner near powers of two: ceil-halving can trigger a final reduction with true reduced argument z just below one, while the implementation clamps the lower reduced log to zero. At the last reduction z>1−1/(2S₂), hence log z>−1/S₂. The implementation's floor of log 2 has at least 15 spare ticks, by a rational truncated atanh series. Since k≥1 here, the full lower k·floor(log 2) still bounds the true logarithm. This also treats the upper endpoint, where clamping upward is conservative. The two-pi interval used with the Stirling remainder has an independent rational Machin check. For the 42/45 return-one-tick shortcuts dominate all shapes by integer shape 3 and use e^{-L}(1+cL+(cL)²/2), with c=1.3 for y≤1/2 and c=4 for y≤2. This decreases for L≥42 and is below 10^{-12} at the respective cutoff. All of these tests used the exact rational domain of the calls.

An additional **450 independent rational price comparisons** used a different beta formula. Write a=r, b=N−r, p=x, q=1−x and A(c,z)=Σ_{i≥0}(N)_i z^i/(c+1)_i. Differentiating p^a q^b A(a,p)/a gives p^{a−1}q^{b−1}; the complementary integral follows by interchanging a,b and p,q. Cancelling the common factor gives the upper-tail identity
```
B_{N,r}(p) = a A(b,q) / (a A(b,q)+b A(a,p)).
```
All coefficients are positive. Beyond a given term, ratios are bounded above by the larger of z and the current ratio, giving a rigorous geometric remainder once this bound falls below one. This computation needs neither gamma values nor quadrature. The stock lower comparisons use the binomial identity. Test sizes include 59/60, 199/200/201, and uniform-clock probes at physical m=200,257,503. These checks independently test selected primitives. They are distinct from the complete finite table replay.

**HOLDS — §6, including unbounded degrees and chain lengths.** For a site of old degree d with l leaves and core receive r, s=(l+r)/[d(d+1)] and a leaf coefficient is t=1/(2d). Its serial credit and all its leaf serial credits give the first two terms of C. If s>t, the adverse host has l−1 surviving original sibling leaves; if t>s, it is a degree-two host with zero siblings. This gives exactly the two directed leaf debits. On a direct core edge only one direction can incur a debit. Charging b₀ plus each site's excess above b₀ covers that edge. Contracting the bare degree-two paths gives |S|−1 site edges and Σ(c−2)=−2. Port fees and matched transfers cancel in the total, leaving a strictly positive surplus 2a+b₀ after the verified local inequalities.

At fixed degree and number of leaves the credit is linear minus nonnegative multiples of positive parts in r. Subtracting the directed excess preserves concavity for fixed direct multiplicities. Low-degree refinement enumerates the exact multiset of direct degrees 2,…,9 and one class ≥10, bounding its aggregate reciprocal receive at both endpoints. The lower coefficient for a neighbour of known degree e is 1/[d e(e+1)]. For high degrees the joint simplex shortcut requires identically zero direct excess. Otherwise the code iterates the chain count and checks the endpoints of the remaining payment count and receive. This avoids a spurious joint-concavity assumption for a product of the chain count and excess.

For contracted paths with endpoint degrees d_x,d_y put q_i=1/(6d_i). A one-vertex chain has coefficient q_x+q_y. For longer chains the end coefficients are 1/12+q_i; all strictly interior coefficients are 1/6. Extending a chain of length ≥3 inserts one equal-coefficient vertex with positive serial credit, giving the length-three upper bound on debt. Convexity of positive parts allows endpoint receive and q bounds. At a large port the forward term sU is at most u_tail/[2(d+1)] for l=0. For l≥1, s≤(d+l)/[2d(d+1)] and U≤min(u_tail,1/[(l+1)M]); hence sU≤min(u_tail/(d+1),1/(4dM)), maximized for the latter expression at l=1. This yields the exact stock-coupled port cap in the code. For d>300, s≤1/(d+1) absorbs direct excess in b₀. Since s−t≤(l−1)/[2d(d+1)], the core-side leaf debit is at most min(1/[2(d+1)M],u_tail/2). Splitting at 2l≤d or 2l≥d+1 gives respectively the large core tariff reserve or positive serial reserve. The integer test covers each case. The finite-m ledger can instead stop at the physical maximum degree m−3. In an n-vertex old tree a site with d neighbours, c of which persist in the core, requires at least 1+d+c distinct vertices. An internal bare chain of length L with two nonbare nonleaf endpoints requires at least L+4 vertices. These give the two further finite-size pruning tests in the code. The additional integer audit explicitly asserted positive bare-node serial credit every time the finite chain calculation was reached.

For a star d≥2, s_center=1/(d+1)>1/(2d)=s_leaf. The normalized total is bounded below by (d/(d+1)+1/2)g_base−[(d−1)/(2(d+1))]U(d+1,d−1). At d>300 the first coefficient exceeds 1.4 and the latter debit is bounded by U(302,300)/2 using the analytic invariant-branch cap. Its j=d−1 is measured *after* the opposite leaf was promoted. It differs by one from the exceptional founder leaf. Finite stars and the two-vertex old edge are checked separately.

**HOLDS — §7 and the high-threshold match.** Every old vertex contributes w/(d+1) to A, so A≥1. The total absolute edge imbalance is at most A, and the leaf-originating share is ≤A/3 by summing the elementary leaf inequalities. Ordered upper shapes 1/2 and 1/3 give debit coefficients 4/3 and 2/3. For s=1−r, B=b+s, v=(1−x)/(Bx), substituting the affine beta coordinate gives the factor v^s/Γ(r) times an expectation against Beta(1,B). The quadratic Taylor upper for (1+z)^{-s} and the two exact beta moments yield 1−d_s v with d_s=(7/9)s(1−(1+s)/6). On 0<v≤1/6 the derivative of v^s(1−d_s v) is positive. Independent rational integration gives a scalar margin ≥0.005039711407. Thus τ≥6 gives the claimed strict gain exceeding 0.0019(1−x)^b/[n(n+1)]. The analytic proof includes equality τ=6.

**HOLDS — §8 and every threshold and label seam.** In the T₃ branch where 3 attaches to 1, labels 2 and 3 are symmetric leaves. In the other half, the next birth produces a center-versus-leaf star term with total weight (1/2)(2/3)=1/3, a middle-versus-leaf path term with weight (1/2)(1/6)=1/12, and a symmetric path with zero gap. The founder star's tag has k=3,j=2 and requires c₂. Its center can track two of three leaves. The two positive size-four comparisons are tested separately. At y≥2, the size-four analytic price comparison with τ=3y≥6 gives both signs directly. At size n=2 in the conditional theorem, the equal s coefficients on the old edge cancel its debit.

For m≥200 the accepted integers cover τ from 0 to 600000/100000 exactly. For 4≤m≤199 the accepted y boxes cover 0 to 39800000/19900000=2. They meet at the integer numerator 600000 corresponding to y=6/199. If a>0 is a lower endpoint numerator with denominator D=19900000, the checked sizes are exactly those with (m−1)a<6D, i.e. m≤floor((6D−1)/a)+1, capped at 199. At a=0 all sizes are checked. Omitted sizes satisfy the high-tail inequality throughout the box. This gives full coverage, including b=198/199, y=2, x=0 and every limit x↑1. Since x↦−b log(1−x) and x↦bx/(1−x) each map [0,1) onto [0,∞), there is no residual interval in either system of units.

**HOLDS — conditional and marginal conclusions separately.** The signed lower is pointwise positive on each old tree; averaging preserves strictness for a fixed b and x. The founder decomposition is a marginal argument. Transitivity at fixed threshold covers every pair of different nonfounder labels and then the founders versus later labels. More generally, for k<i and 3≤i<j, taking the tower conditional on a common old T_k (k≥2) preserves strictness for every adjacent step i,…,j−1. Summation gives strict order for the full pair conditional on that earlier tree. Conditioning additionally on the placement of either new label changes the claim. On a completed twin placement the tails of the two new labels even coincide. The proof has no such quantifier.

## 4. Evidence and disposition

**Prior versus new.** The pinned public baseline already contains positive-mass/moment ordering, finitely many full comparisons, a Poisson-stock endpoint lower, pre-placement signed cancellation, adverse sibling races and a geometric ledger. Those components recur in both approaches. Neon's two-clock affine bridge, physical stopping prices and exhaustive finite-size joining argument provide its route from the scoped prior results to an all-threshold conditional claim. Astra uses persistent physical senders, a different pending price, separate threshold conversion and early rooted transport. This review audits the Neon chain itself. Re-executing its binaries gives reproducibility evidence. Human verification remains outstanding.

The following are the outputs of reviewer-written assertions, using the producer's tables as inputs where indicated. Each command exited 0.

```
# Run in the reviewer's q3-audit suite (not exported in this minimal checkpoint):
python3 -IB check_tables.py
PASS all exact contiguous intervals=52; endpoint numerators 0..600000 denominator 100000
PASS all generator checks (continuing lower, of which next-frozen, stopped lower, upper) = 50195933 79929 4519459 53268480
PASS small exact contiguous intervals=84; endpoint numerators 0..39800000 denominator 19900000
PASS small generator checks (continuing lower, of which next-frozen, stopped lower, upper) = 48033670 66919 4181114 51138048
PASS finite accepted (box,m) pairs=11307; attempted including subdivided parent boxes=12105; difference=798

python3 -IB check_finite.py
PASS all 84 founder path inequalities; minimum raw margin 779708298/10^12 on [34900000,39800000]/19900000
PASS all 84 founder star inequalities; minimum raw margin 20252629370/10^12 on [34900000,39800000]/19900000
PASS integer-size physical generator checks (lower, upper, static stopping checks) = 14635151 12655282 1551686

python3 -IB check_enclosures.py
PASS integer-size beta upper: 216 one-sided rational comparisons, 208 nontrivial
PASS integer-size stock lower: 108 one-sided rational comparisons, 78 nontrivial
PASS uniform beta upper: 72 one-sided rational comparisons, 72 nontrivial
PASS uniform stock lower: 54 one-sided rational comparisons, 54 nontrivial

python3 -IB check_constants.py
PASS hardcoded two-pi interval by rational Machin bounds
PASS 42/45 early cutoffs through integer-shape-three binomial majorant
PASS all log-two floor 693147180559930/1000000000000000, guaranteed spare ticks >= 15
PASS small log-two floor 693147180559930/1000000000000000, guaranteed spare ticks >= 15
```

The generator counts include attempted boxes later subdivided. The two initializer implementations also produce overlapping mathematical cases. These counts report executed inequalities. The 11,307 accepted finite pairs count each accepted box and eligible size exactly once. The producer's printed 12,105 is a correct count of all attempted ledger calls. The 798 extra attempts cannot be counted as additional accepted coverage. The smallest **certified** founder margins on [0,2], in the last box, yield the explicit marginal lower
```
P(Z₂>x) − P(Z₃>x) ≥ (4·20252629370 + 779708298)/(12·10¹²)
                  = 40895112889/6000000000000
```
for 0≤−log(1−x)≤2. Outside this interval the analytic size-four price comparison applies. This is a reviewer-extracted consequence of the accepted arrays.

The previous reviewer code and producer replay remain separately available. `rival_checks.py` (reviewer continuation suite, not exported) produced the independent high-tail lower, 720 affine overshoot identities and 12,480 cap drifts. `review_neon.sh` (reviewer continuation suite, not exported) previously reproduced both producer binaries under undefined-behaviour sanitization, both of their own primitive audits (1,226 comparisons), and three planted defects. The full receipt is in `NEON-CXX-REPLAY.txt` (internal receipt, not exported). The independent rational calculations and generator assertions have their own receipts under `q3-audit/` (internal receipts, not exported).

**Independent direct placement sum.** I generated every unlabelled tree through twelve old vertices by all one-leaf extensions and center-rooted canonical forms. For each shape, I computed both actual ordered placement masses from its reciprocal neighbour degrees and counted surviving sibling leaves in the completed geometry. I applied the producer's raw H, R and U prices to the uncancelled signed terms using exact rational coefficients, with no tariffs, degree cutoffs or chain contraction. This checks the geometry against the observed arrays for all applicable finite boxes and these sizes. The global ledger proof above still supplies the unbounded-size step.

```
# Run in the reviewer's q3-audit suite (not exported in this minimal checkpoint):
python3 -IB check_trees.py
PASS 986 unlabelled old shapes, sizes 2..12
PASS direct uncharged signed sums for 72725 shape/box pairs on 856 applicable size/box pairs
```

The complete output includes counts and exact minimum conditional lower bounds at each size. `q3-audit/TREES.txt` (internal receipt, not exported)

## 5. Reviewer corollary for the founders

**HOLDS by the following separate argument.** With q=1−x and 0≤x<1, the marginal founder gap satisfies
```
(15/272)q < P(Z₂>x) − P(Z₃>x) ≤ q.
```
In particular the upper-endpoint order is linear. This corollary was derived during the audit. It uses the physical stock and pending bounds, with fourteen finite rational rows and no threshold partition certificate. It is stronger at the founder seam than the qualitative claim under review. I make no priority claim for the estimate.

After the third birth the branch with 3 attached to 1 makes 2 and 3 symmetric. In the other half, label 2 is the degree-two center and label 3 a leaf on the same size-three path. Thus the marginal gap is exactly half their center-minus-leaf tail difference. On the center mark both leaves. The stock lower gives
```
H_{3,2}(x) = [B_{3,2}(x)+B_{3,1}(x)]/e = [1−x²+q²]/e = 2q/e.
```
Let h be an upper bound on the probability that this marked leaf ever receives a child. Apply the static recursion reconstructed in §1 through k=5, stopping at k=6 with c₁ or c₂ as appropriate. The two unpromoted continuation states at size four have U(2,0)=483157/777600 and U(3,2)=422819/972000. At size three the exact weights are 1/2 for promotion, 1/2 for sibling loss, and 2 for host growth, all divided by three. Thus
```
h = [1+U(2,0)+4U(3,2)]/6 = 13068889/23328000 ≤ 9/16.
```
Calling these two U values uses Lemma 3.3 at N=4 as stated. The step from three to four is explicit.

If q≥9/10 then h≤(5/8)q. If q≤9/10 use the physical beta upper with shapes (1/2,5/2), since the marked leaf's only sender has degree two. This tail is bounded by that of Beta(1/2,2): use the same numerator gamma variable and add an independent Gamma(1/2,1) to the denominator's second colour. The latter law has density (3/4)t^{-1/2}(1−t). Its tail divided by q equals (1−s)(1+s/2)/(1+s), where s=√x≥3/10. Both positive factors 1−s and (1+s/2)/(1+s) decrease on [0,1]. At 3/10 the product is 161/260<5/8. Thus at every threshold the size-three leaf has tail at most (5/8)q, retaining its actual possibility of no promotion even at x=0. The marginal gap is at least (1/e−5/16)q>(15/272)q because e<68/25. From the initial edge the founder itself is dominated by Beta(1,1), giving the upper bound q.

```
# Run in the reviewer's q3-audit suite (not exported in this minimal checkpoint):
python3 -IB founder.py
PASS 14 rational pending rows; U(2,0)=483157/777600, U(3,2)=422819/972000, U(2,1)=13068889/23328000 <= 9/16
PASS beta-tail/remaining-fraction upper <= 5/8 in both ranges (large-x corner 161/260)
PASS all thresholds: (15/272)(1-x) < P(Z2>x)-P(Z3>x) <= 1-x
```

Integrating the two pointwise bounds against t↦p t^{p−1} gives (15/272)/(p+1)<E[Z₂^p−Z₃^p]≤1/(p+1) for every real p>0. In particular the gap has order 1/p at large p. The pinned all-powers founder argument already establishes a 1/[2(p+1)(p+2)] lower for p≥1 using the same size-three center/leaf symmetry. The sharper large-p order here follows from the pointwise linear bound.

The universal conditional claim and the full marginal claim each receive **HOLDS** on the identified source bytes. The two false auxiliary descriptions, the omitted large-k convention, and the attempted-versus-accepted distinction have the repairs given above. I found no headline obstruction. The founder corollary has author-side mathematical and exact arithmetic checks. It was derived after the frozen proof pin and awaits independent review.
