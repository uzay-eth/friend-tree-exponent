> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes; earlier growth, spatial, spectrum, birth-label, and fixed-shape results. See the [release index](../../../../../README.md#claim-23) for exact scope and citations.

# Fixed shapes at a proportional degree cutoff

**The marked Poisson conjecture holds.** Let the original friend tree start from the edge, and let the deterministic cutoffs satisfy `N_n/n → c > 0`. Delete vertices with **original** degree `> N_n`. Count only components containing an original core (original degree at least two), and root each at its highest vertex in the ordinary birth orientation (root the initial edge at label 1). Write `Z_v=lim_n D_v(n)/n` and `H_c={v: Z_v>c}`. For any finite set of distinct rooted ordinary shapes `R` of size at least two, the marked counts converge **stably relative to the full physical history**: on an extension, conditional on that history, take independent

```math
P_{h,R}\sim\mathop{\mathrm{Pois}}\nolimits (I_R),\qquad
I_R=\int_0^\infty \mathbb P(T^0_u=R)\,du,\qquad h\in H_c.
```

Here `T⁰` starts at the inert virgin (root plus one real leaf, one exterior slot in the root degree), and an atom of multiplicity `P_{h,R}` has parent label `h` and parent mark `Z_h/c`. Thus the unmarked totals are conditionally independent `Pois(|H_c| I_R)`. This holds also for the full size-`|R|` component counts: eventually every counted occurrence of such a fixed shape is intact below its missing parent edge. For example `I_pair=2/3`. If `c≥1`, the limiting list of hubs and all the fixed-shape limits are empty.

The singleton-core count is different: it is eventually exactly the finite number of bounded original cores all of whose **final** neighbours lie in `H_c`. It can be carried jointly with the Poisson limit; a nonfounder singleton has parent label `h` and limiting parent mark `Z_h/c`. For the possible orientation-root singleton, use no parent mark. Ordinary degree-one leaves isolated by deletion are excluded throughout.

## 1. Every positive deterministic threshold is regular

Use the physical rate-one vertex clocks. The total size `M(t)` is Yule, `e^{-t}M(t)→Y>0`, and its size-`n` epoch satisfies `τ_n=log(n)−log(Y)+o(1)` almost surely. Then `e^{-t}D_v(t)→YZ_v`. We use the [hub dichotomy and fixed-hub offspring theorem](../../../../../../astra/source/round2/fringe/hub-children.md): unbounded degree implies positive mass and every positive-mass vertex has a positive-mass child. The second input follows in particular from the positive logarithmic offspring coefficient in (5) of that note. The rate of arrows from a hub into its existing core neighbours is `c_v(t)/D_v(t)`, where `c_v(t)=O(1+t)` almost surely by the personal clock bound. This rate is integrable on the hub event. The compensated counting-process criterion (stop at accumulated intensity `m`, then take the union over integer `m`) gives only finitely many such arrows there. In particular the edge between two hubs eventually carries no arrow in either direction.

We need atomlessness for the **embedded fraction**, not merely for `YZ_v`. Fix a label `v` and `x>0`. On `{Z_v=x}`, select a hub-child edge of `v` and some finite integer size `m` at which both endpoints are cores and after which the edge carries no arrow. Take the countable union over `m` and such labelled edges. Conditional on a finite size-`m` history, cut the specified edge and retain each endpoint's slot. Generate independent inert continuations on the two sides `A,B` (`v∈A`). Up to the first cross arrow these have the same internal trajectories as the physical continuation: add the two suppressed directions by extra thinned clocks with rates `1/D_boundary`. The marginal inert trajectories are unconditional. It suffices to show that the inert ratio is atomless at `x>0`, since on no cross the physical and inert ratios coincide.

For an inert side with `m_i` initial ordinary vertices, introduce a proposal Yule process `Q_i` starting at `m_i`. At a proposal choose a uniform slot in `Q_i`; if it represents a real vertex, take its physical neighbour choice with the exterior slot retained, otherwise do nothing. Increment `Q_i` at every proposal and the true size only at genuine adds. This projects to the inert rates and makes the proposal **embedded chain** independent of proposal holding times. For `G_i=lim e^{-t}Q_i(t)`, the law is `Gamma(m_i,1)`. The limits `b_i=lim M_i(t)/Q_i(t)>0` and `a_i=lim D_v(t)/Q_i(t)≥0` (the latter on `A`) exist: degree has the usual supermartingale limit, and the [inert size limit](../inert/README.md#3-ordinary-size-and-the-high-cutoff-joint-law) is strictly positive. Along the proposal epochs these ratios are functions solely of the proposal embedded chain. Thus, conditional on the two proposal chains, `(G_A,G_B)/(G_A+G_B)` is Dirichlet with parameters `(m_A,m_B)`. The desired ratio is

```math
\frac{a_A G_A}{b_A G_A+b_B G_B}.
```

It is zero when `a_A=0` and otherwise strictly increasing in `G_A/(G_A+G_B)`. Hence it has no atom at `x>0`, even after intersecting with the no-cross event. This proves fixed-label positive atomlessness without treating the physical size and degree masses as independent.

At ordinary birth size `v≥3`, dominate `D_v(n)` by a two-colour Pólya count with one red out of `v`, using the physical gain probability bound `D_v(n)/n`. The urn proportion is a bounded martingale with `Beta(1,v−1)` terminal value. For any `p>1` the nonnegative submartingale maximal bound on its `p`th power gives

```math
\mathbb P\big(\sup_{n\ge v}D_v(n)/n\ge a\big)
\le a^{-p}\frac{\Gamma(1+p)\Gamma(v)}{\Gamma(v+p)}.
```

For `p=2` this is `2/[a² v(v+1)]`. Borel–Cantelli with `a=c/2` leaves only finitely many candidates ever, to which apply the pointwise limits and atomlessness. Almost surely the deleted labels equal `H_c` for all sufficiently large `n`; at late physical times the same is true when using `N_{M(t)}`. No independence across label urns is needed.

## 2. The old components

Work pathwise once the deletion set is the fixed finite `H_c`. The component of any fixed retained vertex grows monotonically from the time all deleted vertices exist. If it has finite final size at least two, two adjacent members would both have bounded original degrees: their only other possible neighbours are in `H_c`. The arrow rate across that edge would stay bounded below, contradicting bounded degrees. Thus any eventually bounded component is a singleton.

A component containing a core and not lying strictly below a high parent is the orientation-root component. The other components that have a cut somewhere below their top have tops on the finite paths from the orientation root to `H_c`. These fixed exceptional components therefore give no size-`s` count eventually for `s≥2`. Every remaining component containing a core is an intact descendant branch of an actual promotion by a high parent, with the potential initial-edge branch handled as another fixed exception. A core singleton has at least two high neighbours, so it belongs to the finite subtree spanning `H_c`: it is the middle of a length-two path between its high neighbours. Once all `H_c` are present, these possible singleton cores belong to a fixed finite set. Their indicators are eventually exactly the required final-neighbour indicators; a future gain to a new vertex is a gain to a retained vertex and ends singleton status. This proves the finite, eventually constant singleton formula and also handles `H_c=∅`.

## 3. Late branching with a finite-history forecast

Fix a finite list of shapes of size at most `s`, a horizon `l>0`, and `0<ε<1`. Let `F_b` be the physical history at an early deterministic time `b`. Use the proxy set `H^{b,K}` of **cores** with labels at most `K` and degree fraction `>c` at time `b`, and forecast

```
Y_b = e^(-b) M(b),    t_b = log(n) - log(Y_b),
a = t_b-l-ε,          u_- = t_b-ε,       u_+ = t_b+ε.
```

The window is `F_b`-measurable and can be used after `b` (for fixed `b,l,ε`, `a≥b` eventually). First localize at a deterministic label cap `K` with `P(H_c ⊄ [K])` small; as `b→∞` the proxy equals `H_c∩[K]` with probability tending to one. Also `lim_{b→∞} limsup_n P(τ_n∉[u_-,u_+])=0` for each `ε>0`. The predicted time need not be independent of the past.

On the window `[a,u_+]`, use each proxy's fresh rate-one initiating clock as its immigration proposal clock. At every point independently preassign a virgin trajectory, even if the physical ring will fail to promote. At a real promotion, use this trajectory's fresh drivers for the promoted leaf and its new leaf and for every later arrival inside the branch. The outgoing choice to the parent is null internally, but the send to the real parent is still made. This gluing is legitimate at promotion: a leaf has not seeded a core branch and its unused post-promotion drivers can be refreshed. Distinct selected parents are already cores at `b` and cannot be selected as promoted leaves. Future new child branches of these parents are geometrically disjoint; rule out using any young port twice by the rule that only degree-one vertices are eligible. One may preassign all inert paths using independent drivers and in the live process reveal only the clocks and choices as needed in causal order. Conditional on `F_b`, at `a` the proxy clocks are independent rate-one Poisson clocks, with independent inert marks for every point. This is a statement about the comparison law conditional on **early history**, not the law of future innovations given `H_c` or given that the coupling succeeds.

A proxy ring choosing any existing core is a failure; this includes every unwanted incoming arrow into its matched old or newly created branch. Degree and constant-rate clock bounds give, for each real hub,

```math
\int_a^{u_+} \frac{c_h(r)}{D_h(r)}\,dr\longrightarrow0\quad\text{a.s.}
```

To apply compensation, integrate the **predictable** rate over the `F_b`-chosen parents and interval without a future-history gate. The integral is bounded by `K(l+2ε)`; on the event that the proxy equals `H_c∩[K]` it converges to zero. Thus its limsup in expectation is at most `K(l+2ε)P(H^{b,K}≠H_c∩[K])`, which vanishes with `b` at fixed `K,l,ε`. Separately discard paths with a hub beyond `K`, at their probability cost. Stabilization ensures there is no high descendant in any branch promoted in the window at late enough `n`: such branches start with just the promoted leaf and its newborn, and all high labels are already present as cores. On success the late physical counts are the inert marked counts evaluated at `τ_n`.

### Bounds for omitted ages and for the random observation time

An intact branch of age `r` and size at most `s` has had at most `s−2` intrinsic changes, with a genuine rate at least `3/2` throughout. The upper probability `B_s(r)=P(Pois(3r/2)≤s−2)` holds from any actual virgin promotion even with adaptive live input. Before `b` there are finitely many already promoted branches of the first `K` parents, and each grows past size `s` in finite time almost surely. For later promotions before `a`, a contribution at any point of `[u_-,u_+]` requires size at most `s` already at `u_-` and age at least `l` there. Conditional on `F_b`, use the uniform bound at each promotion and then its rate at most one to bound their expected count by `K ∫_l^∞ B_s(r)dr`. The time endpoints here are known at `b`; no conditional bound at a future-selected birth is used. The integral tends to zero as `l→∞`.

For the ideal branches starting at `a`, the number with size at most `s` at any time in the window has conditional mean at most `K ∫_0^∞ B_s(r)dr`. The total intensity of genuine jumps in such branches is at most `s` per branch, and new immigrants arrive at rate at most `K`. Hence, even at the random observation `τ_n`, the probability of any fixed-shape count differing from that at `t_b` while `τ_n∈[u_-,u_+]` is at most

```math
2\varepsilon K\left(1+s\int_0^\infty B_s(r)\,dr\right).
```

This bounds **changes of the observable**, not the number of ordinary tree jumps (the latter is of order `nε`). No independence between the size epoch and the ideal observations has been used.

## 4. Stable limit and parent marks

At `t_b`, the conditional Poisson marking law for the ideal system gives independent counts for each proxy and each distinct shape, with means `I_R(l+ε)=∫_0^{l+ε} P(T⁰_r=R)dr`; a mark can be assigned the `F_b`-known value `D_h(b)/(c M(b))`. Each trajectory is classified by exactly one shape at the observation, even if the same trajectory visits several of the requested shapes at other times. At size `n` the physical mark of any contributing hub converges almost surely to `Z_h/c`. Under the above coupling the counts are tight; the finite label cap and a uniform-continuity bound on the mark interval show that replacing the physical mark by its early forecast has vanishing effect as `b→∞`.

Explicitly, test the marked vector by bounded Laplace tests and multiply by any bounded `F∞`-measurable variable `U`. Replace `U` by `E[U|F_b]` with an `L¹` error tending to zero. The ideal conditional Laplace functional has factors

```math
\exp\!\left(-\sum_{h\in H^{b,K}}\sum_R I_R(l+\varepsilon)
 \left[1-e^{-f_R(h,D_h(b)/(c M(b)))}\right]\right)
```

for bounded nonnegative mark tests continuous in the relative degree. Take `n→∞` first, then `b→∞` at fixed `K,l,ε`. Let `ε↓0` and `l→∞`, then remove the label cap; its error is at most a constant times `P(H_c ⊄ [K])`. These operations prove stable convergence to the stated conditional kernel. In particular this is conditional independence **in the limit**, not a finite-time claim of independent immigrant streams under conditioning on the future high set. Since the physical singleton counts stabilize and their finite marks converge, they may be included in the stable test. The `I_R` can also be computed by the [finite shape recursion](../fixed-shapes/README.md): the integral equals the hit probability times `1/(|R|−1/D_root)`.

## 5. Relation to the slow diagonal

The Poisson marginal does not by itself prove infinitely many empty samples. The separate [slow-diagonal argument](../../../../../../../reviews/cold/astra/checks/forest-fine-diagonal-sol.md#3-a-deterministic-slow-diagonal-that-fails-as) uses pathwise recurring empty intervals of a dominating queue built from finitely many **deterministic-label** clocks before restriction to `H_c`. The fixed-cutoff regularity here removes its atom-exclusion restriction for `c>0`. No previous diagonal file is edited.

## Exact finite fixture

From the repository root, `python3 -B round2/open/forest-age/macroscopic/check.py` reports:

```text
2532 exact original-degree masks; orientation-root singleton verified
210 exact marked Poisson factorial cancellations; 18 clock-shift budgets
19 exact inert Yule slot generators and strict Dirichlet ratio crossings
```

The masks are enumerated on chronological trees and checked by a direct connectivity search after deletion. The Poisson controls sum over separate categories, including the unused category, and the split controls check that a rejected boundary arrow changes the auxiliary clock only. The fixture checks these finite identities; stability, threshold atomlessness and the random-size transfer are proved above.
