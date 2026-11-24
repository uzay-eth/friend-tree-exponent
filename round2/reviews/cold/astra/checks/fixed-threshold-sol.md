# Original Q3: cold audit of the fixed-threshold prefactors

**Verdict: both bounds hold for each fixed `0<x<1`.** Moreover, the same no-cross argument gives the stronger lower bound **`c_x v^(μ-2)(1-x)^v`**, using the physical mean core/birth balance. No sharp prefactor or adjacent sign at a fixed positive threshold follows. Reviewed source (`open/age-order/fixed-threshold/README.md`), SHA-256 `b94829768c867437941d6e0879d8bb4c4821259ca652a02c9105bdadcdbc4bd0`.

## Upper bound: the half-unit saving is real

At birth the physical target weight is `r_u/(v-1)`, and while the newborn remains a leaf its unique incoming arrow has rate `1/D_parent` (divide by the current size for a step). The original target degree can only grow. The single-stage promotion bound is therefore the displayed sum with weights `r_u/(D_u+1)`. For core targets use `r_u<=D_u`; for old leaf targets the sum is exactly `Σ_core L/(2D)`. This argument is unconditional over birth targets and gives no independence between the promotion size and the environment afterward.

At first promotion the old parent is permanently a core even if it was a leaf before the marked birth. Its contribution is at most `1/2`; all others contribute at most one. The chain for `D-1/2` is dominated by a real-mass unit-reinforcement urn with red `3/2`, total `m`. Counts have the same fractional part, so a one-unit gain below equality cannot overtake. The beta law for this comparison follows by the independent continuous-time gamma representation of the **urn**, or its rising factorial moments. It is not a physical conditional degree law. The comparison holds from every promotion state, including parents of high degree and degree-two marked vertices with adverse future promotion history. The density bound yields `C√m (1-x)^(m-3/2)`; with `E C_(v-1)=O(v^μ)`, summing the single-stage probabilities over `m>=v+k` costs the geometric factor `1/x`. Exponent and indexing match the advertised delay bound.

A useful endpoint refinement follows from writing the beta tail integral with `t=x+(1-x)y`, so the remaining normalized weight of `y` is `Beta(1,b)`. Jensen gives `E√t<=√(x+(1-x)/(b+1))`. Thus for an absolute `C` and `q=1-x`,

```text
Pr(M_v>=v+k, Z_v>x)
  <= C v^(μ-2) q^(v+k-3/2) (1+√(v x))/x,  k=1.
```

For general `k` the same bound holds with `v x` inside the square root replaced by `(v+k)x` (and one may replace the leading `1/v` in the proof by `1/(v+k-1)`). In detail, at each `m` the product of its delay bound and beta cap is at most `C v^(μ-1) q^(m-3/2)(1+√(m x))/m`; the last quotient decreases in `m`, so sum the geometric series. At `k=1`, the algebraic loss near one is only `q^-1/2` relative to `q^v` on this upper bound. This does not control exact near-one asymptotics at a fixed label.

## Lower bound: conditioned family and actual denominator

The two setup jumps in the source always give an inert exterior with a tagged leaf, even on the hostile placement onto a bare core, or onto an old leaf whose only remaining neighbour is a core. Keep both exterior slots in the degrees, but not in ordinary size. Poisson no-cross clocks must be **separate** from the occupied-slot Yule proposal stream: construct the inert marginal with that stream, and then add independent thinned cut-arrow clocks. Conditional survival on the inert paths is exactly `exp[-∫(1/D_mark+1/D_ext)]`. In particular the method does not assert an unconditional surviving fraction obtained by declaring a cross arrow impossible on the auxiliary path.

The tagged Yule family can be built autonomously with a unit birth rate per real or phantom tag and a uniform initiator per birth; a removed tag remains in the same family as a phantom. It can never become a real stock member again. When a real leaf is promoted, future actual off-stock initiations use separate unit clocks, not the same phantom birth clocks. This separation is needed when conditioning on the family amplitude. Conditional on `U=u` the schedule runs at intensity `u e^t`, by the gamma posterior/negative-binomial Bayes identity; choices of initiator stay uniform and independent of chronology. The focal rate-one clock is also independent. Between focal marks, a scheduled birth preserves the expected phantom fraction; each focal mark increases it by at most `1/J`. Under silence to time one, averaging over schedules and marks bounds its limiting expected fraction by `[-log(1-e^-1)]/u`. For the maximal inequality one can condition on the family schedule (or also on the focal ring times) but **not** on the initiator choices that drive the fraction. On finite horizons the fraction is a nonnegative submartingale under that conditioning; integrate the resulting bound, then increase the horizon, to get failure `O(1/u)`. A separate compensated Poisson bound over dyadic means makes `J(t)>=u(e^t-1)/2` for all `t>=1` with the same high probability. These two events bound the reciprocal integral without applying an unconditioned Yule estimate to a tilted real tree.

The marked numerator `W` need not be `<=U`: new untracked leaves, if any, can only help. The statement needed, and obtained by the phantom count, is a bound on the **positive shortfall** `(U-W)_+`. Off-stock excess `N-J_physical-1` has per-member growth at most one and jumps of size at most two at focal rings. Real leaf births preserve this excess, and the focal immigration and actual off-stock clocks retain their rates under the auxiliary-family conditioning. Starting from a clean star at time one gives the uniform `E(S_mark-W)<=2e^-1`; the normalized limits and `S_mark>=W` are physical. A union bound on the two losses and reciprocal failure gives the simultaneous event with absolute thresholds and uniform positive conditional probability for `u>=u0`. The donor-side construction only needs the exterior to have a retained tagged leaf; its other neighbours may be hostile. All cuts use a common restart clock.

The exterior Yule proposal limit `S_+` bounds the actual exterior limit. Its exact Laplace martingale multiplies the proposal clock by `q(t)=1/(1+α e^-t)` under the tilt `exp(-α S_+)`, with no change in conditional marks or occupied-slot selection. Hence it multiplies accepted inert rates by that same factor. Under the tilt `S_+~Gamma(N_ext,1+α)`. Use the marked-side test `U>α S_++(1+α)D_*` together with `W>=U-D_*` and `S_mark<=W+D_*`; then `W>α(S_ext+D_*)` implies the **actual** ratio exceeds `x`. In particular conditioning on `S_+` does not replace the marked-side size by an independent gamma. For fixed `x` the exterior bounded-reciprocal event has uniformly positive tilted probability from the retained leaf (embedded birth–death survival plus a fixed finite reciprocal budget). Intersect with `S_+` large by subtracting its vanishing tilted lower-tail probability, rather than claiming independence. The normalizer is exactly `(1-x)^N_ext`.

## Stronger lower bound by averaging immediate entrances

Let `n=v-1` and choose only targets `u` with **an old leaf before the birth**. Then mark `v` is born on `u`, and `u` promotes `v` on the very next jump. The marked side is again a virgin with two ordinary vertices, and the exterior has `n` ordinary vertices, with at least one retained leaf. Conditional on `T_n`, this setup costs exactly

```text
Σ_(u:L_u>=1) r_u/(D_u+1) / [n(n+1)]
  >= (2/3) Σ_u L_u/D_u / [n(n+1)],
```

where the inequality uses `r_u>=L_u` and `D_u>=2`. Thus setup probability is at least order `v^(μ-2)`: the physical uniform-core mean limit and birth balance (`occupation-neon/flow/README.md`) give `E Σ L/D ~ μ E C_n` and the Claim 6 mean amplitude is positive. This input is needed for the sharper power; even without it the sum is at least one on every tree with a core, giving the uniform `v^-2` entrance bound. The subsequent tilted construction is uniform over the successful setups and costs `c_x q^n` with `q=1-x`. This proves the improvement and shrinks the gap between the powers to `3/2`. It does not assume each possible birth target offers an old leaf; the exact fixtures exhibit failures of that stronger assertion. On a conditional episode that has already fixed a hostile leaf or bare target, use the source's two-step repair instead.

## Endpoints and adjacent sign

The fixed-`x` powers are not a uniform asymptotic as `x` approaches either endpoint. One can make the restriction on the lower proof more transparent: test `U>u0+α S_++(1+α)D_*`, so no assumption that the tilted `S_+` itself be large is necessary. To avoid a reciprocal budget blowing up as `1/q`, split exterior time at `t0=log(1/q)`. The clock multiplier is at least `1/2` thereafter, and before that the boundary reciprocal is at most `1/2`. The embedded one-leaf family survives to operational time at `t0` with probability at least its eventual survival probability `1/2` (the operational time there is deterministic, at most `log 2`). Restart a tagged family from one leaf then; the remaining reciprocal budget is bounded on an event of uniform positive probability, using `dt<=2 d(clock)`. Thus no-cross costs at worst an extra `√q` from this initial exterior delay. This gives a constructive lower constant of the form `c√q exp(-D_*/q)` multiplying `v^(μ-2) q^(v-1)`, valid without taking `x` fixed first, with absolute `c>0`; it does **not** identify the true near-one scale. As `x↓0` the lower constant need not vanish, but these bounds do not recover the full birth-scale behavior (the delay sum bound has a `1/x` loss). At `x=0` one is comparing hub probabilities, and at `x=1` both strict tails vanish. A bounded oscillatory multiplier can keep the same exponential rate and polynomial envelopes while reversing adjacent signs; that is not a counterexample for the actual growth process.

## Exact fixture

```sh
nice -n 15 taskset -c 3 timeout 12s python3 -Bu round2/checks/fixed-threshold-sol.py
```

Passed 436 labelled birth/promotion rows (including adverse placements and the factor `2/3`), 336 rational conditional Yule transitions, 464 phantom-fraction rows, and exact ratio, gamma-tilt and endpoint/oscillation controls. The limiting comparisons are analytic; the fixture does not estimate a rare-event probability numerically.

## 42b5 follow-up: fixed-label ordering near one

The local checkout does not resolve object `42b5`; the reviewed near-one note has SHA-256 `3ffe7e91ed5c59a1236668573131a4b048f24f4ab58109be1abc13c9e2934b78`.

**Holds.** Near-one note (`open/age-order/preplacement-transfer/positive/README.md`). The first promotion at size `v+1` has probability at least `1/[v(v-1)]` conditional on **each** birth tree. Then the leaf stock at size `n` on the append event is `n-v`, so the exponent from the finite append product is exactly `v`, not `v-1`. Conditioning from size `n` on the complete chosen history changes neither the subsequent physical transition law nor the future independent size clocks. The gamma-minus-loss expectation yields `E(S-W|T_n)<=v+1`; the union bound in the submitted proof is valid for arbitrary dependence between `S` and `W`. For the next label, earliest promotion is at `m=v+2`, so the shifted beta exponent is `m-3/2=v+1/2`. The beta cap is conditional on **any** promotion state, and can also be used given the prebirth history after averaging over future promotions. No assertion compares two post-placement conditional tails with different fixed environments by exchangeability.

Consequently the near-one rank of the original unconditional tails for fixed labels is `1=2>3>4>…` in the pairwise eventual sense. Every fixed finite list admits a common interval; there is no uniform-in-label interval from this proof. For any two fixed labels `3<=v<w`, `P_w(1-q)/P_v(1-q)=O_v,w(q^(w-v-1/2))`. Labels 1 and 2 have tails of order `q`, by the size-three centre event and the ordinary beta `(1,1)` cap. For `v>=3` the log tail power is sandwiched between `v-1/2` and `v` (liminf and limsup respectively); an exact log power is not shown. In the random-rank side result, differentiation really averages adjacent differences with Bernstein weights; one starting block cannot recover the sign of each coefficient. For an explicit admissible synthetic list, start at `1/2` and subtract the successive differences `(5,-3,5)/16`; the four values are `(8,3,6,1)/16` and all lie in `[0,1]`.

### Wider rigorous joint region

The exact size-clock independence gives a better constant than the gamma/Markov split: at deterministic restart size `n`, `S` is independent of **the embedded fraction** `Z_v`, though not of `W=S Z_v`. Therefore `E(1-Z_v|T_n)<= (v+1)/n` on the append event. Put `A_v=(v+1)^2/v` and `n=ceil(A_v/q)` with `q=1-x`. Markov gives a success probability at least `1/(v+1)`. The append denominator is at most `(n-1)^v < (A_v/q)^v` (also when `A_v/q` is an integer). This yields the coefficient below. It optimizes the bound `A^-v[1-(v+1)/A]` over `A>v+1`; this is a statement about this elementary Markov/append comparison, not global optimality of the region.

Define the deterministic averaged entrance coefficient and a later-label bound by

```text
e_v = E Σ_u r_u/(D_u+1) on T_(v-1),
b_v = e_v/[v(v-1)],
L_v = b_v v! / [(v+1) A_v^v],
U_v = [e_(v+1)/(v(v+1))] Γ(v+2)/[Γ(3/2) Γ(v+3/2)].
```

`b_v` is the exact immediate-promotion probability, counting **all** targets, since the append argument needs no retained exterior leaf. Given any birth tree a uniform version of `L_v` follows by replacing `b_v` with `1/[v(v-1)]`. For the later label at promotion size `m`, its mass is bounded by `e_(v+1)/[v(m-1)]`. In summing the beta bound, the coefficient `Γ(m)/[(m-1)Γ(m-1/2)]` decreases (successive ratio `(m-1)/(m-1/2)`), so

```text
P_v(1-q) >= L_v q^v,
P_(v+1)(1-q) <= U_v q^(v+1/2)/(1-q).
```

The best width from these two displayed majorants alone is `q<Q_v = [2r_v/(1+sqrt(1+4r_v²))]^2`, `r_v=L_v/U_v`. Taking the minimum of the later upper estimate and the promotion-independent `u_v q^(v+1/2)` from 42b5 also admits the alternative width `(L_v/u_v)^2`; use the larger of the two widths, capped by one. The formula with `e_v` can use exact finite expectations, but does not require enumerating them. In fact `e_v` is of order `v^μ`: its upper bound is `(3/2) E C_(v-1)`, while the lower bound is `(2/3) E Σ L/D` as above. Thus averaging both sides gives `r_v` of order `e^(-v)/v` by Stirling, and the guaranteed near-one width is of order **`e^(-2v)/v²`**. To obtain a wholly elementary positive width at *every* size without asymptotic constants, in `L_v` use the conditional `b_v>=1/[v(v-1)]`, and in `U_v` use `e_(v+1)/v<=1/2` (given any birth tree the permanent parent has degree at least two). The gamma quotient is `2·4^(v+1) ((v+1)!)²/[π(2v+2)!]`, so replacing π by 3 gives a rational upper coefficient. A convenient rational sufficient width is `min(1/2,(L_v/(2U_v))²)`.

The improved fixed-positive-threshold no-cross bound above does not extend the large-label sign region: comparing that lower envelope with the later beta/promotion upper envelope leaves a factor of order at most `v^-3/2` at fixed `x`. Even the stated endpoint-dependent construction cannot close the gap where those majorants apply: the extra `q` powers available near one cannot overcome its `exp(-D_*/q)` factor. Near zero the delay-sum majorant deteriorates, and the compact-scale proof is the relevant one. For each **fixed** `K` and all sufficiently large `v`, it certifies `0<=x<=K/(v-1)` by taking `n=v-1`; it gives no quantitative diverging `K(v)`. The combined certified slices are that compact birth-scale interval and `1-Q_v<x<1`, disjoint for all sufficiently large labels at each fixed `K`. This leaves every fixed `x∈(0,1)` eventually outside the joint certified region. No new adjacent sign is inferred there.

The CPU 3 replay also checked conservative rational widths for labels `3,4,8,16,32,64` (for example `7.3098e-6` at 3 and `3.2890e-12` at 8, rounded down here), including the ceiling in `n` and the lower success factor. Those explicit widths use the worst-case entrance coefficients rather than the averaged `e_v`.
