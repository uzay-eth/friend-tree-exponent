# Independent check: parity influence along old physical edges

**Holds for real common edges already present at the transmitting vertex's first infection.** At discount `s=1/2` the uniform expected discounted excess intensity through the eligible outward edges is at most **0.992198062**; at `.565` it is at most **0.838093404**. This does not contract the number of all future infections: an explicit counterexample is below. The check is on the physical or inert-slot intrinsic coupling, not on the marked residual branching family.

## Coupling and conditional hypotheses

Start with two copies of one finite ordinary tree, with one extra leaf at a paired vertex in copy A. Every initial ordinary vertex, including every degree-one leaf, is paired. A fixed inert slot, if present, is common and counted in the endpoint's degree, but neither sends nor receives a changing event. Color the common tree by parity from the initial discrepancy. Even vertices are favored in A, odd in B. For each oriented real common edge synchronize at the minimum of the two rates `1/D_source`; give the surplus to the smaller-degree copy of the source. The latter is the favored copy of the *target*, whose color is opposite. A simultaneous event attaches a fresh paired leaf of the opposite color to its target. A single-copy event into a common target attaches a private leaf. A paired source selecting a private neighbour causes growth solely within that private branch. A private source can send into the common tree only over its unique attachment edge; that target is already favored in the same copy. Private branches cannot reconnect elsewhere because all events add one new vertex and one edge. Thus the common tree stays connected, and at each common site all private neighbours are on its favored side. The event rates are the true directed marginals, including for paired leaves that become cores on different histories. Finite-time nonexplosion follows from the ordinary size birth rates, no greater than size.

A common site is infected when its two degrees first differ. Extra neighbours are never removed or re-paired in this construction, so infection is irreversible. At a first infection away from the initial seed, the degree gap is exactly one and the infector is a previously infected real common neighbour. Neither a private send nor a common simultaneous event can infect a clean vertex. The infected set is connected and the common graph is a tree, so this infector is the **only** previously infected common neighbour. These statements apply conditionally on the full filtration at infection time, including to initially clean sites that grew before being infected. No conditional independence of later rings is needed. Put `d` equal to the lower degree at that stopping time, including the inert slot if any. At most `d−1` real outward neighbours then exist after excluding the infector. At the initial discrepancy there is no infector: use at most `d` real neighbours instead, and the same per-edge estimate; the nonseed contraction need not hold at the seed.

Let `I(t)` count ever-infected common neighbours of this site `x` after the stopping time, including any newly appended paired neighbours. Any increment has to arise from a differential send by `x`; the component beyond an uninfected neighbouring edge contains no infected site (tree connectivity). The combined intensity of these sends is bounded by the total of the lower-degree source copy's arrows, at most one. Consequently `E I(t) <= 1+t` conditionally. A common neighbour can send repeated excess returns with no new infection, but their combined rate into `x` is at most `I(t)`. Private neighbours can send repeatedly into `x` at rate at most one each; their number is precisely the gap `delta`. All other arrivals leave the gap unchanged. Therefore the conditional mean of `delta` satisfies `u' <= u+1+t`, with `u(0)=1`, giving `u(t)<=3 exp(t)-t-2`. This uses a stopped drift inequality and the integrable Yule domination before removing the stop. Including a silent slot only decreases the rates.

For an individual already-common outward edge the excess intensity is `1/D_low − 1/(D_low+delta)`. Lower degree never decreases, and `delta/(d(d+delta))` is increasing and concave. Jensen gives the second envelope below. For the first, the total instantaneous gain rate of `D_high` is at most `D_high` (one inverse-degree arrow per real incident edge), so `E D_high <= (d+1)exp(t)`. Convexity of reciprocal gives `E[1/D_high] >= exp(-t)/(d+1)`. Thus

```
E(excess intensity) <= min( 1/d − exp(-t)/(d+1),
                            B(t)/(d(d+B(t))) ),  B(t)=3 exp(t)-t-2.
```

The envelope concerns all excess sends, including after the receiving site is infected; retaining those sends is a valid upper bound on the first-infection measure by the predictable intensity identity. In particular conditional expectations of the sum of `exp(-s*(T_y-T_x))` over first infections through eligible edges are bounded by the stated coefficient. The discount is measured from the actual infection time of the sender.

## Exact arithmetic and limits

`bound.py` uses rational Taylor sums through degree 12 with geometric upper remainder for each step exponential; it then maintains enclosing integer grids scaled by `2^64`. At a dyadic mesh point `t` the `-t` term is represented exactly. Both candidate gap envelopes, and their pointwise minimum, are nondecreasing in `t`; the discount is nonincreasing. Each integration cell uses the appropriate opposite corners for its lower and upper products. The omitted time integral is at most `2 exp(-T/2)` for `s>=1/2`. For **all** `d>=128`, split the integral at 4: the degree-multiplied intensity is at most `3 exp(t)/128` until then, and at most one thereafter; this gives `(6/128)(exp(2)-1)+2 exp(-2) <= 0.570157572`. `d=1` has zero old outward contribution. The lower numerical endpoints below are for the **analytical envelope**, not for actual influence.

| Discount | Gap envelope | Certified interval for its supremum over integer degrees |
|---|---|---:|
| `.5` | tight, `1/1024` grid | `[0.991355217, 0.992198062]` |
| `.525` | crude, `B=3 exp(t)-2`, `1/256` grid | `[0.960821624, 0.964260985]` |
| `.565` | tight, `1/1024` grid | `[0.837294684, 0.838093404]` |
| `.5` | crude, `1/256` grid | `[1.025972861, 1.029478862]` |

**Concrete limitations.** `coupling.py` produces a degree-one newly infected paired leaf, with zero old outward neighbours, that later infects a paired leaf appended after its infection. Hence replacing “already-common outward edges” by “all outward edges ever” is false, even though the mean-gap estimate itself counts and pays returns from such edges. The checker also exhibits strictly increasing discrepancies from a repeated common-neighbour return and then a private-neighbour return, so a single-gap state approximation is not closed. For the **fixed initial common skeleton**, newer growth cannot reconnect distinct old vertices, and the bound can be iterated along the members of this skeleton using the conditional inequality at each infection. It cannot on its own be iterated over arbitrary newly created physical branches, nor does paired core status or equality of permanent birth marks follow from parity. A transfer to Claims 2 and 5 must separately pay mark changes, killing and signed copy differences, or give a physical-to-residual identity that accounts for them.

More explicitly, let `T_x` be the infection time of an initial paired vertex, `T_seed=0`, and `H_k=E sum_{dist(seed,x)=k} exp(-s T_x)` with the sum restricted to the initial skeleton and infinite infection times paying zero. Along that skeleton every first infection necessarily uses the inward edge. Conditional on each infection time, all skeleton children are already-common neighbours; summing over vertices yields `H_{k+1}<=q(s) H_k` for `k>=1`, with the certified `q(s)` above. For the seed use `H_1<=1/s` by the total differential-send intensity bound. This corollary counts infection events on the fixed skeleton; it does not claim a bound for the state or future descendants of each infected vertex.

## Replay from repository root

```
taskset -c 13 timeout 14s python3 -Bu round2/influence/check/coupling.py
taskset -c 13 timeout 14s python3 -Bu round2/influence/check/bound.py
taskset -c 13 timeout 25s python3 -Bu round2/influence/check/bound.py --refined
taskset -c 13 timeout 14s python3 -Bu round2/influence/check/bound.py --challenge
```

The finite audit reconstructs the two exact single-copy rate maps, checks their coupled marginal rates and every outgoing transition on 432 coupled states each for no inert slot, an inert slot at the seed, and one at its neighbour. It rejects swapped surplus ownership and active-slot rate inflation. This is an independent audit of the conditional physical estimate; it supplies no eight-digit exponent claim.
