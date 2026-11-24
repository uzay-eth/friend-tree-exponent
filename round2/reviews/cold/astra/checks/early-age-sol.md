# Cold audit: early adjacent ultimate hub probabilities

**Verdict: holds.** The early-age proof (`open/age-order/early/README.md`) certifies `p1=p2>p3>p4>p5`. I independently reproduced the zero-horizon intervals and certified the three strict differences with four extra births. The finite computations give rational enclosures for the *ultimate* probabilities; they do not identify individual hub probabilities as exact rational numbers. Equalities below are identities of the true probabilities.

## No-send bound, including stock zero and the receiver type

Fix `x→u` with sender `x` already a core, and track `j` other leaves of `x`. Initially set `B=D_x−j`; keep this constant. Growth of tracked stock by one adds a new permanent edge; losing a tag by promotion or by discarding it does not delete the edge. Thus `D_x≥max(2,j+B)` at future stock `j`, even though the full environment changes. Before the first send into `u`, when `u` was a leaf it is still a leaf: only `x` can promote it. Its rate-one initiations therefore always create a tracked immigrant at `x`. For `u` a core, ignore all such immigration. The monitored stock gain rate is `j+a` with `a=1` or `0`, its loss rate is at most `j/max(2,j+B)`, and the send rate is at most `1/max(2,j+B)`. These bounds are predictable in the full filtration and do not require an autonomous physical stock law.

For the decreasing barrier `f(j)=min(1,2/(j+1))`, set payoff one at the send. Its drift is at most zero, using `a=0` and the less favorable denominator `j+1` at `j≥1`: the resulting drift is `−1/3` at `j=1`, and `−(j−2)/((j+1)(j+2))` at `j≥2`. At zero it is immediate from `f(0)=f(1)=1`. The frozen-at-send barrier bounds the probability of a send by each fixed time; sending that time to infinity yields the all-time bound. For the sharper bound, use the increasing harmonic payoff for a finite surrogate with boundary values `0` on send and `1−f(m)` at stock `m`. With zero immigration, stock zero is also a boundary of value zero, even though the true process can still evolve there. The harmonic recurrence at a nonboundary stock `i` is

```
(i+a)(v[i+1]−v[i]) = ((i+1)v[i]−i v[i−1]) / max(2,i+B).
```

Starting with a positive arbitrary scale at the lowest active state and using the recurrence forward gives increasing values, then normalization at `m` gives the same Dirichlet payoff as the writer's product of upward-passage ratios. Actual gains contribute nonnegatively to its drift; lowering losses or sends only increases it. Finite exit or a send is almost sure for this monitored problem: on each unit interval, conditional on still being active, a prescribed run of up to `m` gains and no competing monitored events has uniformly positive probability. This follows from the lower gain rate one and upper bounded loss/send intensities on the band, e.g. by splitting off rate-one gain proposals and dominating the monitored hazards. Zero stock with no immigration is already stopped. Bounded optional stopping then gives a lower bound on the probability of no send forever by applying `1−f(m)` after a hit of `m`. This is not just the probability of surviving until the finite stock cutoff.

For a receiving leaf, promotion is necessary for hubhood, so the complement is a hub upper bound. For a bare receiving core, hubhood requires a next gain, and until then every possible source is an old core neighbour. The upper bound is a union bound with `a=0` for each neighbour. Using `a=1` in this sector can give a false upper bound: for zero initial tracked stock it would certify a positive chance of no send from an immigration mechanism that is not guaranteed. The fixture checks that the zero-stock payoff is exactly zero for `a=0` and positive for `a=1` at the same finite cutoff.

## Exact small cases and the old-edge denominator

For a focal core, thin independent per-tag rate-`1/d` death proposals to the true hazard `1/D`; retire the tag even on rejection. Preserve the rates of untracked promotions and other arrows; the total physical generator therefore stays unchanged. The autonomous conservative pair `(d,l)` has an upward degree jump with conditional probability `d/(d+1)` and a tag loss with probability `1/(d+1)` at each of its own genuine jumps. The physical degree is never below `d`. At the cutoff the rate-`(1,1/d)` families give survival probability `1−d^(−l)` and survival implies a hub by the stock and mass comparison (`fringe/hub-children.md#1-a-vertex-with-unbounded-degree-has-positive-mass`). The fixture propagates mass forward through degrees rather than using the writer's recursion, obtaining `437/512` for the center of `T3` at cutoff four. This is a lower bound on the true center probability. The simpler bound is `3/4`. Independently the spotted leaf's parent has `(j,B,a)=(1,1,1)`; the Dirichlet payoff at stock cutoff four is `108/305`. Thus `p2−p3=(H_center−H_leaf)/2≥127/2440`. The factor `1/2` is from the original two labels being equally likely to be the center. The two leaves of a fixed `T3` have exactly equal true hub probabilities, even at `h=0`, when a non-signed subtraction of their marginal bounds would generally give a nonzero-width interval.

For the old-edge calculation, note that the orientation change comes from the *first target's* reciprocal degree. Set `r_u=Σ_(w~u)1/D_w` in the old tree. For distinct old adjacent targets `u,v`, the second target intensity for `v` is `r_v−1/[D_u(D_u+1)]`; for nonadjacent old targets it is just `r_v`. Combining the two orders leaves `(r_v δ_v−r_u δ_u)Δ(u,v)/[n(n+1)]`, where `Δ` refers to the placement with first target `u`. The denominator is exactly `n(n+1)` because a unit-rate initiator is chosen among all vertices at each of the two genuine jumps. No factor of two, normalization over old targets, or reweighting by the final degree is allowed. Repeated targets give twins; a chain is unpaired and has normalized coefficient `r_u/[(D_u+1)n(n+1)]`.

As exact identities of *ultimate* probabilities, write `H(v)` for the hub probability of the named vertex in the named final shape. From `T2` the old-edge coefficient is zero and the two chain placements both give `P4`, so

```
p3−p4 = (H(P4, middle)−H(P4, tip))/6.
```

From `T3` let `F` be the five-vertex tree with a three-arm junction, two short arms and one two-edge arm. The exact identity is

```
p4−p5 = (H(F, long-arm middle)−H(F, long-arm tip))/18
       + (H(P5, next to tip)−H(P5, tip))/24
       − (H(F, short tip)−H(F, long-arm tip))/72.
```

Indeed the old center has `rδ=1/3`, each old leaf `rδ=1/4`; the two old edges contribute `2(1/4−1/3)/(3·4)=−1/72`. The two chain types give the positive coefficients `1/18` and `1/24`. The fixture enumerates both genuine target distributions before pairing, checks the oriented formulas against the enumerated weights, and separately aggregates the endpoint orbits. Exact zero at `h=0` follows whenever the two rooted trees are isomorphic (e.g. the two end vertices of `P5` under reflection). The writer's transposition criterion is sound and sufficient for twins; it can miss a genuine equality requiring other vertices to move. This costs precision only. Zero-horizon intervals that straddle zero for the last two gaps are not equality certificates.

## Replay and enclosures

The independent fixture early-age-sol.py reads no writer code. Its signed propagation uses the same physical stopping identity, but tests rooted isomorphism rather than only a transposition. This makes the four-birth enclosures slightly tighter. Rounding here is outward to twelve decimal places.

| Gap | Independent `h=0` | Independent `h=4` |
|---|---|---|
| `p2−p3` | `[.150183708474,.500000000000]` | `[.212301826490,.366757155952]` |
| `p3−p4` | `[-.005070360056,.166666666667]` | `[.024677810172,.132386004504]` |
| `p4−p5` | `[-.009421948742,.106161344325]` | `[.006926642666,.086086384333]` |

Direct enumeration to size seven, with the all-time terminal probabilities bounded there, also gives these exact-arithmetic enclosures of the individual ultimate probabilities:

| Label | Lower | Upper |
|---|---:|---:|
| `1,2` | .552778647902 | .758823481295 |
| `3` | .214896184166 | .517646962589 |
| `4` | .123944502075 | .450756097659 |
| `5` | .069412591316 | .411113744083 |

The separate marginal enclosures overlap for the last two gaps, consistent with using signed pullback for those comparisons. From the repository root, these commands exited zero:

```sh
timeout 25s taskset -c 13 nice -n 10 python3 -Bu round2/checks/early-age-sol.py
timeout 25s taskset -c 13 nice -n 10 python3 -Bu round2/open/age-order/early/check.py --horizon 0
timeout 25s taskset -c 13 nice -n 10 python3 -Bu round2/open/age-order/early/check.py
```

The first printed `exact small certificates: 108/305, 127/2440, 437/512; swap controls: edge, cherry, reflected P5` and `nonzero orbit coefficient magnitudes: T2 [1/6]; T3 [1/72, 1/24, 1/18]`, followed by the tables above. The writer printed identical `h=0` intervals and, at `h=4`, `[.186507074525,.392551907917]`, `[.022645238637,.134418576039]`, and `[.005316116377,.087696910622]` in the order shown. All writer files were left untouched.
