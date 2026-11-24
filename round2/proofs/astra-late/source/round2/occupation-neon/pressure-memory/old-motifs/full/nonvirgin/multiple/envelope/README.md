> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../../../../README.md#claim-21) for exact scope and citations.

# Shared all-copy exit envelope and D18 cut obstruction

**The direct D18 transport fails pointwise:** the exact upper witness, split by root ownership and normalized at the virgin, has a positive residual `0.606151101368…` on a reachable fringe with heterogeneous old grandchildren. A separate genuine marked-fringe weight works for all copy counts and all remaining horizons. Its scale is 13.1% smaller than the plain degree-weight scale, but it still does not resolve the short return. A prior-assisted scalar comparison at `s=.56572365` gives **`F(s)∈[.9691216940,1.0189758872]`**; it determines no sign.

## Exact transport formula and obstruction

Write a marked fringe `t` with mark index `i`, root degree `D`, and children `y`. In the physical process the child `y` retains its birth mark index `j_y`; the true parent send is `p=1/D`, while `a_y=1/j_y≥p`. Let `P(t)` be the original upper D18 potential of the planted tree, using the exact physical census and planted parent slot. Define the root-owned allocation `A(t)=P(t)−Σ_y P(y)`. A motif entirely inside a child has the same physical labels before and after this cut, so it cancels exactly. A leaf gain of `y` leaves its child list unchanged; a promotion at the root makes one new `J_min(D,20)` regardless of whether the corresponding separate production is a stopped trigger.

For the intrinsic generator `L⁰` (no exterior arrival on the planted root), this gives the exact cut formula

```
L⁰ A(t) = L⁰ P(t) − Σ_y L⁰ P(y)
          − p Σ_y [P(y⁺)−P(y)] − b(t) P(J_min(D,20)).
```

The last two terms are the actual sends into the child parts and the new part created at promotion. For the true stopped-degree-two generator one must then add

```
(R̂−s) A(t) = L⁰ A(t) − (s+1/i) A(t)
            + Σ_y [(a_y−p) A(y)+p A(y⁺)]
            + b(t) 1{D>2} A(J_min(D,20)).
```

The sum of the physical global inequalities does not supply the **difference** of generator defects in the first formula. Moreover the cut terms in the second formula are computed at the actual marks, even though `A` itself is mark-blind. On an old state the same test applies to the first-virgin equation with terminal values equal to the virgin trace of `A`: a copy landing on a virgin or a nontrigger newborn then takes exactly the value used in the formula. Dividing by `A(J₂)>0` normalizes all virgin marks to one, with no change in the pointwise residual's sign.

Write `(i,L,children)` for the marked tuples. The exact bad row at the trial is

```
(2,1, ((2,4, ((3,0,((2,1,()),)), (4,3,()))),
       (3,2, ((2,1,()),))))
```

This is reachable: the root promotes first at degree two and next at degree three, followed by the final gain to degree four; each side arm can gain and promote in birth-mark order, and the bare middle of degree two promotes at degree two. A finite prescribed order of intrinsic arrows has positive probability. On this row the normalized value is `.85590171088`, intrinsic term `.19090551634`, copy term `1.07740028039`, direct term `.25`, and death term `−.42795085544`. After subtracting the discount the exact residual is `1591146641093019431/2625000000000000000>0`. Five of five heterogeneous-old test rows fail this supersolution direction. Changing its exterior mark index from 2 to 20 without changing its tree adds exactly `(9/20) A(t)/A(J₂)=.385155769896` to the residual. The birth mark of the root can differ through its original exterior parent, with this same intrinsic geometry subsequently formed. This rejects the natural root-owned transport with its own terminal trace, including an overall positive rescaling. It does not reject all mark-dependent correctors or all possible choices of terminal vector. Replay uses the original eight-family evaluator, not just `f/g` rows.

## One certified positive envelope for all-copy first-virgin action

For the stopped-at-next-virgin excursion of the positive comparison (`round2/sol-moments/positive-comparison/README.md`, historical reference), the pre-hit particle generator is `R`, without nontrigger newborns in the generator. Those emissions instead earn terminal vector `w`, as do old-child copies hitting a virgin. Use one nonnegative `w` with `w₂=1`; write `wmax=max w`. Let `W` be the published decreasing-degree marked weight with `d=.19`, `q=.272`, `ρ=.462`, `m=.765931373`. The replay reconstructs all its degree pairs up to 20, checks the intrinsic gain and promotion inequalities, and uses the saturated tail. In its recursive tree form it is

```
W(t) = [h(D,c) + Σ_y a_y W(y)] / (1/i+d),   c=1+#children.
```

Thus it retains heterogeneous pre-existing arms and grandchildren. Off the virgin it obeys `(s−R)W≥m+(s−ρ)W`. Since `W≥m/(.5+d)`, one shared scale is

```
B_s(w,t) = wmax W(t) / [m (1+(s−ρ)/(.5+d))].
```

Its drift pays the weighted nontrigger toll `b w_min(D,20)≤wmax`, and its value pays `w_i` at every virgin hit. Stopped Dynkin gives the same pointwise upper on every truncated excursion action for remaining time `h≤280`, and on the infinite excursion. In numerical comparison the scale multiplying `wmax W` is `1.134984448895`, versus `1.305599999232` from `1/m`. The exact heterogeneous fixture prints the envelope `3.8331803531` per unit `wmax`. The omitted one-row mass after 280 from a virgin is bounded by `wmax·1.152731590e−12`, using this particular `.462` cone. The existing `.448` duration proof is sharper on the extreme-time tail; this new scale is a comparison to the local `.462` exit envelope, not an improvement to the best aggregate `6.202e−12` truncation error. The short action itself still needs a lower and tighter upper evaluation. The high-root seam (`round2/sol-moments/positive-comparison/seam/README.md`, historical reference) can use this shared anchor, but does not price low-root/high-descendant exits by itself.

The filtered two-copy table in [filtered](../filtered/README.md) measures a different absolute observable (trigger reward with a specified total copy count, including direct resets). Adding its upper as if it measured the first-virgin exit would double count some continuations and omit many other ones. No such composition is used here.

## Numerical source bound and what remains

Separately, use the audited post-EDGE prior `[.565678056632067,.565770469400619]`. Log convexity from the certified `r(K_.4481)<4.062` and the upper exponent prior gives the trial `r(K_s)≤R_s` computed by an outward rational exponential. The inherited all-label comparison gives a lower `g₂` at this same Perron multiplier, and the immediate trigger floor removes at least `g₂/[2(s+.05+1.5)]` uniformly. Thus, as in the positive comparison, `F(s)≤1+(R_s−1)/(1−θ_s)`. This upper is evaluated directly at the target (not copied from the lower trial). For the lower, at an earlier anchor `t` with `F(t)≤B_t`, monotonicity and convex secant extrapolation from the true root imply

```
F(s) ≥ 1 − (s−L)(B_t−1)/(L−t)
```

if `s≥μ`, and otherwise `F(s)≥1`. The best of the replayed anchors is `.565`. This produces the interval above; the lower is expressly **prior-assisted**, not an independent short-kernel sign. It narrows the previous coarse interval at this trial, but cannot deliver a new exponent endpoint. To use the pointwise route for a new sign one still needs a mark-sensitive approximate short action and a bound on its **signed** residual or adjoint pairing, including the cut increment exposed by the D18 counterexample.

## Replay

```sh
taskset -c 8 nice -n 17 timeout 85s env OPENBLAS_NUM_THREADS=1 ./round2/.venv/bin/python -Bu round2/occupation-neon/pressure-memory/old-motifs/full/nonvirgin/multiple/envelope/replay.py b3f18_lo.bin b3f18_hi.bin
```

The replay passed the exact cut formula on the logged test rows, isolated the same-geometry mark increment, reconstructed the weight, and produced `[.9691216940922066,1.0189758871875265]` before outward display rounding. It uses the pinned original payload layout checks and inherited kernel anchor; it does not rerun GMP or the full threshold-280 integral.
