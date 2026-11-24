> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Unused first-body slack of repaired 80fe026

**Authenticated outward interval: `[0.565676955958, 0.565771614545]`.**
The previous corrected release was `[0.565676930574, 0.565771634792]`;
the new endpoints move inward by `2.5384e−8` and `2.0247e−8`.

## What is left after the correction

Fix the same `s, ε, δ=ε(s−λ), M, t, z` and corrected `Φnew` as in
80fe026. Write `Q_x=q(x)+t(k−h)` for the all-eight *potential* charge,
so `Φ≤ΣQ_x≤MC` and `Σ t(k−h)=0`. Let `F(x,ρ)` be the original first-body
*drift* charge at `λ`, and `A_ρ` the two physical census transitions.
For nonfounders the unused row is

```
R(x,ρ) = F(x,ρ) + δ(M−Q_x) + (A_ρ−s)z(x) ≥ 0.
```

The row identity uses the same `ρ` in the physical gain and in `F`.
With `V=z−z(J)` and `m=s z(J)−δM`, the corrected signed drift satisfies

```
ε(ℒ−s)Φnew ≥ m(C−1) + Σ_nonfounder R(x,ρ) + Σ old non-f rows.
```

In detail, the shift contributes `−δΦ`, charged conservatively to
`−δΣQ_x`; replace the founder charge by `−δM` and pay it with the
signed constant in `Φnew`. The nonfounder terms have the exact expansion
`F−δ Q_x+(A_ρ−s)V=m+R`. The unchanged newborn correction is zero.
The unspent nonnegative difference `δ(ΣQ_x−Φ)` is not booked. All seven
non-`f` charges stay at the original `λ`; **no original `f` occupation
floor is added**. The base margin `m` is counted once, as it was in the
previous release. The new denominator remains the same `πΦnew` ceiling.

Through degree 64 the replay re-extracts the same authenticated affine
rows and checks both exact `ρ` corners of `R`; the smaller corner is at
the scale of directed Bellman rounding. For degree above 64 it uses the
verified all-stock floors and the same arity cap: gain cancels exactly
there and the promotion corner inequalities certify nonnegativity even
in the saturated overflow arity. The priced objective retains only
`D≤18`; all excess elsewhere is safely dropped.

## Matched pricing and comparison

The LP keeps exact signed charges as rationals and certifies rounded
floating dual proposals by checking every column, subtracting the worst
positive reduced-cost defect per mass class. The two LPs have exactly the
same 7,640 columns and objectives; switching memory on adds 17 balance
rows only. One floating upper-side proposal failed; the fallback proposal
passed the same exact rational check. Prices below are approximate
displays; the exact directed fractions are in [`full.log`](full.log).

| Unused row expectation | Lower witness | Upper witness |
|---|---:|---:|
| Direct physical census lower bound | 1.806e−12 | 1.959e−12 |
| Endpoint-flow lower bound | 3.1583e−8 | 2.5247e−8 |
| Endpoint flow plus receiver-memory lower bound | 3.4442e−8 | 2.7681e−8 |

For comparison, repricing original `f` with **the same mask** and same
endpoint-flow constraints, then using the same seven non-`f` rewards and
the original witness denominators, yields
`[0.565676491598,0.565771989808]`. This is a comparator, not a term in the
corrected result. The repaired release and the new computation transfer
the original seven disjoint rewards unchanged.

For the memory lift, put `q=1/D_y`, `w=q r_y`; the source contribution is
`w q/(1+q)`. An exact chord and tangent enclose the latter rational
function, including the `q=0` tail. The common oriented edge table keeps
both endpoint reciprocal moments and the child's gain-normalized donor
moment. Quadratic receiver terms use rational Bézier triangles enclosing
`(B,B²)`. Saturated source arities restore excess stubs and aggregate
upward pressure separately, and tail gain-normalized rates are enclosed
by `[b,(1+b)/2]`. All these controls are present in both LPs. See
`pressure-memory/marginal/replay.py` (`../../../../../../pressure-memory/marginal/replay.py`, historical reference)
for their independent matched first-body/flux pilot.

## One-command replay

From the repository root, using the two authenticated production payloads:

```sh
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 15 /usr/bin/time -f 'wall=%e cpu=%U+%S exit=%x' round2/.work/venv/bin/python -Bu round2/occupation-neon/interpolate/compose/residual/replay.py --lower b3f18_lo.bin --upper b3f18_hi.bin --cut 18
```

The full run regenerates the baseline composed ledger, authenticates both
binary hashes and the prior receipts, checks the repaired and planted
physical regressions, then regenerates the residual LP prices. It passed
in `155.80 s` elapsed, `150.92+2.41` CPU seconds, exit `0`. This is an
incremental certified gain, not an eight-decimal enclosure.
