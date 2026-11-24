> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../README.md#claim-21) for exact scope and citations.

# Matched edge correction, composed with the repaired D18 release

Against frozen commit `80fe02660b9f0a0d4e7fdc8240b0815615ff5896`, exact rational replay gives **`[0.565677980686114, 0.565770513279953]`** outward. This is an inward move of `1.050112114e−6` lower and `1.121512047e−6` upper compared with the frozen outward `[0.565676930574, 0.565771634792]`. This is a composition of a pointwise residual inequality with the frozen stationary reserve; the two new trial points do not themselves have positive uniform pointwise margins.

| | Lower | Upper |
|---|---:|---:|
| Trial `s` | `.565677930574` | `.565770634792` |
| Edge cutoff (each endpoint `D,k`) | `40,8` | `32,7` |
| Exact Bellman states | 91,732 | 44,492 |
| `f`-only margin at same `s` | `−3.753084e−6` | `−3.585438e−6` |
| Matched signed margin `m` | `−2.350525e−6` | `−2.081259e−6` |
| Total correction absolute norm/core, rounded up | `.000676630` | `.000200471` |
| Edge value sup, rounded up | `.000648892` | `.000176328` |

## Local adjoint and matched birth

Each counted edge goes from a **nonfounder** `p` to its actual core child `w`. The edge state is `(D_p,k_p,D_w,k_w)`. In the physical generator on an existing edge, the two gains have rates `L_p+1/D_w+u`, `L_w+1/D_p+v`, where `u∈[0,k_p/2]` includes the parent's own parent and the other children, and `v∈[0,k_w/2]` includes the child's children. Promotions at either end have rates `L_i/D_i`. The two external sums may vary jointly and adaptively. The existing edge survives each promotion and takes that end's post-promotion descriptor.

Write the available production `g` charge as `a+b u+c v` after the reserve subtraction below. Solve `(s−A_edge)V ≤ max(0,a+b u+c v)` with zero continuation beyond the cut. The same `u,v` enter the reward, transitions and holding denominator. The physical charge is nonnegative after reserving; the hinge is only for a relaxed box. Verify the rounded inequality at the four corners and any hinge intersections. Cross multiplication makes the comparison affine on each polygon. The bounded monotone census DAG permits exact descending evaluation. On the unretained part `V=0` and no edge can enter the retained part starting from a cut violation. The upward-rounded recurrence also encloses the largest robust margin in this fixed box/cut class; it does not bound actual physical residuals of other correctors.

On a promotion at parent `(D,k)` add the exact new edge value `V((D,k+1),(2,0))` in the parent's census Bellman right-hand side, at rate `L/D`. The full signed correction is `Σ_nonfounder(z−z(J)) + Σ_edges V − δM/s`; the physical newborn census correction is zero. Edge births from the founder are excluded on both sides. A promotion changes all old incident edge values as well as creating the new one, as the exact whole-tree difference check verifies. The global constant covers the founder's original denominator cap. Census comparisons use the repaired five stock floors and both exact gain-rate endpoints, with arity-only continuation above degree 64. Norm bounds are much smaller than the original strict pointwise positivity constant.

## Reserve and composition

The original `λ` and charge decomposition are kept. Spend charged `f` **only** in its original census equation. The extra corrector spends `g` **above** two pointwise reserves: the frozen uniform `m_inner` on every both-inner edge and the frozen mixed virgin `g` floor for parent degree 18–191 and arity ≤10. The replay subtracts these before forming Bellman rows; the high-parent floor on the intersection is independently recomputed from the frozen high-parent extractor and checked exactly against the new extractor. It leaves all six higher-body charge families intact. Consequently the three blocks of the frozen all-seven occupation reserve remain valid without overlap. Pressure-memory has not been combined.

The resulting inequality is `ε(ℒ−s)Φ_new ≥ m(C−1) + Q_reserved(T)` on finite trees. Normalize and integrate as in the frozen release; the founder and global constant vanish. The numerator `m + π Q_reserved` is strictly positive at both final trials, even though `m<0`. Divide its certified lower bound by the old all-degree stationary denominator plus the one-sided uniform correction sup (the upper sign can drop `−Σ_edges V≤0`). This conservative denominator suffices. Frozen reserve components in the code are rounded down from the authenticated `compose/{lower,upper}.log` receipts. At degree 36, arity 8 the lower target still missed by `3.2254e−8` on the robust margin even after granting the frozen ledger's rounded-up numerator; widening to 40 crossed it. This was a cutoff limitation, not a physical impossibility result.

## Replay

From `.`, run the sides serially. No reuse of the withdrawn first-body endpoint is involved.

```sh
ulimit -v 10485760
export OPENBLAS_NUM_THREADS=1
P=./round2/.venv/bin/python
R=round2/occupation-neon/interpolate/edge/replay.py
nice -n 10 taskset -c 7 timeout 160s "$P" -Bu "$R" b3f18_lo.bin --degree 40 --arity 8 --challenge --mask disjoint --s .565677930574
nice -n 10 taskset -c 7 timeout 200s "$P" -Bu "$R" b3f18_hi.bin --degree 32 --arity 7 --challenge --mask disjoint --s .565770634792
```

Passed in 73.60 s and 47.09 s. Each side authenticated the original full payload, matched all 2,079 frozen first-body rows against hash-pinned data, used the five pinned stock floors, checked the Bellman inequalities in rational arithmetic, and passed 119 exact finite-tree generator comparisons with 391 nonzero post-promotion edge values. The smallest whole-tree surplus over the stated signed margin was `5.135075426193229e−6` lower and `4.264437811200185e−6` upper. The seven-core regression of the withdrawn source is in `regression.py`. No independent reviewer has signed off on this new edge composition.
