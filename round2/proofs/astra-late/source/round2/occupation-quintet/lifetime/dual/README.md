> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# Original quintet singleton Bellman dual

The frozen `D=12000,c<=30` Bellman subsolution is exported exactly for both witnesses. The sided birth-law LP at degree 24 certifies additional stationary singleton slack of at least **93,573.157772619572338 / 170,995.906971746003733 ticks** beyond the original `LO*V(J)`. This produces total singleton floors **391,079.994162541737385 / 731,984.671964200118015 ticks**, rounded down. No new denominator or higher-host cap is used here.

## Integration contract

[`lower.json`](lower.json) and [`upper.json`](upper.json) use schema `quintet-f-dual-v1`. Metadata pins the source commit, receipt and production payload hashes, the original tick scale, signed `cap_f_ticks`, the nominal `lambda` and exact original quintet `LO` and `HI`. `degree=12000`, `arity=30`, `cut=64` describe the solve and exported mask.

- `values`: `[D,c,integer_V_subticks]`, every `2<=D<=65` and `1<=c<=min(D,30)`. Divide the third field by `V_grid=10^12` to get ticks. The retained state uses **physical core degree** `c`, including the mother. `V_J_ticks` is the exact rational value at `(2,1)`.
- `rows`: `[D,c,A,B]` for the same arity range through `D=64`. Exact rational strings `A,B` represent the reduced reward `A+B*rho` in ticks. Use `Fraction` on the strings. Missing rows have **zero masked reward**. Missing rows outside the export are not assertions that the global Bellman value or reduced reward vanishes there. The overflow continuation for `c>30` is zero by construction.
- Each sided occupation column `(D,k,b,p,u,v)` uses `c=k+1` and objective `A+B*(u+v)` when listed. `u` is the actual summed child receive and `v` the actual maternal receive; these are the same column coordinates in the gain `D-k-1+u+v` and the reciprocal edge marginals. Do not independently minimize the objective over `rho` while using a different receive for the transitions. Use zero cost for overflow/nonoccupation columns. This is a score per unit occupation, **not** a birth score; there is no `b` factor in the objective.

## Identity and verification

Let `F` be original signed `f_slack` in ticks (including `cap_f=0/-1053`). Let `A_phys` be the focal old-vertex generator, without newborn source, and write `R=F+(A_phys-HI)V`. Its intercept uses `L*(V_gain-V)+(L/D)*(V_bud-V)-HI*V`; its receive slope is the original reward slope plus `V_gain-V`. The gain successor at `D=64` is the actual value at `D=65`; the bud rate vanishes at `c=D`. The bounded original subsolution has `R>=0` globally. The export reconstructs the original sweep with the same analytic continuation and directed grid rounding, verifies every integer Bellman inequality, and independently checks exact affine nonnegativity against the unrounded rows at both receive corners in the exported slice. `V(J)` exactly matches the published values `262953589366589949/500000000000` and `991668031155825523/1000000000000` ticks.

For the physical stable-age law, discounted Dynkin from birth, or stationarity of the sum over nonfounders including the newborn source, yields `pi(A_phys V)=mu*pi V-mu*V(J)`. The bounded value and bounded jump drift suffice: all tabular degrees are finite, and `V=h[c]/D` beyond the degree cut has gain drift bounded using rate `<=D`, plus promotion rate at most one. Hence

```
pi F = mu*V(J) + (HI-mu)*pi V + pi R
     >= LO*V(J) + pi R_masked.
```

This is a strengthening of the old singleton payment, not a separate reuse of the same reward. Founder omission is inherited. The sided LP works in units `eta=pi/mu`, so its certified nonnegative lower bound is multiplied by `LO` once. The native exact repair checks every rational reduced cost with the updated Quintet `LO` and `HI`; HiGHS only supplies candidate duals. The inherited birth-class kernel is a physical lower bound that remains valid because its construction uses an upper discount at least the exact quintet `HI`. No D18 numerical rewards or denominators are used.

## Evidence

From `.`, run:

```sh
taskset -c 8 nice -n 17 python3 -Bu round2/occupation-quintet/lifetime/dual/export.py --replay
taskset -c 8 nice -n 17 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 ./round2/.venv/bin/python -Bu round2/occupation-quintet/lifetime/dual/sided.py --limit 24
```

The export authenticated both original blobs and verified **1,514 values, 1,484 affine rows and 2,968 exact corner tests per side**. Minimum corners were `0` lower and `1238687088504128371/19016504500787500000000000000000` upper. Export SHA-256:

```text
d698df2bdd20fd6d5a802f84d18124a1fcd6f491bb4dde8ae5038ca5d27f6169  lower.json
a23fe15ac518c35f9e7dcc92bc1fa497253f16e7677a9d94e572815682253fa8  upper.json
```

The degree-24 LP passed with 613 rows and 5,554 columns. The exact additional stationary rewards, before rounding the display, are:

```text
1276783006042809110549818991696409322453507199781/13644757069600662706896800000000000000000000
212108873684137139278951921799899102453819851509/1240432460872787518808800000000000000000000
```

The degree-64 LP built 2,393 rows and 19,634 columns but timed out at 140 seconds; there is no degree-64 LP price in this result. All degree-64 affine rows are independently verified and available for integration. The original lifetime replay file was left unchanged.
