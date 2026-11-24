> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../../README.md#claim-25) for exact scope and citations.

# Conditioned quintet reward on the common sided census

The frozen `fd63517` D360/c18 lower continuation is exported for combination as a **reduced score**, without modifying the freeze. The projection is positive in every retained cell in the sense that its affine score is nonnegative and not identically zero. It is **not a material new occupation payment**: its maximum over the full census boxes is below `2.708e-9 / 6.517e-9` **ticks per core** (lower/upper witness). The one-dimensional source marginal remains useful for joint pricing.

## Contract

[`lower.json`](lower.json) and [`upper.json`](upper.json) use `quintet-conditioned-f-reduced-v1`. Every row is `[D,k,A,B,C]`, with exact `Fraction` strings and objective `A+B*u+C*v` in original ticks. Here `k=c-1` counts real children, `u` is total receive from those children, and `v` is receive from the mother. Use the **same** `u,v` as for the gain transition and the edge-flux marginals. Unlisted cells contribute zero reduced score. All 6,326 cells with `2<=D<=360` and `1<=c<=min(D,18)` are present.

There are two ways to integrate this export:

1. Add `birth_floor_ticks` to a certified stationary lower bound on the occupation of the exported score. `birth_floor_ticks` is the exact old conditioned floor; it already includes the physical birth factor.
2. In a joint stationary objective, include `b*birth_V_ticks[D=2 or D>=3]` in the **promoting mother's** census, with `b=(D-k-1)/D`, in addition to the reduced occupation score. This source term applies at **all** maternal degrees, including outside the residual mask and on overflow. For an `eta=pi/mu` normalization, both terms are in the same units and the resulting nonnegative lower bound is multiplied by `LO` once. Do not add the separate source floor in this version.

For the optional `beta_D2_bounds` cut in an `eta` LP, bound `sum eta_column * b * 1{D=2}` between the two stored endpoints, since `eta(b)=1`. In `pi` units the same cut is relative to `pi(b)=mu`. The ambiguity of founder-mother births vanishes in normalized occupation, since those occur at mean rate at most one from the unique founder. Objective `A+B*u+C*v` is an occupation score and takes **no additional** factor of `b`.

## Exact construction

`export.py` hash-pins the frozen conditioned replay and adds only an in-memory observer after each stored value. It runs one original lower sweep per witness, with the original unconditioned analytic continuation. No other solver is rerun. The captured virgin values (grid `10^12`) are exactly those of the frozen proof:

| Witness | Mother degree 2 | Mother degree ≥3 |
|---|---:|---:|
| Lower | `1322979352773609470` | `1020282677750521357` |
| Upper | `1755696543328942415` | `1457540240862028342` |

For each hidden state, recompute `F+(A_phys-HI)V` from the unrounded original row and exact stored successor values. The mother/child **own gain** rate boxes can be minimized on their corresponding `ΔV` before projecting. The focal gain and original reward both retain the same actual `rho=u+v`; neither minimizes it separately. Hidden states supply the following admissible boxes in the external common coordinates:

- Mother exact degree 2: `v=1/2`. Mother overflow: `v in [0,1/3]`.
- Zero children: `u=0`. If the latest child's exact degree is `j`, `u in [1/j,1/j+(k-1)/2]`. For child degree ≥4, `u in [0,1/4+(k-1)/2]`.

The residual after hidden rate minimization is affine on each such box. All hidden corners are checked nonnegative first. For each common census, take the minimum of all constraints from the same `(u,v)` corner, then select an affine minorant that is also nonnegative on the entire common box `[0,k/2]×[0,1/2]`. Selection is by box-center value among a finite menu of the constant floor and scaled ramps from box corners; this is a reproducible valid minorant, not a claim of the optimal possible affine minorant. Each selected line is verified against **every applicable hidden corner**, hence the entire hidden boxes by affinity. The broad fallback outside the retained mask has nonnegative reduced reward by its original Bellman subsolution and monotonicity of degree/arity. Exact initial values, the discounted probability interval for physical degree-2 mothers, and the exact previous conditioned source floor are checked in the export.

The stationary balance is `pi F = pi(b*V_birth(D)) + (HI-mu)*pi V + pi R_full`, with `V_birth(D)` selected by the actual mother's pre-promotion degree. Promotion leaves her degree unchanged. The newest child's observable clock starts only on that real event, current child replacement is in the old-vertex generator, and the focal newborn source is inserted only once. Since `V>=0`, dropping the discount correction and replacing `R_full` with the exported minorant gives the claimed joint lower objective. For the degree-2 source cut the underlying newborn-at-`(2,1)` has birth hazard `1/2` and competing gain in `[1,3/2]`, so compensation gives the two stored probability endpoints `(1/2)/(HI+2)` and `(1/2)/(LO+3/2)`.

## Evidence

From `.`:

```sh
taskset -c 8 nice -n 17 timeout 160s python3 -Bu round2/occupation-quintet/lifetime/conditioned/reduced/export.py --replay
```

The exporter authenticates both original payloads. It captured 72,322 frozen values and checked 215,889 hidden residual corners per witness. Of 6,326 exported rows, 6,092 lower / 6,233 upper have a nonzero receive coefficient. Export SHA-256:

```text
d9afcb58793f49c6e5c13d8e60c0e7a89ccd3d70b0a973a47d7eb373305b97c0  lower.json
c4cfe3ad3b4a6ab49be79831994742701c2c4bb5e6b4f61009d056b72f3808fa  upper.json
```

The near-zero projected residual is a measured limitation of the export: free mixing over the hidden latest-child census still admits near-Bellman-tight cases. No joint solver result or arbitrary hidden-source marginal is claimed.
