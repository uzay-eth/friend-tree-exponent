> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../README.md#claim-25) for exact scope and citations.

# Gold1 negative H3/P4 occupation

**Result.** Exact all-age floors for two disjoint strict nonfounder chain classes sum to **21,591.6797303 lower / 20,790.0758831 upper ticks per core** (`1.96375e−8 / 1.89085e−8` in potential units). These are genuine negative original charged rows, composable with the parent's nonnegative host-cap reserve. The `1e−6` per-core target was not reached. No exponent endpoint or cap booking is asserted here.

Inputs are the private `quintet_{lower,upper}.bin` copies under `round2/.work/replay-inputs/` and the receipt from its clean detached source `b04cc50c1e391e0420eb7b006ac781d484c68543`. The reader reuses the earlier independent `transfer/check.py` binary layout on read-only mmaps; it does not independently rehash or run the production scan. Capping and arithmetic follow the pinned `PROOF.md` and `verify_quintet.cpp`, not the builder's floating-point optimizer.

## Physical rows and screen

Use outward sign `σ=−1` on lower and `+1` on upper, and keep the side's production trial `λ`. An existing occurrence of higher potential `h` has

```
outward row = σ[-λh + Σ_x (L_x + Σ_internal i_y + Σ_external i_w) Δ^L_x h
                    + Σ_x b_x (Δ^B_x h + newborn-next-motif weight)
                    + incoming exact transfers − outgoing exact transfers].
```

Here `i_y=1/D_y`, `b_x=(D_x−c_x)/D_x`, and `c_x` counts the parent for a nonfounder. The newborn is `(2,1)` and the old vertex uses its **post-promotion** arity but unchanged degree. The incoming terms are exactly one per deleted tip (edge→H3 at its middle or H3→P4 at the new internal vertex). For outgoing terms subtract `(i_w−r1_x)Δ^L_x h` for each exterior neighbor whose addition produces a supported P4/S3 (from H3), or P5/B5 (from P4). Coverage depends on the source internal labels, not the donor degree; covered gains retain `r1_x` per exterior donor. Uncovered gains retain actual reciprocals and are affine in the exterior sum at every role. Pinned C++ §§5–6 give these branches; genuine exact physical rows use the lower-level gain difference instead of the projected center/radius upper bound.

For a production envelope, evaluate `σ(i_w−r1_x)m+|i_w−r1_x|r` on the coarse product law (donor independent of the other endpoint coordinates). `screen.py` tracks the excess over the true incoming transfer; subtracting this excess recovers the physical row with the **same** outward convention. H3 and P4 themselves receive one-coordinate envelopes (`23`, `34p`); two- and three-coordinate downstream controls exercise `34s,45bp,45bs,45s4` without booking these rows. `45p` is also checked. A separate concrete event/deletion oracle agrees on **27 physical H3/P4 rows per side**, including **22 nonzero outgoing transfer ledgers** per side. That oracle enumerates physical vertex sets after each literal gain and promotion; it does not use the row evaluator's extension formulas.

For all degrees `≤6`, seven nonfounder arity strata gave `2,075` states per side. H3 `(2,2,1)` had `46/125` lower and `55/125` upper robustly negative; P4 `(2,2,2,1)` had `316/625` and `339/625`. The largest negative robust magnitudes in those seven strata were only `40.24` and `100.06` ticks (upward rounded). The exact endpoints in this cut had zero envelope excess. The production maximum H3 caps are `225/408`, and the P4 maxima are `42/41` ticks; these are family maxima, not newly inferred production block caps.

The sparse boundary screen illustrates why that low cut cannot answer the real question. Numbers below display outward row **upper** bounds over the unknown exterior reciprocal boxes, before/after the production incoming envelope (evaluated exactly before C++'s per-term and final ceiling operations).

| Side, shape, degree and arity | Physical upper bound (million ticks) | With envelope (million ticks) |
|---|---:|---:|
| lower H3 `D=(48,2,48), c=(2,2,1)` | −246.642 | −206.657 |
| upper H3 `D=(60,2,60), c=(2,2,1)` | −109.818 | −102.542 |
| lower P4 `D=(27,2,2,27), c=(2,2,2,1)` | −181.097 | −164.086 |
| upper P4 `D=(18,2,2,18), c=(2,2,2,1)` | −385.924 | −332.712 |

Large row values need a residence price. For example, high endpoint arity `(16,2,1)` has even larger values in the sparse scan, but this is not a proof of appreciable occurrence mass.

## Registration and charge compatibility

`price.py` pays `R=max(0,−outward row)` on a complete terminal class and zero outside it. In the standard chain cut only the registered children are present at each role: the apex has one exterior parent and a directed chain has terminal `c=(2,…,2,1)`. Prefix role degrees and core census are the actual physical labels. A focal newborn starts the length-one prefix. Every promotion that creates the designated next child seeds one next prefix, at the old vertices' **post-event** states. A promotion also updates or exits any existing prefix according to its own arity mask. No cap-clamping or time truncation is used. Degrees and arities only increase, so a fixed prefix that leaves its mask cannot reenter. Likewise a completed embedding in the terminal mask must have had every prefix in its own mask at each insertion; deleting the latest birth is its unique parent prefix. Fork H3 and the three bent P4 arrival orders are separately selectable and disjoint from the strict chains. The registration oracle checked these orders against connected physical subsets with repeated destinations.

For a prefix value `V`, use the physical leaf-gain rates and promotion rates. At a promotion of the designated next parent, add the value of the longer prefix once, in addition to continuation of the old prefix if retained. At a gain, couple that vertex's exterior sum **in both the reward and the gain hazard**. A Bellman corner inequality at discount `H=HI` is

```
(H + total event rate) V
  ≤ reward + Σ_i (gain rate_i) V_after_gain_i
           + Σ_i (b_i) V_after_promotion_i + b_designated V_new_prefix.
```

Within each half-box cut by `a+b·rho=0`, the reward is affine. Checking all box corners plus edge intersections with that plane verifies the inequality everywhere. Values are rounded down on a `2^20` grid of fractional ticks and every inequality is checked exactly. Reverse degree/arity order and decreasing prefix length make the recursion triangular. The ceiling control reverses the inequality, rounds up, and uses `LO` as discount.

Equivalently, let `U(T)` be the sum of the retained prefix values. Each promotion creates one new length-one prefix independent of which other prefixes it continues. Summing the local inequalities gives the finite-tree identity/inequality `L U ≥ HI U − R_total + V(J) L C`. Hence stable-age balance gives `π[R_total] ≥ μ V(J)+(HI−μ)π[U] ≥ LO V(J)`, with the reverse bound `π[R_total]≤HI V_upper(J)` for the same mask. The finite-tree inequality also avoids assuming a normalized stationary law if used directly in an exponent proof. Prefix/row counts per core and nonzero jump rates per affected prefix are bounded in these finite masks. The founder is excluded as a registration root and creates no length-one value for itself, but each new core vertex has a core parent, even when the promoter is the founder.

For every occurrence with original nonnegative cap `K`, `h_out+max(0,−h_out)≤K`. Thus **all** original nonnegative higher-row charges are still paid by the parent cap ledger; this lane credits only the negative part. H3 is hosted at its middle. P4 ownership is half at each internal vertex, so `max(0,−h_out/2)+max(0,−h_out/2)=max(0,−h_out)`; no extra factor of two enters. There is no overlap between H3 and P4 row rewards, even when their underlying trees and promotions overlap. No singleton/edge reward, positive block-cap reserve, or new denominator is used. In particular the auxiliary `ticks/m_max` display is scale context, not a new endpoint.

## Executed receipts

From Astra root, set `R=round2/occupation-quintet/higher`. Commands were run with `python3 -Bu`; longer commands used `taskset -c 11 nice -n 19 /usr/bin/time -f 'user=%U sys=%S'`.

| `price.py` arguments | Lower ticks/core | Upper ticks/core | Terminal states/side | CPU user+sys (s) |
|---|---:|---:|---:|---:|
| `--caps 65 5 65` | 8795.217915702293 | 10453.159609321352 | 16384 | 20.73+0.16 |
| `--caps 36 5 5 36` | 12796.461814642478 | 10336.916273749122 | 19600 | 33.31+0.42 |
| `--caps 65 5 65 --upper-control` | 205007.1115989722 | 207924.06492124725 | 16384 | 20.83+0.11 |
| `--caps 36 5 5 36 --upper-control` | 43597.09042843258 | 39736.747186092136 | 19600 | 33.30+0.19 |
| `--caps 60 4 60` | 7021.672791997933 | 8515.680286487264 | 10443 | 12.93+0.07 |
| `--caps 48 4 48 --projected` | 677.0750531486226 | 633.9929211729021 | 6627 | 6.68+0.37 |
| `--caps 48 4 48 --extra 1` | 5302.625530123369 | 6350.555555881816 | 43710 | 59.82+7.55 |
| `--shape fork --caps 4 60 60` | 49.988930467066595 | 43.317574346009096 | 6962 | 8.27+0.32 |
| `--shape bent-older --caps 4 4 28 28` | 18.206295089246463 | 21.36503665367688 | 4374 | 6.74+0.04 |

The main H3 lower/upper terminal Bellman checks were `37145/37643`; robust positive terminal states `10196/9535`. P4 checked `40161/39952`, with `17622/18005` robust positive. The two main **ceilings sum to less than 248605/247661 ticks**, or `<2.27e−7` potential units on both sides. Thus refining the same masks cannot reach `1e−6`; that is not an exclusion for other higher masks. Do not add overlapping variants in this table.

Other passing commands: `screen.py` default (`1.90+0.07` CPU s), `screen.py --sparse` (`40.22+2.89`, 26992 selected states/side), `screen.py --literal` (`1.61+0.05`, 27/22 rows/ledgers per side), `screen.py --envelopes` (`0.42+0.06`; `15,15,75,15,75,75,375` checks per side for `23,34p,34s,45p,45bp,45bs,45s4`), and `price.py --check-registration` (`fork=6,bent-older=5,bent-younger=5,bent-early=3,chain3=8,chain4=3`, no duplicate physical embedding). The CLI prints exact rational prices alongside decimal displays.
