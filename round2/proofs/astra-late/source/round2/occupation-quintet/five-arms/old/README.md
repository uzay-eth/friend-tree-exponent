> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# Gold1 B5: strictly registered exterior child, live through occupation

The low receiver/donor degree mask pays **at least `2.83447e−9 / 1.90388e−9` potential units/core** (lower/upper witnesses). Its rigorous ceilings are **`<1.270e−8 / <1.071e−8`**, so this scope cannot reach `1e−6/core`. The strongest receiver in the strict mask is the short tip. Giving its sister tip a separate `D=38…55` window did not recover the signal: the two ceilings for that variant are below `3.979e−11 / 5.023e−11`.

## Reservation and scope

Use the original Gold1 physical B5 row, production trial `λ` and outward sign from [the seam replay](../README.md), rewarding only `max(0,−row)`. Root a nonfounder at hub `0`, short tips `1,2`, long-arm interior `3`, far tip `4`. Track one **actual** exterior child `5` of the receiver `x∈{0,1,3,4}`. The six tracked vertices have exactly the children in the augmented tree; their only remaining neighbor is the hub's actual parent. The tagged donor is a core leaf, so additional children of the donor exit the mask too. Post-promotion `c` is the actual core arity, including the tracked donor in its receiver's census. At payment require supported B5, receiver and donor degree `≤5`, and all other degrees `≤8`. The older-sister variant uses receiver `1`, degree caps `(7,5,55,7,5,5)`, and requires `D₂≥38` at payment.

All chronological orders consistent with ancestry are admitted. For a short-tip receiver, call it `1` regardless of whether it or its sister is older. Otherwise require short tip `1` to be born before short tip `2`; this removes just the symmetric duplicate. Roles, extensions and rewards are attached to physical vertex sets, not labels. A birth can legitimately seed multiple prefixes for different role assignments, and the old prefix ceases to be strict. For a paid row the unique exterior child and unique role assignment prevent duplicate credit, including across receiver masks. Every paid row has an exterior child, which excludes the earlier strict B5 reservation; the older-sister variant is disjoint from the low-degree union by degree. 835c confirmed its current H3 work does not pay B5 rows.

The strict scope excludes alternative root orientations, extra participant children and older tips outside the explicit variant. Also, the measured B5 **uncovered donor-gap** `0.978/0.741e−6/core` in the pooled receiver/donor D=2–5 bin is an allocated diagnostic component, not a measurement of this mask's negative-row reward. The older-tip diagnostic was for all B5 orientations and included much wider degree bands. No inference about their total occupation follows from the ceilings here.

## Exact live-donor comparison

At each chronological prefix state, a physical newborn starts at `(2,1)` and its parent's degree is unchanged; that parent's arity increases by one. Leaf stock is `D−c`; promotion occurs at `(D−c)/D`. Every tracked core gains at the leaf-stock rate plus the sum of actual reciprocal degrees of its tracked neighbors. In particular the tagged donor contributes **the same current `1/D₅` to the B5 receiver's gain hazard and its signed outward row**. Donor gains update that reciprocal through a real `D₅→D₅+1` transition; there is no birth mark, independent donor draw, or reward equal to the difference from a pointwise box bound. `g`, FY, and newborn motif values are not paid here.

The hub's own parent is the only untracked contact. Its *same live reciprocal* appears in the gain hazard and reward, with the Bellman inequality checked uniformly over `[0,1/2]`. This comparison is valid without a law for the parent, rather than assigning a box-gap reward to it. The exact signed row comes from `Top.row`, not a projected center/radius bound. For terminal states the sign breakpoint, if interior, is checked as well as endpoints. With each substep increasing degree or prefix size, exact directed rounding on a `2^20` fractional-tick grid verifies the triangular comparison. Degree cap exits and extra births are irreversible; failure of complete internal support cannot reverse. The older-tip lower threshold gives zero **reward**, not zero continuation, when it has yet to be crossed. Summation over actual birth registrations and the Gold1 discount bracket gives `LO·V_floor(2,1) ≤ π[R] ≤ HI·V_ceiling(2,1)`; see [the preceding comparison](../README.md).

## Replayed results

| Mask and witness | Floor (ticks/core) | Ceiling (ticks/core) |
|---|---:|---:|
| low-degree union, lower | `3116.54214226939` | `13960.9923844050` |
| low-degree union, upper | `2093.34494551223` | `11766.1682247918` |
| older sister, lower | `3.42633517725443` | `43.73878277682451` |
| older sister, upper | not run | `55.22610257464006` |

Exact low-degree union floor fractions (lower/upper) are `2040053084480070672876140613/654588640663961600000000` and `5481119289294508924951213593/2618354562655846400000000`; exact ceiling fractions are `2272536175610329161519138181/162777552844226560000000` and `4788170174962767378218601341/406943882110566400000000` ticks/core.

Run from `.` with `R=round2/occupation-quintet/five-arms/old/live.py`, on CPU 13, nice 19, wall deadlines. Every listed completed job exited zero. The CLI emits exact rationals and the counts of terminal states and checked inequalities.

- `python3 -Bu "$R" --side {lower,upper} --fixture` (separate invocations): each passed 15 exact physical row checks, 15 live degree-change checks, and 100 chronology checks; walls `3.08/2.13s`. Literal census uses explicit external child and parent and individually computes gain, promotion and incoming deleted-tip transfer terms.
- `python3 -Bu "$R" --side {lower,upper}`: low-degree floors, walls `51.91/64.06s`.
- `python3 -Bu "$R" --side {lower,upper} --ceiling`: low-degree ceilings, walls `57.60/46.19s`; `9,604/25,088/20,580/25,088` terminal states per receiver `0/1/3/4` on each witness.
- `python3 -Bu "$R" --side lower --receiver 1 --internal-cap 7 --tip-cap 5 --old-role 2 --old-cap 55 --tip-min 38`: older-sister floor, wall `31.26s`.
- Same older-sister command with `--ceiling` on lower and upper: walls `45.47/32.72s`, each `82,944` terminal states, `165,888` terminal corners. Exact ceilings `8899615030995871703871299/203471941055283200000000` and `22473924575560371455085827/406943882110566400000000` ticks/core.
