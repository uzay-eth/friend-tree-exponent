> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../../README.md#claim-21) for exact scope and citations.

# Native high-apex P4 and high-terminal Yb transport

**A two-sided action below `1e−6` per core is verified on the specified P4 support.** At multiplier `θ=1.2415`, the lower-minus-`θ`-upper contrast of the terminal directed-donor excess lies in **`[3.234195491229276e−7, 8.271905907218992e−7]`** (see the exact fractions below). This is one contribution to the signed action, not an enclosure of the entire P4 charge or the global signed action. A larger lower-only run certifies `6.9841289325e−7` per core, short of `1e−6`. No exponent endpoint is booked.

| Objective | Cut `Da ka Db kb Dc kc Dt kt Ddon kdon` | Signed stationary interval/core | Native rows | Wall / peak RSS |
|---|---|---|---:|---:|
| P4 contrast | `48 5 5 2 5 2 7 2 8 2` | `[3.2341954912e−7, 8.2719059073e−7]` | 2,996,822 | 78.61 s / 279,856 KiB |
| P4 lower-witness donor excess | `80 6 6 3 5 2 8 2 10 2` | `[6.9841289325e−7, 1.6374865174e−6]` | 15,375,661 | 191.97 s / 278,460 KiB |
| Yb contrast | `8 4 5 2 5 2 40 2 7 2` | `[−7.0709790046e−12, 2.5480314188e−11]` | 1,914,645 | 55.82 s / 316,352 KiB |

P4 exact first-row endpoints are `161709774561463793/500000000000000000000000` and `3308762362887597/4000000000000000000000`. Bounds include the broad discount bracket `[.565673,.565775]`. Timing includes authenticated payload reads and exact reward extraction (not just the native sweep). Native value slabs use 9.50 MB, 25.18 MB and 47.18 MB respectively. No process was stopped at a row quota.

## What is counted

For each original witness and a terminal external child whose reciprocal is `r`, use its nonnegative directional excess `max(0,−b)/2 + b r`, with `b` the exact charged-row terminal child coefficient. It is the part left by the terminal child contribution to the independent-port affine minimum. The contrast uses the *same* physical donor reciprocal for both witnesses and combines coefficients before optimizing; its reward can be signed. The signed interval is not obtained by subtracting two independently priced positive components.

* **P4:** a nonfounder apex with current `18≤D≤Da`, `2≤k≤ka`, and the three descendant path roles in the stated finite degree/arity boxes. Sum over every physical terminal child donor of every such oriented path. The distinguished donor is born at an actual terminal promotion. Its arity remains capped even after leaving the exact-degree box. Path caps are monotone; donor degree overflow is a live aggregate, not a death.
* **Yb:** the D18 table does not support a high apex for this family. Instead use the short apex child as the terminal (`18≤D≤Dt`), and a distinct long apex arm plus its tip. Follow promotions in the order long arm, long-arm tip, short arm, then short-arm donor. Thus this result covers only Yb pairs with the prescribed actual birth order, not all possible ordered embeddings of the shape. Reward extraction reorders these roles into `Yb`'s production tuple.
* The high-apex mask is disjoint from the finite low-apex `pressure-memory/donor_lift.py` P4 mask, and both transport masks miss the frozen size-12 and `D≤16` higher reserve. They may intersect full higher charges in the staged hierarchy, old-flow, other donor bounds, or earlier high-apex payment: compare on an overlapping objective or reserve this exact directional piece separately. The wider positive run and the signed run are not additive. Mixed/gauge tests are not used.

## Physical recurrence and finite verification

Each stage tracks a particular physical connected prefix. At promotion, the existing prefix successor is always counted, and if the promoter is the next tagged parent, also count the new prefix at its **post-promotion census** with a virgin child. Every completed tagged configuration has one such birth history, so the linear sum over all births preserves actual embedding multiplicity. The two ages are not averaged independently. The gain law keeps all actual internal reciprocals; untracked donor sums range over `[0,n/2]` with known exterior count `n`. The tagged donor send is a distinct live channel into its receiver: the very same rate multiplies the payoff slope and its gain destination. All non-tagged promotions still act on physical `(D,k)`.

At a donor degree gain past `Ddon`, keep `(tail,k)` with outgoing reciprocal `[0,1/(Ddon+1)]`, no visible gains, and promotion rate `[max(0,1−(k+1)/(Ddon+1)),1]`. Thus future donor gains continue and hide only their effect on the bounded send; promotions continue to alter its retained arity. This uniform interval is valid under arbitrary adapted histories, including unbounded degree. Outside the paid arity cap, monotonicity permits zero continuation. We did not certify an unbounded apex or terminal degree payoff; they are explicitly finite payment masks.

The skip-free gain sweep eliminates apex degree slabs in descending `D`, using only the current slab and the next one. Within each slab, states descend in stage length and the other monotone census coordinates, so each listed continuation is computed before use. A bound on a holding quotient is integer-verified by optimizing its cross-multiplied linear residual over all admissible rate corners; every original rational rate is outward enclosed (`10^-12` scale) and value and charge intervals are outward rounded (`10^-18`). Integer cross products use 128-bit arithmetic. Multiply the newborn interval by `[LO,HI]` with endpoint signs respected, using the physical stable-age law. Bounded local reward and monotone cut ensure this applies to the signed difference as well as each positive component.

`check.py` independently enumerates all rational corners and recursive physical prefixes for both graph shapes, including signed tolls, genuine newborn terms and the unbounded donor aggregate. It passed 6,610 P4 and 7,020 Yb reference states; the native enclosure's excess rounding error at the root was `<2.6e−17` in these fixtures. The authenticated source and payload hashes are inherited from the frozen EDGE source. The native check is not an independent cold review of the production operands.

## Replay

Run serially from the Astra checkout:

```sh
ulimit -v 10485760
export OPENBLAS_NUM_THREADS=1
P=./round2/.venv/bin/python
R=round2/occupation-neon/interpolate/edge/higher/transport
B=<scratch-path>
g++ -std=c++17 -O3 -Wall -Wextra -Werror "$R/sweep.cpp" -o "$B"
nice -n 10 taskset -c 7 "$P" -Bu "$R/check.py" "$B"
nice -n 10 taskset -c 7 "$P" -Bu "$R/replay.py" b3f18_lo.bin b3f18_hi.bin --native "$B" --config 48 5 5 2 5 2 7 2 8 2
nice -n 10 taskset -c 7 "$P" -Bu "$R/replay.py" b3f18_lo.bin b3f18_hi.bin --native "$B" --config 80 6 6 3 5 2 8 2 10 2 --objective lower
nice -n 10 taskset -c 7 "$P" -Bu "$R/replay.py" b3f18_lo.bin b3f18_hi.bin --native "$B" --family Yb --config 8 4 5 2 5 2 40 2 7 2
```
