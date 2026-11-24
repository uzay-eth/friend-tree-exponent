# Cold audit: original Gold2 lower singleton, conditioned through D430

**Holds** for the original singleton lower price **207990.215716345 ticks/core**, rounded downward. Reviewed proof (`occupation-quintet-gold2/lifetime/README.md`) `cc6c8b3c3c9c5d52917e41b17c18fad2018bb57e`, checker (`occupation-quintet-gold2/lifetime/check.py`) `4a62abe3065ff0d08d1eb5d229a1c5e9e996799d` and conditioned kernel (`occupation-quintet/lifetime/conditioned/replay.py`) `8d647aacd635a01bbf0decacebdd18c62078e69d`. The existing Gold2 payload authentication, original pointwise caps and stable-age occupation law are inputs. Both replays mmap the lower payload read-only; this audit does not rehash the private payload. This verifies the singleton input to a restoration, not the sibling H3 price or the combined denominator claim. No pressure LP or withdrawn H3 payment enters.

## Signs, physical rates, and tick units

Let `L=D−c`, `ρ=Σ 1/D_neighbor` over the core slots, `f_0=f(D,c)`, `Δ=f(D+1,c)−f_0` and `J=f(D,c+1)−f_0+f(2,1)+g(postmother,(2,1))`. In the lower witness `σ=−1`. In original potential ticks (`2^40` ticks per unit), the singleton running price in an uncentered row is exactly

```
55 + (L+ρ)Δ + (L/D)J − λ f_0.
```

Thus the reward slope is `+Δ` per unit of reciprocal receive, not `−Δ`, not `Δ/D`, and not a second factor of the total receive. The gain continuation uses physical rate `L+ρ` in **both** the numerator and holding denominator. On `D<252,c≤22` the charge instead substitutes `c·stored_center/2^30` for `ρ`, so the running slope is zero there while the **gain rate retains physical `ρ`**. This is the original centered singleton account; the other original families handle the centering transfer. At `D=251,c=23` one is already outside the centered class; at `D=252` all rows are uncentered. At `D=3999`, the pre-jump singleton is tabled and its gain successor at `4000` is analytic. The promotion uses the **post-promotion arity** with the original overflow convention and the proper linear tail interpolation across the two edge labels. The new core's `f+g` is purchased once in the running reward at `L/D`; promotion also changes the *old* focal state and changes the latest-child coordinates to `(2,1)`. Those are different objects and neither is an additional source credit.

The broad lower quotient adds `ρ·(Δ+V_gain)` to its numerator relative to receive zero, and `ρ` to its denominator. In the numerical recurrence both the reward intercept and reward slope are floored to `10^-12` of a potential tick. Floored continuation values are also in lifetime ticks on this auxiliary grid. Since all physical rate coefficients are nonnegative, flooring the slope is safe regardless of its sign. The common integer denominator in the code clears time rates, not a change of units. In the conditioned row, the known mother/latest-child reciprocal sum enters the reward with coefficient `z` and the gain with coefficient `W_gain`; the unresolved focal port adds the same joint coefficient `z+W_gain`. The other two uncertain channels are **mother gain** and **child gain**; their ports use only the corresponding continuation value, not the focal reward slope. Latest-child leaf promotion uses its separate fixed rate. The fractional-linear minimum over a box can be sorted by each port's successor-plus-reward slope: starting from all lower rates, include a port while its slope is below the current quotient, then no later port can improve it. Endpoints include the known physical rate, and boxed rates are permissible even when the actual rates are dependent. The max with the broad lower at a row is valid as a choice of two available physical lower continuations.

## Truncations and stationary price

The conditional row includes `D≤430,c≤18`. At its degree or arity exit the handoff goes to the broad value at the reached **true** census, not the state before the boundary jump. A live mother is dropped into the `≥3` bucket only upon her true gain; a latest child is reset only on an actual focal promotion. Tracking does not remove the old donor's contribution; after loss of precision its reciprocal remains in the proper unresolved box. The broad solve truncates at `D=7999,c=30`; exiting the arity mask uses zero, safe for the nonnegative original running price. At `D=8000`, the lower continuation is `h_c/D`. The suffix proof uses the same uncentered price, with the analytic residual polynomial checked at both receive endpoints. Dividing by `1/D` uses `P(0)≥0`, replacing `P(0)D` by its minimum at the cut before applying the Bernstein bound. The relative decrease in `h_c/D` on a gain contributes at most one per unit time, while the promotion rate is in `[1−c/8000,1]`; the suffix's two endpoint tests therefore give a lower discounted potential. The label count and arity bounds used in these comparisons are physical, not a forced evolution of the actual outside tree.

At birth a nonfounder is truly `(2,1)`. Its mother is either exactly degree two or at least three at that birth; the daughter has no earlier child. A nonfounder has at most one degree-two maternity event, while in its first virgin state. There the desired rate is `1/2` and the competing gain rate lies in `[1,3/2]`; discount at the true exponent and bound with the original Gold2 endpoints. The stationary birth measure (equivalently discounted birth compensation, ignoring founder entrance) then gives the fraction range used to mix the two roots. This does not use an independent physical clock under a sampled future birth. For the lower lifetime the discount is the original upper endpoint, and the stationary price multiplies by the original lower birth intensity **once**. The mixed lifetime before this factor is `106709227048431615194065682065030463/290238887678415021500000000000` ticks. Its outward product is the checked price; it is not already net of `19217` ticks/core. For a combined repair, the original per-core error must be subtracted only once from the sum of separately validated disjoint original rewards, before division by the original denominator. The singleton alone yields the old lower endpoint reported by the replay, below `99/175`.

## Independent seams and real output

The cold checker reconstructs the charge directly from pre/post singleton values and the original newborn edge labels. It compares 88 distinct rows at three receives, finding 73 nonzero slopes and 15 zero slopes. It then evaluates 75 unrounded full broad quotients across the layout and analytic degree seams and the degree-8000 handoff, including arity-30 exits. It instruments only an in-memory copy of the conditioned loop to observe the candidate **before** the broad fallback, and compares five direct physical quotients on compatible old-donor completions. The degree-430, arity-18 rows exit to broad values on either focal gain or promotion; the live-child fixture also crosses to the degree-4 child bucket. The two boundary rows share actual receive between reward and gain, with nonzero slope. The independent row check does not claim that arbitrary intermediate tables or the private finite cap were reproduced from scratch.

```text
taskset -c 8 nice -n 17 timeout 120s python3 -Bu round2/occupation-quintet-gold2/lifetime/check.py --side lower --conditioned 430
lower conditioned 430 states/corners 86742 750 birth roots (Fraction(225111373432610969, 500000000000), Fraction(173838070999720941, 500000000000)) mixed 106709227048431615194065682065030463/290238887678415021500000000000
lower lower price n 239534 V ticks 106709227048431615194065682065030463/290238887678415021500000000000 reward ticks 207990.215716345 den 1808603193693 endpoint 0.565713784375142 inward 0.000000298041786

taskset -c 8 nice -n 17 timeout 120s python3 -IBu round2/checks/gold2-singleton-repair-sol.py
PASS direct charge vs writer intercept/slope at seams: 73 0 15
PASS 75 full broad quotients at layout, analytic and terminal seams
PASS five direct physical Bellman quotients; d,c,mother,child,reward-slope:
430 18 -3 (4, 0) 67383747
430 18 2 (3, 1) 67383747
430 1 2 (0, 0) 5117669
252 18 2 (3, 2) 201001239
251 18 2 (3, 1) 0
PASS exact conditioned mixed value and outward 207990.215716345 ticks/core 86742 750
```
