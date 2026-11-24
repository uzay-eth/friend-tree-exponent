> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../../../README.md#claim-21) for exact scope and citations.

# Birth-conditioned edge occupation cuts

**Matched result:** In the existing degree-4, order-2 embedding LP with shared degree through order 2 and the same identification mass cuts on both sides, these time cuts narrow the virgin-pair interval by **8.25%**. Two held-out shell sums narrow by **7.35%** and **6.41%**. No exponent or D18 endpoint was computed. The implementation appends inequalities to the existing model; it introduces no new column, no alternate generator, and no duplicate graded model.

## Exact physical birth identity

Write `η=π/µ`. Track a uniform physical newborn `x` and its permanent mother `p` at an independent exponential age of rate `µ`, using the [stable-age identity](../../../../../../../astra/source/round2/degree-law.md#birth-times). For a nonnegative edge census payoff `R`, orient the stationary edge count from the mother, then reverse it to the unique child by physical mass transport. Thus the same objective is

```
η Σ_(x child of p) R(p,x) = E_β ∫₀∞ exp(−µu) R(p(u),x(u)) du.
```

For the newborn start, if the mother had pre-promotion census `t=(D,k)`, her post-promotion census is `(D,k+1)` and the tracked child starts at `(2,0)`. The mass of births with that pre-promotion census and grandmother census `g` is **exactly** `b(t) m(g→t)`, with `b(t)=(D−k−1)/D` and `m` the pre-existing `η` labelled edge columns. This is the same formula `β=π B/µ` on a mark of the promoter's old mother. The corresponding numerator for the unobserved grandmother-degree tail is `b(t)[m(t)−Σ_(D_g≤cap) m(g→t)]`. It is nonnegative by the mother projection. The tracked edge is permanent; later births to the same mother are not re-counts of this edge. Neither the pre-promotion degree nor a birth-class marginal is drawn independently of the exterior.

## Rational time-window prices

Solve a discounted Bellman enclosure on the two **actual endpoint** censuses. Stop on exiting the greatest mother or child degree in the reward support (the degree is monotone and the reward cannot return). At all orders the actual receive rates fall into these intervals:

```
mother: [L_p+1/D_x , L_p+1/D_x+(k_p−1)/2+u_g]
child:  [L_x+1/D_p , L_x+1/D_p+k_x/2]
```

Here `u_g=1/D_g(0)` on a known grandmother birth-degree class, and `1/(cap+1)` on its omitted tail. Her degree never decreases, so the bound holds **at every age**, including when the exterior makes adaptive choices. `u_g=1/2` recovers the unstratified price. Both included-edge inverse degrees change with the actual endpoint gains. Mothers' and children's ordinary promotions occur at their exact `L/D`; they create further core children and decrease leaf stock without changing degree. Exterior changes alone leave the two tracked censuses unchanged. No future survival or independence assumption is made.

For each discount endpoint and each local census, go backwards in decreasing degree and arity. Divide `R+Σ r V(next)` by `s+Σ r`. The extrema over exterior rates occur at the corners of the two intervals since the ratio is separately fractional linear. Upper and lower values are rounded outward and rechecked with exact rationals. These inequalities apply to the **true generator on the whole adaptive tree**, not a claim that the projected censuses are Markov. Stopped discounted Dynkin proves the two integral bounds, with discounts `low` for upper and `high` for lower since `R≥0`, and zero after exit. Multiplying these prices by the exact birth source gives two linear LP rows per reward. In particular for the `(2,1)→(2,0)` residence with grandmother degree `d`, the price before multiplication by the birth hazard is bracketed by `1/(high+5/2+1/d)` and `1/(low+5/2)`; this checks exit by either endpoint's gain or the child's promotion. The general backward solve also includes histories leading **into** a census and uses the same included-edge product at every step.

These are necessary physical inequalities, even though the finite formulation uses no exact age process and no time-one evaluation of a whole large ball. The unstratified version barely moves the matched model (at the five listed objectives the displayed differences were all within `4e−12` in relative width). The addition that matters is the constraint on **future mother input from the grandmother's degree at birth**. The LP still allows a coarse continuum of side gain rates; none of those extrema assume independent exterior donors, and no payoff from the D18 slack ledger is used.

## Replay and scope

```sh
OPENBLAS_NUM_THREADS=1 ./round2/.venv/bin/python -Bu \
  identification/time-cuts/check.py --held-out
```

From `embedding-flow/` the run used 718 columns and 803 baseline inequality rows, adding 108 rational rows for individual edge censuses only. Shell objectives were not added as constraints. The existing solver proposes both duals with HiGHS and repairs them in rational arithmetic; the table displays rounded values of these certified outer bounds on the **unnormalized** masses (`η`, not `π`).

| Edge test | Baseline outer interval (display) | With time cuts (display) | Width reduction |
|---|---|---|---:|
| `(2,1)→(2,0)` | `.06356008 – .07442998` | `.06444678 – .07441970` | 8.25% |
| `D_p=3`, virgin child (held-out sum) | `.03826728 – .05291505` | `.03887381 – .05244484` | 7.35% |
| `D_p=4`, virgin child (held-out sum) | `.02676279 – .03907520` | `.02718549 – .03870913` | 6.41% |
| `(4,2)→(2,0)` | `.00748733 – .01338629` | `.00756772 – .01335757` | 1.85% |
| `(4,2)→(3,1)` | `.00072103 – .00303196` | `.00075465 – .00294246` | 5.33% |

The cap-5 attempt (`--cap 5 --held-out`) hit `HiGHS Status 0: Error` on the baseline optimization and provides no larger-cap result. These are improvements on physical census objectives, not evidence that the paired production contrast tightens by the same percentage. The next useful test is a matched signed edge-weight objective; for signed rewards the `low`/`high` shortcuts above cannot be used without separating the signs or checking the discount corners directly.
