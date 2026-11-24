> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../../README.md#claim-25) for exact scope and citations.

# Marked-neighbor entry through degree 1999

The signed seam entry `π I` in the [joined-edge balance](../README.md) is in **`[-8.5300112e−3, -3.6674464e−4]` per core**, on the original support. Width drops from `2.8126131e−2` to less than `8.163267e−3`: **3.445× tighter**. This bounds the physical seam directly. It uses the exponent endpoints only to enclose the exponential age discount and the birth-rate multiplier, not the whole-score identity or the original bracket as a price for the positive sector. There are no H3 rewards in the calculation.

## Exact target and change of clock

Let `u` be the focal core, `τ` its first gain to ordinary degree 2000, and `W(u,v)=−g_pair((2000,c_u),(D_v,c_v))` at the crossing. `W≥0` on the entire original face (checked); it is zero when the other endpoint has `D_v≥180`. The target magnitude is

```
−π I = μ E_birth[ exp(−μ τ) Σ_{v~u} W(u,v) 1{c_u,c_v≤106}; τ<∞ ].
```

Before crossing, the root is at `D=1999` and the gain has intensity `1999−c_u+Σ_core 1/D_v`. This expression uses all true incoming arrows, including the marked neighbor. Apply the stable-age identity and compensate that actual first crossing to get the display; the root arity and every neighbor mark are unchanged by the gain. This is the same full-input seam as in the parent edge balance, including its arity handling. The target pays the *new selected edge weight*, not the old-to-new raw difference of the low edge. There is no synthetic draw of a low neighbor conditional on the root census.

## Actual child edge as a mark

Each physical promotion before `τ` registers one distinct child, which starts at `(q,c)=(2,1)` with its real host and with the host's actual degree `d` at that promotion. Follow that same edge through the remaining race to `τ`. At every retained pair state the two gains include `1/q` into the host and `1/d` into the child. The remaining child input is in `[0,(c−1)/2]`, and its own promotion rate is `(q−c)/q`. The host has at least two core neighbors (its mother and this child); over surviving `2≤c_u≤min(d,106)` the actual host gain is between `d−min(d,106)+1/q` and `d−2+1/q+1/2`. Both tagged degrees and the child arity are retained as evolving coordinates; the host arity is relaxed within this interval, rather than replaced by a typical value. Reciprocal sends alter only the receiver. The residual rates may be adaptive to the full environment.

For the upper solve, price the terminal corner at the **actual host crossing**, maximizing over its clamped arity at that instant. If the marked child exceeds degree 24, the upper continuation carries the largest admitted weight above degree 24 (`.025237214` rounded up), discounted through the remaining host gains at rate at most `d`. This tail is not paid early. Exiting the original host support has true payoff zero and the upper bound can overpay by continuing. For the lower solve, minimize the terminal corner over the same host arities, kill on marked child degree exit, and kill on any subsequent host promotion (hazard at most one). This last stopping rule keeps the post-registration host arity admitted and leaves a genuine lower bound for each physical child; it does not suppress birth of later children in the outer counting process. The upper and lower barrier inequalities use min/max of a linear-fractional ratio over the two gain-rate endpoints after denominators are cleared. The replay checks exact integer inequalities at each corner, rounding outward. The rows apply conditionally at actual registration, including arbitrary history and competition among siblings.

Next, compensate all actual host promotions before the crossing from physical birth `(2,1)`. At each root promotion the recurrence adds the mark price at the real current degree and continues with the host arity incremented. Conditional bounds make the expected sum valid without independent child lifetimes. The outer recurrence stops after the registration at old `c=100` (post `c=101`, still admitted). For the upper bound, unregistered children number at most `c(τ)1{c(τ)>100}` and each is worth at most `CAP`; use monotonicity from crossing to the same exponential observation to integrate the geometric arity envelope. The lower bound can discard the missed nonnegative contributions. The unique mother is not a child registration: her mark may be aged at the host's own birth. She is kept in the target, with zero lower price and upper price `CAP·μ_upper·P_π(D≥2000)_upper`; the degree hitting bound is independently replayed from the physical gain and promotion corners. This is a bound on the unmapped mother contribution, not an assertion that it is irreducible.

| Seam magnitude component | Certified interval per core, rounded outward |
|---|---:|
| Child subtotal, including upper registration overflow | `[3.6674464e−4, 6.0747133e−3]` |
| Mother | `[0, 2.4552980e−3]` |
| Late registrations, already in child upper | `[0, 6.259e−19]` |
| Total magnitude `−π I` | `[3.6674464e−4, 8.5300112e−3]` |

These bounds can replace the seam interval in the joined-edge identity. The mother starting-mark law, the slack from relaxed host arity while tracking a child, and the constant child tail price remain possible sources of improvement. The original width comparison uses exactly the preceding all-neighbor interval, not a different seam mask.

## Replay

```sh
taskset -c 5 nice -n 19 python3 -IBu round2/occupation-quintet/true-drift/edges/entrance/check.py
```

The script verifies the original detached source pin and reads only the sparse high-edge face and small controls, with no whole-table scan. The five ordinary-tree fixtures check both marked sends against the actual donor slots, the direct crossing price, the child degree cutoff, and zero on the double-tail face.
