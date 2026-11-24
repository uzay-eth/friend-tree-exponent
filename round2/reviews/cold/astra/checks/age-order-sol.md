# Q3 age order: independent cold review

**Verdict:** The partial result (`open/age-order/README.md`) holds in the original edge-seeded tree. For each fixed deterministic `u`, a single explicit `v0(u)` works for all `v >= v0(u)` and all `0 <= x < 1`. The size-204 example gives the opposite limiting order at the hub threshold **conditional on one labelled tree**. Neither unconditional adjacent-label question in §7(3) is settled.

## Physical clocks and the stock estimate

In continuous time a directed arrow `a -> b` has predictable intensity `1/D_a` when the edge is present; the newborn attaches to `b`. At size `j`, divide by `j` for the one-step probability. An unpromoted leaf has just one incoming direction, from its permanent parent. Its own ring gives its parent a new leaf; it does not promote itself. A tagged leaf's ring creates a fresh tag at its parent. Tag loss is possible only on the parent's ring selecting that tag, and the total loss rate is `B_t/D_parent <= 1`. These statements hold with the focal leaf excluded from the tags, including after it promotes. Later changes in the exterior tree can affect the denominators without changing these rules.

To make the lower stock coupling explicit, use a rate-one loss proposal and thin for physical tag removals; the lower stock `U` loses one at every such proposal before absorption. When `U <= B_t`, retain a physical tag birth as a birth of `U` with probability `U/B_t`. This preserves the physical rates and gives lower-stock birth rate `U`. The stopped martingale `exp(-(t∧τ))(U_(t∧τ)-1)` has mean `k-1` and squared deviation at most `k+1/2`, since `E U_t <= k exp(t)` and the squared jump sizes have factor `exp(-2t)`. A first violation of `U_t >= (k/2) exp(t)` requires excursion at least `k/2-1`; optional stopping or the squared-martingale maximum inequality gives the claimed exceptional probability. Up to that violation, the incoming arrow intensity to the excluded leaf is at most `(2/k)exp(-t)`. Integrate the **stopped predictable** intensity and add the violation probability. This avoids conditioning on an entire future good-stock path, which would not in general preserve the physical hazard. For `k >= 8` the result is less than `16/k`: `k/2-1 >= 3k/8` and `k+1/2 <= 17k/16` suffice. The preliminary discussion of stopping on first hit of `m` tags in the author note is unnecessary; the actual crossing is against the moving envelope.

For the gamma comparison, replace every removed tag with a phantom that keeps growing at the tag birth rate. The combined population is a Yule process of initial size `k`. Each loss at time `t` contributes a conditional mean of `exp(-t)` to the scaled phantom count at a later observation time. The compensator therefore bounds the total scaled expected loss by one. Since the live tags are bounded by the physical degree, Fatou gives `E[(G_k-W)_+] <= 1`, hence `E W >= k-1`. No independence between the lost mass and `G_k`, or between physical `W` and the size limit `S`, is used. This comparison also gives `P(W=0) <= P(G_k<k/2)+2/k` at large leaf count. The lower birth-death comparison in hub-children (`fringe/hub-children.md`) is valid by giving each tag loss proposals at rate `1/d` once the focal degree is `d`, and removing the tag from the lower population at every proposal. Survival gives unbounded leaf count. The large-stock failure bound applied at its first hit of `k` then proves positive mass on survival almost surely. Likewise, at each first degree hit `d`, the new leaf starts a comparison that survives with probability `1-1/d`; letting `d` grow rules out unbounded degree at zero mass.

## Uniformity at the two endpoints and an explicit cutoff

**Near 1.** For `u=1,2`, prescribe the size-three star centred at `u`, with probability `π=1/2`. For `u>=3`, first build the star at vertex 1 through size `u`, then attach label `u+1` to leaf `u`. This has probability `π=1/[u(u-1)^2]`, and `u` now has one leaf. Thus take `m=3` or `u+1`, respectively, and `B=m-1`. Insisting that every subsequent birth through size `n` attaches to `u` has total probability at least `π (m-1)! n^(-B)`. Indeed at size `j` the `j-B` leaves supply that much incoming rate, and the product telescopes to `(m-1)!` divided by `B` factors each less than `n`.

At that final state restart the clock; `S` has mean and variance `n`, while `0 <= W=S Z_u <= S` and `E W >= n-B-1`. With `ε=1-x`, choose `n=max(m,16,ceil(8(B+1)/ε))`. Chebyshev and Markov bound the union of `S<n/2` and `S-W>=ε n/2` by `1/2`. On its complement the strict inequality `Z_u>x` follows, including the case `S=n/2`. With `A=m+17+8(B+1)` this gives `P(Z_u>x) >= c ε^B` for `0<x<1`, with `c=π(m-1)!/(2 A^B)`. The degree process of a newborn `v>=3` admits the comparison urn from red count 1 at size `v`: at equality its gain probability is at most the urn's, and a one-unit gain from a strictly lower integer count cannot overtake. Its limit bounds `Z_v` with upper tail `ε^(v-1)`. Let `p` be the least nonnegative integer with `2^p c >= 1`. For `v>=B+2+p` and `1/2<=x<1` the ratio of the upper to lower bounds is at most `2^(-(p+1))/c <= 1/2`. At `x=1` both strict tails vanish, as required by support.

**At and near 0.** At size `j=v-1>=3`, condition on the tree **before** the birth. All leaf targets combined have weight at most `C_j`; core targets with degree `<K` have combined weight at most `K C_j`. For the other bad core sector, `L<D/2` implies `D<2c_core`; the sum of core degrees in the induced core forest is at most `2 C_j`. Its combined weight is thus at most `4 C_j`. The remaining targets carry at least `K/2` old leaves, explicitly excluding the new vertex. For `K>=16` the stock bound gives

```text
P(Z_v>0) <= (K+5) E[C_j]/j + 32/K.
```

This uses the implication that a hub must first promote, but does not equate promotion with hub status. For an elementary fully explicit mean bound, let `H_j` count vertices of degree at least three. At size `j>=4`, no degree-two vertex has two leaf neighbours. The rate of leaf targets is at most `(C_j+H_j)/2`; the rate of degree-two targets is at most `3(C_j-H_j)/2` since at least one of their two neighbours is a core. Thus for `P_j=C_j+(2/9)H_j`,

```text
j E[ΔP_j | T_j] <= (5/6) C_j + (1/6) H_j <= (5/6) P_j.
```

Since `P_4<=3` and `sum_(i=4)^(j-1) 1/i <= log(j/3)`, iteration gives `E C_j <= 2 j^(5/6)` for `j>=4`. This supplies effectiveness without assigning any quantitative rate to dominated convergence or rerunning the numerical upper certificate.

Use `L=c/2^B`, `K=max(16,ceil(128/L))` and `T=ceil(8(K+5)/L)`. For `v-1>=T^6`, the bound on `P(Z_v>0)` is at most `L/2`, while `P(Z_u>x) >= P(Z_u>1/2) >= L` for `0<=x<=1/2`. In particular this covers the atom at zero; the comparison beta urn alone has no useful bound there. A sufficient integer cutoff is

```text
v0(u) = max(u+1, 5, 1+T^6, B+2+p).
```

It is deliberately coarse and does not imply adjacent-label monotonicity.

## Size 204: exact hazards and what conditioning buys

In the prescribed shape `D_1=202`, `D_2=2`, `u=203`, `v=204`. The initial promotion hazards are exactly `1/202` and `1/2`, respectively. The former varies as `D_1` grows; it is not frozen at its initial value. The old leaf has exactly 200 other tags at 1, so the stopped bound gives `P(Z_203>0 | shape) <= 29851/980100`.

Until the first of `2 -> v`, `v -> 2` and `1 -> 2`, vertex 2 keeps degree two and `v` stays a leaf. The first two channels can be run as independent fixed-rate clocks of rates `1/2` and `1`, using the two independent half-rate splits of vertex 2's clock while its degree is two. Use a separate proposal clock for `1 -> 2`: its actual initial rate is `1/202` and it never exceeds that value, since degree never decreases. The residual rate for vertex 1's other directions is assigned at their physical rates; this is a generator-preserving thinning, not an independence claim about actual accepted arrows. Winning the race against the proposal ensures promotion with probability at least `101/304`. The fresh leaf then gives conditional hub probability at least `1/2` by the birth-death coupling. Thus `P(Z_204>0 | shape) >= 101/608 > 1/8`, sharpening the author's valid bound.

The exact shape probability is `1/(203·202²)`: the star must hold through size 203 (probability `1/202`), then label 204 attaches to leaf 2 (probability `1/(203·202)`). By continuity of strict tails as `x` decreases to zero, the conditional gap also persists for all sufficiently small positive `x`. Decomposing the unconditional difference by this event leaves an uncontrolled complementary contribution. The exact event weight and conditional gap certify a conditional reversal only. They give no sign for the unconditional difference, and the initial hazard comparison alone gives no limiting conclusion.

## Targeted replay

`age_order_sol.py` reconstructs the directed-arrow rates and full shape probability without importing the author checker; verifies exact rational endpoint budgets and stock bounds; and checks the elementary core drift on every chronological tree at sizes 4–7. The replay checks the finite arithmetic and generator, while the stopping and all-threshold assertions are analytic.

```sh
nice -n 15 taskset -c 3 python3 -B round2/checks/age_order_sol.py
```
