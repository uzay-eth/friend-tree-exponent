# Independent audit: serial-middle Mellin lower bound

**Final pinned verdict: holds.** Reviewed the writer proof and replay (`open/age-order/hub-ranking/mellin/README.md`) at proof/check/data blobs `9e75bd41e8a204b100a086d21e9fea76ad844fb0` / `9786bbdf7230ee86b9aa518d421d9fb471866130` / `defa8793b129e4b3c3f5ef6f23fd4033fc4ff8fc`. The pinned proof gives strict signed conditional order for every real `p≥3` on every old tree; the low-power table supplies only the serial inputs to the separate ledger. The earlier `1/256` signed grid failed its first box; the final `1/512` grid clears all 2,048 boxes and reproduces the same 256 low-power rows byte for byte.

For the signed cost, a degree-`d` core with `L≥1` old leaves has `s_core≥L/[d(d+1)]`, while each leaf has `s_leaf=1/(2d)`. The sum of `(s_leaf−s_core)_+` at the core is at most `A_l(d−1)/(d+1)` for `A_l=L/(2d)`. Three times the latter is at most `A_l+w_core`, using `w_core≥2d A_l/(d+1)`. The all-edge absolute debit mass is at most `ΣD_us_u=A`; the seed has no edge coefficient. Together with the serial tip debit this yields `4U_2/3+2U_3/3` per unit `A`, since the tip bounds are ordered.

**Unrefined tip recurrence.** At promotion the old-parent slot contributes at most `1/k` forever, so the tag's effective affine birth seed from degree two is `1+1/k`, not `2−1/k`. Prior to promotion the parent cannot shrink, the promotion hazard is at most `1/k` and the promotion payoff cap is nonincreasing with the actual parent degree. This proves the constant-hazard boundary upper at parent degree `48` or higher. Below the cutoff the tip's own rate-one ring always increases the parent degree by one; the promotion arrow and this ring are disjoint. With constant exponent `p` and a fixed parent degree `k`, the normalized first-event upper numerator is the promoted payoff divided by `k`, plus the next bound. The checker uses the smaller denominator obtained at `lo`, rounds up, and takes the maximum with the next bound. This last enlargement makes all additional parent degree gains nonpositive generator terms, even when the individually rounded interval caps have no certified monotonicity. Thus the resulting test is a supersolution on the original physical chain. The pointwise promotion ceiling follows from the decreasing, convex `Beta(r,1−r)` moment at `lo`, its integer chord, and multiplication by `1+hi/r`. The grids align with integers, so no chord is extended outside its enclosing unit interval. The large-power locked-tip upper uses a different starting degree, namely one, hence effective seed `1/k`.

Final pinned replay (`taskset -c 8 nice -n 10 timeout 75s python3 -IBu round2/open/age-order/hub-ranking/mellin/check.py`):

```text
serial middle 1 563194756241683/1000000000000000
serial middle 2 127084241249381/250000000000000
serial middle 3 96109692461153/200000000000000
PASS 256 lower-power boxes and 2048 signed boxes on [3,7], minimum 1892360505093/1000000000000000
PASS p>=7 normalized serial lower 103/280 signed margin at least 22890923/1058158080
```

**The physical lower comparison holds.** The discounted continuation and gain seams are verified below on the entire real interval `[1,7]` for the stated finite cutoffs. At the writer's `D0=30`, return cap 40, the independent lower values at `p=1,2,3` are greater than `.563194756241683`, `.508336964997524`, `.480548462305765`, respectively (the right sides are outward lower ticks). The independent real-power tail argument also works for `p≥7`. The finite serial calculation by itself does not certify the signed old-edge ledger at small exponents.

## Terminal lower comparison really uses physical stock

Track actual leaves of the mark. Every selected leaf gives a new tracked leaf at rate one; more gains may be added. Loss of tracked leaves has rate at most one, since all such promotions are choices of the marked initiator. For `j≥2` set `F_p(j)=Γ(j−1+p)/[Γ(j−1)Γ(p+1)]` and set `F_p(0)=F_p(1)=0`. For real `p≥1` this is simply the positive polynomial `∏_(k=1)^(j−2)(p+k)/(j−2)!` at each integer `j≥2`, with empty product one. It is nondecreasing in stock, and the worst generator gives

```
(A−p) F_p(j) ≥ p F_p(j) [1/(j−1) − 1/(j+p−2)] ≥ 0.
```

The formula includes the loss to `j=1`, exactly zero when `j=2`; at `j=1` the birth drift is positive and at zero the inequality is trivial. In particular at `p=1` equality in the worst generator at positive stock `j≥2` is harmless. Do not use a Gamma formula at its pole or make an absorbing-state convention that retains positive stock reward after total extinction.

Thus the discounted test is a submartingale after localization. To remove localization or pass to infinity, the *tag degree* (not the entire exterior) is bounded by a rate-one per-member Yule process starting from the current degree, via the physical gain bound `r≤D`. Its time-normalized moments of any fixed finite order are uniformly bounded, so the discounted polynomial tests are uniformly integrable (use an order strictly higher than `p`). The Gamma ratio at large **stock** gives `F_p(j)~j^p/Γ(p+1)`. It suffices that the tag's normalized degree converges to `W`: since `j≤D`, the `limsup` of the discounted stock test is at most `W^p/Γ(p+1)` even if the tracked stock itself has no known limit. Reverse Fatou under uniform integrability then gives the lower bound. In the full-stock choice there is also direct convergence since the mark's core-neighbour count increases at most by its own rate-one rings. The horizon and moment limits are at a fixed physical seed. To translate to a fixed-size share, use `W=Y Z` from the *common completion size*, with `Y~Gamma(m,1)` independent of its graph path under the original free law. No such independence is inserted at a clock-defined return time.

## Finite continuation and return after actual last loss

At focal degree `d<D0` and positive stock `j`, the monitored rates are `j` and `j/d` for birth and actual promotion, and a birth changes the state to `(d+1,j+1)`. Let `V` be the normalized discounted lower test. Add discount rate `p`, yielding

```
V(d,j) = [d V(d+1,j+1) + V(d,j−1)]/[d+1+pd/j],
V(d,0) = R_(p,d)(2) V(d+1,1).
```

The zero state is admissible only immediately after the last **tracked physical leaf** was promoted (or in the subsequent designated return phase). It is not a bound for an arbitrary bare core. At that promotion the designated neighbour has degree two. While the focal mark has no next gain, its degree stays `d`; at neighbour degree `k`, its directed return hazard is `1/k`, and its own growth hazard is at most `k−1+1/d`. Any other focal gain is a successful reset with a new actual leaf. Define the discounted return value with terminal zero at `K` and

```
R(k) = [1/k + (k−1+1/d) R(k+1)]/[p + 1/k + k−1+1/d].
```

The needed pessimistic-growth sign is `R(k)≥R(k+1)`; it follows analytically by backwards induction from `R(k+1)≤1/[1+p(k+1)]<1/(1+pk)`, since the recursion is a weighted mean of `1/(1+pk)` and `R(k+1)`. Also `R≤1`, so other focal gain arrows only help. The reset at focal gain uses a conditional physical continuation from the attained state, no independence of return time and that state. The terminal stock value at `D0` and the seam `V(d+1,j+1)≥V(d,j)` complete the comparison against arbitrary extra gains. The monitored exits are finite, because the finite degree return phase has genuine hazard at least `1/(K−1)` until exit and each positive-stock stage has a tracked birth clock. Degree increases at each successful return and stock decreases at fixed degree after promotions; no infinite series of retries is hidden before the finite cutoff.

**Analytic gain seam, for every real `p≥1` and any finite cutoffs.** Induct downward from the nondecreasing terminal row. In row `d` write `v_j=V(d,j)` and `G_j=V(d+1,j+1)`. The child row is nondecreasing, hence so is `G`. At `j=1`, `(1+p)v_0≤V(d+1,1)≤G_1` by `R(2)≤1/(1+2p)`. At `j≥2`, if `v_(j−2)≤v_(j−1)`, the preceding recurrence gives `G_(j−1)≥(1+p/(j−1))v_(j−1)`, so `G_j≥(1+p/j)v_(j−1)`. In either case the latter inequality is exactly the criterion for `v_j≥v_(j−1)`. It also gives `v_j≤G_j` by substitution into the current recurrence. This closes both inductions and removes the mesh from the physical-comparison argument. Rational boxes remain useful to enclose numerical values uniformly in `p`.

## Entire real parameter interval

The replay does not infer intermediate signs from sampled integer powers or from convexity of the moment. It uses rational interval reconstruction: enclose each terminal polynomial at the two rational ends, round down/up to integer ticks of scale `10^15`, and propagate positive recurrence numerators outward with the denominator extrema. On the same interval check each gain seam against the lower next-state tick and upper current-state tick. Return monotonicity is also tested via the undiscounted numerator of `R(k)−R(k+1)`, `1/k−(p+1/k)R(k+1)`. The boxes cover all real `p∈[1,7]` including both ends. This proves the seam even if the normalized serial moment or its individual state values fail to be monotone in `p`.

## The unit-loss argument for all real `p≥7`

Use a scalar stock process started from one, with rate-`j` births and rate-one losses when `j>0`, absorbing at zero. It embeds below the physical stock: use a subset of the actual leaves, compensate at most one actual tracked loss with the remaining rate as a discarded track (never a fictitious physical promotion). At any extra physical gain ignore it. The monotone stock-rate comparison, or this direct tracking construction, preserves the stock floor under arbitrary exterior input. Use the **unshifted** rising factorial `A_s(j)=(j)_s/Γ(s+1)` for `j>0`, with value zero at the cemetery state for **all** `s≥0`, including `s=0`. For `s≥1`, the gain term cancels the discount and the loss term equals `−A_(s−1)(j)`. This identity includes `j=1`, where the loss jumps to zero. Consequently for `r_s(u)=E[e^(−st)A_s(J_t)]`, `u=1−e^−t`, one has `r'_s=−r_(s−1)`. For the final index `s≥0` used in the remainder, `0≤r_s≤1` by domination by the no-loss Yule family (or direct pure-birth eigenfunction compensation). Normalized tests are uniformly integrable by the same Yule moments.

For real `p≥7`, integrate this derivative identity seven times, with last integrand `r_(p−7)≤1` to bound the negative remainder from below. At `u=1` the resulting lower bound is `Σ_(j=0)^7 (−1)^j/j!=103/280`. Using fewer than seven differentiations at `p=7` is unnecessary; no derivative at an index in `(0,1)` is invoked, and at `s=0` the quantity is survival probability, not the constant one. This floor strictly exceeds `[4(1/2)_7+2(1/3)_7]/[3Γ(8)]`. Each term of that debit is decreasing for real `p≥7`, since for `a<1` the ratio `Γ(p+a)/Γ(p+1)` has log derivative `ψ(p+a)−ψ(p+1)<0`. This comparison is valid independently of the low-power old-edge ledger. It applies to any ledger whose adverse normalized moment is bounded by that displayed debit, as with the one-third split and the locked leaf caps. It does not splice a signed comparison at `p<7` without matching adverse bounds there.

## Replay

```text
taskset -c 8 nice -n 10 timeout 65s python3 -IBu round2/checks/preplacement-mellin-sol.py
PASS 110 exact terminal drift rows, including equality at p=1
PASS p>=7 unit-loss Taylor floor strictly beats the descending Gamma debit at p=7
p=1: baseline lower tick 563194756241683/10^15 (D0=30,K=40)
p=2: baseline lower tick 508336964997524/10^15 (D0=30,K=40)
p=3: baseline lower tick 480548462305765/10^15 (D0=30,K=40)
PASS 192 rational boxes on [1,7] at D0=8, K=40; uniform lower tick 331354814525722/10^15
PASS 192 rational boxes on [1,7] at D0=12, K=40; uniform lower tick 374654701010217/10^15
PASS 384 rational boxes on [1,7] at D0=20, K=40; uniform lower tick 408966217885242/10^15
```
