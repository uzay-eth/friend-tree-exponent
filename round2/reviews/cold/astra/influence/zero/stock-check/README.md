# Independent signed large-leaf expansion

`check.py` is a rational checker with its own labelled directed-arrow enumerator. It uses the true silent-parent dynamics on every old arm, including sends from an arm root into its children and sends back up from those children. An arm is `(leaves, core_children)`; equal arms still have distinct incidences. The evaluator accepts a bounded signed payoff `h` as a pointwise interval oracle with global limits, works at positive real discount `s`, and returns outward intervals for `F_h` and the **added-root-leaf** difference. The default fixtures use `s=.56572365`.

## Exact coefficient equations

The physical generator on `(D,Y)` is `G=A+(D+P−k−1)Δ_D+E/D`. A root leaf gains at rate `L`; the old children send up at rate `P`. A root promotion, at rate `L/D`, appends `J` **at the same degree** (so the leaf count falls by one in that same event). A root arrow into each old arm adds one leaf to that arm at rate `1/D`. These give the stated `A,E,P` and the toll `1+(Σh(Y_i)−k−1)/D`.

Write `U_m=1/s+Σ_(j=1)^m B_j f_j` and `R_b=(b−A)^−1`. Then

```
f1 = R_(s+1)[Σh(Y_i)−k−1]
f2 = R_(s+2)[(k+2−P) f1+E f1].
```

At order `j≥2` the input to `R_(s+j)` is `(j−1)(j+k−P)f_(j−1)+Σ_(i<j) (i)^rising(j−1−i) E f_i`. Start from `Δ B_i=−i B_(i+1)`, `DΔ B_i=−i B_i+i²B_(i+1)` and repeatedly expand `B_i/D=B_(i+1)+i B_(i+1)/D`. The exact uncanceled residual, with sign **reward minus row**, is

```
r_m = B_(m+1) m(m+k+1−P) f_m
      + (B_m/D) Σ_(i=1)^m (i)^rising(m−i) E f_i.
```

The checker also tests this formula through order four with arbitrary rational trial functions and a literal full physical generator. Those identities do not assume the trial function is a resolvent solution.

## Independent finite interval coefficient evaluator

Let `a(Y)=(s+1−L^0_Y)^−1(h−1)(Y)`, `c=(a(J)−1)/(s+1)`, and `v(Y)=1−1/D_Y`. Then `f1=c+Σa(Y_i)` and `E f1=Σ[a(Y_i^+)−a(Y_i)]−(k+1)a(J)`. In particular for **any signed constant** `h=q`, `a=(q−1)/(s+1)` exactly. For `h=1`, `f1=−1/(s+1)` and `f2=f1 R_(s+2)(2+k−P)`. Its separate numerical control uses `R_(s+2)P=Σ R^0_(s+2)(1/D)(Y_i)+R^0_(s+2)(1/D)(J)/(s+2)`, the last term from convolution with the virgin arrival stream.

Expanding the `f2` forcing as a constant plus one-arm and unordered distinct-pair sums gives

```
g0 = 2c−a(J)
g1(y) = (1+v(y))a(y)+v(y)c+a(y^+)−a(J)
g2(y,z) = v(y)a(z)+v(z)a(y).
```

Set `b=s+2`. Compute `b2=R_b^(old pair) g2`, `b1=R_b^(old single)[g1+b2(·,J)]`, and `b0=(g0+b1(J))/b`. Then `f2=b0+Σb1+Σ_(i<j)b2`. The old pair has independent old-arrow clocks. Immigration contributes `b1(J)+Σ_i b2(Y_i,J)` to `A f2`: one rate-one Poisson arrival appends one arm to the entire list. In particular it neither appends independently to each old arm nor creates a new pair consisting of two immigrants in one jump. The arithmetic uses signed interval sums and multiplies endpoints in the correct order. The recursion stops after a finite **number of old-arm jumps**, pricing the remaining resolvent by global conditional payoff limits divided by its discount; every finite arithmetic rounding is outward. This does not kill the infinite-horizon coefficient. The reciprocal fixture uses that degree is nondecreasing for a sharper upper terminal. Coefficient intervals may be wide for outward `H`, particularly when signs are mixed.

## Full-horizon absolute remainder

This bound is independent of the intervals just computed. Put `M=max|h−1|`, `a*=M/(s+1)`, `C1=max(a*,(1+a*)/(s+1))`, and `b=s+2`. Then `|f1|≤C1(k+1)` and `|E f1|≤(3k+1)a*`. The source for `f2` is bounded by `(C1+3a*)(k+1)(k+2)`. In `A`, arity is exactly `k+Pois(t)`, so integrating the linear and quadratic Poisson moments gives

```
C2 = (C1+3a*)(1/b + 2/b² + 1/b³)
|f2| ≤ C2(k+1)(k+2)
|E f2| ≤ 2(2k+1) C2(k+2)(k+3).
```

Use `B_(p)≤D^(-p)`. The resulting residual polynomials `G_m`, with `|r_m|≤D^(−m−1)G_m(k)`, are implemented in `absolute_residual`. The physical root arity grows only by promotions, thinned from a rate-one clock even when conditional on the current tree. Thus `k_t≤k+N_t` pathwise and `D_t≥D`. For the nondecreasing polynomials here the absolute convolution is bounded by `D^(−m−1) ∫e^(−st) E[G_m(k+N_t)]dt`. Newton differences express the integral exactly: `Σ_j Δ^j G_m(k)/s^(j+1)` since `E[(N_t)_j]=t^j`. No independence between the physical degree and the clock, short-time approximation, or local one-ring bound is used.

The physical process is nonexplosive (total jump rate is at most total size, each event adds one vertex). Stopped Dynkin and passage to the limit give `F_h−U_m=V_s r_m`: for `m=1,2` the growth bounds make `U_m` uniformly bounded, using `k+1≤D`, and `r_m` is integrable by the displayed Poisson comparison. More generally the triangular recursion inductively bounds `f_j` by a constant times `(1+k)^j`: `E` costs at most one extra factor of `1+k`, and the `A` resolvent integrates polynomial moments. The same domination yields order `(1+k)^(m+1)/D^(m+1)` for the value remainder at fixed order and `s>0`.

For the root derivative use `Δ_D U_m=−Σ_j j B_(j+1)f_j`, then add the **two marginal absolute remainder radii** at `D` and `D+1` with the same old arm list. In particular the order-two value formula certifies a derivative after its leading term with an all-degree `O((1+k)^3/D³)` remainder. The included `−2B_3 f2` is exact for the trial expansion but this absolute comparison does **not** certify an `O(D^-4)` remainder after that term; using a third value coefficient would give such an absolute bound for fixed `k`. The cheap Poisson bound can be very wide for small `D`; it asserts no sign there.

## Replay

```sh
nice -n 15 taskset -c 3 timeout 8s python3 -B round2/influence/zero/stock-check/check.py
nice -n 15 taskset -c 3 timeout 8s python3 -B round2/influence/zero/stock-check/check.py --tests-only --plant same-event
nice -n 15 taskset -c 3 timeout 8s python3 -B round2/influence/zero/stock-check/check.py --tests-only --plant shared-immigration
```

The plants respectively drop the simultaneous promotion cross term in the `((k+1)L)` drift and multiply the one append clock by the current number of arms in the ordered-pair statistic. Each must exit 1. The numerical replay includes `0`, `1`, `−2/3`, `1/D_Y`, `L_Y/D_Y`, and `H-pair-TV`. The last uses the paired-TV contraction (`pair/README.md`), `q=.8600511260` for `s≥.565658`. Since `0≤V_s b≤1/s`, the root-leaf difference has norm at most `1/s`; applying the inverse bound gives `|H|≤1/[s(1−q)]`.

The earlier cap-48 signed surrogate profile and stock-check numbers derived from it are withdrawn. They used a reversed signed leaf reward on the paired `v` class (A leaf / B core), which understated the reward by `v(1/A+1/B)`. The checker no longer accepts that profile or uses its invariant box. The constant and specified-payoff checks are unaffected.

## Cold reference: normalized moments and corrected-v2 arm sign

The stock derivation (`stock/README.md#coefficients-and-residual`) writes the recurrence in raw coefficients. Put `u_n=f_n/(n−1)!`. Independently dividing its physical row coefficient gives

```
u1 = R_(s+1) Q
u_n = R_(s+n)[(n+k−P_o) u_(n−1) + (1/(n−1)) Σ_(j<n) E u_j]
r_m = (m−1)! [m(m+k+1−P_o) B_(m+1) u_m + (B_m/D) Σ_(j≤m) E u_j].
```

The rising factors cancel because `j^rising(n−1−j)=(n−2)!/(j−1)!`. `low_stock.py` checks the raw and normalized polynomial recursions agree through order five. The normalization does not itself shrink the physical remainder: the factorial returns when evaluating the residual. On all 200 physical `(D,k)` rows for `5≤D≤20`, the full Poisson majorant with `D^(-m-1)` is minimized at `m=1` over **all** positive orders. To certify the infinite comparison, let `p_m` majorize `u_m`. Positivity gives `p_(m+1)≥[(m+1)/(s+m+1)]p_m` coefficientwise; the first residual factor grows enough to offset this fraction since `s<1`, and the cumulative `E` majorant only grows. Hence `G_(m+1)≥m G_m` and the coarse radius cannot decrease once `m≥D`. Keeping both exact rising-factorial prefactors gives a modest data-dependent improvement (order two for two low-arity rows at `D=20`) over the tested finite bank through order 21. This certifies the stated radius comparisons for the true-`H` enclosure; it does not identify the optimal order for the unknown actual error. The coarse order-one radius already exceeds `0.54` at `D=20,k=0`. For an independently verified uniform approximation `h` with `|H−h|≤ε`, the value discrepancy adds at most `ε V_s C1`, or `qε` for the root derivative; the first coefficient discrepancy is at most `ε[k/(s+1)+1/(s+1)^2]`. A proposal without a certified error cannot substitute for `H`.

At `s=.56572365` the exact arm check uses the **corrected** mixed-leaf-v2 (`surrogate/V2.md`) bounds for the true `H`. The only accepted profile is pinned to that replay. The independent one-arm comparison uses physical promotion `(D,L)→(D,L−1)` at rate `L/D` and gain `(D,L)→(D+1,L+1)` at rate `L+Σ_i 1/D_i ∈ [L,L+k/2]`. No external root send occurs in the silent arm. Both discounted Bellman inequalities are checked with rational arithmetic at both rate endpoints, including the entire v2 tail at cap 48. For `J` the checker further uses four old-arm jumps with literal labelled arrows, then prices the full continuation at the corresponding census endpoints. This gives `V_(s+1)(H−1)(J)≤0.1991078904`. For the leaf-only old arms with `D=3,4,5,10,20` the upper bounds are all strictly negative. On 97 of 1,127 retained censuses the resolvent upper bound is strictly negative. The available bounds do not prove global `H≤1`, nor `V_(s+1)(H−1)(J)≤0`.

Writing `w=V_(s+1)(H−1)`, the leading high-hub constant is `−f1=(1−w(J))/(s+1)−Σ_i w(Y_i)`. It is **strictly positive** for any fixed list of the certified nonpositive arms with at most two additional `J` arms: even with two `J` arms, `f1≤−0.1132998431`. This is a sign of the true leading coefficient at the stated discount and therefore an eventual fixed-arm sign; it is not a claim about the derivative at `D=5…20` or the sign of the newborn scalar.

The smoothed producer's pinned post-refactor replay (`smoothed/envelope/README.md#exact-virgin-certificates`) gives `g(J)=V_(s+1)H(J)∈[.068345438944,.304160175903]` on the original corrected v2 profile. The stock checker verifies its arithmetic consequence without rerunning the producer solve. Subtracting `20000000/31314473` gives `w(J)≤−10475384384010255881/31314473000000000000<0`. With this input, any fixed finite list of nonpositive arms may contain arbitrarily many copies of `J`; for `j` such copies its leading constant is at least `(1−w_*)/(s+1)−j w_*`, where `w_*` is that rational upper bound. Positivity of `g(J)` is not needed for this conclusion. The independently checked older bound `w(J)≤.1991078904` alone proves the at-most-two result, not negativity of `w(J)`.

For all sufficiently large **arm** degrees, use the physical pair coupon with one extra leaf. The true reward `b+CH` is a mixture of leaf value one, core values in the global v2 interval and the silent slot value zero. Match the common slots at mass `1/A` each. The unmatched probability mass on either side is at most `min(1,m/B)` and the reward difference is at most that mass times the global reward width. Applying the paired drift bound (`pair/README.md`) `E m_t≤exp((1+1/d)t)` with initial degree `d` and integrating by an outward rational `1/64` grid proves `H≤.9817086146` for every initial arm degree `d≥192`. Degree is nondecreasing, so `w≤−.0116823843` there. Feeding the tail and the intermediate coupon rows back through the single-arm Bellman comparison certifies `w≤0` on 4,750 of the 18,335 finite censuses below 192, without the provisional `g(J)` input. The bare `D=2,L=0` single-core-child census still has upper bound `1.1387135830`; this comparison does not certify the sign for every arm or the global leading sign.

The smoothed producer's sharper oriented comparison (`smoothed/envelope/faces.py`) was also replayed against the rematched signed LC profile (SHA-256 `8512947b923e4c2ac2877f964ff39866a4d0b8057292bbd9a6ad061c7bcefc82`). It gives an analytic `g≤.531408704001` for every arm degree `D≥128`, then checks the remaining finite census rows back through **degree 94**. All rows at degrees 94–127 have `g<1/(s+1)` (the worst census at 96 has upper `.627626659651`). Thus all arms of degree at least 94 can be used as nonpositive arms in the leading-sign certificate, conditional on the rematched profile enclosure. This replay does not use the provisional `g(J)` bound or establish the result below degree 94.

```sh
nice -n 15 taskset -c 3 timeout 12s python3 -B round2/influence/zero/stock-check/low_stock.py \
  ft-zero-envelope-v2-s56572365-cap48-tight.txt
```
