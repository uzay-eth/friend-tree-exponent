> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** supporting Astra work; the cited inputs and any recorded review are scoped in the linked claim index and historical record; no blanket independent review is asserted. See the [release index](../../../../../README.md) for exact scope and citations.

# Arithmetic at the physical crossing: degree–height exclusions and a truncation test

**Result.** Using the exact original Gold2 + three funded Astra floors, an algebraic exponent of degree 2, 3, 4 or 5 must have primitive polynomial height at least **63, 15, 6 or 3**, respectively. These are sharp for this interval. A rational specialization with square discriminant in Long's first-promotion transform would require denominator at least **94,251**; without that extra condition the lower bound remains **426**. None of these bounds proves irrationality.

For an all-rational exclusion, positivity and complex-line control do not supply the missing denominator estimate. Below is an actual physical truncation lower bound, a reduced-denominator test that would suffice with an additional arithmetic input, and a rational-root countermodel in which the denominator growth can be computed exactly. This does not replace the physical transform by a toy transform when drawing conclusions about the exponent.

## 1. Exact exclusions from the funded physical interval

Take the [published composition of the original Gold2 interval and three occupation floors](../../../occupation-quintet-gold2/rationals/README.md) as input. The endpoints used here, without rounding, are

```
L = 126905728962161948969119589499374909399 / 224328301684105796116575000000000000000
U = 127993535256830043 / 226242120050000000.
```

For degree `d`, height means `max_i |a_i|` of the primitive integer minimal polynomial `Σ a_i x^i`. The following table has both an exclusion below the threshold and an irreducible polynomial attaining it with one real root strictly inside `[L,U]`.

| Degree | Minimum height in the enclosure | Witness | Irreducible modulo |
|---:|---:|---|---:|
| 2 | 63 | `2x²−63x+35` | 3 |
| 3 | 15 | `13x³+13x²+15x−15` | 11 |
| 4 | 6 | `x⁴+3x³+3x²+6x−5` | 17 |
| 5 | 3 | `x⁵+3x⁴+3x³+3x²+2x−3` | 2 |

For each height strictly below the table, enumerate the nonconstant coefficients with positive leading coefficient. Integer interval arithmetic encloses the range of their sum on `[L,U]` and determines every eligible constant term. For survivors, either a strict enclosure excludes zero or the derivative has one sign and both endpoint evaluations have the same strict sign. The checker permits bisection with a hard depth limit and fails unless every candidate is excluded. Each witness changes sign and has no monic divisor of degree at most half its degree over the indicated finite field. In particular the old quadratic height-seven witness `3x²−7x+3` no longer fits. The tabulated exclusions apply to nonmonic as well as monic minimal polynomials.

The square-discriminant question is equivalent to `s=u/v−v/u` for reduced positive integers `u>v`. The reduced denominator is **exactly** `uv`: neither `u` nor `v` has a prime in common with `u²−v²`. In this interval `u<2v`, since `t−1/t` increases and is `3/2` at `t=2`. Thus stopping at `uv≤94251` uses only `v≤floor(sqrt(94251))` and `v<u≤min(2v,floor(94251/v))`. The sole hit is `(u,v)=(353,267)`, giving `53320/94251`. This improves the prior `57200`; it does not assert a special-value result for the physical trigger at that subclass.

## 2. What the physical complex line proves

Write `F` for the actual [positive degree-two trigger](../../../../../../astra/source/round2/finite/gauge/physical-trigger.md#4-positive-chronological-series-and-uniform-error), normalized by `F(μ)=1`. Split off the initial, zero-transition reward `1/[2(s+2)]`. At `s=μ+iy` the rest is the Fourier transform of a positive discounted measure of mass `1−m`, where `m=1/[2(μ+2)]`. Thus

```
1 − |F(μ+iy)| ≥ m [1 − (1 + y²/(μ+2)²)^(-1/2)],     y ∈ R.       (1)
```

This proves a strict gap off the real crossing and a quantitative quadratic lower bound near it. It also explains why positivity does not give a uniform gap as `y→0`. In fact the full chronological expansion shares the same initial holding factor: stripping it, use absolute domination on the critical line to get `|F(μ+iy)|≤(μ+2)/|μ+2+iy|`. Sharpness near zero of the *order* in (1) also follows from the finite second moment of the discounted time law (use the positive analytic convergence to the left): its characteristic-function magnitude is `1−Var(t)y²/2+o(y²)`, with `Var(t)>0` because the initial contribution has positive density on an interval. The real derivative is strictly negative. Accordingly the renewal response `1/(1−F)` has a simple positive-real residue and no other pole on the critical line. The analytic continuation and this isolation locate the pole in the complex plane but do not impose an arithmetic type on its real coordinate.

## 3. Growing physical truncations: a sharper necessary scale

Let `F_N` be the actual size-killed resolvent: kill upon the first step to `N+1`, and sum all lengths under that size cap. The chronological expansion has positive coefficients, so all paths that hit size `m>N` contribute to `F−F_N`. For each `m>N≥20`, take the path growing the initial degree-two, mark-`1/2` star to degree `m`, making one retained direct reset to `J_(1/20)`, and ending there on its trigger reward. These are distinct paths for each `m`, including after the reset, and cannot occur in `F_N`. At real discount `s∈[L,U]`, their first product is

```
Π_(d=2)^(m−1) (d−1)/(d+s+1/2−1/d)
 ≥ Γ(m−1) Γ(s+5/2) / Γ(m+s+1/2)
 ≥ c m^(−s−3/2).
```

The last bound is uniform on the compact discount interval by the Gamma ratio (or by summing logarithms); the retained reset rate divided by its holding factor is at least `c/m`, and the terminal reward at the new virgin is uniformly bounded below. Summing **all** those omitted paths gives

```
F(s) − F_N(s) ≥ c' (N+1)^(−s−3/2)      (N≥20, s∈[L,U]).           (2)
```

In particular the physical size tail cannot decrease geometrically. This replaces the older single-path `N^(−4)` lower estimate by a power about `2.066` at `μ`. A coarser explicit check needs no Gamma functions: the individual path pays at least `1/[m²(m+1)(m+2)]` for `s≤1`, so summing only `N+1≤m≤2N` gives at least `1/[16(N+1)^3]`.

For a fixed trial `r=p/q∈[L,U]`, rational in lowest terms, let `D_N(r)` be the **actual reduced positive denominator** of `1−F_N(r)` (equivalently of `F_N(r)`). Suppose `F(r)=1`. Then `D_N(r)(1−F_N(r))` is a positive integer for every `N`. With the [published size-tail bound](../../../../../../astra/source/round2/finite/gauge/physical-trigger.md#41-a-single-finite-rational-solve-with-an-error-tending-to-zero), this implies

```
D_N(r) ≥ (1/2837) ((N+1)/2)^0.11764.                            (3)
```

In particular eventual denominator boundedness would suffice, and more generally **any unbounded sequence** of cutoffs with denominators strictly below the right side of (3) excludes that trial. To exclude all rationals by this route one needs a statement of this kind for every candidate (with constants allowed to depend on the candidate), or a stronger nonzero-integer-forms construction. Conversely, by (2) a sequence with `D_N(r)/(N+1)^(r+3/2)→∞` cannot produce small cleared values at a hypothetical rational crossing even after exact value reduction. These statements involve reduced denominators at evaluation, after both function-level and specialization-level cancellation. The raw determinant denominator is an upper bound, not evidence of that growth after cancellation. Its resultant lower bound with growing matrix dimension does not become a small-value proof. The chronological primitive-polynomial obstruction in the [parent note](../README.md) is for a different family. Prime appearances in coefficients at infinity are likewise not valuations at `p/q`; in a live degree-`p` row the holding factor can cancel the same `1/p` that created a formal coefficient denominator (see [the earlier specialization analysis](../critical-value/README.md)). No bound proving (3) fails is available for the actual rational specializations here.

## 4. Sharp limitation of using complex-line and positivity data alone

For any rational `r∈[L,U]`, including `241/426`, set `c=1−1/[2(r+2)]>0` and use the positive rational-history response

```
G_r(s) = 1/[2(s+2)] + c Σ_(n≥1) [1/n²−1/(n+1)²] (r+n+2)/(s+n+2).
```

The residues are positive and summable; the tail converges locally uniformly off the negative integers. Thus it is meromorphic with infinitely many genuine poles and transcendental over `C(s)`. It is a positive Laplace transform, strictly decreasing on the nonnegative real axis, and `G_r(r)=1` is a simple crossing by telescoping. Its initial term is identical to the physical direct term, hence (1) holds with `r` in place of `μ`. Moreover `G_r=O(1/|y|)` uniformly along vertical lines by summable residues; compactness off the crossing and the local simple zero give a strip immediately to its left in which `1/(1−G_r)` has the single crossing pole. If a common first holding factor is desired, use instead `J_r(s)=1/[2(s+2)]+(r+2)/(s+2)·(G_r(s)−1/[2(s+2)])`, a positive sum of serial two-hold paths. It satisfies the stronger shared-factor complex-line bound just as the physical trigger does. Its extra poles at `−n−2` have nonzero residues, and it has the same critical value and finite-prefix errors. Both have positivity, rational histories, the meromorphic crossing, and quantitative vertical separation.

At that exact rational crossing, the finite prefixes of either model satisfy `1−G_{r,N}(r)=c/(N+1)^2`. If `c=A/B` is reduced, the reduced denominator is **exactly** `B(N+1)^2/gcd(A,(N+1)^2)`, and the cleared numerator is `A/gcd(A,(N+1)^2)≥1`. Thus a genuinely polynomially shrinking tail, of a scale compatible with (2)–(3), coexists with a rational critical point and growing reduced denominators; there is no hidden passage to bounded denominators. This model does not reproduce the full physical transition identities or the formal coefficient congruences. Its purpose is to separate the analytic information in (1)–(2) from the arithmetic information that would have to come from the actual generator. Matching complex-line decay and a mere order of convergence would not close that gap.

## Replay

`nice -n 10 timeout 12s python3 -IB round2/open/arithmetic/transform/check.py` uses exact rational/integer arithmetic for the height exclusions, finite-field witnesses, square discriminant, rational-root truncations and explicit physical path lower envelopes. The official endpoints and occupation receipts are inherited from the linked audit, not regenerated here.
