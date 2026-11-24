> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** supporting Astra work; the cited inputs and any recorded review are scoped in the linked claim index and historical record; no blanket independent review is asserted. See the [release index](../../../../README.md) for exact scope and citations.

# Arithmetic of the original friend-tree exponent

**Results.** The reviewed physical interval excludes rational denominators
below 76, primitive quadratic polynomials of height below 7, and rational
square-discriminant specializations of Long's first-promotion transform
with denominator below 57,200. Separately, the fully reduced primitive
numerator of the **actual positive degree-two trigger**, truncated by
chronological lineage length, grows at least exponentially at its critical
point. That is an obstruction to one arithmetic proof mechanism, not a
claim that the critical point is irrational.

## Attribution and which equation is which

Christopher D. Long's pinned third manuscript, *Transcendence targets for
the friend-tree exponent: exact first-promotion reduction and arithmetic
obstructions* ([source at `5bd25a7`](https://github.com/long-mathematics/friend-tree-growth/blob/5bd25a7d506cb64a8444b86334ab531d295f6c78/friend_tree_transcendence_progress.tex)),
uses Uzay's original-model exponent theorem. Long supplies the inert
physical-generation pressure crossing, and separately a full-environment
spectral-radius crossing, with critical essential radius one. His third
manuscript gives the first inert-star promotion transform as a
`3F2`, proves its positive residual and rational-telescope obstruction,
and gives primitive denominator/height laws and full-model noncancellation
for specified finite-population certificate families. Its final
small-value argument requires much more than noncancellation. In
particular his `primitive_population_relations.tex` proves a matching
order for **his** normalized family; the trigger family below is different.

Long's **modified VAN** (including PAN and DAN) has a local independent
reproduction reduction and an equation involving
`F(x)=∫_0^1 u^x exp(u−1) du`. The companion irrationality theorem for finite
algebraic activation parameter uses Siegel–Shidlovskii at rational trial
exponents; it does not prove transcendence even for modified VAN. The third
manuscript explains why the quoted `E*` finite-exception theorem does not
certify the required value at `−1`. None of these facts transfers to the
original process, where old-core sends remain and the first inert promotion
transform alone is strictly less than one at every positive discount.

For the **original** model, Astra's [19-label kernel](../../../../../astra/source/round2/fringe/finite-kernel.md)
retains permanent capped *birth marks* but uses true reciprocal physical
send rates. Its entries are infinite path sums and `r(K_mu)=1`.
The [degree-two trigger](../../../../../astra/source/round2/finite/gauge/physical-trigger.md) removes
only direct births from degree-two promoters: if `K=C+k e_2^T`, then
`T(s)=e_2^T(I−C)^−1 k`, and `T(mu)=1`. `T` is not the first-promotion
transform, nor the [dormant auxiliary scalar](../../../../../astra/source/round2/finite/gauge/isolated-transfer.md).
The auxiliary prime-denominator proof uses terminal `b=L/D`; its unique
large prime comes from a terminal large star. The trigger's terminal
reward is `h=1/2` at virgins and zero elsewhere. That large star contributes
zero to `h`; the cited proof does not show that `T` is transcendental.

The [exact physical-generation marking](../../../../../../reviews/cold/astra/influence/generation/check/README.md)
has `R_z=L^0−a+zQ` and `K(s,z)=zJ(s−R_z)^−1 bL`. Both old-copy transfers
and direct newborns carry one `z`; the inert starting mark is zero.
The occupation source is `s·1`. Thus the exact coefficient equals Long's
discounted physical-generation count. It is not a coefficient of powers
of `K(s,1)` alone. The compact-pencil critical gap, where available,
gives exponential convergence in **generation**, not rationality of those
infinite-state coefficients or an arithmetic rate for a size cutoff.
Long's signed renewal germ has a real critical boundary value; that
result gives neither rational coefficients at rational discount nor a
Taylor-radius bound to the boundary. His degree-weighted comparison
also preserves the bipartite product identities and arbitrary finite
population prefixes with rational exponent one. These are distinct
obstructions: neither is the trigger pole argument below.

## Primitive trigger numerators grow, even after cancellation

Let `T_m` sum all chronological trigger histories with at most `m`
**transitions** before the terminal reward (equation (4) of the trigger
note). Each history has a positive constant numerator and one holding
factor `s+alpha(S)` at every visited row, including the terminal, where
`alpha(S)=a_S+n(S)−1/D_o`. This prefix is a finite rational function.
Let `p_m` denote the primitive integer numerator, with positive leading
coefficient, after rational-function cancellation.

**Pole survival.** The smallest holding parameter is `31/20`, achieved
at `J_20` only. Indeed all virgins have size and degree two, and every
other row has size at least three, reciprocal at most one half, and
mark at least one twentieth. Let `r_m` be the maximum number of visits
to `J_20` on a retained path with positive terminal reward. At
`s=−31/20` every other factor is strictly positive. Therefore all paths
with `r_m` visits have strictly positive leading Laurent coefficients;
there is no cancellation of this pole in `T_m`, or in `1−T_m`. This
uses the rational continuations of finite sums, not negative-discount
Laplace integrals. To make visits, grow the starting `J_2` to degree 20
using 18 leaf gains, then take its nontriggered direct birth to `J_20`.
The same 19-transition cycle can be repeated from `J_20`. That virgin
already earns terminal reward; no exit to `J_2` is required for this pole
bound. Thus `r_m >= floor(m/19)` for `m>=19`.

**No large integer content can defeat the pole factor.** In general,
for a rational function `1−T` with `T` strictly proper, take the reduced
denominator `B` monic over `Q` and write `1−T=A/B`, with `A` monic of the
same degree. The primitive integer numerator is just `cA`, where `c`
is its positive integer leading coefficient. It follows exactly that
`p(s)=c B(s)(1−T(s))`. All zeros of `B` here are negative real holding
poles; primitive normalization can remove integer content but cannot
replace the factor `B(s)` by a small constant. In the trigger family
every monic factor at `mu` exceeds one. Therefore

```
p_m(mu) >= (mu+31/20)^floor(m/19) · (1−T_m(mu)).
```

**Omitted physical path.** Set `N=max(20,m+2)` and grow a single star
from size two to degree `N`, retaining input label 2. Its gain factors
are `(d−1)/(s+1/2+d−1/d)` for `2<=d<N`. At `s=mu<1`, each is at least
`(d−1)/(d+2)`; their product is `6/[(N−1)N(N+1)]`. Use a nontriggered
direct birth to `J_20` at degree `N`, with discounted factor at least
`(N−1)/[N(N+2)]`. The path length already exceeds `m`, and the new
virgin earns terminal reward `1/[2(s+31/20)]>1/6`. There is no need to
leave that virgin. Hence

```
1−T_m(mu) >= 1/[N² (N+1)(N+2)],
p_m(mu) >= (mu+31/20)^floor(m/19) / [N² (N+1)(N+2)].
```

In particular `liminf log(p_m(mu))/m >= log(mu+31/20)/19 >0`.
This applies to primitive numerators after *all* rational function
cancellations, with no height upper-bound estimate or arithmetic
hypothesis on `mu`. It rules out using these single primitive polynomials
as small values for irrationality or transcendence. It applies also to
chronological prefixes with an additional size cap at least 20: the
reset cycles survive, and imposing a second cutoff cannot recover the
omitted path. It does
not apply automatically to unrestricted-length size-killed resolvents:
their finite-matrix poles are eigenvalues, and resummation can change or
remove the repeated holding pole. Nor does it exclude signed combinations
of several prefixes, Padé forms, or other normalizations that change the
integer polynomial.

## Error, height and what finite matrices can certify

For a fixed rational candidate `p/q`, the test at that candidate is
an enclosure of the **unscaled** `T(p/q)−1` excluding zero. Alternatively,
for a nonzero integer polynomial value at a hypothetical `mu=p/q`, the
lower bound is `q^(−deg)`, without any height premise. For an algebraic
candidate, the resultant also brings in height and conjugates, as in
Long's final section. A rational matrix approximation, even with a
strict lower/upper order and convergence, does not give either small-value
test without a quantitative arithmetic comparison. For example the
published size rate for the trigger is `O(N^−.11764)`; the generation
tail is geometric in the number of untriggered newborns **after**
infinite-residence `C,k` are supplied. Neither is a clearing-denominator
bound. The preceding pole/omission result is stronger than comparing
these error upper bounds to an upper bound on height, which alone could
not prove failure of all possible cancellations.

For a concrete available height bound on the rational **size-killed**
system, let `M` be its number of retained types, `H=lcm(1,…,max(N,20))`.
At size `<=N`, every individual rate is a sum of reciprocals with
denominators dividing `H`. The augmented triggered system with feedback
`h ev_J` has an integer characteristic polynomial after row scaling by
`H`, degree `M`. Each row has coefficient-vector absolute sum (summing
entry `l1` norms first) at most `H(4N+4)`: intrinsic mass is at most `N`,
diagonal rate at most `N+1/2`, all old copies together at most `N/2`,
and a direct newborn rate at most one. The determinant expansion is
bounded by the product of these row sums, so its coefficient height is
at most `[H(4N+4)]^M`. The same bound applies with the triggered column
removed. The matrix determinant lemma gives their quotient `1−T_N` on
the domain where the killed trigger resolvent exists. This is a valid
literal height bound, **not** a lower bound on primitive height or a
useful resultant certificate at the observed size-error scale.

Even function transcendence plus full matrix positivity, a simple
crossing, and a positive trigger cannot settle special-value arithmetic
as abstract hypotheses. For example take rational `r>0` and
`G(s)=Σ_(j>=1) 2^(-j)(r+j)/(s+j)`. This is a positive Laplace response
with rational finite histories, and extends meromorphically with infinitely
many genuine poles `−j`, hence is transcendental over `C(s)`. Yet `G(r)=1`
exactly. Take the positive 19-by-19 matrix whose entries are `G(s)/19`.
Its radius is `G(s)` and its trigger scalar is `G(s)/(19−18 G(s))`
near the crossing, also transcendental since the Möbius map is invertible.
The root is exactly rational. This countermodel is independent of all
friend-tree assumptions and does not claim the physical entry functions
have that form.

## Certified finite exclusions and replay

Use exactly the reviewed fractions from
[`occupation-affine-verified.txt`](../../../../../astra/source/round2/outputs/occupation-affine-verified.txt),
lines 30 and 34. The exclusion conclusions are **conditional only on
that physical enclosure**; the replay checks arithmetic and does not
rebuild the enclosure. `43/76` is the first rational in the interval
ordered by denominator. For irreducible primitive quadratic polynomials,
the smallest height is 7, witnessed by `3X²−7X+3`; this does not exclude
quadratic irrationals without a height constraint. For rational square
discriminant, parameterize `s=t−1/t` with reduced rational `t=u/v>1`.
The reduced denominator is exactly `uv`. Enumeration by `uv<=57200`
shows the first possible ratio is `32361/57200`; Long's illustrative
`2603/4602` from the wider enclosure is strictly below the reviewed lower
endpoint. No wider square-discriminant class is excluded: it is dense.

```
timeout 8s python3 -B round2/open/arithmetic/check.py
```

The replay checks exact exclusions, the 19-step reset cycle against the
actual capped generator, four omitted histories (including the optional
return to label two), and the rational-crossing countermodel. The analytic
proof, not a finite replay at selected negative points, establishes all
pole orders needed above. None of these statements settles irrationality
of the original exponent.
