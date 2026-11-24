> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes; Claims 6, 10, 11, and 18. See the [release index](../../../../../README.md#claim-22) for exact scope and citations.

# Birth labels and the small-mass hub spectrum

**Result.** Write `K` for the hub-spectrum constant (the `A` of [Claim 11 / hub spectrum](../../../../../../astra/source/round2/fringe/hub-spectrum.md)), to distinguish it from the scaled endpoint tail `A(y)` of [scaled/](../scaled/README.md). Let `X=lim C_n/n^μ`. Then

```math
 \boxed{K\,\mathbb E X=\int_0^\infty y^{\mu-1}A(y)\,dy.}
```

The argument below derives the joint law for **annealed count-weighted sampling**. The [quenched proof](quenched/README.md), independently [cold-audited](../../../../../../../reviews/cold/astra/checks/marked-spectrum-sol.md), gives the same joint limit for uniform sampling within almost every realized tree, without count-biasing the tree law. The lead supplied the degree-two Pólya comparison after first promotion and the promotion-time truncation; these give the domination and integrability needed to match Claim 11 on the original discrete law.

## 1. Rare-label envelope and the Mellin normalization

Let `M_v` be the discrete size at which vertex `v` first promotes, with `∞` for no promotion. For `v≥3`, on `{M_v=m}` its degree is two and the conditional tail is bounded by the `Beta(2,m−2)` tail. Indeed the successive one-step degree gains are dominated by the two-colour Pólya count started from degree two at *that* size, because the physical attachment rate is at most the degree. Thus for `0<x<1`,

```math
 \Pr(Z_v>x\mid\mathcal F_{M_v})\le (1-x)^{m-2}[1+(m-2)x] \quad(M_v=m<\infty).
```

The right side decreases with `m` (subtract consecutive terms, obtaining the nonnegative difference `(m−1)x²(1−x)^(m−2)`). Since `M_v≥v` and positive mass requires finite promotion, the [late-birth bound](../late/README.md) gives

```math
 \Pr(Z_v>x)\le \frac{7\mathbb E C_{v-1}}{v-1}\, (1-x)^{v-2}[1+(v-2)x]
 \le C v^{\mu-1}(1+vx)e^{-vx} . \qquad\text{(1)}
```

After bounding finitely many small sizes, this holds with a single constant for all `v≥3`. For fixed `y≥0`, `v^(1−μ)Pr(vZ_v>y)→A(y)` uniformly on compact `y`-sets by [scaled/](../scaled/README.md). Applying (1) at `x=y/v` shows `A(y)≤C(1+y)e^(−y)`. Put `N_x=# {v:Z_v>x}`. The original-law Riemann sum gives

```math
 x^\mu\mathbb E N_x\longrightarrow J:=\int_0^\infty y^{\mu-1}A(y)\,dy.
```

On `ε≤xv≤M` use the uniform finite-label limit; below `ε` the normalized expected count is `O(ε^μ)+O(x^μ)` by (1), and above `M` it is bounded by the tail of the integrable envelope `C y^(μ−1)(1+y)e^(−y)` (or a decreasing-envelope Riemann bound for large `M`). All divergences at zero are integrable since `μ>0`.

For uniform integrability as `x↓0`, first discard `v>M/x` by the same bound. Among the remaining labels, those promoted by size `⌊RM/x⌋` are bounded in count by `C_⌊RM/x⌋`, apart from a harmless two initial labels. After multiplying by `x^μ`, these truncated counts have bounded second moments for fixed `M,R` by the core-count `L²` bound. Labels `v≤M/x` with finite promotion later than `⌊RM/x⌋` have expected normalized count at most `C M^μ(1+log R)/√R + o(1)`, by the uniform late finite-promotion bound applied at size `⌊Rv⌋`; finitely many small labels cost `o(1)`. Choose `M` then `R` to make the discarded expected mass small. Hence the almost-sure Claim 11 limit `x^μ N_x→KX` also holds in `L¹`. The Claim 11 non-strict cutoff gives the strict version by squeezing with `(1+ε)x`; no no-atom assertion for quenched vertex masses is involved. Comparing expectations proves the boxed identity.

For a local version with no core amplitude, let `ν` be the finite measure on `(0,∞)` with tail `ν((y,∞))=A(y)/E X=π(F_y)`. It is the uniform-core birth intensity restricted to positive local endpoint mass: no finite total mass is asserted for the measure that includes all zero-mass births. Then

```math
 K=\frac{1}{\mu}\int_0^\infty u^\mu\,\nu(du).
```

Equivalently, if `U_0` has probability law `ν/ν((0,∞))` (deterministic late birth conditioned to be a hub in the limit), then `K=(c/(μ E X)) E[U_0^μ]`, with `c=A(0)`.

## 2. Joint law for a hub sampled from the annealed count

For a hub above `x`, put `Y=xv`, `R=Z_v/x` and `U=vZ_v=YR`. The annealed sampling law weights each tree by `N_x/E N_x` and then chooses a uniform one of its `N_x` hubs. Equivalently, the expectation of a test is the expected sum over all above-threshold hubs divided by `E N_x`. For `0≤a<b≤∞` and `r≥1`, the same Riemann argument, at threshold `rx` but with label coordinate still `xv`, gives

```math
 \Pr_{\mathrm{count},x}(a<Y\le b, R>r)\longrightarrow
 \frac{1}{K}\int_a^b y^{\mu-1}\nu((ry,\infty))\,dy. \qquad\text{(2)}
```

Equivalently the limiting joint law of `(Y,U)` is `K⁻¹ 1_(0<y<u) y^(μ−1)dy ν(du)`. This is a probability by the Mellin identity. Its boundary lines `Y=a` and `R=r` have zero mass by integration in `y`; both coordinate tails are tight by (2) and its normalized marginals. The change of variable `y=u/r` factorizes the joint law of `(U,R)`:

```math
 U\sim\frac{u^\mu\nu(du)}{\mu K},\qquad
 R\sim \mu r^{-1-\mu}\,dr\quad(r>1),\qquad U\perp R,\qquad Y=U/R.
```

In particular the scaled birth label has density `y^(μ−1) A(y)/(K E X)`. The Pareto marginal follows directly within almost every realized tree, since `N_{rx}/N_x→r^(−μ)` by Claim 11 and `X>0`. The [quenched marked-count limit](quenched/README.md) gives the full joint law within almost every realized tree: the factor `X` in the random intensity cancels upon normalization by the realized count. The sampling schemes still differ at the tree level: the annealed count-weighted tree law has a density that converges in total variation to the law with density `X/E X` relative to the original tree law, while uniform in-tree sampling leaves the tree law unchanged.

## 3. The adjacent profile and a second Mellin identity

For `y>0` the scaled telescoping identity becomes

```math
 y^{\mu-1}A(y)-z^{\mu-1}A(z)=\int_y^z u^{\mu-2}B(u)\,du,\qquad z>y.
```

Continuity of `B` therefore gives `A∈C¹(0,∞)` and `B(y)=(1−μ)A(y)−yA'(y)`. It follows that `ν` has the continuous density `−A'(u)/E X` on `(0,∞)`. Integration by parts first on a compact interval, then using continuity at zero and the exponential bound at infinity, gives

```math
 \int_0^\infty y^{\mu-1} B(y)\,dy = \int_0^\infty y^{\mu-1} A(y)\,dy=K\,\mathbb E X.
```

No differentiability at zero or exponential envelope for the individual finite-label adjacent differences is needed here.

## Exact arithmetic fixture

`taskset -c 13 nice -n 10 timeout 10s python3 -uB check.py` checked 475 degree-two beta-tail and monotonicity rows by direct polynomial integration. It also checked 20 exact rectangle/tail rows in an atomic Mellin fixture, including the Pareto exponent and its independence from the power-biased mark. The original-law Riemann limit, truncation and spectrum comparison above are analytic.
