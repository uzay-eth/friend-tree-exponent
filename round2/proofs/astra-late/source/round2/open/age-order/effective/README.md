> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results. See the [release index](../../../../../README.md#claim-18) for exact scope and citations.

# A computable late-label cutoff

**Numerical follow-up.** The [explicit cylinder certificate](cylinder/README.md) closes the numerical cutoff question below: it proves `πL>0.120996` and `p_(n+1)>p_(n+2)` for every `n≥2^(2^396)`. It uses a uniform bound on finite late hits and stationary moment bounds on excluded balls. The comparisons between the early prefix and that cutoff remain open.

The existing [adjacent-gap theorem](../adjacent/README.md) can be made effective for the original edge-seeded chain. Below is a terminating rational certificate search for a starting label. It does not require a uniform time-to-hub estimate over all environments. The search and the resulting bound are designed for computability; no numerical starting label has been evaluated. The separate [early certificate](../early/further/README.md) covers through label 11.

Write `l=113/200`, `u=1133/2000`, `rho=51/100`, `q=9/25`, `A=20/3`, and let `a_n=E C_n`. In the notation of the adjacent-gap theorem,

```math
n(n+1)(p_{n+1}-p_{n+2})=E\sum_{o\text{ core}}G(T_n,o),\qquad
|G|\le3,\qquad \pi G=(1-\mu)\pi F\ge \beta:=(1-u)l/2>0.
```

The bound on the stationary mean uses the [birth-to-promotion balance](../late/README.md) and the identification of the two limiting coefficients. It does not replace the deterministic prebirth law by the uniform-core law.

## 1. Enclose the hub score by finite cylinders

For a placed tag of degree one or two, take an integer `d≥3` and a time `T`. A lower endpoint is `(1−1/d) P(hit d by T)`. An upper endpoint is the expectation of one on the hit event and, otherwise, the capped sum of the neighbours' first-send bounds at time `T`. For a core neighbour with `j` other leaves one can use `min(1,2/(j+1))` from the [first-send barrier](../early/README.md); for a leaf neighbour use one. The upper endpoint counts the first future gain after `T`, on the unchanged neighbour set up to that gain. The lower endpoint uses a fresh tracked leaf at the hit. At fixed `d`, the asymptotic interval width is at most `1/d`: on nonhub paths with no hit, the upper charge tends to zero, and on hub paths the hit occurs eventually. In fact the limiting upper endpoint is `P(ever hit d)` and the limiting lower endpoint is `(1−1/d)` times that probability. This argument applies to the degree-two middle as well as tips.

Localize all such tests to an old core ball of radius `R≥6` about the original score root, retaining the physical boundary degrees and making the exterior slots inert. Insert the two prescribed vertices before running the clocks. One can work with the post-insertion cells: the tag and the initial neighbours whose stocks can be queried are contained in cells based within distance five of the original root (a leaf starts in the cell of its core parent). Additional neighbours of the tag arise from its own cell or the cell of its initially adjacent core; in particular an initially adjacent leaf promoted before the test stays in its cell. This uses that the only edges between cells are the old core edges. For any one placement, if `V` is the sum of the initial core degrees in that post-insertion ball of radius five, the [ordered-proposal coupling](../../../../../../astra/source/round2/operator.md) bounds the failure probability by

```math
P_R=\min(1,V T^{R-5}/(R-5)!).
```

Indeed every omitted old core is at distance at least `R+1` from the original root; by the time an influence reaches a queried cell it has crossed an ordered path with a suffix of length at least `R−5`. Reverse the suffix as in the coupling proof and union bound over queried cells. Use the unconditioned truncated endpoint expectations with additive `±P_R`; do not condition either process on successful coupling. At fixed input and `T` this error vanishes as `R` increases.

Impose a total-degree cap `D` throughout the input core ball, reading the boundary core degrees as well as leaf marks. There are finitely many admitted balls. In each truncated system with `M` initial ordinary vertices, the genuine change rate, excluding null slots, is at most its current size. Stop at size `B`; the probability cost by time `T` is at most `M e^T/B`. Finite killed chains and the hit flag are evaluated with outward rational uniformization. Sum the endpoint intervals with their original signed coefficients. The sum of the absolute coefficients on **differences** is at most three per core, exactly by the charge calculation in the adjacent theorem. Intersect with `[-3,3]`. On excluded balls use `[-3,3]`; denote the resulting cylinder endpoints by `L,U`.

These intervals admit arbitrarily small `pi(U−L)`. First choose `d` large, then `T` large using bounded convergence of the endpoint widths under `pi` and the coefficient bound (a difference uses two endpoints). For that `T` choose `R` large by the same argument, then `D` large at that fixed radius, and finally the size and arithmetic cutoffs. This is an existence proof for success of a dovetailed **finite rational test**, not a request for an oracle for any of those limiting integrals. When a jump estimate is used below, write `L=L_0−3` and `U=U_0+3`, where `0≤L_0≤6` and `−6≤U_0≤0` are supported on the admitted balls.

There are also explicit core moment cutoffs. From `pi Q_t=e^{mu t}pi`, start at `pi`, ignore the single old root after discounting (first truncate the test at a fixed value), and sum over promotions in its cell. The promotion intensity in expectation is at most `e^{mu s}`. A new core starts at total degree two and core degree one, and has age-`v` moments bounded by `2^{1/4} e^{v/4}` for `D_o^{1/4}` and by `(5/4)e^{v/4}` for `(5/4)^{c_o}`. Thus `pi D_o^{1/4}≤2/(l−1/4)` and `pi (5/4)^{c_o}≤(5/4)/(l−1/4)`. For example `c^2≤40(5/4)^c` for nonnegative integer `c`: compare the binomial terms of orders zero, one and two. Hence `pi(c_o 1_{D_o>J})≤B_tail J^{-1/8}` with `B_tail=10/(l−1/4)`. Reverse a simple path by [mass transport](../../../../../../astra/source/round2/operator.md), using successive mark cutoffs `J_i` at depth `i`. The union bound is the root moment tail plus `sum_{i=1}^R (J_1…J_{i−1}) B_tail J_i^{-1/8}`. One may choose the `J_i` successively, then use their maximum as `D`. At deterministic size, excluded-ball mass under the **expected count divided by `a_n`** is also bounded using §3 on the exclusion indicator. No higher uniform-in-size degree moment of the physical tree is assumed.

## 2. Certify the spatial semigroup constants

The required estimate is [spatial-limit (1)](../../../../../../astra/source/round2/fringe/spatial-limit.md), with a certified `C, sigma`:

```math
|e^{-\mu t}Q_t f-h\pi f|\le C4^k e^{-(\mu-\sigma)t}\|f\|_\infty .
```

Here `k` includes an extra level if needed to read boundary core degrees. The following winding-number search makes the spectral step in the [renewal proof](../../../../../../astra/source/round2/fringe/almost-sure.md) effective. All matrix norms below are maximum absolute row sum.

* For the nineteen-label kernel `k(t)=J R_t B_0` (here `B_0` is the renewal birth operator), `||k(t)||≤A e^{rho t}` and `||k'(t)||≤64A e^{rho t}` almost everywhere. Indeed `||R B_0 g||_W≤64` for `||g||≤1`: root gains below the mark cap cost at most 40; above the cap the gain increment of `L/D` costs at most one after multiplying by the rate; root promotions cost at most two; killing costs at most one; and the child-weight identity bounds the absolute copy contribution by `(1/2+d_weight)W`, with `d_weight=3/20`. Thus on `Re z≥sigma>rho`, `||K(z)||≤A/(sigma−rho)` and `||K(z)||≤B_1/|z|`, where `B_1=A+64A/(sigma−rho)`. Also `||K'(z)||≤A/(sigma−rho)^2`.
* Enclose complex `K(z)` with the finite resolvent of [finite-kernel (3)](../../../../../../astra/source/round2/fringe/finite-kernel.md), using the same physical root-slot correction `r_0=n−1/D_root`. Entrywise errors are dominated by the real first-exit bound at `Re z`. Uniformly in the half-strip one can choose a positive rational `delta≤1` with `rho+delta(1+q/2)≤sigma`; the error is at most `A (2/(N+1))^delta/(sigma−rho)` per entry. This follows by taking the absolute value of the post-exit reward before integration. No complex positivity is used.
* Search rational `sigma∈(rho,mu)` using real enclosures of `mu`. Take an integer `H>2 B_1+1`; use the rectangle with real sides `sigma,H` and imaginary sides `−H,H`. Enclose `det(I−K)` on subdivisions of its perimeter using the resolvent error and derivative bound. Certify nonvanishing and winding one (e.g. compare with a rational polygon contained segment by segment in zero-free disks). Also enclose `||(I−K)^{-1}||` there by adjugate and determinant bounds; use a rational upper bound `M_*`. A successful subdivision exists for some `sigma`, by the simple Perron pole and the nonlattice argument in the renewal proof. Outside this rectangle in the half-plane `Re z≥sigma`, the Neumann bound is at most two. Winding one and the known real root certify the whole strip. Search also for a rational `eta>0` with `eta≤mu−sigma` using real brackets.
* The residue norm has a rational upper bound `D_*=(4H) M_*` by the rectangular residue integral. Inverting `(I−K)^{-1}K^2` on the vertical line gives the renewal density remainder `||v(t)−e^{mu t}D_res||≤C_v e^{sigma t}`, with `C_v=A+2H M_*(A/(sigma−rho))^2+4B_1^2/H`. Here `D_res` is the residue itself, and the leading `A` accounts for `k(t)`. The constants drop the harmless factor `1/(2 pi)`.

For reproduction from a newborn the convolution remainder constant is bounded by `C_J=A+A C_v/(sigma−rho)+A D_* /(l−rho)`. For arbitrary planted types the weighted constant in renewal estimate (3) is bounded by `C_P=1+C_J/(sigma−rho)+A D_* /(l−rho)^2`. These follow by integrating `e^{rho s}` against the remainder and separately bounding the missing tail of the leading convolution; the discounted forcing from newborns has norm at most `A/(l−rho)`. The weighted norm of the rank-one leading operator for arbitrary planted types is at most `A D_* /(l−rho)^2`. Inserting the immigration bound `A+2e^{qt}` and the shallow-layer error of [spatial-limit §2](../../../../../../astra/source/round2/fringe/spatial-limit.md) gives, for example,

```
C = 1 + C_P (A/sigma + 2/(sigma-q))
      + A D_* /(l-rho)^2 (A/l + 2/(l-q)).
```

To evaluate `pi f` for any fixed bounded cylinder, start at the single core of `T3`, compute `Q_t f` and `Q_t1` with a size cutoff, then divide. The neglected expected population when the total Yule size from three reaches `B` is at most `12 e^{2t}/B`, since on cutoff the final size is at least `B`. The semigroup estimate gives the interval by separately enclosing `e^{-mu t}Q_t f` and `e^{-mu t}Q_t1`; divide the resulting leading-coefficient intervals once the latter has positive lower endpoint. It has limit `h(T3)>0`, so this also supplies a rational `h_*>0` below `h(T3)`. Exponentials and powers in the procedure can be enclosed to any precision with rational arithmetic. In particular there is an accepted candidate from §1 with `pi L>beta/2`; test this strict inequality by dovetailing finite interval evaluations. One can also certify `pi(U−L)<epsilon` for any requested positive rational `epsilon`.

## 3. Return to the deterministic original labels

Let `f` be supported on capped radius-`k` balls with maximum degree cap `D`, or differ from such a function `f0` by a constant `b0`. For the sum over cores a valid one-birth jump bound is

```math
J_f=16(k+1)(D+1)^{k+1}\|f_0\|_\infty+|b_0|.
```

Indeed only the gain vertex changes a mark, or in a promotion the old parent and the newly promoted core affect balls. For a root whose summand is nonzero on either side, charge in the relevant graph to a changed vertex in its ball; when the new core was absent, use the old parent for that side. Reverse paths within a capped ball. There are at most `sum_{i=0}^k D^i` roots for each charged vertex, per side, and a summand changes by at most `2||f0||`. This counts disappearance of a ball from the cap as well as entry. The constant summand changes only on a promotion.

For `n` large use `m=floor(n^(19/20))≥3`, `t=log(n/m)`. From the exact tree `T_m` run an independent Yule size clock and the original jump chain for elapsed time `t`. Its size `N` has mean `n`, variance `n(n−m)/m`, and is independent of the entire discrete graph sequence. Couple the evaluation with size `n` on that same sequence. Then

```math
\left|E\sum_{o\in T_n}f(T_n,o)- E\sum_{o\in T_m}Q_t f(T_m,o)\right|\le J_f n/\sqrt m.
```

For a ratio, apply this to `f−pi f`, whose semigroup leading term vanishes. If `J` bounds its jumps and `M` bounds its sup norm, this yields the fixed-label error

```math
\left|a_n^{-1} E\sum_{o\in T_n} f-\pi f\right|
\le C4^k M n^{-\eta/20}\frac{e^{\mu t}a_m}{a_n}+ \frac{J n}{\sqrt m\,a_n}.
\qquad\text{(*)}
```

Both remaining ratios have explicit bounds. Use the clock from `T3`: at elapsed time `log n` the size has mean `3n`, and `a_j≤max(1,j/n)a_n` for `j≥3` by monotonicity and the promotion bound `a_{j+1}≤(1+1/j)a_j`. Hence `Q_{log n}1(T3)≤4a_n`. Conversely the size at elapsed time `log n` from `T3` reaches `n` with probability at least `1/4`: one of its three independent Yule families has geometric size with tail `(1−1/n)^(n−1)≥e^−1>1/4`. Hence, from the certified semigroup bound,

```
a_n ≥ h_* n^mu/8         if C n^(-eta) ≤ h_*/2,
a_m ≤ 4 (H_*+C) m^mu     where H_* ≥ h(T3).
```

Consequently `e^(mu t)a_m/a_n≤32(H_*+C)/h_*`, while `n/(sqrt(m) a_n)≤16 n^{-(l−21/40)}/h_*` using `m≥ n^(19/20)/2`. The exponent `l−21/40=1/25` is strictly positive. Neither the clock independence nor monotonicity of the **core count** gives monotonicity of the signed score; that is why the bounded jump estimate is needed.

Use `f=L` in (*) with `M=6` and a bound `J` for `L−pi L` (for instance add three to the jump bound for `L`). Pick an integer `N_0` beyond which the two decreasing rationally bounded error terms sum to less than `beta/4`, with the preceding core lower bound valid there. Then all `n≥N_0` satisfy

```math
p_{n+1}-p_{n+2}\ \ge\ \frac{\beta h_*}{32}\frac{n^\mu}{n(n+1)}>0.
```

For a two-sided normalized error, apply (*) to both `L` and `U` and add the certified width `pi(U−L)`. The monotone analytic envelopes above may be verified at a single integer `N_0` by rational bounds on powers. This certificate search supplies a computable cutoff, not a proof of positivity for the intervening labels. The general finite pullback can certify a remaining comparison if its gap is strict, by the two endpoint tests, but cannot certify equality by shrinking intervals alone.

`threshold.py` implements the final arithmetic given certified inputs. It returns the cutoff as a power of two without expanding the integer, which can be enormous even for optimistic input constants. Its self-check verifies the arithmetic inequalities and tests rejection of a cutoff whose clock error is too large. It does not implement the matrix winding search or the finite-cylinder evaluation.
