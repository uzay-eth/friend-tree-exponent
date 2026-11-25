> **Public checkpoint export note:** Edited export of `mixed-astra/ELASTICITY-MATCH.md` (frozen original SHA-256 `6ae3a3be86530a53831c922ae3e6cfaa194e4004d7d256968b12c24155aabe26`). Changes are limited to GitHub-compatible GFM math formatting (`` $`…`$ `` inline and fenced `math` blocks; `\tag` rendered as `\qquad (...)`, `\operatorname`/`\hbox` normalized) and local path/session/link redaction. All mathematical content, formulas, labels, errata, and scoped audit verdicts are preserved. **Line numbers cited in errata and audits refer to the frozen original source bytes listed in [`../PROVENANCE.md`](../PROVENANCE.md), not line numbers in this reformatted export.**

# A sharper stopping inequality for the physical renewal mean

For the actual tree, $`\beta(p)=0`$ for every fixed $`927/1000\le p<1`$. The earlier reviewed upper interval began at $`187/200`$. The lower interval $`\beta(p)>0`$ for $`p/q<R`$, $`q=1-p`$, remains as in CRITICAL.md, where

```math
B(a)=\pi(1+a)/\sin(\pi a),\qquad R=\min_{0 < a < 1} B(a).
```

The same function $`B`$ appears on the upper side by preserving the age distribution in the stopped convolution. At five sixths, strict physical supercriticality requires the finite-age ratio $`y(t)k(t)/q`$ to exceed $`5/R`$ somewhere. A sign at five sixths still requires physical input. These bounds hold pointwise on their stated parameter intervals.

## Stochastic ordering at a first zero

Use original-clock time. The physical virgin kernel is $`k(t)=q\,\mathbb E[1/D(t)]`$. Put $`Y=y-1`$, where $`Y(0)=0`$ and $`v=Y'=p+k*v`$. Then $`Y>0`$ off zero. Fix $`0<a<1`$. If $`H(t)=tv(t)-(1+a)Y(t)`$ has a first positive zero $`T`$, then

```math
a pT\le\int_0^T (T-r) k(T-r)\,dY(r),\qquad
Y(r)/Y(T)\ge (r/T)^{1+a}\quad(0\le r\le T). \qquad (1)
```

Indeed $`H=-apt+k*H+(sk(s))*v`$ by subtracting the integrated renewal equation. At small positive ages $`H=-apt+O(t^2)`$. Up to the first zero $`Y(r)/r^{1+a}`$ decreases. In particular the probability measure $`dY(r)/Y(T)`$ is stochastically smaller than the measure on $`[0,T]`$ having cumulative $`(r/T)^{1+a}`$. For every nondecreasing nonnegative test function, its integral against the former is at most its integral against the latter. Truncation justifies integrable singular test functions.

If $`k(s)\le Cq/y(s)`$ up to $`T`$, then, strictly for $`s>0`$,

```math
s k(s)< C A s^{-a},\qquad A=qT^{1+a}/Y(T).
```

Apply (1) to the increasing function $`(T-r)^{-a}`$. The resulting integral against the power density uses

```math
\int_0^T s^{-a}(T-s)^a ds=T\,\Gamma(1-a)\Gamma(1+a)
=T\pi a/\sin(\pi a).
```

It follows that $`p/q<C B(a)`$ is necessary at that first zero. If no first zero exists, $`Y(t)=O(t^{1+a})`$ and $`v(t)=O(t^a)`$ at infinity. The Laplace transform of $`v=p+k*v`$ then gives $`\widehat k(s)<1`$ for every $`s>0`$. Since physical $`\int k\ge1`$, its mass equals one and $`\beta=0`$. In particular if $`\beta(5/6)>0`$, choose a minimizer of $`B`$; some finite age satisfies $`y(t)k(t)/q>5/R`$. This gives an explicit size of the discrepancy that a sufficient upper comparison would need to control.

**Exact limit of the scalar comparison.** For the deterministic equality $`k=q/y`$ used in the reviewed HIGH-P.md, the preceding test gives $`\int k=1`$ for every $`p/q\ge R`$, including equality. Indeed $`\int k<1`$ would bound $`v`$ by $`p/(1-\int k)`$ and make $`q/y`$ nonintegrable. Conversely the necessary condition from CRITICAL.md gives $`\int k>1`$ for $`p/q<R`$: its proof uses only $`k\ge q/y`$ and $`\widehat y(s)=1/s+p/[s^2(1-\widehat k(s))]`$ under unit mass, both exact here. Thus the scalar comparison already allows unit mass at the specific ratio 5, since $`B(1/2)=3\pi/2<5`$. The physical implication at that parameter is the ratio test above.

For $`p/q>R`$ the scalar equality even has $`Y(t)\sim c_p t^{1+a_-}`$ for a finite positive constant, where $`a_-`$ is the smaller root of $`B(a)=p/q`$. Use $`a=a_-`$ in (1), so $`Y(t)/t^{1+a_-}`$ is decreasing. It is bounded below away from zero by the lower-root estimate (3) of ZERO-SCALE.md with $`W=y`$ and $`d=1+a_-`$. That estimate needs a probability renewal kernel, $`W=1+pJ`$, $`k\ge q/W`$, and its renewal convolution identity, all holding here. The constant follows by monotone convergence of the ratio. This asymptotic concerns the equality comparison alone.

## Keeping the early physical bound

The early interval of the physical proof has a different kernel bound. Here are the exact hypotheses needed from Proposition 1, Steps 1–2 of Neon's pinned HIGH-P.md, SHA-256 f5379623485393a366290577dac25d39e5e4744a4ad48307420bfa3262c49c64. Use the fresh killed-family comparison as scoped in REVIEW-HIGH-P.md. If the first four parameter conditions there hold at $`p_0`$, then for every fixed $`p\in[p_0,1)`$ any putative first zero has $`T>t_1(p)`$, with

```math
Y(t_1)<p t_1/(1-i),\quad s k(s)\le q/p\ \ (s\le t_1),
\qquad k(s)\le C_2q/y(s)\ \ (t_1\le s\le T).
```

The constants $`i,C_2`$ are fixed at the certificate parameters. We use this already reviewed stopped estimate as an input. Its previous final contradiction is unnecessary.

Let $`C=\max((1-a)/(1-i),C_2)`$ and abbreviate $`t_1=t_1(p)`$. Define the decreasing majorant in $`s`$ by

```math
w(s)=\begin{cases} C t_1^{-a}/(1-a), & 0 < s \le t_1, \\ C s^{-a}, & t_1 < s \le T . \end{cases}
```

For $`s\le t_1`$, the decreasing $`Y(r)/r^{1+a}`$ gives $`q/p<A t_1^{-a}/(1-i)`$. For $`s>t_1`$, $`s k(s)<C A s^{-a}`$ as before. Thus $`s k(s)\le A w(s)`$ everywhere, with strict inequality on the positive-length late interval. Apply stochastic order in (1) to $`w(T-r)`$. When integrated against the power density, the early constant costs no more than $`C s^{-a}`$. In fact the decreasing functions $`(T-s)^a`$ and $`s^{-a}`$ on $`[0,t_1]`$ have nonnegative covariance under the uniform law there, and the mean of $`s^{-a}`$ is $`t_1^{-a}/(1-a)`$. Integrating the full power majorant as above again gives

```math
p/q<C\,B(a). \qquad (2)
```

Consequently the same first four conditions, with $`p_0/q_0\ge C B(a)`$, suffice uniformly on the entire parameter interval. The uniformity comes from the stopped bounds for each $`p`$ and the explicit scalar inequality $`p/q\ge p_0/q_0`$.

## Rational certificate for the physical interval

Take $`(p_0,a,i,\delta,\eta,\theta)=(927/1000,11/25,159/200,1/4,11/100,13/50)`$. For reproducibility the scalar checks are as follows. Write $`\rho=(1-p_0)/p_0`$, $`X=\exp(i/\rho-4/5)`$, $`b=\rho(1-\delta)/[(1-\eta)\delta]`$, $`c=\eta+(1-\eta)\log(1-\eta)`$. Set

```math
B_1=c\delta X,\quad
L=\frac{\theta^2 e^{-2b}X(1-i)(1-\delta)}
{2(1+a)^2[2+(i+17\rho/10)^2]},
```

```math
C_2=\frac{e^b}{(1-\theta)(1-\delta)^{1+a}}
+\frac{101 X e^{-B_1}}{100(1-i)}
+\frac{101 e^{-L}}{100(1-i)(1-\eta)\delta}.
```

The first four conditions are $`a\ge\rho/(1-i)`$, $`X>100`$, $`(1-\delta)X\ge1`$, $`i p_0>(1-p_0)^2`$, $`B_1>1+a`$, and $`L(1-a-\rho/\sqrt2)>a`$ with positive bracket. The script encloses exponentials and logarithms with the preserved rational Taylor primitives, including strictly positive upper endpoints for both tail terms. It uses $`\pi<22/7`$ and $`\sin(\pi a)=\cos(\pi(1/2-a))\ge 1-(22/7)^2(1/2-a)^2/2`$. The exact author replay output was

    $ python3 -B check_elasticity_match.py
    p=927/1000: B1 > 17, B2*f0 > 5.60, C < 2.751; both tail bounds positive
    weighted stopping cost < 12.68 < 927/73; uniform upper interval starts at 927/1000

This contradicts (2) throughout the new interval. Together with the strict lower interval this narrows the set of undecided physical parameters. The number and ordering of transitions inside that set remain open.
