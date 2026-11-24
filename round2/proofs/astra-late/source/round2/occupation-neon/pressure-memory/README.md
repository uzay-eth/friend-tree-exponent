> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Retaining exterior donor degrees in stationary motif calculations

In the original friend tree, the same old child can send to its parent many times at its current reciprocal degree. Only an arrow received by that child decreases its reciprocal. Following this distinction gives an exact family of stationary constraints and an explicit error bound for replacing high donor degrees by a finite collection of joint moments.

For the four-vertex path with censuses `(L,k,1)=(0,1,1),(0,3,1),(0,1,1),(3,0,1)`, 150 polynomial moments and one tail count suffice to represent the entire exterior-donor contribution to the paired D18 contrast at multiplier one within **$`6\times10^{-9}`$ per uniformly sampled core**. The moments are weighted by incidences of this specific path in the actual tree. The result bounds the error of the representation. Computing the stationary moments and accounting for the other path censuses and motif families are further steps.

## The old-donor balance

Orient the core tree from the first core, and sample a core under its stationary law $`\pi`$. It has one parent in the limiting tree. For a core $`x`$, let $`k_x`$ count its core children, $`L_x`$ its leaf neighbours, $`p_x=1/D_x`$, $`b_x=L_x/D_x`$, and $`g_x=L_x+p_{\mathrm{parent}}+\sum_{y\text{ child}}p_y`$. Put
```math
 R_x=\sum_{y\text{ child}}p_y,\qquad B_x=k_x/2-R_x,\qquad
 M_j(x)=\#\{y\text{ child}:D_y>j\}\quad(j\ge2).
```
The degree of an existing child never decreases. Each newly promoted child starts at degree two and contributes zero to $`B_x`$ and to every $`M_j(x)`$. A gain at $`y`$ increments $`B_x`$ by $`e_y=1/[D_y(D_y+1)]`$, at its full physical rate $`g_y`$. In particular a send from $`y`$ to $`x`$ changes $`D_x`$ alone. For a census test $`h(D,k)`$, define $`A_xh=g_x\Delta_D h+b_x\Delta_k h`$, where the second difference is at fixed degree. A single physical birth cannot change both the old receiver census and its deficit. Conditional on the full current tree, the generator on the old root therefore gives
```math
 L(h f(B_x))=f(B_x) A_x h+
 h\sum_y g_y\{f(B_x+e_y)-f(B_x)\}.
```
For the time derivative of a sum over all cores one includes the separate newborn term $`b_x h(2,0)f(0)`$. In particular
```math
 \pi[(\mu h-A_x h)B_x]=\pi\!\left[h\sum_y g_y e_y\right],\qquad
 \pi[(\mu h-A_x h)M_j(x)]=\pi\!\left[h\sum_{D_y=j}g_y\right].
```
One may first take compact degree support for $`h`$. The birth-age formula and discounted Dynkin give the equations, and the arity has a uniform exponential moment for any $`t`$ with $`e^t-1<\mu`$. The last bound follows by applying the generator to the sum of $`e^{t k}`$ and using $`b\le1`$. Thus all the displayed expressions for finitely supported $`h`$ are integrable. The same proof applies to $`h=1`$, where the census differences vanish. For other tests, integrability of their gain differences and the boundary terms must be checked. [Stable age](../renewal/README.md) · [Oriented limit and arity moments](../flow/README.md)

The threshold sources use the same marked parent–child edge as their tests. Within a retained census, $`g_x B_x=(L_x+p_{\mathrm{parent}}+k_x/2)B_x-B_x^2`$, which supplies a second moment constraint. If only $`R_x`$ and $`k_x>0`$ are known, convexity of $`p^2/(1+p)`$, concavity of $`p/(1+p)`$, and $`p_x\le g_y\le D_y`$ give
```math
 \frac{p_x R_x^2}{k_x+R_x}\le \sum_y g_y e_y
 \le \frac{k_xR_x}{k_x+R_x}.
```
For $`k_x=0`$ all terms vanish. Revealing the degree of an old donor makes its jump size exact. Revealing its other senders further constrains its rate. A fixed tagged motif uses the same identities between incidence births; its entry and exit measures must also be included. [Labelled generator and compensated births](../../../../../../reviews/cold/astra/checks/pressure-memory-sol.md)

## An effective joint degree and child-count tail

In the physical birth-Palm law a core starts at $`(D,k)=(2,0)`$, in a random live environment. The law $`\pi`$ observes it after an independent exponential age of rate $`\mu`$. Let $`Y`$ be a pure-birth process starting at two, increasing from $`d`$ at rate $`d-1/2`$. Let $`N`$ be a unit-rate Poisson process independent of $`Y`$. Uniformly over birth environments, there is a coupling
```math
 D(t)\le Y(t),\qquad k(t)\le N(t).
```
Indeed the true gain rate is at most $`L+(k+1)/2\le D-1/2`$. Use the jumps of $`Y`$ as opportunities for a true incoming arrow, accepting with probability $`g/(Y-1/2)`$ and selecting its sender at its physical conditional rate. Independently supply every ring of the tracked core from $`N`$. At those rings sample the true uniform neighbour. All other arrows have their original conditional rates. An incoming arrow and the tracked core’s own ring are disjoint transition classes. This constructs independent *dominating* processes at fixed elapsed time in the adaptive exterior. Mixing over the common observation age preserves the joint upper comparison.

Use any rational $`0<s\le\mu`$; the code takes $`s=565673/10^6`$ from the original double-sided D18 verification. Couple the physical age below $`\mathrm{Exp}(s)`$. For $`d\ge2`$ put
```math
 P_d=\prod_{j=2}^d\frac{j-1/2}{j-1/2+s},\qquad
 H_d=1/s+\sum_{j=2}^d\frac1{j-1/2+s}.
```
Then
```math
 \pi(D>d)\le P_d,\qquad \pi[k\mathbf1_{D>d}]\le P_d H_d.
```
To see the second bound, write $`\tau`$ for the hitting time of $`d+1`$ by $`Y`$. Conditionally on observation at age $`t`$, $`N(t)`$ has mean $`t`$ and is independent of $`\tau`$. The integral from $`\tau`$ to infinity of $`s t e^{-st}`$ is $`(\tau+1/s)e^{-s\tau}`$. Exponential tilting of the independent holding times increases each rate by $`s`$, giving the stated expression. The same calculation yields a bound for factorial arity moments of order $`m`$ by integrating $`s t^m e^{-st}`$ from $`\tau`$ to infinity. [Original input and verification](../../../../../../reviews/cold/astra/attacks/neon1/REHEARSAL.md)

The exact layer expansion is $`B_x=\sum_{d\ge2} M_d(x)/[d(d+1)]`$. Truncating at $`J`$ replaces each reciprocal of degree above $`J`$ by $`1/(J+1)`$. The per-donor excess is $`w_J(D)=(1/(J+1)-1/D)_+`$. Put $`i=J+1`$ and $`a_i=1+si/(i+1/2+s)`$. The unweighted and arity-weighted errors satisfy
```math
 \pi w_J\le \frac{P_i}{i a_i},\qquad
 \pi[k w_J]\le \frac{P_i H_i}{i(a_i-1/H_i)}
 \quad(a_i>1/H_i).
```
For the first inequality telescope $`P_d/d-P_{d+1}/(d+1)=P_d a_d/[d(d+1)]`$ and use monotonicity of $`a_d`$. For the second use $`H_{d+1}-H_d=(d+1/2+s)^{-1}`$. Multiplying the corresponding difference by $`d(d+1)/P_d`$ gives at least $`a_d H_d-1`$, which is at least $`(a_i-1/H_i)H_d`$ for $`d\ge i`$. Both majorants decay, respectively, as $`O(J^{-1-s})`$ and $`O(J^{-1-s}\log J)`$. For bounded weights on a child occurrence, reversing the physical parent–child sum uses the unweighted donor marginal. For a parent that can serve several children, the arity-weighted estimate retains that multiplicity.

## A physical D18 path with strongly donor-dependent reward

Write $`A,B,C,D`$ in order down the four-vertex path whose states are specified above. For this occurrence the lower-minus-upper charged row has exterior coefficients $`c_p=0.00588543562`$ for the parent of $`A`$, and $`c_c=0.00701260718`$ for the sum of the two other child reciprocals at $`B`$. They are the exact differences of the two original payload slope coefficients. A paired physical fixture at these labels changes the ratio of the two slack scores from less than one to over eleven million when only the exterior donors change. Exact row evaluator (`round2/occupation-neon/all-contrast/scan.py`, historical reference) · [Independent labelled checks](../../../../../../reviews/cold/astra/checks/stationary-p4-sol.md)

Fix a hidden side child $`y`$ of $`B`$. Its parent and grandparent determine $`B,A`$, and at most two of its siblings can be the path’s $`C`$. Each $`C`$ has a unique core child. A donor serving as the parent of $`A`$ can be used by at most three paths per core child of that donor. Both multiplicities are attained. Reversing the sums on finite trees and then passing to the uniform-core law proves that the stationary cost of truncating these two reciprocals is at most
```math
 2 |c_c|\pi w_J+3|c_p|\pi[k w_J].
```
There is no conditional stationarity assumption about the donor given the particular path. Arity integrability justifies the reversal.

For $`2\le d\le J=65536`$, partition the degrees into 15 disjoint dyadic slabs with closures $`[a,2a]`$. For order $`m`$, on each slab use the polynomial
```math
 Q_{m-1}(d)=\frac{1-T_m(3-2d/a)/T_m(3)}{d}.
```
Here $`T_m`$ is the ordinary Chebyshev polynomial. The numerator vanishes at zero so the quotient is a polynomial in $`d`$ of degree $`m-1`$. Throughout the slab the relative reciprocal error is at most $`1/T_m(3)`$. Use the constant $`1/(J+1)`$ above the last slab. For $`m=10`$, $`T_{10}(3)=22619537`$. The low-degree approximation error is at most $`(|c_c|+3|c_p|(1-s))/T_{10}(3)`$. Indeed each child incidence has been charged at most twice, $`p\le1/2`$, and $`\pi(k/D)\le1-\pi b=1-\mu`$. One may take the moments of the joint positive incidence measure with weight $`c_p`$ for an exterior parent and $`c_c`$ for an exterior child. This needs 10 moments on each slab and the weight of the remaining tail. The constant part of the charged row multiplies its ordinary occurrence mass separately.

| High-degree cutoff | Absolute tail charge per core, rounded up |
|---:|---:|
| 1,024 | $`2.463\times10^{-6}`$ |
| 8,192 | $`1.153\times10^{-7}`$ |
| 32,768 | $`1.470\times10^{-8}`$ |
| 65,536 | $`5.225\times10^{-9}`$ |

Including the polynomial error at cutoff 65,536 gives a total below $`5.8740\times10^{-9}`$. This budget concerns approximation of a running reward under the physical stationary law. Changing the process’s jump rates to match a polynomial approximation would require its own error argument.

Run `python3 -B round2/occupation-neon/pressure-memory/check.py` from the repository root. The standard-library replay checks the scalar majorants, the path multiplicities, and polynomial degree and relative error. It computes the tail products and tilted means with directed rounding at $`10^{-24}`$ precision, and checks the final inequality in rational arithmetic.
