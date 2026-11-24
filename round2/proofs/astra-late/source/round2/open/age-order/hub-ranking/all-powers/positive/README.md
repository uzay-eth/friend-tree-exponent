> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits, scoped Claude derivations and replay, and Gemini proof review and independent checks; prior two-birth identities and finite box certificates as cited. See the [release index](../../../../../../../README.md#claim-27) for exact scope and citations.

# Strict moment order of every positive order in the random friend tree

Start the tree from the edge with labels 1 and 2. Each new vertex attaches to a uniformly chosen neighbour of a uniformly chosen existing vertex. Write $`Z_i=\lim_{m\to\infty}D_i(m)/m`$ and let $`\mathcal F_n`$ be the observed history at size $`n`$. For every real $`p>0`$,
```math
  \mathbb E Z_1^p=\mathbb E Z_2^p>\mathbb E Z_3^p>\mathbb E Z_4^p>\cdots .
```
For every $`n\ge2`$ and every history up to that size, with the same constant $`c=9/5000`$ at every positive power,
```math
 \mathbb E[Z_{n+1}^{p}-Z_{n+2}^{p}\mid\mathcal F_n]
 \ \ge\ c\,\frac{\Gamma(1+p)\Gamma(n)}{\Gamma(n+p+2)}. \qquad\text{(1)}
```
At the founder seam, $`\mathbb E[Z_2^p-Z_3^p]\ge c/[(p+1)(p+2)]`$. The labels 1 and 2 have equal laws by symmetry of the edge. [Audit of the fractional range](../../../../../../../../../reviews/cold/astra/checks/preplacement-all-fractional-sol.md) · [Audit of the unbounded range](../../../../../../../../../reviews/cold/astra/checks/preplacement-power-order-sol.md)

## Fractional powers on every physical tree

After the next two births the size is $`n+2`$. Conditional on each completed tree, the true Yule size amplitude has law $`\mathrm{Gamma}(n+2,1)`$ independently of the future embedded tree. The physical degree amplitude $`W`$ is the product of the size amplitude and $`Z`$. Thus the two-placement identity for $`\mathbb E W^p/\Gamma(1+p)`$ multiplies its signed numerator by $`\Gamma(1+p)/[n(n+1)(n+2)_p]`$ to get the degree-share gap. Here $`(m)_p=\Gamma(m+p)/\Gamma(m)`$. One uses the common size $`n+2`$ before pairing the placements. At a later random promotion the remaining-time physical process retains its actual conditional law. [Placement identity](../../../adjacent/README.md) · [Moment normalization](../README.md)

| Power | Trees | Lower bound on the normalized numerator |
|---|---|---:|
| $`0<p\le473/512`$ | At least two old cores | $`33/1000`$ |
| $`473/512\le p\le1`$ | At least two old cores | $`1/40`$ |
| $`0<p\le1`$ | Seed, finite star profiles, unbounded star tail | $`>0.13`$ in each of the three cases |
| $`1\le p<\infty`$ | Every old tree | $`9/5000`$ |

A core has ordinary degree at least two. The first 473 boxes retain the adverse tip's actual leaf siblings at the common completion, include the serial host's degree in its favourable continuation at small powers, and charge both endpoints of each suppressed core path. The 39 boxes through one use the discounted stock harmonic and the edge ledger. These are closed rational enclosures in the exponent, with the first box used for strictly positive $`p`$. The two families agree in normalization at $`473/512`$ and at one. [473-box physical proof](../../mean-ledger/fractional/README.md) · [Fractional stock harmonic and final 39 boxes](../../fractional-middle/README.md)

On the seed the serial weight is one. For the other exceptional trees evaluate the signed centre plus its old leaves at zero core-neighbour count. Above the finite degree cutoff use the lower bound on the centre's serial contribution and the upper bound on its total old-leaf debit. The minima over all 512 boxes for the seed, the exact finite star profiles and the analytic high-degree line are respectively
```math
 \frac{406542520567}{2500000000000},\qquad
 \frac{2085176012479}{12000000000000},\qquad
 \frac{265041348327659}{2012097656250000}.
```
The last two minima occur in the box starting at $`473/512`$. For $`p\ge1`$ the earlier all-tree proof gives (1) directly. [Exact star replay](../../mean-ledger/fractional/star_bound.py) · [Powers from one to infinity](../README.md)

## The first strict pair

At size three label 3 is a leaf. With probability one half label 2 is the centre, and in the other half it is exchangeable with label 3. For $`0<p\le1`$ the centre's normalized stock harmonic from two leaves is at least $`(2+p)/e`$. For an exponential variable $`E`$ of mean one, $`\mathbb E\log E\ge-1`$, so $`\Gamma(1+p)=\mathbb E E^p\ge e^{-p}`$. Therefore the centre's raw physical moment is at least $`C_p=(2+p)e^{-1-p}`$. The physical leaf has both a hub-probability upper $`a=197/305`$ and a mean upper $`b=387/1000`$. Hölder's inequality gives the raw upper $`L_p=a^{1-p}b^p`$. [Founder bounds](../../small-powers/README.md) · [Actual degree-two tip race](../../tip-moments/README.md)

Since $`\log(C_p/L_p)`$ is concave, its minimum is attained at an endpoint. Even using $`e<11/4`$, both $`2/(ea)`$ and $`3/(e^2b)`$ exceed $`51/50`$. Hence $`C_p-L_p\ge b/50`$. Returning from size three to limiting shares multiplies the centre–leaf difference by $`1/[2(3)_p]`$. Since $`\Gamma(1+p)\le1`$ by log convexity for $`0\le p\le1`$, this is a gap of at least
```math
 \frac{b/50}{(p+1)(p+2)}\ \ge\ \frac{c}{(p+1)(p+2)}.
```
For $`p\ge1`$ the size-three comparison already gives $`1/[2(p+1)(p+2)]`$. [Size-three proof for larger powers](../README.md)

## Comparing more distant labels and mixtures of powers

For $`2\le m\le i-1`$, $`i<j`$, sum (1) after conditioning every intermediate pair on the same $`\mathcal F_m`$. The beta recurrence gives
```math
 \mathbb E[Z_i^p-Z_j^p\mid\mathcal F_m]
 \ \ge\ \frac{c}{p+1}
       \{\mathrm B(i-1,p+1)-\mathrm B(j-1,p+1)\}. \qquad\text{(2)}
```
Unconditionally the bound also holds for $`i=2`$, using the founder bound for the first term. The same formula holds on the history at an almost surely finite stopping size $`m\ge2`$ for two fixed positive label offsets, by partitioning on $`m`$.

Taking $`p\downarrow0`$ in (1) gives the hub indicators by dominated convergence with $`0^p=0`$. The direct hub theorem gives the larger constant $`371/12500`$ at this endpoint. Integrating (1) or (2) against any nonzero finite positive measure on $`(0,\infty)`$ gives strict order for the corresponding mixtures of powers. For example, for any $`a>0`$, the bounded function $`g_a(z)=1/(a-\log z)`$, extended by $`g_a(0)=0`$, satisfies the same strict expectation order because $`g_a(z)=\int_0^\infty e^{-ap}z^p\,dp`$. Equivalently, the extended variable $`-\log Z`$, with value $`+\infty`$ on nonhubs, is ordered by its Laplace transforms. Strict first-order stochastic dominance at every positive threshold is a separate unresolved question. [Direct hub-probability ranking](../../README.md) · [Limits of moment information for sharp thresholds](../../beta-tests/README.md)

From the Astra root run `python3 -IB round2/open/age-order/hub-ranking/all-powers/positive/check.py`. It authenticates nine dependencies and checks the constants, exponent seams and beta telescopes. The [star checker](../../mean-ledger/fractional/star_bound.py), [473-box checker](../../mean-ledger/fractional/check.py) and [39-box checker](../../fractional-middle/check.py) reproduce the underlying fractional bounds.
