> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude consistency checks; earlier spectral constants and, for the named early-label splice, Fable transports. See the [release index](../../../../../README.md#claim-24) for exact scope and citations.

# Eventual strict stochastic order at all thresholds

Let $`Z_v`$ denote the limiting degree fraction of birth label $`v`$ in the original edge-seeded random friend tree. Write $`p_v(x)=\Pr(Z_v>x)`$ and $`d_v(x)=p_v(x)-p_{v+1}(x)`$.

**There is a finite $`v_0`$ such that every $`v\ge v_0`$ satisfies $`d_v(x)>0`$ simultaneously for every $`0\le x<1`$.** At zero this compares hub probabilities. At one both tails vanish. The proof gives an existential cutoff. The exact early certificates and this eventual result leave finitely many adjacent pairs unresolved. Fable Claim 5 proved all-threshold domination through pair $`(19,20)`$ using two programs, with a further single-program receipt for $`(20,21)`$. Its near-one order for each fixed birth label supplied the original comparison urn. [Fable 5 audit](../../../attacks/fable5/README.md)

## The three threshold regimes

Let $`A(y)`$ and $`B(y)`$ be the positive continuous scaled tail and scaled signed gap of Claim 18. For bounded $`y=vx`$, uniformly on each fixed compact interval including zero,
```math
 v^{1-\mu}p_v(y/v)\longrightarrow A(y),\qquad
 v^{2-\mu}d_v(y/v)\longrightarrow B(y),\qquad B(y)\ge(1-\mu)A(y)>0.
```
The original statement uses $`n=v-1`$. Uniform convergence and continuity on compact intervals allow the displayed normalization and argument $`y/v`$. [Scaled tails and signed differences](../scaled/README.md)

For every fixed $`0<a<1`$, uniformly on $`a\le x<1`$,
```math
 p_v(x)\sim \kappa v^{\mu-3/2}x^{-1/2}(1-x)^{v-1/2},\qquad
 \frac{p_{v+1}(x)}{p_v(x)}=(1-x)(1+o(1)).
```
Here $`\kappa=(\mathbb E X)/\sqrt\pi`$ times the stationary mean of the physical rate of a birth onto a leaf multiplied by its limiting path likelihood. The proof keeps the conditional outside-arrow law and its clock. In the comparison, a Pólya endpoint fixes the Bernoulli parameter for attachments to the tag. Counting outside jumps gives uniform estimates even as this parameter tends to one. [Interior equivalent](../fixed-threshold/uniform/README.md) · [Upper-edge proof](../fixed-threshold/uniform/upper-edge/README.md) · [Independent upper-edge audit](../../../../../../../reviews/cold/astra/checks/fixed-upper-edge-sol.md)

For every sequence with $`x\to0`$ and $`vx\to\infty`$, the same tail equivalent holds and $`d_v(x)\sim x p_v(x)`$. This signed estimate uses both birth orders. The middle vertex of a new two-edge chain, grown from an old leaf, sees the same limiting exterior likelihood as the one-birth tag. Its beta density is exactly $`2v\theta`$ times the single-birth density at endpoint $`\theta`$. Its placement rate supplies the factor $`1/2`$. The other placements are negligible compared with the gap itself. [Transition proof](../fixed-threshold/transition/README.md) · [Exact two-birth audit](../../../../../../../reviews/cold/astra/checks/fixed-transition-sol.md) · [Second hostile review](../../../../../../../reviews/cold/astra/checks/fixed-transition-second.md)

If there were failures at arbitrarily late labels, take a sequence of thresholds realizing them and pass to a subsequence on which $`vx`$ is bounded or tends to infinity. The bounded case contradicts compact-scale positivity. In the unbounded case pass to a further subsequence with $`x`$ convergent in $`[0,1]`$. A zero limit uses the signed transition estimate. A positive limit uses the upper-edge ratio, uniformly including limits at one. This proves the common cutoff. The same argument gives the first-order gap $`d_v(x)\sim x p_v(x)`$ along **every** sequence $`vx\to\infty`$. Strict ordering also holds for expectations of all positive powers simultaneously beyond the same cutoff, by integration of the tails.

## The tail profiles and one uniform formula

The overlap of the two scales determines the large-$`y`$ profiles, including their constants,
```math
 A(y)\sim\kappa y^{-1/2}e^{-y},\qquad
 B(y)\sim\kappa y^{1/2}e^{-y},\qquad
 -A'(y)\sim\kappa y^{-1/2}e^{-y}.
```
For any sequence $`y_j\to\infty`$, choose labels $`v_j`$ large enough to approximate each positive fixed value $`A(y_j),B(y_j)`$ to relative error $`1/j`$ with the compact-scale limits. Require also $`y_j^2/v_j\to0`$. Then $`x_j=y_j/v_j\to0`$, $`v_jx_j\to\infty`$, and $`(1-y_j/v_j)^{v_j-1/2}\sim e^{-y_j}`$. Substitution in the two transition equivalents gives the first two formulas for every sequence $`y_j`$. The derivative follows from the separately proved identity $`B=(1-\mu)A-yA'`$. [Profile identity and hub mass intensity](../spectrum/README.md)

In particular the following two relative equivalents hold **uniformly on the whole interval $`0\le x<1`$**. With $`y=vx`$ and $`E_v(x)=e^{vx}(1-x)^{v-1/2}`$,
```math
 p_v(x)\sim v^{\mu-1}A(y)E_v(x),\qquad
 d_v(x)\sim v^{\mu-2}B(y)E_v(x).
```
On bounded $`y`$, $`E_v\to1`$ uniformly and both profiles have a strictly positive minimum. On sequences with $`y\to\infty`$, insert the large-$`y`$ equivalents into the right sides. These expressions agree with the one-label and signed-gap estimates above, including sequences $`x\to1`$. Every sequence of thresholds has a subsequence covered by one of these cases. Very small probabilities near one enter through the exact factor $`(1-x)^{v-1/2}`$.

## Consequences for moments

For fixed $`\alpha>0`$, the known moment coefficient is $`b_\alpha=\alpha\int_0^\infty y^{\alpha-1}A(y)\,dy`$. Thus $`b_\alpha\sim\alpha\kappa\Gamma(\alpha-1/2)`$ as $`\alpha\to\infty`$. Indeed the integral over a bounded interval divided by $`\Gamma(\alpha-1/2)`$ tends to zero, and the relative error in the integrand on the rest can be made arbitrarily small. The stationary positive birth-mass intensity $`\nu`$ has tail $`A/\mathbb E X`$. Its density therefore has equivalent $`\kappa y^{-1/2}e^{-y}/\mathbb E X`$. Its Laplace integral satisfies $`\int e^{t y}\nu(dy)\sim(\kappa\sqrt{\pi}/\mathbb E X)(1-t)^{-1/2}`$ as $`t\uparrow1`$, by integrating the positive density between two relative bounds above a fixed cutoff. [Fixed-power moment laws](../moments/README.md)

For **every** pair of sequences $`v\to\infty`$, $`\alpha=\alpha_v\to\infty`$, the two equivalents take a beta-integral form,
```math
 \mathbb E[Z_v^\alpha]\sim\alpha\kappa v^{\mu-3/2}
       \mathrm B(\alpha-1/2,v+1/2),\qquad
 \mathbb E[Z_v^\alpha-Z_{v+1}^\alpha]\sim\alpha\kappa v^{\mu-3/2}
       \mathrm B(\alpha+1/2,v+1/2).
```
In the second formula both labels use the **same power** $`\alpha_v`$. Dividing the respective beta integrals gives a relative moment gap asymptotic to $`(\alpha-1/2)/(\alpha+v)`$.

To check the passage to moving powers, apply layer cake to the two uniform formulas. Divide the first integrand by its large-$`vx`$ form, whose normalized weight is the beta law $`T\sim\mathrm{Beta}(\alpha-1/2,v+1/2)`$. The remaining profile ratio $`A(y)e^y\sqrt y/\kappa`$ is bounded on $`[0,\infty)`$ and tends to one. For the gap use $`T\sim\mathrm{Beta}(\alpha+1/2,v+1/2)`$. The remaining ratio $`B(y)e^y/(\kappa\sqrt y)`$ tends to one and is bounded by $`C(1+y^{-1/2})`$. In either case, for beta parameters $`a,b`$, $`\mathbb E[(vT)^{-1}]=(a+b-1)/(v(a-1))\to0`$. Hence $`vT\to\infty`$ in probability and $`\mathbb E[(vT)^{-1/2}]\to0`$. These bounds justify both relative integrations without a prescribed relation between the two diverging indices.

The compact-interior proof is pinned at `389eeff`, its upper-edge and small-threshold extensions at `4492e16`, and the earlier compact birth-scale theorem at `25a2a2f`. The reviewed fixed-label near-one equivalent at `82864bc` uses a different order of limits. The numerical cutoff $`2^{2^{396}}`$ in the earlier hub-probability paper applies to the threshold $`x=0`$. Computing a cutoff that works for every $`x`$ remains a separate problem.
