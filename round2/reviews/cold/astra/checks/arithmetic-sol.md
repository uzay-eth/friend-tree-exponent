# D audit — original trigger arithmetic

**Verdict: CORRECT in scope.** The primitive-numerator obstruction holds for the positive chronological prefixes of the original cap-20 trigger, including prefixes with an additional size cap at least 20. The three exact exclusions hold conditional on the archived physical interval. No blocking gap found. Neither result establishes irrationality, a clearing-denominator bound for evaluated fractions, or a pole of the infinite trigger. This audit uses the physical first-moment identity and bracket as inputs; it does not reprove the bracket or review Long's manuscripts.

## Rows, poles and the visit count

Use physical-trigger §4, equation (4) (`finite/gauge/physical-trigger.md#4-positive-chronological-series-and-uniform-error`). A live fringe of size $`n`$ has intrinsic jump rate $`n-1/D`$: the only unused vertex-clock choice is the root sending to its inert parent slot. The holding denominator is $`s+\alpha`$, $`\alpha=1/i+n-1/D`$, where $`i`$ is the permanent root label. Killing adds $`1/i`$. Old-child copies and untriggered direct newborns belong in transition numerators, not the holding denominator: they do not end the original particle's residence. A chronological path follows the new copy or newborn when it takes such an edge. Its terminal virgin earns $`1/2`$ in the numerator at **any** mark. A path of length $`\ell`$ has $`\ell+1`$ holding factors, not $`\ell`$. Every row has finitely many successors and intrinsic size can rise by only one per transition, so each prefix is finite and strictly proper.

On a virgin, $`\alpha=3/2+1/i`$; off the virgin set $`n\ge3`$ and $`D\ge2`$, so $`\alpha\ge51/20`$. Only $`J_{20}`$ attains $`a=31/20`$. Let $`r`$ be the largest multiplicity of $`s+a`$ among positive terminal paths. For $`r>0`$, multiplying the prefix by $`(s+a)^r`$ and taking $`s\to-a`$ kills every lower-order term and leaves a strictly positive sum: all other factors are positive **at that point**. For $`1-T_m`$ the coefficient has opposite sign, but the same order. This is a statement about rational continuations of finite sums. At another holding value $`b>a`$ the remaining factors can have either sign; their orders need not survive. As a positive-history control, $`1/[(s+2)(s+3)]+1/(s+3)=1/(s+2)`$ cancels the nonminimal pole. Any remaining poles still come from the original linear holding factors with multiplicity at most their maximal path multiplicity.

In fact $`r=\lfloor m/19\rfloor`$. The lower bound concatenates 18 star gains and the degree-20 direct birth, starting at $`J_2`$ and then at $`J_{20}`$ for every repeat. The direct birth has rate $`19/20`$ and the terminal virgin has rate $`1/2`$. For the reverse bound, restart at the initial virgin or at each visit to $`J_{20}`$. There is no internal mark-20 vertex at a restart, and the current root mark does not transfer to a child. The first subsequent acquisition of a mark-20 **child** (by intrinsic promotion) or a mark-20 **newborn** (by direct birth) requires a promoter of degree at least 20. If an old-copy edge is to bring the lineage to a virgin of mark 20, the relevant interior mark must have been acquired previously in this segment. The segment starts at size two, every transition increases size by at most one, and every fringe root or interior degree is at most the fringe size, even with the root slot. At least 18 transitions are needed before any promoter can have degree 20, and at least one more to act. Thus visits, including the first, cost at least 19 each. The explicit witness uses size at most 20, so the same order holds with the additional indicated size cap. This sharpening is not needed by the growth proof.

## Primitive normalization and omitted mass

Write $`1-T_m=A_m/B_m`$ in lowest terms over $`\mathbb Q[s]`$ with $`B_m`$ monic. Since $`T_m(\infty)=0`$, $`A_m`$ is also monic of the same degree. Clear the coefficient denominators of $`A_m`$ with their least common multiple $`c_m`$. If $`c_m A_m`$ had integer content $`g>1`$, then $`g\mid c_m`$ by its leading coefficient, and $`c_m/g`$ would already clear every coefficient. Thus $`p_m=c_m A_m`$ is primitive with leading coefficient $`c_m\ge1`$. Clearing a common denominator for the rational function and then extracting the numerator's content gives exactly the same polynomial. Every factor of the reduced **monic** $`B_m`$ is $`s+b`$ with $`b\ge31/20`$; at $`s=\mu>0`$ all are greater than one. None of the potentially cancelled nonminimal factors is needed for $`B_m(\mu)\ge(\mu+31/20)^r`$.

For $`N=\max(20, m+2)`$ the star path keeps label 2 for $`N-2`$ gains, resetting to $`J_{20}`$ on the next edge. Its $`N-1>m`$ transitions exclude it. Each gain weight at $`s=\mu<1`$ is at least $`(d-1)/(d+2)`$, whose product for $`d=2, \dots, N-1`$ is $`6/[(N-1)N(N+1)]`$. The direct weight is at least $`(N-1)/[N(N+2)]`$. The terminal holding payment is $`1/[2(\mu+31/20)]>1/6`$. By positivity and $`T(\mu)=1`$ this gives

```math
1-T_m(\mu)\ge [N^2(N+1)(N+2)]^{-1}, \qquad
p_m(\mu)\ge \frac{(\mu+31/20)^{\lfloor m/19\rfloor}}{N^2(N+1)(N+2)}.
```

An additional size cap only removes paths, so it cannot recover the omitted mass. The liminf exponential rate in the source follows without assuming $`\mu`$ rational or algebraic, and without a height estimate. Summing unrestricted lengths after killing by size is a different operation: it inverts a matrix with off-diagonal transitions, so its poles are governed by that matrix rather than by maximum visit counts. For the raw size-killed height bound, $`H=\mathrm{lcm}(1, \dots, \max(N,20))`$ clears each rate (a sum of reciprocals or repeated reciprocal terms). The row coefficient norm after scaling is at most $`H[1+(N+1/2)+N+N/2+1]`$, bounded by $`H(4N+4)`$. The determinant expansion is bounded in coefficient $`\ell^1`$ by the product of these row norms. The rank-one determinant identity yields $`1-T_N`$ as the augmented/unaugmented quotient wherever the latter matrix is invertible. This is a raw determinant upper bound, not a lower bound on a reduced polynomial or on the denominator at a trial value.

## Rational evaluation is a different reduction

The lower bound just proved is for the **polynomial value** $`p_m(\mu)`$ at the unknown critical point. At a hypothetical $`\mu=p/q`$ in lowest terms, if $`d_m=\deg p_m`$ then $`q^{d_m}p_m(p/q)`$ is a positive integer; the elementary lower bound $`q^{-d_m}`$ supplies no contradiction to exponential **growth**. One may instead rule out a fixed $`p/q`$ by an enclosure of the **unscaled** $`T(p/q)-1`$ excluding zero. A positive prefix smaller than one at a candidate does not do that on its own; at the true point every finite prefix is strictly below one. The positive analytic crossing is compatible with a rational root in the source's countermodel: $`G(r)`$ is the geometric sum one and the rank-one matrix gives $`G/(19-18G)`$. The locally uniform meromorphic sum away from the negative integers has nonzero residue $`2^{-j}(r+j)`$ at each $`-j`$; infinitely many poles imply nonalgebraicity over $`\mathbb C(s)`$ (a polynomial equation can have finite-plane poles only at zeros of its leading coefficient). The invertible Möbius response preserves that property, not special-value irrationality.

At a known rational input the scalar fraction in lowest terms need not be obtained simply by homogenizing primitive polynomials. For instance take positive histories
```math
T_M(s)=\frac{1/M}{s+2}+\frac{9/2-3/M}{(s+2)^2}, \qquad M\ge2\text{ even}.
```
The reduced monic denominator is $`(s+2)^2`$ and the primitive numerator of $`1-T_M`$ is $`M(s+2)^2-(s+2)-(9M/2-3)`$, primitive since its leading and linear coefficients are $`M`$ and $`4M-1`$. At $`s=1`$ the ratio is $`(9M/2)/(9M)=1/2`$ regardless of $`M`$. This is an adversarial algebraic control, not a physical trigger or a claim about its evaluated gcd. For any one coprime integer polynomial pair a fixed nonzero resultant bounds common divisors of the homogenized values of coprime inputs (use the dehomogenized Bézout identity on each of the two charts modulo each prime power). Across these prefixes that bound itself varies; no uniform scalar cancellation bound follows from primitivity or positivity. The original obstruction is correctly limited to using the stated single primitive polynomials as small values; it does not rule out alternative arithmetic constructions or evaluated-fraction methods.

## Exact exclusions and replay

- **Denominator 76:** the exact closed physical interval lies strictly between Farey neighbours $`13/23`$ and $`30/53`$ (determinant one) and contains their mediant $`43/76`$. For any reduced $`p/q`$ between them, $`q=53(23p-13q)+23(30q-53p)\ge76`$. Presence of the mediant in the interval is not equality to $`\mu`$.
- **Quadratic height 7:** for positive leading coefficient $`a`$ enumerate $`1\le a\le7`$, $`|b|,|c|\le7`$, keep primitive triples of positive nonsquare discriminant and test for a zero on either monotone segment split at $`-b/(2a)`$. The unique survivor is $`(3,-7,3)`$. This covers irreducible degree two and uses exact endpoint signs, not floating root estimates.
- **Square-discriminant denominator 57,200:** for $`s>0`$, $`s^2+4`$ is a rational square exactly when $`t=(s+\sqrt{s^2+4})/2=u/v>1`$ is rational, with $`s=(u^2-v^2)/(uv)`$. For coprime $`u, v`$, no prime of $`uv`$ divides $`u^2-v^2`$. Independently of the author's $`u, v`$ search, enumerate reduced $`p/q`$ in the interval for $`q\le57200`$ and integer-test whether $`p^2+4q^2`$ is square. Only $`32361/57200`$ survives, from $`u/v=275/208`$. The wider-interval example $`2603/4602`$ has square discriminant but lies below the exact lower endpoint. The unrestricted specialization class is dense by continuity of $`t-1/t`$ on $`t>1`$.

Run from the repository root:

```sh
nice -n 15 taskset -c 8 timeout 8s python3 -B round2/checks/arithmetic_sol.py
nice -n 15 taskset -c 8 timeout 8s python3 -B round2/open/arithmetic/check.py
```

Both commands exited zero. The cold checker reports:

```text
PASS: exact archived endpoints; Farey denominator 76, quadratic height 7, direct discriminant denominator 57200
PASS: capped-generator reset edges, positive witness Laurent coefficients, 7 exact omitted star paths
PASS: 5 physical prefix reductions; nonminimal-pole and scalar-gcd adversarial examples
```

The author checker independently reported `43/76`, `32361/57200`, and `(7, 3, -7, 3)`, and passed its path and rational-crossing controls. The cold checker uses the same capped intrinsic generator for the reset edges, but builds and cancels its own small-depth rational functions. Those five physical prefixes do not reach $`J_{20}`$; the general pole order and growth bound are analytic. Neither replay reconstructs the physical enclosure or tests the value of the infinite trigger at an unknown argument.
