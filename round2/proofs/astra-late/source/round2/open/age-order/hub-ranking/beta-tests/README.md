> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# Conditional binomial-tail order, and the limit at an interior threshold

Write `B_(k,r)(z)=P(Bin(k,z)≥r)` for `1≤r≤k`. Sampling attachment indicators here means independent Bernoulli samples **conditional on the limiting share** of the tag. It does not mean the consecutive physical attachments of the growing tree are iid. For the older and newer of the next two labels, the following new comparisons hold **conditional on every old tree** in the stated size range:

* **All old sizes:** strict order of `E B_(k,r)(Z)` when `r≥400(k−r+6)`. Examples with fewer required hits: one allowed miss (`k=r+1`) for every `r≥1650`; two allowed misses for every `r≥1800`.
* **Fixed hit count, broader samples:** for `r=3,…,8`, strict order for every `r≤k≤10r` at the respective old sizes `n≥97000,21000,23000,30000,39000,50000`. For example this covers the two tests `(k,r)=(9,3),(12,4)`, whose nominal hit fractions are both `1/3`. The size bounds are deliberately sufficient rather than optimized.

These are positive-kernel comparisons, including polynomials with **negative** monomial coefficients. They do not certify sharp order at `x=1/3` for the missing pairs. In fact the all-power inequalities by themselves, even with the quantitative gap and the simple individual attachment caps, cannot do so. Below is an exact obstruction at both that threshold and an explicit beta-smoothed test concentrated there.

## Common setup

The [two-birth identity](../../adjacent/README.md) applies to any one-variable continuation test. For an old tree of size `n≥2`, put `b=n+1`, `m=n+2`, and `A=Σ_u r_u/(D_u+1)`, using the old attachment intensities `r_u`. Each serial middle has its positive weight (`A` altogether); the serial tip debit is also weighted by `A`, and the absolute signed parallel coefficients sum to at most `A`. As in the [Mellin ledger](../mellin/README.md), the adverse parallel leaves with degree-two parents have coefficient sum at most `A/3`; for the rest the parent has degree at least three. No condition is being imposed *after* the two births while still claiming conditional order of the labels.

The binomial tail has a **single positive** kernel

```
K_(k,r)(x)=r binom(k,r) x^(r−1)(1−x)^(k−r),
E B_(k,r)(Z)=∫₀¹ K_(k,r)(x) P(Z>x) dx.
```

This is just the density of `Beta(r,k−r+1)` and follows by differentiating the binomial sum and telescoping. Conditioning and integrating the difference is legitimate because the two expectations are bounded. All estimates below apply to individual physical continuations before recombining the two-birth identity.

## 1. Few-miss tests on every old tree

Use the [physical initiator-colour comparison](../../global-transport/sharp-large/README.md#protect-the-family-for-a-constant-probability-prefix) on a serial middle. The original protection event on the first `2b+1` blue initiations and the conditional audited red-family loss give

```
P_middle(Z>x) ≥ C_b (1−x)^b,
C_b = (b−1)/(3b) · (1−log(3/2)) · (1−1/b)^b.
```

This input holds for `b≥3`: none of the colour construction or loss estimate imposes a lower bound on `b x`; the `256` cutoff in that proof is used later to compare the two tails pointwise. Each `b`-dependent factor increases with `b` (`log(1+t)<t` verifies the derivative of the last factor). The seven-term alternating log upper gives `C_b≥28/729=:C`. For any adverse leaf, the [locked urn cap](../../global-transport/sharp-large/README.md#all-adverse-arrows-including-zero-stock) yields `P_tip(Z>x)≤(1−x)^b/√(π b x)` from the common completion size. This bound need not be sharp or below one at small `x`.

Integrate both sides against the **same** beta kernel. Let `l=k−r+1` and `Q=∫K_(k,r)(x)(1−x)^b dx>0`. The serial lower is `CQ`. For the tip, retain the sharper cap's factor `sqrt(1−x)` and denominator `sqrt(π(b+1/2)x)` (equation (61) in the source). After tilting the kernel by `(1−x)^b`, the new variable is `X~Beta(r,l+b)`. Cauchy–Schwarz and `E[(1−X)/X]=(l+b)/(r−1)` give the signed numerator lower bound

```
A Q [ C − 2 sqrt((l+b)/((r−1)π (b+1/2))) ].
```

Use `π>25/8` and `(l+b)/b≤(l+3)/3` to get the exact sufficient inequality `4(8/25)(l+3)/(3(r−1))<C²`. The replay checks the stated corners and the bound `r−1≥400(l+3)` in the general family. Both the small-`x` adverse part and the large-`x` favourable part are included; a test concentrating high is not simply treated as a sharp high-threshold test.

## 2. Fixed hits using unsigned Mellin bounds

Let `j=k−r`, and let `U` have law `Beta(r,1)`. The positive derivative representation gives, for any `z∈[0,1]`,

```
B_(k,r)(z) / binom(k,r) = z^r E[(1−zU)^j].
```

The factor is in `[0,1]`, and its loss from one is at most `j r z/(r+1)` by Bernoulli's inequality and `E U=r/(r+1)`. Use this loss only for an **unsigned serial middle**. At the common fixed completion the continuous degree mass is `W=S Z`, with `S` independent of the embedded tree path and Gamma shape `m`. The bound `E W^(r+1)≤(2)_(r+1)=(r+2)!` is pure Yule domination from degree two, not a signed Mellin gap. Thus after dividing the beta-test reward by `binom(k,r) r!/(m)_r`, the serial lower is at least

```
M_r − (k−r) r(r+2)/(m+r).
```

Each tip upper is its genuine normalized `r`-moment upper, since the factor is at most one. Use the [directed physical endpoint values](../mellin/README.md) `M_r,U_(2,r),U_(3,r)` and pay the old-edge terms with `(4 U_(2,r)+2 U_(3,r))/3`. The common completion factors cancel in this proof, even though they could not be reused at random promotion times. In particular this is not an alternating subtraction of controlled and uncontrolled signed moments. The replay calls the frozen physical endpoint functions with exact `p=r`; for each fixed `r` it calculates the rational size needed to leave positive slack at `k=10r`, which covers the smaller `k` as well.

## 3. What transfers back to a sharp threshold

For any law on `[0,1]`, write `S(x)=P(Z>x)` and `M_(k,r)=E B_(k,r)(Z)`. For `r_-≤r_+` and `x∈(0,1)`, a two-test direct sandwich is

```
S(x) ≥ [M_(k,r_+) − B_(k,r_+)(x)] / [1−B_(k,r_+)(x)],
S(x) ≤ M_(k,r_-) / B_(k,r_-)(x).
```

The first follows by bounding the test by its value at `x` below the threshold and by one above; the second uses its value at `x` above. To certify a sharp comparison with finite tests, compare an actual lower bound for the older expectation at `r_+` to an actual upper bound for the newer at `r_-` through these formulas. Merely knowing two same-test signs is insufficient. For example `k=90`, `r_-=26`, `r_+=35` give errors `B_(k,r_+)(1/3)<1/5` and `1−B_(k,r_-)(1/3)<1/5`, but this is a filter calibration, not a sign for that pair of expectations.

For a limiting transfer at fixed old size, choose `r_k=ceil(kx+k^(2/3))` (admitted for all sufficiently large `k`). The binomial tail tends to zero at every `z≤x`, including `z=x`, and to one at every `z>x`, by the binomial variance bound. Bounded convergence now works for both laws without needing atomlessness. Nonnegative signed tests along that whole diagonal would imply weak sharp order at `x`; a uniform positive margin or a strict finite sandwich is needed for strictness. The proved families do not provide such a diagonal at `x=1/3`: few misses concentrate near 1 (the stated diagonal family only permits very high ratios), and the fixed-hit family requires taking old size large relative to the sample bound.

## 4. Rigorous obstruction to a moment-only deduction at `1/3`

The following are **abstract probability laws**, not claimed to arise from the physical placements. They obey the quantitative all-power floor for old size `n=20`, the quantitative hub gap there, and even the individual simple birth-size Pólya tail caps for labels 21 and 22. Set `q=1/800`, `w=2^−21`. The older candidate law `P` has a component of mass `q` with `Beta(1,21)` law and atoms of mass `w` at each of `1/4,1/2`. The newer candidate law `Q` has mass `2w` at `3/8`. Complete both to probabilities at zero.

For every real `p≥1` the two atoms of the older dominate the central atom in `p`-th moment by convexity. The positive continuous component alone gives `q B(p+1,21)·21` in the moment difference, equal to `420q Γ(p+1)Γ(20)/Γ(22+p)`; `420q=0.525` exceeds even the sharper constants in the all-power claim. The hub gap is `q`. The individual tail caps `P(Z>x)≤(1−x)^20`, `Q(Z>x)≤(1−x)^21` follow by checking the upper edges of each atom interval (exact rational checks), with the continuous component compared pointwise. Nonetheless at `1/3` the tail difference is `q(2/3)^21−w<0`. Moreover the explicit beta test `(k,r)=(600,201)` has negative gap: its point-mass contribution divided by `w` is less than `−9/10`, and its continuous contribution is exactly `q (k−r+1)_21/(k+1)_21`. The checker verifies that their sum is negative. The obstruction concerns deduction from those summaries, including an arbitrarily long list of additional real powers `p≥1`, not the truth of the missing physical comparison. The next ingredient for such a comparison must use additional physical placement information or a genuinely localized family of certified tests.

## Replay

```
python3 -IB round2/open/age-order/hub-ranking/beta-tests/check.py
```

The replay checks the exact constants, six finite-sample size cutoffs, the explicit counterexample, and 234 kernel/mixture identities at rational inputs. It imports the existing Mellin endpoint evaluator read-only and does not change any other proof lane.
