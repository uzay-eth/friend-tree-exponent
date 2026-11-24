> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../README.md#claim-26) for exact scope and citations.

# Fable 5(d,e): independent endpoint challenge

**HOLDS in scope**, at `88e13f36b39162f6ee60561a2b465a982b34af02`, `round2/newground/stochastic-domination/COUPLING.md`, §5. For (d) the certified monotone coupling in (c) is an input, not replayed here. The fixed-label tail exponent in (e) is a statement about the actual physical tree, not just the comparison urn. The lower bound and strictness include the delicate distinctions between `x=0`, `0<x<1`, and `x=1`.

## (e): the physical deficit, fractional urn, and matching lower event

At a fixed reachable size-`n0` tree, for any existing focal vertex `v` let

```text
r = Σ_(w~v) 1/D_w,        δ = D_v − r,       a = r,       b = n0 − r.
```

Both `a,b` are positive. The deficit is a sum of `1−1/D_w` over neighbours; it is nondecreasing since old degrees never fall and every new neighbour initially contributes zero. Thus at all later sizes the physical probability of gaining an edge is at most `(D_v(n)−δ)/n`. Couple this shifted degree below a two-colour fractional Pólya red count from `(a,b)` by the same uniform for the two gain indicators, drawing the rest of the tree update conditionally on gain or no gain. The order persists even when one red count is strictly larger. The standard urn Beta law holds for *positive real* initial weights: the exact probability of a specific red/white word is `(a)_R(b)_W/(a+b)_(R+W)`, the Beta mixture of independent Bernoulli trials. The fixed shift vanishes in embedded fraction, giving conditional stochastic upper domination by `Beta(a,b)`.

At the birth of a nonfounder `v≥3`, its sole neighbour is the attachment recipient, of **post**-birth degree `d≥2`. So `a=1/d≤1/2`, `b=v−1/d`. Using the weaker fixed shift `1/2` directly at birth gives Fable's unconditional upper comparison, and also that upper comparison conditional on *any* size-`v` birth history. For `0<ε≤1/2`, the beta upper tail is bounded by a constant times `ε^b` (the factor `t^(a−1)` is bounded on `[1/2,1]`). Thus there is no atom at one; there can be an atom at zero, consistent with the upper comparison. For `x<0` both tails are one; at `x=1` both are zero.

For the lower bound, require **every** attachment from size `n0` through size `m` to be at `v` (the final step is from `m−1`). Old neighbour degrees then stay exactly fixed, including the old parent, and every extra neighbour is a leaf. The *physical* path probability is exactly

```text
Π_(j=n0)^(m−1) (j−n0+a)/j = (a)_(m−n0)/(n0)_(m−n0) ≍ m^(−b).
```

The comparison urn is not being used as an order reversal for this lower bound. If the recipient at `v`'s birth is a leaf just before birth, then `d=2`; at size `j` on the forced path the degree of `v` is `j−v+1`, its parent contribution is `1/2`, and its new-leaf contribution is `j−v`. The chance of this birth recipient event is strictly positive: a finite star prefix exists and its centre can then send to a leaf. This achieves the least possible conditional exponent `v−1/2` and recovers the claimed unconditional one.

Here is a direct conditional continuation check. Write `L_n` for the number of leaf neighbours of `v`. Once `v` is a core, its positive change comes from a gain at rate `r_T(v)/n`; its negative change is a promotion at rate `L_n/[nD_v(n)]`. Hence `(L_n−1)/n` is a bounded submartingale since `r_T(v)−L_n/D_v(n)≥L_n−1`. Also `(D_v(n)−L_n)/n→0`: the number of core neighbours is at most its initial value plus the number of times `v` has been selected as initiator, which is `o(n)` a.s. (independent indicators with means `1/n`, or a dyadic Markov/Borel–Cantelli estimate). The submartingale limit is therefore `Z_v`; given the forced history through size `m`, `E[1−Z_v]≤(n0+1)/m` because `L_m≥m−n0`. Choose `m=ceil(2(n0+1)/ε)`; Markov on `1−Z_v≥ε` gives at least `1/2` for the **strict** tail `Z_v>1−ε`. The ceiling and the complement with `≥` handle the integer/threshold equality seams. Neither this continuation bound nor the upper bound conditions on future urn success.

This in fact proves conditional exponent **`n0−r` for every fixed physical restart tree**. In particular, the unconditional exponent should not be asserted for every birth history: at size four, attachment to an old degree-two centre produces post-birth degree three and exponent `11/3`, whereas attachment to an old leaf produces exponent `7/2`. For initial labels 1 and 2, apply the same restart calculation to the initial edge (`a=b=1`): both tails are of order `ε`, and their laws are equal by symmetry. Fable excludes them correctly from the `v−1/2` formula. Different exponents imply strict near-one ordering for any two distinct labels with the earlier one `≥2`, or for label 1 against any label `≥3`. They do not imply stochastic domination at any other thresholds or prove an asymptotic constant for the normalized endpoint tail. Fable's noncomputed near-one crossover is sufficient; it need not be uniform in the labels.

## (d): strictness from a residual cylinder, not just weak order

At the birth size of `w=u+1`, the `w`-root is always a leaf. The `u`-root is a core with positive marginal probability: at the preceding size it is a leaf, and its neighbour can choose it to receive `w`. This includes `u=2` at the original edge. Consequently the monotone transport in (c) has a positive set of paired seeds with A's root degree at least two and B's root a leaf. The matching has one matched child at the roots. The matched A-child's rate includes the A-root reciprocal, at most `1/2`; the corresponding B-parent receives `1` from its leaf root. Its other contributions dominate the paired A-child contributions by the odd-pair rate comparison. Therefore A's unmatched-child contribution gives strictly positive residual rate at the A-root, and the `1−1/D_A(root)≥1/2` difference gives strictly positive residual rate at the B-parent. The residual coupling can choose its two sites independently conditional on a residual step, without changing either marginal. Force these sites for finitely many consecutive steps: every new leaf is unmatched, the matching persists, and the same strict slack remains. This constructs a positive-probability **coupled** finite cylinder; it is not a simultaneous assertion about the two degrees in one realization of the labelled tree.

After `k` forced steps at common size `m`, A has at least `k` leaves, so `E[1−Z_A|history]≤1−(k−1)/m→0`. B's root remains a leaf, so `E[Z_B|history]≤1/m` from the ordinary bounded degree-fraction supermartingale. For each fixed `0<x<1`, Markov plus a union bound yields positive conditional probability of `Z_A>x≥Z_B` for large enough `k`, without presuming independence of the post-cylinder continuation. The residual construction gives physical marginals even conditional on the **joint finite** history because the one-step marginal probabilities are the physical Markov probabilities at each paired state. The process order forces the reverse event to have probability zero. This proves the tail inequality strictly.

At `x=0`, the low expectation of B alone is insufficient. Let `j` be the number of **other** leaves at B's parent while B's tagged root remains a leaf. For `j≥2` stop on tagged promotion or on `j=1`, assigning value one at either stop and use `1/j` before stopping. At parent degree `D≥j+1`, the `j+1` leaf clocks all add a leaf there; other arrivals only improve the bound. The parent promotes any other leaf at total rate `j/D`, and the tagged leaf at rate `1/D`. The continuous-time drift (or the embedded drift multiplied by size) is at most

```text
−1/j + 1/[D(j−1)] + (1−1/j)/D ≤ 0 .
```

Hence the probability of ever promoting the tagged leaf is at most `1/j` from the forced history. If it never promotes its degree is one forever and `Z_B=0`. Combine that bound with A's high conditional expectation by the same union bound. At `x=1` no strict inequality is possible, and none is claimed. The proof gives a separate positive gap for each `x`, not a uniform positive gap as `x↑1`. Outside the certified pairs the near-one result is independent of (c), while the all-threshold strictness argument needs the ordered coupling; it supplies no hypothetical transport certificate. It also does not turn an unconditional transport into one valid given an earlier tree.

## Exact replay

```text
taskset -c 8 nice -n 10 timeout 8s python3 -I -B round2/attacks/fable5/endpoint/check.py
PASS 40318 physical deficit transitions; 2802 leaf-stock drift rows; 536 tagged-leaf barrier rows
PASS 2555 fractional Pólya paths; all-red path matches actual all-attachment product
PASS 24 positive residual forcing steps with a fixed matching and protected B leaf
```

The companion enumerates ordinary labelled parent arrays and computes actual reciprocal send rates. It tests the endpoint-support and barrier identities independently of the rival code. The asymptotic exponent and strictness use the analytical arguments above; no empirical endpoint-tail extrapolation is used.
