> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT asymptotic audits and scoped Claude founder and reversal checks; earlier degree, hub, and spatial results. See the [release index](../../../../../../README.md#claim-18) for exact scope and citations.

# An explicit lower cylinder and a numerical late-label cutoff

For the original edge seed, put `p_v = P(Z_v>0)`. Using the [explicit spectral constants](../../../../../../../../reviews/cold/astra/checks/age-effective-sol.md), the construction below gives a bounded cylinder `L≤G` with **`πL>0.120996`**, and proves

```
N0 = 2^(2^396),     p_(n+1) > p_(n+2) for all n >= N0.
```

The certificate is a priori: it pays all untested cylinders using stationary age moments, and does not evaluate stationary weights or enumerate trees at label `n`. Its finite local transition formula and numerical parameters are specified below; the replay checks the analytic budgets in compact arithmetic and does not expand the enormous local kernel. The [exact prefix](../../early/further/README.md) reaches label 11. This leaves `11≤n+1≤N0` outside the two proved ranges; in particular the cutoff is not a sign proof for a specific such middle gap.

## 1. A uniform modulus for a finite next gain

The aim is to replace the unevaluated no-hit error in the [effective search](../README.md) by a numerical bound, including bare and live inputs. Consider a marked tag of initial degree one or two. Stop at its degree-`d` hit, and consider one inter-gain interval before that hit. Given the state at the start of the interval, the tag's degree is fixed and less than `d` until the next gain. Its next gain, if any, must come from one of its current neighbours. For one such neighbour `x`, work until the interval ends (by any source), or until `x` has at least `M=2^56` other leaves, excluding the tag. The latter event will only be used after `x` is core.

In each block of length `2d+128`, the conditional probability of one of these two stops is at least `1/32`, no matter what the exterior does:

1. In the first `d` units, if `x` is still an ordinary leaf, its rate-one donor clock would give the tag a gain on its very first ring while it is a leaf. Therefore either the interval ends or `x` is core at the end of these units, except with probability at most `e^(−d)≤1/8`. If already core this restriction costs nothing.
2. If `x` is core but has no other leaf, the tag's channel to `x` has rate at least `1/d` until the interval ends. This remains true if the tag itself is a leaf, when its choice is deterministic. A proposal at that rate either produces one other leaf at `x` or the prior stop has already occurred. Allow `d` units, for conditional success at least `1−e^(−1)>5/8`. If there is already another leaf, start the next stage immediately. A stage can finish early; only its maximum duration is charged.
3. Track one such other leaf and its descendants at `x`. Retain all its leaf births and use independent rate-`1/2` removal proposals for each tracked leaf, covering its actual removal rate `1/D_x≤1/2`. This embeds a linear birth–death process of rates `1,1/2` per member up to the stop. It can be run after stopping for the comparison. From one individual its first two moments at time `128` are `e^64` and `4e^128−3e^64`, so Paley–Zygmund gives probability at least `1/16` of reaching `M` during this stage (`e^64/2 > 2^56`, using `e>5/2`). There is no conditioning on the graph after stopping or on the tag remaining unselected through the block. The product of the three success floors exceeds `1/32`.

On hitting the other-stock cutoff before the gap ends, the [first-send barrier](../../early/README.md) bounds the probability of a later send from `x` to the tag, before any intervening tag gain, by `2/(M+1)`. It uses the fixed target slot and fits either a core target or an ordinary leaf. The comparison can be stopped on the first tag gain; before then the target remains in the specified status. Therefore the probability that the next gain is from `x` and occurs at a finite age greater than `t` is at most

```
(31/32)^floor(t/(2d+128)) + 2/(M+1).
```

This does not bound the probability of an *infinite* next gap, which may be large. Take `d=2^12`, `t=2^24`, and `T=d*t=2^36`. A hit after `T` but at finite age requires at least one completed gap longer than `t`; there are fewer than `d` gaps and fewer than `d` candidates for the incoming edge at each gap's start. A union bound at the successive stopping times gives an unconditional finite-late-hit error `δ≤d²[2^−60+2/(M+1)]<2^−28`, uniformly in the initial graph. This argument also covers an infinite locally finite exterior using the directed channels in the old-cell construction.

Let `P=P(hit degree d by time T)` on the full physical process. Then the ultimate hub probability is between `(1−1/d) P` and `P+δ`: a fresh leaf at the hit has survival probability at least `1−1/d`, and becoming a hub requires a finite hit. These bounds need neither global tree size nor a uniform conditional delay to hubhood itself.

## 2. Prescribed cylinder and finite arithmetic

Use old-core radius `R=2^40`. Admit an input ball only if:

* up through depth six, core degrees are at most `2^10` and total degrees at most `2^640`;
* throughout the radius-`R` ball, core degrees are at most `2^60` and total degrees at most `D=2^E`, with `E=8(60R+81)`.

Degrees on the boundary count all slots; the cylinder read radius is `k=R+1`. On each admitted input, for every prescribed serial and adjacent-parallel placement in `G`, insert both vertices and observe the tag's hit flag by time `T`. Run inside the original radius-`R` old-core ball, retaining exterior slots in the degrees and turning exterior arrows off. Cells involved in each tag's finite hit test are based within distance five of the original score root after the prescribed insertions. An old leaf promoted in an insertion adds at most one core at the placement site; the two insertions increase the local vertex or degree counts by at most a constant. Thus the sum `V` of the initial degrees of the queried cells is below `2^800` by the depth-six caps. By the [ordered-path bound](../../../../../../../astra/source/round2/operator.md), the error in the hit probability is at most `V T^(R−5)/(R−5)! ≤ 2^−30`. Indeed `j!≥(j/e)^j`, `j=R−5≥8T` and `e<4`. The error is unconditional, without filtering trajectories by whether the coupling succeeds.

An explicit rational endpoint formula is as follows. Count ordinary vertices in the truncated process, excluding inert exterior slots. Initially their number is at most `2^a`, where `a=E+60(R+1)+10` by the branching cap. Freeze at the hit (payoff one), or kill at size `B=2^b` (payoff zero), with `b=a+2T+32`. Before freezing/killing, total birth intensity is at most the current size; extending without absorption and using Yule domination shows the size-kill hit error is at most `2^a e^T/B<2^−32`. Let `A` be the rational generator of this killed, hit-flag chain read directly from the capped graph: sum `1/D_x` over the physical donor slots producing each successor, omit exterior arrows, and subtract the sum of the active rates on the diagonal. For the starting state `z`, use the rational number

```
P0 = [(Id + (T/U) A)^U 1_hit](z),
U = 2^u,  u = 32+2(b+36).
```

This expression also defines a finite sum by multiplying the one-step probabilities for flag histories, or it can be evaluated by repeated squaring. No values of any unknown occupation probabilities enter it. The transition matrix is nonnegative and stochastic (include the two absorbing outcomes). At step length `T/U`, comparison with the exact Poissonized kernel of rate `B` gives TV error at most `(BT/U)^2` in one step (compare a Bernoulli count of mean `BT/U` with the Poisson count of that mean). By contraction and telescoping, the error over `U` steps is at most `(BT)^2/U = 2^−32`. Thus `|P0−P|≤ε=2^−30+2·2^−32`.

Use leaf/middle hub endpoints `max(0,(1−1/d)(P0−ε))` and `min(1,P0+ε+δ)`. Their width is at most `1/d+δ+2ε`. Form each signed difference with the two correct endpoints and sum with the exact placement weights of [the adjacent score](../../adjacent/README.md), reversing endpoints for a negative coefficient. The sum of the absolute weights of differences is at most three per core (charge each old leaf to its parent). Define `L` to be the resulting lower endpoint clipped below at `−3`, and use `−3` on every excluded cylinder. This is a single specified finite rational cylinder, not a search for a favourable mask or a list of assumed-sign cylinders. The closed transition formula is a compact specification: only the error and cutoff arithmetic are executed in the replay. Explicitly tabulating all the admitted balls or evaluating all their endpoint values would be impractical at these parameters and is unnecessary for the sign certificate.

## 3. Weighted loss of every excluded cylinder

For a stationary uniform core, the [promotion-age moments](../README.md) give `πD^(1/4)≤2/(l−1/4)` and `π(5/4)^c≤(5/4)/(l−1/4)`, with `l=113/200`. Here `c` is core degree, `D` total degree. They come from a fresh core's degree-two Yule bound, its rate-one promotion clock and the stationary intensity bounded by one per cell member. Using `c²≤40(5/4)^c` and Cauchy–Schwarz on truncated tails gives

```
π[c 1_(D>J)] ≤ 32 J^(−1/8),
π[c 1_(c>K)] ≤ 32 (5/4)^(−K/2) ≤ 512·2^(−floor(K/8)).
```

Both constants follow with `l=113/200`; the second dyadic inequality follows from `(5/4)^4>2`. The same estimates bound root failure probabilities because the stationary one-ended core has `c≥1`. For a first capped-path failure at depth `i≥1`, reverse the path by mass transport. The first step has at most `c` choices at the bad endpoint and the remaining `i−1` previous vertices were capped. Account for each depth by a union bound. In particular, with uniform intermediate core-degree cap `K`, the exclusion mass through depth `r` is at most `(1+r)K^(r−1)` times the sum of the two displayed tails (for `r≥1`).

For the six-level prefilter the code bounds the total by `7·2^50 (512·2^−128+32·2^−80)<2^−22`. For the outer cap, `(1+R)(2^60)^(R−1)≤2^(41+60R)` and each of the two tail contributions after multiplication is below `2^−30`. These are bounds on the excluded measure under `π`, not assumption of independence of the marks and not discarded probability at deterministic size. Since both `G` and `L` are bounded by three, the excluded contribution to `π(G−L)` costs at most six times this measure. On the accepted balls the cost is at most six times the single-tag width. Therefore

```
π(G−L) ≤ 6 (prefilter + outer + 1/d + δ + 2ε) < 0.001467.
```

The stationary flux identity gives `πG≥β=97971/800000` independently of the cylinder, so `πL>0.120996>β/2`. This establishes the needed sign without evaluating either `π` or a single astronomical transition matrix. In particular no estimate of the apparent reversal in a sample is used.

## 4. Cutoff arithmetic and practical fixture

The audit's certified inputs are `η=2^−348`, `C=2^453`, `h_* =2^−60`, `H_*=2^68`. Its three exponent tests for `N0=2^s` are `s≥515·2^348`, `s≥20·2^348(982+2k)`, and `s≥25(71+ceil(log2 J))`, with `J=96(k+1)(D+1)^(k+1)+6`. The replay bounds `ceil(log2 J)` by `8+42+(k+1)(E+1)<2^92`, never constructing `D`, `J` or `N0`. All three tests pass at `s=2^396`; the signed de-Poissonization of the effective argument applies with `L` above.

As a size diagnostic, a much cheaper scalar test uses the [exact permanent-slot barriers](../../preplacement-transfer/README.md) on the three radius-one local profiles in `check.py`. Two adjacent old cores have `(D,L)=(2,0)` and `(2,1)` in each profile. The root has core degree two, with 0, 1, or 8 leaves. After combining both orientations on each of its old core–leaf edges, the outward lower floors are `−0.151601366`, `−0.097185100`, and `+0.000300006`. To make the first two local floors positive requires about `0.151601366` and `0.097185100` additional lower-score units, respectively. Negative lower endpoints are not negative actual scores. Even the positive eight-leaf fixture does not by itself certify the stationary sign; it supplies neither coverage nor a good bound on the other profiles. The full cylinder above pays that gap with a uniform endpoint modulus and certified moment tails.

## Replay

From `.` (CPU13, niceness 10):

```sh
timeout 20s taskset -c 13 nice -n 10 python3 -Bu round2/open/age-order/effective/cylinder/check.py
timeout 15s taskset -c 13 nice -n 10 python3 -Bu round2/checks/age-effective-sol.py
```

The first replay checks the rational error and exponent budgets, verifies the hit-flag kernel on small exact tests (including inert denominators and a hit coinciding with the size cutoff), and evaluates the three modest scalar fixtures. The second verifies the analytic spectral, normalization and de-Poissonization inputs numerically in rational arithmetic.

```text
weighted stationary error < 0.001467; required < 97971/1600000
pi L > 0.120996; certified N0 = 2^(2^396)
finite Euler kernel: hit flag, inert denominators and cutoff tie checked
bare/live: scalar lower -0.151601366; additional lower score to reach zero <= 0.151601366
one-leaf/live: scalar lower -0.097185100; additional lower score to reach zero <= 0.097185100
eight-leaf/live: scalar lower 0.000300006; additional lower score to reach zero <= 0.000000000
```
