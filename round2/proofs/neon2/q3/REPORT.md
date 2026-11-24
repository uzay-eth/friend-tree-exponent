# Addario-Berry et al., Open Question 3: exact adjacent comparisons

**Pinned claim commit:** `6cc2c8873885770d2868fad543596760d2b54e55`.  
**Review and dependencies:** Gemini 4 Neon proof/code audit within Neon-2; independent Claude exact reimplementation and same-family GPT replays. See [Claim 2 dependencies](../README.md#inherited-proof-inputs).

[Release replay and dependencies](../README.md) · [Astra review](../../../reviews/attacks/astra-on-neon2-2.md) and [Fable review](../../../reviews/attacks/fable-on-neon2-2.md).

**Author:** GPT-6.2 Sol Ultrafast; proof and code audited by Gemini 4 Neon.

Write `m_u = E[Z_u]` and `h_u = P(Z_u > 0)`. For the original edge seed, the exact finite-state certificates prove

```
m_1 = m_2 > m_3 > m_4 > m_5 > m_6 > m_7 > m_8,
h_1 = h_2 > h_3 > h_4 > h_5 > h_6 > h_7 > h_8.
```

They use tree sizes at most 16. The full stochastic-order question remains open here, including for these adjacent pairs. The equalities follow from exchangeability of the seed labels.

## Signed rooted-tree reductions

Let `K(R)` be the law of the limiting root mass from the rooted, unlabelled tree `R`. Each birth attaches at `x` with probability `a_T(x)/n`, where `a_T(x)=sum_{y~x} 1/D_y`. Taking conditional future laws at the birth of the younger label gives a signed measure on rooted trees, and equal rooted positions cancel. Rooted isomorphism is detected by the balanced-parenthesis string around the lexicographically sorted child strings.

- At size 3, label 3 is an end and label 2 is the center with probability `1/2`. Thus `Law(Z_2)−Law(Z_3)=(K(T3c)−K(T3l))/2`.
- At the next birth, the size-3 center receives the new leaf with probability `2/3`; labels 3 and 4 are then both star leaves. The two old ends each receive it with probability `1/6`. Only attachment at label 3 gives different orbits: 3 is interior and 4 is an end. Thus `Law(Z_3)−Law(Z_4)=(K(P4i)−K(P4e))/6`.
- At size 5, the exact reduction has five nonzero terms. At size 6, it has nine; neither reduces to a single orbit contrast by rooted isomorphism alone.

Here are the last two reductions, in units of `1/72` and `1/120`, respectively:

| Size | Root string | Coefficient units |
|---:|:---|---:|
|5|`((((()))))`|−3|
|5|`(((()())))`|−3|
|5|`(((())()))`|−1|
|5|`(((()))())`|3|
|5|`((()())())`|4|
|6|`(((((())))))`|−1|
|6|`((((()()))))`|−2|
|6|`((((())())))`|−2|
|6|`((((())))())`|1|
|6|`(((()()())))`|−2|
|6|`(((()()))())`|2|
|6|`(((())()()))`|−1|
|6|`(((())())())`|2|
|6|`((()()())())`|3|

Analytic computation of the table: start at size `m=3` with the unique path, or at size `m=4` with the star/path mixture of weights `2/3,1/3`. For each first target `x` append leaf `b` to get `T'`; for each second target `y` in `T'` append leaf `c` to get `T''`. Add

```
(a_T(x) a_T'(y) / (m(m+1))) [K(T'',b) − K(T'',c)].
```

By grouping equal targets under the starting tree symmetries, the size-3 path contributes the size-5 row directly. Before applying the size-4 mixture, the star contributes the size-6 strings in table positions `2,5,6,7,9` with coefficients `(-2,-2,2,-1,3)/80`; the path contributes positions `1,3,4,8` with coefficients `(-1,-2,1,2)/40`. This finite two-birth identity is exact for limiting laws; propagation to size 9 is an independent check, not a limiting argument.

## Mean and hub bounds

**Mean.** At size `n≥3` let the root have degree `d` and `k` adjacent leaves. Its degree gain probability `a/n` is at most `d/n`. Every hit on the root adds one leaf; a promotion of one of its leaves is a separate event of probability `k/(nd)≤1/n`. Hence `E[L_(n+1)−1 | T_n]≥(1+1/n)(k−1)`. The positive part of `(L_n−1)/n` is also a submartingale by conditional Jensen. Its bounded limit is at most `Z`, since `L_n≤D_n`; this gives bounds `max(k−1,0)/n ≤ E[Z|T_n] ≤ d/n`. In particular the size-3 center has mean at least `1/3`, while the size-3 end already has mean at most `(1+1/6)/4=7/24` at size 4, yielding the standalone analytic gap `m_2−m_3≥1/48`.

**Hub lower.** Tag the root's `k` current leaf neighbours and all leaf offspring produced by their own rings. If there are `j` living tags, births have rate `j`, deaths at most `j/max(d,j+c)` with `c=d−k` (initial untagged neighbours are permanent). Couple below the tags a birth-death chain with these extremal rates and absorption at zero. Its survival implies the physical root becomes a hub: the chain reaches arbitrarily large tag counts on survival (it cannot stay forever in a finite positive band); from `j` physical tags the Yule-minus-loss bound of the age-order README makes nonhub probability tend to zero as `j` grows. For the chain, set `p_0=1` and `p_i=prod_{j=1}^i 1/max(d,j+c)`. The probability of hitting level `b` before zero is `sum_{i<k} p_i / sum_{i<b} p_i`, from the first-step recurrence for adjacent increments of the hitting function. Let `b→∞`. Bound the denominator above by `sum_{i=0}^{28} p_i + p_28/(max(d,29+c)−1)` using a geometric tail. Also take the maximum with the elementary fixed-degree bound `1−d^(−k)` for `d≥2`. At `k=0` use zero.

**Hub upper.** For each neighbour `x` charge the probability it ever sends a birth into the root. Before the first such arrow, tag all its `r` other leaf neighbours and their leaf-born offspring. With `j` surviving tags there are `c=D_n(x)−r−1` permanent other untagged neighbours. Thus tag births have rate `j`, tag deaths occur at rate at most `j/(j+1+c)`, and the charged arrow has rate at most `1/(j+1+c)`. A decreasing nonnegative test function `f` with success value 1 is therefore an upper envelope if its drift under these maximal rates is nonpositive. The simple `f(j)=min(1,2/(j+1))` has this property: at `j=1` only the birth term contributes, and for `j≥2` the extremal drift is `−2j/((j+1)(j+2))+1/(j+1+c)≤0`. Stopping at the arrow and using boundedness gives the stated probability bound. To sharpen, solve the exact first-step system on `1≤j<30`, with boundary `f(0)=1`, `f(30)=2/31`, using integer-scaled rates `j(j+1+c), j, 1`. The code checks that this solution is decreasing and satisfies every equation; the boundary values bound future success by the simple envelope. Thus the finite solution also bounds success, by stopping at the arrow or a boundary. A hub must gain again and its first future gain comes from an existing neighbour; sum these arrow charges and clip at 1.

**Exact sign certificates.** At a given size, propagate signed orbit weights `q_R` under the same rooted transition kernel. If `g(R)` is a future mean or hub probability with bounds `[l_R,u_R]`, use the lower endpoint `sum_{q_R≥0}q_R l_R + sum_{q_R<0}q_R u_R`, and swap endpoints for the upper. This uses all equal-orbit cancellations; simply subtracting independent full-label intervals can fail even when the signed bound certifies. The lower certificates below are exact; the hub column is conservatively floored to a millionth using rational arithmetic.

| Pair | Size | Certified mean-gap lower bound | Certified hub-gap lower bound |
|---|---:|---:|---:|
|2,3|10|`572536956559579/3950456832000000`|`221023/1000000`|
|3,4|10|`16129411373/1881169920000`|`22803/1000000`|
|4,5|10|`32213487647/33861058560000`|`6177/1000000`|
|5,6|12|`8105547789593251/21901332676608000000`|`3955/1000000`|
|6,7|14|`113768087106544029413632801/664465320523586469888000000000`|`2655/1000000`|
|7,8|16|`685232928195421574000390165654450453/8169185426332976091577398067200000000000`|`1826/1000000`|

The size-12 individual mean bounds requested are:

| Label | Lower | Upper |
|---|---:|---:|
|1,2|`452018868722384593765169/1931697542076825600000000`|`1528596577808524482571/3942239881789440000000`|
|3|`11154305836385049599861/241462192759603200000000`|`674821379387643087407/3942239881789440000000`|
|4|`20769131259594331473719/965848771038412800000000`|`521735854646219931077/3942239881789440000000`|
|5|`1841579750375855635699/160974795173068800000000`|`25077924621779917369/219013326766080000000`|

The individual intervals already separate 2 from 3, but overlap for 3 and 4. For that pair, at size 10 the `P4` interior lower bound is `0.165136851…`, versus end upper bound `0.162545644…`; signed propagation adds further cancellation. Full exact intervals at sizes 6, 10, 12, 14 and 16 are replayable; pass `--full-hub-fractions` to print the unreduced-to-millionths rational hub endpoints.

## Monte Carlo (diagnostic)

`mc_t13.csv` contains all 25 rows, seven thresholds, and paired negative differences for all `u<v≤25`. The simulator samples discrete births and exponential waiting times until the first 25 labels exist, then schedules independent rate-one clocks only for cores. The leaf stock at each core grows as a Yule process between interventions; lazy sampling uses geometric/Yule increments (negative binomial for a stock). Every core clock selects uniformly among its core neighbours and marked or unmarked leaves. A hit to a leaf turns it into a new core with one attached leaf; a hit to a core increments its leaf stock. This accounts for the untracked population exactly in distribution. At time 13 the observed mean size was `885894.432 ± 1402.170` SE, versus the exact Yule mean `2e^13 ≈ 884826.784`.

At 200,000 independent runs, no mean or `D≥20` proxy point estimates reversed for any pair. `D≥20` estimates for labels 2, 3, 4, and 25 were `0.699115`, `0.428390`, `0.354740`, and `0.142970`. The `D≥50`/`D≥100` columns are threshold sensitivity diagnostics. None estimates the limiting hub probability with a controlled finite-time bias. The population mean column uses `D/N`; tail columns use `1[D/N>x]`. SEs are the sample-SD-over-root-n values (equivalent to delete-one jackknife for these sample means), using within-run contrasts for pairwise checks. The largest negative tail contrast was at `19,20`, `x=.30`: `−0.000030 ± 0.00001225` from just six discordant events, two-sided exact sign-test `p=1/32` before accounting for the full scan. No compelling reversal.

## Replay and checks

From `round2/proofs/neon2/q3`:

```sh
python3 -B symmetries.py
python3 -B cert_q3.py --size 16 > replay_n16.txt
cmp cert_n16.txt replay_n16.txt
g++ -O3 -std=c++17 -Wall -Wextra -pedantic sim_q3.cpp -o sim_q3
./sim_q3 200000 13 330283 > replay_mc_t13.csv
```

Symmetry enumeration cross-checks direct labelled arrow tallies against rooted propagation through `n=9`; sign swaps fail for each nonzero reduction. The certificate independently checks raw labelled-history and rooted sums for all four observables at `n=8`; local one-step inequalities are audited through `n=7`, with planted invalid `L/n` and reverse `D/n` claims rejected (1,342 and 5,900 counterexamples). The certificate also checks that reversing a signed contrast negates and reverses its bounds. Wall times: symmetry 3.74 s, size-16 certificate 140.31 s, production simulation 205.31 s.
