# Cold audit: tiny powers and the founder pair

**Holds** for the small-powers proof (`open/age-order/hub-ranking/small-powers/README.md`) blob `c2ee87a3a4c3ec15595063abde1ec810261f0b9a` and checker blob `40659d732606294bb3a10608cbfd69bfc0201c5d`, read in ee0c's directory. The all-tree result audited here is strict conditional order at every old size `n≥2` and real `0<p≤10⁻⁹`; the two-core-or-more raw Mellin numerator retains `371/12500`. The independent founder result is unconditional strict order of labels 2 and 3 throughout `0<p≤1`; labels 1 and 2 have the same law by seed symmetry. The near-one perturbation also ran in the producer replay but is outside this cold verdict.

## Tiny powers: pay only for paths the frozen proof purchases

At common size `m=n+2`, all exponential holding times have size-only rates and are independent of the embedded tree chain. Thus the moment conversion is `E W^p=(m)_p E Z^p` conditional on each fixed completion, with the same factor for both placements; these are raw moments, with no division by `Γ(1+p)`. At a later stopping time it is the physical residual amplitude that restarts, and the moment gets the elapsed-time discount. The leaf upper `P(W>0)+p` follows by splitting at `W=1`, using concavity on `W≥1`, and using the degree-one Yule bound `E W≤1`. This takes expectation under the physical completion; it uses no return-time independence or new small-mass assumption for adverse tips.

The birth-window construction (`open/age-order/global-transport/birth-window/README.md`) implements the frozen lower as the chance to earn the terminal survival payoff. The serial improvement's first race stops at the old-parent change, monitored gain or loss, or the desired old-parent arrow; on an ordinary exit its physical continuation is the baseline, or the actual rescue after last loss. The favourable parallel tip starts with its own preliminary return. These two preliminary races are on different endpoint continuations. From focal degree one or two, at most 19 gains and at most 20 positive-stock losses occur before focal degree 20, since each loss spends an existing tracked leaf. Each positive-stock stage has tracked gain hazard at least one until its end, including when some exterior arrow brings the stage to an earlier end. In a pending return the focal degree stays fixed, the specified sender has degree at most 31, and the intensity of the specified send is at least `1/31` *throughout the race*, including after any sender-degree increase. Consequently the conditional survival of the entire timeout-free, capped race is at most `exp(−t/31)`; the gain or cap can end it sooner. One must not multiply this mean by a second count of the degree changes inside a race. At most 20 whole races and one initial serial race suffice, even with dummy zero waits, yielding `Eτ≤40+20·31+1=661`. This is an original-law bound on the stopped time with failures included. It is not a conditional mean bound after selecting success.

At a genuine degree hit the lower linear families have death rate `1/20` and terminal survival `1−20^(−j)`. Age 128 is measured from the stopping time. The chance a surviving family has fewer than `H=2^32` particles then is bounded by the *unconditional* probability of a positive count below `H`, at most `2H exp(−(19/20)128)` per family. At stock at least `H`, reset the actual mass bound at this age. Use the physical gamma-minus-exponential lower with the history at that age: for stock `l≥H`, `P(Gamma(l,1)≤H/2)≤4/H`, and `P(Exp(1)≥H/4)≤4/H`. These marginal estimates are union-bounded; survival need not be independent of either draw, and we do not use the gamma bound conditional on eventual survival. Thus the probability that the frozen survival path fails to yield residual mass above one is less than `2⁻²⁷`. On the good event its original contribution is at least `exp(−p(τ+128))`; the elementary exponential loss bound uses `E[(τ+128)1_good]≤789`. Adding the one-sided upper loss gives `η=2⁻²⁷+790p`.

Pay the signed debit only at its adverse owner, using `|s_u−s_v|≤s_adverse`. At a retained core the serial weight and old-leaf incident edge weight each total at most `3/2`, while directly paid core edges total at most one. Hence `4η` fits the frozen retained surplus. The producer's `6η` chain allowance for lengths through three also fits `δ=1/200000`; each later bare interior site has positive serial credit and an equal-coefficient new edge. For stars the per-core loss fits the frozen bounds (the low-star normalized credit has serial weight at least `1/2`); the seed has no edge debit and retains `G(2)−η>0`. Thus the old skeleton reserve is preserved without a loss multiplied by tree size.

**Separate improvement relayed to 1374.** For lengths through three the serial weight is at most `n_chain/3≤1`. At each boundary both endpoints are cores, so the absolute difference of their coefficients `s=w/d` is at most `1/3`. An internal edge joins bare degree-two sites with coefficients in `[0,1/6]`; there are at most two such edges. The chain loss is therefore at most `2η`, not `6η`. Keeping the same frozen budgets and the unchanged `4η` retained allowance extends this elementary transfer to `0<p≤3/(2·10^9)`; the retained slack at that endpoint is `3017247/13107200000000`. This is a review-side extension of the same proof, separate from the writer's pinned `10⁻⁹` claim. It does not certify intermediate fractional powers.

## Founder pair: independent fractional argument

At size three the shape is always the three-vertex path. By the seed symmetry label 2 is its centre with probability `1/2` and otherwise an exchangeable leaf; label 3 is always an exchangeable leaf. The conditional size-three amplitude factor is the same across both configurations. It is enough to show the raw centre moment is greater than the raw leaf moment at this one shape. Use the exact fractional harmonic boundary (`open/age-order/hub-ranking/fractional-middle/README.md`) with two tracked leaves: the centre lower is `Γ(1+p)e⁻¹(2+p)`. For an exponential `X` of mean one, `E log X≥−1` because the negative contribution on `(0,1)` has absolute value at most `∫_0^1−log x dx=1`. Hence Jensen gives `Γ(1+p)≥exp(−p)`.

For the leaf, the independently audited no-send probability `108/305` gives `a=197/305` as an upper for `P(W>0)`. The physical shifted-race mean bound gives `b<.387`. Hölder on `{W>0}` gives leaf moment at most `a^(1−p)b^p`; the endpoint `p=1` is the mean bound itself. The log ratio of `(2+p)exp(−1−p)` to `a^(1−p)b^p` is concave. With `e<11/4` and even the larger ceiling `b=.387`, both endpoint ratios exceed `16000/15609>1`. This argument uses zero only as an endpoint of the analytic ratio, not as a convention for `W^0`. It needs neither the tiny-power transfer nor the old-edge ledger. The founder gap in `Z^p` is half the centre–leaf raw gap divided by `(3)_p`.

## Replays

`taskset -c 8 nice -n 10 timeout 90s python3 -IBu round2/open/age-order/hub-ranking/small-powers/check.py` passed the frozen high/low, star, time, and mean-chain budgets, with these result lines:

```text
PASS small powers through 1/1000000000: endpoint difference loss < 7.97450580597e-07
PASS near-one powers from 14999/15000: retained loss <= 9/20000, chain loss <= 11/40000
PASS founder comparison for all 0<p<=1
OPEN p=1/2: degree-only fractional ledger does not close
```

```text
python3 -IB round2/checks/preplacement-small-powers-sol.py
PASS 661 mean retry clock; original maximum loss 41809377/52428800000000; 4η and 6η fit
PASS refined chain budget 2η; extended retained slack 3017247/13107200000000
PASS founder endpoint ratios; log concavity gives uniform ratio > 16000/15609
```
