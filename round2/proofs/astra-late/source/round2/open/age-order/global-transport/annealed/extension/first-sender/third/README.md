> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../../../README.md#claim-26) for exact scope and citations.

# At one third: an upper certificate for the first-family score

At the first uncovered older label `b=21`, the one-family conditional-Markov score of the [sub-half method](../joint/README.md), even evaluated with its exact marginal integrals and a generous physical first-miss credit, is **less than `.272` at `x=1/3`**. The chosen paired beta debit exceeds `.28`. This rules out closing that pair by improving the integration or inserting the first miss probability into the same score. It does **not** disprove pre-placement order or rule out a sharper residual, additional exposed families, actual-stock debits, or a different coupling. The positive-threshold coverage from the frozen notes is unchanged. The hubhood result at `b450a31` is a separate zero-threshold comparison.

## Exposure and the average over physical misses

Use the serial-middle colour setup. Let `θ` be the red endpoint, `u=1−θ`, `F=f` the first marked blue sender index, and `R` the red census just before it. Given the entire colour word and `F=f`, each of the preceding `f` blue choices draws its sender uniformly among blue vertices **other than the mark**. No red neighbour can be spoiled before the first marked send. Reds only add red neighbours there. All `R` reds are leaves and `R−1` is the sum of `f+1` geometrically distributed red waits conditional on `θ,F=f`.

Let `K=D−R` be the number of non-red neighbours and `K_h` its value just before the blue choice of index `h<f`. An original parent of degree at least two is permanent. Of the other `K_h−1` non-red neighbours, each can send to the mark at rate weight at most one; the old parent contributes at most `1/2`. Thus, even conditional on the full word and the pre-choice history,

```
E[K_(h+1)−K_h | word, F=f, past] ≤ (K_h−1/2)/(b+h−1).
```

Only adding a new blue neighbour can increase this count. Starting with `K_0=1`, the product telescopes and `E[K | word,F=f]≤1+f/[2(b−1)]`. This is an upper bound **conditional on the word**, not a claim of independence between stock and degree. At the exposed physical arrow the miss probability is exactly `K/(R+K)≤K/(R+1)`. Taking expectations in the correct order gives

```
P(miss | θ,F=f) ≤ (1+f/[2(b−1)]) E[1/(R+1) | θ,F=f]
                  ≤ min(1,(1+f/[2(b−1)]) min(1/2,2/f)),   θ≥1/3,
```

with `2/f=∞` at zero. The last step uses `1/(R+1)≤1/R`, the geometric reciprocal formula `E[1/R]=u(1−u^f)/(θf)` at `f≥1`, and `u/θ≤2`. This is an integrated physical upper credit for miss opportunities, not the worst-case `K≤f+1` substitution. A uniform lower credit exists only as a function of stock, as in the joint note.

## Which score is bounded above

Condition `θ>x`. The score uses `t=b(θ−x)/u` with survival `(1+t/b)^(-b)`, and uses `F=floor((b−1)(1−q)/q)` for independent uniform `q`. The one-family lower-payoff formula is

```
G_f(t,a) = (1−C_f(1+a/b)/(t−a))_+  for a<t, and zero otherwise,
C_f = b Σ_(j>f) 1/[(b+j)(j+1)] = b/(b−1) Σ_(k=f+2)^(f+b) 1/k.
```

We use the exact harmonic coefficient instead of its larger log upper bound. Assign a uniform dummy founder when the first arrow misses, using the real founder on a hit. Conditional on the pre-arrow history and destination, the exposed founder share has `Beta(1,R−1)` law (unit atom at `R=1`), so the dummy-family marginal `a=bθW/u` has uncapped survival `(1+a/b)^(-(f+1))` for purposes of this payoff, as before. The score examined here is `E[G_f(t,a) 1_hit+G_f(t,0) 1_miss]`. Whether or not the residual estimate can be sharpened on the selected subpopulations, this particular score is bounded above by

```
E[G_f(t,a)] + E[G_f(t,0) P(miss | θ,F)].
```

This inequality drops the overlap entirely, so it does not multiply independent stock and family averages to model the physical miss branch. Nor is it a bound on a scheme that integrates a new, stronger residual conditional on the full census or on the destination. The score itself is not asserted as a lower bound after selecting its hit and miss rows, which would require rechecking the residual conditional on that selection.

## Exact replay

Upper rectangles use the two finite `F` endpoints of each `q` bin, the right `t` endpoint with exact rational bin mass, and the left family-loss endpoint. The larger `F` simultaneously gives a smaller charge and a stochastically smaller loss. The conservative miss factor uses larger `F` in the non-red-stock factor and smaller `F` in the reciprocal factor. Family CDFs are rounded **up** at grid endpoints, with off-grid mass charged at most its appropriate envelope (`a≥12` pays zero on the retained `t≤12` bins). The `q≤1/512` and `t>12` tails cost their entire probability because the original score lies in `[0,1]`. The exact replay compares against the chosen debit using `π<22/7` and `(2/21)^(1/6)>2/3`.

```
nice -n 10 python3 -IB round2/open/age-order/global-transport/annealed/extension/first-sender/third/check.py
PASS exact upper averaging: hit < .232, miss credit < .038, discarded bins < .0021
PASS one-family certificate < .272 versus chosen beta debit > .28 at (21,1/3)
```

The practical next comparison is to tighten the residual using the **post-arrow census** `(p,r)`, or expose later victims and retain a physical or dominating stock variable at each new checkpoint. The fixed one-family score has an absolute deficit of more than `.008` even with this generous first-miss credit; repeating the same quadrature with more bins cannot close it.
