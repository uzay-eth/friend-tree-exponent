> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# A conditional cutoff of 256 by tracking one leaf family

For every old tree of size `n`, the next two labels are strictly ordered at every `0<x<1` with **`(n+1)x≥256`**. The argument uses neither a good event for an infinite likelihood product nor a bound on the old parent's degree or stock. It concerns the upper range of the conditional problem only.

## An exact beta word on the physical tree

Start from any serial middle completion of size `m=n+2`. The tag `o` has one fresh leaf besides its old parent. Colour that fresh leaf red and every other vertex blue, including the tag. A newborn inherits its **initiator's** colour, irrespective of its attachment target. Colours never change. The sender is uniformly sampled in the physical jump chain. The colour word is therefore exactly Pólya `(1,m−1)`, with endpoint `Θ~Beta(1,m−1)`. Conditional on `Θ=θ` the word is iid Bernoulli(`θ`), and each step still chooses a fresh uniform sender *within the indicated colour*, then a uniform actual neighbour. This disintegration follows by the Pólya word probabilities: the rates of the two colours depend only on their counts. It does not tilt an artificial clock or discard any old-parent arrows.

Call a red vertex good if it is still a leaf of `o`, and bad otherwise. A red newborn is good if its sender is good (such a sender has only the tag as neighbour); any birth from a bad sender can only improve the good count relative to treating all its offspring as bad. A good leaf can be lost only by a send from `o`, and that send loses at most one. Thus bad fraction `F` has conditional drift at red steps at most zero. At a blue step, if `B` is the current blue count and `R` the red count, the chance to select `o` is `1/B`, and a good leaf is selected at most with probability one. The increment of `F` has conditional mean at most `1/(BR)`. Equivalently one can permanently mark lost families and all their descendants bad, in which case their fraction is a bounded submartingale. In either convention the following finite-time expectation bound suffices. Since the share of red good leaves is at most the tag degree fraction and `R/j→θ`, Fatou implies

```
E_θ[(θ−Z)_+] ≤ θ Σ_(f≥0) 1/(m−1+f) E_θ[1/R_f],                 (1)
```

where `R_f` is the red count **after** the red string just before blue jump `f` (index starts at zero). Indeed `(R/j−D_o/j)_+≤(R/j)F`; use bounded convergence for the red factor and Fatou for the deficit. The physical degree share has its usual a.s. limit under the unconditioned law, hence at beta-a.e. `θ`. The proof does not need convergence of an auxiliary bad fraction or independence of `F` and the old tree.

## Exact integral and the endpoint shift

Set `b=m−1`, `u=1−θ`, `z=bθ`. Then `R_f=1+Σ_(i≤f)G_i` for independent geometric waits with survival `θ^k`. In particular the loss charges include the reds **preceding** the first blue jump. By `1/(1+k)=∫₀¹ t^k dt` and geometric summation,

```
C(θ) := b θ Σ_(f≥0) E[1/R_f]/(b+f)
      = b u ∫₀¹ t^(b−1) log[(1−ut)/(u(1−t))] dt.              (2)
```

Insert `1/(b+f)=∫₀¹ t^(b+f−1)dt`: summing in `f` at a fixed generating-function argument `s` gives `u/(1−θs−ut)`; integrate `s` first. Every exchange precedes subtraction and is justified by positivity. If `T~Beta(b,1)`, then `E[-log(1−T)]=H_b`, the usual harmonic number, and concavity gives `E log(1−uT)≤log(θ+u/(b+1))`. Hence

```
C(θ) ≤ u [log(z+1)+1−log u].
```

Using this unconditional loss bound with a single deterministic endpoint shift would introduce a `log z` in the exponent. The needed improvement is a constant-probability prefix measured in **blue choices**, not in physical time or total jumps.

## Protect the family for a constant-probability prefix

Let `k=2b+1` and let `E` forbid choosing the tag as sender on the first `k` blue choices. Conditional on the *entire* colour word, including every geometric string, its probability is exactly

```
w = ∏_(0≤f<k) (1−1/(b+f)) = (b−1)/(3b).
```

Thus imposing `E` leaves the law of the word unchanged, even given `θ`. All reds are still good after these blue choices. Future senders given their colour and past retain their physical kernels. The bound (1) now starts at `f=k`. For `f≥1`, the generating function of the geometric sum gives

```
E_θ[1/R_f] = ∫₀¹ [u/(1−θs)]^(f+1) ds = u(1−u^f)/(θ f).
```

Consequently `E_θ[(θ−Z)_+ | E] ≤ u log(3/2)/b`: bound the remaining decreasing sum `Σ_(f≥k)1/[f(b+f)]` by its integral from `k−1=2b`. This also verifies why the first blue wait in (2) must be included. At `θ≥ x+(1−x)/b`, the gap `θ−x≥(1−x)/b≥u/b`. Markov gives `P_θ(Z>x | E) ≥ 1−log(3/2)`. The strict success event uses the failure event `Z≤x`, so equality at a threshold causes no issue. Integrating the exact beta tail gives the following lower bound for **every** serial middle, with no dependence on its old parent:

```
B = [(b−1)/(3b)] [1−log(3/2)] (1−1/b)^b (1−x)^b.
```

For `b≥256` the three bracketed factors are greater than or equal to `83/250`, `297/500`, and `9/25`, respectively. The second follows by truncating the alternating logarithm series after seven terms. For the third, `(1−1/b)^b≥exp(−b/(b−1))`, and `e<11/4`, `exp(1/255)≤255/254`. In particular their product `c` satisfies `800c²>4` with these rational floors. The imposed event has probability roughly one third, even as the endpoint goes to one, regardless of whether the old parent gains immediately.

## All adverse arrows, including zero stock

Use the [signed two-birth identity](../README.md#signed-comparison). With `A=Σ r_u/(D_u+1)>0`, the serial coefficients all buy `B`; serial tips cost at most `AU`; the absolute parallel edge coefficients sum to at most `A`, so only the adverse endpoint on each edge costs at most another `AU`. Here `U=Q_(m,1/2)(x)` follows from the locked urn at a leaf with permanent core parent. It needs **no other leaf siblings**. With `b=m−1`, bounding the beta density's `t^-1/2` by `x^-1/2` and applying Jensen to a `Gamma(m−1/2,1)` variable gives

```
U ≤ (1−x)^b sqrt(1−x)/sqrt(π (m−1/2)x)
  ≤ (1−x)^b /sqrt(π b x).
```

Hence `2U/B<1` for `bx≥256` by `π>25/8` and the rational bound `800c²>4`. This proves the claimed cutoff via the exact signed numerator lower bound `A(B−2U)`. It does not imply a single-digit sufficient threshold or an obstruction to one. It bypasses rather than strengthens the earlier all-parent endpoint likelihood estimate: the beta variable here is the exact **initiator-colour** proportion, not a locked attachment endpoint. In particular the argument loses the old parent's extra initial reciprocal mass and does not claim a sharp tail prefactor.

The companion checks physical red/blue kernels on old stars and on a tree with a nonzero adverse edge whose tip has zero other leaf siblings. It also checks geometric reciprocal sums and the explicit rational cutoff budget. Finite controls supplement the analytical argument; they do not estimate limiting tails by simulation.

```text
taskset -c 13 nice -n 10 timeout 12s python3 -IB round2/open/age-order/global-transport/sharp-large/check.py
PASS 47520 exact physical arrow histories: stars and nonzero zero-stock adverse fixture
PASS 32 reciprocal controls and exact cutoff budget Y=256
```
