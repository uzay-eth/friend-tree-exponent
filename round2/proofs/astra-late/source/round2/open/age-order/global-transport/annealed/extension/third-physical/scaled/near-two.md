> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../../../../README.md#claim-26) for exact scope and citations.

# A size window for the near-two comparison

The target is a uniform conditional comparison near a fixed value of `b x/(1−x)` below 2 as `b→∞`. This is separate from the finite tenth-threshold splice.

## A nonvanishing serial improvement

At a serial middle choose the tracked stock initially equal to one. At every later size `N`, the old host is still a core, so the physical degree of the tag is at least `j+1` at stock `j`. Its loss probability is therefore at most `j/[(j+1)N]`. Its gain probability is at least `j/N`, even if the old host contribution is ignored.

Couple the physical stock to an auxiliary stock with birth `j/N` and death `j/[(j+1)N]` up to reaching `J`; upon reaching `0` or `J` kill the auxiliary experiment and count zero. Before killing, the ordered coupling exists: at equal stock the physical tails dominate, and at adjacent stock the smaller auxiliary upward plus the larger physical downward probabilities are bounded by `(J+1)/N<1`. At larger distances they cannot cross. At the terminal size `⌈A m⌉`, pay the Poisson envelope to the surviving auxiliary stock. This gives a valid lower, even though the killed auxiliary payoff is not monotone at its cutoff.

For fixed `A,J` and `m x→y`, define the subgenerator `Q` on `1,…,J−1`, with off-diagonal birth `j`, death `j/(j+1)`, and the diagonal including both rates even at exits. The product of the auxiliary size kernels converges to `exp(log(A) Q)`: the sum of `1/N` over the window tends to `log(A)`, and the sum of `1/N²` tends to zero. The terminal binomial tails in `H` converge uniformly for `y` in a compact interval to the Poisson CDF at `Ay`. After normalization by `q^(m−1)`, the test becomes

```
e^y [exp(log(A) Q) h(Ay)]_1,
h_j(v) = e^(−1) Σ_{k<j} P(Poisson(v)≤j−k−1)/k!.
```

Unlike the one-jump gain, this improves the limiting coefficient by a fixed amount at `y=2`. The same argument allows any certified host contribution throughout the window, by adding the host degree to the local state. In particular, for small host degree this would be an additional credit rather than a prerequisite for the stock-window lower bound.

**Exact lower bound.** Take `A=4`, `J=40`. The [rational replay](near_two.py) gives a strict limiting coefficient above `.462` throughout `1.99≤y≤2`. The terminal Poisson CDF is bounded below at `4y=8`, and the outside factor is bounded below at `y=1.99`. Uniformization uses the nonnegative substochastic matrix `I+Q/40` and drops all terms after 160. Each row is rounded down. The positive-power exponential factors use Taylor bounds with a geometric tail, and `log(4)` is enclosed by the positive atanh series at `1/3`. Uniform convergence and the strict margin imply the same `.462` lower for all sufficiently large `m`, uniformly on the displayed interval; no explicit finite seam is asserted. For example `b x/(1−x)=399/200` puts `m x` in the interior for all sufficiently large `b`.

```
nice -n 15 python3 -IB round2/open/age-order/global-transport/annealed/extension/third-physical/scaled/near_two.py
PASS limiting serial coefficient > .462 uniformly on 1.99 <= y <= 2
```

## Pending-tip reduction

For an upper computation set size ratio `a=N/m`, and use time `t=log a`. Before promotion retain `(k,j)` with parent degree `k` and other leaf stock `j`. The physical weights in this time variable are promotion `1/k`, sibling loss `j/k`, and gain in `[j+1,(k+j+1)/2]`. The two gain endpoints remain valid without assuming that stock and host degree are independent. At promotion the upper boundary is the gamma survival of shape `1+1/k` at `a y`; this is the limiting beta survival with the degree at the actual promotion. At a terminal `a=A` or a degree cutoff the shape is `1/k`. For a uniform discrete proof, these gamma payments need outward bounds on a compact size/threshold box, and the physical recursion needs a uniform truncation/error bound. Taking a fixed number of births while letting `m→∞` cannot capture this effect: the chance of any state change in those births is `O(1/m)` for fixed `k`.

The formulas above do not yet prove a signed comparison below 2. The remaining work is an outward replay for both continuations and the signed ledger on a common interval. Exploratory floating-point values are not used as certificate inputs.
