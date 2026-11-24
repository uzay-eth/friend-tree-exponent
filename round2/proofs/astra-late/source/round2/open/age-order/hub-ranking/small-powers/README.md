> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits, scoped Claude derivations and replay, and Gemini proof review and independent checks; prior two-birth identities and finite box certificates as cited. See the [release index](../../../../../../README.md#claim-27) for exact scope and citations.

# All-tree adjacent Mellin order below one

For the physical edge-seeded tree, the older of the next two labels has strictly larger conditional moment of its limiting degree share, on **every** old tree at size `n≥2`, for real exponents in either interval

```math
 0<p\le 10^{-9},\qquad \frac{14999}{15000}\le p\le1.
```

The strict unconditional comparison of labels 2 and 3 holds for **every `0<p≤1`**; the two founders have identical laws. The first all-tree interval extends the hub ledger to positive Mellin powers. The second is a direct transfer from the mean ledger; the independent [fractional-middle certificate](../fractional-middle/README.md) gives a wider signed interval near one.

## Common physical amplitude

At the common completion size `m=n+2`, restart the physical rate-one vertex clocks. For either endpoint let `W=lim e^(−t)D(t)` and `S=lim e^(−t)N(t)`, so `W=SZ`. Conditional on the starting completion, the size holding times have rates `j` depending only on size and are independent of the embedded graph sequence. Therefore `S~Gamma(m,1)` is independent of `Z` under that conditioning, and

```math
 \mathbb E W^p = \frac{\Gamma(m+p)}{\Gamma(m)}\,\mathbb E Z^p.
```

Use the same factor for every completion in the signed two-birth identity. Work with the **unnormalized** `W` moments, so the hub-probability comparison is their limit as `p↓0` and the mean ledger has the same units at `p=1`. No gamma factorization is used at a later random promotion, return, or stock hit. At such a time we use the physical remaining-time amplitude and its clock discount.

## Near zero: pay only the successful retry paths

For any leaf endpoint, `E W≤1` by degree-Yule domination from degree one. The elementary inequality `w^p≤1+p(w−1)` on `w≥1` gives

```math
 \mathbb E W^p\le \Pr(W>0)+p\,\mathbb E(W-1)_+\le \Pr(W>0)+p.
```

It is important that the expectation is under the original completion, with zeros omitted only in the indicator term. Thus every adverse frozen leaf upper can be increased by `p` without assuming anything about how late its first promotion occurs.

For the lower endpoint use the monitored construction from the [birth-window proof](../../global-transport/birth-window/README.md), including the initial favourable-tip return and the serial old-neighbour credit. It reaches degree 20 with a tracked stock, or fails. At most 40 positive-stock decisions occur (each gain advances degree, and each loss spends tracked stock); each wait has conditional mean at most one by the tracked birth intensity. A **whole** pending return race, including all sender-degree increases until cutoff, has conditional mean at most 31: while it is pending, the specified send has intensity at least `1/31`. There are at most 20 such races including the possible initial favourable-leaf promotion. The initial serial old-neighbour race, when used, has conditional mean at most one by the tracked leaf's rate-one clock. Thus elapsed time `τ` through monitored decisions satisfies `E τ≤40+20·31+1=661`, with zero waits padding terminated paths.

At a successful degree hit, the survival event of the terminal tracked families implements the frozen `1−20^(−l)` payoff. At age 128 afterwards, outside error less than `2^-27`, survival gives actual remaining mass greater than one. In detail, the probability of a surviving family remaining below stock `2^32` is at most `40·2^32·2^-121` across the at most 20 families, by the linear birth–death count. At stock at least `2^32`, the physical gamma-minus-exponential coupling costs at most `8/2^32` to ensure remaining mass greater than `2^32/4`. Use a union bound, with no independence of the success event from the final mass. On the good event, the original mass satisfies `W^p≥e^(−p(τ+128))`. Since `1−e^(−t)≤t`, the error against the frozen successful-event probability is at most `2^-27+789p`. This pays the logarithmic clock cost in expectation rather than requiring every retry to finish by a deterministic time.

Consequently the serial gap and either signed parallel difference lose at most `η=2^-27+790p` from their frozen one-sided bounds, including the adverse leaf upper. Lower bounds defined as zero require no tracked event. At a retained core, the serial weights of the core and its old leaves sum to at most `3/2`. The old-leaf edge weights charged there also sum to at most `3/2`, since `|s−t|≤s+t`; direct core edges can cost at most one further unit of weight at their adverse owner. Hence the error at each retained vertex is at most `4η`. For suppressed chains of lengths one through three use the conservative `6η` from the birth-window ledger: the sum of serial weights is at most one and at most four edge weights each are at most one. Longer chains add homogeneous interior sites with equal `s`, zero extra edge debit and positive serial credit. The exact [slack replay](../../global-transport/birth-window/check.py) gives residual at least `1/200000` at every retained core and greater slack on suppressed chains. At `p≤10^-9` both losses fit. Stars use the same per-core error against their earlier positive floors; the edge seed has only serial terms. The reduced-tree baseline `B+2A=.02968` is therefore preserved when there are at least two old cores, now on the raw Mellin numerator scale.

For founders, at size three the centre has two tracked leaves and survival probability at least `3/4` under the per-leaf death rate `1/2` comparison. Run directly to age 128: the error toward remaining mass greater than one is at most `4·2^32·2^-64+8/2^32<2^-27`. Against the elementary leaf upper `1−108/305` from the [founder proof](../../early/README.md), the centre still has advantage at least `127/1220−η`. Label 2 is half centre and half an exchangeable leaf, whereas label 3 is always an exchangeable leaf in this shape. Their moments use the same true factor at size three.

## From the mean side

For `p=1−q`, `0≤q≤1/15000`, a lower one-sided perturbation from the mean follows from `w^p≥w−q w log(w)` on `w≥1` and `w^p≥w` on `[0,1]`. Use `w log(w)≤w²/e` on `w≥1`. The physical degree-Yule domination gives `E W²≤6` for a serial middle or the seed centre and `E W²≤2` for a leaf, with the comparison based on their actual initial degree. Hence the lower losses are at most `(6/e)q` and `(2/e)q` respectively. On the upper side, for `0<w<1` integrate the derivative in the exponent to get `w^p−w≤q w^p(−log w)≤q/(ep)`; on `w≥1` there is no increase. Since `e>65/24` and `p≥14999/15000`, both `1/e` and `1/(ep)` are less than `3/8`. These bounds operate on the same physical continuations as the mean ledger, not on a newly normalized endpoint law.

It follows that the serial gap lower can decrease by at most `λ_s=(21/8)q`; each adverse parallel difference can increase by at most `λ_e=(9/8)q`. In the [mean ledger](../mean-ledger/README.md), at a retained core the total error is at most `(3/2)(λ_s+λ_e)+λ_e=(27/4)q`, by the same owner accounting as above. For each suppressed chain through length three the serial weights sum to at most one and each of its at most four edge weights is at most `1/3`, since both endpoints are old cores. Its loss is at most `λ_s+(4/3)λ_e=(33/8)q`. Longer chains again add positive homogeneous credit. At the lower exponent endpoint these losses are `0.00045` and `0.000275`, below the actual first-moment margins `0.0005` and `0.000291477`. Stars and the seed have ample surplus by the same estimate. For the founder comparison use centre `E W≥1` by the physical gamma-minus-exponential bound from stock two, and the exact mean first-send upper `U(2)` for the leaf. The centre minus leaf remains at least `1−U(2)−λ_s>0`.

The near-one argument is a uniform perturbation of the **physical first-send mean estimates**; it does not interpolate endpoint ordering in `p`. With the independent fractional-middle checker, the signed comparison still fails to close at `p=1/2`: its serial lower exceeds `.6240`, its degree-two-parent adverse upper exceeds `.4878`, and its all-degree tail precondition fails. This is a failure of the certificate, not a counterexample. The intervening all-tree range requires sharper signed comparisons, in particular ones that keep more information about endpoint leaf stock.

## Founder comparison for all fractional powers

Apply the [exact stock boundary](../fractional-middle/README.md#exact-stock-boundary-without-a-moving-cutoff) to the size-three centre with two tracked leaves. Its normalized lower is `e^-1(2+p)`. If `X~Exp(1)`, then `E log X≥−1` by bounding the negative integral on `(0,1)` by `∫_0^1 −log x dx=1`. Jensen gives `Γ(1+p)=E X^p≥e^-p`. Thus the raw centre lower is `(2+p)e^(-1-p)`. For the exchangeable leaf let `a=1−108/305` and `b=U(2)` be its probability and mean uppers. Hölder gives `E W^p≤a^(1−p)b^p`, since `W^p=1_{W>0}^{1−p}W^p` for `0<p<1`; at `p=1` use the mean bound itself. The log ratio of the centre lower to this leaf upper is concave in `p`, since only `log(2+p)` is nonlinear. At `p=0,1` the two centre bounds `2/e,3/e²` exceed `a,b`, respectively, already using `e<11/4`. The ratio is therefore strictly above one throughout. The size-three centre/leaf mixture for label 2 and common amplitude then give the claimed founder gap.

## Replay

`python3 -IB round2/open/age-order/hub-ranking/small-powers/check.py` replays the frozen endpoint and first-moment budgets read-only, then checks the rational loss and founder inequalities. No old-tree enumeration or expansion of a label cutoff is needed.
