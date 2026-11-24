> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes for hubhood, with GPT-only audits of the accepted positive-threshold proof; threshold addenda use Claim 24 or pinned Fable finite transports where stated. See the [release index](../../../../../../README.md#claim-26) for exact scope and citations.

# Annealed promotion balance without a stationary law

**What closes.** For every `n≥2`, the probability of promotion on the very next jump after birth is strictly larger for label `n+1` than for label `n+2`, in the actual edge-seeded chain. This is an all-label result for that finite-delay statistic, not for hubhood. There is an exact stationary-free signed cohort criterion for hubhood and threshold tails. A natural strategy for proving it by nonnegative synchronized delay slices fails: the first pending-promotion slice after both births is strictly **negative on every old tree**, hence also after any old-jump or cohort averaging. Even cumulative first-promotion order can reverse conditional on a reachable completed placement. These obstructions do not exhibit an unconditional reversal of hub probabilities. Fable Claim 5’s finite full domination and the existing finite hub receipts remain separate coverage; this note does not extend their label range for hubhood. See the [transport audit](../../../../attacks/fable5/flow/README.md) and [base pre-placement argument](../README.md).

## 1. A strictly decreasing immediate-promotion statistic

On a physical old tree `T` of size `n`, let `s_u=r_u/[D_u(D_u+1)]`, `w_u=D_u s_u`, and `A(T)=Σ_u w_u`. Label `n+1` promotes on the following jump with conditional probability `A(T)/[n(n+1)]`: the birth rate at `u` is `r_u`, and its new leaf can then only be promoted by `u`, at rate `1/(D_u+1)`.

If the next leaf attaches to `u` of degree `d`, write `T^u` for the new tree. At `u` the rate becomes `r_u+1`; at every old neighbour `v` it falls by `1/[d(d+1)]`. The new leaf's own incoming rate is `1/(d+1)`. Thus

```math
 A(T^u)-A(T)=\frac{3d+4-2r_u}{2(d+1)(d+2)}
              -\frac1{d(d+1)}\sum_{v\sim u}\frac1{D_v+1}.\qquad\text{(1)}
```

The coefficient includes the new leaf's contribution `1/[2(d+1)]`. Since `r_u>0`, the increment is strictly below `3/[2(d+1)]`. With `P` the actual single old jump, `P A(T)<(1+3/(2n))A(T)`. The conditional probability for label `n+2` to promote on its next jump is `P A(T)/[(n+1)(n+2)]`, strictly below `(n+3/2)/(n+2)` times the preceding conditional probability for `n+1`. This uses the exact old-target law, and the conditional statement evaluates both forecasts from the same old tree. It still applies when `n=2` and both original vertices are leaves.

## 2. Exact synchronized first-promotion decomposition

Compare labels `a=n+1` and `b=n+2`. For a leaf `j` at size `k`, write `λ_j(T_k)=1/D_parent`; set it to zero after its first promotion. For each reachable `T_k` and pending leaf `j`, write `b_j(T_k^j;x)` for the future probability of `Z_j>x` when the next jump attaches a leaf to `j`. The physical dynamics thereafter use the actual postpromotion state; beta caps can bound this continuation but are not its value. The first-promotion decomposition, conditional on `T_n`, is

```math
 p_{a\mid T_n}(x)-p_{b\mid T_n}(x)
 = \frac{\sum_u w_u b_{a}(T^{u,a};x)}{n(n+1)}
 +\sum_{k\ge n+2}\frac1k\,\mathbb E\big[
 \lambda_a(T_k)b_a(T_k^a;x)-\lambda_b(T_k)b_b(T_k^b;x)\mid T_n\big].\qquad\text{(2)}
```

When a tag is no longer pending, its product in the sum is understood to be zero. `T^{u,a}` denotes the two-birth chain. For each summand, condition on the actual past and the next attached recipient. The two separate unsigned sums each converge by disjointness of the first-promotion times. This is also a derivation of the exact conditional rearrangement of discounted first-promotion forecasts: insert a deterministic discount for the destination size, or keep a conditional continuation `b_j`, but neither can be moved outside the expectation if it depends on the realized future state. The leading term of (2) is **all** serial middles, including those from first attachment at an old core.

For any deterministic cutoff size `m≥n+2`, an explicit stationary-free certification criterion is to compute the right-hand side through `k=m−1` using lower/upper intervals on the two continuation factors, then subtract the expected younger pending upper endpoint at `T_m`. For a pending younger leaf, an upper endpoint is `min(1,2/(J+1)) Q_{m+1,3/2}(x)`, with `J` other leaves at its actual parent and `Q` the locked beta cap (at `x=0`, use one for `Q`). This follows from the first-send bound and the earliest possible promotion size. The older residual can be dropped or retained using a lower charge. A positive certificate proves the physical adjacent gap at that `n,x`; a universal bound of this form would prove the all-label assertion. No such positive uniform certificate is claimed here. The formula also holds for finite cohort sums by summing (2) with nonnegative deterministic cohort weights. We use it next to isolate a contribution of the physical law on which cohort averaging cannot help.

## 3. The first common pending-promotion slice is strictly negative

For this calculation replace the future payoff `b_j` by one: count **first promotions**, not hubs. Look at the jump from size `n+2` to size `n+3`. If both births attach to the same old target, both pending rates agree. In a parallel placement at distinct old targets, the pending-rate difference is

```
H_0(u,v) = q_u−q_v,   q_u=1/(D_u+1).
```

In a serial placement the older tag is already promoted, and the younger pending rate is `1/2`; the difference is `−1/2`. Exact reversal of the birth order on distinct old targets therefore gives the *unnormalized* pending slice

```math
 -A(T)/2+E(T),\qquad
 E(T)=\sum_{uv\in E(T)}(s_v-s_u)(q_u-q_v).\qquad\text{(3)}
```

The denominator is `n(n+1)(n+2)`. Orienting the edge sum also writes `E(T)=Σ_u s_u Σ_(v~u)(q_v−q_u)`. Because `q_v≤1/2`,

```math
 E(T)\le A(T)/2-\sum_u w_u q_u,
 \qquad -A(T)/2+E(T)\le-\sum_u w_u q_u<0 .\qquad\text{(4)}
```

In contrast, the first slice (promotion during the second birth) is positive, `A(T)/[n(n+1)]`. Both signs remain strict after any old-tree averaging, in particular after a real old jump or any finite growing cohort of deterministic adjacent pairs. So nonnegative domination at **each common exact destination size**, or termwise paying adverse edges out of later pending promotion increments, cannot be justified even with the correct annealed law. This leaves open repayment using the earlier serial term, a signed combination of slices, or a proof on hub-weighted promotions. Worth distinguishing: (4) does not say the edge contribution alone is negative, and replacing `b_j` by one changes the statistic.

For decreasing deterministic destination weights `q_k`, summation by parts says the discounted difference equals the sum of cumulative first-promotion differences times the nonnegative drops in `q_k`, plus any terminal weight times the lifetime difference. Thus a negative exact slice is not an obstruction to *all* decreasing discounts. The cumulative difference through the first shared jump is in fact positive: the two contributions sum to `[(n+2−1/2)A+E]/[n(n+1)(n+2)]`, and `|E|≤A/2`. Cumulative order at all destination sizes requires another argument.

## 4. Why a positive repayment path cannot be imposed after placement

In the [reachable size-204 fixture](../../../../../../../../reviews/cold/astra/checks/age-order-sol.md#size-204-exact-hazards-and-what-conditioning-buys), the older and younger leaf are attached respectively to a large star and its degree-two arm. The older eventual first-promotion probability is bounded above by the same first-send bound `29851/980100`, whereas the younger has promotion probability at least `101/304`, by the three-clock race. This proves that the **post-placement** cumulative order reverses at a finite destination size (take the finite cumulative sums to their strictly ordered limits). It also reverses for all sufficiently shallow positive exponential discounts on that pending future, by bounded convergence. These are actual promotion probabilities, not beta proxies. They cannot be averaged by simply retaining a positive delay order for every realized post-placement history.

This last example conditions on the *completed* two-birth tree; it gives **no reversal** after averaging both placements from an old tree. Likewise there is no claim of negative cumulative order in the unconditional birth law. In (2) the necessary missing control is an upper bound on the annealed later deficit *relative to the actual all-middle continuation credit*. A negative slice does not quantify that total, and first promotion and eventual hubhood are different marks. A proof using a finite growing cohort would still have to control that unpaid total at the individual entrance rank rather than infer it from the sign of a weighted average.

## 5. An all-label result for genuinely discounted physical marks

Let `M_j` be the destination size of first promotion, with `∞` if there is none. Discount by **absolute destination size** (the same discount for the two labels). For nonnegative `q_l` with `q_(n+2)>0`, consider `E[q_(M_j) 1_(Z_j>x,M_j<∞) | T_n]`. The younger pending rate at size `n+2`, averaged over the two births, is exactly `P A(T_n)/(n+1)`. At all later sizes its rate with pending indicator can only decrease on every trajectory: its permanent parent's degree cannot decrease before promotion, and the indicator then kills the rate. Therefore, dropping the older positive later credits in (2) and capping the younger continuation by one gives the strict comparison whenever

```math
 \sum_{i\ge1} \frac{n+3/2}{n+1+i}\frac{q_{n+2+i}}{q_{n+2}}\le b_* ,\qquad\text{(5)}
```

where `b_*` is any positive uniform lower bound for the serial-middle continuation at the target threshold and on the given old tree. To see the normalization, divide the upper charge `P A(T_n)/(n+1)·Σ_(k≥n+2) q_(k+1)/k` by the serial unit `A(T_n) q_(n+2)/[n(n+1)]` and use the strict `P A(T_n)<(n+3/2)A(T_n)/n`. No empirical input, joint independence, or dropped serial middle is needed. If counting first promotions without hubhood, take `b_*=1`, which proves strict comparison on every old tree for geometric discounts `q_(l+1)/q_l≤1/2`. With hubhood, the fresh-leaf survival lower bound alone gives `b_*=1/2`, so the **hub-weighted** mark is strictly ordered for every old tree and ratio at most `1/3`. The stronger serial-middle bound from ee0c's separate conditional analysis can be substituted once certified; for example any bound at least `3/4` allows `3/7`. None of these proves strict ordering of the undiscounted hub mark.

At positive thresholds the younger upper charge in (5) can instead be multiplied termwise by `Q_(n+2+i),3/2(x)` (all future promotions happen on those destination sizes); the middle lower bound must still be at the actual threshold. With `q≡1` and `x=0`, the simple sum in (5) diverges, though the true pending terms are summable. In the bounded birth-scale range `x=y/n` the beta caps alone cannot give a usable uniform version of (5): for `n+3≤l≤2n` and any fixed `0≤y≤K`, `Q_(l,3/2)(y/n)≥(1-y/n)^(l−1)`, by comparison with `Beta(1,l−1)`. For `n≥2K` this is at least `e^(-4K)` throughout the interval, so the weighted sum still grows at least linearly in `n` for fixed `K`. A proof must use cancellation against the older pending cohort, or a sharper bound on the integrated younger residual than degree monotonicity and the beta discount. The synchronous negative slice (4) shows why a termwise sign argument cannot provide that cancellation.

There is another limitation of trying to pass from the discounted comparison to hubhood. A discount depending on **each tag's own age** has a different coefficient for the same destination size, so it cannot be substituted into the signed synchronized sum. Even for the common absolute discount, the inequality for ratios up to `1/3` supplies no sign near ratio one; the signed generating function is not a nonnegative series, by (4) already for the unmarked promotion test. Thus this is an all-label statement only for the explicitly discounted marks, not an Abel-limit proof of unweighted dominance.

## Exact replay

```
python3 -IB round2/open/age-order/global-transport/annealed/check.py
PASS 6 physical seeds, 51 true old-jump successors, both birth orders and discounted budgets
```

The direct side constructs directed donor rates from adjacency, checks both births with their original conditional probabilities, and reads the pending hazards on the resulting tree. The fixtures include the initial edge, paths, a star, a comb and a three-ray tree. This is not a sweep over old birth histories or a test of the ultimate mark by simulation.
