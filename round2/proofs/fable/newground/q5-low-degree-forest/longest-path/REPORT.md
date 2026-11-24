LONGEST DEGREE-2 PATH: L_n ~ c log n/log log n WITH c = mu PROVED

**Author:** Fable's GPT-6.2 Sol Ultrafast proof author. **Pinned source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. **Review and dependencies:** Same-family GPT cold referee; Claude lead checked chain rates and first moments only. Uses Astra's mean bound and positive almost-sure limit; see [F7](../../../README.md#f7).

## Audit table

| Claim | Verdict / scope | Writer and checker |
|---|---|---|
| Longest degree-2 path, counted in vertices | Almost surely, `L_n log log n / log n -> mu`; more sharply, `L_n = u(log n) + O(1)` a.s., where `u log(2u/e) = mu log n` at `u=u(log n)`. Inputs: `E C_n = O(n^mu)` and almost-sure positive core growth from Astra | GPT family (this child), self-checked only; no independent family review yet |
| First moment | Exact birth-flux count for downward runs, exact clean-chain lifetime/race and matching factorial bounds up to polynomial factors on the expected number of downward length-k paths at size `n` on the extreme scale | GPT family, self-checked only |
| Largest retained component M_N, fixed cutoff | `M_2 = L_n + O(1)`; for every fixed integer `N >= 3`, `M_N = Theta(log n)` almost surely (order, no sharp prefactor) | GPT family, self-checked only |
| Monte Carlo | 500, 500, 100 observations at the three requested sizes; finite-size diagnostics | GPT family; exact labelled-tree oracle for implementation |
| Astra inputs | Sharp mean upper bound (Claims 2/3), positive fixed-size a.s. limit (Claim 6); cited as results, not reproved | Written by Astra (GPT); previous Claude-family and Gemini-family audits are recorded on `astra.md`; see also `attacks/neon-on-astra-6.md` |

I am a GPT-family model. No other model reviewed these new claims in this child thread. All logs below are natural. Time is continuous unless a size subscript `n` is used.

## 1. Heuristic and exact first-moment structure

The right starting population is the **promotions**, not the ordinary births. A first degree-2 vertex is made by a promotion (a leaf becomes a core), and the number of promotions is the number of cores, up to an additive seed term. After rooting at the centre of the size-three tree, a nonfounder degree-2 vertex has exactly one child. A downward degree-2 run must have passed through a clean pendant chain of every smaller length, even if the last child of its final vertex is no longer a leaf at observation.

Suppose the clean pendant chain has `r` degree-2 vertices, ordered from its top to its tip, and then one leaf. It extends at rate `1/2` (tip core rings and chooses the leaf). The top has an exterior parent of degree `d >= 2`. Directed sends along the `r-1` internal edges have total rate `r-1`; the tip leaf feeds its parent at rate `1`; the exterior parent feeds the chain top at rate `1/d`. These are exactly the spoilers; hence their total rate is `r+1/d`. A successful extension turns the old leaf into a degree-2 core and creates the next leaf, keeping the first `r` degrees unchanged. After the `k`th vertex of the run has been made, future extensions are harmless, but further receipts at any of the first `k` vertices are forbidden.

For the planted chain with inert parent, or with a fixed incoming boundary rate `a`, the probability that a newborn top ever produces a clean run of `k` is exactly

```math
 q_k(a)=\prod_{r=1}^{k-1}\frac{1/2}{r+1/2+a},\qquad 0\le a\le1/2.
```

For a physical parent one must account for the feedback on its degree. One equivalent construction gives each **directed edge** an independent rate-one marked Poisson process, and accepts a proposal if its uniform mark is at most the reciprocal of the donor degree. The total accepted rate at a donor of degree `d` is `d/d=1` and the chosen neighbour is uniform. Proposals on the parent-to-top edge with marks at most `1/d_0`, where `d_0` is the degree at top promotion, dominate actual incoming arrows, since degrees only increase. Removing all incoming arrows bounds the physical run probability above by the inert chain; treating every boundary mark at most `1/2` as a spoiler bounds it below by the `a=1/2` race (the sufficient event is defined from independent interior and boundary clocks). In particular

```math
 2^{-(k-1)}/k!\ \le\ q_k^{\rm physical}\ \le\ q_k(0)
 \le\ 2^{-(k-1)}/(k-1)!.
```

These statements are for ever producing a run, not for survival to the observation time. For the latter, while the first `k` degrees remain two, all sends along the internal edges spoil at rate `k-1`; the current bottom child sends into the run at rate `1/D_child` and the parent at rate `1/D_parent`. The exact conditional hazard is therefore `k-1+1/D_child+1/D_parent`, bounded between `k-1` and `k+1/2` (the parent is a core). In particular, for a post-formation wait `w`, conditional survival is between `exp(-(k+1/2)w)` and `exp(-(k-1)w)`, uniformly in the formation history. More exactly, suppress every directed arrow whose recipient is one of the first `k` vertices, and keep all other transitions with their physical rates. Then

```math
 S_k(w;T,u)=e^{-(k-1)w}\,\mathbb E^{\mathrm{supp}}_{T,u}
 \exp\!\left[-\int_0^w\left(\frac1{D_{\rm child}(s)}+\frac1{D_{\rm parent}(s)}\right)ds\right].
```

This is the usual killed-chain formula: give each ordered edge its own independent unit-rate proposal stream and thin by the reciprocal donor degree. All suppressed streams are then independent of the suppressed trajectory. The expectation is not in general an exponential of a constant, because either exterior endpoint can gain edges without breaking the run.

For a more explicit exact first moment, for `k>=2` let `R_k(t)` be the number of downward paths of `k` degree-2 vertices whose top is not the founder. Let `A_r(T)` be the set of nonfounder tops of clean downward pendant chains of `r` degree-2 vertices ending at a leaf. Then

```math
 \mathbb E R_k(t)=\frac12\int_0^t \mathbb E\!\sum_{u\in A_{k-1}(T_s)}
            S_k(t-s;T_s,u^{+})\,ds.
```

Here `u^+` denotes the configuration just after extension, `S_k` is survival of the first `k` degrees in that configuration. This is the compensator of the unique threshold-creation event for each `u`; an earlier spoiler cannot be undone. In discrete time the exact analogous expression is the sum for `m=3,...,n-1`, replacing `ds/2` by `1/(2m)`, `T_s` by `T_m`, and `S_k` by the `(n-m-1)`-step survival probability in the post-extension chain. For the total path count, add the founder term: for each split `i=0,...,k-1` count the indicator that the founder and the prefixes of lengths `i` and `k-1-i` on its two original sides all have degree two (an absent prefix fails its indicator). The prefixes, when they qualify, are unique.

For lower bounds in a unit window require exactly `k-1` extensions, no other rings by any current fringe vertex (including outgoing sends from the root to its exterior parent), and no boundary proposals. The probability density for this interior event at ordered times `0<s_1<...<s_{k-1}<1` is

```math
  (1/2)^{k-1} \exp\!\left[-\sum_{r=1}^k(r+1)(s_r-s_{r-1})\right],
 \qquad s_0=0,\ s_k=1.
```

No boundary proposals independently has probability at least `e^{-1/2}`. The density integrates over a simplex of volume `1/(k-1)!`. Thus, unlike allowing age to grow with `n` in a no-spoil calculation, the fixed unit window gives `exp(-O(k))/(k-1)!`; the ever-hit race gives a matching upper bound up to exponential factors. A path can straddle the founder, giving two downward pieces in that single case; this has no effect on the constant. Putting `k=gamma t/log t` yields an expected count on the scale `exp((mu-gamma+o(1))t)` and predicts `gamma=mu`.

## 2. Almost-sure upper bound

Start at size three (re-rooting physical time there). Degrees only increase. For a path through the founder that uses `i` vertices on one initial side and `k-1-i` on the other, each side must have ever hit its required depth while the founder was still degree two. Freeze the founder at degree two by suppressing sends into it. Its marked stream splits into independent rate-`1/2` streams for the two initial sides. These sides are independent when the founder is frozen. A side begins as a leaf, which is eventually promoted by its founder stream; after that promotion suppress all further incoming sends from the founder on that side. This only removes spoilers before the first hit of a required clean length: the interior trajectory agrees up to that first hit whenever the physical run succeeds. The side's probability of ever reaching `j` is at most `q_j(0)` (`q_0=1`). Thus the sum of founder split probabilities, bounding both the ever-hit probability and the expected number of founder paths at any observation, is at most `k^2 2^(k-1)/(k-1)!`: use `q_j(0)<= (j+1)/j!` and the binomial formula. Counting both sides is necessary; a spoiler can shorten one side while growth on the other continues, so a single uninterrupted race in total founder length is not the right count.

For every nonfounder top, `q_k(0)` bounds the probability of ever attaining the threshold conditional on its promotion history, even if subsequently observed at a stopping time. Rooting at the founder, any all-degree-2 path not through it must run straight down: a turning top would have two children and a parent, hence degree at least three. At the deterministic size horizon `m`, conditional expectation at each promotion time therefore gives

```math
 \Pr(\exists n\le m:L_n\ge k)\le O(m^\mu)\,q_k(0)
                   +\exp(O(k))/((k-1)!).
```

One can use indicators for tops promoted by `m` and ever successful, since `n_prom <= m` is known at promotion. For `m=2^j`, `k=ceil((mu+eps) log m/log log m)`, the bound is `exp(-eps log m + o(log m))` and is summable. Using the next dyadic horizon deals with nonmonotonicity of `L_n`. This proves the a.s. upper bound.

## 3. Almost-sure lower bound: disjoint planted trials

At time `t` write `H(t)` for the number of nonfounder terminal cores with one leaf and one core parent (so degree two). We need exponentially many, and get them from recent promotions without a terminal-core frequency theorem. In a unit window starting at `b`, give every freshly promoted vertex a token as long as it stays terminal of degree two. A token is lost if its core receives an attachment (rate `1+1/D_parent<=3/2`) or its leaf is promoted (rate `1/2`). No event destroys two old tokens: in the former case only the recipient leaves this state, and in the latter only the recipient's parent does. A promotion creates a new token even when the same event destroys the parent's token.

Realize each token's losses by thinning its own independent rate-two proposal stream started at birth; label streams by order of promotion. Split accepted proposals into the possible loss events proportionally to their rates, and run the remaining events at their remaining rates. The rates are exactly the physical transition rates (partitioned by lost token, or none), so the tree has the correct law. Every token whose first proposal wait exceeds one must remain terminal until the window ends. The first `M` proposal waits in the preassigned sequence are independent `Exp(2)` variables unconditionally, including unused labels. If `J` promotions occur in the window, then

```math
 \Pr\{J\ge M,\ H(b+1)<M e^{-2}/2\}\le e^{-b_0 M}
```

for an absolute `b_0>0`, by the binomial lower-tail bound. This is not a binomial claim conditional on `J`.

The size process is a rate-one Yule process per vertex; `e^{-t}|T(t)| -> Y>0` a.s., hence by Claim 6 `e^{-mu t}C(t) -> X Y^mu>0` a.s. It follows that `J=C(t)-C(t-1)>=floor(exp((mu-delta)t))` eventually a.s. for `0<delta<mu`. Taking this deterministic threshold and summing the preceding error proves `H(t)>=c_0 exp((mu-delta)t)` eventually a.s. along integer times.

At time `t`, take the terminal cores and their unique leaves, with distinct, disjoint fringes. No tested core is the parent of another tested core, though several fringes may share a boundary parent. Use fresh increments of the directed-edge clocks. On every tested incoming edge retain the proposals with mark at most `1/d_0`, where `d_0` is the donor degree at time `t`. They form independent rate-`1/d_0` streams on different edges, including edges with a common donor, and dominate all actual incoming attachments. Prescribe the sufficient events within each fringe using separate fresh streams on new edges, indexed by local order of birth. Require no boundary proposals for the whole interval `[t,t+2]`. In each fringe require `k-1` clean extensions in `[t,t+1]` and require no other outgoing arrows from its current vertices up to `t+2`. The density from §1, with horizon 2 and ordered creation times in `[0,1]`, yields a uniform conditional success probability at least

```math
 p_k:=e^{-1-2(k+1)}\,2^{-(k-1)}/(k-1)!.
```

These sufficient events are independent conditional on the initial tree; outside changes can alter the effective thinning but cannot affect a sufficient event. The conditional probability of no success is at most `exp(-H(t)p_k)`. Set `k=floor((mu-eps)t/log t)` and choose `0<delta<eps`; on the good trial-count event this conditional error is summable in integer `t`. A success preserves the run throughout `[t+1,t+2]`. Borel–Cantelli thus covers **all** sufficiently late times, rather than just a lattice (the success event need not be independent across overlapping windows). At the size epoch `tau_n=log n - log Y + o(1)` a.s. this is the claimed a.s. lower bound. Countably many positive epsilons suffice. The earliest finite seed times do not matter.

### Matching first-moment scale

At deterministic physical time `t`, run the terminal-chain trial from `t-1`, requiring success through `t+1` as before. The token bound and positivity of `X Y^mu` imply `E H(t-1)>=c e^{mu t}` for a constant `c>0`: take `M=floor(a e^{mu t})` with `a>0` small enough that the promotion count exceeds it with probability bounded below. The expectation upper bound uses the count version of §2 (sum the same ever-success indicators) and `E C(T(t))=O(e^{mu t})`: Yule holding times are independent of the size-indexed tree sequence, so Jensen's inequality and `E C_n=O(n^mu)` suffice. Thus with `k~gamma t/log t` the logarithmic path expectation divided by `t` tends to `mu-gamma`.

The corresponding discrete bounds also hold on the extreme scale. Start at size `floor(m/2)` and use the loss proposals from there. The time to size `m` is at most 1 with probability `1-e^{-Omega(m)}`, since the future size at time 1 is a sum of `floor(m/2)` independent geometric Yule sizes with mean `e` (Chernoff). Claim 6 gives at least `floor(a m^mu)` promotions by size `m` with probability bounded below for suitably small `a>0`. On that event and a wait of at most 1, all the first `floor(a m^mu)` tokens whose waits exceed 1 survive until size `m`; the same binomial bound, without conditioning on promotion count, gives `E H(T_m)>=c m^mu`. Now take `m=floor(n e^{-3/2})`. The elapsed time between sizes `m` and `n` lies in `[1,2]` except with probability `e^{-Omega(n)}`, by the two geometric Chernoff bounds. Conditional trials given `T_m` give expected success count at least `p_k E H(T_m)`. On the exceptional timing event subtract at most `m` successes; on its complement all successful distinct tops give distinct paths at size `n`. In particular for `k~gamma log n/log log n`, the expected total length-`k` path count has logarithm `(mu-gamma+o(1)) log n`. More explicitly it is between `c n^mu p_k - n e^{-Omega(n)}` and `O(n^mu)q_k(0)+k^2 2^{k-1}/(k-1)!`. These bounds do not assert a sharp prefactor or an extreme-value limit law beyond the leading ratio.

### Sharper localization: bounded additive error

Let `F(u)=u log(2u/e)` and let `u(x)` be the large solution of `F(u(x))=mu x`. Then

```math
 L_n=u(\log n)+O(1)\quad\text{a.s.}
```

In particular, with `x=log n` and `ell=log x`, this implies

```math
 L_n=\frac{\mu x}{\ell}
     +\frac{\mu x}{\ell^2}(\log\ell+1-\log(2\mu))
     +o(x/\ell^2)\quad\text{a.s.}
```

**Short observation window that costs only polynomially.** Start from a terminal fringe and suppress its incoming boundary arrows. The probability of ever producing `k` is exactly `q_k(0)`. Given success, the waiting times at lengths `r=1,...,k-1` are independent `Exp(r+1/2)` variables, since in each exponential race the waiting time and the winning event are independent. Their sum has mean at most `1+log k`. With conditional probability at least `1/2` the formation time is in `[0,2(1+log k)]`. Partition this range into bins of width `1/k`, extending the last bin to full width. Some deterministic bin `[a_k-1/k,a_k]` has conditional probability at least `1/(2 ceil(2 k(1+log k)))`. Require formation in this bin and no accepted outgoing arrows from the `k+1` current fringe vertices for the next `2/k` after formation. The latter conditional probability is `exp(-2(k+1)/k)`. The run then survives throughout `[a_k,a_k+1/k]`. In the physical tree additionally require no incoming boundary proposals with marks at most `1/d_0` up to `a_k+1/k`. This independently costs a factor at least `exp(-(a_k+1/k)/2) >= c/k`. Thus for all `k>=2` the sufficient event has conditional probability at least

```math
    p'_k \ge c_1 q_k(0) k^{-3},\qquad a_k\le 2(1+\log k)+1/k,
```

and the sufficient events for disjoint terminal fringes are independent as above. To see why outgoing boundary sends during formation are harmless: their only effect is on the exterior parent, and all incoming changes during the trial are ruled out by the boundary event. The interior trajectory ignores outgoing sends except in the final short silence requirement.

**Upper bound.** The elementary factorial inequalities in §1 and Stirling give `log q_k(0)=-F(k)+O(log k)`; more precisely the discrepancies are between `-(1/2)log k+O(1)` and `(1/2)log k+O(1)`. At dyadic horizon `m` use `k=ceil(u(log m))+20`. By `F'(u)=log(2u)`, the nonfounder term of §2 is `O(k^(-19))`; the founder term decays as `exp(-mu log m+O(k))`. Both are summable in the dyadic index. Transfer to all sizes using the next horizon. Implicit differentiation gives `u'(x)=mu/log(2u(x))`, so the horizon change costs `o(1)` in this comparison.

**Lower bound.** In each integer interval `[j,j+1]` set `k=k_j=floor(u(j))-20`, for large `j`. Use `k` windows of length `1/k` with left endpoints `j+i/k`, `0<=i<k`. For each use a trial start time `s=j+i/k-a_k`. The token estimate applies on `[s-1,s]` with threshold `M_s=floor(e^{mu s}/s^2)`. The exceptional probabilities `exp(-b_0 M_s)` are summable over these polynomially many start times per `j`. Also `C(s)-C(s-1)` divided by `e^{mu s}` has a strictly positive limit as *real* `s->infty`, so eventually every such window has at least `M_s` promotions. Thus eventually `H(s)>=c_0 e^{mu s}/s^2` at all these starts a.s. Using `s>=j-2(1+log k)-1/k`, the trial count is at least `c_2 e^{mu j} k^{-2mu}/j^2`. By Stirling and `F(u(j))-F(k)>=20 log(2k)`, the product of the good trial count and `p'_k` is at least `c_3 k^{20-4-2mu}/j^2`. Since `j=O(k log k)`, this grows faster than `k^2`. Hence the conditional probabilities of failure on the good-count event, bounded by `exp(-H(s)p'_k)`, are summable over the mesh. Independence between different start times is not required. The event at each start covers its full observation window, giving `L(T(t))>=k_j` for all `t` in all sufficiently late intervals. At the size epochs `t=log n+O(1)` a.s., this differs from `u(log n)` by only a bounded additive amount.

For the expansion use `log u(x)=ell-log ell+log mu+o(1)` from the first-order asymptotic and divide `mu x` by `log(2u(x)/e)`. No claim is made about convergence of the bounded residual.

**Sharper first moment at size `n`.** This also replaces the crude fixed-window lower bound by a polynomial-loss version. Let `m=floor(n exp(-a_k-1/(2k)))`; at the extreme scale `m>=c n/k^2`. Conditional on `T_m`, the sufficient trials give expected count at least `p'_k H(T_m)`. The holding-time sum from size `m` to size `n` is within `[a_k,a_k+1/k]` except with probability `O(exp(-c m/k^2))`. Indeed its mean differs from `a_k+1/(2k)` by `O(1/m)`, and the sum of the squared scales of its independent exponentials is `sum_{r>=m} r^-2=O(1/m)`. For the centred sum its log-mgf is bounded by `C z^2/m` for `|z|<=m/2` by the expansion of `-log(1-z/r)-z/r`; choose `z` of order `m/k` with a small constant. Do not condition on timing to assert independence of successful trials; just subtract `m` times the exceptional probability. Using `E H(T_m)>=c m^mu` already proved above, this gives the bounds

```math
 c n^\mu q_k(0) k^{-3-2\mu}-O(n e^{-c' n/k^4})
 \ \le\ \mathbb E R_k(T_n)\ \le\ C n^\mu q_k(0).
```

On the extreme scale the timing remainder is negligible and the founder term, if included, is smaller than the lower bound. Thus the logarithm of the expected total length-`k` path count is `mu log n - F(k) + O(log k)` for `k` of order `log n/log log n`.

## 4. Largest component at fixed cutoff

Here `M_N(n)` counts ordinary vertices in the largest connected component of the induced subgraph on original degree at most `N`. For `N=1`, it is 1 after size two. For `N=2`, each nonsingleton component is a run of degree-2 vertices plus at most two original leaves. Thus `L_n <= M_2(n) <= L_n+2` (with the singleton case allowed).

### Upper bound for `N>=3` (also holds at `N=2`)

We show an exponentially small conditional bound for creation of a large component in a fixed time block, uniform in the old state apart from its vertex count. Fix `N>=2` and put `h=16N^2`. At a deterministic start time `s`, call all current ordinary vertices old. Work at the **end** of the block `[s,s+h]` first. Every component of retained vertices at that time either contains a connected set `S` of old vertices or lies entirely in one new branch from an old anchor (follow the unique path toward the old tree). If it contains old vertices, the old intersection is connected, and its new vertices all descend within the union of new branches attached to `S`.

Use independent rate-`1/N` mandatory proposal clocks on each directed old edge. They must be realized as attachments while the donor's degree is at most `N`. For every fixed connected old set `S` of size `j` with original degrees at most `N` at time `s`, the `2(j-1)` oriented internal clocks are independent. If `S` is retained at time `s+h`, their total marks (each raises a degree in `S`) cannot exceed `Nj`. For `j>=2` their Poisson mean is `2(j-1)h/N >= jh/N=16Nj`; using `P(Pois(u)<=u/2)<=exp(-u/8)` gives probability at most `e^{-2Nj}`. For `j=1` just use 1.

For the new volume use Cauchy–Schwarz and **marginal** bounds; old-degree survival and new-branch growth need not be independent. The retained component with old vertices `S` can receive at most `Nj` new branches directly rooted off `S`, because each injection of a new branch adds one to the degree of a vertex in `S`. Track only the first `Nj` injections at vertices in `S`, and give each potential branch an independent dominating Yule trajectory, preassigned in order of injection. Each such new branch has just one old boundary vertex. When its size is `m`, all donors in it contribute at most `m` to its birth rate; the boundary donor contributes at most 1. Thus the branch, which has size 1 at injection, is dominated by a rate-one Yule process begun with **two** individuals. For a simultaneous coupling, run each independent Yule trajectory on its own age clock from the corresponding injection. At each proposal in that trajectory, thin by the ratio of the branch's actual birth rate to the current Yule population, and choose the actual recipient proportionally to its attachment rate. This is valid inductively: the Yule population is always at least the current branch size plus one. Distinct branches have disjoint recipient sets, so all other births can run at their remaining rates, restoring the exact physical rates. Random injection times do not affect the pathwise bound by the sum of the first `Nj` independent *full-horizon* Yule populations (including unused trajectories). No independence from the old mandatory count is required; each marginal bound can be established in a separate realization of the same physical law.

Let `F(λ)=E exp(λ Y_h^{(2)})`, finite for small positive `λ` since `Y_h^{(2)}` is a sum of two geometric variables. Choose `λ>0` small enough that `λ+N log F(λ) < 1` (possible by continuity at zero). For each fixed `S`, the probability it is retained and has at least `m` vertices including all its new branches is bounded, by the geometric mean of the two marginal bounds, by

```math
  e^{N} e^{-Nj} \exp[-\lambda m/2 + j(\lambda+N\log F(\lambda))/2].
```

For `j=1` the prefactor compensates for the absent old-edge bound. The number of size-`j` connected sets of old vertices of degree at most `N` at the start is at most `|T(s)| N^{2(j-1)}` by encoding each in a depth-first walk. We have `2 log N < N-1/2` for all `N>=2`: it holds at 2 and `N-1/2-2 log N` is nondecreasing thereafter. Therefore summing over `j` gives an upper bound `C_N|T(s)|exp(-λm/2)` for the old-containing component event.

For components disjoint from the old tree, an anchor can inject many distinct branches (if its degree is high). Count each new branch at its birth in the block, then condition on its injection history: with remaining time at most `h`, its full size is dominated by `Y_h^{(2)}`. The conditional expected number of injections at any old anchor or all old anchors is bounded by the expected total births in the block, `(e^h-1)|T(s)|`. Compensation or summing over birth epochs and using the uniform conditional tail therefore bounds the all-new component event by `|T(s)|(e^h-1) F(λ)e^{-λm}`. Together, conditional on the old configuration, the two bounds imply `P(M_N(s+h)>=m | T(s)) <= C_N |T(s)| exp(-λm/2)`.

For a simultaneous-time upper bound use windows with `s=j-h` at integer `j>=h`, ending with observation times in `[j,j+1]`. The old mandatory survival bound still uses `[s,j]`; new branch injections and growth can last `h+1`. Choose `λ` with the latter horizon in `F`. The bound `Nj` on injections on the old set uses retained degrees at the **actual** observation time; to avoid confusion, call the cardinality `r`. Keep the first `Nr` injected branches and omit later ones for that set; at any observation when that set is retained, every relevant injection is among the first `Nr`. The all-new component case uses domination of the increasing total branch volume at `j+1`. Thus the same exponential bound applies to `sup_[j,j+1] M_N`. Since `E|T(j-h)|=3e^{j-h}`, choose threshold `m=K_N j` with `λK_N/2>1`, and apply Borel–Cantelli. Transfer by `tau_n=log n+O(1)` to get `M_N(n)=O(log n)` a.s.

### Lower bound for `N>=3`

In any terminal fringe from §3, produce many vertices all of degree at most three by the following branching recipe. A virgin core `x` has a single leaf `y`. Require that `y` ring once and feed `x`, and then that `x` promote each of its two leaves, in a specified order. This makes exactly two virgin core children (each with one leaf and a parent), and leaves `x` of degree three thereafter. All other intra-fringe rings are forbidden. Full binary trees of these cores are connected and survive the `N>=3` cut, whether or not the exterior parent does.

There is a common unit time window for a full binary core tree of depth `d` with at least `2^d` vertices, with success probability `>=exp(-A 2^d)` for an absolute `A`. Allocate successive phase lengths `ell_i=2^{i-d}` for `i=0,...,d-1`, which sum to less than 1. For each of the `2^i` virgins at level `i`, schedule its three required rings in successive thirds of the phase, and forbid all other intra-fringe rings until time 2. The required selection probabilities at the parent in the second and third thirds are `1/3` each. On these trajectories the fringe has at most `3*2^d` vertices (there are `2^{d+1}-1` cores and `2^d` final leaves). Thus using the union of Poisson clock densities and integrating over all allowed triplets for each virgin gives a lower bound

```math
 e^{-1-6\cdot 2^d}\prod_{i=0}^{d-1}\big[(\ell_i/3)^3/9\big]^{2^i},
```

where `e^-1` is a lower bound for no boundary proposals for time 2. To justify multiplication within a phase, prescribe the ordered sequence on each individual virgin (different virgins plus their leaves are disjoint), integrate over each virgin's three disjoint thirds, and allow every interleaving across virgins. The absence of all other rings contributes at least the displayed total hazard factor, independently of interleaving. The displayed product is `>=exp(-A 2^d)` since `sum_{i<d} (d-i)2^i <= 2*2^d`. This removes the factorial penalty that was unavoidable for a single pendant chain.

Take the largest `d` with `2^d<=a t`, for fixed `a>0` small enough that `Aa<mu/2`. The same conditional independence across `H(t)` terminal fringes as in §3 and its a.s. lower bound show that, except on a summable set of failures on the good trial-count event, one success occupies `[t+1,t+2]` with at least `at/2` vertices. This yields `M_N(n)=Omega(log n)` a.s. for `N>=3`, with the same lower bound for every such `N`. Countable intersections give the assertions for all fixed integer cutoffs on a single probability-one event.

## 5. Reproduction and measurements

Run from `round2/children/q5x`:

```sh
taskset -c 13 nice -n 5 g++ -std=c++17 -O3 -Wall -Wextra -Wpedantic sim.cpp -o sim
taskset -c 13 nice -n 5 ./sim check
taskset -c 13 nice -n 5 timeout 550s bash -c './sim 7 100 6201172026 > sample7.csv && ./sim 6 400 1217489 > sample6.csv'
taskset -c 13 nice -n 5 round2/.venv/bin/python -B summarize.py
taskset -c 13 nice -n 5 g++ -std=c++17 -O1 -g -fsanitize=undefined -fno-sanitize-recover=undefined sim.cpp -o sim-ubsan
taskset -c 13 nice -n 5 ./sim-ubsan check
taskset -c 13 nice -n 5 timeout 550s bash -o pipefail -c './sim 7 100 6201172026 | cmp -s - sample7.csv && ./sim 6 400 1217489 | cmp -s - sample6.csv'
```

The compiler issued no diagnostics. Both the regular and UBSan-labelled-tree oracles printed `PASS: 12800 explicit labelled-tree states; donor/neighbor ranks, leaf multiplicities and all four component scans agree`. The two raw CSVs had 301 and 801 lines respectively, including headers. Both seeded datasets replayed byte for byte (`cmp` exited 0 with no output). The first validation attempt did catch an order mismatch in the oracle: it originally sorted core children by their adjacency-list order, whereas the compressed tree stores them in order of promotions. Sorting the oracle's core-neighbour list by core index fixed the mismatch; simulation logic was not modified. The summary checks sample counts, all observation nodes, and the deterministic `L <= M2 <= L+2` relation.

The simulator stores the core tree and weighted leaf counts in a dynamically appended Fenwick tree. Choosing an ordinary leaf donor chooses its parent with probability proportional to the parent's leaf count. Choosing a core donor uses a uniform index over its original degree; if the index is a leaf, that leaf becomes a new core with one newly born leaf. The initial edge produces the unique size-three unlabelled tree, which is where the simulator starts. Dynamic programs on the resulting rooted core tree measure the longest path in **vertices** and ordinary volume of retained components at three cutoffs.

| n | R | L histogram (length: count) | Mean L ± MC SE | Median | `mu log n/log log n` | Mean / leading | `log log n` |
|---:|---:|---|---:|---:|---:|---:|---:|
| 100,000 | 500 | 1:1, 2:81, 3:314, 4:96, 5:8 | 3.058 ± 0.029 | 3 | 2.666 | 1.147 | 2.443 |
| 1,000,000 | 500 | 2:13, 3:231, 4:218, 5:37, 6:1 | 3.564 ± 0.030 | 4 | 2.977 | 1.197 | 2.626 |
| 10,000,000 | 100 | 3:14, 4:70, 5:15, 6:1 | 4.030 ± 0.058 | 4 | 3.280 | 1.229 | 2.780 |

The means of `M_2,M_3,M_4` were `(4.052,7.160,10.180)`, `(4.556,8.360,12.374)`, and `(5.020,9.560,14.510)` respectively. Measurements across sizes share trajectories, so the table is not three independent ensembles. The means are 15–23% above the leading approximant, and the ratio moves farther above 1 across this size range. `log log n` rises only from 2.443 to 2.780 between the smallest and largest sample. Evaluating the sharper deterministic centre `u(log n)` gives `5.000, 5.553, 6.084`, about two vertices above the sample means at these sizes; this illustrates how large a bounded additive offset can be relative to the path length. No finite-size correction was fitted.

```sh
taskset -c 13 nice -n 5 round2/.venv/bin/python -B -c 'from math import log,e; from scipy.special import lambertw; mu=.56572366; print([(n, round(mu*log(n)/lambertw(2*mu*log(n)/e).real, 3)) for n in (10**5,10**6,10**7)])'
```

Output: `[(100000, np.float64(5.0)), (1000000, np.float64(5.553)), (10000000, np.float64(6.084))]`.

## 6. Limits and dead ends

* Using `n` ordinary sites counts many leaves attached to hubs as though each could supply a fixed-rate first promotion. The correct total count of realized starts up to the horizon is `C_n+O(1)`.
* Ignoring killing and integrating sequential promotions over a time window of length comparable to the path length discards the factorial upper bound. The competing internal sends supply a length-dependent rate; killing at receipt, not at every ring, is the right exact count.
* A one-fringe coupling error of order `1/d` is far larger than the success probability at the target scale unless the parent degree is enormous; subtracting that error does not establish the lower bound. Using no independently allocated incoming proposals is a multiplicative factor bounded below by a constant.
* An observation-time upper bound alone does not interpolate by path or component monotonicity. The upper bounds above cover entire size horizons or physical-time blocks, and lower successes persist over a unit interval.
* No centred extreme-value law, sharp `N>=3` logarithmic prefactor, or uniformity when `N=N(n)` grows is proved.
