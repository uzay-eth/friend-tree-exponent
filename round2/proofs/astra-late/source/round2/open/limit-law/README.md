> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude and Gemini reviews; Claim 6, inert limits, and attachment influence; the cusp also uses repaired Fable Claim 6. See the [release index](../../../../README.md#claim-17) for exact scope and citations.

# Law of the Claim 6 limit: uniform jumps, tails and support

Let `X=lim C(T_n)/n^μ` from the edge (equivalently, from `T₃`). The following strengthens the L² and support-endpoint conclusions of [Claim 6](../../../../../astra/source/round2/standalone/claim-6.md): the topological support of `X` is `[0,∞)` and there are positive constants `c,C,x₀` such that

```
exp(−C x^(1/(1−μ))) ≤ P(X≥x) ≤ exp(−c x^(1/(1−μ)))       (x≥x₀).
```

In particular the tail has a sharp double-log exponent, `log(−log P(X≥x))/log x → 1/(1−μ)`, and all positive moments exist. This does not identify an asymptotic constant in the exponent. The proof concerns the embedded jump-chain limit, not the continuous-time eigenmartingale divided by a random size on a concentration event.

## 1. Uniform gradient in a free tree

Normalize `Φ(T)=lim_(t→∞) exp(−μt) E_T C(T_t)=d₂ Ψ(T)` in the notation of Claim 6. The limit is independent of orientation. Every finite free seed with a core is covered by the physical limit; the edge is covered by conditioning on its first birth (rate two, yielding `T₃`), which gives `Φ(edge)=2Φ(T₃)/(2+μ)`. The same eigenfunction and comparison estimates extend to the edge with the additive constant in `C≤K(1+Φ)`.

For an arbitrary one-leaf extension, root the **ordinary** common tree at the attachment vertex, without adding any stub. Couple each ordered common edge at the minimum of the two rates and assign the excess to its true marginal. The roots' common degrees are real degrees on both sides. Pair simultaneous new leaves, color the common graph by parity from the attachment, and keep private branches attached to only one common vertex. At a common site all private neighbours are on its parity-favoured side: a surplus send selects the lower-degree source copy, which is the higher-degree target copy. All common vertices, including degree-one leaves, are paired. These are the free marginals of the [first-infection coupling](../../../../../../reviews/cold/astra/influence/check/README.md).

In particular, if the attachment is to a genuine leaf, the discrepancy is **at that same paired vertex**: its lower degree is one and its upper degree two. No extra discrepancy is seeded at its old parent. Its core indicator differs from time zero and is counted as one retained singleton. Initially infected vertices need only the total surplus intensity bound one, valid even when their smaller degree is one. At all later infections the single arriving private leaf makes the initial degree gap exactly one; the only already infected common neighbour is the incoming one (connected infected set in a tree). The conditional outgoing infection bound over its previously existing real neighbours thus applies unchanged without a slot. The degree-one case has no eligible outgoing edges. Private branches cannot reconnect elsewhere; returns from them or from already infected common neighbours are included in the conditional rate estimate.

For clarity the planted gradient's (`round2/influence/gradient/README.md`, historical reference) partition works without any core parent assumption. At the first infection of a retained anchor `x`, leave the already-common outward edges in the infection recursion. Charge every later ordinary branch attached directly at `x`, including the incoming private leaf itself, to a disjoint account on its own side. Attachments targeting vertices in an account stay there forever. A clean unowned vertex that first gains a private leaf must receive it from a retained anchor along an already-existing common edge, since branches have only their edge back to an infected anchor. The coupled core-count difference therefore consists only of the retained singletons and the signed branch accounts. Already-present structural leaves remain in the recursion, even if they promote on only one side.

A branch account starts as one leaf. Until its first promotion it costs no core; that promotion has to be caused by its anchor. Any single anchor spends at most rate one on these first branch promotions in each marginal, with no bound needed on how fast ordinary leaves are appended there. At the promotion the branch has one core, and its later core population is exactly the spatial cell of that core when labels are reset at that time. Its attachment vertex is a core on that side (it has the edge into the branch and at least one other neighbour). The full exterior stays active. The [free spatial estimate](../../../../../astra/source/round2/fringe/spatial-limit.md) gives its conditional expected count `≤C_Q exp(μu)` after elapsed time `u`; no comparison with a planted or silent parent is made. Nonexplosion and strong Markov hold in the paired filtration, whose conditional marginal rates equal the individual physical rates.

Take `s_b=.525`, `κ=.97`, `b=1/s_b`. The initial seed charges at most `b` discounted outward infections; each noninitial retained vertex charges at most `κ` by the conditional first-infection bound. Thus `E Σ_x exp(−s_b τ_x)≤D_b=1+b/(1−κ)`. At time `t` the discounted expected bound on the retained singleton differences is at most `D_b exp(−(μ−s_b)t)`. The two branch-promotion intensities bound their normalized contribution by `(2C_Q/μ) E Σ_x exp(−μτ_x)`. In particular

```
|Φ(T⁺)−Φ(T)| ≤ (2 C_Q/μ) [1+ (1/μ)/(1−.84)] =: L.
```

For the sharper constant we used the exact rational [.565 certificate](../../../../../../reviews/cold/astra/influence/check/README.md) at `μ≥.565` for branch contributions; the `.525` discount still makes the singleton term vanish strictly, including if `μ=.565`. Division by `d₂>0` gives the requested uniform bound for `Ψ`. The argument also gives a uniform finite-time bound `|E_T C_t−E_(T⁺) C_t|≤e^(μt) D_b(1+2C_Q/μ)`; it is not inferred by interchanging a supremum and a limit. No gradient is assumed for Long's full-environment genealogy. Long's independent inert and free eigenfunction constructions (`round2/long-comparison/README.md`, historical reference) provide related core comparability, not this attachment bound.

## 2. Conditional concentration and upper tail

Set `γ_n=Γ(n)/Γ(n+μ)` and `Y_n=γ_n Φ(T_n)`. The conditional drift from the continuous-time eigenfunction and the total rate `n` is `E[Φ(T_(n+1))|F_n]=(1+μ/n)Φ(T_n)`. The gamma recurrence therefore centers the jumps exactly:

```
Y_(n+1)−Y_n = γ_(n+1)[Φ(T_(n+1))−Φ(T_n)−(μ/n)Φ(T_n)].
```

Because the centered term is a one-step change minus its conditional mean, its absolute value is at most `2L` without any size estimate on `Φ`. Thus `|Y_(n+1)−Y_n|≤A n^(−μ)` for a finite deterministic `A`. Claim 6 identifies its limit with `X`; alternatively square-summability gives L² convergence of this already identified martingale. For every reachable history and every `t>0`, conditional exponential supermartingales (or conditional Azuma–Hoeffding) followed by the limit give

```
P(|X−Y_n|≥t | F_n) ≤ 2 exp(−t²/[2 A² Σ_(j≥n) j^(−2μ)]),
Var(X|F_n) ≤ A² Σ_(j≥n) j^(−2μ) ≤ V n^(1−2μ).
```

One-sided tails omit the factor two. This yields unconditional sub-Gaussian concentration by starting at `n=3`, and uniformly vanishing conditional variance over *all* reachable states. For the sharper right tail use `Φ≤K C≤K n`, so `Y_n≤B n^(1−μ)` deterministically. Choose `n=floor((x/(2B))^(1/(1−μ)))` for large `x`. Then `Y_n≤x/2` and `n` is bounded below by a constant times `x^(1/(1−μ))`. The conditional upper bound is `exp(−c x² n^(2μ−1))≤exp(−c' x^(1/(1−μ)))`. In particular `E exp(ε X^(1/(1−μ)))<∞` for sufficiently small positive `ε`.

## 3. Matching lower right-tail rate: keep degree at most three

The path history from Claim 6 costs order `exp(−n log n)` and misses the desired rate. Instead use a history that keeps **maximum degree at most three**, allowing either leaf promotion or attachment to a degree-two vertex. For such a free tree of size `j≥3`, with `L` degree-one vertices and `N₂` degree-two vertices, the tree degree sum gives `N₃=L−2`, `N₂=j−2L+2`. If `L≥j/4`, the probability that the next event targets some leaf is at least `L/(3j)≥1/12`, since its source has degree at most three. Otherwise `N₂>j/2`; the total probability of targeting a degree-two vertex is at least `2N₂/(3j)>1/3` by adding its two directed incoming channels. These events preserve the cap. The actual success probability may be larger; even the adaptive rule that uses just the indicated class survives through size `n` with probability at least `12^(−(n−3))`.

Every degree-three-capped tree satisfies `C≥(n−2)/2`: this follows from `L=N₃+2≤C+2`. By the reverse `Φ` comparison, uniformly on those histories `Y_n≥b₀ n^(1−μ)` for a fixed `b₀>0` and all sufficiently large `n`. Choose `n=ceil((2x/b₀)^(1/(1−μ)))`. Then `Y_n≥2x`, and the full conditional martingale variance bound gives `P(X≥x|F_n)≥1/2` uniformly there for large `x` (its exceptional probability is at most `V n^(1−2μ)/x²`). Multiplying proves the matching lower bound. It also implies `E exp(ε X^(1/(1−μ)))=∞` for sufficiently large `ε`; neither bound asserts the threshold value.

## 4. No gaps in support; a near-zero lower estimate

For any `x>0` and open neighbourhood of `x`, choose `N` large and prescribe the star through size `N`. This event has positive probability; `Y_N≤O(N^(−μ))` because the star has one core. From this star, elongate a single arm by repeated attachment to its endpoint. Each choice has positive probability because the endpoint has a neighbour of finite degree. Along this particular infinite *sequence of shapes*, not an event asserted to have positive probability forever, `C(T_j)=1+j−N` and hence `Y_j→∞` by the reverse comparison. Between consecutive states `|Y_(j+1)−Y_j|≤A j^(−μ)`. So at a first crossing of `x`, `Y_j` differs from `x` by at most `A N^(−μ)`. The whole finite prefix has positive probability. Taking `N` still larger makes the uniform conditional variance `V N^(1−2μ)` small enough that the chosen neighbourhood has positive conditional probability for `X`. This proves every `x>0` lies in the support (and the closure contains zero); positivity almost surely is the original Claim 6 conclusion. No monotonicity of the eigenfunction along the prescribed prefix is needed. Full support does not imply a density or absence of atoms.

One can also quantify a lower bound near zero. On the size-`n` star event, of probability `2/(n−1)`, its conditional mean `Y_n` is at most a constant times `n^(−μ)`. Choose `n` of order `ε^(−1/μ)` with `Y_n≤ε/2`; Markov's inequality gives `P(X≤ε)≥c₀ ε^(1/μ)` for small `ε>0`. No matching near-zero upper rate is established by the uniform variance bound.

## Short controls

```
taskset -c 8 timeout 8s python3 -B round2/open/limit-law/check.py
taskset -c 8 timeout 8s python3 -B round2/influence/check/bound.py --s .565
```

The first command checked 268 exact paired free states from an initially genuine leaf, 1,672 degree-capped labelled histories, and 18 delayed-star arm steps. The second gave the rationally outward displayed first-infection budget `.839292428 < .84`. These replays test the physical coupling and reachable scripts; the tail and support conclusions use the analytic martingale estimates.
