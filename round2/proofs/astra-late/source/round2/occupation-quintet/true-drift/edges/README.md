> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# Matched raw high-core and joined-edge scalar

At the same exact trial `s=.565724`, the [negative high-core block](../moments/README.md) becomes **positive by at least `9.1908e−4` per core** when the original old-edge residual is included on every edge with a `D≥2000` endpoint. This includes both physical directions of core gain on those edges and gains from exterior donors. Adding high sends into virgin low singletons gives a broader scalar in `[9.1908e−4, 3.00631e−2]`. The remaining sector is negative, in `[-3.01349e−2, -8.3879e−4]`. These are signed integrals in potential units, not shifts in the exponent. Their cancellation is forced to within `8.8%` of the positive scalar by the already-certified bracket; the sign of the full raw score is unresolved.

## Support and exact balance

Use `F=f_lower+f_upper`, `g=g_lower+g_upper`, the original `2^40` scaling, and the same pre-event cutoff `K=106` as the parent work. Let `U` be the sum of old `g` weights of edges with at least one `D≥2000` endpoint and both endpoint arities `≤K`. These are unoriented edge copies. Exact reads of the endpoint face show that **all paired `g` words for two interpolating endpoints are zero**, so an edge can contribute only with its other endpoint below `d_g=180`. This is a fact about the shipped pair, not an assumption about how often two high cores are adjacent. It also implies `U≤0` using the originally verified far-edge range `[-CAP,0]` for the two weights, rounded conservatively to `CAP=350000000000/2^40` for their sum.

Define `E` as the true old-edge discount and every old selected edge's gain and promotion jump on this **pre-event** support. Discount and old promotion are counted once; every gain uses the true receiver, and retains a core donor only when its arity is `≤K`. Thus for a mixed edge `u–v` its two internal arrows pay `(1/D_u) Δ_D^v g` and `(1/D_v) Δ_D^u g`; a donor outside the edge is appended once for each occurrence affected. `E` evaluates the actual successor, including arity `K+1`. As usual, an endpoint gain or promotion changes no other old mark. The distribution of these arrows is not replaced by product marginals.

Write `B` for the high promotion's newborn edge, already in `H`, with its full post-promotion mark; write `I` for the entry of an old low-low edge at a gain `1999→2000`. **`I` is not a raw jump row:** it is the new selected weight in the change of the indicator for `U`. The raw change at that event is in the low pre-event sector. A low promotion never creates a selected edge; a high promotion never creates a larger supported motif. Compensation on `U` gives the exact bookkeeping identity

```
π E = (μ−s) π U − π B − π I + ε_K,
π(H+E) = π(H−B) + (μ−s) π U − π I + ε_K.
```

Here `I` uses full incoming at the seam, and `ε_K` is the signed difference from the arity boundaries: omitted donors on selected gains, the old potential removed at an endpoint promotion `K→K+1`, and newborn edges counted in `B` at old `c=K` but excluded from capped `U`. The bound on its absolute integral is `CAP[K π(c>K)+(K+1)π(c≥K)]`. The first term uses the overflow donor's total outgoing rate to cores `c/D≤1` and at most `K` affected edges at a surviving receiver. The other two use the promoting endpoint's rate `≤1`. The subtracted seams are evaluated at their **new** edge value, including the actual low neighbor. There is no dropped old-neighbor term in `E`: preexisting old edges at the same event can contribute both through ordinary drift and through entry, depending on which edges were already selected.

The balance follows by applying the physical mean-law compensation to the signed capped edge potential, using fixed finite-census approximations if needed. All ordinary jumps of selected edge rows are integrable by the [parent raw-drift proof](../README.md); the additional single-degree seam is bounded by the finite flux below. The capped potential has a per-core incidence bound. The old-edge identity uses the actual common `μ` throughout; no stationary-zero substitution at trial `s` is made.

## Quantitative inputs

The implemented rational two-state-coordinate solve for `π(D≥2000)` uses the first hit of the seam from physical birth `(2,1)`. At each `D<2000,c≤34`, bound the Laplace hitting transform with both physical gain corners `D−c,D−c/2`, the true promotion rate, and both discount corners enclosing `μ`. Use value 1 at the degree hit, and values 0 and 1 at arity exit for the lower and upper solve. Since the stable-age observation is exponential, these hitting-transform values bound the high mass with **no extra factor of μ**. To bound high promotion intensity deduct at most `π[c;D≥2000]/2000` from the lower bound. The common-age proposal moments from the preceding block give the high arity moment bounds `c₁,c₂` for `π[c;D≥2000], π[c²;D≥2000]`.

* **Birth routing:** sparse knot/INF interpolation of `g(post,J)` for *all* clamped arities and `D≥2000` gives a uniformly negative coefficient. Multiplying it by the signed-safe promotion intensity interval, adjusting for `c>K`, bounds `π B`. A separate virgin-residence floor for `−π U` counts each retained birth once, and pays its signed weight until the child first changes or the parent exceeds `K`. The uniform kill rate is at most `7/2`, and `−π U≤CAP·c₁` by incidence reversal.
* **Degree entry:** each core has exactly one first seam crossing. Its numbers of low neighbors at that crossing are bounded by its core arity there. By monotonicity this discounted crossing arity is at most `c₁`: couple the actual root with its independent exponential observation and count the core arity at observation on the event of an already completed crossing. Normalize crossing flux by multiplying by `μ`. This gives `π I≥−CAP·HI·c₁`. For a genuine negative floor, start from a first entry to `1999` with `c≤34` and require a promotion before the next root gain, then a root gain before the virgin changes. The first event has discounted probability at least `(1−34/1999)/(HI+2000)` using the rounded upper discount in code. In the second race, the root gain is at least `1998−34`, and competing hazards are at most `7/2` even if another root promotion is treated as failure. The crossing is then on an actual virgin edge with `g≤−.01`. A lower bound on the entry mass is the high-core mass floor minus `P(c≥34) upper`; this also safely covers the arity requirement. Nothing in this construction prices that cylinder twice as a high promotion.
* **Joined directions:** the replay checks the high-end interpolation slope, post-promotion slope change, and reverse gain jump on the high-sender/low-receiver face, allowing the low endpoint to cross its interpolation seam. Beyond degree 180 on both endpoints the paired jump is zero. For all high-receiver gains (including leaf, internal and exterior core donors), use `G≤D`; for the core-only low-sender part use the actual sum bounded by `c/2`. The reverse high-to-low arrow is weighted by `1/D_high`. The reported directional intervals overlap their containing subtotal, so they are diagnostics of the same certified decomposition, not quantities to add again.
* **Remaining sender split:** on a virgin low singleton `(2,1)`, a high donor contributes exactly `(F(3,1)−F(2,1))/D_high>0`; there is no old low-low edge through this receiver. Add this slice to `H+E`. High sends into other low singletons and appended to low-low old edges are kept in the remaining sector. An independent absolute bound uses `64` per raw singleton/edge jump, `K` affected edges per surviving receiver and the high-donor reciprocal moment `c₁/2000`. That standalone bound is very loose (`.535` per core); the signed aggregate comparison below avoids treating such a sender bound as an independent additive error.

The same-support bound for `H−B` intersects the prior certified `H` interval minus `B` with a direct expansion that removes `g` from the birth polynomial. In the latter the leading constant is `F(J)−s a`, with linear coefficient `−k+(1+s)h+c[(1+s)k−F(J)]` in `1/D`, in the notation of the [high-core proof](../moments/README.md). Enlarge the linear coefficient of the quadratic envelope by `max(0, k−max|z_c−k|)` to cover the shortened promotion polynomial. The same common-age degree/arity and overflow norms apply. No charged `f/g/FY` reward is used.

| Physical row, per core | Outward bounds |
|---|---:|
| High birth edge `π B`, already in `H` | `[-1.894151e−3, -9.208628e−4]` |
| Old edge seam-entry flux `π I` | `[-2.812615e−2, -1.864199e−8]` |
| Selected old-edge subtotal `π E` | `[9.196680e−4, 3.002139e−2]` |
| Gains at high endpoint, all senders | `[-1.033143e−4, 1.711118e−10]` |
| Promotions at high endpoint | `[-3.077580e−8, 2.427473e−5]` |
| Internal high-to-low arrow | `[-8.147936e−11, 1.614107e−5]` |
| Low core sends into high receiver, internal plus exterior | `[-2.887924e−7, 4.783053e−13]` |
| Gains at low endpoint, all senders | `[-8.473853e−6, 3.228214e−2]` |
| Promotions at low endpoint | `[-1.400723e−2, 1.174979e−4]` |
| Added high send into virgin low singleton | `[0, 4.181342e−5]` |

## Comparing with the remaining low-action sector

Let `A=H+E+virgin sends`. The exact full raw score is `π A + π C`, where `C` includes **every other** old or newborn raw summand with its original sign: the low-action singleton block, old low-low edges including the actual `1999→2000` raw jump, low promotions and their newborns, the remaining reciprocal sender corrections, higher motifs (including high tips and exterior donors), and uncapped arity rows not in `A`. Thus `C` is the residual after this low/high split, not a claim about a strictly all-low neighborhood. The whole raw identity is `(μ−s)π(Φ_lower+Φ_upper)`. The inherited all-tree comparison and original gold bracket imply `[-7.178818e−5,+8.029322e−5]` for this full score. Subtract the matched positive interval to bound the **signed residual**. In particular `|π A+π C|/π A < .088` on the same physical law. This is a quantitative cancellation using the pre-existing bracket, not an independent pricing of `C` or a new endpoint. The dominant width in the positive scalar is the seam-occupation envelope. A substantially tighter interval needs joint low-neighbor occupations at the first seam crossing, rather than merely `c₁` there.

835c independently measured the **charged** original family spectra. Its held-out priority region for H3 (middle `c=2–4`, max tip `D=19–1998`, other tip unrestricted, outside prior strict masks and founder sets) measured `5.727/7.009e−6` per core at the original lower/upper trials. That is a diagnostic pointer for further work, not a raw low-sector bound or additive evidence here.

## Serial replay

```sh
taskset -c 5 nice -n 19 python3 -IBu round2/occupation-quintet/true-drift/edges/check.py
```

The replay checks the original detached source pin, every cross-multiplied Bellman corner, and five direct ordinary-tree fixtures (including the exact degree seam, high-high weight zero and two distinct high donors). It reads 125,888 sparse words per original payload, chiefly the high `g` face, with no whole-table scan. The checks are entirely rational before displaying decimals.
