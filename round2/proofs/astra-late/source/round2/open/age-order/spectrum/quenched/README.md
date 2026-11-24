> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT, Claude, and Gemini within their recorded scopes; Claims 6, 10, 11, and 18. See the [release index](../../../../../../README.md#claim-22) for exact scope and citations.

# The empirical birth-label and mass process

In the original edge-seeded friend tree let `Z_v` be the limiting degree fraction of label `v`, and `X=lim C_n/n^mu`. Use the finite **intensity** measure `ν` on `(0,∞)` from the [annealed spectrum](../README.md): `ν((z,∞))=π(F_z)`, where `F_z` tests continuous mass measured **from the instant of birth**. Then, on one probability-one event, for all `0≤a<b≤∞` and `r>0`,

```math
x^\mu\sum_v\mathbf1\{a<xv\le b,\ Z_v>rx\}
\longrightarrow X\int_a^b y^{\mu-1}\nu((ry,\infty))\,dy. \qquad\text{(1)}
```

The [independent cold audit](../../../../../../../../reviews/cold/astra/checks/marked-spectrum-sol.md) confirmed the result under the inputs cited below. The almost-sure assertion is **simultaneous** over all these rectangles, including label endpoints zero and infinity. The `L¹` assertion holds for **each fixed deterministic rectangle** with `r>0` as `x↓0` through arbitrary deterministic values, and for each fixed deterministic continuous compactly supported test on `(0,∞)^2`. It does not assert an `L¹` bound for the supremum over all rectangles or uniform control as `r↓0`.

Thus the empirical measures of `(xv,Z_v/x)` converge vaguely on `(0,∞)^2` to `X` times the pushforward of `y^{mu−1}dy ν(du)` under `(y,u)↦(y,u/y)`. Conditioning to `Z_v>x` and choosing uniformly among those hubs **in the realized tree** gives the joint law in the annealed note, without count-biasing the tree. In particular the asymptotic independence of the power-biased `vZ_v` and the Pareto `Z_v/x` is a statement about that sampling law, not independence of the marked vertices' futures.

## 1. Actual clock and local cohort mean

Start the continuous clock at the three-vertex tree. Write `B_v` for the actual birth time of subsequent vertex `v`, `Y=lim e^{-s}|T_s|>0`, and `V_v=YZ_v` for the mass in this clock. Then `v e^{-B_v}→Y` a.s., and `e^{-mu s}C(T_s)→H:=Y^mu X`. For real `l<d` and `r>0` first prove the time-cohort limit

```math
e^{-\mu t}\#\{v:t+l<B_v\le t+d,\ V_v>r e^{-t}\}
\longrightarrow H\int_l^d e^{\mu z}\nu((r e^z,\infty))\,dz . \qquad\text{(2)}
```

Let `F_q^u` be the per-core birth intensity score for the newly placed tip satisfying `e^{-u}D(u)>q` at age `u`. Measure both the degree and the elapsed time from the actual insertion; the old leaf targeted by a birth is never counted as the tip of that birth. For `q≥q0>0` and all sufficiently large `u` (so `q0 e^u≥2`), `0≤F_q^u≤7`: the event requires promotion, and the all-time promotion intensity bound applies to both the core source and the old leaf sources. Likewise `0≤F_q≤7` for the infinite-age mass. The finite-horizon history with marks after a given placement uses the old-cell coupling; the continuity and atomlessness for the infinite-age mass are proved in the [scaled note](../../scaled/README.md).

In fact if inputs approach any fixed locally finite marked core, `q_i→q>0` and `u_i→∞`, then `F_{q_i}^{u_i}` on those inputs tends to `F_q` on the limit. For a fixed placement, the uniform `L¹` terminal error from the scaled note is at most `e^{-u/2}+(c(0)+u+2)e^{-u}` for a degree-one tip, including any initial-neighbour correction. The same bound works on any locally convergent input sequence (the initial degree and relevant initial neighbour category are eventually fixed). Finite-horizon law continuity, this uniform approximation, and the absence of a positive mass atom give the claim for each placement. There are only finitely many placements, counting repeated old-leaf sources with their eventually fixed multiplicity, at the limiting input. By Claim 10 local empirical weak convergence and bounded sequential continuity (apply Skorokhod representation on the rooted marked-tree space), for every `z∈[l,d]` on the local-LLN event,

```math
\frac{1}{C(T_{t+z})}\sum_{o\,\mathrm{core}} F_{r e^z}^{T-t-z}(T_{t+z},o)\longrightarrow \nu((r e^z,\infty)), \qquad\text{(3)}
```

whenever `T−t→∞`. Here the convergence of the empirical law holds at *all* real times by the nonexplosive jump chain and the discrete local LLN. Values at jump times in the following integral may equally be taken just before the jump.

Let `A_t` count tips born in `(t+l,t+d]` with `D_v(T)>r e^{T-t}`, for a deterministic observation time `T>t+d`. At an instant `s` of this window, the intensity of a birth times the conditional probability of its terminal success, integrated over the possible initiating arrows, is

```math
g_t(s)=\sum_{o\,\mathrm{core}} F_{r e^{s-t}}^{T-s}(T_{s-},o).
```

The threshold is correct because `e^{-(T-s)}D(T)>r e^{s-t}` is equivalent to the common terminal degree cutoff. `g_t` is predictable: it averages over the prescribed insertion and its unconditioned physical continuation. But it is **not** asserted to be a compensator of an adapted process revealing eventual success at birth. The next step supplies the error between `A_t` and its time integral. From (3), the a.s. boundedness of `e^{-mu s}C(T_s)` at large `s`, and dominated convergence on the fixed window,

```math
e^{-\mu t}\int_{t+l}^{t+d}g_t(s)\,ds
\longrightarrow H\int_l^d e^{\mu z}\nu((r e^z,\infty))\,dz. \qquad\text{(4)}
```

## 2. Cohort fluctuation, without independent terminal marks

Set `s0=t+l`. Use the backward martingale `M_s=E[A_t|F_s]`, `s0≤s≤T`, in the actual birth-time-marked process. Its intensity of jumps is the ordinary size. A one-leaf perturbation, keeping all existing birth times fixed and giving new matched vertices the common time of the matched jump, changes the conditional mean of the terminal cohort count by at most `3 e^{2(T-s)}`. To see this, use the dirty-core coupling of [the core-count proof](../../../../../../../astra/source/round2/fringe/almost-sure.md), whose clean paired cores retain identical degrees, birth marks and promotion status. All terminal degrees counted are at least two. Inserting at a core initially dirties that core in both trees; inserting at a leaf additionally creates the promoted dirty core and dirties its old parent in both copies. Clean-core rings synchronize promotions of paired leaves, and paired leaf clocks add matched leaves; a dirty-core ring creates at most one dirty core or dirties one pair. Dirty leaf initiations add private leaves at a dirty parent without creating a core. Thus the expected combined dirty-core count is bounded by `3 e^{2u}`; discard its contributions and all remaining eligible endpoints cancel. No control of the number of dirty leaves is required. The insertion changes no old birth mark, and the new leaf and any descendants are counted on their true sides. The same coupling covers the future of not-yet-born vertices. This is the bound for a jump of `M_s`, including the new vertex if eligible.

Write `I_t=∫_{s0}^{t+d}g_t(s)ds` and let `P_s=E[I_t|F_s]`. At a fixed `s` the accumulated past integral has no jump. To bound the jump of the remaining expectation, use immediate-birth compensation and the tower property on the remaining part of the window: its conditional expectation equals the conditional mean of the terminal cohort count restricted to births strictly after `s`. Again the same coupling gives `3 e^{2(T-s)}` (do not include the inserted vertex in this remaining cohort). Finally `E[A_t|F_s0]=E[I_t|F_s0]` by the insertion intensity calculation and tower; there is no deterministic-time birth at `s0`. Use the physical jump filtration, including birth times but not future holding times or drivers. The first-jump equation gives the compensated-jump representations (between jumps the forecasts evolve by the deterministic time variable and the continuous accumulated integral). Square integrability follows from finite-horizon Yule moments, and stopped jump isometry removes finite-size cutoffs. Hence

```math
E\big[(e^{-\mu t}(A_t-I_t))^2\big]
\le 36 e^{-2\mu t}\int_{s0}^{T} e^{4(T-s)} E|T_s|\,ds
\le K_l\exp\{-(2\mu-1)t+4(T-t)\}. \qquad\text{(5)}
```

One may obtain (5) directly by subtracting the martingales, as their jumps differ by at most `6 e^{2(T-s)}`; the looser separate isometries give the same constant. This does not involve any covariance assumption about two hubs.

## 3. Terminal degrees track physical limiting masses

Take `T=t+2 log(t+2)`, so `k=e^{T-t}=(t+2)^2`. The [hub-spectrum estimates](../../../../../../../astra/source/round2/fringe/hub-spectrum.md) give maximum core degree `O(T)` a.s. and the conditional gamma comparisons for `e^T V_v` given the history at `T`. For fixed rational `0<ε<1/4` and `r>0`, conditional on `max_v c_v(T)≤ε r k/4`, the probability of *any* existing vertex violating either implication

```
D_v(T) > (1+ε) r k  =>  V_v > r e^(-t),
V_v > r e^(-t)       =>  D_v(T) > (1−ε) r k
```

is at most `K E|T_T| exp(−c ε² r k)` after taking expectations with the indicated history gate (adjust constants, or slightly shrink the internal comparison margin, for strict endpoints/integer degrees). The first implication uses `Gamma(L,1)−Exp(1)` in the local comparison, the second the degree-based upper gamma comparison. All vertices in the birth cohort are present by `T`. This bound is summable on each fixed time mesh; the history gate holds eventually a.s. there, without needing a probability bound on the gate to get that conclusion. Equivalently one may use the independent two-colour fractional Pólya upper comparison after promotion for deterministic-size upper tails; its auxiliary fraction is not a common mark for different vertices and does not give a lower concentration estimate. For the simultaneous sandwich we use both local gamma comparisons in the common physical clock.

Equation (5) is now `O_l((t+2)^8 e^{-(2mu−1)t})`. This is a quantitative centered second-moment estimate. For example the centered normalized error is a.s. `O(e^{-γt})` along any fixed mesh for `0<γ<(2mu−1)/2`, by Chebyshev and Borel–Cantelli. It is not a rate for (1) or (4): (3) uses the qualitative shifted local LLN, and the gamma sandwich permits fixed then vanishing relative margins. Intersect the good events for rational endpoints, thresholds and margins, and all meshes `j/m` with positive integer `m`. At mesh points use (4) at thresholds `(1±ε)r` and send `ε↓0`. The limiting integral is continuous for finite windows and `r>0` (positive atomlessness of `ν`, or integration in `z` and dominated convergence). This proves (2) on the meshes.

For arbitrary `t` between neighbouring mesh points `s` and `s+1/m`, there is no monotonicity of a translated window. An inclusion instead bounds its count below by that at `s` with window `(l+1/m,d]` and threshold `r`, and above by that at `s` with window `(l,d+1/m]` and threshold `r e^{-1/m}`. Use rational bounds on the latter threshold if necessary, and the factor `e^{-mu(t-s)}` between the normalizations. Sending `m→∞` gives (2) for all real `t` with rational parameters, then widen/narrow with rationals to obtain all real finite parameters on the same event.

## 4. Return to labels, extend the windows, and compare sampling laws

For `x=e^{-t}`, a label in a compact interval of `xv` bounded away from zero has `B_v-t=log(xv/Y)+o(1)` uniformly along that interval, pathwise, by `v e^{-B_v}→Y`. Conversely the same relation is uniform on the corresponding compact time windows. The cutoff for `V_v` is exactly `rY e^{-t}`. Apply (2) with the sample-dependent endpoints `log(a/Y),log(b/Y)` and threshold `rY`, permitted by its simultaneous pathwise form (or sandwich with rational parameters). Changing variables `y=Ye^z` cancels the clock factor since `H Y^{-mu}=X`. This gives (1) for `0<a<b<∞`. It uses neither independence from `Y` nor substitution of a deterministic birth time for an individual mark inside the conditional score.

For zero or infinite endpoints use the [unmarked almost-sure spectrum](../../../../../../../astra/source/round2/fringe/hub-spectrum.md), which gives the total count above `rx` with normalization tending to `KX r^{-mu}`. By the [Mellin identity](../README.md), this equals the right side of (1) on the full label axis. Thus the a.s. limsup of the normalized discarded count outside `(a,b]` is the matching integral outside that interval (indeed the difference of total and compact counts has that limit). Exhaust the axis and squeeze. The initial labels handle any strict/inclusive convention at zero with normalized error `O(x^mu)`. The annealed note also proves the total-count limit in `L¹`; hence its normalized counts are uniformly integrable as `x↓0`. Every rectangle count at fixed `r>0` is bounded by that total count for threshold `rx`, times the fixed normalization factor. This gives `L¹` in (1), with no new moment assumption. Quantitatively the annealed degree-two Pólya bound at the actual promotion size gives the expected normalized omitted-label bounds `O(a^mu)+O(x^mu)` at the left and `O(∫_b^∞ y^{mu−1}(1+ry)e^{-ry}dy)` at the right in the limiting upper envelope, for fixed `r>0` and large `b`. The auxiliary beta tail is conditioned on the full history at promotion before averaging; sharing it across marks would be unjustified.

Specifically, at the promotion size `m` generate subsequent fresh uniforms independent of that past. At size `j`, use the same uniform to test the physical target probability for this vertex and the auxiliary probability `P_j/j`, using a separate conditional target draw for the complementary physical outcomes. The auxiliary chain starts at `P_m=2` and dominates the physical degree by induction. Its limiting fraction has the `Beta(2,m−2)` law conditionally on the promotion history. This is a one-tag coupling; the beta variable need not be independent of the future of the physical tag coupled to it. Conditioning additionally on mass success or taking products of beta tails for several vertices is not part of the argument.

To read out the sampling law, divide (1) by the total normalized count above `x`, whose limit is `KX>0`. This gives weak convergence of the realized conditional probability measures (rectangle convergence, tightness from the limiting probability and exhaustion). In variables `(Y_label,U)=(xv,vZ_v)` the law is `K^{-1} 1_{0<y<u} y^{mu−1}dy ν(du)`; its factorization is the same as in the parent note. In contrast, the annealed count-weighted **tree** law carries the limiting density `X/E X`; uniform in-tree sampling requires no change in the tree law. The random empirical intensity retains the random factor `X` before normalization.

## Finite check

`python3 -B round2/open/age-order/spectrum/quenched/check.py` passed 108 exact conditional forecast rows and caught 102 controls that incorrectly included the already-inserted tip in the remaining forecast. The check uses the exact physical jump chain with a bounded birth-label window and finite terminal size; the moving-time LLN and clock transfer are analytic.

The [cold audit's independent fixture](../../../../checks/marked-spectrum-sol.py) checked the backward-forecast jump variance in three physical finite-size windows, found nonzero cross terms under separate birth centering, and detected the wrong drift from re-including the inserted tip. It also checked 150 exact Beta moment/tail rows and 144 strict translated-window inclusions. These are bookkeeping and dependence checks; they do not give a rate for the full almost-sure or `L¹` limit.
