> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and scoped Claude consistency checks; earlier spectral constants and, for the named early-label splice, Fable transports. See the [release index](../../../../../README.md#claim-24) for exact scope and citations.

# Fixed thresholds: polynomial bounds and the adjacent-ratio gap

For the original edge seed, put `P_v(x)=Pr(Z_v>x)`. For every fixed `0<x<1` the bounds below give

```math
 c_x v^{-3}(1-x)^v\ \le P_v(x)\le C_x v^{\mu-1/2}(1-x)^v. \qquad\text{(1)}
```

These improve the exponential-rate statement to a logarithmic error in `log P_v(x)` but do not identify a sharp prefactor or determine the sign between adjacent labels. The bounds are on the actual deterministic birth label. No assumption that the random birth parent has a typical degree is used.

## 1. The upper bound: discount promotions by jump count

At the birth of label `v`, conditional on the size-`v−1` tree, the target `u` is selected with probability `r_u/(v−1)`, where `r_u=Σ_(z~u)1/D_z`. Let `M_v` be the size at first promotion, with `∞` if there is none. Given the birth target, the permanent parent's degree is at least `D_u+1`. Before promotion the discrete hazard at size `j≥v` is at most `1/[j(D_u+1)]`. Averaging over the actual target law gives, for `m≥v+1`,

```math
 \Pr(M_v=m)\le\frac{1}{(m-1)(v-1)}\mathbb E\sum_u\frac{r_u}{D_u+1}
 \le\frac{3\mathbb E C_{v-1}}{2(m-1)(v-1)}. \qquad\text{(2)}
```

For the last bound sum core terms using `r_u≤D_u`, and sum old leaf terms by charging `1/(2D_o)` to their core parents. This is an upper bound on a **single** promotion size, distinct from the eventual promotion bound.

There is an urn improvement after promotion. The marked vertex has a permanent core parent, so its conditional attachment rate is at most `D−1/2` for the rest of its life (all other neighbours contribute at most one). At first promotion the marked degree is two at size `m`. Compare `D−1/2` with a Pólya urn with red mass `3/2`, total mass `m`; for equality of red masses couple the increments using the common uniform, and otherwise the unit increments cannot cross. The limiting cap is `Beta(3/2,m−3/2)`; subtracting the constant `1/2` does not alter the degree fraction. For `b=m−3/2`, integrate the beta density using `t^(1/2)≤1` on `[x,1]` to obtain

```math
 \Pr(\mathop{\mathrm{Beta}}\nolimits (3/2,b)>x)
 \le \frac{\Gamma(m)}{\Gamma(3/2)\Gamma(b+1)}(1-x)^b
 \le C\sqrt m\,(1-x)^{m-3/2}.
```

Strong Markov at `M_v`, (2), and `E C_n=O(n^μ)` bound the contribution from `m≥v+k`, for `k≥1`, by

```math
 C_x v^{\mu-1/2}(1-x)^{v+k}. \qquad\text{(3)}
```

This is a geometric cutoff in *number of births*, not elapsed continuous time. It controls discarded promotion sizes relative to the upper scale in (1), not yet relative to an identified equivalent for `P_v(x)`.

## 2. The lower bound: an exact exponential leaf-family tail

Use the physical cut setup of [scaled/](../scaled/README.md): immediately after label `v` is born, require that it initiate the next jump (probability `1/v`, putting another leaf at its parent), and then that its parent select it at the following jump (probability at least `1/(v+1)^2`). The cut between the marked vertex and its old parent now gives an inert marked side of two ordinary vertices and an inert exterior of `v` ordinary vertices. Both boundary degrees retain their exterior slots. The exterior parent has at least one leaf. Suppress cross arrows using independent thinning clocks kept in reserve; conditional on inert paths, restoring both channels and using no cross costs `exp(−∫(D_mark⁻¹+D_ext⁻¹))`. Time is measured from the cut.

### Marked side: use the Yule limit directly

For this lower event also suppress **all** focal initiations, keeping the rate-one focal clock independent. The side is then a pure star; its leaf count `J(t)` is a rate-one Yule process from one, with `U=lim e^(−t)J(t)∼Exp(1)`. To restore the focal clock it is enough to forbid focal initiations only until a suitable fixed time **if** there is a controlled later evolution. A stronger construction for the polynomial lower bound suppresses only the promotions that affect a designated successful family, as follows.

Run the original inert star-side physical dynamics, and count descendants of its first leaf in an auxiliary Yule tree. A physical parent initiation may remove a family member; if it does, continue that family member's Yule descendants as phantom. For the lower bound, instead realize an auxiliary process in which all focal initiations are disabled until time `h`, and then allow the physical dynamics. Conditional on the star at `h`, the degree mass after `h` has a gamma-minus-loss bound and the side size a gamma upper bound, but their difference at the tail scale matters. The following separate constant-loss estimate supplies the needed precision.

At any inert marked side, let `J` be its leaf stock and let `N` be its ordinary size. The quantity `N−J` counts the focal vertex and all ordinary vertices off its leaf stock. Starting at a star this quantity is one. A focal promotion moves one leaf off stock and places one more vertex off stock; changes in `N−J` due to other initiators are bounded in rate by `N−J−1`. Consequently `N−J` is dominated by a process of per-member birth rate one on its excess above one and immigration of two at rate one. It has normalized limit with expectation at most two. This bounds the loss `S_mark−W_mark` without any independence from the marked mass. To get a lower estimate **conditional on a high Yule amplitude**, use the next conditional construction.

### A conditional high-amplitude Yule comparison

Use a common Yule leaf family from one, which continues independently after any physical loss by giving phantom descendants to removed tags. Conditional on its limit `U=u`, the Yule population is `1+Pois(u(e^t−1))` at time `t`, and across times these increments form a Poisson process of intensity `ue^t`. Indeed the posterior for `U` given the Yule history at time `t` is `Gamma(J(t),e^t)`; applying Bayes to the transition probabilities gives the stated independent increments. All focal ring clocks retain rate one under this conditioning: the auxiliary Yule trajectory has autonomous birth clocks independent of focal rings (use the same clocks for both real and phantom births). Conditional on a schedule of Yule births, the parent of each new auxiliary member is uniform among the current members, with independent choices. These choices likewise do not depend on the focal clock.

Let `L_ph(t)` count phantom family members. At any focal ring, at most one additional auxiliary member is designated phantom. At a scheduled auxiliary birth, the expected **fraction** of phantom members is unchanged, since the initiator is uniform. Therefore the expected final phantom fraction conditional on the auxiliary birth trajectory and focal ring times is at most the sum over focal rings of `1/J(t)` (the bound is valid even for an adapted choice of the lost tag). It follows by Fatou that

```math
 \mathbb E[U-W_{\rm mark}\mid U=u]\le u\int_0^\infty\mathbb E[(1+\mathop{\mathrm{Pois}}\nolimits (u(e^t-1)))^{-1}]\,dt.
```

The right side is `O(1+log(1+u))`, **not** bounded uniformly in `u`. The integral is of order `(log u)/u`, because the period just after the cut has only `O(u t)` auxiliary particles. This is the obstruction to achieving (1)'s lower exponent with a direct constant-loss argument. It can, however, be removed at constant cost by insisting that the focal clock remain silent for one unit of time. Conditional on `U=u`, this event still has probability `e⁻¹`. After time one the expectation bound on the phantom mass from rings is at most

```math
 u\int_1^\infty \frac{dt}{u(e^t-1)}=-\log(1-e^{-1}). \qquad\text{(4)}
```

On the same silence event the expectation bound for `S_mark−W_mark` from the previous paragraph is at most `2e⁻¹`, even conditional on `U=u`. Indeed off-stock initiations are ordinary independent rate-one clocks, separate from the real/phantom family birth schedule; focal clocks after time one are independent of that schedule. Normalize and apply the expectation bound to the dominating process from time one. New untracked leaves, if any, only improve the lower bound `W_mark≥U−limsup e^(−t)L_ph(t)`; off-stock accounting always uses the full actual stock. Thus there exist absolute constants `D_*,L_*` such that with conditional probability at least `1/2` for all sufficiently large `u`, the marked side simultaneously satisfies

```
W_mark ≥ u−D_* ,   S_mark−W_mark ≤ D_* ,   ∫ dt/D_mark(t) ≤ L_* .
```

For the reciprocal assertion, before time one the integral is at most one. After time one use independent rate-one births on the *unconditioned* physical path would be invalid here. Instead on the conditioned Yule path, the phantom fraction is a nonnegative submartingale conditional on the auxiliary schedule, with expectation bounded as in (4). The maximal inequality shows that the fraction never exceeds `1/2` with probability `1−O(1/u)`. Also the count trajectory satisfies `J(t)≥u(e^t−1)/2` for every `t≥1` with probability `1−O(1/u)` by the Poisson process maximal bound (union over dyadic bins in its mean parameter, using its square-integrable compensated martingale). These events imply the reciprocal bound. This argument conditions only on the auxiliary tagged-family schedule; external initiators on the inert side can be run with independent drivers throughout by splitting the directed birth channels.

### Exterior: integrate the exponential tail against the tilted comparison

Couple the inert exterior size below a Yule size from `v`, with limit `S_+∼Gamma(v,1)` as in [scaled/](../scaled/README.md). Let `α=x/(1−x)`. Given its path with `S_+` large, the marked side can achieve success by `U>α S_++(1+α)D_*` and the three simultaneous bounds above: then `W_mark>x(S_mark+S_ext)` with `S_ext≤S_+`. Integrating over independent marked-side `U` gives a lower factor `c_x e^(−α S_+)`. The exterior event must also have bounded reciprocal. Tilt by `e^(−α S_+)`, whose normalization is exactly `(1+α)^−v=(1−x)^v`. Under the tilt, all inert exterior rates are multiplied by `q(t)=(1+α e^(−t))⁻¹≥1−x`, by the proposal-Yule Laplace martingale. An embedded tracked birth–death family with rates `1,1/2` in that time-changed clock, started from the exterior parent's retained leaf, has probability bounded away from zero of a bounded reciprocal integral (use survival and then a fixed finite bound). Since tilted `S_+∼Gamma(v,1+α)` tends to infinity in probability, intersecting with `S_+` large preserves a uniform positive probability. Undo the tilt, multiply by the bounded no-cross factor and the setup probability, and obtain the lower bound in (1).

## 3. What still prevents a sharp prefactor or adjacent ratio

The postpromotion beta distribution bounds the physical endpoint; it is not an exact conditional kernel. Even near time zero after promotion, `D−r` depends on the permanent core parent and which leaves promote, and the deficit depends on the history generating a large endpoint. In the rare event the distribution of the postpromotion environment is therefore further tilted by its own endpoint tail. The exterior size tilt above is exact for the **auxiliary proposal size**, while the no-cross histories furnish only a lower sublaw. Neither can identify the physical Radon–Nikodym weight at the relevant endpoint. A sharp endpoint singularity of a dominating or dominated beta mixture is not a singularity calculation for the actual law.

Also, even an upper/lower polynomial sandwich of the form (1), together with exponential decay of each finite promotion-delay weight, does not rule out arithmetic oscillations in the coefficients. For instance multiplying `(1−x)^v` by a bounded positive two-periodic sequence preserves those scales but can reverse the adjacent sign infinitely often. This is a logical obstruction to inferring the sign from envelopes, not a counterexample for the physical graph process. A proof by summing over fixed promotion delays would need convergence of the birth-labelled, physically tilted endpoint kernel at each delay (and a summable error on the *sharp* scale); no such convergence is proved here. The conditional graph law at birth cannot be replaced by an iid sample of an unweighted local graph.
