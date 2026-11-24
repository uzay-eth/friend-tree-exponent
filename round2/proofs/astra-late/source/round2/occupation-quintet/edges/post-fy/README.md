> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT proof audits and Claude replay of earlier freezes; original Neon-2 Gold1 witnesses and inherited stationary results. See the [release index](../../../../../README.md#claim-25) for exact scope and citations.

# Quintet edge slack after all physical FY gaps

**Both original signs give positive new stationary floors, far below the target.** On edges with `max(D)>18`, the certified residual after removing every applicable FY gap is at least **`4.13508081414e−11` lower / `9.82714999612e−11` upper per core**, in potential units. These are about **`8.27e−6` / `1.97e−5` of the `5e−6` target**, respectively. These are lower bounds on a covered part of the residual, not upper bounds on the full account. No endpoint or board payment is applied.

| Disjoint paid edge states | Lower ticks/core | Upper ticks/core |
|---|---:|---:|
| Both `D<=26`, or one `D<=4` and both `D<180`; always `max(D)>18` | 37.6057966843 | 90.1621611172 |
| One `D>=180`, other `D<=12` | 7.8597599681 | 17.8881708807 |
| Both `D>=180` (entry route restricted) | 0.0001377170 | 0.0003248884 |
| Total (approximately) | **45.465694369** | **108.050656886** |

The target is `5,497,558.13888` ticks/core at `S=2^40`. The total receipts below, rather than rounded display values, specify the exact floors.

## Residual and disjoint accounting

Use the original outward verifier sign `sigma = −1` lower, `+1` upper. Start with the [original physical edge allocation](../README.md), including the original centered `f→g` transfer and the post-promotion newborn H3 term. At each applicable non-M3 exact G2 endpoint `x`, with opposite physical mark `y`, let `j=Δ_x^L g(x,y)` and `A_x` be its original FY interpolation. Replace the signed donor contribution for **each** external incidence by

```
M_x(j) = max_knots sigma * (xi*j - A_x(xi)).
```

The signed primal `sigma*(c_x−1)*A_x(1/D_y)` is unchanged. The maximization uses the actual physical gain jump (including the 179→180 seam) and includes the sign **inside** the maximum. Repeat independently at the other endpoint. The actual reciprocal donor belongs to a knot segment on which the expression is affine, so `G_relaxed = G + sum(all physical FY gaps)` identically. Thus `cap_g−G_relaxed` uses precisely the residual after paying even FY gaps that earlier masks have not priced. With the `max(D)>18` guard, this is disjoint from the original degree-18 `g` slack and from any non-overlapping collection of old and new physical FY atom floors. It takes no higher-row, singleton, or extra Jensen reserve. This does not establish disjointness *among* earlier FY ledgers.

Outside FY, genuine external receive is left live (except the original M3 `g→H3` transfer); the same outside donor coordinate enters the residual and the physical transition. Arity `c=1` has no external donor; bare states have zero promotion. Inner overflow uses the actual `D,c` in all rates and the capped payload ID only in the weights. No state or edge is counted twice when its two endpoint labels agree.

## Large degrees, infinity and the cap

For a mixed edge with tail degree `D>=180`, use `t=180/D`. The physical jump of the finite endpoint is `(1−t) j_inf+t j_knot`, and FY primal interpolation is exact on `[0,1/180]`. Convexity of the maximum of signed affine functions gives `M(j_phys) <= (1−t) M(j_inf)+t M(j_knot)`. The finite endpoint's adjacent-tail send contributes `t*sigma*j_phys/180`, bounded above by `t*max(sigma*j_inf,sigma*j_knot)/180`, just as in the verifier's tail row. The centered `f` transfer is affine; the projected H3 birth is affine as well since its endpoint seam is below 180. The tail's physical gain and promotion become `t*gamma*Δg` and `t*beta*Δg`, with `gamma=r/(D+1)` and `beta=(D−c)/D`. For exact tail arity, `gamma` belongs to `[(180−c+1/D_finite)/181,1]` and `beta` to `[(180−c)/180,1]`; saturated arity `c>=27` uses `gamma∈[0,1]` and a zero promotion jump. These are the authenticated original verifier ranges. Each lifted corner is bounded above by the signed `g` cap, including the INF corner where the finite endpoint can still evolve. Hence `G_relaxed<=cap_g` even on mixed tails. No maximizer is fixed across different physical states.

The mixed-tail **payment** uses the safe projected upper bound on `G_relaxed` just described: discard the nonnegative INF residual, and retain `t` times the certified tail-knot residual with live `gamma`, `beta`, and finite outside receive. That choice discards any extra savings from the convexity or adjacent-send envelope. For two tails, all **784** G2 tail/INF weight pairs are zero on each payload. Endpoint changes cannot leave the tail block; neither endpoint has FY, centered `f` transfers or supported newborn H3. Thus the exact row is zero and its whole post-FY residual is `cap_g` at every physical double-tail degree and arity.

## Signed physical balance

For each paid mask build a nonnegative old-edge subsolution `V` at the original **exact** `HI` and a promoter subsolution `W`. The two Bellman inequalities are `HI*V <= R_mask + A_edge V` and `HI*W <= b*V(post-promoter,J)+A_vertex W`; the second uses the new edge exactly once and evolves the original promoter. On summing edge and core clocks, stable age yields `pi(R_mask) >= pi(b V_birth) >= LO*W(J)` at the original exact `LO`: the nonnegative `(HI−mu)` terms can be dropped. Zero extensions leave a nonnegative incoming generator term wherever a test resumes from an unpaid state. All paid scores and test potentials have bounded local drift on their supported degree/arity regimes. In the unbounded tail, `t=180/D` bounds the drift contributions of tail gain, and the other endpoint stays finite. The double-tail test is constant after entry. The founder has no limiting age mass.

* **Finite solve.** The mask contains every physical arity, including exact-degree inner overflow in the wing. Visit unordered endpoint pairs in reverse gain/promotion order; test all four receive corners with exact rational arithmetic. The initial `max(D)<=18` corridor earns zero but propagates future outside credit. The reward uses the exact physical FY maximum, not a projected upper bound.
* **Mixed solve.** Use `V(D,c;y)=t*U(min(c,27),y)` with all 27 tail arity classes and finite tip degree `<=12`, all tip arities. Dividing the Bellman residual by `t` converts tail gain into `−gamma*U`; tail promotion is `beta*(U_next−U)`. The finite endpoint's adjacent gain is `t/180`; its whole transition uses `D_finite−c_finite+(c_finite−1)u+t/180`. Testing both endpoints of each rate interval and `t∈{0,1}` is sufficient by separate affinity. The edge bridge propagates the continuation across the exact `179→180` transition and down to `J`; its larger finite endpoint is restricted to arity at most 10 before the seam, with full arity on the short side. No bridge reward is paid. Tail promoter births are also covered via a `t*W(c)` suffix through arity 26, with the promoter overflow continuation safely set to zero.
* **Double-tail solve.** An already-tailed pair keeps constant value `cap_g/HI`. Register from births whose promoter is already at `D>=180`; the new tip follows its minimum degree-gain rate while at arity one and is killed at its first promotion before reaching the tail. The corresponding promoter entry follows the same killed route from `J`. The quotient is exact, and increasing gain helps because the successor is always greater. This counts only a conservative entry route into the full positive double-tail reward.

These masks are disjoint **at the time of reward**; sharing tracks, bridge states, or physical birth clocks does not duplicate reward. The uncovered medium–medium finite region and mixed edges with larger finite tips were not priced here. No all-support upper control is claimed.

## Replay

```sh
taskset -c 8 nice -n 17 /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  python3 -Bu round2/occupation-quintet/edges/post-fy/replay.py \
  --finite 26 --wing 4 --tip 12 --bridge-stock 10
```

The successful replay exited 0 (`wall=143.67 rss=329752`), authenticating the original receipt and both snapshots. On each sign it checked reversal, full gap removal and cap domination on 1,296 exact pair/receive fixtures (540 finite, 315 mixed, 441 double-tail), and checked exact Bellman corners for 203,256 finite pairs, 2,079 mixed cells and 131,593 bridge pairs. Earlier exploration with wing 8 finished the lower side but hit a 240-second limit during the upper side; it is not a matched result.

Exact total floors in ticks/core (double-tail product rounded down at `10^-12` tick before addition):

```text
lower 5676522650555158448184312449/124852874882500000000000000
upper 53961740580850154833509463977/499411499530000000000000000
```
