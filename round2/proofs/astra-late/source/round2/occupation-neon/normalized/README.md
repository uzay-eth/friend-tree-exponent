> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Normalized physical occupation for paired D18 slack

**Result.** A 20-state rational occupation relaxation prices the actual D18 first-body slack on degrees 2–5, leaves the omitted age occupation in explicit overflow states, and imposes the exact physical balance `pi(b)=mu`. In the next omitted degree class (`D=6`), the paired contrast has both signs row by row, yet its mean is strictly positive. The prototype does not improve the full exponent bracket: with only these rows, the unpriced all-motif slack budget is still about `152e-6` per core. No multigigabyte artifact is opened.

## 1. Opposing slacks: share the budget before pricing the difference

Write `a=.565673`, `b_*=.565775`, and put `delta_- = (L-a) Phi_-`, `delta_+ = (b_*-L) Phi_+` for the two physical tree drifts. Use `pi` also for the per-core mean of a sum of local charges. The accepted D18 decomposition assigns nonnegative terms to both drifts. With `P_±=pi Phi_±`, the two identities are

```
S = pi(delta_-+delta_+) = (mu-a)P_- + (b_*-mu)P_+,
E = pi(delta_--gamma delta_+) = (mu-a)P_- - gamma(b_*-mu)P_+.
```

These are averaged drift identities, not differentiation of a real-axis asymptotic. Integrate the finite-tree mean potential derivative on `[t,t+1]` before taking the normalized limit, as in the [renewal occupation note](../renewal/README.md). One may also define the per-core slack limits by this argument. The nonnegative budget bounds on residuals need only a retained-charge limit, even if one has not separately made every omitted charge into an integrable stationary observable. A founder-exception contribution can stay in the omitted ledger. Use the denominator ceilings from the renewal note as inputs; this replay does not compute new denominator or ratio estimates.

Partition by **charge and class**, including multiplicity: `S=S_H+S_M`, `E=E_H+E_M`. Stop reading known rows once `S_H >= L_H` and `e_l<=E_H<=e_u` have been certified. If `U_S` is an upper bound on the displayed total budget, then `0<=S_M<=U_S-L_H`. For each omitted disjoint class `i`, suppose `d_i` is its sum of nonnegative opposing charges and `e_i` their contrast. Then always `-gamma d_i<=e_i<=d_i`. A proved narrower row cone `q_i^- d_i<=e_i<=q_i^+ d_i`, a physical occupation enclosure for `pi(e_i)`, or an independently paid pointwise envelope `|e_i|<=epsilon_i h_i` can be added as linear constraints for the omitted classes only. For example, if `q^-<=e_i/d_i<=q^+` on all missing rows, one can use `min(0,q^-)(U_S-L_H)<=E_M<=max(0,q^+)(U_S-L_H)` without any missing-row occupation model. An average contrast of definite sign in one class constrains that class even when its row cone straddles zero.

For a trial `s`, substitute `s` for `mu` on the two right sides, use sound bounds on `P_±`, and impose the common remaining sum budget, not a copy of it for each class. Equivalently, recovered per-side totals are `(gamma S+E)/(1+gamma)` and `(S-E)/(1+gamma)`. Omitted entries cannot be made nearly proportional by normalizing their occupation: assigning all remaining positive slack to either one side realizes the two extremes of the bare cone. This explains why a tiny *retained* contrast alone gives no useful final interval. Conversely, an allocation with calibrated missing-class occupation, absolute contrast bounds, or stronger cones has a direct path to an interval exclusion test.

## 2. Normalize the *occupation*, not the observed prefix

For the physical birth-time measure `beta`, the [stable-age identity](../../../../../astra/source/round2/degree-law.md#birth-times) reads `pi f = mu E_beta int exp(-mu t) f(X_t) dt`. On a true core, let `D` be current degree, `k` core children, `L=D-k-1`, `B=L/D` the promotion rate. From the newborn `(2,0)`, upward moves are `(D+1,k)` at `G=L+rho` and `(D,k+1)` at `B`. If the parent birth degree is `d_p`, `rho` lies in `[0,1/d_p+k/2]` forever because the parent's degree is nondecreasing. The actual incoming hazard can be adaptive; no frozen or independent exterior is assumed. For a degree cut `D<=n`, continue the mass above the cut in states `(tail,k)` up to `k=n-1` and absorb all larger `k` in `(tail,n)`. There, internal degree gains are invisible and the visible promotion rate is between `max(0,1-(k+1)/(n+1))` and `1`. Keep `B in [0,1]` at the absorbing aggregate too; only its transitions become invisible. In particular the entire occupation after exit still pays the test `1` and may pay its share of `pi(B)`.

For an affine endpoint action with rates `Q_a` and age rate `s`, normalized occupation variables `x_(j,a)>=0` obey, for every state `i`,

```
sum_(j,a) x_(j,a) [1_(j=i) - Q_a(j,i)/s] = 1_(i=virgin).
sum_(j,a) x_(j,a) = 1,
LO <= sum_(j,a) B(j,a) x_(j,a) <= HI.
```

The second equality follows by summing the first, including every overflow state. The last is the exact physical `pi(B)=mu` relaxed to the certified incoming interval. A row with an actual adaptive `rho` uses an occupation-weighted mixture of its endpoints. Allowing both global endpoints for `s` independently in each state makes this a relaxation, not a claimed optimal physical policy; `1/mu` lies between the two reciprocals, so split the occupation variables accordingly. The same argument handles the adaptive invisible changes in the overflow states. A joint LP can use separate blocks for birth strata and sum their mass-balance constraints with right-hand sides equal to the stratum weights. In particular this uses the known total `1` to *assign* omitted mass; dividing a suboccupation by its retained mass would give a selected-state law and is invalid for `pi`.

The replay uses the dual, not a floating LP solver. For any rational `t`, obtain `V` by backwards minimization with downward rounding and check every endpoint inequality

```
(s + total_visible_rate) V(j) <= s (f(j)+t B(j,a)) + sum_i Q_a(j,i) V(i),
pi f >= V(virgin) - max(t LO,t HI).
```

Here the sum on the first line uses only off-diagonal successors. Apply occupation balance to the inequality; `-f` gives the opposite bound. The dual search only proposes multipliers; the final inequalities use exact fractions. Constants ±1 replay exactly, including overflow. More moment or conditional-mass constraints can be added by including their payoffs with signed multipliers; the exact constant bookkeeping is the same.

## 3. Pay only for missing classes

An exact-birth-degree mass floor follows by selecting at most one child per newborn, requiring gains from its current leaves up to target parent degree `d`, then the first promotion. Each prescribed nonterminal gain has rate `l` and the terminal promotion rate is `l/(l+1)`. The concurrent adverse hazard is bounded by the other intrinsic transition and exterior `1/2`. With discount `HI`, multiply the factors `wanted/(HI+l+l/(l+1)+1/2)`. This bounds the discounted probability of a selected child from each physical newborn. Counting the selected children by birth compensation, then dividing by total discounted births and taking the stable-age limit, gives a lower bound on `beta(D_parent at birth=d)`. This is the same analytically justified selected-birth input as the [physical pricing note](../../../../../../reviews/cold/astra/occupation-neon/check/README.md); no observed histogram is treated as a bound. Exact `d` bins are disjoint. Conditional on a bin, the virgin occupancy alone is at least `LO/(LO+3/2+1/d)` by the first-residence hazard. More generally recompute the dual value with the narrower parent-rate interval; weight just the improvement over the universal value by the bin floor and keep the remaining birth mass at the universal value. Apply the `pi(B)` dual offset once globally, not as a constraint `E[B|bin]=mu`. As an independent mass bound, the comparison with a rate-one Poisson birth clock observed at an independent exponential age gives `pi(k>=q) <= (1+LO)^(-q)`.

All numbers below refer to the pinned charged first-body rows, regenerated only for the requested degree classes from the sparse JSON. The `D=6` example uses a pointwise mixed-sign payoff, not an assumed rowwise proportionality. Slack columns are in units of `1e-6` per core, with `gamma` the exact opposing virgin-charge ratio (about `1.241277477486`).

| Charge/mass | Flow alone | With promotion balance and birth strata where shown |
|---|---:|---:|
| `pi(D>5)` | `[.2445016101,.4255618213]` | `[.2756378078,.4255618213]` |
| `D<=5` first-body sum lower | `2.7847442376` | `2.8005560613` |
| `D<=5` paired contrast | `[-.0001597383,-.0000126049]` | `[-.0001308304,-.0000166930]` |
| `D>5` first-body floor lower (only `k<=4`) | `.6216839717` | `.6811266279` |
| `D=6` first-body contrast | `[.0000254990,.0000507082]` | same |
| `D=6` first-body sum lower | `.1615586683` | `.1846283073` |

Every `D=6` contrast row is between `-.0016327640` and `.0019353022`, yet its **average** has the positive sign shown. Solve the now-retained `D<=6` contrast jointly for `[-.0000930759,.0000330224]`, and reduce the all-class sum cap to `151.8951272131`. At the earlier cut, the remainder sum is enclosed by `[.6811266279,152.0797555205]` with the first-body floor; the full paired contrast with only the generic missing cone is still `[-188.7733061395,152.0797388274]`. The tail floor and `D=6` sum overlap and are **not** both subtracted from the budget. The large remaining cap identifies the next needed input: certify narrow contrast envelopes or charge-weighted cones on the omitted higher motifs and mixed/tail classes. It cannot be removed by more precision on the retained first-body rows alone.

One can instead pair the two all-`L` floors on the disjoint suffix `D>6, k<=4`. Solve a single joint occupation test for the exact charges through degree 6 and those two floors above it. The resulting paid sum is at least `3.7595087005`, and its contrast is `[-.0002739598,-.0001912629]`. Assign only the **excess over each floor** to the residual on the suffix. The residual cap becomes `151.1208028813`, and the full contrast is in `[-187.5831229559,151.1206116183]` with the generic residual cone. This tightening does not change the [renewal exponent bracket](../renewal/README.md).

## Replay

```sh
/usr/bin/time -f 'wall=%e rss=%M exit=%x' timeout 12s taskset -c 13 nice -n 16 python3 -B round2/occupation-neon/normalized/check.py
```

The default degree-five example completed with `wall=6.43 rss=25880 exit=0`; the same command with `--cap 6` completed with `wall=3.03 rss=25888 exit=0`. The latter also certified positive first-body contrast on the next omitted class `D=7`: `[.0000222243,.0000419394]`. The replay checks the paired row operands against the pinned sparse signed corners and verifies every dual inequality used for the printed values. No GMP verification, full payload hashing, or denominator search was run.
