# Pressure memory, threshold flux, and reaction time — cold review

**Verdict: holds for the labelled mechanics, marginal lift, and the completed parent proof's fixed-P4 representation error `<5.874e-9` per core.** The old-motif improvements replay exactly as conditional stopped-entry comparisons. The P4 result bounds the error of a moment representation of one stationary incident score, not its unknown stationary value or a change in lifetime kinetics.

## Completed parent proof: adaptive domination and the fixed P4 budget

Reviewed source (`occupation-neon/pressure-memory/README.md`) SHA-256 `329c722931f9fcab8a0d05b4d875f4cfc192c69da06b06736d6d253cac24ac79`, checker `766b83aaa8d85d4d41b3c79388fb16b45118f6dc3edd1d535958aba050a4f22e`. The coupled physical process may depend on *both* dominating clocks. The clocks themselves do not: pre-sample an autonomous pure-birth `Y` of rate `Y-1/2` and an independent unit Poisson `N`, with fresh predictable-thinning marks and independent randomization for the other arrows. At a `Y` event, accept a physical arrow into the tracked core with probability `g/(Y-1/2)` and choose its sender proportional to the exact arrow rates. At an `N` event perform the tracked core's real ring, with uniform neighbour. These are disjoint arrow classes even for mutually adjacent cores: the tracked core never sends to itself, and a send cannot increase both its source and receiver degrees. Reconstruct the remaining arrows at their conditional physical rates; the unaccepted `Y` jumps are auxiliary only. Conditional on the **past** of the joint process, the two proposals retain their prescribed intensities; revealing their entire future at time zero would not justify that assertion. The true marginal generator is recovered, hence feedback can correlate physical `D,k` and the acceptance history without correlating the autonomous `Y,N`. Initially `D=Y=2`, `k=N=0`; each degree gain uses a `Y` jump and each promotion uses an `N` ring, so pathwise domination persists. The comparison clock at the common observation age need not give independent `Y(A),N(A)`; independence is used at deterministic time or between the Yule hitting time and the Poisson process before age mixing.

For `τ_d=inf{t:Y(t)=d+1}`, `P_d=E exp(-sτ_d)`, and differentiating the finite product gives `-P'_d=P_d Σ_(j=2)^d (j-1/2+s)^(-1)`. At independent `Exp(s)` age, the expected Poisson count on `{age>=τ_d}` is `E[(τ_d+1/s)exp(-sτ_d)] = P_d H_d`. Couples at `s<=μ` are legitimate for these increasing positive tests, including arity weighting. In the telescope let `r_d=P_(d+1)/P_d`, `i=J+1`. Multiplying the weighted decrement by `d(d+1)/P_d` gives exactly `a_d H_d - d*r_d/(d+1/2+s)`; the subtracted quantity is less than one, and both `a_d` and `H_d` increase. This proves the denominators `a_i` and `a_i-1/H_i` in the two tail bounds. The first nonzero layer is `d=i`, not `J`: the donor of degree `J+1` has zero excess, and degree `J+2` has exactly `1/[(J+1)(J+2)]`. Monotone positive layer summation is enough; no hidden boundary cancellation is needed. Directional rounding is safe: round `P_i,H_i` up in the numerator and `H_i` down inside the correction denominator.

Independently read the **raw P4 potential gain differences** in the two original D18 payloads. At the prescribed inner labels, all exterior slopes are just these gain differences with the witness drift sign (received transfers depend on the included slots, not on exterior reciprocals). The signed lower slopes are `.00374338137,.00447940521`, signed upper slopes `-.00214205425,-.00253320197`, in parent/hidden-child order. Their differences are exactly `.00588543562,.00701260718`; these are at multiplier one, not the stopped-flow `1.2415`. Fix a hidden side child: its parent determines B, the grandparent determines A, and at most two *other* children can be C, each with a unique D. Fix an exterior parent: each A is one of its `k` core children, has a unique B, and B has three choices of C. The independent fixture attains both bounds simultaneously with two matching A children, giving six parent uses and two uses per hidden side child. The orientation excludes reversal and counts each occurrence once, without dividing by three or multiplying by the path length. Reversal on finite trees followed by the existing arity integrability proves the unweighted side term and the `k`-weighted parent term; no distribution is conditionalized on the occurrence.

Use slabs `[2^r,2^(r+1)-1]` for `r=1,…,14`, and `[32768,65536]` for the fifteenth. Each integer `2,…,J` is assigned exactly once. At shared endpoints the two even-order polynomials even agree, though the incidence measure must still assign the endpoint to only one slab. The response is a degree-nine polynomial because the numerator vanishes at zero; `T_10(3)=22619537` and relative error is at most its reciprocal on the full real slab. The tail constant consumes one extra incidence-weighted count. For low degrees the absolute weighted reciprocal budget is `c_c+3 c_p(1-s)`, by `2p_y<=1` on the side role and by `π(k/D)<=1-μ` on the parent role. These same two weights may be pooled in one positive measure at the stated multiplier. Four-census row intercept times occurrence mass is separate. An independent 30-decimal directed enclosure gives `5.873994280771333e-9`; the exact rational assertion is `<5.874e-9`. The original checker asserts the looser `<6e-9` but its computed number also lies below the sharper bound. The result does not cover different censuses, slopes, or approximating the transition generator. No lifetime DP was reproduced.

## Physical generator: two distinct kinds of receive rate

For nonfounder `x`, write `p_x=1/D_x`, `R_x=Σ_(y child) 1/D_y`, `B_x=k_x/2-R_x`, `g_x=L_x+p_parent+R_x`, `b_x=L_x/D_x`. The founder omits the parent term. An incoming arrow to `x` changes `(D,k)` to `(D+1,k)` and leaves `B` fixed, even if its source is one of the children. A promotion at `x` changes `(D,k)` to `(D,k+1)` and registers a new child of degree two, with **zero** contribution to `B`. The child's new vertex test uses its actual newborn state; in particular `f(B)h(D,k)` has newborn value `f(0)h(2,0)`, not generally zero. A gain at an already registered child `y` increases `B` by `e_y=1/[D_y(D_y+1)]` and leaves `x`'s census fixed. Its total gain rate is `g_y=p_x+Q_y` with `Q_y=L_y+R_y`, all live outside contributions, not merely the arrow from `x`. Thus

```text
L(f(B)h) = f(B)(g_x Δ_D h + b_x Δ_k h)
           + h Σ_y (p_x+Q_y)[f(B+e_y)-f(B)].
```

Parent gains and child promotions can change later hazards but do not jump this statistic. No product cross-increment exists for the two factors: the receiver of a birth is unique, and the newborn's registration adds zero deficit. The standard-library fixture enumerates 436 labelled physical trees and checks 25,218 exact product/threshold evaluations; all of `drop-Q`, `wrong-newborn`, and `fake-simultaneous` are rejected. The formula allows nonlinear `f`; the marginal closure below uses the linear `f(B)=B` case only.

## Exact layers and joint edge cuts

For threshold `j>=2`, put `Y_j=# {child D_y>j}` and `Z_j=k-Y_j`. Then

```text
L(hY_j) = Y_j (g_x Δ_D h+b_x Δ_k h) + h Σ_(D_y=j) g_y,
L(hZ_j) = Z_j g_x Δ_D h + b_x[(Z_j+1)h(D,k+1)-Z_j h(D,k)]
          - h Σ_(D_y=j) g_y.
```

These are **not** interchangeable registration rules: `Z_j` receives one at birth and `Y_j` does not. Add the rows to obtain the exact generator of `hk`; taking `h=1` gives `μ π[Z_j]=μ-E_edge[1{D_y=j}g_y]` and `μ π[Y_j]=E_edge[1{D_y=j}g_y]`. Exploit the whole family jointly: for a test cell, allocate the source as `E_edge[h(parent) 1{D_child=j} j w_child]` with `w_child=g_child/D_child`; share both this flux and threshold counts with the **same** edge measure. Cut child bins at the exact threshold and its predecessor so that `D_child=j` is an actual atom, not a guessed fraction of a bin. Products `g_x Y_j` need the matching receiver exposure (`Y_j` and its product with receive rate); adding only a free threshold source does not close those rows. Complementary and nested counts satisfy `0<=Z_j<=Z_(j+1)<=k`. The exact identity `B=Σ_(j>=2)Y_j/[j(j+1)]` follows by telescoping. Retaining `j<=J` makes the error `Σ_(D_y>J) [1/(J+1)-1/D_y]`, with **true** `1/D_y` in the tail. Hence `0<=B-B^(J)<=Y_J/(J+1)`.

For the base global minimum of a pressure-rate factor `T_x=LB_x=Σ_y g_y/[D_y(D_y+1)]`, the only unrestricted lower bound is zero, attained at `k=0`. Even at fixed `k>=1`, take all child degrees to infinity and use `g_y<=D_y` to see no positive uniform bound over those trees. The physical bound `p_x Σ_y 1/[D_y(D_y+1)]<=T_x` is pointwise and can yield a positive lower price on a finite child mask. For the omitted source `T_x-T_x^(J)`, the crossing layers give exactly `Σ_(D_y>J) g_y/[D_y(D_y+1)]` (the gain at `D_y=J` is already paid in layer `J`). If `|h|<=H`, edge reversal and the one-parent child marginal give

```text
π(|h| (B-B^(J))) = E_edge[|h(parent)| 1{D_y>J}(1/(J+1)-1/D_y)],
π(|h| (T-T^(J))) <= H π[1{D>J}/(D+1)].
```

Both omissions are nonnegative for `h>=0`; for signed tests use positive/negative parts or a uniform absolute bound. In particular the source error is at most `H π(D>J)/(J+2)`. The degree-gain Yule comparison `g<=D` and stable-age integration give the elementary bound `π(D>J)<=Γ(2+μ) Γ(J+1)/Γ(J+1+μ)` for `J>=2`, using the Laplace transform of the degree-two pure birth hitting time to `J+1`; this is `O(J^-μ)`. No unweighted degree moment is needed. For bounded parent masks the stock error is then `O(J^(-1-μ))` as well. For paid terms carrying growing factors in `k` or `g_x` these estimates need the corresponding weighted endpoint moment; they are not blanket error bounds for unbounded reward. A bounded finite census `h` with capped differences needs only bounded rates on the active gain sites (and `b<=1`). More generally `π` has exponential `k` moments by the prior edge-flow proof; using both the parent moment and the degree tail requires a justified weighted bound (e.g. Hölder), not independence of endpoints.

## Marginal lift and old motifs

Reviewed `pressure-memory/marginal/replay.py` in the writer checkout, with 2f74. Its linear-memory source is `w F(q)`, `q=1/D_child`, `w=q g_y`, `F(q)=q/(1+q)`. On a fixed reciprocal interval the chord and tangent give lower/upper affine envelopes for `F`; lift `q` with common endpoint weights and `w` with interval endpoint weights. At the given actual pair `(q,w)`, tensor interpolation preserves `q,w` and both products with the affine envelopes; convexly mixing those two lifted rows gives `wF(q)` exactly. Independence here is only a convex representation of a specified pair, not a statement about physical occupations. Project `w` back to the child census with the same edge mass, including the zero reciprocal endpoint. Bounds `b<=w<=(1+b)/2` use `g=L+ρ` and `ρ<=(k+1)/2`; there is no substitution `w=0` at the infinite-degree endpoint. In saturation the first moment `k` is `SAT` plus a separate stub count, the outgoing weighted parent slots are `kp=1-b-p`, and the **entire** upward pressure is carried by the separate nonnegative `tailup` column. The respective total masses are bounded by `π(k)=1` and `π(ρ)=1-μ`; these bounds justify the dual repair, not an upper bound on each local unseen incoming rate. The `B²` triangle represents the exact joint moment with `B` and hence `(L+p_parent+k/2)B-B²=g_x B`. Since dropping the 21 memory equalities leaves the identical columns at `D24`, the mathematical feasible interval can only shrink; separately repaired printed endpoints need not nest at noise scale. The bounded marginal run reached both column counts, but did not finish its full on/off output, so this review makes no new numerical dual improvement claim.

The old-motif replay passed all 54 comparisons, 51 strictly inward. The known degree bin defines both the reward and the receive intensity; a change at the donor is a different receiver event with rate in `[1/D_x,1/D_x+m-1]`. Tail reciprocal bins never regenerate a finite known degree, and registering a fresh side child enters at degree two. The two-tag-to-four-tag source and continuation are paid at the same promotion event: an old cohort continues and can emit an additional tagged prefix, but the action clock is counted once. `P4` through-terminal at degree 18 changes the lower price from `-9.7773179284e-9` to `+1.3234289140e-9`; `Yt` through-terminal cuts the upper from `2.6928114459e-8` to `1.3475663149e-8`. These are per-entry values; a first-entry intensity and disjoint owner mask are still necessary to pay a stationary ledger. The broad promoted-donor gain cap there is safe, but when the monitored donor is degree two the following physical comparison is sharper.

## Uniform degree-two reaction time

Track a nonfounder from its promotion to a core, and let `p0` be the current reciprocal degree of its parent (or any upper cap). The parent degree never falls. In state `(D,k)=(2,0)` the donor has a single ordinary leaf: degree gain is at most `1+p0`; its exact rate `1/2` promotion moves it to `(2,1)` **without** changing degree. In `(2,1)` its other neighbour is a core of degree at least two, so its exit rate from degree two is at most `1/2+p0`. There is no return after the first gain. For the residence reward `1{D=2}`, comparison with the two-state killed chain gives the uniform conditional discounted lower price

```text
V_bare = 1/(s+1/2+p0),
V_birth = [1+(1/2)V_bare]/(s+3/2+p0).
```

The signs can be checked directly in the two generator inequalities at the maximal gain rates, using nonnegative values and exit value zero; this is valid under adaptive exterior evolution. At `s=0` the comparison gives mean wait at least 1 in the worst parent class `p0=1/2`; integrating the chain at age `Exp(μ)` gives normalized physical degree-two occupation at least `LO*V_birth(HI,1/2)>0.2908713732985`. Conditional on a Palm birth with parent degree at least 18 the corresponding bound is `>0.3855630891860`. That is a **birth-mark** conditional statement, not the same as conditioning a stationary tree on the current parent's degree. No positive reaction **speed** uniformly across all old hosts follows: bare local ports can have arbitrarily small incoming reciprocals. The result is a minimum residence/delay, not an upper bound for the waiting time or a mixing estimate.

## Replay and coordination

```sh
nice -n 15 taskset -c 3 timeout 15s python3 -Bu round2/checks/pressure-memory-sol.py
```

Passed: 436 labelled trees, 25,218 exact generator/threshold evaluations, rational edge-lift and reaction checks. The original `old-motifs/replay.py` passed on CPU 8 under its documented 30-second command. I sent 2f74 the threshold source rule, marginal/saturation review, and verified reaction bounds. This review is separate from the completed Neon-2 response-transcendence attack.

```sh
nice -n 15 taskset -c 3 timeout 15s round2/.venv/bin/python -Bu \
  round2/checks/pressure-memory-sol.py --parent-tail b3f18_lo.bin b3f18_hi.bin
nice -n 15 taskset -c 8 timeout 15s python3 -Bu round2/occupation-neon/pressure-memory/check.py
```

Both new runs passed. The independent mode reads only sparse D18 cells with header checks and inherits the prior payload pins; it does not rehash or rerun GMP. The parent checker reported 234 telescopes, 12 partial-tail comparisons, 153 core shapes, and 33 reciprocal-polynomial controls.
