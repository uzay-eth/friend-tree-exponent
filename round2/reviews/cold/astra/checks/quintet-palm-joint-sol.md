# Quintet conditioned birth-Palm ratio: cold audit

**Holds on both signs.** The pinned producer vectors certify **`[0.565703539413800, 0.565746271349400]`** against the original sided cut-24 matrix. The exact signed repairs give positive margins of about **0.750302733676 / 0.559140454945 eta ticks** (lower/upper). Both use pure conditioned roots, with no reduced-`f` or old-donor price.

## Physical ratio

Let `s=+1` for the lower witness and `−1` for the upper, `x=s*(mu−lambda)`, and `M` the stationary mean original potential in ticks. Denote by `H(D,c)` the current nonnegative higher-family host cap, bounded by the old global cap. Put `e_fg=cap_f+max(0,cap_g)`, retaining the negative upper `cap_f`. Write `G` for the sum of eligible edge floors, in stationary ticks per core: the exact degree-18 edge floor, the original outside-FY floor, and the lower-only extended FY floor. The signed bound is

```
x*M >= pi(F) + G − e_fg − pi(H).
```

Here `F` is the original singleton slack. Its conditioned Bellman values use the mother of the actual newborn: `V2` if the mother has degree two at the promotion, `V3` if she has degree at least three. The mother's degree is unchanged by that promotion; her arity is not. The source is therefore `b(D,c)*V(D)` on the actual mother's occurrence, with `b=(D−c)/D`. This is not the current degree of an independently sampled mother, the focal newborn degree, nor a projection of the current edge law. Discounting the nonnegative lifetime at the accepted exact `HI>=mu` and using stable age gives

```
pi(F) >= pi(b V(D)),   eta = pi/mu,   eta(b) = 1,
x*eta(W) >= eta[N],    N = b*V(D) + G − e_fg − H.
```

There is no additional `LO` multiplying `b*V` in this formulation. The edge floors are already per core and enter as constants on each occupation column, hence contribute `G*eta(1)=G/mu`. Dividing those floors by `LO` before insertion is not justified. Lower-only extended FY logs are in potential units per core and are multiplied by `2^40` exactly once. The original FY `(175,5,8)` has paying-neighbor degree `<=5`; the extended cut requires paying-neighbor degree `>5`; both require `max(D_host,D_receiver)>18`, so neither uses the degree-18 edge account. The third vertex's degree does not define the exclusion. Higher savings replace the global higher cap once. No post-FY or other outside-edge floor is included in this checker.

The default uses no reduced `f` reward. If two proofs both bound the **same** `pi(F)`, any fixed convex blend is valid. For joint weight `theta` on the original unconditioned dual, replace `b*V(D)` by `(1−theta)*b*V(D)+theta*(b*V0+R_mask)` with `0<=theta<=1`. Here `R_mask` is the original nonnegative affine singleton residual in the shared `u_child+v_mother`, with zero on the degree-overflow sentinel. `b*V0` still applies there. This is a replacement, not a second payment or a cellwise maximum of separate averaged inequalities.

## Signed trial and exact certificate

For a proposed endpoint `z`, use `a=s*(z−lambda)`. On each occupation column compile `N−a*u`: use `u=W_hi` for `a>=0` and `u=W_lo` for `a<0`. The envelopes compare the **same original stationary denominator** after the lower-only degree recharge telescopes; its finite founder/constant correction does not survive normalization. The denominator is positive in the true lift. A nonnegative lower bound on `eta[N−a*u]` therefore proves `x>=a` and the endpoint in the correct direction. A negative lower bound does not disprove the endpoint, and must never be clipped. The tail degree index 25 is a class for `D>24`: retain its true lifted birth coordinate and use `V3` for it, with no exact degree-25 `R`. All nonoccupation columns have zero objective; they still need dual inequalities.

The independent verifier uses the original sided physical cut-24 matrix with saturation 7 and the 12-class Palm kernel. That kernel uses discount `0.565775>=HI`, exact mother degree classes 2 through 12 and `>=13`, and true daughter promotion rewards through degree 34. `beta>=beta*K` concerns physical birth fractions; a Palm slack coefficient in its own row is `−1`. The total lifted Palm slack is `1−beta*K*1<=1`. On the other column classes (occurrences, edges, extra stubs, tail upward receive), the total mass per class is at most `1/LO`. In particular the saturated tail child receive is not an occupation reward but is still on its true shared marginal. No sign assumption on the objective, birth multiplier, or dual bound is needed.

For exact multipliers `y`, with `A z=b_rhs` on the true lift, compute `delta_k=max(0,max_{j of kind k}[(A^T y)_j−c_j])` across **all** columns. Then

```
eta[N−a*u] >= y_birth − delta_palm − (sum of other deltas)/LO.
```

If every exact inequality holds, all defects vanish. The checker accepts `--multipliers` as JSON exact rational strings (or exact integers) in the pinned builder's row insertion order, in ticks. If none are supplied in a custom trial, HiGHS only proposes a dual; rational repair certifies it. The physical kernel and host envelopes are independently rebuilt using the prior cold-check routines, without importing the conditioned-joint producer.

## Same-law perspective and virgin row

For the positive release trials use the single upper host objective vector `U`; the physical ratio relaxation seeks `N·z/(U·z)` on one nonnegative lift `z` with `Az=b_rhs`. In particular it does not separately optimize the birth mixture and denominator. Use `tau=S/(U·z)` and `w=tau*z`, so `A w=tau*b_rhs` and `(U/S)·w=1`. Conversely every feasible perspective point has `tau>0`: if `tau=0`, the mass row together with the birth row gives `sum(s_j*w_j)=0` over occurrence columns with `s_j>=LO>0`, contradicting `(U/S)·w=1`. Dividing by `tau` recovers a feasible original lift, and `N·w/S` is the original ratio. The independently rebuilt floating perspective proposed approximately `1.169414101615e−6 / 2.008650784091e−6`; those floats are diagnostics only. The stored offsets are `5847069/5000000000000` and `10043253/5000000000000`. The strictly positive exact fixed-trial inequalities guarantee the sign needed to use the upper denominator, regardless of the sign of individual numerator columns.

The original virgin age row already bounds the degree-two fraction. For the four `(D,k)=(2,0)` columns its coefficient is `3/2+v`, with `v∈{0,1/2}`; on every other occurrence it is `−s_j`. The mass and birth rows give `sum s_j*z_j=1`. It follows that `sum_virgin (3/2+v+s_j)*z_j=1`. Since the exact degree-two birth flux is `z_virgin/2`, it lies in `[1/(2*(HI+2)),1/(2*(LO+3/2))]`. The checker verifies the column identities and matrix size; there is no extra virgin constraint or substituted birth-law approximation. The earlier [conditioned-root audit](quintet-conditioned-sol.md) supplies the physical lifetime proof, while the attached root log binds the exact class-specific values here.

## Preflight executed

Prefix for each command:

```sh
taskset -c 8 nice -n 17 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 \
  round2/.venv/bin/python -Bu round2/checks/quintet-palm-joint-sol.py
```

* `--side lower --trial 0`: exit 0; `613 rows / 5554 columns`, signed repaired margin approximately `2811152.169353585` eta ticks, `fixed-trial sign certified`.
* `--side upper --trial=-1/1000000 --reduced-blend 1/2`: exit 0; same matrix, repaired margin approximately `6689368.947603217` eta ticks, `fixed-trial sign certified`. This is a control of the negative-trial and replacement branches, not a choice to add reduced `f` to the release recipe.
* Both runs also control the signed verifier on the real matrix with objective `−3*b`, whose exact bound must be `−3`, and a planted Palm multiplier with a two-tick Palm-slack violation.

## Pinned release replay

The producer source, roots and both JSON files match the supplied SHA-256s. Using the command prefix above with `--side lower --stored` and `--side upper --stored` exited 0 on one niced core. Each replay decodes the 613 integer multipliers as `integer*10^6/10^12` ticks in original row order and checks all 5,554 columns. Both the producer's uniform `1/LO` repair and the independent tighter Palm-specific repair give exactly the same margin here because there is no Palm-slack defect. Stub and tail-up defects are also zero. The exact margins are:

```text
lower 246796378766836906157775417183400184888681/328929067814462075504886439151250000000000
upper 337626417760847229096216385970226521841377/603831138983070476772255298925000000000000
```

Exact nominal plus the lower offset is `2828517697069/5000000000000`; nominal minus the upper offset is `2828731356747/5000000000000`. The replay uses only the three specified edge accounts (two on upper); all earlier input receipts remain pinned. Only the two audit files were changed.
