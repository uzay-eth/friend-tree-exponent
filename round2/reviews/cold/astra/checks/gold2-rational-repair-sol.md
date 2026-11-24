# Cold audit: repaired Claim 28 (original Gold2, no LP)

**HOLDS** on `rationals/repaired/README.md` blob `43f1bc6dd9574096b60cfad82936d26f44a43d60` and `check.py` blob `5f40e60547dd4e4d07e368126748e45d1e16eec1`, conditional on the original Gold2 certificate and physical stable-age law. The repaired enclosure is

```math
 \left[\frac{10043824822832044427}{17754235186087500000},
       \frac{127993535256830043}{226242120050000000}\right].
```

The final writer text correctly uses the signed adverse H3 baseline `P_0` and `F(z)=[z−P_0]_+−[−P_0]_+`. On `z≥0` this is nonnegative, convex, and zero at zero for **either** sign of `P_0`. Hence `Σ_i F(z_i)≤F(Σ_i z_i)≤[−P_physical]_+`. This matches the actual sibling solver and its independent audit; no payment or endpoint changes.

## Four charges, one potential

The absolute D430 singleton floor is 207,990 ticks/core. Its lower replay gives `207990.215716345`; this is not an increment over another singleton bound. Small and older FY contribute 395,261 and 52,341 in disjoint ordered payer slots (payer `D≤5` versus `6≤D≤34`, exact G2 outside M3). The actual-sibling lower is 251,635 ticks, rounded down from `16363762149897743918379249111/65029541606195200000000`. That charge is on the original H3 row after transfers; it spends at most the negative part of that row, leaving the original nonnegative cap valid. The three insertion orders and the multi-sibling bound are checked by the independent sibling replay. Physical vertices can overlap across families; none of the credited signed slots does. No threshold H3, old 699,435 credit, upstream, seam, D18 edge, or higher cap rebate is used. The four floors give `907227−19217=888010` net ticks/core after the original error, once.

The upper comparison is one host allocation of the **original** potential. The singleton goes to its core, G2 to the child with zero allowed for the absent incoming edge, H3/P5 to their middles, B5 and S3/S4 to their hubs, and each P4 splits equally between its two internals. Independent reconstruction reads the actual original weight slabs following the pinned wire layout. It sums the singleton, the nonnegative G2 endpoint maximum, and nonnegative motif maxima with multiplicities `C(c,2)` for H3/P4/P5, `C(c,3)`, `C(c,4)`, and `c C(c−1,2)` for B5. P4/B5 maxima account for the partner's `c_partner−1`; P5 accounts for both partners. P4's half cancels the factor of two in `c(c−1)`. The maxima are of weights, not signed generator caps; the lower comparison's `alpha(c−2)` recharge is not included. The 299 `c≥1` entries match the frozen ceiling export exactly. The isolated `c=0` case is just an original singleton and has no newborn source. A founder at `c≥1` is covered even without an incoming edge. Apply the bound at every host, and use the original pointwise `m_max` outside the support of the saving; in particular there is no degree restriction on the original tree or the original drift.

## Physical drift and birth normalization

At an existing core, a birth to the core changes `(D,c)` to `(D+1,c)` at rate `g=L+Σ_{core neighbours} 1/D_z`, between `L` and `L+c/2`. A birth to a current leaf of that core is its promotion, moves the existing core to `(D,c+1)` at rate `p=L/D`, and creates a new core `(2,1)`. These events are disjoint; unchanged projected coordinates have zero generator increment. At `L=0` the promotion term is exactly zero. At the stopping boundary the gain continuation is zero, and degree monotonicity prevents re-entry. For `B=m_max−U` on the retained degrees and zero beyond, the floored recurrence verifies `(HI−A)V≤B` for the generator `A` of the **existing** core, including non-Markov exterior environments, because this residual is affine in the one real ranged gain rate. The nonnegative `V` also makes this a lower-occupation test at the actual discount `μ≤HI`.

One can derive the per-core normalization directly, without a population-amplitude factorization. Write `F_V(t)=E Σ_existing cores V` and `K(t)=E C`. If `b(t)=E Σ L/D`, then `K′=b(t)` because each promotion creates exactly one core. Dynkin gives `F′_V=E Σ A V + b(t)V(2,1)`. Integrate each identity against `e^(−μt)` over `[0,T]` and divide by `T`. Boundedness and `e^(−μt)K(t)→k>0` kill both boundary terms, giving `π(A V)+μ V(2,1)=μ πV` by local mean convergence, and hence `πB≥μ V(2,1)+(HI−μ)πV≥LO V(2,1)`. Equivalently the stationary age integral uses **one** factor `μ` for physical births per counted core and discount `μ`. There is no extra Gamma factor in the CT comparison. Passing through a fixed-size Yule normalization multiplies the common leading population amplitude by the same factor for numerator and denominator; it does not change this per-core law. Localized expectations are justified by finite-time Yule domination and bounded tests. The singleton and higher rewards use their existing physical occupation bounds; the small FY donor's unbounded degree contributes via a bounded reciprocal with uniform truncation error per host.

For the original lower potential `P`, use its generator inequality after restoring the four payments. Multiply the expectation by `e^(−μt)`, integrate, and divide by `T`. The boundary goes to zero by the original comparability; positive net first forces `μ>λ`. Only then substitute the upper limit `k(m_max−πB)` for the weighted potential average. Thus the ceiling `M=ceil(m_max−LO·V(2,1))` gives `μ≥λ+888010/M`, with no need for a limit of the signed total potential. All units remain ticks, with `2^40` used only to convert the FY receipt fractions to ticks.

## Arithmetic and replay

| Cut | `V(2,1)` | `M` | Exact trial margin at `99/175` |
|---:|---:|---:|---:|
| 6 | 619269501408 | 1458274085072 | `3681682957/781250` |
| 18 | 680528936614 | 1423618796411 | `321296970031/12500000` |
| 24 | 686326891945 | 1420338814887 | `2422917810989/87500000` |

The (closed) repaired interval stays strictly between `99/175` and `142/251`, whose determinant is one. For any interior reduced `a/b`, the positive integers `175a−99b` and `142b−251a` give `b=251(175a−99b)+175(142b−251a)≥426`. Direct reduced-fraction enumeration finds none through 425 and only `241/426` through 426, for each cut. The frozen height scan at the exact D18 lower endpoint reports degrees 2–5 sharp heights `63,15,6,3`, unrestricted denominator `426`, and square-discriminant denominator `94251` (`53320/94251`).

- `round2/.work/venv/bin/python -IBu round2/checks/gold2-rational-repair-sol.py`: writer hashes authenticated, 2,792 physical fixtures and 32,580 directed arrows; original lower payload freshly SHA-256 checked; 299 matched `U` entries; all 40/340/598 exact Bellman corners passed; exact endpoints and denominator scans agree.
- `python3 -IB round2/occupation-quintet-gold2/rationals/repaired/check.py`: all three cutoffs and Farey exclusion passed on the pinned writer checker.
- `nice -n 17 round2/.work/venv/bin/python -IBu round2/occupation-quintet-gold2/lifetime/check.py --conditioned 430 --side lower`: `86742` conditioned states, `750` corners, `207990.215716345` ticks/core.
- `nice -n 15 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 round2/.work/venv/bin/python -IBu round2/checks/gold2-lower-sibling-sol.py`: matched the exact sibling-only receipt, including `2048` labelled two-sibling rows.
- `python3 -IB round2/checks/gold2-fy-sol.py`: exact FY floors and 30 sparse cases/24 registrations passed; that old script's separate H3 composition output is withdrawn and unused.
- `python3 -IB round2/open/arithmetic/transform/sensitivity/check.py --interval 10067018953060470031/17795234955137500000 127993535256830043/226242120050000000`: reported the height and denominator results above.
