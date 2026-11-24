> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Four signed original potentials at one physical trial

**Result.** At `s=.56572365`, four unrestricted signed coefficients cannot make the raw residual uniformly accurate to `10^-8` on physical trees, even after removing an arbitrary common per-core offset. The exact finite obstruction is at least **`.001675953898` per unit of positive stationary denominator**. Even ordinary chains alone force at least **`1.199755e-6`**, and small branching trees alone at least **`1.510010e-6`**. These are uniform finite-control obstructions, not bounds on the stationary average or on the true exponent.

## Raw rows and normalization

The columns in order are original D18 lower and upper, original D20 lower and upper. “Upper” identifies the payload; its column is still the raw `(ℒ−s)Φ` with no sign flip. D18 is the original Neon pin `cec3855cc52ab4265ce3e36b8a7a38f1ae257f77`; D20 is the Claim 3 pin `51f4cfadac675fdf1cef84d44935eebb5456aabb`. The sparse replay checks sizes, headers and pre/post file stamps; whole payload authentication is inherited from the pinned input replays. It uses no routed charges or estimates from sampled trees.

For each actual oriented physical tree `T`, let `r_T=((ℒ−s)Φ_i(T)/C(T))_i`, where `C(T)` counts its cores. The compressed generator uses the actual directed donor reciprocals, adds one leaf to its receiver, and distinguishes a gain by a core from the promotion of one of its leaves. In the latter case it re-evaluates the whole post-event potential including the newborn core. Direct checks enumerate **ordinary** directed arrows on complete labelled trees and use the separate `audit.Cert.potential` implementation. That evaluator does not use the compressed motif enumerator. Rooted stars include the actual founder census, not a planted surrogate.

Put `d=(πΦ_i)_i`. Each component lies in the closed comparison interval below because the original finite-tree comparisons survive stationary normalization. D18 lower endpoints are deliberately weakened to zero. No product-marginal estimate or value of `μ` is needed for this box.

| Column | Lower on `d_i` | Upper on `d_i` |
|---|---:|---:|
| D18 lower | 0 | 1.84124847315 |
| D18 upper | 0 | 1.95820675182 |
| D20 lower | .28368159134 | 1.84219674916 |
| D20 upper | .28863648651 | 2.00854054323 |

For every box vertex `v`, replay an exact representation `(v,0)=Σ_T y_T(r_T,1)`. HiGHS only proposes the support. Fraction Gaussian elimination reconstructs the coefficients and verifies the equality exactly. Interpolate the 16 representations with the usual box barycentric weights. This gives a representation for **the actual, unknown joint denominator vector**, with `Σ|y_T|` at most the worst verified cost `B`. For any signed `w` and offset `a`,

```
|d · w| = |Σ_T y_T (r_T · w − a)| ≤ B max_T |r_T · w − a|.
```

This is homogeneous, allows arbitrary signed coefficients, and assumes neither `sum(w)=1` nor pointwise positivity of the mixture. A mixture with `d·w>0` can be rescaled to one. The absolute test drops the constant column and fixes `a=0`. The proof also applies to any comparison-normalized mixture whose positive stationary denominator is inherited from a finite uniform comparison.

## Controls and costs

The deterministic catalog has 541 entries: homogeneous and alternating-stock chains/stars, smaller asymmetric branches, stocks across the degree-18/20 seams and high leaf plateau, and high-arity stars across both arity seams. No candidate relies on a sampled endpoint. The displayed bounds are rounded **down** by rational arithmetic; the worst cost basis need not be exactly optimal.

| Rows allowed | Common offset? | Exact replay floor, 541-row discovery | Fixed-support exact-only replay floor |
|---|---|---:|---:|
| All | no | .002174579989 | .001121567735 |
| All | yes | .001675953898 | .001228676614 |
| Chains, leaf stock ≤64 | yes | .000001199755 | .000000748925 |
| Small branches | yes | .000001510010 | .000001510010 |

The strong control uses high-arity stars (separate cutoff seams). It is not just a star/overflow artifact: the fixed chain support uses three-, six- and nine-core chains, including a degree-two alternating-stock chain and mixed stocks 17/19. The small-branch support uses four- and seven-core physical trees and no high arity. The 16 frozen support trees can be evaluated in `--fixed` mode with no floating-point basis selection or LP tolerance in the proof. Both paths additionally check complete labelled ordinary trees.

As a diagnostic only, fixing `sum(w)=1` and minimizing the common-offset cost gives the following numerical results. All fits allow negative coefficients. The chain fit extrapolates strongly: coefficients approximately `(18.3685, −.1617, −17.7372, .5305)`, yet improves the best pair by only about 10%. The small-branch fit improves by about 0.3%. The global fit improves about 22% over D20 alone. These values compare forcing contrasts, not certified stationary denominators; the exact box obstruction above does not use this gauge.

| Rows | D18 pair | D20 pair | All four |
|---|---:|---:|---:|
| All | 7.2514e-3 | 4.3929e-3 | 3.4383e-3 |
| Chains | 5.8471e-5 | 5.8507e-5 | 5.2652e-5 |
| Small branches | 9.7175e-5 | 9.7479e-5 | 9.6927e-5 |

This rules out the uniform raw-fit route at this trial, including fits up to an additive per-core offset. It leaves open signed stationary cancellation. Any continuation must use actual occupation or new corrections; exempting only the rare high-arity seam controls will not suffice, because low-arity chains and branches produce independent obstructions. There is no routed-charge ownership claim.

## Replay

From `.` (default inputs are `b3f{18,20}_{lo,hi}.bin`; override with `--blobs` in column order):

```sh
P=./round2/.venv/bin/python
R=round2/occupation-neon20/four-potential/replay.py
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 taskset -c 8 nice -n 15 /usr/bin/time -p timeout 160s "$P" -Bu "$R"
OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/time -p "$P" -Bu "$R" --fixed
```

Full discovery verified 64 rational box-corner representations across the four tests and 168 independent ordinary-arrow whole-tree matches. The fixed support also verified all 64 representations and 92 independent matches, and completed in 9.69 seconds. No giant table scan or input rehash was run.
