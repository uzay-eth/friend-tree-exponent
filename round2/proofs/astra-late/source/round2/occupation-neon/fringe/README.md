> **Author/model:** Astra (GPT family). **Pinned source commit:** `2afedb1e25a6b95be30752c67022f792730d6f21`. **Review scope and dependencies:** GPT physical pricing audits and scoped Claude archive replays and arithmetic; the applicable Neon witness, stationary growth, and physical occupation. See the [release index](../../../../README.md#claim-21) for exact scope and citations.

# Physical 19-mark pricing of the Neon D18 residual

**Outward bracket `[0.565676250330, 0.565772131313]`.** This pays **6.00%**
of the original published Neon gap of `0.000102`, including the previous
disjoint payments. Relative to the existing size-12 higher-motif bracket,
the inward moves are `6.52672e-7` and `5.46556e-7`. The new numerator
payments are at least `6.20406268e-7` and `3.95200834e-7` per core, in
potential units, for the lower and upper witnesses respectively. The
large majority of the gap remains unpriced.

## Actual rooted tolls and the parent split

Use the original D18 binary potentials, their nonnegative **charged**
decomposition, and the two recorded full GMP verifications. `charges.py`
reconstructs all eight rooted families `f,g,hV,hC,P4,Yt,Yb,S3`, not
unrecharged whole-tree differences or a pointwise slack ratio. The
potential is contracted with the pinned Bernstein controls at the exact
physical degrees. On each old embedding the source from a newborn
extension is assigned before adding transfers; for each enlarged
embedding its incoming terms are exactly the connected tip deletions.
Transfers from a tail receiver use `ψ` for `f→g`, and the actual
`T2` and `T3/Yref` combination for `g→hV/hC`. Level-3 to level-4
transfers use the original inner-support rules. First-body and all other
families are separate channels. Equal tips contribute once per actual
embedding, with the correct multiplicity. `check.py` independently
reconstructs complete potential-before/after drifts. On 120 trees per
side, including both-tail and leaf/arity-cap fixtures, the sum agrees
exactly, all eight families are exercised with nonzero charge, and a
one-sided omitted tail transfer fails on four fixtures per side.

After evaluating all in-crown donor degrees, a sum assigned to an apex
has form `q(S,p)=a(S)+p b(S)` where `p` is the *live exterior parent*
reciprocal. Pay `m(S)=a+min(0,b/2)` on this crown, and pay
`p b(S)−min(0,b(S)/2)` one level up on its actual parent's crown. Both
are nonnegative. The latter sees the recipient crown and the donor degree
simultaneously. Summing over physical cores gives back the actual toll,
with at most the negligible founder boundary. In stationary occupation
this is the depth-one distinguished-position reconstruction, not an
assumption about independent marks or a substitution of birth degree for
current degree. Bounded restricted tolls suffice; the construction does
not need a moment bound for the omitted residual. For the disjoint
payment the size gate is on the **charged** crown in both terms, including
when the second term is paid at its parent.

## Positive 19-mark occupation and stationary mixing

Apply the true killed residual inverse to the nonnegative crown tolls
and to the 19 positive direct-birth sources. The denominator is the
actual `s+a+n−1/D`; old/plus copies have rates `a_Y−1/D` and `1/D` and
can reset to retained smaller crowns. The finite solve retains size at
most 8 plus a selected root-degree/young-arm strip (`D≤22`, total size
of descendant arms at most four). Each retained intrinsic move grows
size; its triangular block can be eliminated by backward substitution,
leaving a positive copy-feedback Schur system. The pilot solves the
equivalent sparse full system; it does not form the reduced dense
matrix. After rounding each normalized rate down, `positive.cpp`
checks exact nonnegative finite residuals in every fitted row and
channel. The resulting columns are lower bounds for the **physical**
resolvent at the certified upper discount, and therefore also at `μ`.
Neither the LP nor the finite truncation is used as a model for the
unknown true reproduction matrix.

Normalize the true left Perron vector by `ell·1=1`. Then
`πf=μ ell J G_μ f` for the crown/parent allocation. For the computed
lower newborn matrix `M≤K_μ`, the independent birth-column computations
give `ell≥ell M`. For any nonnegative vector `t` and a computed lower
occupation column `u`,

```
πf ≥ LO ell u ≥ LO min_i [u_i − t_i + (M t)_i].
```

Use the LP only to propose `t`, then evaluate the displayed dual bound
exactly and take the maximum with the worst-newborn result. In the
ungated non-first-body test, this improves the lower-side paid mean
`1.47502740e-7 → 1.69000953e-7` and the upper-side paid mean
`4.65412992e-7 → 5.28792328e-7`. Neither beats the previous size-12
occupation on its own. Restricting to charged crowns larger than 12
and excluding the high-host families `g,hC,S3` reserved below makes
the added fringe payments disjoint. They contribute `1.18312e-10`
and `1.15523e-10`. This identifies a genuine transfer, but also shows
where further marked-state enumeration would have little leverage.

## Whole-crown-cap-free high-root/young-arm payment

Use the stable-age identity on the actual physical newborn, with no
worst-newborn residual-killing penalty. Track only root degree `D`, its
number `k` of core children, and counts `m,n` of marked-for-observation
children in states `(D_child,k_child)=(2,0),(3,0)`. Mark a child for
observation at its real birth if fewer than `M` are observed; never
recruit older children. This observation flag is not a permanent
renewal birth mark. The exact loss/transition rates for these counts:

* An observed degree-two child promotes at `1/2` and exits the observed
  pair of states; its degree gain at `1+1/D` moves it to degree three.
* An observed degree-three child exits at `2+2/3+1/D`.
* The root promotes at `L/D`; the event changes `k` and recruits a
  new degree-two child if room remains. The root gains at
  `L+m/2+n/3+u`, where `0≤u≤(1+k−m−n)/2` is the joint live rate
  from the parent and unobserved children.

Pay the true selected `g,hC,S3` embeddings of the observed children.
For each embedding, all remaining observed senders use their own rates;
all other senders occur in the signed charge through the **same** `u`
as in the root gain. The paid sum is `max(0,a+b u)`, a lower bound
because the underlying selected rows are nonnegative. Each first-step
ratio is fractional linear between its possible hinge and the two
endpoints, so the exact minimum over those points proves the discounted
Bellman inequality for arbitrary adapted old arms. There is no bound
on the sizes or degrees of those arms, and no assumed independence
between the parent's degree path and a child's survival. The reward is
only paid for hosts with `D≥18`; its physical crown necessarily exceeds
12 there. The numerical comparison stops at `D>120` or `k>14`, and at
most eight simultaneous children are observed. The integral is over
all ages (no time cap), paying zero in omitted degree/arity phases. It is
cap-free in total crown geometry, not in the projected host domain.

| Positive occupation lower payment | Lower witness | Upper witness |
|---|---:|---:|
| `D≤65,k≤8,M=4`, degree-two only | `5.281301589e-7` | `2.871559990e-7` |
| Same limits, retaining degree-three phase | `5.967553797e-7` | `3.630872502e-7` |
| `D≤120,k≤14,M=8`, retaining degree-three phase | `6.202879563e-7` | `3.950853114e-7` |

The degree-three phase gives a larger gain than raising the projected
caps afterward. Unpriced older/non-star arms and other high-degree
motif families remain the real target. The existing size-12 bound only
pays the six higher families and never pays any `D≥18` host. The new
high-family payment and selected large-crown fringe payment are
mutually disjoint. The first-body channel is compared by **maximum** to
the older all-time first-body payment. Combining the three disjoint
non-first-body parts with the first-body payment gives mean numerator
floors `4.410155511e-6` and `3.921982939e-6`. Divide by the independently
checked [full stationary denominator ceilings](../denom/README.md)
`1.356833021311` and `1.367169995210`. This is the earlier mean-drift
identity at the original Neon endpoints, with no change to the potentials.

## Replay

From the project root:

```sh
g++ -O3 -std=c++17 round2/occupation-neon/fringe/positive.cpp -o <scratch-path>
OPENBLAS_NUM_THREADS=1 nice -n 10 taskset -c 7 ./round2/.venv/bin/python -Bu round2/occupation-neon/fringe/check.py b3f18_lo.bin b3f18_hi.bin
ulimit -v 10485760
OPENBLAS_NUM_THREADS=1 timeout 105s nice -n 10 taskset -c 7 ./round2/.venv/bin/python -Bu round2/occupation-neon/fringe/transfer.py <scratch-path> b3f18_lo.bin b3f18_hi.bin --verify-sha
```

The transfer replay passed both full pinned SHA-256 checks and all
residual checks on 79,211 marked rows, 344 charge geometries, and
27 channels, in 55.02 seconds at 473,756 KiB RSS. The independent physical
fixture ran in 1.20 seconds. The virtual-memory limit accommodates both
3.50-GB read-only maps; the input GMP proofs were reused, not rerun.
