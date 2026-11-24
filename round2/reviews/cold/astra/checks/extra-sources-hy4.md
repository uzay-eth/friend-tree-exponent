CORRECT

# Cold audit: eligibility of the additional source APIs

Scope. This audits the *inputs* my tail review conditioned on:
`reserve/README.md` (`motif-cert/reserve/README.md`), `reserve.hpp` (`motif-cert/reserve/reserve.hpp`),
`source.md` (`motif-high/occupation/source.md`) and `credit.hpp` (`motif-high/occupation/credit.hpp`).
It is not a new tail-source derivation and creates no generic tail-sources folder; the affine calculation that
consumes these credits was confirmed separately. Only the two frozen snapshots were read; no author file was changed.

## 1. Frozen positive allocations bound the maximized residual

`certificate.json` (`motif-cert/certificate.json`) gives singleton positive allocations 154,457 (upper) and
95,343 (lower); `reserve.hpp` uses exactly these as `fcap` and asserts the computed floor is nonnegative. My probe
returns 79 singleton floors in `[155337,156965]` (upper) and `[2226,72369]` (lower), all nonnegative, and every
`edge_unused_floor` equals `max(0,−residual_ceiling)`. The exhaustive census log
(`mc-reserve-all.log` (`outputs/mc-reserve-all.log`)) records maximum physical ceilings 520/10,646 inside the
archived caps 524/10,651, with all 979,299 singletons and 6,714,612 ordered pairs checked. The higher-row intrinsic
credits are likewise bounded pointwise by the frozen positive corner parts (`source.md` §Compatibility): capped
intrinsic totals at most 403/S upper and 7,184/S per core lower, never exceeding the archived charge c_M.

## 2. Paid U−ψ and −Rmax cannot overlap

`compute_unused` (`reserve.hpp:208-215`) forms the residual with `maximum(z,mode)` at every slot, i.e. with ψ replaced
by its **supremum** U; `edge_donor_floor` separately charges `U−ψ(1/D_y)` (`known_fy_gap`). Thus the row credit
`floor(max(0,−R_max))` and the donor credit `floor(U−ψ)` are disjoint pieces of the decomposition
R_actual = R_max − Σ(U−ψ), and their sum against the signed row is at most max(0,R_max). My negative control
recomputes a hypothetical row from the *actual* donor gap (`ce−df`): on the sampled queries this would double-claim
476,182,202 units (upper) and 6,209,623 (lower), so using the true donor in the row is exactly the forbidden overlap;
the shipped code never does it. The shipped held-out check already rejects `pool-as-exact`, `double-edge` and
`double-singleton` defects (`mc-reserve-verified.log` (`outputs/mc-reserve-verified.log`): 2,594/190/237 upper and
4,818/474/238 lower planted rejections), confirming distinct physical incidences.

## 3. Discrete FY grid vs continuum, exceptional knots

The frozen header has Q_K=10080 and 15 knots: 0, 1/72, 1/48, 1/32, 1/24, 1/18, 1/14, 1/10, 1/7, 1/5, 1/4, 1/3,
**5/12**, **11/24**, 1/2. Exactly 5/12 and 11/24 are non-real (Q_K not divisible by the numerator), matching
`reserve/README.md:41`. The first segment is exact at 1/72 (`D0·knot_num[1]=Q_K`, also asserted by the constructor).
My checker reproduces the constructor's interval scan: for all 59 degrees in 2..72 that are not themselves knots,
both bracketing knots are real, so continuous and physical optimization agree on every affine segment and a maximum
over feasible knots computes the physical supremum. Pooled donors use the documented endpoint/knot set and clamp
to zero where needed.

## 4. Degree-72 true seam

`difference()` uses gain denominator `D+1` at D=D0=72 (pre-gain tail mass 1), not D. A patched-header control
(replacing `x.D+1` by `x.D` in a temporary copy) leaves the D=71 edge ceiling unchanged (75 upper, −890 lower) but
changes the D=72 ceiling from −7,475,439 to −20,483,660 (upper) and from −28,148,921 to −15,376,892 (lower).
The +1 is therefore load-bearing exactly at the true seam; the shipped code uses it.

## 5. High-role permutations

`Shape::from_parents` enumerates input permutations to match each frozen rooted shape, taking the first matching
order for symmetric children (`credit.hpp:24-40`). My probe builds canonical and role-swapped shapes at identical
physical states: swapping the two symmetric children of the cherry C, the two symmetric children of the Y, and any
two of the trident T gives identical `row` values on both witnesses. Role maps are thus a canonical bookkeeping
choice, not a source of spurious credit, and `validate_shape` rejects non-permutation role arrays.

## 6. Physical gain and integer floors; representative oracle

My probe recomputes the signed physical gain directly from the frozen `mh::Snapshot` tables (gain = table after a
leaf minus table before, using the level's projection at the seam), independent of `Source`'s own row engine, then
compares the API credit with `floor(J(d,1/7))`, J(d,p)=max(d,0)/2−dp, and with the pooled floor
`floor(min_{p∈[0,1/7]} J)`. On representative active triple states all 8 exact and 2/7 pooled cases match on
upper/lower, including negative gains where the exact floor is `floor(−d/7)` and the pooled floor is zero. The
integer floor is `mpz_fdiv_q` (downward); a round-toward-zero control differs whenever J is nonintegral. Fast/GMP
agreement is additionally enforced by `prepare` (`hi==lo`, else exact recomputation with lo≤exact≤hi).

## 7. Verdict and limitations

Commands, from the repository root on one CPU:

```
taskset -c 0 python3 -B round2/checks/extra_sources_hy4.py upper
taskset -c 0 python3 -B round2/checks/extra_sources_hy4.py lower
```

Both print `EXTRA SOURCES COLD AUDIT: PASS` with the evidence above. I did not rerun the exhaustive census, the
sweeps or the held-out tree measurements; those are cited from the frozen logs. The raw signed inequality and the
affine calculation remain separately refereed. No eligibility defect was found.
