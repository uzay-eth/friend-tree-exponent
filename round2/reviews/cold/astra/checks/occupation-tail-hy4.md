CORRECT

# Independent referee review: reciprocal-degree occupation tail

Scope. This reviews `continuation.md` (`occupation-tail/continuation.md`),
`solver.cpp` (`occupation-tail/solver.cpp`), `physical.py` (`occupation-tail/physical.py`) and
`replay.py` (`occupation-tail/replay.py`), conditional on the separately reviewed raw certificate
(`reserve/README.md` (`motif-cert/reserve/README.md`), `reserve.hpp` (`motif-cert/reserve/reserve.hpp`)),
the higher-row sources (`source.md` (`motif-high/occupation/source.md`), `credit.hpp` (`motif-high/occupation/credit.hpp`))
and the first calculation (`inside.md` (`occupation/inside.md`)). No author file was changed. The raw tables
`mc-upper.bin` (`outputs/mc-upper.bin`) and `mc-lower.bin` (`outputs/mc-lower.bin`) were opened read-only only.

## 1. Source disjointness

The five added families are gated by distinct mask bits in `credit()` (`solver.cpp:165-196`) and call different
frozen APIs: singleton unused row (`singleton_floor`), pair unused row (`edge_unused_floor`), internal-edge donor
on a triple (`edge_donor_floor`), tightened intrinsic row of a triple/quartet (`row`), and an uncovered donor of a
triple booked on a four-set (`triple_donors_on_four` / `triple_donor`). Bit 32 selects an *alternative*
implementation of the fifth family (`if (!(mask_&32)) … else …`), so the two are mutually exclusive, not additive.

These are disjoint from the baseline credits by the prerequisite proofs and by the explicit statement in
`continuation.md:11`: the baseline five-vertex credit uses four-set external donors and the baseline edge credit
uses the singleton donor gap; neither is among the new sources. Consistent with
`reserve/README.md:53` (`motif-cert/reserve/README.md`) (both singleton credits retained together) and
`source.md:49` (four-set donor, triple donor and high-degree singleton donor have
distinct physical origins). I found no summand claimed twice.

## 2. Physical lower affine extension at all degrees

I re-derived the tail rate bounds. For a tail receiver at knot d with physical degree D≥d, other coordinates fixed,
H = t·H_d+(1−t)·H_inf with t=d/D. A leaf gain gives t′=d/(D+1) and hence

  ΔH = (t′−t)(H_d−H_inf) = d/(D(D+1))·(H_inf−H_d) = t·[1/(D+1)]·(H_inf−H_d).

The physical rate a satisfies a≥D−c≥d−c and a≤D−c+c/2≤D+1, so a/(D+1)∈[(d−c)/(d+1),1]; the solver's
gain range `{d−c, d+1, d+1, old(infinity)}` (`solver.cpp:241-243`) encloses it. With the corner endpoint chosen by
sign(H_inf−H_d), the corner term is min over that interval, and the physical gain equals t·(a/(D+1))·(H_inf−H_d) ≥
t·(corner term). The promotion interval `[1−c/d,1]` (`solver.cpp:252-253`) is valid since b=1−c/D≥1−c/d; the
infinity corner carries zero gain. The same argument covers an exact receiver whose neighbour is interpolated: the
rate interval is fixed and contains the actual rate at every corner, so the min endpoint bounds the multiaffine jump.

Projection composes: (d/D)·(d′/d)=d′/D, implemented in `Axis::project` (`solver.cpp:50-61`) and the tensor product in
`Grid::interpolate` (`solver.cpp:87-113`), with one downward rounding of a nonnegative contribution. This is the
correct affine lower extension for every D≥d, including D unbounded.

## 3. Chronological deletion / old-state semantics

Connected sets are ordered oldest-first; a promotion at member i creates the unique enlarged set by adjoining the
newborn in state (2,1) with `parent[n]=i` (`solver.cpp:208-210, 248-251`), and deletion of the youngest member of that
enlarged set recovers the old set — exactly the bijection in `inside.md:29`. `quartet_donors` (`solver.cpp:132-147`)
enumerates every leaf-deletion of a quartet, renumbers the retained triple in chronological order, and matches the
frozen rooted shape. `old()` (`solver.cpp:224-228`) asserts every in-table successor has a strictly larger label
index; rows are visited in reverse (`solver.cpp:211`), so successors are already solved. Out-of-cutoff successors
return value 0, which is the favourable sign. Nonnested finer/coarser cutoffs are checked in the constructor
(`solver.cpp:299-300`). My aggregated tree test (below) exercises these semantics end to end.

## 4. Row-min signs

For each range the code sets the rate to `low` when `destination>=z` and `high` otherwise (`solver.cpp:262`, and
`Range::term` at `solver.cpp:119-122` uses the same sign of `destination−v` in the exact residual). That is exactly
min over a nonnegative interval of r·(V′−V), the correct lower bound for (3) of `inside.md`. The acceptance is a
final exact integer residual (`solver.cpp:271-278`), so the floating proposal cannot leak into the result.

## 5. Implementation and exact arithmetic

Signed 128-bit intermediates, `floor_ratio` and `narrow` are bounded as argued in `source.md:77-79` and
`inside.md:51`; I found no overflow path in the interpolation masses, the discount term −4096·num·v/den (num≤10^6),
or the per-range terms. I compiled `solver.cpp` with `-fsanitize=undefined -fno-sanitize-recover=all` and ran it on a
small grid on one CPU without any UB report. A test-only reconstruction of the single-vertex row from the dump
shows that flipping the endpoint rule would wrongly accept an oversized score (e.g. upper row (31,9,0): 5639→5640),
confirming the sign rule is load-bearing.

## 6. Both fractions and outward rounding

Exact `Fraction` recomputation from the unchanged baseline endpoints and ceilings reproduces both claims and both
outward ticks (the lower floor and upper ceiling of `end·10^10`):

```
lower  end = 565662790716854568223/1000008252270910000000   equal=True  tick 5656581227
upper  end = 1698046074520368376961/3001152517940340000000   equal=True  tick 5657979941
```

with v_lower=15,289,028 / y_lower=2,654,889,165,321 and v_upper=17,458,365 / y_upper=2,655,887,184,018, added on
the lower side and subtracted on the upper side, matching `replay.py:59-63`. The snapshot digests match the frozen
certificate (`sha256sum` of both raw tables equals the constants in `reserve.hpp:229-231`).

## 7. Independent exact test with planted errors

`occupation_tail_hy4.py` is my own exact-Fraction test on a grid **different** from the
author's `physical.py` grid, D=[[32],[18,14],[12,9,7],[10,7,5,4]], C=[[9],[6,5],[4,4,3],[4,3,2,2]]. It re-derives
the interpolation and the aggregated identity source ≥ u·T − G + v·B from scratch, enumerates 269 real small trees
(reachable states plus decorated cores up to degree 10^9+3), and asserts the inequality at every one. Results:

```
[upper] aggregated inequality PASS on 269 real trees; tail identity 30000/30000 exact; flat 30000/30000 rejected, shift 30000/30000 rejected; newborn fault rejected on 252 trees; INDEPENDENT REFEREE TEST: PASS
[lower] aggregated inequality PASS on 269 real trees; tail identity 30004/30004 exact; flat 30004/30004 rejected, shift 30004/30004 rejected; newborn fault rejected on 191 trees; INDEPENDENT REFEREE TEST: PASS
```

Two specific injected errors are demonstrated to be rejected: (i) flattening the interpolation weight (t=1) and
shifting it to (d−1)/D both break the exact tail identity on every sampled tail row; (ii) inflating the newborn score
breaks the aggregated physical inequality on many trees. The planted wrong-endpoint rule is rejected at the row level
as noted in §5. Commands (from the repository root, one CPU):

```
taskset -c 0 python3 -B round2/checks/occupation_tail_hy4.py upper
taskset -c 0 python3 -B round2/checks/occupation_tail_hy4.py lower
```

## 8. Limitations

I did not rerun the full production tables or the author's full UBSan replay (76 s CPU, ~1.4 GB RSS per side); per
instructions those long runs were generated by the author. I relied on the recorded production outputs
`occupation-tail-affine-upper.txt` (`outputs/occupation-tail-affine-upper.txt`) and
`occupation-tail-affine-lower.txt` (`outputs/occupation-tail-affine-lower.txt`) for the row counts and the two v
values, and verified all downstream arithmetic exactly. The prerequisite headers, the raw certificate and
`inside.md` are taken as separately reviewed. No gap was found in the material under review.
