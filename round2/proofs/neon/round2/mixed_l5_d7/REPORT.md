# Coarse L5: splitting BIG and moving the exact cutoff to seven

**Authors:** Gemini 4 Neon (lead) and GPT-6.2 Sol (prototype and checks). **Archive pin:** `c2498cc482bfef861644092d33f727cbbc93b280`. **Review and dependencies:** GPT exact tree fixtures and Gemini lead review cover this extension of the [mixed-row prototype](../mixed_l5_study/REPORT.md); the separate [Neon Claim 4](../../README.md) carries production review and release scope.

**Recommendation:** Pilot `D5=Dc=7` with two endpoint buckets in the GMP verifier at `D0=20` first. It had the narrowest float bracket and adds 122.69 MB per endpoint at `D0=20/24`. Keep `6/2` as the 25.98 MB fallback; defer `8/8`, which costs over 0.52 GB per endpoint and has only been counted. Measure the certified gain at `D0=20` before choosing a `D0=24` production run.

## Float experiment

Fixed `D0=12`, `KT=11`, `LF=6000`, `KF=100`, L4 `(10/8,9/7,10/8,8/6)`, 18 sweeps, A+B mask, one thread. Width is the exact decimal difference of the two final binary64 values. The values are float diagnostics, not certificates.

| Dc / buckets | λ lower | λ upper | Exact binary64 width | Wall lower/upper (s) | Peak lower/upper (KiB) | Last shift |
|---|---:|---:|---:|---:|---:|---:|
| 6 / 1 (prior study) | 0.56566461861668471 | 0.56578671991690521 | 0.00012210130022050069698025254183448851108551025390625 | 301.59/253.94 | 684516/684516 | 8.4e-11 |
| 6 / 2 | 0.56566661824833919 | 0.56578490667614756 | 0.0001182884278083662366043427027761936187744140625000 | 377.36/367.26 | 686384/686532 | 8.0e-11 |
| 7 / 1 | 0.56566704788809485 | 0.56578445411904543 | 0.00011740623095057856062339851632714271545410156250000 | 595.62/482.26 | 769548/769940 | 7.97e-11 |
| 7 / 2 | 0.56566877353270828 | 0.5657830484317693 | 0.00011427489906101762073831196175888180732727050781250 | 689.87/690.66 | 781332/780692 | 8.34e-11 |

At six, the split narrows the bracket by exactly `0.00000381287241213446037590983905829489231109619140625` (3.122712% of the unsplit width). It closes **20.981268%** of the prior `1.81727458e-5` gap to the optimistic A+B clip, leaving `1.43598734e-5`. This is the net split effect; the experiment does not isolate fiber spread from the two changes to gain-rate terms. A fiber-only control would retain the split state space and all split gain rates while feeding the **union** of the two fibers into both split sources. The D6 net gain is enough to justify measuring the split at seven.

With one bucket, raising the cutoff to seven narrows the six-cutoff unsplit reference by `4.69506927e-6` (3.845225%), closing 25.835773% of the optimistic-clip gap. Relative to the much smaller D6 split, it narrows the bracket by only `8.82196858e-7` for an extra 85.07 MB of payload.

At seven, splitting narrows the one-bucket bracket by `3.13133189e-6` (2.6671%) for an extra 11.65 MB. The winning layout narrows the D6 unsplit width by `7.82640116e-6` (6.409761%), closes **43.066696%** of the original clip gap, and narrows the historical no-L5 width by **14.469262%**. It is `4.01352875e-6` narrower than D6 split. Raw sweeps and timing are in `d{6,7}b{1,2}m{2,1}.log` for the three fitted configurations; timings include shared-host load. No width was measured at `D0=20/24`.

## Construction and Fraction audit

Only exact states may occupy an L5 interior. Large endpoints are either one `BIG` bucket or two buckets, `BIG_INNER` (`Dc ≤ D < D0`) and `TAIL` (`D ≥ D0`). The apex still has one state per bucket, shared across both physical `h` values. At a supported L4 row, only exact A/B receivers have their external credit capped at `1/Dc`; each deleted exact L5 leaf of degree `d` contributes its own `(1/d − 1/Dc)` times the marked L4 difference. The spectator bucket limits which L4 controls enter the min/max. The same weights multiply the before and after controls at each spectator, so take extrema of **differences**. Promotions at an exact L4 receiver add the newborn L5 value at rate `L/D` even with large spectators; any large promoter produces an unsupported interior, including overflow. `TAIL` stays `TAIL`, and contributes at most `1/D0` as an incoming rate. `BIG_INNER` contributes at most `1/Dc` as an incoming rate and can gain across the split only at `D0−1`. Its own crossing term lies between the smaller and larger of `0` and `(D0−1)·(q5(TAIL)−q5(BIG_INNER))`, because at the boundary `L+ρ ≤ L+(k+h)/2 ≤ D0−1`. This bound jointly includes noncrossing states and both signs of the jump. The split comparison measures the combined effect of changed fibers, tighter TAIL neighbor credit, and the new crossing slack. It does not isolate the fiber term alone.

[`check_trees.py`](check_trees.py) verifies both drift directions on 320 concrete trees per configuration. It uses signed L4 fixture products with exact Bernstein interpolation, rooted-symmetric signed L5 fixtures, independent forward/reverse edge-transfer enumeration, and independent newborn growth. It asserts both termwise and summed inequalities; it does not replay fitted C tables. Raw output: [`audit.log`](audit.log).

| `Dc` / buckets | Transfers | Strict source intervals | Crossings (negative/positive/zero jump) | Stays inner / stays TAIL | Births containing BIG(_INNER) / TAIL | Overflow births excluded |
|---|---:|---:|---:|---:|---:|---:|
| 6 / 2 | 1,865 | 619 | 322 (158/164/0) | 415 / 453 | 1,021 / 606 | 2,708 |
| 7 / 1 | 2,737 | 1,181 | — | 1,649 stays BIG | 2,097 / — | 2,708 |
| 7 / 2 | 2,737 | 880 | 500 (239/252/9) | 458 / 691 | 1,353 / 957 | 2,708 |

Planted overflow newborns failed ownership in all three runs. In both split runs the checker also rejected a missing crossing term, wrong TAIL fiber, zero TAIL incoming credit, and a missing TAIL-spectator newborn.

## Storage coverage and sizing

Coverage counts L4 storage row/receiver incidences with at least one external core edge; two symmetric receivers count twice. These counts are neither tree frequencies nor contributions to λ. Keeping the D5=6 all-small exclusion fixes the denominator for comparison. The last row gives D5=7 mixed-only coverage using its own all-small exclusion.

| Cutoff / exclude all slots below | A eligible / total | B eligible / total |
|---|---:|---:|
| 6 / 6 | 9,132,669 / 101,920,154 (8.960611%) | 3,002,870 / 98,749,300 (3.040903%) |
| 7 / 6 | 20,875,764 / 101,920,154 (20.482469%) | 8,498,815 / 98,749,300 (8.606456%) |
| 7 / 7 | 20,476,760 / 101,521,150 (20.169945%) | 7,933,800 / 98,184,285 (8.080519%) |

At 7/7, 63,713,090 A and 71,920,800 B incidences still have an oversized receiver; another 17,331,300 A and 18,329,685 B have a small receiver but a large interior spectator. Bucket splitting does not change coverage.

| `Dc` / buckets | L5 rows | Row bytes | Incremental payload bytes at `D0=12 / 20 / 24` |
|---|---:|---:|---:|
| 6 / 1 | 1,831,740 | 14,653,920 | 23,612,156 / 23,615,804 / 23,618,204 |
| 6 / 2 | 2,127,088 | 17,016,704 | 25,974,940 / 25,978,588 / 25,980,988 |
| 7 / 1 | 12,760,335 | 102,082,680 | 111,040,916 / 111,044,564 / 111,046,964 |
| 7 / 2 | 14,215,993 | 113,727,944 | 122,686,180 / 122,689,828 / 122,692,228 |
| 8 / 1 | 64,354,262 | 514,834,096 | 523,792,332 / 523,795,980 / 523,798,380 |
| 8 / 2 | 69,821,821 | 558,574,568 | 567,532,804 / 567,536,452 / 567,538,852 |

Row count is independent of `D0`; payload also includes layout, `2^19` cache slots, and `3*NG` fiber IDs. Existing L4 tables dominate total memory at production sizes. D5=8 is count-only, needs an S3 child cutoff of at least 7, and costs over 0.52 GB incremental per endpoint. Exact counts use the nine rooted shapes and multiset ranks for identical children; the C allocator reports the same D5=8 counts at all three `D0` values.

## Production GMP integration blueprint

1. **File contract.** Keep the v1 magic and 96-byte fixed little-endian header, set `version=2`; immediately after the fixed header insert `u32 Dc, u32 buckets, u32 receiver_mask, u32 reserved_zero, u64 nL5` (24 bytes, serialized individually). Keep all v1 arrays in order; append `nL5` signed little-endian int64 entries at the existing scale `S`. Enumerate shapes in `l5_shapes` order for size 5; within shape use `h=0` then `h=1`, ascending dense rank: ascending base global state IDs filtered by shape and h, then bucket IDs in BIG/TAIL order at endpoints, child groups sorted by shape and multiset colex rank. No per-row keys or trusted extrema. Recompute all ranks, sizes, cutoff compatibility, and `nL5`, checking arithmetic and exact file length before allocation. Reject unknown flags, out-of-range values, absent domains, and a nonzero reserved field. Keep v1 on its original path.
2. **One L5 row inequality.** Work on scaled values `I=S·q`. For every physical leaf position of the decoded shape with exact degree `d<Dc`, delete it and keep the marked neighbor. For eligible receivers enumerate the Cartesian product of compatible controls at every coarse spectator: inner `Dc≤D<D0` for `BIG_INNER`; INF and T1/T2/T3 for all allowed capped `k` (including overflow) for `TAIL`; union for unsplit BIG. Use both physical `h` at apex, only `h=1` elsewhere. Compute `ΔI4=I4(after)−I4(before)` with GMP integers (zero for an inactive after row), select min for lower/max for upper, and add `(Dc−d)·extΔI4/(Dc·d)`. Use a full-key cache only after computing the extrema from the current certificate; no float or independently extremized before/after. Add known degree reciprocals and `L` times each exact-slot gain difference, and `L/D` times its promotion difference. Add the direction-selected zero-or-upper-rate product for exact-slot slack, with upper rate `c_ext/2 + sum(1/Dc for BIG or BIG_INNER neighbors, 1/D0 for TAIL neighbors)`. For each `BIG_INNER` slot add the analogous corner from the `(D0−1)` crossing bound; `TAIL` has no own change. Subtract `(P/Q)·I5`; verify the sign in every row with existing GMP rational accumulators or an `mpz` numerator over a positive common denominator (e.g. `lcm(1,…,Dc,D0)·Q` for these terms). Check both directions via separate certificates.
3. **Couple to existing L4 checks.** In the four L4 verifier loops, project every old control tuple. Replace only eligible exact receiver external intervals `[0,c_ext/2]` with `[0,c_ext/Dc]`; use exact rational scaling without changing the separate `ext_half` rates for genuinely tail receivers. Add the supported q5 newborn to the **promotion delta** for each exact promoter, independently of receiver mask, with unchanged `L/D`. Keep the same newborn value at each control for a tail spectator, including INF and overflow, so interpolation gives weight one. Leave tail-promoter and overflow bud equations unchanged because L5 support excludes any large promoter. Extend the exact concrete-tree final check to enumerate L5 motifs and newborns independently; reject wrong bucket, inverted source direction, missing crossing, and missing newborn with negative controls. Fit at fixed rational λ and replay after quantization; float convergence is not a substitute for row checks.
4. **Exact interior-anchored comparison.** Charge each supported L5 motif to its unique topmost interior (interiors induce a connected subtree). For a charge of degree `D<Dc`, grow an ambient tree rooted there, stopping expansion at any large endpoint. The root has at most `D` branches, every other vertex at most `Dc−2` forward branches. Every charged motif embeds in this tree, since expansion through a large vertex would make it an unsupported interior. Hence `B(D)=[x^4](1+F)^D = D·binom(D−1+4(Dc−2),3)/4`, with `F=x(1+F)^(Dc−2)`, bounds the number of charges by counting connected subsets containing the root. This gives uniform maxima 1,425 at six and 3,450 at seven. Scan L5 integer extrema including zero: `m=min(0,min I5)`, `M=max(0,max I5)`. In `verify_comparison_bounds`, for each exact charge state add `B(D)*m` to `loz` and `B(D)*M` to `hiz` in `mpz_t` **before** multiplying `loz` by four and applying the existing zero-sum recharging. Other states get zero. Then apply the existing global min/max and require upper `4S*c0>0`, both modes `S*C0>0` and `S*Phi(T3)>0`. `T3` cannot contain an L5 motif. If needed, scan extrema by the exact charge state instead of globally.

## Reproduction

From `round2/mixed_l5_d7/`:

```sh
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 TMPDIR="${TMPDIR:-$PWD}"
ulimit -c 0
python3 -B build_sources.py
gcc -O3 -march=native -ffp-contract=off -fopenmp diag_slack.c -lm -o "$TMPDIR/neon-mixed-d7"
export SOLVER="$TMPDIR/neon-mixed-d7"
bash run.sh 6 2
bash run.sh 7 1
bash run.sh 7 2
python3 -B build_checker.py
(python3 -B check_trees.py --d 6 --buckets 2 &&
 python3 -B check_trees.py --d 7 --buckets 1 &&
 python3 -B check_trees.py --d 7 --buckets 2) > audit.log 2>&1
python3 -B size_and_coverage.py
for b in 1 2; do
  for d0 in 12 20 24; do
    bash run.sh 8 "$b" --count --D0 "$d0" --KT "$((d0-1))" --KS3 7
  done
done
python3 -B summarize.py d6b2 d7b1 d7b2
rm "$TMPDIR/neon-mixed-d7"
```
