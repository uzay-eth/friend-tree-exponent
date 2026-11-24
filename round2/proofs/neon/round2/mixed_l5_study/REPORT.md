# Mixed-row L4 slack and coarse-endpoint L5

**Authors:** Gemini 4 Neon (lead) and GPT-6.2 Sol (prototype and checks). **Archive pin:** `c2498cc482bfef861644092d33f727cbbc93b280`. **Review and dependencies:** GPT exact fixtures and Gemini lead review cover diagnostics based on the [pure-inner Level-5 prototype](../level5_study/REPORT.md). See [auxiliary evidence](../../README.md#auxiliary-evidence).

## Reproduction

From `round2/mixed_l5_study/`:

```sh
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OPENBLAS_NUM_THREADS=1 TMPDIR="$PWD"
python3 -B make_diag.py
gcc -O3 -march=native -ffp-contract=off -fopenmp diag_slack.c -lm -o "${TMPDIR:-.}/neon-coarse"
SOLVER="${TMPDIR:-.}/neon-coarse" bash run.sh grid
SOLVER="${TMPDIR:-.}/neon-coarse" bash run.sh coupled
python3 -B summarize.py
python3 -B summarize.py --report
python3 -B check_coarse.py
python3 -B count_variants.py
/usr/bin/time -f 'TIME wall=%e user=%U sys=%S peak_KiB=%M exit=%x' \
    "${TMPDIR:-.}/neon-coarse" --threads 1 --l5 6 --coarse 6 --l5-check
"${TMPDIR:-.}/neon-coarse" --threads 1 --l5 6 --coarse 6 --count
"${TMPDIR:-.}/neon-coarse" --threads 1 --l5 6 --coarse 6 --count \
    --D0 20 --KT 19 --LF 16000 --KF 200 \
    --DP4 12 --KP4 9 --DYT 11 --KYT 8 --DYB 12 --KYB 9 --DS3 10 --KS3 7
```

The solver command in `run.sh` uses `--D0 12 --KT 11 --LF 6000 --KF 100 --DP4 10 --KP4 8 --DYT 9 --KYT 7 --DYB 10 --KYB 8 --DS3 8 --KS3 6 --nit 18 --threads 1 --l5 6` in modes 2 (lower) and 1 (upper). The recorded grid was launched first with the initial attribution-only source, built using the same flags to `${TMPDIR:-.}/neon-mixed-diag`. The later source retains those equations behind `--coarse 0` (default); it adds the coupled experiment `--coarse 6 --coarse-mask 3`. The queue runs A/B/C before the coupled pair and remaining grid. Each solve is single-threaded; the queue and verification commands were also pinned to CPU 9 during the run.

## Fine-grained attribution

Each ablation reruns the coupled **pure-inner** D5=6 solver and changes only the remaining mixed L4 external-credit interval. A = interior L4 receivers; B = inner-degree L4 endpoint receivers; C = known-`k` tail L4 endpoint receivers (external rate extrema in T2/T3). T1 has no such external-rate dependence. The overflow hull is unchanged. `zero` replaces the selected interval by `[0,0]`; `clip` replaces it by `[0,c_ext/6]`. Both are diagnostic interventions rather than valid rate bounds without a transfer. Compare against `reference` and the joint intervention `zeroABC`; individual changes need not add. The earlier float reference width before pure-inner L5 was `0.00013360681999630142`; the earlier coupled pure-inner and fully suppressed mixed widths were `0.00012800250560285331` and `0.000085083282862830067`, at 16 sweeps. The local grid is 18 sweeps. Widths from `summarize.py` are exact decimal representations of the difference between the two round-tripped binary64 endpoints; this does not make the endpoints rigorous.

**Inner endpoints contribute the largest measured mixed-rate sensitivity, followed by interiors.** Using the earlier pure-inner L5 run as reference, independently suppressing B reduces width by 18.10%, A by 13.33%, and C by just 0.73%. The coupled experiment covers eligible A+B receivers and leaves C conservative.

<!-- BEGIN RESULTS -->
| Case | λ lo | λ hi | Exact binary64 width | Wall lo/hi (s) | Peak lo/hi (KiB) | Last shift |
|:--|--:|--:|--:|--:|--:|--:|
| reference | 0.56566050377360144 | 0.56578850684117155 | 0.0001280030675701038234137740801088511943817138671875 | 225.16/225.64 | 674276/674280 | 7.2e-11 |
| zeroA | 0.56566903239809441 | 0.56577996836772781 | 0.0001109359696334077938217888004146516323089599609375 | 231.42/242.54 | 674280/674152 | 5.6e-11 |
| zeroB | 0.56567468670125687 | 0.56577951792288395 | 0.00010483122162707747548893166822381317615509033203125 | 198.80/194.03 | 674280/674276 | 4.5e-11 |
| zeroC | 0.5656610973891516 | 0.56578815942505201 | 0.000127062035900404879384950618259608745574951171875 | 261.47/234.34 | 674280/674280 | 7.1e-11 |
| clipA | 0.56566524850087485 | 0.56578433368773517 | 0.00011908518686032376621142248040996491909027099609375 | 220.40/240.61 | 674284/674280 | 6.3e-11 |
| clipB | 0.56566936739544915 | 0.56578279936192899 | 0.00011343196647983777580748210311867296695709228515625 | 222.66/201.93 | 674280/674280 | 5.5e-11 |
| clipAB | 0.56567436682391092 | 0.56577829537831392 | 0.00010392855440299886282673469395376741886138916015625 | 258.60/231.69 | 674280/674280 | 4.5e-11 |
| clipC | 0.56566089296696176 | 0.56578827924003283 | 0.0001273862730710728641270179650746285915374755859375 | 186.62/196.56 | 674284/674280 | 7.1e-11 |
| zeroABC | 0.56568463870540375 | 0.56576972228192535 | 0.000085083576521594039832052658312022686004638671875 | 229.21/194.28 | 674284/674280 | 2.6e-11 |
| coupled | 0.56566461861668471 | 0.56578671991690521 | 0.00012210130022050069698025254183448851108551025390625 | 301.59/253.94 | 684516/684516 | 8.4e-11 |

`zeroA` narrows the reference width by 13.333%.

`zeroB` narrows the reference width by 18.103%.

`zeroC` narrows the reference width by 0.735%.

`clipA` narrows the reference width by 6.967%.

`clipB` narrows the reference width by 11.383%.

`clipAB` narrows the reference width by 18.808%.

`clipC` narrows the reference width by 0.482%.

`zeroABC` narrows the reference width by 33.530%.

**Coupled endpoint-only BIG: 8.611% width reduction versus the earlier no-L5 baseline.** The absolute reduction is 2.05× that of the earlier pure-inner L5 run. The measured width is at or above `1.00e-4` at D0=12. This result does not establish a D0=20/22 width.

**Coupled minus optimistic A+B width: 0.0000181727458.** The optimistic A+B run narrows the no-L5 baseline by 22.213%.
<!-- END RESULTS -->

### Why the real transfer trails the optimistic A+B clip

1. **Coverage:** `clipAB` reduces the external rate on *every* mixed A/B row. Coarse-BIG requires every L4 interior and the receiver to have `D<6`. `count_variants.py` counts storage-row/receiver incidences with at least one external core edge, excluding the pure-inner D5=6 block. These are exact layout counts, not occurrence weights or λ sensitivities:

   | Receiver | Eligible in coarse-BIG | Receiver is too large | Small receiver, large spectator interior | Total |
   |:--|--:|--:|--:|--:|
   | A | 9,132,669 (8.96%) | 79,915,955 | 12,871,530 | 101,920,154 |
   | B | 3,002,870 (3.04%) | 80,910,900 | 14,835,530 | 98,749,300 |

   The receiver cutoff leaves `6 <= D_u < D_cap` untouched at interiors (`D_cap=8,9,10` by shape), and `6 <= D_u < D0` untouched at inner endpoints. A small receiver can also be ineligible because another interior is too large. Rows can have multiple receivers; the three breakdown columns partition incidences, not motifs.
2. **Fiber spread:** where the credit is transferred, an exact four-state difference becomes the min/max over BIG spectator controls. Signed differences and multiple independent control axes can make that interval wide, especially when BIG includes both tail and inner degrees. The exact checker exercises strictly interior cases of this enclosure (647), but uses fixture functions and does not quantify its impact on converged λ. A two-bucket split attacks this approximation directly.
3. **L5 generator slack:** the transfer retains its contribution in L5, where each exact receiver still sees `[0,c_ext/2]` for its own external neighbors and an additional interval for any BIG neighbor. These rates multiply `Δ q5`, not `Δ q4`. The oracle creates no L5 correction for the extra clipped rows and has neither of these new terms. No individual λ share is identified for fiber spread or L5 generator slack by the existing grid, and the three effects cannot be added as width losses.

## Mathematical construction and exact checks

The C self-check visited all 1,831,740 coarse rows in both modes; zero source gave exact zero response, and a varying P4 source gave 1,203,500 nonzero rows in each mode. Time was 14.40 s wall, 7.06 s user, 0.13 s sys, 230,656 KiB peak, exit 0 (includes time-sharing with the attribution solve on the same CPU).

<!-- BEGIN DERIVATION -->
## Exact endpoint-only BIG construction

Write `h` for presence of a core parent, `k` for number of core children, `L` for number of leaf children, and `D=L+k+h`. Use `D_c=D_5=6`. Keep exact `(L,k,h)` at degree below six. Permit a `BIG` vertex (physical degree at least six) only at *undirected* L5 endpoints. Rooted shape carries the apex and all edge directions. BIG at the apex deliberately forgets `h` as well as `D,L,k`; at a nonapex endpoint `h=1` is implicit. Potentials are zero on unsupported motifs. At L4 an external credit is eligible for clipping if its receiver is exact, all L4 interiors are exact, and the receiver belongs to the configured role mask. Large L4 endpoints *other than the receiver* are allowed, whether inner or tail. The C receiver is never clipped by this design: attaching the deleted leaf would make it a BIG interior. A/B attribution applies to all mixed rows; the coupled design only clips their supported subset. Large L4 interiors are outside that subset even if the receiver is small.

### Source and coarse fiber

For an eligible incidence `(M,u,z)` use `r=1/D_z`, `r'=min(r,1/6)`, `x=r-r'`. Then `r Δ_u q4 = r' Δ_u q4 + x Δ_u q4` for either sign. For ineligible incidences retain `r` intact. For a receiving physical L5 motif `N`, each leaf `z` with `D_z<6` can give one pair `M=N-z`, `u=neighbor(z)`; if the eligibility predicates hold, that is exactly the inverse of the original incidence mapping. This includes an apex leaf when present, and includes each distinct physical leaf even if labels or rooted subtrees coincide. No automorphism factors enter the sum. An L4 domain exit at the increment uses zero for the after-value; never assume `Δ_u q4` has one sign. In this configuration all old L4 preimages of supported L5 deletions are active under their original degree/child-count cutoffs; the implementation also zero-extends the lookup.

For a given marked coarse preimage, hold all its exact coordinates fixed. For each BIG endpoint, take every compatible old inner control with `6 <= D < D0`, plus INF and T1/T2/T3 for all allowed capped `k` including overflow. For a BIG apex include both `h` values. For a descendant use `h=1`. The physical L4 before and after use the **same** probability weights at these endpoints since the gain is at the exact receiver. Thus the difference itself is a convex combination of control-point differences, not an independent subtraction of two extrema. The minimum and maximum over the Cartesian product enclose every *physical* fiber member. The bounds can be wider than the physical fiber extrema because the Bernstein controls need not individually occur on a physical tree. Multiply each by `x >= 0` and sum marked deletions. Correlation between the deletions cannot invalidate this sum; it can make it loose. The source cache stores extrema by the full ordered four-vertex pattern and marked receiver, is reset every sweep, and checks the full key on collisions.

### Existing motifs and promotion ownership

For an existing supported L5 motif, a BIG vertex stays BIG on gain and on promotion, with own change exactly zero. It contributes an incoming `1/D` at its unique neighbor in `[0,1/6]`. An exact vertex sees its `L`, exact reciprocal degrees from exact motif neighbors, and an additional interval `[0, c_ext/2 + (# BIG motif neighbors)/6]` for the remaining gain rate. The additive rate uncertainty enters linearly with coefficient `q5(after gain)-q5(now)`; take the sign-selected corner in the requested direction. This interval extremization is exact over the box, not necessarily attainable across neighboring rows of a physical tree. A gain at degree five enters BIG if the motif slot is an endpoint, or zero if it is an interior. Promotion preserves degree; it changes exact `(L,k,h)` to `(L-1,k+1,h)` at exact rate `L/D`. Unsupported motifs with a large interior cannot later enter support as degrees are nondecreasing and shape is fixed.

Every newborn L5 motif contains the promoted leaf and has exactly one old preimage after deleting it. Its old promoter must be an L5 interior: a connected old four-set uses at least one motif edge at the promoter and the newborn adds another. Therefore a BIG promoter creates **no supported newborn**, including any `k>KT` overflow. A promoter with `D_c<=D<D0` has an exact physical bud rate, but its newborn term is zero under this support rule. A tail endpoint can be BIG *elsewhere* in a surviving newborn motif; the bud contribution is constant over its Bernstein controls (including INF and overflow), so injecting it into each relevant L4 row for the exact promoter yields exactly `L/D` times that newborn value after interpolation. No `b3_tail_const/var` bud argument needs modification.

For a broader projection that allows a BIG promoter to become interior, degree monotonicity alone does **not** justify omitting its birth. A known-`k` tail birth would have to contribute `B(1-c/D)` with the actual coarse destination `B`. In particular, if the same `B` is added to the T1/T2/T3 bud destinations, their aggregate accounts for `(1-s^3)B(1-c/D)`, leaving the INF contribution `s^3 B(1-c/D)` to account for separately; INF also does not carry `k,h`. In overflow the original bud change on the *old* motif vanishes because the capped label is unchanged. This says nothing about newborns. For the unchanged overflow L4 equations, support must forbid those newborns, as it does here; otherwise supply new bounds for their generator contributions. A separate `TAIL` label does not by itself solve either problem.

### Direction, comparison, and limits

On each tree the original potential drift equals the sum of clipped L4 rows, the exact rerouted L5 sources, the existing L5 changes, and all newborn terms charged to L4. Replace each source and L5 gain-rate expression by a lower or upper bound in the same direction as the L4 rows. Convex interpolation of spectator endpoints preserves the inequalities. This is a pointwise drift reduction and does not assert an eigenvalue bound from a floating fixed-point iteration. An exact certificate needs rational row checks, including the added terms, and an appropriate positive comparison constant.

For comparison charge a supported motif to its topmost *interior* (unique: the interior vertices induce a connected subtree). Start from that vertex and grow by four adjacent vertices along the motif. Every attachment parent is interior: a motif endpoint can never provide a previously unadded neighbor once it has been attached. All interior degrees are at most five. Thus there are at most `4! 5^4 = 15000` candidates per charged core vertex; multiplying this count by the signed extrema of `q5` including zero provides uniform comparison corrections. The bound can be tightened by a shape-by-shape count per exact charge state. A coarse apex may have unbounded degree, so the naive apex-rooted count is invalid. No five-core motif exists at the initial `T3`.

The checker uses deterministic signed test potentials (multilinearly embedded L4; rooted symmetry for L5), independent forward and reverse transfer enumeration, and independent growth from each newborn. It measures the tree inequality for the modified terms, not a GMP replay of fitted float tables. Exact output: `Exact coarse tree inequalities PASS: big_birth_excluded=11832, big_stays_big=1257, birth_with_big=1566, births=3325, credits=45694, gain_crosses_dc=705, l5_rates=6053, mixed_credits=914, multiple_deletions=575, overflow_birth_excluded=2708, strict_source=647, transfers=1896, trees=320`.
<!-- END DERIVATION -->

## Entry counts and memory

The dense L5 layout enumerates all nine rooted shapes, with multiset ranking for identical rooted child shapes. The row count is independent of `D0` as long as `D0>6`. The degree cutoff at an L5 *interior* prevents support for a BIG promoter or BIG receiver, but does not require any full-resolution endpoint indexing in L5. The C allocated payload includes row values, the 524,288-slot direct-mapped source cache, the layout struct, and lists of old Bernstein/inner controls. Payload excludes allocator overhead. Values below are per endpoint; both directions run sequentially.

| D0 | L4 configuration | Base entries | L5 entries | L5 row bytes | Total L5 payload bytes | Whole solver table payload estimate |
|--:|:--|--:|--:|--:|--:|--:|
| 12 | 10/8, 9/7, 10/8, 8/6; LF/KF 6000/100 | 82,520,610 | 1,831,740 | 14,653,920 | 23,612,148 | 693,474,660 B (0.646 GiB) |
| 20 | 12/9, 11/8, 12/9, 10/7; LF/KF 16000/200 | 1,398,838,433 | 1,831,740 | 14,653,920 | 23,615,796 | 11,265,782,492 B (10.492 GiB) |

Whole solver estimate = eight bytes per base entry, a second copy of `f`, plus total L5 payload; excludes all other scratch, rate arrays and overhead. Both incremental budgets pass. The unmodified base production profile at D0=20 already exceeds 10 GiB; the full solver remains above it. A D0=20 count is not a solve or an extrapolation of the coupled width.

### Exact variant sizing (no additional fitted runs)

`python3 -B count_variants.py` uses the rooted-shape multiset recurrence. With `a` motif children at a vertex with fixed `h`, the exact-state factor is `sum_{d=2}^{Dc-1} max(0,d-h-a+1)`. Sum apex factors over `h=0,1` and add one state per BIG bucket only at undirected endpoints (a shared apex state per bucket). An equal-child-shape group of multiplicity `m` and child-state count `n` contributes `choose(n+m-1,m)`. The extra byte estimates keep the same cache and fixed-capacity layout; two fiber lists can partition the existing backing allocation with small fixed metadata overhead if stored separately. The table reports the packed layout estimate with unchanged metadata footprint, excluding allocator overhead.

| Exact boundary / clipped cap | BIG buckets | L5 rows | Row bytes | Projected incremental bytes D0=12 / D0=20 |
|:--|--:|--:|--:|--:|
| 6 / 6 (implemented) | 1 | 1,831,740 | 14,653,920 | 23,612,148 / 23,615,796 |
| 7 / 7 | 1 | 12,760,335 | 102,082,680 | 111,040,908 / 111,044,556 |
| 6 / 6 | 2 | 2,127,088 | 17,016,704 | 25,974,932 / 25,978,580 |
| 7 / 7 | 2 | 14,215,993 | 113,727,944 | 122,686,172 / 122,689,820 |

The D5=7 row raises `Dc` with D5 so all nonzero donor weights remain exactly identified; simply keeping one BIG class at `Dc=6` would mix `D_z=6` positive weights with `D_z>=7` zero weights. For the two-bucket 6/6 variant, a monotone split such as `BIG_INNER` (6 through `D0-1`) and `TAIL` (`>=D0`) gives tighter fibers but changes the generator: gaining at `D=D0-1` crosses the bucket boundary. To keep a pointwise proof with one BIG_INNER row, bound the crossing term jointly over eligible physical rates (e.g. zero or the sign-selected product of a safe crossing-rate bound and the bucket jump). A finite safe crossing rate is `D0-1` since `L+ρ <= L+(k+h)/2 <= D`; use zero for noncrossing states and `[0,1/D0]` for TAIL credit. Promotion leaves degree and bucket unchanged. This split still prohibits BIG interiors and does not increase receiver coverage. If finer fibers alone are the goal, the projected storage increase is 2,362,784 row bytes at 6/6; a numerical benefit remains unmeasured.
