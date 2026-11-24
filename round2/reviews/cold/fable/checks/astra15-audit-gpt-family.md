CLAIM 15: HOLDS

Checked inside the GPT family. This conclusion uses the earlier Claim 9 certificate for the raw signed residual cap, positivity, both receiver reserves, and mask63; I did not redo that work. I reviewed the two new proof steps, checked the exact composition and full finite dual prices for both signs, and independently computed the five-pass geometry and the reduced lower-mesh scores. The claimant's full lower production run and both sets of unsigned input table checks also passed on my copy.

## 1. Composition and source accounting

Run in this directory:

```text
$ taskset -c 13 nice -n 5 python3 -B compose.py
v=63134678 l*v-e=3555601807/100 endpoint=1356995264017948174596364819/2398901654969101844480000000 outward=0.565673570322
v=61537143 l*v-e=6922575559/200 endpoint=2725370105716046104606536597/4816977775848776007680000000 outward=0.565784239110
PASS: inward adjacent tick and substitution of upper prior endpoint rejected
```

Use the nonnegative auxiliary occupation weights at discount `u=1133/2000`. Integrating against `exp(-st)` and sending `s↓mu` gives `(u-mu)[T]≥0` even if `u=mu`, because the normalized `T` is bounded and nonnegative and its initial boundary term disappears. The resulting inequality in ticks per core is `sigma (lambda-mu)m≥mu*v-e`. The LP uses `l≤mu≤u` on nonnegative occupation counts; in the two-sided stationarity rows these give `flux≥l*mass` and `flux≤u*mass`. The benefit is strictly positive because `l*v>e` as shown above. Thus dividing by the positive `m≤S*mu*y` lowers the guaranteed benefit, and `v-e/mu≥v-e/l` uses `e≥0`. Lower sign (`sigma=-1`) adds the shift; upper (`sigma=+1`) subtracts it. The lower outward tick is floored; the upper is ceiled. The upper bound `u` cannot be substituted for `l` in that last denominator.

The physical ledgers distinguish a **receiving potential occurrence and its gain slot**, not just a donor vertex or an event. A single event can therefore finance several rows legitimately. Proof lines (paths relative to Astra):

- `round2/occupation/inside.md:19–23`: the old quartet source is a donor `z` outside the four-set with gain at its contact; “For F1 the identical formula uses the singleton potential f(x) at 72≤D_x<1399 … Each pair is directed towards x.” At :43, “These two tables book distinct directed incidences of F1.” The parent and child labels in the ledger are receivers.
- `round2/five-next/split-audit.md:22–27,33`: “Every leaf deletion is charged once”; each per-deletion candidate “receives that deletion's own source”; on later passes “each previous pass funds its immediate successor using this exact difference for old supports of sizes one through four.” The proof explicitly uses “only the original quartet donor allocation and the resulting auxiliary residuals of supports up to four.” Those later residuals belong to newly introduced auxiliary weights, not to the raw singleton/triple/quartet potentials.
- `round2/occupation-tail/continuation.md:7–11` lists raw singleton cap-minus-maximised-row, raw pair-row slack, edge donor, raw triple/quartet-row slack, and an uncovered **triple** donor. “Within each family the row slack is calculated with its outside inputs still maximized, and the actual donor pays a separate difference.” In particular the intrinsic quartet can share its support with the raw quartet donor gap without sharing its summand: `round2/motif-high/occupation/source.md:28–49` shows `R_actual=R*−sum J`, `I=floor(max(0,-R*))`, and bounds `R_actual+I+selected floored J≤max(0,R*)`. Its line 49 explicitly distinguishes the five-set quartet slot, the four-set uncovered triple slot, and the singleton. `round2/motif-cert/reserve/README.md:53` similarly distinguishes `cap_f−R_f,max` from the further physical-rate difference.
- `round2/joint-enclosure/arithmetic.py:40–61` checks the archived Q4 and increment against the old proof but emits only `quartet-replacement`, `parent`, `child`, `mask63` into the sum. The LP only prices unsigned potential, so this combination books no other signed credit. The raw physical gap is floored toward a lower credit (and the new source first truncates the signed gain change toward zero, which cannot increase `J` at fixed `p∈[0,1/2]`); upper raw residuals, per-state upper maxima and analytic allowances go in the other direction. The K7/144 score reductions are downward bounds on excess allowance, so subtraction leaves an upper comparison.

## 2. Short replay on the copy

Source copy from `git -C ASTRA_REPO archive d87dc4c round2 | tar -x -C copy`. The unmodified command on CPU 13 deliberately fails with `ValueError: run with taskset -c 3`. I changed only the copied wrapper's affinity guard and `MOMENTS_CPU` from 3 to 13. The `round2` archive omits `verify/tree-drift/motif_identity.py`; copied that listed static input and verified its SHA-256 matched `git -C ASTRA_REPO show d87dc4c:verify/tree-drift/motif_identity.py | sha256sum` (`28bff16a65a623f4d750818c645f93b41a3b4d296cf34e20f9ddeee1fb260806`). Used symlinks pointing to the existing snapshots and mbar data for read-only access. With those inputs:

```text
$ taskset -c 13 timeout 60s nice -n 5 python3 -B round2/joint-enclosure/replay.py quick
quick: snapshot and input digests; frozen transcripts; exact joint arithmetic
rejected archived Q4 overlap
rejected duplicate mask63
rejected unreviewed mask64
rejected archived Q4 hidden in replacement
rejected narrowed growth prior
rejected flipped residual
rejected wrong upper shift sign
rejected inward lower tick
rejected inward upper tick
lower: v=63134678, e=115075, y=126536000881/62500000000, lambda=5656453/10000000, prior=113/200..1133/2000, endpoint=1356995264017948174596364819/2398901654969101844480000000
upper: v=61537143, e=155608, y=1016333625559/500000000000, lambda=11316233/20000000, prior=113/200..1133/2000, endpoint=2725370105716046104606536597/4816977775848776007680000000
quick PASS: outward physical interval [0.565673570322, 0.565784239110]
$ taskset -c 13 timeout 60s nice -n 5 python3 -B round2/joint-enclosure/replay.py preflight
mandatory imports PASS: low_oracle, tree_check, physical_check, mutations
preflight: isolated 74 static inputs; no raw, mbar or executable copied
```

`quick` checks frozen manifest and snapshots against digests, unsigned sheet identity/flags/source digests and scale/nominal header, checks the split transcript and pinned mask63/Claim 9 ledger, and reruns exact composition with plants. It does not prove the signed rows, recompute the five-pass production weights, or verify the LP dual. `preflight` checks isolated staging/imports and original prior, nonraw split hashes and sheet identities; it does not independently validate the omitted raw/mbar or production proofs. No 53-minute cold replay run.

## 3. Unsigned denominator

The joint claim uses `moments-v5-unsigned-fine-{lower,upper}-dual.json.gz` (K7/144, **second** numerical row of `moments/RESULTS.md`). The last, `unsigned-birth`, row is a different, slightly tighter model with affine-in-birth score savings and is not used by Claim 15. The K7/144 sheets do still use each tail bin's actual promotion intensity for moment constraints (`degree_tests=true`, `tail_score` absent/false). They prove `D(x)≤y*nu(x)` on the physical relaxation, with `D= mbar − unsigned_score_defect` in unscaled units per core. `D` majorizes the normalized original **unsigned** snapshot potential `m/S` (no multiplication by `sigma`, no absolute value of each term). `nu=[B]=mu` by the Abel integration, so this is `m≤S*mu*y`, with positivity from the older comparison. The other checked objective is `D≤q*[1]`. The model prices no raw residual cap or donor source.

`occupation.py` and `marginal_lift.py` model finite vertex and oriented edge gains with rate `D−c+sum 1/D_neighbor`, promotion at `(D−c)/D`, with exactly one new vertex `(2,1)` and one new edge from the **post-promotion** state. At/above the cutoffs they separate the exact moment generator from inequalities that enclose it: joint degrees >7 are represented with two auxiliary reciprocal colours 0 and 1/8, with correct first moment and independent colours for distinct vertices. The same vertex's colour is reused wherever it recurs. Overlap marginals constrain oriented paths, bare connectors, cherries and the four-site populations. Finite finer vertex states up to 144 use their actual degree and reciprocal flows. Further inequalities bound the count and generator of `(1−T/D)_+ 1{k>=h}` in the tail cells by one-sided monotonicity. The `l/u` on each stationary nonnegative count and the explicit `l≤nu≤u` use only the published broad bracket. The new-edge birth and child-degree changes occur at promotion, rather than counting a leaf gain as a birth.

The box caps used in rational certification are physical too: a particular oriented edge or path count is at most total edges (a bottom fixes its ancestors); cherries with finite apex degree ≤7 are bounded by `K` times total edges using `choose(k,2)≤K*k`. At a finite top or end of a connector test the remaining input is at most `(c−1)/2` per occurrence. Finer-vertex reciprocal sums and all tail parent/child rate sums are bounded by half the same incoming/outgoing edge count; total bins' mass and child-slot rates are bounded by total vertex and edge counts, respectively. Bounded observables/fluxes yield subsequential limits of this finite set of coordinates without requiring pointwise convergence of a full joint tail law.

I used the claimant's primal model builder but neither its optimizer nor either dual verifier. My separate verifier checks sparse row indices, uniqueness, allowed signs and every column, explicitly supplies the nonnegative box price `beta_j=max((Aᵀp+Gᵀw)_j−c_j,0)` for `x_j≤U_j`, verifies `Aᵀp+Gᵀw−beta≤c`, then recomputes `b·p+h·w−U·beta` exactly and matches the stored rational margin. It verifies **dual** feasibility, not a numerical optimum or primal optimizer's vector. Primal feasibility for physical Abel limit points is supplied by the occupation argument, not by price arithmetic.

```text
$ taskset -c 13 nice -n 5 timeout 2400 bash -c 'ulimit -v 8000000; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1; /usr/bin/time -f "wall %e s; peak %M KiB" python3 -B exact_dual.py'
lower columns=73586 equalities=3529 inequalities=50209
M-per-core bound=114689698691/100000000000 margin≈7.578152671508e-11 exact_margin_matches=YES nonzero_box_prices=14403
M-per-birth bound=126536000881/62500000000 margin≈4.149946434999e-11 exact_margin_matches=YES nonzero_box_prices=12694
PASS 2% smaller ceiling rejected by this dual: margin≈-0.501082563
PASS negative inequality price and both out-of-range row indices rejected
upper columns=73586 equalities=3529 inequalities=50209
M-per-core bound=575732818361/500000000000 margin≈8.996226518047e-12 exact_margin_matches=YES nonzero_box_prices=14182
M-per-birth bound=1016333625559/500000000000 margin≈3.46069977408e-10 exact_margin_matches=YES nonzero_box_prices=9223
PASS 2% smaller ceiling rejected by this dual: margin≈-0.503085144
PASS negative inequality price and both out-of-range row indices rejected
wall 25.38 s; peak 479816 KiB
```

**Denominator multiplicities and tails.** In `motif-cert/comparison.cpp:26–49,64–72`, one `f` per vertex; one edge, vertical and chain per bottom; all `k_top−1` possible fork side children per bottom; and `C,Y,T` at the apex with `choose(k,2),choose(k,2),choose(k,3)`. The five unsigned deduction indices do not duplicate the same allowance: incoming edge price (E only); from small-middle triples the vertical and possibly chain/fork; if their parent is high, the extra conditional triple price uses only the chain/fork omitted by the first; a four-site bare-chain price spends the remainder after maximizing over the top; small-apex cherries allocate a pair's C/Y and `1/3` of each star spread over all three constituent pairs; high-middle/apex outgoing-edge prices instead charge the full vertical/chain/fork to the bottom edge and divide C/Y by 2 and T by 3 across their constituent child edges. In the high-apex case this leaves `(k−1)/2,(k−1)/2,(k−1)(k−2)/6` multiplicities. Equal child **states** are not equal physical child **vertices**. The signed `frow/grow` columns have coefficient zero in the score-only objectives/constraints.

The unbounded tail is not silently set to zero: cells cover finite D=145…1399; for D≥1400, the singleton uses the exact affine/quadratic leaf-fraction envelope with the nonpositive reciprocal correction dropped for lower, and each of the four incoming bottom allowances is affine in `1/D` up to upward integer rounding (four extra ticks). At D≥1400 all internal apex supports have already stopped. The caps include the interval `1/D∈[0,1/1400]`, every compatible clamped child count, and leaf-fraction panels. The finite dual does **not** prove these analytic envelopes, the physical applicability of the score deductions, the mean-growth/Abel theorem, or the founder localization. Beyond the joint cutoff the model only keeps selected moments, not the full tail joint law. Founder and adjacent omitted terms need the separate proof: their counts have polynomial expectation in continuous time at fixed radius (`occupation/physical-test.md:5–17`); the actual rate-weighted gain differences stay bounded via the reciprocal interpolation and the analytic singleton, or the complete omitted generator is integrated by parts. Founder-origin births are included in total physical `B`, and neglecting the founder's own rate changes its normalized value by zero. Counting all ordinary births once gives `B/C=1−sum_edges(1/D_a+1/D_b)/C` in the limit; there is no factor two or division by the overall attachment count.

I ran the claimant's independent Python decoder/table checker (the system `python3` first failed on missing NumPy; used the existing project Python only as interpreter):

```text
$ taskset -c 13 nice -n 5 timeout 220 bash -c 'ulimit -v 8000000; export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1; for side in lower upper; do ASTRA_REPO/round2/.work/venv/bin/python -B round2/joint-enclosure/check_tables_k7.py "$PWD" "$side" || exit; done'
lower: independently checked 979299 mbar entries and analytic tail; profile={'v': 10439, 'b': 173, 'e': 4189}, overlap={'t': 12934, 'q': 12348, 'c': 6090}, conditioned=3286, tail-lines=1182
upper: independently checked 979299 mbar entries and analytic tail; profile={'v': 10439, 'b': 173, 'e': 4189}, overlap={'t': 12934, 'q': 12348, 'c': 6090}, conditioned=3286, tail-lines=1182
```

I also wrote `physical_score.py`, which exhaustively enumerates connected 2–4 member subsets of each fixture, independently classifies the shapes, sums actual unsigned weights with the frozen `Fraction` reader, and compares against the sum of all deduced prices and the pointwise vertex cap. This check uses actual bins/caps, including a separate reconstruction of every finite tail-bin range, and runs with degrees crossing the joint, lift, frozen-table and analytic seams (through 2,000,003). It excludes founder-localized terms on both the potential and deduction sides and globally includes the true founder core degree for the birth identity. Across 800 fixtures per side it exercised every family and every group of deductions. The largest actual/majorant ratios cross 0.99, so this is a sensitive multiplicity check, independent of the lead's simulations. Lowering the **pointwise** majorant 2% or double-spending deductions causes definite fixture failures; this is distinct from the 2%-smaller **stationary** bound rejected by the price sheet. The raw snapshot header's built-in proposal lambda differs from the nominal trial exported in the profile; my initial assertion of equality failed. The trial is passed explicitly to signed-row builders (`v5_profile.cpp:27–40`); the unsigned denominator reads the stored weights and does not use either trial's signed caps.

```text
$ taskset -c 13 nice -n 5 timeout 120 bash -c 'ulimit -v 6500000; /usr/bin/time -f "wall %e s; peak %M KiB" python3 -B physical_score.py'
lower: 800 fixtures, active E/V/C/P/F/Y/T={'E': 9313, 'V': 6748, 'C': 3731, 'P': 2221, 'Y': 1385, 'F': 2646, 'T': 758}; price/plant counts={'E': 800, 'vertical': 747, 'cherry': 423, 'bare': 136, 'double-all breaks ceiling': 800, '2% smaller breaks ceiling': 584, 'conditional': 707, 'high-parent': 617}; max phi/majorant=0.992185364; min slack ticks≈9.09675e+10
upper: 800 fixtures, active E/V/C/P/F/Y/T={'E': 9313, 'V': 6748, 'C': 3731, 'P': 2221, 'Y': 1385, 'F': 2646, 'T': 758}; price/plant counts={'E': 800, 'vertical': 747, 'cherry': 423, 'bare': 136, 'double-all breaks ceiling': 799, '2% smaller breaks ceiling': 440, 'conditional': 707, 'high-parent': 617}; max phi/majorant=0.991124993; min slack ticks≈1.34439e+11
wall 20.92 s; peak 256448 KiB
```

`generator_rows.py` exercises the full K7/144 model on locally tree-like cyclic fixtures (to remove the founder seam without creating cycles on a support of size ≤4). It uses the claimant's low-joint colour/one-step oracle, separately collects actual high tail/group variables, tests the actual rates against all degree-test generator enclosures, verifies the nonstationary constraints and `[B]` precursor and plants the three listed mistakes. Stationarity inequalities are deliberately not evaluated as if the finite instantaneous fixtures were stationary. The generator model is common to the two signs; the numerical fixture test ran with the lower profile.

```text
$ taskset -c 13 nice -n 5 timeout 210 bash -c 'ulimit -v 6500000; /usr/bin/time -f "wall %e s; peak %M KiB" python3 -B generator_rows.py'
12 cyclic K7/144 fixtures: 73586 box bounds and 3529 equalities each; 296988 pointwise row checks total
8 cyclic K7/144 oracle fixtures: 12289 exact generator flows each and all pointwise entry bounds
PASS plants (fixture hits): tail B=mass 12, missing newborn intensity 12, half sibling incidence 6
wall 112.04 s; peak 343920 KiB
```

## 4. Five-pass quartet replacement

Each chronological five-set and deletion of an undirected leaf corresponds to a unique raw `(quartet, contact, outside donor)` incidence, including the oldest member when it is an undirected leaf. Multiple deletions correspond to distinct raw summands. For pass 0's two whole-source candidates each inequality uses the same entire source, while the summed deletion candidates split that source. Taking `max(full_1,full_2,sum_d piece_d)` is valid: at a fixed state, take a maximizing candidate; its current value equals the chosen weight and every successor value of the maximum is no smaller. This argument uses nonnegative transition rates and post-promotion old descriptors. The newborn and expansion maps are bijective by deletion of the youngest support member, since only the youngest can be the new child at the birth event.

In particular the new credit is not restricted to the archived small D≤6 window of `inside.md:23`. In `higher-rows.md:45–47,64–90`, every quartet outside slot remains uncentered. Fix the quartet's physical degrees and actual post-gain/post-birth weights. Its true transferred row depends on an exterior reciprocal only through that slot's signed gain coefficient, not the quartet's internal interpolation or incoming triple transfers. The raw numerical upper proof still applies if this unobserved reciprocal is set to either end of `[0,1/2]`, independently at each slot; its tail own-rate bounds continue to hold. This is the all-degree maximized-row argument in `motif-high/occupation/source.md:34–36`. Thus one can claim the raw slot gap for retained rows even with other endpoints interpolating, or if the contact gain exits quartet support, without enlarging the raw cap. One may set the gap of any unretained incidence to zero.

On each exact coordinate the robust row minimizes `rate*d` over the interval with independently bounded exterior reciprocal slots. Its actual-minus-worst difference on one external slot is `p*d−min(0,d/2)=J(−d,p)`, where `d` is the actual old-weight jump. Later passes use this gradient after the previous candidate combination (`reserve.cpp:65–74,92–97`), and spend only that preceding pass's leftover. At a previous omitted zero row where a gain enters support, its positive gain term already supplies this slack; in a stable pooled degree label the old score gain is zero. Missing labels/child entries may contribute zero; entries for ambiguous child states take a lower bound over compatible states. Given a pooled donor, `p∈[0,1/d_min]`; the source floors `J` at the minimizing endpoint. Pooled raw endpoint gradients keep the recipient's true before/after comparisons under the same interpolation of other coordinates; outside support a quartet is zero and the support exit jump is included. The integer row checks use 128-bit signed arithmetic for `2000*2*prod(D)` times the residual; the floating score proposal is accepted only if the exact residual is nonnegative. Summing the five inequalities cancels the first four leftover sources; the fifth leftover need not be spent. The uniform expected core-degree moment estimate controls the integrated auxiliary scores (including multiple tail vertices). See `flow.hpp:66–95,98–143`, `quartet.hpp:55–112,115–142`, and `split-audit.md:21–33`.

Independent mesh census: `ranks.py` enumerates chronological parent lists, uses closed-form finite-axis counts plus pooled axes, independently classifies leaf deletions and adapted meshes, and multiplies axis widths. Rank counts include allocated rows even if their value is zero/inactive. Geometry is identical for both signs; this does not calculate scores.

```text
$ taskset -c 13 nice -n 5 python3 -B ranks.py
pass 0: [(5, 2021486872), (4, 56279691), (3, 70491064), (2, 32031570), (1, 1061159)]
pass 1: [(5, 7312312), (4, 6449232), (3, 2975049), (2, 976140), (1, 52278)]
pass 2: [(5, 1177302), (4, 381038), (3, 2975049), (2, 976140), (1, 52278)]
pass 3: [(5, 1177302), (4, 381038), (3, 2975049), (2, 976140), (1, 52278)]
pass 4: [(5, 1177302), (4, 381038), (3, 2975049), (2, 976140), (1, 52278)]
PASS all 50 entries; total=4431601576
PASS omitting oldest donor rejected: rank-5 pass-0 1796257670 vs 2021486872
```

Compiled `reserve`, `gradient_probe`, `census` from the archived copy with `c++ -std=c++20 -O3 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=undefined` on CPU 13 with 5 GiB VA limit (exit 0). Claimant small replays, CPU 13/nice 5 with `TMPDIR` inside this directory, from `copy`:

```text
$ python3 -B round2/five-next/check_slopes.py "$root/run"
upper 2560 exact quartet slopes, 6720 compatible pooled comparisons, 4908 nonzero pooled bounds; 1920 multiple endpoint cases, 998 nonzero; 6475 infinity-only donor overcharges
lower 2560 exact quartet slopes, 6720 compatible pooled comparisons, 5099 nonzero pooled bounds; 1920 multiple endpoint cases, 1048 nonzero; 6373 infinity-only donor overcharges
upper four physical family totals and four founder-donor totals agree with Fraction
lower four physical family totals and four founder-donor totals agree with Fraction
$ python3 -B round2/five-next/check.py "$root/run" small-split
lower small-split 3572 trees 31935 query configurations {'F': 2338, 'P': 451, 'T': 96, 'Y': 505, 'elementary successors': 34402, 'inflated reserve rejected': 3571, 'outside-current-table': 1403, 'wrong-sign overcharge': 1003}
upper small-split 3572 trees 31935 query configurations {'F': 2401, 'P': 631, 'T': 102, 'Y': 537, 'elementary successors': 34402, 'inflated reserve rejected': 3571, 'outside-current-table': 1403, 'wrong-sign overcharge': 630}
```

**Independent score recomputation on the reduced mesh.** `bellman.py` reimplements the chronological geometry, pooled labels and ambiguous compatible destinations, interpolated signed sources, previous-pass gradients, and all Bellman rows on the `small-split` mesh. It reads physical weights via Astra's frozen snapshot reader, but neither imports nor calls Astra's geometry/solver. Rather than a floating proposal it finds the precise piecewise-affine residual root by increasing trial value through the successors in order, switches endpoints at each crossed breakpoint, and floors the resulting exact ratio. It then verifies its own integer inequalities. The three lower passes agree exactly with the claimant's small-split run; the smallest mesh reproduces real scores, not just dimensions. A separate source-mutation run removes oldest donors and changes the first newborn score substantially.

```text
$ taskset -c 13 nice -n 5 timeout 200 bash -c 'ulimit -v 4000000; /usr/bin/time -f "wall %e s; peak %M KiB" python3 -B bellman.py lower'
lower independent small-split pass 0: 6046616 ticks (matches replay)
lower independent small-split pass 1: 812794 ticks (matches replay)
lower independent small-split pass 2: 75136 ticks (matches replay)
lower total=6934546
wall 126.99 s; peak 424340 KiB
$ taskset -c 13 nice -n 5 timeout 140 bash -c 'ulimit -v 4000000; /usr/bin/time -f "wall %e s; peak %M KiB" python3 -B bellman.py lower omit-oldest'
lower PASS oldest-donor plant changes small-split pass 0: 3560916 vs 6046616
wall 57.17 s; peak 256452 KiB
```

**Full production replay (claimant's code, lower only).** Run from `run/` with its scratch projection subdirectory present, directly reading the frozen lower witness. See `run/lower.log` for all 37 lines, including the 25 rank counts independently checked above. I paused this job while running other CPU-13 checks; the reported wall time includes pauses.

```text
$ taskset -c 13 nice -n 5 timeout 2400 bash -c 'ulimit -v 5242880; /usr/bin/time -f "wall %e s, max RSS %M KiB" ./reserve ASTRA_REPO/round2/outputs/mc-lower.bin split 5'
lower quartet pass 0 newborn 29374764/1099511627776
lower quartet pass 1 newborn 5855805/1099511627776
lower quartet pass 2 newborn 890500/1099511627776
lower quartet pass 3 newborn 146067/1099511627776
lower quartet pass 4 newborn 22982/1099511627776
lower quartet replacement 36290118/1099511627776
wall 1594.08 s, max RSS 2866076 KiB
```

```text
$ taskset -c 13 nice -n 5 python3 -B -c 'from pathlib import Path; a=Path("run/lower.log").read_text().splitlines(); b=Path("copy/round2/outputs/five-next-split.log").read_text().splitlines(); assert a[:36]==b[:36]; print("PASS: lower production lines 1–36 match frozen transcript byte for byte")'
PASS: lower production lines 1–36 match frozen transcript byte for byte
$ sha256sum run/round2/five-next/build/primary-lower.bin
c1f1a1f7de60024992d38aa6295c6aaa26328ba4d544c537917548f4a0681fbf  run/round2/five-next/build/primary-lower.bin
```

The projection hash also matches `five-next-split-cold-receipt.txt:4`. The upper **production scores** are checked via the pinned transcript and prior production evidence only; I ran the upper small-split and slope controls and independently recomputed the production rank counts for both signs. No independent production-mesh score computation on either sign.

## 5. Coverage limits

No independent full production solver or all-row raw scanner. No re-proposal/optimizer or reproduction of the LP's matrix by a second implementation; no independently proved LP optimum. The full numerical unsigned table check is the claimant's Python checker, separate from its C++ extractors; I did not regenerate the C++ profiles in this audit. No independent proof of the published prior, renewal/mean comparison, or the asymptotic exponent result. No all-graph fixture sweep; physical tests are selected configurations and complement the proof. The physical-score fixtures rely on the claimant's frozen weight decoder and do not cover the founder terms pointwise; the row fixtures do not model the founder or empirical convergence to stationarity. No new Claim 9 recomputation, no 53-minute joint replay, no stronger `unsigned-birth` denominator or renewal-kernel result, no eight-decimal claim.
