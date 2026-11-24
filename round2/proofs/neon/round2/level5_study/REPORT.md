# Priority C: B3-FullH width and cost study

**Authors:** Gemini 4 Neon (lead) and GPT-6.2 Sol (prototype and checks). **Archive pin:** `c2498cc482bfef861644092d33f727cbbc93b280`. **Review and dependencies:** the study's GPT exact transfer fixtures and Gemini lead review cover diagnostics using the [Level-4 potential](../proof-b3f-motif4.md), not a separate drift certificate. See [auxiliary evidence](../../README.md#auxiliary-evidence).

## Scope and reproduction

All commands below run from `round2/level5_study/`. Each solve/verify uses one OpenMP thread. Inputs: Neon commit `c2d12bf1b22369bdc8eafe4120ba33831479a2b8`; builder SHA-256 `0942bd2ca029db576cce5b5d817ab729dd738a888436c4afe1f5a4e6cc22c67f`; format SHA-256 `3064dff6311e460291d277a020dd73ebbbbdb32f7f8c0258159a5543d470f180` (GCC 13.3.0).

```sh
export TMPDIR="$PWD" OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1
python3 make_diag.py
gcc -O3 -march=native -ffp-contract=off -fopenmp diag_slack.c -lm -o diag_slack
bash run_diag.sh
python3 summarize.py
python3 l5_count.py
gcc -O2 -Wall -Wextra proto_next.c -lm -o proto_next
./proto_next
```

`run_diag.sh` runs both endpoints of each case sequentially. Logs record the binary64 `lam` at each sweep and final `lam` to 17 significant digits, wall/CPU time, and peak RSS per endpoint. `summarize.py` prints only completed pairs and forms the width by binary64 subtraction. The grid uses `D0=12, KT=11, LF=6000, KF=100`, P4/Yt/Yb/S3 `(10/8,9/7,10/8,8/6)`, and 18 floating sweeps. The D0=14 repeats were removed at the user's direction; their layout counts remain for cost comparison. No final target rounding or fixed-lambda safety sweeps are used in the grid.

### What each ablation measures

- `base`: original float equations and empirical rate sampler. To avoid double allocation, `g`, hV/hC and L4 rebuild in place after their last previous-sweep consumer has finished; `f` keeps its two arrays. Updates and normalizations occur in the original order.
- `center`: only the L4 external core-neighbor rates are fixed to role-wise empirical centers. This applies to inner vertices and the rate arguments of known-`k` T2/T3. T1 has no such rate dependence; the `k>KT` overflow hull uses `alpha`, with unknowable `k` in that row, and is left unchanged. This measures the effect of removing those interval extrema at those chosen centers. It is not a theorem about the largest possible improvement from future L5+ potentials.
- `d5-6`, `d5-7`, `d5-8`: only those same L4 rate intervals are changed to `[0,c_ext/D5]`, with lower endpoint zero and no added centering or solving of L5. This is an optimistic attribution; it is not a valid drift bound on general trees.
- `h-center`: only hV/hC T1/T2/T3 rows use a single representative instead of each relevant tail extremum. T2 uses the midpoint of its total-rate interval; known-`k` T3 evaluates the rational Bernstein row at `t=tmax/2` and the midpoints of both rate intervals; overflow T1/T2 use `alpha=1/2` and overflow T3 uses `alpha=t=1/2` and the midpoint of its affine incoming-credit term. Inner h rows and the g and L4 equations retain their original envelopes. This isolates sensitivity to the h tail envelope for the stated representative; choosing a representative does not produce a rigorous inequality or a mathematical upper bound on a future improvement.

The ablations re-solve all tables, and their effects are not additive. Existing `h → L4` transfers and original `g` Bernstein bounds remain active in every case.

## Completed floating runs

| D0 | Case | λ lo | λ hi | Width | Reduction | Wall lo/hi (s) | Peak RSS lo/hi (KiB) | Max last-two-sweep shift |
|---:|:---|---:|---:|---:|---:|:---|:---|---:|
| 12 | base | 0.56565788754263469 | 0.56579149436263099 | 0.00013360681999630142 | 0.000% | 215.76/168.31 | 661528/661656 | 2.6e-10 |
| 12 | center | 0.56568394486865603 | 0.56576480694192277 | 8.0862073266740886e-05 | 39.478% | 173.29/147.32 | 661780/661652 | 9.4e-11 |
| 12 | d5-6 | 0.56567495571581106 | 0.56577823581536657 | 0.00010328009955551387 | 22.698% | 174.15/156.25 | 661652/661532 | 1.5e-10 |
| 12 | d5-7 | 0.56567645783386855 | 0.56577698989866121 | 0.00010053206479265775 | 24.755% | 182.01/185.27 | 661652/661652 | 1.4e-10 |
| 12 | d5-8 | 0.56567762241469222 | 0.56577601440211078 | 9.8391987418566806e-05 | 26.357% | 224.72/199.80 | 661652/661780 | 1.4e-10 |
| 12 | h-center | 0.56566754344890324 | 0.5657804535430605 | 0.00011291009415725917 | 15.491% | 198.41/203.95 | 661652/661780 | 2.4e-10 |
| 12 | l5-6 | 0.56566050397627932 | 0.56578850648188217 | 0.00012800250560285331 | 4.195% | 198.95/202.56 | 674024/674276 | 1.5e-09 |
| 12 | l5-6-no-inner-slack | 0.56566218554896741 | 0.56578692013676202 | 0.000124734587794606 | 6.641% | 175.48/195.55 | 674024/674024 | 1.5e-09 |
| 12 | l5-6-no-mixed-slack | 0.5656846387582406 | 0.56576972204110343 | 8.5083282862830067e-05 | 36.318% | 130.38/133.92 | 674280/674024 | 5.8e-10 |

`sample_build`: `TIME wall=284.64 user=278.81 sys=1.07 peak_KiB=667092 exit=0`.

`sample_verify`: `TIME wall=110.42 user=109.46 sys=0.45 peak_KiB=647424 exit=0` — `=== VERDICT: PASS (total_fails = 0, exact P/Q = 565655/1000000) ===`.

GMP benchmark: 82,492,914 printed row slots; 109.91 CPU seconds; 1332.4 ns per printed slot (includes startup and comparison).

Approximate CPU hours per endpoint by original table size:
- D0=20 production: 0.52 h
- D0=22 production caps: 0.80 h
- D0=24 base: 0.65 h
- D0=24 large: 1.24 h
- D0=25 large: 1.50 h

Real clipped L5 (D5=6): width change 4.195% narrower versus D0=12 baseline (below the 30% decision threshold).
- Suppressing all-inner residual: width falls by 3.2679178082473115e-06 (2.553% of the real L5 width; non-rigorous sensitivity).
- Suppressing mixed residual: width falls by 4.2919222740023244e-05 (33.530% of the real L5 width; non-rigorous sensitivity).

## Exact table sizes and memory

Potential entries include allocated `g` cells that the verifier does not visit. File size includes the header and int64 rate/reference arrays. The GMP program reads one set of int64 arrays into RAM; its RSS also includes layout, the runtime, scratch and tree data. The unmodified builder allocates two potential copies plus layout/caches/write buffer. The diagnostic avoids a second copy of every table except `f`.

| D0/KT | L4 cutoffs | Potential entries | File bytes per endpoint | Verifier file-equivalent GiB | Original builder potential GiB |
|:--|:--|--:|--:|--:|--:|
| 12/11 | 10/8, 9/7, 10/8, 8/6 | 82,520,610 | 660,209,960 | 0.615 | 1.230 |
| 14/13 | same | 153,652,753 | 1,229,289,152 | 1.145 | 2.290 |
| 16/15 | same | 270,633,133 | 2,165,160,448 | 2.016 | 4.033 |
| 20/19 | 12/9, 11/8, 12/9, 10/7 | 1,393,618,233 | 11,149,119,456 | 10.384 | 20.767 |
| 22/21 | 12/9, 11/8, 12/9, 10/7 (LF/KF=16000/200) | 2,155,055,415 | 17,240,668,400 | 16.057 | 32.113 |
| 24/23 | 12/9, 11/8, 12/9, 10/7 (LF/KF=16000/200) | 3,230,810,807 | 25,846,772,304 | 24.071 | 48.142 |
| 24/23 | baseline cutoffs | 1,748,953,263 | 13,991,911,952 | 13.031 | 26.062 |
| 24/23 | 12/10, 11/9, 12/10, 10/8 | 3,346,876,551 | 26,775,298,256 | 24.937 | 49.872 |
| 25/24 | same larger cutoffs | 4,060,908,685 | 32,487,589,432 | 30.256 | 60.512 |
| 26/25 | same larger cutoffs | 4,899,342,893 | 39,195,099,808 | 36.503 | 73.006 |
| 30/29 | baseline cutoffs | 5,318,553,885 | 42,548,962,592 | 39.627 | 79.252 |

These counts use `./diag_slack --threads 1 --D0 D --KT K --count`, with `--DP4 12 --KP4 9 --DYT 11 --KYT 8 --DYB 12 --KYB 9 --DS3 10 --KS3 7` for D0=20, and `--DP4 12 --KP4 10 --DYT 11 --KYT 9 --DYB 12 --KYB 10 --DS3 10 --KS3 8` for the larger D0=24/25/26 rows. Even if all L4 tables were free, `f+g+hV+hC` at D0=30 would still store 401,107,989 entries (2.988 GiB as int64s). At roughly 25 GiB RAM available at the start of the study, D0=24 with larger L4 would need background allocations to finish before it could run without pressure. A strict **32 GB decimal** budget rules out the D0=25 larger row before overhead; **32 GiB** allows one such table set with under 2 GiB to spare.

## Exact candidate transfer and why the oracle is optimistic

For a physical L4 motif `M`, a slot `u` in it and an external core neighbor `z`, let `Δ_u q(M)=q(M after gain at u)-q(M)` (including zero if a cutoff is exited). At eligible rows (all four degrees `<D5`) use

```
x_z = max(1/D_z - 1/D5, 0)
(1/D_z) Δ_u q(M) = min(1/D_z,1/D5) Δ_u q(M) + x_z Δ_u q(M).
```

Only `D_z<D5` has a nonzero `x_z`. Charge `x_z Δ_u q(M)` to the induced five-vertex tree `N=M∪{z}`. It is connected and the unique extra edge is `uz`, so every term corresponds to deleting a leaf `z` of a five-vertex motif; sum over all deletions that give active L4 contributions. This pairing is exact in signed arithmetic and survives either sign of `Δ_u q`. `proto_next.c` tests the clipped-rate equality over 5,650,239 cases (`max_abs_error=6.9388939039072284e-18`).

An unqualified `[0,c_ext/D5]` in mixed L4 rows is unsound with all-five-inner support. Example: choose P4 with one endpoint of degree `≥D5`, add a degree-2 core child outside it at a different slot, and ensure the latter slot can gain a leaf in the active P4 domain. That extra neighbor contributes `1/2` and no L5 receiving table exists. Keep `[0,c_ext/2]` in mixed rows, or allow some old L4 slots to be outside the pure-inner cutoff in L5. Even a fully tracked L5 reroutes drift; it does not delete it.

| D5 | Pure-inner L5 entries, all nine rooted shapes | int64 bytes | two copies (bytes) |
|--:|--:|--:|--:|
| 6 | 1,564,220 | 12,513,760 | 25,027,520 |
| 7 | 11,404,845 | 91,238,760 | 182,477,520 |
| 8 | 59,172,498 | 473,379,984 | 946,759,968 |

The exact recurrence counts each rooted shape separately, split by apex `h=0/1`. At node with `a` motif children, the number of intrinsic `(D,k,h)` states is `sum_{2≤D<D5} max(0,D-h-a+1)`; equal child shapes use multisets of their recursively enumerated states. Counts include all physical pure-inner rows, before any pruning tied to L4 cutoffs. Costs exclude L5 indexing, comparison, transferred-source and budding arithmetic, and above all any mixed rows needed to approach the optimistic ablations.

### Coupled L5 prototype at D0=12

`proto_l5.c` implements the transfer rather than merely shrinking the interval. For each L5 row, the source is the sum over leaf deletions of `x_z Δ_u q(M)`; an exact degree increment may take the affected L4 slot to `D5` while keeping `q4` active. On existing L5 motifs, gain includes the known `L+sum_inside(1/D)` and an external interval `[0,c_ext/2]`; the state after exiting `D<D5` contributes zero. Budding at `u` acts at rate `L_u/D_u` and swaps `(L,k) → (L-1,k+1)`. Every newborn L5 motif is charged to its unique preimage obtained by deleting the newborn leaf; the L4 bud values include these terms before multiplying by the same rate. All-rooted-shape indexing identifies equal siblings without dropping multiplicities from generator sums. The same exact inequalities can be checked rationally, but this float prototype does not certify numerical endpoints or comparison constants.

For the eventual exact comparison, charge every L5 motif to its apex. If `D_apex >= D5` the charge is zero; otherwise a uniform upper bound on the number of supported 5-motifs rooted there is `4! (D5-1)^4`, by growing a four-edge connected motif one vertex at a time (step `j` offers at most `j(D5-1)` choices). Multiplying this by the extrema of `q5` including zero gives a valid, possibly loose, correction to the existing per-vertex comparison. To avoid a loose correction one can bound charge shape by shape per apex state. `Phi(T3)` is unchanged by L5 terms. A production verifier needs exact rational checks of the modified L4 and new L5 rows plus positive comparison constants where required.

```sh
python3 make_diag.py
gcc -O3 -march=native -ffp-contract=off -fopenmp diag_slack.c -lm -o solve_l5
OMP_THREAD_LIMIT=1 /usr/bin/time -f 'TIME wall=%e user=%U sys=%S peak_KiB=%M exit=%x' ./solve_l5 --threads 1 --l5 6 --l5-check
bash run_diag.sh l5
```

The D5=6 self-check visited all **1,564,220** rows in both modes, gave all-zero output for zero source and **1,083,040** nonzero rows per mode for a varying P4 source (4.83 s; peak 220,288 KiB). At D5=6 only **85,120 of 78,586,725** L4 storage rows are fully eligible (~0.108%); at D5=7 the count is **405,640** (~0.516%). These are row coverage fractions, not drift-slack fractions. The subsequent `l5-6-no-inner-slack` and `l5-6-no-mixed-slack` runs each remove one class of remaining external-rate interval, holding the actual L5 transfers in place. `no-mixed-slack` also removes the known-`k` L4 tail rate extrema, leaving `k>KT` overflow hulls unchanged. These two auxiliary runs are non-rigorous sensitivity tests.

`python3 check_transfer.py` independently checked source and birth ownership using exact `Fraction` arithmetic on physical core trees: **675 tree/cutoff pairs, 38,262 external credits, 29,361 exact transfers, 5,766 cutoff-crossing gains, 26,160 newborn motifs**, all passed.

## Recommendation and verification-cost plan

**Prioritize scaling the existing B3-FullH certificate to D0=22 over implementing a production pure-inner D5=6 verifier.** The real D5=6 gain is only **4.195%**, and the width sensitivity to mixed L4 rows is **13.13×** the sensitivity to all-inner rows. Even the optimistic D5=8 test is below 30%; this is evidence about these configurations, not a bound on every L5 architecture. Separate Neon production logs [`cert_b3f_20_lo.log`](../outputs/cert_b3f_20_lo.log) and [`cert_b3f_20_hi.log`](../outputs/cert_b3f_20_hi.log) report **GMP PASS on both** ends of `[0.565680, 0.565768]` (width **8.8e-5**) with `LF/KF=16000/200`, L4=`(12/9,11/8,12/9,10/7)` and **11,190,881,056 bytes per endpoint**; these are external runs, not rerun here. D0=22 with the same cutoffs costs **17,240,668,400 bytes per endpoint** before overhead, a cheaper next measurement than jumping straight to D0=24. The D0=18 exact reference has width **1.02e-4**. **5e-5 is an experiment target; 3e-5 and 2e-5 lack measured support.** At D0=12, even the fully centered L4-rate run only reaches **8.0862e-5**.

The D0=22/24 production-cap counts used the D0=20 L4 flags above with `--LF 16000 --KF 200 --D0 22 --KT 21` or `--LF 16000 --KF 200 --D0 24 --KT 23`. To investigate coarse L4 endpoints, account exactly for both the `h → L4` transferred source and the jump at the new endpoint boundary. Deleting tail-endpoint rows or substituting a mean would not suffice.

Local one-thread GMP benchmark commands:

```sh
gcc -O3 -march=native -ffp-contract=off -fopenmp diag_slack.c -lm -o diag_exact
gcc -O3 -march=native -fopenmp ../cert/cert_b3f_verify_gmp.c -lgmp -lm -o verify_gmp
bash run_diag.sh extensions
```

The unmodified verifier passed the locally built **D0=12 lower** endpoint `565655/1000000` with zero failures in **110.42 s** wall time and **647,424 KiB** peak RSS. The benchmark and projection in the measured section use the original B3-FullH verifier; multiply the projected time by two for both endpoints, excluding builds. T3 hull checks do more GMP operations than physical inner rows, so extrapolations from entry counts are planning estimates, not guaranteed costs. Pure-inner L5 rows have no Bernstein corners but involve additional transfer and budding arithmetic; no L5 GMP throughput was measured. The exact lower-only benchmark was for calibration; it did not certify the prototype or improve the known rigorous bracket.
