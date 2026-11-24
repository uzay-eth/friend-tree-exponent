# Independent mathematical audit: B3-FullH-20-L4

Scope: Theorem 1 and Corollary 2, the two hash-pinned certificates, and the verifier and candidate builder at `e40c01cb34bb5957d95281d7ab1241f172fa478a` (`exact/neon`). The mathematical enclosure and the literal production metadata are assessed separately from generic build/verification behavior.

## Findings

1. **Corollary 2 gives the wrong number of stored potential entries.** The eight arrays contain `1,398,838,433`, not `1,398,841,597` entries. The exact lengths in file order are `6,432,404; 273,529; 27,746,004; 13,880,592; 354,856,788; 540,925,385; 331,822,656; 122,901,075`. There are `22,687` rate/transfer entries, so the stated aggregate `1,398,860,120` and file size `96 + 8 × 1,398,860,120 = 11,190,881,056` are correct. In particular the stored edge array has `523²` entries; the verifier's `121693` counts parent-child tuples, not all stored edge entries. The discrepancy has no effect on drift or products.
2. **Generic builder pair staging can exceed capacity.** At `cert_b3f_build.c:1115–1117`, the 16,384-element Ybot array is indexed by an unchecked `npairs++`. The header allows `D0=DYB=KYB=40, KT=40, KF=41, LF=40` with otherwise valid settings. At `ksum=24`, the count is exactly `30,899`: sum over k of `(max(0,40−k)+max(0,39−k))*max(0,39−(24−k))`, restricted to `k≥2` and `24−k≥1`. This is a static count of the actual writes, not a claim that a full D0=40 allocation/build was run. The corresponding P4 maximum at `D=40` is `15,755`; the production Ybot/P4 cutoffs are small enough. The production verifier does not use these builder arrays.
3. **Hardening is partial.** The generic verifier still dereferences unchecked allocations in concrete-tree evaluation (`496–499, 635–636, 2033–2051`), comparison (`687–690, 725–745, 757–782`) and `dLf` (`1206`). This can make an out-of-memory run terminate abnormally, not produce a false pass. The native binary reader/GMP scalar interface also assumes little-endian LP64 with the expected 96-byte struct and lacks an ABI guard. On 32-bit `size_t`, the generic table length products and byte counts are not checked (the largest allowed Ytop table already exceeds `SIZE_MAX` entries); `long` must support both the file positions and 64-bit GMP scalar arguments. The audited host is x86-64 LP64. The verifier parses with `strtoll` and a full-consumption check but no `ERANGE` test; the observed oversized inputs still reject because all use sites impose tighter bounds. Override P/Q replaces the original fields before header validation (`940–942`), so an override run tests the replacement exponent, not the integrity of the original two fields. The shipped negative-control script and `run_cert.sh` default to four threads; the executed control suite was adapted to two without modifying it.
4. **Builder CLI is not range-safe in general.** The new `--cushion` uses `atof` without syntax or finiteness checks; auto exponent selection casts `ceil`/`floor` to `int64_t` without a range check and does not revalidate the finalized header. `--nit` is also unchecked (`nit+5` can overflow for `INT_MAX`); other numeric options still use `atoi`/`atoll`. Several builder allocations outside the main table buffers and the final `fclose` are unchecked. These affect generic candidate generation, which is outside the mathematical trust boundary, and were not triggered in the D0=10 smoke runs.

## Mathematical audit

- Rooting at the persistent initial core vertex fixes all parent indicators. An attachment at core `u` has rate `L_u + Σ_{v∼_K u}1/D_v`; attachments at its `L_u` leaves have combined rate `L_u/D_u`, keep `D_u` fixed, and create one terminal core `(1,0,1)`. The sum of rates is `n`, since each vertex's contribution `1/deg(v)` occurs `deg(v)` times; division by `n` is necessary in the one-step expectation.
- The transfers cancel incidence by incidence. For `f→g`, the tail control vector `[0,0,c2/(3D0²),c2/D0²+c3/D0³]` evaluates to `c2/D²+c3/D³`, and `[0,1/(3D0),2/(3D0),1/D0]` evaluates to `1/D`. For `g→hV,hC`, the donor's centered replacement is exactly the receiver's signed rate difference times the same derivative/transfer coefficient; this holds separately for the `T2` and `T3` weights. The L4 donor incidences are: two directed triples for `P4`; the sibling pair and each parent–middle–child triple for `Yt`; the downward path and sibling pair for `Yb`; and the three child pairs for `S3`. Only incidences whose added vertex yields an allowed motif are transferred. Incoming and outgoing transfers use the same pre-event cutoffs.
- On a promotion, only the promoting slot changes and a newly made core can occur only as an endpoint child. Deleting that terminal slot maps newborn `g` to `f`, newborn `hV`/`hC` to `g`, and newborn L4 terms to the indicated triples; this uses the post-bud state at the promoting slot. For `k=KT`, the restricted hV/hC slot disappears on promotion and no such restricted newborn term survives there. The uncut slot transitions to the overflow label; beyond it the uncut projected bud change vanishes, while the newborn f and g still contribute. An inner-to-tail gain lands at `T3` because the new degree is exactly `D0`.
- For leaf gain, put `t'=D0*t/(D0+t)` and expand the cubic as `F0−3at+3(a−b)t²+(−a+2b−p)t³`, with `a=F0−F1`, `b=F1−F2`, `p=F2−F3`. Using `(D+1)(t'−t)=−t` gives the `T1`, `T2`, and degree-2 remainder in §4.2. Multiplying by `1−(c+1−rho)t/(D0+t)` gives the nonoverflow gain; collecting the bud differences, inner-neighbour credits and the two transfer shifts gives §4.1. The degree-4 control ratios bound the rational remainder because its denominator controls are positive, so a value is a weighted average of the ratios. The affine rates attain extrema on corners. For overflow the factor `(D−c+rho)/(D+1)` lies in `[0,1]` for the physical `c`; the proof does not bound it by the capped index. The capped leaf-count checks distinguish `L=LF` (exact budding, bounded transfer) from `L≥LF+1`, where the bud factor and restricted Bernstein controls are bounded separately; VSTAR requires only the newborn terms. Tensor-product weights are nonnegative; their use does not assume stochastic independence of adjacent degrees. The `g(INF,INF)` row has no other contribution and is explicitly required to be zero.
- For comparison, charge `g,hV,P4` to the terminal vertex, each Ybot group to the bottom vertex (its parent/grandparent are unique), and `hC,Yt,S3` to the apex with pair/pair/triple multiplicities. The exact charge envelope includes zero for optional terms and restricted INF. The shift `Σ(k−h)/4` is zero, since both sums count every core edge once. On VSTAR using `KF+1` for the shift still gives a lower bound. `SΦ(T3)` is exactly the root `f(2,0,0)` entry. For lower drift the upper comparison `Φ≤C0 C` suffices even if Φ could be signed on some trees. Induction starts at deterministic `T3`, uses the positive factors `1+λ/j` from `j=3`, and applies comparison at the terminal step. `log ∏(1+λ/j)=λ log n+O(1)` by summability of the quadratic remainders. The mean exponent enclosure follows without an external bracket; the separately cited existence and almost-sure equality were not re-proved here.

The required prefactors are `146041684821/184219674916` and `4×146044513514/115454594604 = 584178054056/115454594604`. Subtracting `0.792758` from the lower prefactor gives exactly `7971742709/23027459364500000`; subtracting the upper prefactor from `5.059808` gives `213734969/901989020343750`. Both are positive. All product factors are positive for integer `n≥3`, including the empty product at `n=3`, so the same strict directions hold after multiplying. The interval is `[7071/12500, 70721/125000]` with width `11/125000 = 0.000088 < 10⁻⁴`.

As an independent check of the comparison sweep, I also adapted the prior referee's row-major reducer. It shares the format header but uses different traversal/aggregation code and does not call the production verifier. It reproduces all six integers (only the lower `SC0`, upper `4Sc0`, and both `SΦ(T3)` are given in Corollary 2; the other two are independently corroborated).

```sh
set -o pipefail
python3 -B round2/checks/b3f20_sol_compare.py | gcc -O3 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all -I. -x c - -o [local path]
OMP_NUM_THREADS=2 [local path] b3f20_lo.bin b3f20_hi.bin
```

```text
b3f20_lo.bin: 4Sc0=113472636536 SC0=184219674916 SPhi(T3)=146041684821
b3f20_hi.bin: 4Sc0=115454594604 SC0=200854054323 SPhi(T3)=146044513514
```

I also exercised the current GMP tail kernel under UBSan. The Python reference obtains the control coefficients by exact interpolation of the direct event formula, without using the verifier's expanded polynomial. The parameterized check covers `D0=20,40`, both certificate directions, overflow/nonoverflow, including nontrivial `t_max<1`, and both choices of transfer where applicable. Commands from `round2/checks`:

```sh
gcc -O2 -Wall -Wextra -fopenmp -fsanitize=undefined -fno-sanitize-recover=all b3f18_sol_tail_probe.c -lgmp -o [local path]
OMP_NUM_THREADS=2 python3 -B -c 'from pathlib import Path; p=Path("b3f18_sol_tail_probe.py"); s=p.read_text().replace("(10,18), (1,2)","(20,40), (1,2)").replace("[local path]","[local path]"); exec(compile(s,str(p),"exec"))'
```

```text
C/GMP tail probe: 864 exact T1/T2/T3 extrema agree with independently recovered Fraction controls
```

## Native arithmetic and single-buffer review

At the current LP64 caps, `NV≤30,061,004`, all state/rank and dense `g` indices fit `int`, and the packed offsets and file size fit `size_t`/`long`: a coarse bound using `NG<2000` and each fixed-slot domain `<2000` puts every length below `2000⁴` and the complete byte count below `8×8×2000⁴` plus the small header/auxiliary arrays. Native `D⁴` is at most `30501⁴ < 2⁶³`; centering sums are at most `65×10¹²` in magnitude and `2*S_r≤2×10¹²`. The Ybot pre-GMP value product is bounded by `38×10¹⁵`; short sums of motif values have ample headroom under `10¹⁵` per entry. Rate denominators in the two-neighbour sum are at most `40²`; the tail `x_den≤3200`, `|x_num|≤66×3200`, and native Bernstein prefactors and coefficient subexpressions fit signed 64-bit by large margins. Larger products with table/transfer values are made in GMP. Header validation precedes layout initialization; arrays are checked for complete reads and inclusive signed entry ranges, without taking the absolute value of `INT64_MIN`. The table reader and layout use allocation-result checks, with the exceptions above elsewhere. Seek/size errors and exact length mismatch reject before table evaluation.

In the builder the seven larger tables use one buffer each; `f`/`fn` remain distinct. References to the previous `g` (`Y_ref` and the f solve) are consumed before `g` is cleared. The subsequent schedule solves `g`, the two triples and then all four L4 tables. Births into higher levels read the previous sweep's still-live higher buffers; transfers from lower levels read the newly solved buffers. Same-table genuine successors precede their sources in descending rank (gain, lower Bernstein order) or descending child count within a rank (bud); unchanged states are handled in their own task. Symmetric packed child slots are traversed canonically. The OpenMP worksharing barriers separate different rank sums; per-table normalization is pointwise. This preserves the old/new dependency structure of the two-buffer scheme on the audited profiles. UBSan is not a race detector; exact acceptance of the output does not depend on the builder being deterministic or correct in general.

## Reproduction

Commands run from the repository root; GCC 13.3.0, `uname -m` = `x86_64`, `getconf LONG_BIT` = `64`. The tracked certificate sources and proof match the pinned commit. The format and verifier also compare byte-identically to their versions at `12202c2`.

```sh
git status --short --branch
git rev-parse HEAD
git diff --exit-code e40c01cb34bb5957d95281d7ab1241f172fa478a -- round2/cert round2/proof-b3f-motif4.md
git diff --exit-code 12202c2 HEAD -- round2/cert/cert_b3f_format.h round2/cert/cert_b3f_verify_gmp.c
gcc -O3 -march=native -fopenmp -Wall -Wextra round2/cert/cert_b3f_verify_gmp.c -lgmp -o [local path]
gcc -O2 -fsanitize=undefined -fno-sanitize-recover=all -fopenmp -Wall -Wextra round2/cert/cert_b3f_verify_gmp.c -lgmp -o [local path]
gcc -O2 -fsanitize=undefined -fno-sanitize-recover=all -fopenmp -Wall -Wextra round2/cert/cert_b3f_build.c -lm -o [local path]
sha256sum b3f20_lo.bin b3f20_hi.bin
```

Both diffs and all three compiles returned zero without output. Hash output:

```text
eb93985b2c90a516b65e85934f4ac847779f42b3fdbfcd2af4dc07fd30c37469  b3f20_lo.bin
55dd6015eca85ca2fba6267d9444ec8807554c456999803d8f9abcc695c090e2  b3f20_hi.bin
```

```sh
export OMP_NUM_THREADS=2
[local path] b3f20_lo.bin --threads 2 --tree-check 30
[local path] b3f20_hi.bin --threads 2 --tree-check 30
[local path] b3f20_lo.bin --threads 2 --tree-check 0 --override-pq 565684 1000000
[local path] b3f20_hi.bin --threads 2 --tree-check 0 --override-pq 565764 1000000
python3 -Bu round2/checks/b3f20_sol_exact.py b3f20_lo.bin b3f20_hi.bin
bash round2/checks/b3f20_sol_smoke.sh
```

The D0=10 script gives the exact UBSan builder and verifier commands, uses `--threads 2`, rewrites only the control script's in-memory thread argument, and asserts that no control returns a UBSan diagnostic. It ran to completion in both modes. Full transcripts from the Python and smoke commands are in `b3f20_sol_python.log` and `b3f20_sol_smoke_complete.log`.

### Python exact checks

The independent decoder validates the complete binary layout and explicit little-endian integers. The shipped Fraction evaluator brute-forces physical attachments on 6 boundary and 30 random trees per side. The separate core-compressed evaluator checks 15 cases per side and verifies each applicable transfer cancellation and newborn incidence; for `n≤100` the compressed generator matches physical enumeration exactly. At high degree grouping identical leaf attachments retains their exact multiplicity. Cases include `k=19→20`, `k=20`, `k=200→201`, `k=202`, `L=16000`, `L=16001`, and a tail with `D=16204`. Each designed residual is strictly positive (lower) or negative (upper). Exact enumeration of the labeled growth law agrees with the product bounds for `n=3..8`, with `E[C_8]=4717423/2177280`.

```text
lo: {'D0': 20, 'KT': 19, 'LF': 16000, 'KF': 200, 'DP4': 12, 'KP4': 9, 'DYT': 11, 'KYT': 8, 'DYB': 12, 'KYB': 9, 'DS3': 10, 'KS3': 7, 'S': 100000000000, 'S_r': 1000000000, 'P': 565680, 'Q': 1000000}; potential lengths=[6432404, 273529, 27746004, 13880592, 354856788, 540925385, 331822656, 122901075]; potential sum=1398838433; total=1398860120
Python exact-Fraction tree check (b3f20_lo.bin): 36 trees checked, fails = 0
hi: {'D0': 20, 'KT': 19, 'LF': 16000, 'KF': 200, 'DP4': 12, 'KP4': 9, 'DYT': 11, 'KYT': 8, 'DYB': 12, 'KYB': 9, 'DS3': 10, 'KS3': 7, 'S': 100000000000, 'S_r': 1000000000, 'P': 565768, 'Q': 1000000}; potential lengths=[6432404, 273529, 27746004, 13880592, 354856788, 540925385, 331822656, 122901075]; potential sum=1398838433; total=1398860120
Python exact-Fraction tree check (b3f20_hi.bin): 36 trees checked, fails = 0
lo=7071/12500; hi=70721/125000; width=11/125000; strict decimal margins: 7971742709/23027459364500000, 213734969/901989020343750
n=3: 146041684821/184219674916 <= 1 <= 146044513514/28863648651
n=4: 2169741311385597/2302745936450000 <= 4/3 <= 32547553303986797/5411934122062500
n=5: 123829306382087406387/115137296822500000000 <= 19/12 <= 18575572169204648770637/2705967061031250000000
n=6: 8614928674308202949749977/7196081051406250000000000 <= 647/360 <= 4307805215043742482452114759/563743137714843750000000000
n=7: 235678603743049508096310120789/179902026285156250000000000000 <= 17189/8640 <= 3535506203895915373940582077121239/422807353286132812500000000000000
n=8: 3184051604940562147168020633305217/2248775328564453125000000000000000 <= 4717423/2177280 <= 477657494664949854765494460365310752617/52850919160766601562500000000000000000
```

### D0=10 UBSan and controls

```text
Wrote b3f20_sol_d10_m2.bin (52653152 bytes), exact P/Q = 565597/1000000 (0.5655970000)
  [COMPARISON] 4*S*c0 = 113027843780, S*C0 = 183278518595, S*Phi(T_3) = 146044432538 (S = 100000000000)
  [TREES] exact-rational GMP drift check on 19 concrete/boundary trees: fails = 0
=== VERDICT: PASS (total_fails = 0, exact P/Q = 565597/1000000) ===
Negative controls summary (b3f20_sol_d10_m2.bin): 23/23 rejected cleanly (failures=0)
Wrote b3f20_sol_d10_m1.bin (52653152 bytes), exact P/Q = 565859/1000000 (0.5658590000)
  [COMPARISON] 4*S*c0 = 115564374748, S*C0 = 194851617082, S*Phi(T_3) = 146049008124 (S = 100000000000)
  [TREES] exact-rational GMP drift check on 19 concrete/boundary trees: fails = 0
=== VERDICT: PASS (total_fails = 0, exact P/Q = 565859/1000000) ===
Negative controls summary (b3f20_sol_d10_m1.bin): 23/23 rejected cleanly (failures=0)
```

Both positive runs printed zero failures in each of the eight tables. Negative-control failure counts, lower/upper respectively: impossible exponent `3427371/4268299`; nonzero INF `1/1`; individual perturbations f `142/143`, g `81/82`, hV `81/82`, hC `101/102`, P4 `1/2`, Yt `1/2`, Yb `1/2`, S3 `1/2`, c2 `1416/1110`, c3 `262/138`, Yp `94/49`, Yw `74/6`. The remaining structural tests produced header errors `-1,-2,-6,-7,-4,-6,-7`, the expected `[f[2]]` magnitude rejection and exact-length mismatch. No runtime sanitizer diagnostic occurred. An initial combined smoke invocation reached the shell’s 120-second timeout after the lower run and part of the upper build; the complete repeat above supersedes it.

Additional CLI/seek queries using `[local path]`: `--threads 2junk` and `--threads 999999999999999999999999` both returned `2 Invalid --threads`; `--threads 2 --tree-check -1` returned `2 Invalid --tree-check`; `--threads 2 --override-pq 999999999999999999999999 1` returned `1 FAIL: invalid header (err=-4)`; a pipe on `/dev/stdin` with `--threads 2` returned `1 FAIL: fseek end failed`.

### Fresh production GMP transcripts

```text
Verifying b3f20_lo.bin | mode=LOWER | lambda = 565680/1000000 | D0=20 KT=19 LF=16000 KF=200 | threads=2
  [f ] checked 6432404 states: fails = 0
  [g ] checked 121693 states: fails = 0
  [hV] checked 27746004 states: fails = 0
  [hC] checked 13880592 states: fails = 0
  [P4] checked 354856788 states: fails = 0
  [Yt] checked 540925385 states: fails = 0
  [Yb] checked 331822656 states: fails = 0
  [S3] checked 122901075 states: fails = 0
  [COMPARISON] 4*S*c0 = 113472636536, S*C0 = 184219674916, S*Phi(T_3) = 146041684821 (S = 100000000000)
  [TREES] exact-rational GMP drift check on 34 concrete/boundary trees: fails = 0
=== VERDICT: PASS (total_fails = 0, exact P/Q = 565680/1000000) ===
exit=0
```

```text
Verifying b3f20_hi.bin | mode=UPPER | lambda = 565768/1000000 | D0=20 KT=19 LF=16000 KF=200 | threads=2
  [f ] checked 6432404 states: fails = 0
  [g ] checked 121693 states: fails = 0
  [hV] checked 27746004 states: fails = 0
  [hC] checked 13880592 states: fails = 0
  [P4] checked 354856788 states: fails = 0
  [Yt] checked 540925385 states: fails = 0
  [Yb] checked 331822656 states: fails = 0
  [S3] checked 122901075 states: fails = 0
  [COMPARISON] 4*S*c0 = 115454594604, S*C0 = 200854054323, S*Phi(T_3) = 146044513514 (S = 100000000000)
  [TREES] exact-rational GMP drift check on 34 concrete/boundary trees: fails = 0
=== VERDICT: PASS (total_fails = 0, exact P/Q = 565768/1000000) ===
exit=0
```

### Sharp inward override controls

The perturbations are `4/1,000,000` in the tightening direction. Rejection means the existing tables do not certify that exponent; it cannot exclude a different certificate at that exponent. `--tree-check 0` disables the supplemental concrete tests, so rejection is attributable to exact certificate-row inequalities.

```text
Verifying b3f20_lo.bin | mode=LOWER | lambda = 565684/1000000 | D0=20 KT=19 LF=16000 KF=200 | threads=2
  [f ] checked 6432404 states: fails = 7975924
  [g ] checked 121693 states: fails = 191
  [hV] checked 27746004 states: fails = 165510
  [hC] checked 13880592 states: fails = 33659
  [P4] checked 354856788 states: fails = 8406648
  [Yt] checked 540925385 states: fails = 6250660
  [Yb] checked 331822656 states: fails = 26154443
  [S3] checked 122901075 states: fails = 523990
  [COMPARISON] 4*S*c0 = 113472636536, S*C0 = 184219674916, S*Phi(T_3) = 146041684821 (S = 100000000000)
=== VERDICT: FAIL (total_fails = 49511025, exact P/Q = 565684/1000000) ===
exit=1
```

```text
Verifying b3f20_hi.bin | mode=UPPER | lambda = 565764/1000000 | D0=20 KT=19 LF=16000 KF=200 | threads=2
  [f ] checked 6432404 states: fails = 7141623
  [g ] checked 121693 states: fails = 203
  [hV] checked 27746004 states: fails = 311559
  [hC] checked 13880592 states: fails = 76840
  [P4] checked 354856788 states: fails = 16538952
  [Yt] checked 540925385 states: fails = 17335204
  [Yb] checked 331822656 states: fails = 42926360
  [S3] checked 122901075 states: fails = 1747823
  [COMPARISON] 4*S*c0 = 115454594604, S*C0 = 200854054323, S*Phi(T_3) = 146044513514 (S = 100000000000)
=== VERDICT: FAIL (total_fails = 86078564, exact P/Q = 565764/1000000) ===
exit=1
```

## Verdict scope

The hash-pinned D0=20 mathematical enclosure **`[0.565680, 0.565768]` is verified**. No gap was found in the production drift inequalities or finite-size product bounds. The overall audit covers Corollary 2 literally and the generic code as requested; the gaps are the incorrect stored-potential count, the Ybot staging-capacity overrun at accepted generic parameters, incomplete verifier allocation/ABI/CLI hardening (and the four-thread script default), and the unchecked builder CLI/conversion/allocation/final-close paths. See Findings 1–4 for exact locations and impact.

VERDICT: GAP
