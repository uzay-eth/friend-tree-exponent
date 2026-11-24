VERDICT: CORRECT

**Scope:** the two hash-pinned D0=18 certificates and the claimed Rung 3 enclosure for the mean growth exponent, with existence as cited in the submission. I found no error in the production drift bounds, comparison, or endpoint arithmetic. The unrestricted binary verifier has real native-overflow bugs; they do not occur on either production input. The theorem's generic verifier wording should be qualified accordingly.

## Mathematical audit

1. **Root and rates.** Root the induced core tree at a core vertex (the original core of `T_3` for the process). The core stays connected; its oriented parent indicator really is `h`. At a core `u`, gain has rate `L+rho`; all `L` promotion targets give the same state and each has rate `1/D`. Promotion replaces `(L,k,h)` by `(L−1,k+1,h)` and creates one core `(1,0,1)`; the promoting parent's `D` is unchanged. Summing gain and promotion rates gives `n`, also directly from the ordered neighbour sampling. The `1/n` in the conditional drift is required and present.
2. **Transfer dictionary.** The positive receiver term and negative donor term multiply the same signed centering difference. Expand `1/D_z` by Lemma 1, with weights summing to one. For inner `f` the coefficient is `Δ_L f`; for tail `f` it is exactly `c2/D²+c3/D³` from the `T2/T3` controls. For `g`, an extra parent or child at the updating endpoint gives `hV` or `hC`: inner uses `Δ_L g`, tail uses the weighted `T2` coefficient `(F0−F1)/D0` and weighted `T3` coefficient `Y_ref`. Incoming terms in the triples (C 1385–1403, 1539–1556) match each direction of the subtractions in `g`, including the parent rate on the upper endpoint and the two symmetric child terms in `hC`. For the four L4 shapes the donors per occurrence are:

   | Receiver | Donors and external vertices |
   |---|---|
   | `P4(a,b,c,d)` | `hV(a,b,c)` gains at `c` against `d`; `hV(b,c,d)` gains at `b` against `a` |
   | `Ytop(p,u;w1,w2)` | `hC(u;w1,w2)` gains at `u` against `p`; each `hV(p,u,wi)` gains at `u` against the other child |
   | `Ybot(u,w2;w1,c)` | `hV(u,w2,c)` gains at `u` against `w1`; `hC(u;w1,w2)` gains at `w2` against `c` |
   | `S3(u;w1,w2,w3)` | each of the three `hC` pairs gains at `u` against its omitted child |

   These are bijections between eligible donor/external-neighbour incidences and the indicated receiver incidences. The donor centers precisely those extra neighbours for which the receiver's strict inner/child-count cutoffs hold, at the pre-event state. Other external inverse rates use `[0,1/2]`. Tail slots may still appear as uncut endpoints in any L4 receiver; this does not center that endpoint's own gain. The argument uses physical distinct vertices; coincident *basis states* for different vertices are permitted by the packed symmetric tables and do not eliminate donor incidences.
3. **Births.** A new core vertex has no core children. In every newborn non-singleton motif it is therefore a terminal child, and deleting it leaves a unique existing lower motif. `f(u)` gets newborn `f` and `g(u,nb)`; existing `g` gets the newborn `hV` (on the incoming edge) or one newborn `hC` per existing child edge. `hV` gets newborn `P4` at its bottom, `Ytop` at its middle, or `Ybot` with the newborn as the sibling at its top; `hC` gets newborn `S3` at its apex or `Ybot` with the newborn below either fork child. L4 has no higher term. Inner cutoffs for the newborn occurrence use the **post-bud** state at the promoting vertex; unmodified slots retain the old state. For tail promotion the `INF` weight is unchanged; its `g` bud difference there vanishes and it cannot carry a newborn triple's designated middle/apex. The `T1/T2/T3` controls take successor `k+1` with the same weights, exactly justifying the added newborn terms in `Fjb` (C 1194–1209, 1284–1286, 1300–1302). At `k=KT`, the bud degree still holds fixed; `g` and other uncut slots go to the overflow label; existing `hV/hC` with the changed restricted slot go to zero (`ua_bud`, `ca_bud`, and `F1b/F2b/F3b`); no newborn restricted triple at that slot exists. At `k>KT` the uncut existing motif's projected bud change is zero, but newborn `f` and `g` still appear in `f`. Inner-to-tail gain goes to the exact `T3` basis point, not to INF; restricted L4 slots exit to zero at their own degree cutoff.
4. **Bernstein tail algebra.** The fresh Python check evaluates the direct gain, bud, neighbour credit, and signed transfer separately, and compares with the C coefficient blocks reduced by their positive common multiplier `M=D0*x_den`. In the nonoverflow case it tests each unit direction of `(F, Fb[1:4], d, Y_ref)`, all production `c=0..18` and `n_in=0,1,2` with `n_in<=c`, both transfer settings, and distinct exact rational `t` values including `t_max`. The symbolic structure is linear in motif/credit values and in the centering or inverse-rate parameters separately, and rational in `t`; after clearing `(D0+t)^3`, the residual is polynomial of degree at most seven. The `F3` neighbour credit in the C `A1` block is `F3_up_sum−n_in*F3`; for restricted centers `F0`, its bud control, and `d0` are zero. Subtracting the actual-to-centered `T2` term replaces `rho` in `T2` by `rho_in+Rbar`; the `T3` numerator still uses the actual rates, and the separate shift is `(Rbar−rho_ext)Y_ref`. The degree-4 control formula independently agrees with C's scaled `c0..c4` and each denominator control is positive. The rational bound is a weighted average of the ratios with weights proportional to the positive denominator controls; the maximum/minimum over rate corners holds by affinity. A type-`INF` slot contributes no own rate term at that basis row, but can still receive sources and contributions from other slots. Skipping `g(INF,INF)` is valid since the separately checked value is zero, neither slot has a transition or receives a transfer there, and `F0b=F0`.
5. **Overflow and caps.** For `k>KT`, physical `c` can exceed its representational cap by any amount. `0<=D−c+rho<=D−c/2<D+1` uses physical `c`, so replacing the multiplicative `alpha` by `[0,1]` requires no phantom bound on physical `c`. Exact re-expansion of `(D+1)ΔB3` gives the C `T1/T2` terms and the degree-2 numerator; its control ratios and the zero corner cover the gain, while the remaining `T3` credit is affine in `t`. At `L=LF` the new `f` index really has `L−1`; `b*B3(g)` uses integer denominator `D⁴`. For `L>LF` and `k<=KF`, the new `f` index is at `LF` (or VSTAR on a `KF` exit); `b∈[(LF+1)/(LF+1+c),1]`. Both cases use a valid enclosure of the `f` transfer from the cubic controls `{0,0,c2*x0²/3,c2*x0²+c3*x0³}` and the two possible clipping extremes for `rho`. `VSTAR` has no change to its old `f` on either event, no applicable transfer, and uses `b∈[0,1]`; the newborn edge uses its correct tail cap and `t<=D0/(KF+1+h)`. The four restricted interval controls are the left de Casteljau boundary; the verifier checks the required corners, including both signs of `f` differences.
6. **State coverage.** Parents use `k>=1` or INF; children use `h=1` or INF; the hV/hC restricted slots use all eligible inner and tail states, with INF as a zero control. The L4 loops cover exactly the fixed inner slots and the corresponding unrestricted endpoints. The omitted `D<2` one-body states are noncore, nonroot `k=0` eligible one-body states still appear, and caps do not lose any actual `D` information where a transfer or tail evaluation needs it. For any single slot, all other slot controls may be fixed before grouping its drift and inner-neighbour credit. Interpolating the verified rows with products of nonnegative slot weights therefore bounds the sum even when several slots are tail and their physical rates are correlated. When an updating inner slot sees a tail inner-motif neighbour, the corresponding credit is instead on that tail slot. There is no double counting and no independence assumption about degrees.
7. **Comparison and induction.** `g,hV,P4` charge once to the bottom core vertex, and each `Ybot` bottom has a unique parent and grandparent and hence at most a unique group of `k_apex−1` siblings. `hC,Ytop,S3` charge to their apex, with exact binomial multiplicities. Missing optional bottom terms and restricted INF terms include zero in their bounds. Convexity bounds tail slots; saturated states retain the same charge envelope and `k-h` for the lower comparison, and VSTAR uses the smallest physical `k` in its class. Each oriented edge contributes `+1/4` at its parent and `−1/4` at its child, so the lower charge shift sums to zero. `T_3` has only the root `f(2,0,0)`. The GMP constants independently reduce to the same integers (below). Multiply the conditional expectation inequality by no random factor, take expectation, and iterate the deterministic positive factors `1+lambda/j` from `j=3`; do the core comparison at time `n`. For lower, `Phi<=C0*C` is sufficient even if pointwise `Phi` is signed elsewhere. `log(product)=lambda log n+O(1)` by summable `O(j⁻²)` remainder, giving the claimed liminf/limsup bounds and thus the enclosure for the stated existing mean exponent. The almost-sure identification is a separate cited theorem and not needed in this inference.

For a compact algebraic cross-check covering all parameters: let `a=F0−F1`, `b=F1−F2`, `p=F2−F3`, `d=D0`, `z=dt/(d+t)`, `s=1−t`, `x=c+1−rho`. The cubic in powers of `t` is `F0−3at+3(a−b)t²+(−a+2b−p)t³`. Since `(D+1)(z−t)=−t`, subtract the `3ts²*a` and `3t²s*2b` terms from the scaled gain; the remainder is `t³*P2/(d+t)²`, with

```
P2 = (3d²p+3d(a−b)) + (3dp+(3−3d)a+(6d−3)b)t + (p−2a+4b)t².
```

Multiply this whole three-term gain identity by `1−xt/(d+t)` for the nonoverflow rate. After keeping `a` on `T1` and `2b−xa/d` on `T2`, the remaining `T3` gain numerator over `(d+t)³` is

```
(d+t−xt)*P2 + s*x*(3a(d+1)/d−6b)*(d+t)².
```

With `delta_j=Fb_j−F_j` (and `delta_0=0`), the bud coefficients are `delta_1`, `delta_2−c*delta_1/d`, and `delta_3*(1−ct/d)+3c*s*(delta_1−delta_2)/d`. Add the neighbour-credit coefficients from §4.2 (valid in both branches) and the signed T2/T3 transfers. Expanding the displayed numerator gives exactly C's three `d1z`, `d2z`, and `f23` blocks; the bud/credit polynomial times `(d+t)³` gives the `A0,A1` block. This establishes the rational rearrangement algebraically in addition to the exact evaluation tests.

## Machine arithmetic and concurrency

I reviewed the whole verifier (2051 lines) and format header. The following are conservative maximum absolute bounds for the native expressions **on both supplied files**, not limits enforced by the parser.

| C / header location | Audit of arithmetic before GMP |
|---|---|
| Header 81–155, 173–351; verifier 663–669, 891–945 | x86-64 Linux LP64, so signed/unsigned `long` accepts these `int64_t`/`uint64_t` values; binary is little endian with 96-byte header. `NG=435`, `NG²=189225`, `NV=5,068,364`; the maximum single table is `175,910,700` entries. All dense `int` indices and `size_t` packed pair/triple offsets, length sums and byte counts fit. Zero-row triangular arithmetic casts `(i1-1)` to unsigned only in a product with `i1=0` in that case. Verified packed dimensions: inner/parents/children/HU/CU = `320/397/210/187/352`; L4 fixed dimensions = `35/20/60/35/24` (P4/Yt/Ybu/Ybw/S3). |
| 59–143 | Positive denominators; signed negation at 77 is safe on all production callers (no `INT64_MIN`). Cross-products and centered source numerators/products are GMP. In rate-sum arithmetic each input inverse denominator is at most 18; two-neighbour denominator ≤324 and numerator ≤36. |
| 180–253, 254–406 | `c<=18` in this branch, ext bounds are nonnegative, `x_den<=648`, `abs(x_num)<=19*648`, `M_s<=11664`, `abs(Mx)<=221616`. Native rho preproducts and both `nrin` cross-products are small. `D0³=5832`; `tn=td=1` on production reachable nonoverflow rows, so the actual `c0..c4` are `12,3,2,3,12`. Even allowing `tn,td<=18` gives all these factors and the scaled native multiples comfortably below `8e6`. A triangle bound on the native coefficient expressions in 347–373 is `513332640`, including `3*D02*(Mx+M_s)` and its larger sixfold relative; large multiplications by motif values and by `M` where appropriate are GMP. `2*S_r=2000000000`; `Y_ref` is multiplied into the shift in GMP. |
| 408–454 | Overflow coefficients use at most `D0²`, `6*D0`, and degree-2 denominator controls; arithmetic on motif values and the neighbour-credit numerator is GMP. No signed/unsigned inversion in selected extrema. |
| 465–660, 1997–2050 | Concrete GMP tree checks have small degrees/indices (default seeded sizes below 60); unsigned xorshift is defined and modulo divisor positive. `P,Q` and all scales are positive, `-P` is safe. Reduction counts fit `long long`; frees happen outside parallel regions. |
| 684–878 | Full scan gave `max\|qYb\|=2871334141`, hence `abs(mult*qYb)<=17228004846` at 707. Pair/triple native multiplicities are ≤136/20 on actual hC/S3 production loops (other pair domains smaller). Charge sums, multiplication by `S`, and the zero-sum adjustment use GMP. |
| 959–1179 | For finite rows `Dv<=14181`, so `D⁴<=40441516074779121`, below both signed and unsigned 64-bit limits. Capped interpolation denominator `dm<=14182`, `dm³<=2852413236568`; `s_n³` and `sn³` are no more than that, as are mixed binomial numerators. `2*D³` and `3*D²` fit; `S_r*D³` can exceed native width but is a GMP multiplication. Rate magnitudes from the full scan are ≤430555556, and every integer number of centered neighbours is ≤18, so weighted centering sums are ≤7750000008. Direct `df_up` and `dLf` differences are in range. |
| 1181–1991 | Every native signed addition/subtraction of table values combines at most four signed entries (including two inner neighbour deltas); each potential-table value has magnitude at most `167836722289`. Newborn additions to `Fb`, residual offsets and signs are therefore safe by a wide margin. Centered child-count multipliers are bounded as above. Other products of value/time or value/rate operate in GMP. |

OpenMP: read-only layout and certificate tables, private per-worker `GmpWork` (also private GMP heap objects), local temporaries, reduction counters, barriers before reporting/freeing; no shared mutable scratch or observable data race found. Native comparisons are on signed values except explicit nonnegative scales, counts, denominators and GMP `_ui` arguments. `zf_add_z` scratch does not alias a caller's numerator/denominator. The UBSan production comparison scan and smoke verifier runs supplement this audit; the production GMP reruns use the normal optimized build.

## Defects and clarifications

- **Generic-input overflow, reproduced.** Setting only `S_r=INT64_MAX` on the valid small certificate passes header validation and reaches signed overflow at C 1024 (`2*S_r`; also the same expression at 388 and 1029). Setting `LF=INT32_MAX` reaches signed overflow at header 175 before the file-size check. Header 165's source expression `KT+1` is also not intrinsically range-safe, although both of my GCC UBSan builds rejected a planted `KT=INT32_MAX` with `err=-7` and no observed sanitizer diagnostic. Other unguarded generic hazards include `KF+1`, `2*NKF*NLF`, layout size/power products, rate/table differences, `mult*qYb`, unsigned-to-signed narrowing, and negating a putative `INT64_MIN` in `zf_add_si`. Unchecked `atoi`/`atoll` inputs, allocation results in layout/comparison/tree code, and `fseek`/`ftell` are additional generic robustness issues. No false acceptance was demonstrated. Generic exactness should require machine-range validation or moving these operations to checked arithmetic/GMP; this finding does not touch the verified hash-pinned range.
- **Wrong negative-test offset.** `negative_controls.py:93` mutates byte 52 (`KF`), not byte 60 (`KT`). The reported `err=-6` is honest output for the wrong test. An actual D0=10 `KT=8` mutation at 60 returned `err=-7`. Also the driver treats every nonzero return, even a crash or OOM, as a clean rejection; the outputs for the 36 runs here show diagnostic rejections, not crashes.
- **Shipped Python docstring overclaims.** `tree_identity_check.py:9–10` advertises a separate telescoping verification; its implementation only evaluates the full potential and brute-force attachment drift. It does that with exact `Fraction`; the floating relative residual is only for display. The new tree script supplies a separate transfer ledger and event decomposition.
- **Proof text.** §3's unchanged Bernstein weights on bud do not imply a surviving newborn triple at `k=KT`, where the successor is out of domain. §4.2's zero bud change is for an existing accepted motif's uncut overflow slot; there are still newborn `f/g` terms. §4.1 should define `rho_T2` as either the centered or the actual rate according to transfer use, state that `F0b=F0`, and give the explicit `T3` numerator or cite the coefficient derivation. §4.3's `L=LF exact` means the bud calculation is exact, while the independent transfer correction is bounded. The “every rooted finite tree” statement needs the convention of rooting **the core** at a core vertex as in §2, not applying the slot convention to a root which is outside the core. The entry count in §1 is rounded: `436,870,169` stored potential-table entries, `436,886,484` including rate and transfer arrays, `3,495,091,968` bytes including header. The verifier's `checked` labels count finite layouts/tuples, not all tested corners, and include three invalid f slots and the special `g(INF,INF)` tuple.

## Commands and evidence

Commands below are from the repository root; temporary stdout logs and binaries are under `[local path]`. I read the proof, the entire header/verifier, both supplied Python files, and `round2/outputs/cert_b3f_18_verified.log`.

```sh
gcc -O3 -std=c11 -Wall -Wextra -fopenmp round2/cert/cert_b3f_verify_gmp.c -lgmp -o [local path]
gcc -O2 -g -std=c11 -fopenmp -fsanitize=undefined -fno-sanitize-recover=undefined round2/cert/cert_b3f_verify_gmp.c -lgmp -o [local path]
sha256sum b3f18_{lo,hi}.bin
for side in lo hi; do [local path] b3f18_$side.bin --threads 2; done
for side in lo hi; do
  [local path] smoke_$side.bin --threads 2
  TMPDIR=[local path] python3 -B round2/cert/negative_controls.py [local path] smoke_$side.bin
done
for side in lo hi; do python3 -B round2/cert/tree_identity_check.py b3f18_$side.bin; done
gcc -O3 -Wall -Wextra round2/checks/b3f18_sol_scan.c -o [local path]
[local path] b3f18_{lo,hi}.bin > sol-b3f-scan.log
gcc -O3 -std=c11 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=undefined round2/checks/b3f18_sol_comparison.c -o [local path]
[local path] b3f18_{lo,hi}.bin
python3 -B round2/checks/b3f18_sol_algebra.py
gcc -O2 -std=c11 -Wall -Wextra -fopenmp -fsanitize=undefined -fno-sanitize-recover=undefined round2/checks/b3f18_sol_tail_probe.c -lgmp -o [local path]
python3 -B round2/checks/b3f18_sol_tail_probe.py
for side in lo hi; do python3 -B round2/checks/b3f18_sol_trees.py b3f18_$side.bin; done
python3 -B round2/checks/b3f18_sol_negative.py
gcc -O0 -g -std=c11 -fopenmp -fsanitize=undefined -fno-sanitize-recover=undefined round2/cert/cert_b3f_verify_gmp.c -lgmp -o [local path]
[local path] sol-huge-kt.bin --threads 2 --tree-check 0
```

GCC 13.3.0; builds completed without compiler diagnostics. Source SHA-256: `f2a1b91575a82e0039b40ac1eb247e86291c803219eb6176306393c70475796b` (C), `05055d3cb789c3ae83ca20321b7e7def5f76b30bf2960f7f60f03c3947cd3729` (header). Both production input digests match the submission:

```
9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec  b3f18_lo.bin
1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41  b3f18_hi.bin
```

**Production re-execution:** see `sol-b3f-prod.log` (two threads, default 25 random plus four hand-built trees per side).

**Independent comparison output:**

```
b3f18_lo.bin: 4Sc0=113465364604 SC0=184124847315 SPhi(T3)=146042017510
b3f18_hi.bin: 4Sc0=115462837236 SC0=195820675182 SPhi(T3)=146044806770
```

The Fraction script established the exact prefactors `146042017510/184124847315 > 0.793168` and `584179227080/115462837236 < 5.059457`. Bounds are `565673/1000000` and `565775/1000000=22631/40000`; difference `102/1000000`.

**UBSan smoke and negative control output** (`sol-b3f-smoke.log`): both smoke certs pass all drift tables and comparison plus 29 concrete trees per side, with no sanitizer diagnostic; exact P/Q are `565597/1000000` and `565859/1000000`. The full negative suite reports `18/18 rejected cleanly (failures=0)` each. Per-control drift-failure counts (lower/upper, in suite order after the four structural controls): nonzero INF `1/1`, f `142/143`, g `81/82`, hV `81/82`, hC `101/102`, P4 `1/2`, Yt `1/2`, Yb `1/2`, S3 `1/2`, c2 `1416/1110`, c3 `262/138`, Yp `94/49`, Yw `74/6`. Override failures: `3427371` at `571252/1000000`; `4268299` at `560200/1000000`. Structural diagnostics: `err=-1`, `err=-2`, `err=-6` (KF mutation), and `file size 52653144 != expected 52653152`. Additional tests: true KT `err=-7`; huge `S_r` and `LF` produce the UBSan errors described above; huge KT `err=-7` at `-O2` and `-O0`.

**Fraction checks:** shipped script `21 trees checked, fails = 0` on each side (`sol-b3f-trees.log`). Fresh algebra output:

```
tail coefficients and hulls: 175104 exact nonoverflow basis identities; overflow and inverse-degree credit passed
invD, f tail transfers, cap/VSTAR Bernstein controls, endpoint width and finite-size constants passed
C/GMP tail probe: 864 exact T1/T2/T3 extrema agree with independently recovered Fraction controls
```

The UBSan tail probe calls the original `add_tail_gmp` directly and uses Newton interpolation from the direct rational generator to independently recover the polynomials, then checks the actual C return fractions in both modes. Synthetic `KT>D0−1` cases also hit the nonunit `t_max` branch absent in the production nonoverflow rows. There were no sanitizer diagnostics.

Fresh tree script (`sol-b3f-fresh-final.log`): 15 per side, including root/nonroot `KT=17` exits with positive L, overflow `k=18`, `KF=180` exit, `k=181,182`, `L=LF`, `L>LF`, and `k>KF` with `L>LF` simultaneously (up to size 14,368). It independently assembles all eight motif values using only the existing Python byte decoder. It checks the vertex charge and edge re-charge identity, checks sender/receiver contributions separately (tensor vs physical formula), assigns newborn drift to the deleted-newborn motif, and compares to physical per-target drift on small trees and exact multiplicity-grouped attachments on large trees. All declared transfer and birth categories occur, including 558 tail `g->hC` and 46 tail `g->hV` transfer incidences per side. The exact measured production drift residual has the expected strict sign on all cases. The separate labelled evolution sums exact probabilities through `n=8` and checks the finite-size product bounds, yielding `E[C_8]=4717423/2177280`.

Administrative checks: `free -h`, `nproc`, `gcc --version`, `df -h [local path]`, `uname -m`, and process/log polling. A NumPy availability probe returned `ModuleNotFoundError`; the scans used C and the algebra used Python's standard library. While testing the additional negatives I initially expected a diagnostic for enormous KT; the actual compiler returned the intended rejection, as reported above. I re-ran the fresh tree and algebra scripts after strengthening their tests.
