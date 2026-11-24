VERDICT: CORRECT

Scope: cold audit of Theorem 2 and Corollary 4 for the hash-pinned `D0=20, Dc=7, B=2, mask=3` pair at `25ce6869abfc5d5379438f5a188829aa361cc565` (`exact/neon`). Both fresh production GMP replays passed. The pointwise drift and comparison yield `5656929/10000000 <= mu <= 565757/1000000` and, for every `n>=3`,

```
(8113447325/10248982441) * prod_{j=3}^{n-1}(1 + (5656929/10000000)/j)
    <= E[C_n] <=
(146044250013/22074015170) * prod_{j=3}^{n-1}(1 + (565757/1000000)/j).
```

## Proof-to-code audit

All line references below are to the audited commit. `format` means `round2/cert_l5/cert_b3f_l5_format.h`; `rows` means `round2/cert_l5/l5_verify_rows.h`; `GMP` means `round2/cert_l5/cert_b3f_l5_verify_gmp.c`.

1. **`proof-b3f-l5.md:7–36`, cutoff and storage.** `format:23–35` validates base v1 magic/mode/scales/caps and subsequently the actual version, extension, `Dc`, buckets, mask and every displayed L4 cutoff. `KT >= D0-1` also eliminates exact-state overflow. For `D<Dc`, `k<=D-h<=Dc-1`. This includes `L=0,k=Dc-1,h=0` in every suitable L4 apex domain, including S3. For a promotion the post-state still satisfies `k+1<=D-h` because the pre-state had `L>=1`. The four L4 shape/state constructors and L5 dispatcher agree: chain `(a,b,c,d)`, unary-with-fork `(p,u,w1,w2)`, binary with one nontrivial branch `(u,w1,w2,c)`, and 3-child apex `(u,w1,w2,w3)`. `format:120–158,182–244` derives the 17 catalog shapes (9 of size five) by unordered child-shape multisets, uses arity-filtered exact states with physical `h`, and admits coarse states only at motif-degree-one slots. Colex multiset rank is `sum_j choose(r_j+j,j+1)` for ascending zero-based child-row indices `r_j`; disjoint mixed-radix blocks plus the apex-`h` offset give every row once. Coarse apex uses the `h=0` block even if its physical `h=1`. Exact row count depends on `(Dc,B)`, not `D0` or the mask.
2. **`:40–52`, exact event and newborn ownership.** Source rings are uniform at total rate `n`; to attach at `u` the generator rate is the sum of neighbor reciprocals, not `1/D_u`. Existing motif support can only shrink on leaf gain; promotion cannot re-support an old motif since all old physical degrees are unchanged. A new supported 5-subtree contains the promoted old leaf `ell`, whose unique attachment is to `u`; deleting `ell` recovers its unique L4 owner. `u` gains a motif edge, so it is interior in the newborn and cannot have `D_u>=Dc`. `format:245–268` requires pre-L4 projected support, an exact promoter and `L_u>0`, buds it and inserts `(1,0,1)`. The four GMP insertions at `1783–1785`, `1875–1877`, `1975–1977`, `2065–2067` add `q5` to the bud target before subtracting the old L4 value and multiplying by `L_u/D_u`. No mask condition on births. All spectator controls carry the same projected `q5` value, including the INF or overflow control for a tail spectator; tensor weights sum to one. The code never puts an L5 birth in a coarse promoter's Bernstein own term.
3. **`:54–67`, zero-sum split and fibers.** For each selected exact receiver in a projected-support L4 row, `format:262–267` replaces the external denominator `2` with `Dc`. The external counts in the four GMP blocks are physical core degree minus motif degree, per receiver; all inner, unsupported and tail cases retain the appropriate original bounds where not transferred. Per external neighbor `z`, the clamped part is exactly `min(1/D_z,1/Dc)` and the positive residual exists iff `D_z<Dc`. A selected exact `u` stays eligible as an interior of `m+z`; conversely `rows:73–82` takes each exact endpoint `z`, deletes it, and tests the receiver in the resulting **4**-motif. Deleting a unary apex re-roots at its unique child without changing that child's physical `h`; original parenthood is preserved in all other cases. This is a bijection on transfer events, with multiple deletions included where applicable. The inverse takes the paired differences `L4(after)-L4(before)` for each complete tuple of unchanged spectator controls. `format:164–176` filters by true L4 apex/child domains and by bucket: `BIG_INNER` is exactly the inner controls with `Dc<=D<D0` in the slot; `TAIL` is INF and T1/T2/T3 with the right role and includes `KT+1`. The apex role permits both physical `h` values; a child permits `h=1` (INF is role-agnostic). For a hypothetical cutoff miss `format:279–311` returns zero. `rows:14–50` takes min/max of the *differences* before multiplication by the positive residual, not independent min/max of the two endpoints. Its cache stores a checked native integer and verifies the complete ordered key (marked receiver, all ordered parent slots, all endpoint states). In the proof's displayed line 66 the middle term should be `S*Delta Phi_4` when the outer terms are `Delta I_4`; code uniformly uses scaled units. This is a notation correction.
4. **`:69–85`, L5 own generator.** `rows:83–113` handles exact leaf gain with the old `L` and exact internal-neighbor rates, the sign-dependent bound on exterior rates, and sign-dependent coarse-neighbor interval (`1/Dc` or `1/D0`); an exact endpoint crossing `Dc` changes bucket while an interior crossing loses its contribution. Exact promotion at `L/D` preserves degree and updates the label. A TAIL endpoint has no own jump. A BIG_INNER endpoint can only jump on gain at `D0-1`; since `r=L+sum_core(1/D_w) <= L+(k+h)/2 <= D0-1`, `rows:85–91` correctly takes the adverse member of `{0,(D0-1)*delta}`. The zero corner covers all interior-of-bucket degrees. `rows:119–139` uses the positive denominator `Q*lcm(1,...,D0)`, initializes with `-P*lcm*v`, and adds the two reciprocal residuals separately. Summing the direction-selected row enclosures is valid even when the chosen corners cannot simultaneously be realized on a concrete tree.
5. **`:89–111`, comparison.** An apex with at least two motif children charges itself, otherwise its unique child is an exact interior and is charged. Let `b(D)` denote the maximum count of rooted connected size-five sets in the full rooted ambient tree with `D` root branches and `Dc-2` forward branches thereafter. Induction on child classes gives the claimed unordered-choice `N(R,m)`: choose each class from remaining branches, multiply recursively for each chosen branch, and use the upper bound `Dc-2` for any non-root *interior*. The unary-apex case uses one predetermined parent and so consumes only the `h` indicator. All shape capacities are simultaneously realized as counts in a filled orientation-specific ambient tree (the upward branch contributes only unary-apex cases). Ignoring orientation, filling to `D` root branches and allowing the upward branch a full forward tree can only increase that count. With `F=x(1+F)^(Dc-2)`, Lagrange inversion gives `[x^4](1+F)^D=(D/4)[t^3](1+t)^(D-1+4(Dc-2))`. The code checks the summed *capacities*, not just observed counts (`GMP:923–937`). Each state/shape extrema includes zero and can safely be multiplied by its capacity regardless of coefficient sign. `GMP:747–854` charges the older tables as in `proof-b3f-motif4.md:149–161`; in particular `Ybot` charges its bottom descendant, which fixes the grandparent and middle and leaves `k_grandparent-1` sibling choices. Comparing independent tail-child pieces against their separate extrema is valid because each uses the same convex child weights; keeping the pieces together yields the sharper inequality `max_j sum_m U_m(j) <= sum_m max_j U_m(j)`. The joint and separate values independently reproduce the exact Upper difference of `134398139`. `GMP:938–967` applies `+(k-h)*S` to the lower charge after multiplication by 4; the exact sum over core edges is zero and substituting `KF+1` in the overflow class is downward. No recharge is used on the upper side. Capped `L` classes in this sweep are valid: beyond `LF` the potential's one-body value is constant while the other coefficients are bounded by the same tail hull, independent of `L`. `T3` has only the root `f(2,0,0)` term.
6. **`:115–123`, implication.** Under the v1 generator decomposition, the only changes in inherited row blocks are the exact L4-to-L5 transfer and L5 newborn terms checked above. The enlarged sum bounds `L Phi - (P/Q) Phi` pointwise. Conditioning at each `n` gives multiplicative factor `1+P/(Qn)>0`; iteration from the deterministic rooted `T3` and the correct comparison side give the exact products for all `n>=3`. Their logarithms divided by `log n` tend to `P/Q` by `log(1+P/(Qj))=P/(Qj)+O(j^-2)`. This bounds liminf and limsup directly; `proofs/existence.md` identifies their common value as the stated mean exponent. No prior numerical bracket for that exponent enters any row check or comparison. At `GMP:975–981` positive `c0` is mandatory only on Upper; Lower needs positive `C0` and `Phi(T3)` (and the verifier also requires the latter two on Upper). Theorem 2's parenthetical suggests the binary enforces strict positivity of all three for both modes; amend it. Both pinned payloads have all three strictly positive.
7. **`:137–161`, counts and fractions.** Independent wire parser rebuilds all table counts; `npars=481`, `nchs=253`, `NG=523`, `INNER_N=396`, `NT=42`, `NHU=228`, `NCU=432`. The eight auxiliary arrays sum to `2*396+4*42+10626+10101=21687`. `1,421,247,831+21,687=1,421,269,518` int64 entries and `120+8*1,421,269,518=11,370,156,264` bytes. `1,421,095,995` is the sum of the verifier's printed state counts. The `f` print includes four `D<2` padding slots that are skipped; the `g` print includes INF/INF, checked separately by exact zero equality, not evaluated by the drift loop. Thus there are `1,421,095,991` actual valid/equality entries counted on that convention (`1,421,095,990` drift entries excluding INF/INF), and some `f` entries require multiple corners. This is a count-description issue, not a coverage issue. The prefactors reduce exactly from `146042051850/184481683938` and `4*146044250013/88296060680`; exponent difference is exactly `641/10000000`.

## Generic implementation observations

- Format parsing: LP64/little-endian is explicit, base and extension fields validated, and the computed row/byte count is cross-checked before payload allocation. `L5` allocation and cumulative sizes use saturating `size_t`; payload elements are rejected outside their native bound before use. For allowed headers `NG<=1993`, `Nv<=30,061,004`, each `NG*NG` index and `Nv` ID fit a signed 32-bit integer. `l5_source_key` needs at most 61 bits and each encoded state fits 13 bits. L5 exact states fit the `64*64` LUTs; every exact receiver strength/own delta product is bounded by `39*2*10^15`. The worst inherited native `D^4` is `30501^4=865478558581622001 < 2^64`; powers of cap endpoints and all rate-mask and native comparison multipliers have ample headroom (e.g. `63*10^15` for a crude weighted Ybot entry bound). GMP is used before the larger denominator and signed coefficient products. `l5_choose` multiplies before division, but repeated child classes of multiplicity two have subtrees of size at most two, and higher multiplicities only repeat leaves; these intermediates are below `1601^4`. Total L5 rows are crudely bounded by `18*1598^5`; total L4 and earlier bytes below `5.05*10^14`; all inherited table-size products and L5 rank arithmetic fit LP64. No arithmetic overflow affecting the accepted mathematical certificate was found.
- A valid generic header can still demand a very large amount of RAM/CPU; size validation is not a resource policy. File reading assumes a stable payload: it does not authenticate bytes itself. The production replay uses separate SHA-256 checks. Native `size_t` index helpers assume the caller supplies valid logical roles and state IDs, not arbitrary integers.
- A one-core size-only probe of the maximum listed header limits (`D0=Dc=40`, `KT=64`, `LF=30000`, `KF=500`, all L4 cutoff pairs `(40,64)`, `B=1`) completed and printed `rows=2512414592774190 file_bytes=20122593030429544` without allocating the payload. Minimum compatible shapes also counted successfully: `D0=Dc=4,B=1` gave `rows=5859 file_bytes=329192`; `D0=5,Dc=4,B=2` gave `rows=8260 file_bytes=722376` (both `LF=KF=D0`, `KT=D0-1`, all L4 pairs `(4,3)`). These are layout probes, not drift tests.
- Builder (`cert_b3f_l5_build.c:316–412,1351–1417`, and `l5_build_rows.h`): numeric CLI arguments are parsed and range-checked, `--size` returns without allocating tables, and quantization rejects nonfinite/out-of-range values before `llround`. It is a floating candidate finder. The L5 recursion increments degree/bucket or child count, never changes rooted shape, and one worker owns each size-five shape across both `h` blocks; later simultaneous `q5` rescaling has a barrier. No cross-shape dependency or data race in the L5 sweep was found. Generic input or fitting settings can fail to converge; a failed serialization can leave a partial output at the user-supplied path (`fopen(...,"wb")`), so a successful independent verifier run is the proof gate. Builder simulation rates and a starting numerical iterate are heuristic choices whose exact serialized results are freely audited; they are not assumed true physical rates or a prior bracket.
- `run_cert_l5.sh:51–57` treats any nonzero override exit as a successful negative control, including resource/environment/usage failure. Use a clean exit code plus the mathematical failure marker (as done by the independent controls). Its D24 option does not enforce production resource limits. The v2 GMP introductory comment still describes v1/eight tables. These are generic test/reporting issues, not failures of either pinned D20 file.

## Independent reduction and physical-tree replay

`b3f20_l5_sol_check.py` imports no project checker. It creates the rooted shape catalog by enumerating unordered forests; opens a read-only map; uses a standalone shape-agnostic C integer extrema scanner over sequential slices, with addressing and all charge formulas in Python; and evaluates each concrete potential by connected induced subset enumeration and literal graph attachments with `Fraction` weights. Sibling attachment events at a common core parent may be combined because the old leaves are exchangeable and the resulting rooted trees are isomorphic. It independently counts actual per-shape charges, checks births against pre-event 4-core sets, and audits the physical own-crossing corner. It also verifies individual clamped/residual splits, matches L4 debits and L5 credits edge by edge, and explicitly enumerates the raw paired-jump control extrema for sampled one-coarse-spectator fibers in each bucket and role. The 52 tests per side include stars at `D0-1`, `D0`, `D0+1`, `109`; endpoint and interior cutoffs on paths/forks/sub-stars; all nine rooted motif shapes; both physical apex `h` roles; both root and non-root hubs with `k=KT+2`; `L=LF-1,LF,LF+1`; and seeded friend trees.

Command and observed output (run from repository root):

```text
$ taskset -c 15 python3 -B round2/checks/b3f20_l5_sol_check.py b3f20_l5_lo.bin b3f20_l5_hi.bin
b3f20_l5_lo.bin: mode=2 P/Q=5656929/10000000 counts=(6432404, 273529, 27746004, 13880592, 354856788, 540925385, 331822656, 131094480, 14215993) aux=21687 bytes=11370156264 logged_states=1421095995
4Sc0=76840591160 SC0=184481683938 SPhi(T3)=146042051850 joint_4Sc0=76840591160 joint_SC0=184481683938
52/52 literal Fraction tree drifts passed: {'apex_cross': 2, 'birth': 225, 'child_cross': 16, 'large_promoter': 82, 'nonroot_apex_cross': 1, 'nonzero_cross': 18, 'shape_10': 0, 'shape_11': 4, 'shape_12': 0, 'shape_13': 3, 'shape_14': 3, 'shape_15': 0, 'shape_16': 5, 'shape_8': 1, 'shape_9': 2}
exact split and paired fiber check: {'apex_inner': 4, 'apex_tail': 6, 'child_inner': 38, 'child_tail': 57, 'transfers': 366}
b3f20_l5_hi.bin: mode=1 P/Q=5657570/10000000 counts=(6432404, 273529, 27746004, 13880592, 354856788, 540925385, 331822656, 131094480, 14215993) aux=21687 bytes=11370156264 logged_states=1421095995
4Sc0=88296060680 SC0=200835990548 SPhi(T3)=146044250013 joint_4Sc0=88296060680 joint_SC0=200701592409
52/52 literal Fraction tree drifts passed: {'apex_cross': 2, 'birth': 225, 'child_cross': 16, 'large_promoter': 82, 'nonroot_apex_cross': 1, 'nonzero_cross': 18, 'shape_10': 0, 'shape_11': 4, 'shape_12': 0, 'shape_13': 3, 'shape_14': 3, 'shape_15': 0, 'shape_16': 5, 'shape_8': 1, 'shape_9': 2}
exact split and paired fiber check: {'apex_inner': 4, 'apex_tail': 6, 'child_inner': 38, 'child_tail': 57, 'transfers': 366}
n=3: histories=1 E[C_n]=1 exact products bracket
n=4: histories=3 E[C_n]=4/3 exact products bracket
n=5: histories=12 E[C_n]=19/12 exact products bracket
n=6: histories=60 E[C_n]=647/360 exact products bracket
n=7: histories=360 E[C_n]=17189/8640 exact products bracket
n=8: histories=2520 E[C_n]=4717423/2177280 exact products bracket
n=9: histories=20160 E[C_n]=101568829/43545600 exact products bracket
```

## Small UBSan controls

`python3 -B round2/checks/b3f20_l5_sol_check.py --controls` built fresh `D0=10, Dc=6, B=2, mask=3` fixtures with the builder and ran each verifier with `--threads 1`, under `-fsanitize=undefined -fno-sanitize-recover=all`. The compiles and builder runs exited 0. The earlier 40-case literal suite passed on each fixture; the expanded 52-case production results are reported above. The valid runs with `--tree-check 30` each passed all nine tables and 34/34 GMP boundary trees, with `total_fails = 0` and exit 0: Lower at `5656316/10000000`, Upper at `5658254/10000000`.

With `--tree-check 0`, the negative controls reported:

| Control | Observed result |
|---|---|
| Lower `--override-pq 5656336 10000000` | `VERDICT: FAIL`, `total_fails = 2042023`, exit 1 |
| Upper `--override-pq 5658234 10000000` | `VERDICT: FAIL`, `total_fails = 2897705`, exit 1 |
| Upper copy, `KS3=4` | `invalid header (err=-21)`, exit 1 |
| Upper copy, last `q5=1000000000000001` | `q5[2127087] = 1000000000000001 exceeds verified native bound 1000000000000000`, exit 1 |
| Upper copy, last `q5=-1000000000000000` | One L5 failure and nonpositive upper `c0`; `total_fails = 2`, exit 1 |
| Upper, verifier built with `-DL5_PLANT_BUCKET` | 137,020 L5 failures; `total_fails = 137020`, exit 1 |
| Literal checker with a zero crossing cap | Rejected at `BIG_INNER crossing corner` |

No UBSan diagnostic was emitted. All payload edits were confined to the temporary fixture copy.

## Fresh production replay

The verifier was freshly compiled from the audited checkout with GCC 13.3.0:

```sh
gcc -std=c11 -O3 -fopenmp -Wall -Wextra -Werror round2/cert_l5/cert_b3f_l5_verify_gmp.c -lgmp -o [local path]
```

Exit 0, no diagnostics; `ldd` showed `libgmp.so.10`, `libgomp.so.1`, and libc. The pre-replay hash command exited successfully and both hashes equal the pins:

```text
$ sha256sum b3f20_l5_lo.bin b3f20_l5_hi.bin
9936d424fb3d87bc9193c79093870f1dd23772c16dd8b27a66233bc42a5a8474  b3f20_l5_lo.bin
cd4553129c758b2c8f74e3f1686146a613e8b4f01f4f720a6a5b86da09b8d317  b3f20_l5_hi.bin
```

Replay started at 16:59:57 UTC on November 24, 2026, after the lead released the four-thread slot. No competing verifier was running; `free -h` showed 29 GiB available. The two commands are chained sequentially under `OMP_THREAD_LIMIT=4`, `OMP_NUM_THREADS=4`, and `ulimit -v 12582912`:

```sh
/usr/bin/time -f 'wall=%e rss_kib=%M exit=%x' [local path] b3f20_l5_lo.bin --threads 4 --tree-check 30
/usr/bin/time -f 'wall=%e rss_kib=%M exit=%x' [local path] b3f20_l5_hi.bin --threads 4 --tree-check 30
```

Lower replay passed: all nine table checks reported zero failures, the comparison scan returned `(76840591160,184481683938,146042051850)` for `(4*S*c0,S*C0,S*Phi(T_3))`, and all 34 exact-rational GMP trees passed. Final lines:

```text
=== VERDICT: PASS (total_fails = 0, exact P/Q = 5656929/10000000) ===
wall=797.37 rss_kib=11138432 exit=0
```

Upper replay started after Lower exited, with 32 GiB available. It passed the same nine table checks and all 34 GMP trees. The comparison scan returned `(88296060680,200835990548,146044250013)`. Final lines:

```text
=== VERDICT: PASS (total_fails = 0, exact P/Q = 5657570/10000000) ===
wall=699.87 rss_kib=11139200 exit=0
```

The printed state counts for both runs were `f=6432404`, `g=121693`, `hV=27746004`, `hC=13880592`, `P4=354856788`, `Yt=540925385`, `Yb=331822656`, `S3=131094480`, `L5=14215993`. After both runs, `sha256sum` on the same two input paths and `diff -u` against the pre-replay hash log returned exit 0, with no differences, at 17:25:25 UTC.
