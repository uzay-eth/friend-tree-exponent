NEON-2 CLAIM 6: HOLDS

The stated gold1 bracket survives this attack, taking existence/identification of the mean exponent as the cited local result. `PROOF.md` omits passage from the edge seed, overstates what vanishes at a far-tail vertex if read literally, and has misleading attribution/dependency wording. The repairs below do not change the endpoint arithmetic or the unconditional bound. This is a GPT-family (`GPT-6.2 Sol Ultrafast`) adversarial review. The board credits Gemini 4 Neon with the proof, builder, 128-bit verifier, tail and combiner, and GPT-6.2 Sol Ultrafast with four prior hostile audits, the Python oracle and UBSan/SHA/binding probes; Gemini independently replayed those checks. My evaluator is new code, independent of that oracle and the verifier's generator/inequalities, but uses the published C++ layout/indexing API.

## Proof audit, in line order

| `PROOF.md` lines | Verdict | Adversarial check |
|---|---|---|
| 1–7 | Follows, subject to the details below | The verified claim concerns the two fixed cert files and mean growth. No orientation is used in the potential; rooting appears only in a counting bound. |
| 9–12 | Follows with an argument | The unique path between two nonleaves has nonleaf intermediate vertices. Thus the induced graph is connected, acyclic and has `C-1` edges, and its only possible zero-degree vertex occurs for `C=1`. |
| 14–22 | Follows | The target `u` is chosen with rate `Σ_{v~u}1/D_v`. At a core this adds 1 to `D` and leaves induced adjacency and `c` alone. At a leaf the only donor is its parent `u`; that leaf becomes `J`, and old `u` changes `c→c+1` but not `D`. The parent has `L` exchangeable leaf targets. When `C=1` the center has gain rate `D` and total promotion rate 1. The edge seed is explicitly outside this formula. |
| 28–43 | Follows | A physical endpoint has `1≤c≤D`. If `D<d`, exact or overflow is disjoint and exhaustive; when `D≥d`, both interpolation weights are nonnegative and sum to 1, including weight 1 at the seam. Taking the mean of `inv` gives `1/D`. |
| 45 | Follows | A fine exact or fine overflow label has exact `D`; monotone cutoffs make its clamped representative sufficient in either coarse case. Fine-tail mass goes from `d_j/D` to `d_j/D × d_{j+1}/d_j`, and fine `INF` maps to coarse `INF`. Core-degree clamps compose. |
| 47–61 | Follows | Singleton, edge, 3-path, 4-path, 3-star, 5-path, broom and 4-star exhaust connected tree types through order 5. Path reversal, star permutations, and the broom's two short arms are quotiented, never counted as separate embeddings. Header inequalities also imply pair-sum support inclusion: `SP5` rows impose both `SP4` bounds, and `SB5≤SP4`. |
| 62 | Follows with an argument | The internal **minimum** on `c` cannot cause an old embedding to enter support: its actual internal incidence in that embedding already forces the minimum. The upper bounds on `c,D` and sums of `D` are monotone under events. Deleting a promoted leaf supplies a supported old parent even though the attachment vertex has one less **physical** `c`; if still internal its remaining motif incidence ensures the minimum; if now an endpoint there is full support. The complete newborn/display accounting is below. |
| 68 | Follows | Signed units are ticks. `add_rat` rounds each signed rational upperwards at `1/1024` tick and the final cap rounds upwards to integer ticks. Thus caps are one-sided enclosures, not claims that all nominal signed residuals are nonpositive. |
| 70–73 | Follows | A term `(inv_v−r1(u))Δ^L f_u` is removed once from `f_u`, added once to its incident edge; `g→H3` likewise once for each old edge and external donor triple, both possible donors represented in the unordered `H3` row. |
| 74–76 | Follows | For fixed `u` the input and output of the transfer have equal multiplicity `c−1` by swapping ordered pairs `(v,w)`. The proof needs piecewise affinity (verified knots from 0 through 1/2) and convexity of a maximum of affine functions, not convexity of the stored primal `A`. Incoming interpolation for a physical tail is exact on the first segment. |
| 77–78 | Follows with an argument | For *each* source embedding and recipient, the eligibility test is constant among all its external donors: only the recipient/internal states determine support, as the extra donor is a new endpoint. Any donor then extends the unique connected support graph with no cycle. Charge the correction once to that extension. Reverse by deleting the donor; where several deletions are possible each corresponds to a distinct original term. All other external sums lie in the indicated interval. The eligibility is a **pre-event** test, unlike the newborn test. |
| 80–88 | Follows | The verifier recomputes fine differences in the old motif, including zero for a post-gain internal state outside support. Convex averaging of the absolute-error inequality gives a bound on the physical source difference conditional on the projected coordinates. The donor lies outside the source and is an independent interpolation coordinate, so multiply by its signed coefficient and use `|coefficient| × radius` before averaging. |
| 90–91 | Follows | At an exact overflow endpoint, only its true `D` matters to the table; promotion keeps its clamped class and gain maps to the next exact overflow or the first tail corner. It cannot be an internal vertex of a supported newborn extension at the same or narrower cap. With one incident motif neighbor removed, the rest is at most `D-c+(c-1)/2≤D-(K+2)/2`. The removed donor is included separately. |
| 92 | Follows with an argument; literal overflow `b=0` is WRONG | `t(D+1)−t(D)=−t/(D+1)`; for exact `c=κ`, `(D−κ+iota_nbr,lo)/(D+1)` is increasing in `D`, so the proof's lower endpoint is valid. For overflow `0≤γ≤1` suffices. The promotion interval comes from monotonicity of `1−κ/D`. Overflow promotion has **zero weight change for that old endpoint**, even if its physical rate `b` is positive (e.g. `D=180,c=27`: `b=153/180`). Read line 92's parenthetical “0 for overflow” as applying to the *contribution*; this is exactly how C++ implements it. A physical tail donating to an inner neighbor creates `t/d × ((1−t)d0+t d1)`; assign all of it to the mass-`t` corner with the max of the signed slopes. Conditional on every other coordinate this remains a valid bound, including multiple tail donors to the same inner hub. |
| 93 | Follows | Two adjacent tail endpoints are possible only on `g`. Both reciprocal donor credits are already inside the respective actual gain rates in the tail formula. `INF,INF` has zero value by the structural test, neither endpoint is exact (so no `f→g` or FY terms), and all tail transitions are charged in their own non-`INF` mass. `INF` on higher rows need not have zero weight or zero transfer. |
| 95–96 | Follows with an argument | The neighboring old core's degree never changes on a promotion. For `D≥DF`, the promoter cannot be internal of a supported newborn motif, so only the newborn singleton and edge are charged to `f`. Existing motifs in which the hub is an outer endpoint can still change; their changes stay in their own rows. Gain gives `(q',x')=((q+x)/(1+x),x/(1+x))`, promotion gives `(q−x,x)`; gain rate `ρ/x`, promotion rate `q`, newborn edge affine in `x` conditional on the class. Taking the signed extremum over the class and `ρ∈[q,(1+q)/2]`, then multiplying by `(1+x)^2`, gives the claimed polynomial, extending to `x=0` after cancellation. Bernstein basis functions are nonnegative and sum to one on each mapped rectangular panel, so the largest tensor coefficient bounds the polynomial. The sign split for division is correct even if the polynomial changes sign within a positive-upper-bound panel. C++ evaluates `D=1999→2000` with the same exact formula. |
| 102–105 | Follows with an argument | Bound positive higher row caps by the host counts below; use `C−1≤C` for edges only after replacing their cap by `max(cap_g,0)`; allow signed `e_f` and then replace the total by its nonnegative part. The proof does not omit `C=1` or `C=2`. |
| 107–111 | Follows with an argument; line 111 is WRONG if read as vanishing of incident motifs | Assign one edge per nonroot; allow zero at the root. Interior host multiplicities are below. The handshaking identity is exact including `C=1`; as `α≥0`, `Σ shifted W≥m_min C` implies the unshifted bound since `SΦ=Σ shifted W+2αS`. Comparability in the singleton case is checked **unshifted** separately, rather than requiring `f(D,0)−2αS≥m_min`. The unshifted `W` supplies the upper bound. `D≥DF` means no *hosted* higher motif; an incident higher motif may have nonzero weight with that vertex as endpoint. The tail `g` bounds use both endpoints of `t∈[0,d2/DF]`, across all parent corners plus zero. The far-tail lower checks separately use `c=0,1,2`, and for `c≥3` use `α(c−2)≥α`. The lower shift is `S/4` ticks, upper shift 0. |
| 113–116 | Follows | Positivity and `C≤SΦ/m_min` absorb `E_tot C`; the lower sign subtracts the penalty and the upper adds it. Constants belong to their own side. |
| 118–126 | Follows | Yule mean and second moment finite, monotonic size bounds every integrand by one fixed integrable random variable on `[0,t]`, and the stopped chain has a finite state set (including the hit of size `M`). Limit by dominated convergence and Fubini gives absolute continuity; comparison of `e^{−λt}EΦ` by its a.e. derivative proves the two Grönwall inequalities. |
| 128–133 | Follows for the 3-vertex seed; GAP in text for original edge seed, filled below | Nonzero prefactors from comparability disappear after `log/t`; no positive seed potential is claimed at size 2. The first-jump argument below is required for the statement about the original edge-seeded process. |
| 134 | Follows using `proofs/existence.md`; WRONG if crediting Addario-Berry et al. with existence | The local cited file proves existence (mean exponents, §§7–8) and says it is Astra's result written up by Claude Opus 5.5. The repo README explicitly says the 2024 paper left existence open. Thus the board's “no explicit Astra dependencies” should be restricted to the unconditional bounds, or the bracket should cite the local existence result as a separate input. |

### Generator: distance and multiplicity

Let `w(M)` be a motif's full rational weight, extended by zero off support. Every old embedding contributes `Σ_{u∈M}(r_u[w(M with D_u+1)−w(M)]+b_u[w(M with c_u+1)−w(M)])`. Recompute the weight at each member, even a physical tail; if a support bound fails, the difference includes the whole loss. A leaf gain outside the motif changes no coordinates in it. A promotion adds exactly one new core; only the old parent changes its decoration, and new motifs all contain the new core as a motif leaf. The graph distance from the attachment target to another member can be four for a newborn `P5`, or five from the physical target to the far vertex of an **old** diameter-four motif containing the parent. This is already covered by updating the old parent's decoration in every old embedding containing it; the ringing source can also be external and is counted in `r_u`.

| Source old motif | Newborn extensions charged to it | Recipient/donor corrections into the indicated larger family |
|---|---|---|
| `f` | `g` at its sole vertex; the separate newborn `f(J)` also charged here | `g`: each endpoint may be the recipient from the other |
| `g` | `H3` at either end | `H3`: center may gain from either end, deleting that donor leaves the old edge |
| `H3` | `P4` at either end; `S3` at center | `P4`: delete either outer donor; `S3`: delete any one of three arms |
| `P4` | `P5` at either end; `B5` at either internal vertex | `P5`: delete either outer donor; `B5`: delete either of two short donors |
| `S3` | `B5` at any of three outer vertices; `S4` at hub | `B5`: delete the long-arm tip; `S4`: delete any one of four arms |
| `P5`, `B5`, `S4` | None (order capped at five) | No outgoing larger-family correction |

Comparability and residual cap host counts for a vertex of core degree `c`: `H3` has `binom(c,2)` choices at its center; `S3,S4` have `binom(c,3),binom(c,4)` at their hubs; a half share of `P4` at internal `u` has `1/2 × Σ_{v~u}(c−1)(c_v−1)` choices; `P5` at center has `Σ_{unordered {v,w}⊂N(u)}(c_v−1)(c_w−1)`; `B5` at hub has `binom(c−1,2) Σ_{v~u}(c_v−1)`. These count *unoriented* vertex sets: no pair of branches can meet in a tree. The verifier combines counts times extrema over allowable internal labels, retaining factors from neighbor degrees in the maxima and including zero. Supports with high-degree endpoints do not force high degree at the host. `verify_quintet.cpp` lines 1220–1355 implement these bounds, with conservative integer rounding.

### Edge seed and identification

Write `F_3(t)` for the mean from the three-vertex tree and `F_2(t)` from one edge. The first ring time is `Exp(2)`, the process then starts from the unique unlabelled three-vertex tree, and the core count is zero until that jump. Therefore `F_2(t)=∫_0^t 2 e^(−2s) F_3(t−s) ds`. With `A e^(ℓt)≤F_3(t)≤B e^(ht)` and the positive gold1 endpoints, the upper bound is `B e^(ht) × 2/(2+h)`; restrict the lower integral to `[0,1]` for `t≥1`, giving a positive constant times `e^(ℓt)`. This completes the stated unconditional edge-seed result. A separate check on discrete size needs no de-Poissonisation: `LΦ=nE[ΔΦ|T_n]`, so iteration from `n=3` gives bounds by products `∏_{j=3}^{n−1}(1+λ*/j)`, whose logarithms are `λ* log n+O(1)`. Thus even without cited existence the same bounds hold for the discrete lower and upper mean exponents; existence is used just to name their common value.

## Replay and independent evaluation

All command paths below are relative to `round2/children/n2c6` unless absolute. CPU bound by `taskset -c 12 nice -n 5`; long component commands used `timeout 590`. Archive inputs remained read-only. Commands (shell variables used for readability):

```sh
git -C NEON2_REPO/ archive b04cc50 round2/neon2/quintet-cert | tar -x -C clean
q=clean/round2/neon2/quintet-cert
b=NEON2_REPO/round2/outputs/quintet-gold1
sha256sum "$b"/quintet_{lower,upper}.bin
g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -Werror -pedantic "$q/verify_quintet.cpp" -o verify_quintet
python3 -B "$q/cert_check.py" --verify-bin ./verify_quintet --lower "$b/quintet_lower.bin" --upper "$b/quintet_upper.bin" --threads 2 > replay.json
cmp replay.json "$q/certificate_gold1.json"
python3 -B -u "$q/audit_independent.py" "$q/certificate_gold1.json" "$b/quintet_lower.bin" "$b/quintet_upper.bin" --negative-controls "$PWD/verify_quintet"
g++ -std=c++17 -O3 -Wall -Wextra -Werror -pedantic eval.cpp -lgmpxx -lgmp -o eval
for side in lower upper; do for group in validation snapshots adversarial; do python3 -B -u run_eval.py "$side" "$group"; done; done
```

Archive/build completed without output. Hash outputs match the board:

```
81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986  quintet_lower.bin
8a36c3afdd2a8f4bbc03461a306af04594c46ada6ef037646261da4d364d9383  quintet_upper.bin
```

`eval.cpp` evaluates the raw motif weights as GMP exact rationals and never reads centering, FY, Bernstein, cap or envelope code. For small trees it explicitly appends the new leaf at each physical vertex and completely recomputes the potential (`brute`); this validates each event delta, not only its sum. For larger trees it computes deltas by visiting every supported old embedding once and dispatching its old and newborn changes to the affected vertices; the final drift sums `Σ_{v~u}1/D_v` times that exact delta over **every** physical attachment vertex, including all expanded leaves. Symmetric leaf attachments can reuse the same delta. `run_eval.py` generates the claimant's 47 and 142 trees plus the lead snapshots and new stress trees; exact comparisons use `Fraction`. In particular the latter do not merely test the local row formulas. Full re-evaluation for every vertex of the largest snapshots would be quadratic, so the checked local dispatch acceleration is used there. Every tree is evaluated with `m_min C≤SΦ≤m_max C`, signed nominal residual cap, and the adjusted exponent bounds. The endpoint ratios shown below are rounded for display only; margins retain exact fractions in `exact.log` for the tightest tree identifiers/ratios.

Clean `--threads 2` replay: `cmp` exited 0, no differences; SHA-256 of `replay.json` was `d26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb`. Exact endpoint fractions were `2825181729593748861/4994114995300000000` and `439124781075521909/776183857175000000`. Each side had 191,375,747 finite rows and 29,879,919 envelope inequalities. Lower `(err,m_min,m_max)=(19186,99882299906,1808571640916)`; upper `(19783,155236771435,1808494729937)`. The original 47-tree oracle's worst signed nominal residual/core was `−801680821011/500000000` lower and `−68525755853/25000000` upper. Supplemental standalone `python3 -B "$q/verify_tail_bernstein.py" "$b/quintet_$side.bin"` for each side confirmed `DF=2000`, 16 panels, exact tail caps `0` and `−868735957253/133466700` (integer upper cap `−6509`) respectively.

`python3 -B -u crosscheck.py lower` and `... upper` compared the already-computed independent results against the claimant's oracle, separately from evaluation. Both printed `47 per-tree exact Phi and drift matches`. They also gave a concrete counterexample to a literal reading of “all higher motifs vanish” at line 111: the 3-core path with degrees `(2000,2,2)` has nonzero `H3` ticks `1649232007178/125` lower and `1660320950412/125` upper. Its middle (not the degree-2000 endpoint) hosts the term, so the actual comparability check includes it.

The `group` variant rebuilt via `g++ -std=c++17 -O3 -Wall -Wextra -Werror -pedantic eval.cpp -lgmpxx -lgmp -o eval_group` also forced complete recomputation for one representative of every leaf orbit at a core vertex, plus each individual core target, on a 2,000-vertex seam star, a 100,001-vertex tail star, and a 2,023-vertex tail broom, both sides. All per-target equalities passed. On the seam star it printed `(SΦ,SℒΦ)=(1808362214356,4089933184598956597/3998000)` (lower) and `(1808494615915,2045282269880375/1999)` (upper).

### Exact tree results

Each row reports min/max drift divided by potential, and minimum margin against the **adjusted exact** endpoint, lower `ratio−λ*` or upper `λ*−ratio`. Certified comparability intervals for `Φ/C` are approximately `[0.090842422565,1.644886325190]` and `[0.141187021140,1.644816375062]`. All numbers below are display-only rounding.

| Side / input | Trees | Expanded vertices | Drift/Phi range | Min margin | Phi/C range |
|---|---:|---:|---|---:|---|
| lower claimant concrete | 47 | 3–2,201 | 0.565702371696–0.566046888729 | 1.93782026260e-7 | 0.683949982953–1.644697508494 |
| lower claimant audit suite | 142 | 3–1,002,033 | 0.565702371961–0.632494970650 | 1.94046857271e-7 | 0.738435995681–1.644398272687 |
| lower `nsnap2_T3.bin` | 800 | 7–217 | 0.565702371941–0.566001251205 | 1.94027301974e-7 | 0.778548374928–1.643059709417 |
| lower `nsnap2_T5.bin` | 600 | 43–1,563 | 0.565702371866–0.565890233908 | 1.93951917162e-7 | 0.874997051270–1.643831281486 |
| lower `nsnap2_T7.bin` | 300 | 379–10,576 | 0.565702371928–0.565775751865 | 1.94013985811e-7 | 0.904220845433–1.644616607488 |
| lower `nsnap2_T8.bin` | 150 | 1,153–32,164 | 0.565704893324–0.565778980832 | 2.71540968341e-6 | 1.019373982217–1.270514968821 |
| lower `nsnap2_T9.bin` | 60 | 2,434–91,508 | 0.565703757048–0.565775593258 | 1.57913409735e-6 | 1.035921813732–1.193462370417 |
| lower `nsnap2_T10.bin` | 20 | 13,879–131,832 | 0.565708539579–0.565732955092 | 6.36166546214e-6 | 1.098518650030–1.199419874255 |
| lower new adversaries | 1,780 | 3–2,200,001 | 0.565702371696–0.736675936811 | 1.93782026260e-7 | 0.547681838292–1.644882516757 |
| upper claimant concrete | 47 | 3–2,201 | 0.563698336019–0.565748276385 | 1.31052417034e-7 | 0.684122058056–1.644816375062 |
| upper claimant audit suite | 142 | 3–1,002,033 | 0.539261331680–0.565748276379 | 1.31058844339e-7 | 0.738533052345–1.644519364681 |
| upper `nsnap2_T3.bin` | 800 | 7–217 | 0.565415535907–0.565748276384 | 1.31053820791e-7 | 0.778965883399–1.643048985115 |
| upper `nsnap2_T5.bin` | 600 | 43–1,563 | 0.565413233182–0.565748276413 | 1.31024992441e-7 | 0.875140017107–1.643829258026 |
| upper `nsnap2_T7.bin` | 300 | 379–10,576 | 0.565652350334–0.565748276345 | 1.31092785850e-7 | 0.904332358199–1.644688512526 |
| upper `nsnap2_T8.bin` | 150 | 1,153–32,164 | 0.565641000089–0.565743617757 | 4.78968062171e-6 | 1.019526099879–1.270617591384 |
| upper `nsnap2_T9.bin` | 60 | 2,434–91,508 | 0.565666030541–0.565746483559 | 1.92387908264e-6 | 1.036095352635–1.193564776472 |
| upper `nsnap2_T10.bin` | 20 | 13,879–131,832 | 0.565709481195–0.565740291036 | 8.11640177128e-6 | 1.098700659238–1.199552633537 |
| upper new adversaries | 1,780 | 3–2,200,001 | 0.485439394003–0.565748276507 | 1.30931005803e-7 | 0.548077310184–1.644816375062 |

Totals per side: 1,930 lead trees, 1,780 new adversarial trees, and 189 claimant validation trees (3,899 evaluations; all four violation counts were zero). Of these, 1,709 per side also had literal full recomputation at each attachment vertex. Adversarial trees include bare paths, decorated caterpillars, brooms, double stars, cherry hubs, spiders (legs 2–4), single stars up to degree 100,000, and 1,500 reproducible random trees biased to cutoff seams. Overall nearest witnesses were the lower degree-1999 star (`4089933184598956597/7229832132995288000`) and upper degree-1000 star (`5115218397283/9041509465780`). The latter also passed the full representative recomputation, printing `(SΦ,SℒΦ)=(1808301893156,5115218397283/5)`. Closest among the lead snapshots only: `T5` index 62 (zero-based) lower and `T5` index 126 upper.

## Published standalone audit and negative controls

The first published standalone-audit attempt used `--negative-controls ./verify_quintet`; `argparse` parsed it as a `Path` and `str(Path('./verify_quintet'))` became `verify_quintet`, which was not on `PATH`. That attempt stopped with `FileNotFoundError: [Errno 2] No such file or directory: 'verify_quintet'` after the lower side's 142 trees. The absolute-path command above is the board-equivalent replay.

The absolute-path replay completed (full stdout in `audit_absolute.log`). Both sides: 114 projection checks; sampled rows in `f,g,H3,P4,S3,P5,B5,S4` order `73,61,549,549,488,549,549,366`; sampled envelope checks in `e23,e34p,e34s,e45p,e45bp,e45bs,e45s4` order `495,495,440,495,495,495,440`; 54 tail checks; 142 physical trees. Worst signed residual/core `−801680821011/500000000` lower and `−369804355371/125000000` upper, matching the independent rerun. On each side `ginf`, `env23`, `env45s4` were rejected with the corresponding diagnostics. Its tighter nominal P/Q control returned finite-only penalties `231371405` lower and `248690151` upper, and best-possible adjusted endpoints `2814113530143278861/4994114995300000000` (lower) and `440290514529804409/776183857175000000` (upper).

Physical file mutation: `python3 -B -u control.py corrupt`, created `mutated_lower.bin` by copying the exact lower file and setting `g(INF,INF)` from 0 to 1; the unmodified verifier returned code 2 with `ERROR: g[INF,INF] != 0`. A subsequent `cmp -l ORIGINAL mutated_lower.bin` found exactly one byte difference (`98043321 0 1`, 1-based offset, octal bytes). This checks the actual rejection, not merely a digest change.

Tightening uses the real `--override-pq P Q` flag plus the shipped analytic combiner, via `control.py`, with a requested target strictly on the forbidden side of the witness ratio and within `10^-12` of it. The C++ verifier is a **cap calculator** for valid tables; it exits 0 at an overly strong nominal value and reports caps. Rejection of that exponent is correctly defined by the **adjusted** endpoint moving to the allowed side of the witness, rather than expecting exit 2 for an otherwise valid table. In both witness directions the script asserts that property.

```sh
python3 -B -u control.py lower 4089933184598956597/7229832132995288000
python3 -B -u control.py upper 5115218397283/9041509465780
python3 -B -u control.py upper 1022061000577/1806564939195  # additional lead-snapshot witness
```

- Lower requested `35356398231/62500000000`, exact violation `3346257754441/56483063539025687500000000`; computed `err_per_core=20114`, adjusted to `1765738557289854933143/3121321872062500000000` (`0.565702170318`). Rejected at requested exponent.
- Upper requested `282874138253/500000000000`, exact violation `132005375883/226037736644500000000000`; computed `err_per_core=21293`, adjusted to `8782495718270710240611/15523677143500000000000` (`0.565748413671`). Rejected at requested exponent.
- Additional upper lead-snapshot request `141437069103/250000000000`, exact violation `53958361583/90328246959750000000000`; computed `err_per_core=21466`, adjusted to `4391247867055742294561/7761838571750000000000` (`0.565748414692`). Rejected at requested exponent.

## Scope

This supports the **fixed gold1 inputs** and the outward bracket `[0.565702177913,0.565748407438]` (width at most `46.2296e-6`), with `μ>0.565700`. I did not independently re-prove existence or convergence to the normalized limit `X`, rebuild the binary tables from the optimization, reimplement exhaustive checking of the 191M row inequalities and 29.9M envelope inequalities, or replay prior UBSan and pathname substitution suites. I audited row coverage and the match of the implemented arithmetic to the inequalities, then replayed the exhaustive verifier and independently tested the physical potential and drift. The new exact evaluations do not prove universality by sampling; the local-row/analytic-tail argument is what does so. The claim does not certify eight digits.
