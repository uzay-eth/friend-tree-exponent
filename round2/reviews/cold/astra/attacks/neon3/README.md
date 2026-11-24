# Astra attack on original Neon Claim 3 — B3-FullH-20-L4

**Verdict: HOLDS.** Scope is original Claim 3 on `neon.md`, exact Neon pin `51f4cfadac675fdf1cef84d44935eebb5456aabb`, and the two D20 payloads identified below. This is a review of the hash-pinned production claim, including the finite-size formulas, with separate checks of the generic repair paths.

## Fresh isolated build and replay

Snapshot and live receipts: `[local path]`, `production.log` and `extras.log`. Created by `git archive` of `round2/cert`, `round2/checks` and `round2/proof-b3f-motif4.md` from the exact commit. Production began at 13:19:34 UTC, on CPUs 2 and 6, with two GMP workers total and no overlap between witnesses. Hashing used CPU 2. Both production verification and hashing used nice +15 and idle-class I/O. The builder was not used to regenerate either production file. The process reads one payload into anonymous memory at a time (10.42 GiB of table data, plus small working state), and no payload copy is made. Depending on shared load, some mapped-in table data can be swapped; live RSS is not the same as total allocated table bytes. A separate post-production queue runs single-worker checks only after both production verdicts have passed.

| Pin | SHA-256 |
|---|---|
| `cert_b3f_verify_gmp.c` | `8b97ecd87176c8709f25c4caac9d300792df4f09a260983f11134a3dbfa0f88e` |
| `cert_b3f_format.h` | `18a9db19f8d15ca07cd3505fd0b52246bc8f18ac0fb864fb37ea2a8bc364190c` |
| `cert_b3f_build.c` | `325ab2bb6aae51a55a9f03a893dafe22a038d2cbaf9f3c77fc1b62e9881516e6` |
| `proof-b3f-motif4.md` | `97197b8e146df91bb014ff89659df44d22e4e59b094965a60101eae21d857f06` |
| Fresh verifier executable (GCC 13.3.0, x86-64 LP64; GMP, GOMP, libc) | `d7c7b7681024cdd5000f0fcd361bb83d0dd4cd576eca7c0fb164c60d8b454ff5` |
| `b3f20_lo.bin` (11,190,881,056 bytes) | `eb93985b2c90a516b65e85934f4ac847779f42b3fdbfcd2af4dc07fd30c37469` |
| `b3f20_hi.bin` (11,190,881,056 bytes) | `55dd6015eca85ca2fba6267d9444ec8807554c456999803d8f9abcc695c090e2` |

Reproduction from the Astra root, with an empty, pre-existing scratch directory `SNAP` and the two read-only external payloads at the paths above:

```sh
CPU1=2 CPU2=6 bash round2/attacks/neon3/replay.sh "$SNAP" > "$SNAP/production.log" 2>&1
CPU1=2 bash round2/attacks/neon3/extras.sh "$SNAP" > "$SNAP/extras.log" 2>&1
```

Core command in the replay: `gcc -O3 -march=native -fopenmp -Wall -Wextra -Werror round2/cert/cert_b3f_verify_gmp.c -lgmp -o verify-gmp`, then `nice -n 15 taskset -c 2,6 ionice -c 3 stdbuf -oL ./verify-gmp BLOB --threads 2 --tree-check 30` for each sign, serially and only after its length and full SHA-256 match.

### Production receipt

Both witnesses reported `VERDICT: PASS (total_fails = 0)` at their respective exact endpoints, with zero failures on 34 exact-rational concrete/boundary trees each. The archived production receipt (`production.log`) is byte-identical to the posted receipt, SHA-256 `1b8f37f67875a75bd4a13dbd40f272d2763280204a2fb6ca70b7cc5d576af514`. This discharges the original signed GMP input for `a874930`, `a2e05e0`, and `ee93320`; it contributes no new occupation payment or numerical credit to those releases.

| Sign | Hash wall (s) | Verify wall (s) | User / sys (s) | Peak RSS (KiB) | Filesystem inputs | Exit |
|---|---:|---:|---:|---:|---:|---:|
| Lower | 16.13 | 2253.17 | 2158.75 / 278.41 | 10921344 | 65093736 | 0 |
| Upper | 11.19 | 1750.52 | 2110.38 / 71.03 | 10931840 | 26396560 | 0 |

Per sign, the verifier reports row counts `6432404, 121693, 27746004, 13880592, 354856788, 540925385, 331822656, 122901075` in `f,g,hV,hC,P4,Yt,Yb,S3` order, all with zero failures. The `g` sweep counts admissible pairs, whereas its compact-profile storage is a dense `523²` matrix. For `(4Sc0, SC0, SΦ(T3))`, the verifier returned lower `(113472636536,184219674916,146041684821)` and upper `(115454594604,200854054323,146044513514)`.

### Independent receipts

- The independent binary reducer matched all eight array lengths and all six integers. In `(L,k,h)` order, the minimum with recharge attained `(0,1,1)` on both signs and the upper charge maximum attained `(15840,18,1)` lower and `(16000,19,1)` upper. As a control, removing the recharge from the lower witness's charge comparison gives a negative raw minimum `−287538520728` on the `4S` scale.
- Exact physical drift passed 84 fixtures per sign through degree 25,020, including all eight motifs and the inner, tail, arity and leaf-stock seams. Per sign, 57 fixtures also matched the uncompressed, directed-arrow calculation component by component. Both deliberately wrong rate/promotion calculations changed the drift. Minimum oriented residuals divided by `S` were approximately `3.762115296e−6` lower and `3.67544908128e−6` upper.
- The independent recursion enumerated every labelled attachment history from the unique rooted `T3` through `n=10`: `1,3,12,60,360,2520,20160,181440`, with every exact mean between the two products. The last two means are `101568829/43545600` and `200684675729/80621568000`.
- The UBSan builds passed all nine invalid-input controls, plus the direct Ybot staging count. With `D0` drawn from `{20,40}` and `k` from `{2,19}`, the C controls agreed with exact interpolation of the physical row in 168 comparisons; a planted transfer-sign reversal gave 84/84 mismatches.
- `extras.log` contains successful repair/tail controls followed by the first oracle attempt, which stopped on Astra’s missing cumulative mass update in the labelled-history recursion. After adding that update, the following rerun returned `oracle wall=143.07 rss=10973176 exit=0`, including both full comparison reductions and both physical suites:

```sh
env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 \
  /usr/bin/time -f 'oracle wall=%e rss=%M exit=%x' \
  nice -n 15 taskset -c 2 ionice -c 3 \
  round2/.venv/bin/python -EBu round2/attacks/neon3/audit.py \
  b3f20_lo.bin b3f20_hi.bin
```

## Mathematical attack surfaces

The isolated release preflight also reran the corrected two-payload oracle, with the same one-worker environment and no CPU pin. See the archived oracle receipt (`oracle.log`): both signs passed, `oracle wall=167.94 rss=10977984 exit=0`.

- **Physical generator.** An ordered donor-to-recipient arrow has unnormalized rate `1/D_donor`, and total rate is `n`, independently of the tree. The core recipient gain rate is `L+Σ_core-neighbor 1/D_neighbor`; aggregate budding rate `L/D` attaches to one of its interchangeable leaves, thus keeping its own `D` fixed and creating a core child `(1,0,1)`. Both first-step outcomes from one edge have `T3` equal to the same rooted path. Our oracle compares full directed arrows to grouped gain and budding, re-evaluates all eight potentials, and tests threshold seam and high-arity fixtures. It deliberately perturbs recipient-rate and budding-degree bookkeeping.
- **All-degree inequalities.** Bernstein degree-three weights are nonnegative and reproduce `1/D`. They are algebraic weights, not an independence premise. At fixed recipient, the `f→g` scalar is `c2/D²+c3/D³` in tail; at an edge the same external donor and same centered T2/T3 multiplier are put back once in the compatible triple. Of the four-vertex motifs, a P4 has two adjacent triple deletions; Ytop has two vertical and one pair deletion; Ybot admits one vertical and one pair deletion; S3 has three pair deletions. Pre-event domain gates are the same on both sides of each transfer. A promotion uses its post-promotion state, never an inner-limit degree one too large. At `D0−1→D0`, an exact gain maps to T3; at `k=KT→KT+1` restricted slots disappear, whereas an unrestricted tail slot moves into its overflow arity. At `k>KT` the true rate factor `(D−c+ρ)/(D+1)` is in `[0,1]` even though represented arity is clipped. At `L=LF`, the bud at that actual leaf stock is evaluated separately before the `L≥LF+1` range; `VSTAR` still pays the newborn terms. The remaining T3 control ratio is bounded because its denominator Bernstein controls are positive, and the extrema in the exterior rates may be taken at the affine corners. The Astra tail probe reconstructs these controls by interpolating the physical event expression, while the full production replay exhausts the claimed table domains.
- **Comparison and products.** Bottom ownership gives one incoming edge, one `hV`, one `P4` and one `Yb` group with `k_apex−1` possible siblings. Apex ownership gives pair/pair/triple multiplicities for `hC/Yt/S3`. The comparison can bound each optional charge with zero included and optimize each tail contribution over its mixture support. The term `(k−h)/4` sums to zero on every rooted core tree, and the clipped `KF+1` is a lower bound for the actual overflow shift. The lower drift and positive seed give a lower bound on expected potential by induction; use `Φ≤C0 C` at the end. The upper drift uses `c0 C≤Φ` with `c0>0`. Both give a factor `1+λ/j` per step starting at `j=3` (empty product at `n=3`). The seed is `f(2,0,0)`, with no higher motifs. The logarithm of the product is `λ log n+O(1)` because the quadratic remainder is summable. This does not require a prior bound on μ; use the existing mean-exponent existence and a.s. equality theorems for interpretation.

## Repair scope at `51f4cfa`

- Corollary 2 now distinguishes potential-table count, auxiliary count, total entries and byte length. The older review's `22,687` auxiliary count is a historical typo; the claim and corrected corollary say `21,687`.
- Pair staging allocates the full Cartesian product, not a fixed 16,384 stack buffer. The problematic generic Ybot `ksum=24` group is counted without allocating any huge generic table in `repair_probe.py`; the logical count is independent of the GMP replay. These builder arrays are outside the production verifier.
- The relevant C allocations now go through checked wrappers; failure aborts and cannot be a certificate pass. ABI guards enforce the required header size/long/size_t widths and little-endian storage. The override path validates the original header first and the replacement second. Builder parsing checks range and finiteness; the finalized auto exponent is validated before integer use, and failed output close is detected. Production uses none of the exponent-override paths. The scripts still have a four-thread default in Neon's code; Astra always passes its own explicit budget. Scope does not include exhaustive arbitrary generic profiles, successful execution after allocation failure, or correctness of the floating candidate generator independently of verification.
