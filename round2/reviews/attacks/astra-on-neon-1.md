VERDICT: HOLDS

**Scope (2026-11-24, updated 10:14 UTC):** Cold attack on Neon’s numbered Claim 1 as posted at 09:42 UTC, proof and verifier pinned to `ft2-neon@cec3855cc52ab4265ce3e36b8a7a38f1ae257f77`, and the two named B3-FullH-18 certificate hashes. The exact enclosure and finite-size product constants survived. **Both original D0=18 production certificates passed the full GMP verifier compiled from a fresh clean detached clone**, in addition to the independent comparison, physical drift and T3 hull checks below. This is not an audit of arbitrary certificate headers or of any subsequent hardening commit. Fable Claim 1’s separate audit still awaits its final freeze.

### Full clean-clone production replay (owner order B)

Ran from `round2/attacks/neon1`:

```sh
REPLAY_OUT="$PWD/outputs/production-final" bash replay-production.sh
taskset -c 3 nice -n 15 bash package-rehearsal.sh
```

The script used `git clone --no-local --no-checkout -- [local path] <scratch>/clone`, detached at the full pin, checked that there were no object alternates, and checked a clean tracked and untracked checkout before and after compiling. For both named blobs, the checked-out bytes matched both the pinned Git blob content and the separately extracted pinned copies. Build command: `gcc -O3 -march=native -fopenmp -Wall -Wextra <clone>/round2/cert/cert_b3f_verify_gmp.c -lgmp -o <scratch>/neon_verify` (GCC 13.3.0). SHA-256 receipts:

| Item | SHA-256 |
|---|---|
| Pinned verifier source `cert_b3f_verify_gmp.c` | `f2a1b91575a82e0039b40ac1eb247e86291c803219eb6176306393c70475796b` |
| Pinned header `cert_b3f_format.h` | `05055d3cb789c3ae83ca20321b7e7def5f76b30bf2960f7f60f03c3947cd3729` |
| Recompiled executable | `6c45a210a42d3fb5d4487215135c826531f1f883043b6947c3febba85aef9b8a` |

Each original 3,495,091,968-byte payload matched its posted SHA-256 before and after its run. The exact runtime command per sign was `taskset -c 3 nice -n 15 /usr/bin/time -v <scratch>/neon_verify b3f18_{lo,hi}.bin --threads 1 --tree-check 30` with `OMP_NUM_THREADS=1` and `OMP_THREAD_LIMIT=1`. The lower process exited before the upper began; there were approximately 23 GiB and 26 GiB available before the respective launches. Neither process swapped.

| Production payload | Full GMP result | Concrete/boundary trees | Wall time | Peak RSS |
|---|---|---:|---:|---:|
| `b3f18_lo.bin` | `PASS (total_fails = 0, exact P/Q = 565673/1000000)` | 34/34 | 14:46.59 | 3,416,064 KiB |
| `b3f18_hi.bin` | `PASS (total_fails = 0, exact P/Q = 565775/1000000)` | 34/34 | 9:20.55 | 3,416,192 KiB |

All eight sections returned zero failures on both signs. The independent oracle reran on the same production files after GMP. Omitting edge recharge in the independent exact comparison again made the upper `4S c0` negative (`−23978775976`). Reversing only the `Y_ref` input in the pinned GMP T3 transfer while retaining the physical interpolation oracle produced **84/84** exact mismatches; the first was `1382688/475` versus `1376019/475`. All **168** unmutated GMP T3 comparisons passed. The production session exited 0 at 10:13:31 UTC.

Full stdout/stderr, exit codes and timing are in `round2/attacks/neon1/outputs/production-final/{session,lo,hi,oracle,tail-oracle}.log`. The compact rehearsal package `round2/attacks/neon1/outputs/neon1-production-rehearsal.tar.gz` (44 KiB on disk; SHA-256 `62f97976ed42c9c8ab195ed5d7b02e08d38c35674b387081dcdae339dc84b1d0`) contains the receipts, pinned source/header, oracle and replay scripts. It uses the two original payloads as external inputs.

### Physical generator and drift certificate

Starting with the one-edge seed, either first attachment gives the path with one core; it has `(L,k,h)=(2,0,0)` and only `f` is present. For every later step, enumerate all physical directed donor→recipient arrows at rate `1/D_donor`; their total rate is `n`, so divide the continuous drift by `n` for the single-step conditional mean. A core recipient gets leaf-donor rate `L` plus the reciprocals of *core donor* degrees. A leaf recipient has only one possible donor, its core parent; aggregating these arrows gives `L/D_parent`. Budding converts an existing physical leaf, leaves its parent’s degree unchanged and creates a core child `(1,0,1)`. This also fixes the Bernstein parameter `D0/D` during budding; only the parent’s `k` category can change. Gains and buds are then replayed by rebuilding all eight motifs, versus direct directional-arrow expansion on the full physical tree.

The independent evaluator uses rational Bernstein weights and does not import Neon’s evaluator or header code. Both D0=18 files pass physical arrows vs local event aggregation and the signed production drift inequality on 83 trees each; all eight motifs were exercised, including exact `D0` and `KT` crossing cases. From the original physical transition law the exact mean at `n=8` is `4717423/2177280`, and all exact means at `n=3..8` satisfy the posted finite-size products. Planted recipient-degree weighting and a spurious parent degree increase on budding each change the computed drift and are caught. These finite trees diagnose topology/rate bookkeeping; they are not a substitute for the unbounded proof.

For the tail transfers, independently checked the convex weights sum to one and reproduce reciprocal degree. On the one-body transfer the T2 and T3 weights recombine to `c2/D²+c3/D³`; on an oriented edge the external sum excludes the included neighbor. Its T2 and T3 weights multiply the same two terms moved into the distinct 3-vertex motifs. Checked 192 exact signed transfer fixtures. The 3→4 transfer similarly assigns an extra eligible neighbor to its unique extended motif; in the `Ybot` case the child continuing the chain and its sibling are different roles. In the nonoverflow T3 calculation I reconstructed the residual from physical gain/bud/inner-neighbor drift, stripped the first two Bernstein terms, multiplied by `(D0+t)^3`, and solved directly for rational Bernstein coefficients. Positive denominator coefficients give a convex combination of coefficient ratios; extrema over credit endpoints are safe because the denominator is independent of credit and the numerator is affine in it. **168** signed D0=10/18 T3 cases match the pinned GMP routine exactly. I did not independently duplicate every motif row or every overflow inequality.

### Independent integer ceilings / comparison

The comparison uses exact signed integer extrema. Incoming edge, vertical path and 4-chain each charge their terminal vertex once (where eligible). A `Ybot` ending there shares the upper two vertices and has `k_apex−1` possible siblings; the other fork/pair and star charges are apex multiplicities `binom(k,2)` and `binom(k,3)`. Optional charges include zero; the child-tail and apex-tail extrema account for the convex mixed states, including the zero apex value at INF. The one-body caps preserve the terms with physical `L>LF`; in the `k>KF` class the lower correction can safely use `KF+1`. Integer scans included every stored motif value used in the extrema, checked the accumulator range for the independent scan, then swept all local states with the exact `4S` lower and `S` upper scales.

| Certificate | `4S c0` | `S C0` | `S Φ(T3)` | Exact product prefactor |
|---|---:|---:|---:|---:|
| lower | 113465364604 | 184124847315 | 146042017510 | `29208403502/36824969463` |
| upper | 115462837236 | 195820675182 | 146044806770 | `146044806770/28865709309` |

The fractions equal the advertised unreduced fractions. Strict outward six-place decimal bounds are `0.793168` below the lower prefactor and `5.059457` above the upper prefactor. Dropping the edge recharge makes the upper comparison fail: `4S c0 = −23978775976` without it, versus `+115462837236` with it. This also independently explains the early upper failure in `cert_b3f_18.log`; `cert_b3f_18_verified.log` is the post-recharge run. The recharge is legitimate because each edge contributes once to `Σk` and once to `Σh`. For the lower growth bound only `Φ≤C0 C` and positive initial `Φ` are needed; the upper growth bound uses `c0 C≤Φ` with `c0>0`. Taking expectations and iterating gives the factors `1+(P/Q)/j`; `Σ log(1+(P/Q)/j)=(P/Q)log n+O(1)` gives the exponent bracket without using a prior exponent interval. Both comparisons to the specified earlier brackets and the claimed width `0.000102` are correct.

The payload has `436,870,169` **potential** entries, `436,886,484` int64 entries including centering/transfer arrays, plus the 96-byte header, for `3,495,091,968` bytes. The proof’s `436,870,000` entry count should be read as rounded; the board has the exact potential entry count. Both production SHA-256 values match the board exactly. The archive command builds from the pinned source; the proof-level use is limited to these production parameter/byte pairs, not a claim that the old native parser/products accept arbitrary input safely.

### Initial independent checks (09:44 UTC)

From `round2/attacks/neon1`, with the source extracted read-only via `git archive cec3855 round2/cert` from Neon’s repository:

```sh
taskset -c 13 nice -n 10 env OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 /usr/bin/time -f 'elapsed=%es' round2/.venv/bin/python -B audit.py --hash b3f18_lo.bin b3f18_hi.bin
```

```text
Bernstein physical credit: 192 exact one-body/T2/T3 transfer fixtures
sha256 b3f18_lo.bin 9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec
physical arrows: 83 trees, all 8 motifs, both planted errors rejected; E[C_8]=4717423/2177280
sha256 b3f18_hi.bin 1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41
comparison without edge recharge: 4*S*c0=-23978775976 at (0, 17, 1); positive=False
physical arrows: 83 trees, all 8 motifs, both planted errors rejected; E[C_8]=4717423/2177280
elapsed=10.35s
```

Further executed on CPU13 at niceness 10, sequentially: `gcc -O3 -fopenmp round2/cert/cert_b3f_verify_gmp.c -lgmp -o neon_verify`; `./neon_verify smoke_hi.bin --threads 1 --tree-check 12` → `PASS (total_fails = 0, exact P/Q = 565859/1000000)`, 6.28 s; corresponding `smoke_lo.bin` → `PASS (total_fails = 0, exact P/Q = 565597/1000000)`, 6.76 s. Planted upper exponent on `smoke_hi.bin --threads 1 --tree-check 0 --override-pq 560000 1000000` → `FAIL (total_fails = 4269587, exact P/Q = 560000/1000000)`, expected exit 1 in 6.25 s. `gcc -O2 -fopenmp tail_probe.c -lgmp -o tail_probe` and `python3 -B tail_probe.py` → `GMP T3 hull: 168 exact comparisons against physical interpolation; all pass`.
