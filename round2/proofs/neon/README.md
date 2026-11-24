# Neon proof and replay index

**Source pin:** `c2498cc482bfef861644092d33f727cbbc93b280`. Gemini 4 Neon authored the certificate approach, with GPT-6.2 Sol implementation and audits; GPT-6.2 Sol drafted Claim 2 and Gemini 4 Neon audited it. The [release ledger](../../LEDGER.md) records claim-specific tier decisions and the exact scopes of GPT, Claude and Gemini reviews.

From the repository root, run `bash round2/verify/neon/run.sh` (Python 3.10+, standard library; about one minute). It checks finite arithmetic and identities, including Bernstein and endpoint algebra, without validating production witness tables. For a compile-only check of both builders and both exact checkers, pass `--compile` (GCC 11+, OpenMP and GMP development headers). Paths below are relative to this page unless specified otherwise. [SOURCES.tsv](SOURCES.tsv) records original source Git objects. Archived receipts show artifact basenames; runner and audit paths are portable. Hashes, results, and checker rules are unchanged by the pathname edits.

| Board claim | Proof and check | Release scope |
|---|---|---|
| 1: B3-FullH-18 + Level 4 | [Theorem 1 and Corollary 1](round2/proof-b3f-motif4.md), [GMP checker](round2/cert/cert_b3f_verify_gmp.c), [hardened receipt](round2/outputs/cert_b3f_18_hardened.log); [Astra](../../reviews/attacks/astra-on-neon-1.md) and [Fable](../../reviews/attacks/fable-on-neon-1.md) | Historical enclosure `[0.565673, 0.565775]` and all-size product bounds. Superseded enclosure; production witness pair not released. |
| 2: Hubs and fixed-label forest | [Proof §§1, 3](round2/open_questions/REPORT.md), [exact enumeration](round2/open_questions/exact_degrees.py); [Astra](../../reviews/attacks/astra-on-neon-2.md) and [Fable](../../reviews/attacks/fable-on-neon-2.md) | Hub dichotomy, strict probability bounds and possible freeze after promotion, exact mean and bias formulae and strict separation of founders from label 3, finite-degree independent set and fixed-cutoff fixed-label limit. The enumeration checks finite fractions, not the infinite-time steps. Fable's Claude review excludes the rejecting-chain proof, which its GPT child checked. |
| 3: B3-FullH-20 + Level 4 | [Theorem 1 and Corollary 2](round2/proof-b3f-motif4.md), [GMP checker](round2/cert/cert_b3f_verify_gmp.c), [lower](round2/outputs/cert_b3f_20_lo.log) and [upper](round2/outputs/cert_b3f_20_hi.log) receipts; [Astra](../../reviews/attacks/astra-on-neon-3.md) and [Fable](../../reviews/attacks/fable-on-neon-3.md) | Historical enclosure `[0.565680, 0.565768]` and all-size product bounds. Superseded enclosure; no released production witness pair. The GPT cold review's entry-count and generic-input fixes were applied in `51f4cfa`; Fable replayed earlier unchanged row logic and tested the corrected checker with negative controls. |
| 4: B3-FullH-20 + coarse endpoint-only Level 5 | [Theorem 2, Proposition 1 and Corollary 4](round2/cert_l5/proof-b3f-l5.md), [GMP checker](round2/cert_l5/cert_b3f_l5_verify_gmp.c), [independent raw and tree checker](round2/checks/b3f20_l5_sol_check.py), [lower](round2/outputs/cert_b3f_20_l5_lo.log) and [upper](round2/outputs/cert_b3f_20_l5_hi.log) receipts | Historical enclosure `[0.5656929, 0.5657570]` and all-size product bounds, superseded as an enclosure. GPT cold replay and proof audit; no rival-family attack before close. No released production witness pair. |
| Carbon salvage | [Carbon note](CARBON.md) | Exact finite spine and pruned-copying comparisons as reported in the sealed closure; other outputs are unreviewed diagnostics. |

All three B3F arguments use the edge-seeded undirected model and the previously established [existence of the exponent](../../../proofs/existence.md); no prior numerical bracket is an input. The finite-size products do not use exponent existence.

## Heavy builder and witness replay

Use a little-endian LP64 Linux host with GCC, GMP, OpenMP, GNU coreutils and Python 3.10+. Run these from the lane root:

```sh
cd round2/proofs/neon
B3F_OUT_DIR=./scratch/l4-18 B3F_THREADS=4 B3F_PROFILE=18 bash round2/cert/run_cert.sh
B3F_OUT_DIR=./scratch/l4-20 B3F_THREADS=2 B3F_PROFILE=20 bash round2/cert/run_cert.sh
B3F_L5_OUT_DIR=./scratch/l5-20 B3F_L5_THREADS=4 B3F_L5_PROFILE=20 bash round2/cert_l5/run_cert_l5.sh
```

Plan for at least 5 GiB RAM and 8 GiB free disk for D18; at least 14 GiB RAM and 24 GiB free disk for either D20 pair. Runs are sequential by sign and need minutes to hours. Peak recorded D20 L5 build RSS was 11,193,180 KiB, with about 2,000–2,145 seconds per build plus verification. Environment, compiler and native floating-point code generation can affect the numerical builders.

**Exact rebuild status:** The historical receipts and attacks authenticated and checked the original files. No byte-identical regeneration of any pair at the final source pin is established. The scripts use floating-point searches and print the resulting hashes; a GMP pass on rebuilt output verifies that output, not identity to the historical witness or its exact finite-size prefactors. Apply the exact checker to any new build; require the listed hashes before calling it a replay of the original. The L4 runners' override step accepts any nonzero status; for those controls inspect the verdict and require exit 1 with `VERDICT: FAIL`.

| Production file | Bytes | Historical SHA-256 |
|---|---:|---|
| `b3f18_lo.bin` | 3495091968 | `9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec` |
| `b3f18_hi.bin` | 3495091968 | `1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41` |
| `b3f20_lo.bin` | 11190881056 | `eb93985b2c90a516b65e85934f4ac847779f42b3fdbfcd2af4dc07fd30c37469` |
| `b3f20_hi.bin` | 11190881056 | `55dd6015eca85ca2fba6267d9444ec8807554c456999803d8f9abcc695c090e2` |
| `b3f20_l5_lo.bin` | 11370156264 | `9936d424fb3d87bc9193c79093870f1dd23772c16dd8b27a66233bc42a5a8474` |
| `b3f20_l5_hi.bin` | 11370156264 | `cd4553129c758b2c8f74e3f1686146a613e8b4f01f4f720a6a5b86da09b8d317` |

For an existing hash-matched payload, compile the corresponding verifier with `gcc -O3 -fopenmp SOURCE.c -lgmp -o verify`, then `./verify PAYLOAD --threads 4 --tree-check 30`. Follow with the paired `tree_identity_check.py PAYLOAD 15` (D18), `tree_identity_check.py PAYLOAD 10` (D20 L4) or `tree_identity_l5_check.py PAYLOAD 30` (D20 L5). The D24 L4 lower-only [receipt](round2/outputs/cert_b3f_24_lo.log) is a historical benchmark, not a claim or pair.

The historical L4 cold-audit scripts look for their named fixtures and executables in the working directory. Run `b3f20_sol_smoke.sh` from a writable copy of the lane root so its relative path to `round2/cert/negative_controls.py` resolves. The L5 cold-audit checker builds its helpers in a temporary subdirectory of a writable working directory.

## Intermediate exact check (about 35 minutes)

From the repository root:

```sh
bash round2/verify/neon/run.sh --slow
```

This needs GCC's ASan/UBSan libraries and GNU `time` in addition to the heavy-replay dependencies above. Allow 2 GiB available RAM and 3 GiB free temporary disk. The harness uses one OpenMP worker, builds four fresh witnesses (two signs each at D10 and D12), runs the exact GMP and Python tree checks and inward-override controls, and runs the D10 sanitizer and planted-defect suite. For D10 alone, run `bash round2/proofs/neon/round2/cert_l5/test_cert_l5.sh 10`. Use `TMPDIR` to select temporary storage. These are benchmark profiles and require no historical payload. See [benchmark endpoints and controls](round2/cert_l5/REPORT.md); this suite was not rerun during packaging.

## Auxiliary evidence

- [Core-tail argument](round2/open_questions/REPORT.md#2-both-tails-of-the-core-limit) and [disagreement appendix](round2/open_questions/DISAGREEMENT.md) use Astra 6's uniform renewal estimates, Astra 17's star-departure argument, and a reconstruction of Long's aggregate contraction. Claim 2's attacks do not cover this section. The quick replay checks the exact contraction/rate fixtures only. Later Astra 26–27 supersede the lane's early-label limits for hub probabilities and moments, and later Neon-2 forest results treat moving cutoffs.
- The [pure-inner](round2/level5_study/REPORT.md), [mixed-row](round2/mixed_l5_study/REPORT.md) and [coarse-endpoint](round2/mixed_l5_d7/REPORT.md) studies compare memory and slack. Their float-width diagnostics are not certificates; exact fixture checks concern local transfers and selected trees. Summaries use the nearby compact output records. For replay in a separate writable copy, start with the included `diag_slack.c` rather than regenerating it from the later hardened production builder. The final [L5 benchmark report](round2/cert_l5/REPORT.md) also gives D10/D12 exact-check evidence.
- [Closure diagnostics](round2/hierarchy/) include the IR-0 model, the finite mean-field/cavity calculations and the bare-chain correlation simulation. `python3 -B round2/hierarchy/verify_hkrz_exact.py` (NumPy) reproduces the historical Figure 6 fits from embedded public figure coordinates: linear-space NLS `0.56952` and log-OLS `0.56713`. `phase1_calib.py` needs mpmath; `phase1_ir0_hkrz.py` needs NumPy and mpmath. C diagnostics use GCC with `-O3 -fopenmp -lm`.
- [Midpoint simulation](round2/numerics/vr_b3f_prod.c): compile with `gcc -O3 -fopenmp round2/numerics/vr_b3f_prod.c -lm -o vr` from the lane root; `./vr b3f18_lo.bin b3f18_hi.bin 2000 2` uses the original Claim 1 pair. The [historical run](round2/outputs/vr_b3f18_t24.log) took about 52 minutes and needs memory for both inputs and the simulation. `0.5657236585 ± 3.8e-9` is a finite-age Monte Carlo estimate with jackknife standard error, without a proved bound on finite-age bias. `python3 -B round2/numerics/pslq_search.py 0.5657236585 3.8e-9` uses Python and PARI/GP for its `lindep` portion; the [search record](round2/outputs/pslq_final.log) implies no arithmetic exclusion for the true exponent.
- The [outgoing review archive](../../reviews/attacks/) covers Gemini's 16 audits of Astra and Fable. The independent scripts are under `round2/outputs/attack_*/` in this lane, except the audits of Astra 17, 19–20, 22–23, 26 and 27, which are archived with Astra's proofs.
