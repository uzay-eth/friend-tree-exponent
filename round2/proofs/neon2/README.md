# Neon-2 proofs and replay

**Authors:** Neon-2 lane (Gemini 4 Neon lead with GPT-6.2 Sol Ultrafast contributions).  
**Pinned source:** `edd297185538b5ec399c04cb7fc80ebe2c07310a`.  
**Review and dependencies:** The [ledger](../../LEDGER.md) records actual GPT, Claude and Gemini review scope by claim and links to all 23 rival attacks. The [inherited proof inputs](#inherited-proof-inputs) are listed below by claim.

[SOURCES.tsv](SOURCES.tsv) records each imported file's original path and Git blob. Code, exact inputs, and receipts are byte-for-byte imports; proof copies have attribution, link, math-rendering, and replay edits.

From the repository root run `bash round2/verify/neon2/run.sh`. It needs Bash, Python 3.10+ (standard library only), a 64-bit little-endian environment, GCC with C++17, `__int128` and OpenMP, and standard utilities (`mktemp`, `cmp`, `cat`, `rm`). `PYTHON`, `CXX`, `TMPDIR` and `NEON2_THREADS` can be set; the last defaults to two and also caps OpenMP unless `OMP_THREAD_LIMIT` is set explicitly. The default run checks the finite fixtures for the analytic proofs, matches the size-16 label-order certificate, and builds and verifies tiny proposals for Claims 6 and 10 with mutations. Allow 1 GiB RAM and several minutes. The tiny proposals are smoke tests, not the published bounds.

## Claim index

In the checker column, paths are relative to this directory. Where cited, inherited results are premises; see the [Astra proof index](../astra/README.md) and [late Astra proof index](../astra-late/README.md). Finite test suites exercise generator identities, arithmetic and controls. The written arguments supply the unbounded and asymptotic steps.

| Claim | Proof | Checker and scope |
|---|---|---|
| 1 | [Absolute continuity and positive lower semicontinuous density](limit-density/README.md); [proof details](limit-density/REPORT.md) | `limit-density/check_density.py`: exact conditional clocks and split fixtures, plus a Beta simulation. The sender-rate control in the original checker is vacuous; see the independent review for a working asymmetric test. The no-kill weighting sends inert null sets to surviving null sets; the converse can fail. The extended density may be infinite on a null set and no upper regularity follows. |
| 2 | [Means and hub probabilities through label 8](q3/REPORT.md) | `q3/symmetries.py` and the size-16 exact receipt check (also available by itself with `bash round2/verify/neon2/run.sh q3-full`; roughly 2–6 minutes, allow 1 GiB RAM). The size-14 receipt is historical, with a known diagnostic-line difference on replay. `sim_q3.cpp` and `mc_t13.csv` are finite-time diagnostics. |
| 3 | [Explicit fringe law, rerooted eigenfunction, count limits and summability](claim14-proofs/README.md) | `claim14-proofs/verify_track3.py`: signed rows, rerooting, bounded-tree attachments and stopped-prefix fixtures. Uses repaired `ea6a093c`, not the original contour-independence assertion. Only intrinsic count/mean limits are contour-free; identifying the named spatial marginal uses the spatial theorem. The hub-mass summability result is self-contained. |
| 4 | [Response transcendence](obstruction/REPORT.md) §B | `obstruction/verify_transcendence.py`: finite coefficients and prime-valuations, not an exhaustive kernel census. The path arguments address the general coefficients. Transcendence of $`U_s(J)`$, $`F(s)`$ and every $`K_{ij}(s)`$ says nothing about the arithmetic of $`\mu`$. |
| 5 | [Simultaneous forest cutoff and separated TV regimes](q5/README.md) | `q5/verify_q5_exact.py`: exact tests for marked scores, holding times and boundary attachments. `q5/sim_q5.cpp` is a diagnostic. No sharp TV transition, fixed-shape intensity or TV theorem for all intermediate cutoffs. |
| 6 | [Unoriented quintet growth-rate proof](quintet-cert/PROOF.md); [exact bounds and digest table](quintet-cert/README.md) | `quintet-cert/verify_quintet.cpp`, `verify_tail_bernstein.py`, `cert_check.py`, `audit_independent.py`. Gold1 and gold2 are **separate historical witnesses**. Gold2 excludes $`142/251`$, not $`99/175`$. Exponent existence is a separate premise. See production commands below. |
| 7 | [Elementary density regularity](density-lp/README.md) | `density-lp/verify_density_lp.py`: exact packing, exponents, rates and barrier; floating-point probes are diagnostic. In the time-shift notation $`F_1(r)`$ is the first side's core count with the first wait removed. Global $`L^{19/18}`$ uses Astra's left tail; the elementary route does not prove boundedness or smoothness. |
| 8 | [All-label mean ordering](mean-order/README.md) | `mean-order/cert_mean_order.py`: 13,299 exact rows, 24 trees and ten controls. No stochastic dominance. |
| 9 | [All-label real moments](all-moments-order/README.md), $`p\ge1`$ | `all-moments-order/cert_all_moments.py`: all 19 boxes on $`[1,12]`$ and analytic tail checks. The proof handles the unbounded range. |
| 10 | [Eight stationary marginals and corollaries](stationary-marginals/PROOF.md); [bounds and digest table](stationary-marginals/README.md) | `stationary-marginals/verify_marginal.cpp`, `verify_marginal_tail.py`, `cert_marginal.py`, `verify_all_marginals.py`. **Historical** production receipt for sixteen exact witnesses. Both endpoints use original gold2 bracket and the spatial identification; neither the builder nor a tiny replay verifies the historical intervals. |
| 11 | [Height, diameter and core-height tail](height-diameter/README.md) | `height-diameter/cert_height.py`: exact upper drift and embedded CMJ rows with controls. Lower bound uses Biggins' branching-random-walk speed theorem and the repaired two-branch survival argument. No single almost-sure limit constant is asserted. |
| 12 | [All-label fractional moments](fractional-moments-order/README.md), $`1/2\le p\le1`$ | `fractional-moments-order/cert_fractional_moments.py`: all 14 boxes, 664,349 exact rows. Independent of Astra 27, which covered all positive powers earlier on 24 November 2026 (16:15 UTC versus Neon-2's 17:06 UTC commit; priority corrected 17:45 UTC). Together with Claim 9 this lane covers $`p\ge1/2`$. |

## Inherited proof inputs

Paths below are the original citation paths, located through the Astra indexes above. The linked Astra archive uses a later `round2/operator.md` that supplies the locally finite exterior approximation in its one-edge estimate and corrects the orientation notation of a banded transfer. The other listed inherited proofs trace to the same source blobs. These inputs require their own proofs and checks; a replay here does not re-certify them.

| Claims | Input and exact scope |
|---|---|
| 1 | Addario-Berry et al. Thm. 3.1; Astra 6 in `round2/standalone/claim-6.md` for the positive finite count limit; `round2/fringe/hub-children.md` §§1, 3 for the cover. The submitted positive-density route additionally uses the inert mean/stock/gradient in `round2/open/forest-age/inert/README.md` §§1–3 (or the Neon-2 grafted gradient under the corresponding normalization). |
| 2 | Addario-Berry et al. for existence of the limiting label shares. No Astra inputs. |
| 3 | `proofs/existence.md` Thm. 1 and Lemma 2.1, gold1 bracket from Claim 6, cap-48 signed contraction in `round2/influence/zero/pair/README.md`, `round2/fringe/finite-kernel.md` §§1–2 and `round2/finite/gauge/physical-trigger.md` §§1, 4 for the real positive semigroup and trigger, and `round2/standalone/claim-6.md` §§2, 4 plus `round2/open/forest-age/inert/README.md` §2 for the stated positive and separation inputs. `round2/operator.md` and `round2/fringe/spatial-limit.md` enter only the naming of the spatial fringe marginal in Theorem 1(2). Theorem 3 requires no Astra result. |
| 4 | The right-half-plane resolvents: `round2/influence/zero/critical.md` (Astra 14), `round2/finite/gauge/physical-trigger.md` (Astra 12) and `round2/fringe/finite-kernel.md` (Astra 5). |
| 5 | `round2/open/forest-age/README.md` (sequential limit); `round2/fringe/{almost-sure,finite-kernel,renewal,hub-spectrum,hub-children}.md` (core/size asymptotics, spectral gap, dirty-core bounds and hub comparisons). |
| 6 | [Exponent existence](../../../proofs/existence.md) Thm. 1 connects the independently bounded lower/upper logarithmic growth rates to the same $`\mu`$. |
| 7 | Repaired Neon-2 Claim 3 intrinsic limit, gold1, `proofs/existence.md` §§1, 3, and the Claim 1 clock construction. Corollary 2 adds `round2/open/limit-law/left/README.md`, which inherits a spatial/contour estimate. |
| 10 | Original gold2 bracket; repaired Claim 3 coboundary and positive mean growth; Astra 4/10 in `round2/operator.md` and `round2/fringe/spatial-limit.md` for the identification with the uniform-core fringe law and the spatial conditional expectations. The one-sided fringe by itself does not give the exterior parent's degree. The Cesàro alternative still requires positive limiting discounted mean growth. |

Claims 8, 9 and 12 use no inherited results. Claim 11 needs no Astra result; see [Biggins, *Branching out*, Thm. 3.1](https://arxiv.org/abs/1003.4715) for the external speed theorem.

## Production-sized certificates

The historical Git receipts record exact output and SHA-256 input digests. The multi-gigabyte inputs are absent from this archive. Builders use a fixed-seed empirical centering step and floating-point Bellman/optimizer sweeps, followed by integer quantization and envelope export. Compiler, OpenMP scheduling and reduction order can change the quantized values. A successful build gives only a candidate; a successful exact check on different bytes establishes only its own resulting endpoints. A byte-for-byte comparison with the historical receipt is meaningful only after actual verification. The original gold1 numerical check covered 191,375,747 finite rows and 29,879,919 envelope inequalities per side; gold2 covered 767,588,844 and 106,953,150 per side. See the separate gold1 and gold2 sections of the [GPT](../../reviews/attacks/astra-on-neon2-6.md) and [Claude](../../reviews/attacks/fable-on-neon2-6.md) attack reports.

Run one production command at a time, with output directories on a volume with sufficient space. Commands are from the repository root. For a check-only run with snapshots obtained separately, supply their directory instead and confirm that the receipt matches.

```sh
runner=round2/verify/neon2/run.sh
bash "$runner" gold1-build "$HOME/neon2-work/gold1"
bash "$runner" gold1-check "$HOME/neon2-work/gold1" > gold1-replay.out
bash "$runner" gold2-build "$HOME/neon2-work/gold2"
bash "$runner" gold2-check "$HOME/neon2-work/gold2" > gold2-replay.out
```

Per side, gold1 writes 1.59 GB and took about 8 minutes to build on the original setup; gold2 writes 6.35 GB and took about 25 minutes per side. The combiner can need about 3.5 GB (gold1) or 12.5 GB (gold2) peak RAM, plus enough disk for both files. A rebuild and a check both run lower and upper sequentially. For the published bracket use the exact rational endpoints in the matched receipt, not the nominal header exponents. The aborted gold3/gold4 investigations are not certificates.

```sh
runner=round2/verify/neon2/run.sh
bash "$runner" marginals-build "$HOME/neon2-work/marginals"
bash "$runner" marginals-check "$HOME/neon2-work/marginals" > marginals-replay.out
```

The marginal builder takes `high ID SIDE 565713486333 1000000000000 565736986679 1000000000000 OUT` for `ID=1..8` and `SIDE=lower,upper`; the wrapper builds all 16 without changing the original exponent interval. Plan for at least 6.7 GB of output and roughly 1.5 GB peak RAM at replay; the original build took roughly 1.5 minutes per file. The builder's secant/optimization loop chooses each nominal observable threshold; the verifier, not the optimizer, supplies rigorous margins across both endpoints. If a fresh gold2 result does not imply the original outward interval, it does not by itself discharge the numerical premise of the listed marginal bounds. The default quick run exercises a tiny pair for observable 1; this verifies the build-to-check route and passes the bracket as a premise, but does not establish that the pinned `high` builder reproduces the sixteen historical digests or intervals.

## Diagnostic dead ends

- [Stationary LP pilot](stationary-lp/REPORT.md): the optional simulations and numerical solves use NumPy 2.3.5, SciPy 1.16.3 and GCC C++17. The standard-library `stationary-lp/verify_farkas.py` checks the two `cert-control-*.json` witnesses by regenerating the sparse matrices from the model code. They exclude only the test interval $`[1/10,1/5]`$; small-cutoff LPs did not improve the prior.
- [All-diameter local eigenfunction obstruction attempts](obstruction/REPORT.md) §A: `obstruction/verify_obstruction.py` exercises the listed failures of the induction. The theorem remains open for diameter at least two.
- [Planted Bellman quotient diagnostics](planted-enclosure/README.md): build each C++ program with `g++ -O3 -std=c++17`; `planted-enclosure/test_transfer_clean.py` uses NumPy and SciPy. These are floating-point probes and extrapolations, not certified exponent bounds.
