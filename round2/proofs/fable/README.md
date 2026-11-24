# Fable: proofs, certificates and numerical evidence

**Edition:** Release packaging index for Fable authors Claude Fable 5.2, Claude Opus 5.5, GPT-6.2 Sol Ultrafast, and Gemini 4 Neon. **Pinned Fable source:** `a18825ccdecfdb8d8465a89a7d99fc33734fa8e9`. Claim-specific reviewer families, check scope, and dependencies are listed in the table and the [Fable ledger](../../LEDGER.md#fable).

The [round-two ledger](../../LEDGER.md) records review scope and grades. The imported proofs and checker inputs are pinned in [SOURCES.tsv](SOURCES.tsv); proof headers name the author families. Public review files are linked below. Start from the repository root:

```sh
bash round2/verify/fable/run.sh
```

The default run uses Bash, Python 3.10+ (standard library only; tested on 3.12), GCC with C++17 and GMP/GMPXX headers and libraries (`libgmp-dev` on Debian/Ubuntu), plus standard shell utilities and SHA-256. It runs serially and uses temporary files under `${TMPDIR:-/tmp}`. The small-tree and statistical checks took 16 s and 290 MiB peak RSS including compilation; the IR-0 check adds roughly one minute depending on CPU load. Set `CXX` or `PYTHON` to select executables. No dependency on Git metadata.

| Result | Proof or evidence | Replay and limits |
|---|---|---|
| <a name="f1"></a>F1: final finite-window statistics | [Claim and conditional-bias algebra](CLAIM.md); [final receipt](data/final_output.txt); [Astra review](../../reviews/attacks/astra-on-fable-1.md) | Quick run reproduces the primary and late ratios, both independent rungs, run and leave-chunk errors, the four conditional bounds and two hold-out ANPQ controls. It also exercises the production generator on small trees. It does not re-run production or infer the bias hypotheses from agreement. |
| <a name="f2"></a>F2: law of the core-count limit | [Proof](newground/law-of-X/REPORT.md) and [required self-contained waiting argument](newground/law-of-X/referee/many-site-repair.md); [Astra](../../reviews/attacks/astra-on-fable-2.md), [Neon](../../reviews/attacks/neon-on-fable-2.md) | All six conclusions use Astra 6's hypotheses (L1)–(L4), including the refined one-leaf estimate in [the probabilistic reconstruction](checks/astra6-reconstruction/blind-proof-probabilistic-half.md). Lower-tail equivalent for the CDF; no limiting Weibull rate constant. |
| <a name="f3"></a>F3: hub probabilities by label | [Proof and first certificate](newground/q3-hub-probability/REPORT.md), [large extension](newground/q3-hub-probability/EXTENSION.md), [large program notes](newground/q3-hub-probability/ext-REPORT.md); [Astra](../../reviews/attacks/astra-on-fable-3.md), [Neon](../../reviews/attacks/neon-on-fable-3.md) | Quick replays the exact size-11 enumeration through label 8 and a separate smaller rooted DP. The two large exact programs agree through label 23; the final pair (23,24) is a directed-grid result from one program. Exact-order scaling uses the core-count mean theorem. Neon 2 first posted the degree dichotomy; F3 also proves it. |
| <a name="f4"></a>F4: promotion age and young fringes | [Proof §§1–2](newground/q5-low-degree-forest/REPORT.md); [Astra](../../reviews/attacks/astra-on-fable-4.md) | The age law uses positive a.s. core growth from Astra 6. The two fringe error bounds use post-promotion parent data and do not require that limit. Scope does not include the component law of the whole forest. |
| <a name="f5"></a>F5: stochastic domination | [Coupling, strictness, tails and local criterion](newground/stochastic-domination/COUPLING.md); [Astra](../../reviews/attacks/astra-on-fable-5.md), [Neon](../../reviews/attacks/neon-on-fable-5.md) | Quick exercises dense, sparse and local criteria at smaller sizes. Full transport range (2,3) to (19,20), one sparse result for (20,21), and local extension through (25,26) are distinct; full local range shares shape laws with its GPT check. No all-label result follows. |
| <a name="f6"></a>F6: density and endpoint | [Proof](newground/law-of-X/density/REPORT.md) with the [required bare-triple repair](newground/law-of-X/density/bare-triple-repair.md); [Astra](../../reviews/attacks/astra-on-fable-6.md) | Quick exercises exact star-exit and moment arithmetic and the original flag counterexample. Boundedness, smoothness, sharp Hölder order and the two endpoint equivalents use Astra 6, 17 and its additions, Neon-2 1 and F2. Stronger regularity has no full outside-family check. |
| <a name="f7"></a>F7: longest degree-two path and maximal low-degree component | [Proof](newground/q5-low-degree-forest/longest-path/REPORT.md) | Quick checks finite-size path laws. The almost-sure `O(1)` error and fixed-cutoff `Θ(log n)` result use Astra 6 and the analytic trial arguments; finite enumerations do not verify those limits. There is no full outside-family review. |

## Astra 6-7 reconstructions

The Claude-family lead's related write-ups are [uniform means (analytic)](proofs-by-lead/uniform-means-from-the-renewal-kernel.md), [almost-sure limit from uniform means](proofs-by-lead/almost-sure-limit-from-uniform-means.md), and [hub mass](proofs-by-lead/hub-mass-theorem.md). The analytic write-up follows reading the earlier arguments and uses a contour shift. Christopher D. Long's [public drafts](https://github.com/long-mathematics/friend-tree-growth) have priority for the critical generation-sum characterization, computability and the public proof of `E C_n = Θ(n^μ)`. The signed scalar agrees with Astra's above `μ`; the operators differ. Long's drafts do not contain the sharp mean equivalent, positive almost-sure amplitude or hub-mass theorem.

## Frozen numerical data

[DATA.sha256](DATA.sha256) authenticates the published text artifacts. `data/f1-moments.tsv` is a scope-limited, lossless projection of the pinned `moments.npz` (source SHA-256 `94f550688381015973d10e083cedd3546476d21d824d4d9e9e23841a839e3109`). For the five named cohorts it keeps the exact chunk counts, reference exponent, selected source column indices (zero-based), selected sums and corresponding full row-major cross-products. Each decimal is the shortest representation that round-trips to its original binary64 value. Chunk order and floating-point values are unchanged. It omits unused observables, raw paths, per-source filenames and the file-level block permutation data. The complete historical printed receipt is preserved, but only the statistics named above are replayed from the projection.

The five `theta_*.txt` files losslessly represent the little-endian binary64 frozen weights, in source order. To reconstruct a weight file for the historical simulator, run for each weight:

```sh
python3 -c 'import struct,sys; sys.stdout.buffer.write(b"".join(struct.pack("<d",float(x)) for x in open(sys.argv[1])))' INPUT.txt > OUTPUT.f64
```

The quick integrity check authenticates the text and tests the derived values. The terminal hub receipts `ext-{lastpair,merged23,grid24}.log` are the original text after decompression. Proof transformations add attribution and required cross-references, and scope F4 to its two published sections. `SOURCES.tsv` records the original unchanged blob hashes, including for transformed files.

## Opt-in computations

Set `P` to the absolute path of `round2/proofs/fable`, and run commands emitting files from a writable scratch directory. Do not compile exact checkers with `NDEBUG`.

**F3 large certificates** (C++17 with GMP; the orbit checker needs C++20). The last lead run historically took about 23 min and 8.2 GB. The four orbit shards were run serially (roughly 35 min total; peak 6.9 GB). Directed-grid size 24 is separate and should be budgeted for a large-memory machine; the historical run used a virtual-memory limit of about 11.2 GiB.

```sh
H="$P/newground/q3-hub-probability"
g++ -O2 -std=c++17 "$H/lead-signed_orbit_dp.cpp" -lgmpxx -lgmp -o dp
./dp 2 18 1
IMPROVED=1 ./dp 19 21 0
IMPROVED=1 SHARDS=4 ./dp 22 22 0
g++ -O3 -std=c++20 "$H/ext-dp.cpp" -lgmpxx -lgmp -o orbit
for i in 0 1 2 3; do ./orbit 23 23 1 22 23 4 "$i" > "sh23-$i.log" || exit; done
python3 "$H/ext-merge.py" sh23-0 sh23-1 sh23-2 sh23-3 > merged23.log
cp "$H/ext-dp.cpp" dp.cpp
g++ -O3 -std=c++20 -I . "$H/ext-grid.cpp" -lgmpxx -lgmp -o grid
./grid 24 > grid24.log
```

**F5 full certificates.** Sparse through (19,20) historically took about 2 min and 3.7 GB; (20,21) alone took 4 min and 10.8 GB. Local through (25,26) needs about 8 min and 6.3 GB. Check that each expected row reports `PASS`, since the local checker can report failure while returning exit zero. The original second transport program shares a pinned law implementation; the copy below only restores the relative include filename, as done in the cross-review. Compare its exact integers with the first sparse transport rows.

```sh
S="$P/newground/stochastic-domination"
g++ -O2 -std=c++17 "$S/lead-transport_sparse.cpp" -lgmpxx -lgmp -o transport
./transport 2 19
./transport 20 20
g++ -O2 -std=c++17 "$S/lead-local_scalar.cpp" -lgmpxx -lgmp -o local
./local 4 25
mkdir -p newground/stochastic-domination attacks/q3
cp "$H/lead-signed_orbit_dp.cpp" attacks/q3/signed_orbit_dp.cpp
cp "$S/child-transport.cpp" newground/stochastic-domination/
g++ -O2 -std=c++17 newground/stochastic-domination/child-transport.cpp -lgmpxx -lgmp -o transport2
./transport2 2 19
```

**Diagnostics and historical reproduction.** These commands use the original numerical methods, and none establishes the F1 bias premises. Python dependencies for the analysis scripts: NumPy 2.1.3, SciPy 1.14.1, mpmath 1.3.0 and SymPy 1.13.3 (Python 3.12).

- To run just the IR-0 interval argument-principle counts and negative control, use `bash round2/verify/fable/run.sh ir0`. This is the solvable cousin, not the friend tree.
- Final inverse-symbolic searches: `python3 "$P/src/isc_search.py" 0.5657236606 6.15e-09` and `python3 "$P/isc2/search.py" 0.5657236606 6.15e-09`. Compare [lead](data/isc_final_lead.txt) and [independent](data/isc_final_independent.txt) receipts. This is a finite family/height search. The [friend-tree kernel grid](data/kernel_scan_T17.txt) is a Monte Carlo diagnostic, not an enclosure.
- <a name="fable-fact-14"></a>Fig. 6 sampling-noise diagnosis: [protocol, simulation budgets and results](fit569/REPORT.md). Refit the included replicate means with `python3 "$P/fit569/analyze.py" "$P/fit569"` (the historical simulation took about 65 minutes on two cores).
- For historical F1 raw regeneration (about 45 core-hours total), build `src/ftprod.cpp` with GCC 13 `-O2 -std=c++17`, and reconstruct the five input weights as above. The exact per-chunk command is `ftprod STREAM FIRST_RUN RUNS OUT THETAS T1 T2 ... Rsum=3 Rpow=3 Dp=7 Dk=8 L1cap=8`. The [claim](CLAIM.md) records streams and checkpoints, the [historical drivers](src/prod.sh) and [late driver](src/late.sh) give chunk IDs and seeds, and [the reducer invocation](src/final.sh) gives sample selectors. As a two-row spot check: stream 9100, first run 50000, count 2, checkpoints 12 16 20 24, weight order sane/full/AN/ANPQ/AR; compare the prefix against [replay_reference.txt](data/replay_reference.txt). The raw files and historical binary are not included.
- Simulation diagnostics for F2, F6 and F7 use the compile/run commands in their reports. They regenerate sample files in scratch; no raw sample dumps are part of the proofs.
