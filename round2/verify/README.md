# Round-two verification

From the repository root on 64-bit little-endian Linux:

```sh
bash round2/verify/run_all.sh
bash round2/verify/run_all.sh --lane neon2  # one suite only
bash round2/verify/run_all.sh --slow        # Neon D10/D12 only; about 35 min
```

Use **CPython 3.12** (`PYTHON` can name its executable), Bash, GCC/G++ with `__int128`, C++20 and OpenMP, GMP development headers, and GNU command-line utilities. GCC 13 is a baseline; some optional replays use `taskset`, `time`, and sanitizer runtimes. On Ubuntu 24.04:

```sh
sudo apt-get install build-essential libgmp-dev python3.12-venv coreutils diffutils gzip time util-linux
```

The aggregate creates a fresh standard-library-only venv under `TMPDIR` (default `/tmp`), runs suites one at a time and caps OpenMP and BLAS at one thread. Use a writable temporary volume that permits executables. The default run installs no packages and fetches no inputs. The runner prefixes checker stdout with the suite name, stops on the first failed suite, and emits one individually labelled `PASS` per completed result in the success summary. The inherited `PYTHONOPTIMIZE` setting is cleared because some checkers use assertions; do not add `-O` or `-OO` to manual Python commands that rely on assertions. If a checker supplies its own optimized-Python test, keep its documented flags.

| Suite | What the default run checks | Replay guide |
|---|---|---|
| Astra | Claims 1–16 finite fixtures, including regenerated excursion rows and the native cap-48 scan; exact arithmetic over frozen ledgers for Claims 9 and 15. | [Astra](../proofs/astra/README.md) |
| Late Astra | Claims 17–24 and 26–28 exact fixtures, including Claim 21 arithmetic over frozen D18 receipts, archived lifetime comparisons, and independent checks. | [Late Astra](../proofs/astra-late/README.md) |
| Neon | Exact degree and disagreement transcripts, and finite Bernstein/endpoint algebra. The latter took about 64 s with one worker and < 13 MiB peak RSS. | [Neon](../proofs/neon/README.md) |
| Neon-2 | Finite identities and exact order, moment, and height fixtures; size-16 label-order receipt and tiny quintet/marginal checks. The size-16 enumeration takes roughly 2–6 min; budget 1 GiB RAM and minimal disk. | [Neon-2](../proofs/neon2/README.md) |
| Fable | Selected simulation-moment reduction, finite exact hub/transport/path fixtures, a density-repair fixture, and certified root counts for the IR-0 cousin (~43 s for IR-0). | [Fable](../proofs/fable/README.md) |

These checks do not replace the accompanying proofs, verify every inherited premise, or replay production rows absent from the clone. In particular, a successful tiny quintet or marginal run verifies the tiny candidate only. Selected simulation moments are numerical evidence, not an exact certificate for the exponent.

## Larger replays and historical receipts

Run the commands in each linked guide separately. Do not run production jobs concurrently. Times and resource figures below are historical measurements or planning budgets, not guarantees for a new host. Reserve additional scratch space for compilation and logs.

| Replay | Time, memory and disk | Extra requirements and scope |
|---|---|---|
| [Astra physical and joint enclosure](../proofs/astra/README.md#full-historical-replays) | The historical joint cold run took 53 min; its five-pass stage peaked at 2.74 GiB RSS under a 5 GiB address-space cap. Allow several GB of scratch, including two raw outputs of about 612 MB each. | GCC C++20 with UBSan, GMP, GNU `time`, Linux CPU affinity and NumPy for joint cold. The separate `spatial-full` replay adds the wider spatial cutoff and semigroup checks. |
| [Late Astra threshold](../proofs/astra-late/README.md) | The complete sweep's time and peak RAM have not been measured for this release. Run one label shard at a time as documented; no generated files are needed. | Python standard library. Opt-in, separate from the independent checks in the default run. The cold review reran only the first shard and the separate label-20 extension. |
| [Neon L4/L5](../proofs/neon/README.md#heavy-builder-and-witness-replay) | D18: budget 5 GiB RAM and 8 GiB disk. D20: 14 GiB RAM and 24 GiB disk; D20 L5 took 33–36 min per build plus verification. `--slow` runs the smaller D10/D12 controls: historically ~35 min, under 2 GiB RAM; allow 3 GiB scratch for four generated endpoints, copies and binaries. | GCC/OpenMP/GMP; the control benchmark also uses ASan, UBSan, and GNU `time`. The original production witness pairs are absent; builders produce candidates. |
| [Neon-2 quintet and stationary marginals](../proofs/neon2/README.md#production-sized-certificates) | Gold1: ~8 min build per side, ~3.5 GB check RAM, 3.2 GB pair. Gold2: ~25 min per side, ~12.5 GB check RAM, 12.7 GB pair. Marginals: ~1.5 min per file, ~1.5 GB check RAM, ~6.7 GB output. | GCC C++17/OpenMP. The [Gold rebuild guide](gold/README.md) covers digest comparisons. For new bytes, verify the freshly produced endpoints; the old receipt does not certify them. |
| [Fable large enumerations and simulations](../proofs/fable/README.md#opt-in-computations) | Size-23 hub enumeration: ~23 min and ~8.2 GB RAM in four shards. A sparse transport run for labels 20 and 21 took 4 min and 10.8 GB RAM. Most exact enumerations emit text logs; the historical raw simulation budget was ~45 core-hours, with disk depending on requested samples and features. | C++17/GMP (C++20 for the larger orbit checker). Reaggregating packaged moments does not reconstruct the unshipped raw runs. |

**Matching a historical receipt:** Verify the input digest and the exact checker verdict. A builder may return successfully with a different candidate; platform, compiler and thread schedule can change floating-point proposals. A new exact pass establishes only the bounds reported on those new bytes. For a negative control, require the documented rejection status and diagnostic: an unrelated file, compiler, or invocation error is not a rejected mathematical counterexample.

## Optional Python dependencies

For the optional commands that import NumPy, SciPy, mpmath, or SymPy, create a separate CPython 3.12 environment and install the [pinned requirements](requirements.txt). For example:

```sh
work=$(mktemp -d)
python3.12 -m venv "$work/venv"
"$work/venv/bin/python" -m pip install --only-binary=:all: -r round2/verify/requirements.txt
export PATH="$work/venv/bin:$PATH" PYTHON="$work/venv/bin/python"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
```

These pins use the NumPy/SciPy versions from the Neon-2 stationary LP pilot. For historical Fable numerical comparisons, use its [recorded versions](../proofs/fable/README.md#opt-in-computations) in a separate venv. Install PARI/GP (`pari-gp` on Ubuntu) only for the Neon inverse-symbolic diagnostic that calls `gp`. Remove the optional environment when finished.
