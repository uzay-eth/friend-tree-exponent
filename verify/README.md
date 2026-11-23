# Verification Suite

Run every finite exact-arithmetic check from the repository root:

```bash
bash verify/run_all.sh
```

The script creates an isolated virtual environment at `verify/.venv` from the pinned Python dependencies in [`requirements.txt`](requirements.txt) (`numpy==2.1.3`, `mpmath==1.3.0`, `scipy==1.14.1`), compiles the C/C++ checkers into a temporary directory, runs all seven verification suites (including negative controls), cleans up temporary build artifacts, and prints one `PASS` line per result.

## Tested toolchain

- **OS:** Ubuntu Linux (`6.8.0-149-generic`, x86_64)
- **C / C++ compiler:** GCC / G++ `13.3.0` (with OpenMP and `-fsanitize=undefined`)
- **GMP:** `libgmp-dev 2:6.3.0+dfsg-2ubuntu6.1`
- **Python:** `Python 3.12.3` with `venv` and `pip`
- **Compression:** GNU `gzip 1.12`

On Debian/Ubuntu, the required system packages are `build-essential`, `libgmp-dev`, `python3-venv`, `python3-pip`, and `gzip`.

## Contents

| Directory | Result verified | Exact arithmetic |
|---|---|---|
| `motif-bracket/` | Sharp two-sided bracket `[113/200, 1133/2000]` on `E[C_n]` (`< 17/30` by `1/6000`), exact `n=3..21` cross-check, 47 rooted trees, 6 negative controls | GMP (`cert_verify_gmp.c`) + Python `Fraction` |
| `banded-upper/` | Independent upper bound `5657/9983` (`< 17/30` by `1/299490`), `226,545,617` rows + `D >= 2400` tail | 128-bit integer arithmetic (`check.cpp`, UBSan) |
| `banded-upper/independent/` | Second verifier written from the specification (`verify_banded.cpp`), smoke `3/5`, 436 exhaustive histories, 50 adversarial trees, 14 negative controls | 128-bit integer arithmetic (`verify_banded.cpp`, UBSan) + Python `Fraction` |
| `transfer-bounds/` | Transfer bracket `[563/1000, 568/1000]` (`9,639,742` states + `D >= 208` tail) and hub bracket `[56/100, 57/100]` | 128-bit integer arithmetic (`transfer.cpp`, `hub.cpp`, UBSan) |
| `baseline/` | Earlier baseline certificates `543/1000` (`1,205,767` checks) and `293/500` (`319,086` checks) used by the almost-sure lower and polynomial-scale lower proofs | Python `Fraction` (`certs9.py`, `certs6.py`, `find_upper_lin.py`) |
| `cross-charging/` | Almost-sure lower hypotheses (H1)–(H5) on `543/1000`, cross-charged drift identity on 200 trees, and `[11/20, 29/50]` certificates (`5,853`/`5,851` checks) + negative controls | Python `Fraction` (`as_lower_check.py`, `xcharge.py`, `xverify.py`) |
| `polynomial-lower/` | Exact local drift, exterior invariance on 84 saturated branches, 102 trees, 80 grafts, and independent cold audit with 14 negative controls | Python `Fraction` (`verify_as_lower.py`, `check_as_lower.py`) |

## Exhaustive drift checks on actual trees

From the repository root run the following command.

```bash
bash verify/tree-drift/run.sh
```

The motif test constructs each non-isomorphic rooted tree of sizes three through 14 with a core root. It attaches a leaf at every possible target, re-evaluates the potential on the resulting tree and sums the exact rational increments with the target's attachment probability. The two certificate directions are tested on 32,972 trees.

For the banded potential a rooted shape does not determine the labels. Starting from the three-vertex path, the test branches over every attachment target at every step through size nine, evaluating the true one-step drift at each of the 23,116 histories. The independent history implementation tracks the labels. A read-only observer retrieves the resulting weights from the production verifier. All finite rows must still pass that verifier before it returns the weights. Neither test uses the certified drift decomposition to calculate the true drift. The two optional arguments select the maximum motif size and maximum banded history size, respectively.
