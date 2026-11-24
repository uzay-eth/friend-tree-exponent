# Neon Claim 1: pinned production rehearsal

Inputs are the two original `B3-FullH-18` payloads. Supply the lower file first and upper file second. The script checks their byte counts and SHA-256 digests before and after each full verification; it does not rebuild or edit them. The local Neon Git repository is read-only to this procedure.

Requirements: Linux x86-64, Git, GCC with OpenMP, GMP development headers/library, GNU time and taskset, Python 3 with NumPy, and ≥8 GiB of currently available memory at each verifier launch. Use an assigned idle CPU. `CPU` defaults to 3, `NEON_SOURCE` to `[local path]`, and `AUDIT_PYTHON` to the local Python used in the audit. `REPLAY_OUT` must be a fresh directory.

```sh
CPU=3 AUDIT_PYTHON=/path/to/python REPLAY_OUT="$PWD/outputs/rehearsal" \
  bash replay-production.sh b3f18_lo.bin b3f18_hi.bin
```

The build uses a detached `git clone --no-local --no-checkout` pinned to `cec3855cc52ab4265ce3e36b8a7a38f1ae257f77`. The clean clone's verifier and header are compared to their Git blob bodies and to the vendored pinned copy used only by the independent probe. Build products live outside the checkout. Both payloads are verified serially with `--threads 1 --tree-check 30`; each verifier process exits and releases its memory before the next sign starts.

Receipts: `session.log` records source and executable hashes, environment and full output; `lo.log` and `hi.log` preserve every checker row count, comparison, concrete-tree result and `/usr/bin/time -v` component per sign; `oracle.log` and `tail-oracle.log` preserve the independent checks. `audit.py` recomputes exact comparison constants and detects omitted recharge; `tail_probe.py` compares rational interpolation to GMP and detects a planted reversal of the T3 transfer coefficient. The package omits the bulky payloads, scratch clone and rebuildable executables.
