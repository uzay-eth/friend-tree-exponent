# Gold2 rebound numerical negative control: wrong env23 join

**Result:** The unmodified `b435b4c824e2ba5835f58b85e67a27de7d6a0d49` C++ verifier rejected a hash-identified, range-valid, size-preserving two-word mutation of the original lower Gold2 snapshot. It reached the numerical `env23` projection check and returned exit **2**, empty stdout, and exactly `ERROR: env23 verification failed` on stderr. No `--perturb` or header override was used.

## Fault and independent arithmetic

Swap the centers of two exact endpoint slots in the first `env23` block, for the H3 receiver `(D,c)=(2,2)`. This simulates joining envelope values to the wrong endpoint IDs. Both radii are zero and left unchanged. Both in-range centers remain present, their signed edits cancel, and the per-block sum, multiset, maximum magnitude and all range constraints are unchanged.

All offsets are zero-based bytes from the file start, including the 512-byte header. Words are signed little-endian 64-bit integers.

| Exact endpoint `(D,c)` | Center byte offset | Original | Mutant |
|---|---:|---:|---:|
| `(2,1)` | `6141599992` | `−93711837683` | `−175951565211` |
| `(2,2)` | `6141600000` | `−175951565211` | `−93711837683` |

`plan.cpp` uses the pinned layout and sparse reads of the unchanged `g` coefficients to calculate the two source differences `g((3,2),w)−g((2,2),w)`: `−93711837683` and `−175951565211` ticks respectively. Each fine endpoint projects entirely onto its corresponding coarse slot. Both original inequalities have `abs(error)=radius=0`; after the wrong join, **both** have `abs(error)=82,239,727,528 > 0`. This local arithmetic does not require a table scan or a physical-tree fixture.

## Authentication and replay receipt

Input directory `round2/.work/neon2-gold2-GE7v6RUF`; `hashes.log` reported both `quintet_lower.bin: OK` and `quintet_upper.bin: OK` before the run. The clean detached checkout's tracked status was empty afterwards.

| Artifact | SHA-256 |
|---|---|
| Original lower (upstream authenticated receipt) | `41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413` |
| Original upper (upstream authenticated receipt; not mutated) | `f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876` |
| Private lower mutant (freshly computed before verification) | `38f90ddc3940aa131312cc9f7685033bf054a49016501b4580aec28fd8f97e9d` |
| Stock `verify_quintet` executable | `6ca090355b6b250672a7fc9fcce824a62263e5c901cc99ba636bd39197b5c12e` |

Only a separate Btrfs reflink inode was writable; it was removed with its temporary directory after the run. File size stayed `6,344,902,168` bytes. The stock verifier computes SHA over the actual input buffer before structural/range and numerical checks; its error path emits no success JSON or digest. The digest above is the independent SHA of the exact mutant passed to that invocation. Therefore this is not an old-digest mismatch or a file-format rejection. It tests this numerical wrong join on the lower witness, not every transfer family or an independent upper mutation.

From `[local path]`, with `B=round2/.work/neon2-gold2-GE7v6RUF` and `C=round2/attacks/neon2-quintet/gold2/controls`:

```sh
taskset -c 6 nice -n 19 timeout 18s g++ -std=c++17 -O0 \
  -I "$B/source/round2/neon2/quintet-cert" "$C/plan.cpp" \
  -o [local path]
taskset -c 6 nice -n 19 timeout 190s /usr/bin/time \
  -f 'elapsed=%e cpu=%U+%S peak_rss_kib=%M' \
  python3 -Bu "$C/replay.py" "$B" [local path] [local path]
```

Actual replay output (planner's diagnostic folded onto one line by the terminal):

```text
available_kib 16746724
env23 block 0 receiver 2,2 slots 0,1 labels 2,1 and 2,2 physical source deltas -93711837683,-175951565211 swapped abs errors 82239727528,82239727528 vs radii 0,0
offset 6141599992 old -93711837683 new -175951565211
offset 6141600000 old -175951565211 new -93711837683
mutant_sha256 38f90ddc3940aa131312cc9f7685033bf054a49016501b4580aec28fd8f97e9d
verifier_exit 2 stdout '' stderr 'ERROR: env23 verification failed\n'
elapsed=83.68 cpu=37.53+13.75 peak_rss_kib=6351152
```

`df -h [local path]` showed 621 GiB free; the script required ≥12 GiB available RAM at entry and before verifier launch, with inherited affinity pinned to one CPU and `--threads 1`. One Gold2 payload reached the stock verifier, and it stopped at the first envelope family, before drift/cap checks. The initial attempt found the copied mode was read-only, exited before mutation/hash/verification, and was corrected by changing mode on the confirmed private inode only.
