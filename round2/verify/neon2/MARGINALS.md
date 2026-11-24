# Stationary marginals: clean rebuild and exact replay

All sixteen `high` snapshots rebuilt successfully and matched the historical SHA-256 digests. The strict release checker exited 0. Independent `cmp` against both the worktree receipt and the release Git blob returned 0. The complete 57,761-byte replay is therefore identical to the [checked-in exact JSON](../../proofs/neon2/stationary-marginals/certificates_all.json), including all eight exact intervals and six corollary records. Its SHA-256 is `cc8658e2e6354577716b44ff7f57bf6e3bb9831a1c7a7ac0ecf5c9bf717dbf9c`.

## Released source identity

The measured build compiled the release worktree's `build_marginal.cpp` with its two local headers. The builder inputs, verifier inputs and receipts below were byte-compared with the release Git objects in `c97ea780442e98396a62b810af4c1e4f681804f5` and matched the source-pin blob IDs in [`SOURCES.tsv`](../../proofs/neon2/SOURCES.tsv) for `edd297185538b5ec399c04cb7fc80ebe2c07310a`. The same blobs are in `cd2f3a83156ae5fac47a6321b241a41e48643a98`, which supplies the [strict release wrapper](run.sh) (SHA-256 `067ad87698d13cd2becc8c1eb98cdc2ba026611b4cdbeb2d3225c399d1597d46`).

| File | Git blob |
|---|---|
| `build_marginal.cpp` | `41034c1224aea59851235efde5e99f5c9408d57b` |
| `marginal_format.hpp` | `519f829d94e327fb8ced5b29bd73dc93a704dcc8` |
| `quintet_format.hpp` | `647a892cfbcbea11efd8c2bc5dc9898f8d564e06` |
| `verify_marginal.cpp` | `c3cd0f5c107a485114a1979ef28f04ec2cc81201` |
| `verify_marginal_tail.py` | `29b2bcce9481d7a2b741cb8e955eba5757939746` |
| `verify_tail_bernstein.py` | `461f48541fe968083ddf0bad9031ca63bf88d256` |
| `cert_marginal.py` | `6f46870b15518e8cc9a8e5b2b416951d919863b0` |
| `verify_all_marginals.py` | `cb9cf9f252de75796624eb737074b24d0e7dc0cc` |
| `certificates_all.json` | `0546888217fceddffe706f227b34f96400b430d5` |
| `certificate_gold2.json` | `56793d82ae8e2860f898a30194b86519d048507a` |

## Clean build and exact check

Run from the release repository root on Linux with Bash, GNU `time`, GCC with OpenMP, and Python 3. The measured run used x86-64 (AMD Ryzen 7 5800X), GCC 13.3.0 and Python 3.12.3. The build and check were each pinned to a single logical CPU with `nice -n 10`; CPU 15 was used in this run. Replace `cpu` with an allowed CPU number on your machine. The build used a fresh output directory; the commands below use portable output names. Run one production command at a time. Use a filesystem with at least 6.7 GB free for the snapshots plus logs and temporary verifier files.

```bash
set -euo pipefail
cpu=15
m=round2/proofs/neon2/stationary-marginals
out=$(mktemp -d "$PWD/marginals.XXXXXXXX")
export LANG=C.UTF-8
export OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 OMP_WAIT_POLICY=PASSIVE
unset PYTHONOPTIMIZE

taskset -c "$cpu" nice -n 10 g++ -std=c++17 -O3 -fopenmp \
  -Wall -Wextra -Werror -pedantic "$m/build_marginal.cpp" -o "$out/build_marginal"
for id in {1..8}; do
  for side in lower upper; do
    stem="obs${id}_${side}"
    taskset -c "$cpu" nice -n 10 /usr/bin/time \
      -f 'elapsed=%e user=%U max_rss_kib=%M status=%x' \
      "$out/build_marginal" high "$id" "$side" \
      565713486333 1000000000000 565736986679 1000000000000 \
      "$out/$stem.part.bin" > "$out/$stem.log" 2>&1
    mv -- "$out/$stem.part.bin" "$out/$stem.bin"
    sha256sum -- "$out/$stem.bin"
  done
done

TMPDIR="$out" NEON2_THREADS=1 taskset -c "$cpu" nice -n 10 \
  /usr/bin/time -f 'elapsed=%e user=%U max_rss_kib=%M status=%x' \
  bash round2/verify/neon2/run.sh marginals-check "$out" \
  > "$out/marginals-replay.json" 2> "$out/marginals-replay.log"
cmp -- "$m/certificates_all.json" "$out/marginals-replay.json"
```

The exact check uses the default rational exponent interval `[565713486333/10^12, 565736986679/10^12]`, with no threshold override or skipped tree oracles. Using exact fractions, the [pinned Gold2 receipt](../../proofs/neon2/quintet-cert/certificate_gold2.json) gives `[35083855339718511/62017003637500000, 127993535256830043/226242120050000000]`, contained in that interval. This check uses the exponent bracket as a premise; it does not rerun Gold2 or independently establish the spatial identification. One thread is a scheduling choice, not a change to the checked inequalities. The verifier checks both exponent endpoints, the finite motif rows and envelopes, the rational tail with 16 Bernstein panels, and 35 exact concrete-tree oracle cases per snapshot. For observable 8's upper bound it adds the analytic bare-core tail at cutoff 1600.

## Measured snapshot resources and digests

Per-file times come from GNU `time` around the builder alone; they exclude compilation, rename and hashing. The replay measurement includes compiling the C++ verifier, checking all snapshots, and the wrapper's receipt comparison. Expected SHA-256 values come from the historical Git receipt; actual values were recomputed from the final files with `sha256sum`, independently of the builder logs and verifier. The historical snapshots are not bundled, so snapshot identity here is checked by digest, not a direct binary `cmp` against historical copies.

| Snapshot | Exit | Wall (s) | User (s) | Peak RSS (KiB) | Bytes | Expected SHA-256 | Actual SHA-256 |
|---|---:|---:|---:|---:|---:|---|---|
| `obs1_lower.bin` | 0 | 125.90 | 119.25 | 626,104 | 413,129,728 | `cdac40bcfd62bac4e75be4790838b7be014fa8aba0ca5f983ccd0b233f36897a` | `cdac40bcfd62bac4e75be4790838b7be014fa8aba0ca5f983ccd0b233f36897a` |
| `obs1_upper.bin` | 0 | 113.33 | 110.82 | 626,108 | 413,129,728 | `97f0151a8ef79cbf2f294a51e728c79a9e2b5e5045cd878b794789c794556ed1` | `97f0151a8ef79cbf2f294a51e728c79a9e2b5e5045cd878b794789c794556ed1` |
| `obs2_lower.bin` | 0 | 121.02 | 119.12 | 625,980 | 413,129,728 | `9c8c5fc1c29c8e2552b6b8cb8677a07e70e255ba19d0a01b497901ff946c3024` | `9c8c5fc1c29c8e2552b6b8cb8677a07e70e255ba19d0a01b497901ff946c3024` |
| `obs2_upper.bin` | 0 | 115.60 | 114.66 | 626,108 | 413,129,728 | `52dd0dddcc373f57d4e5fb65b3e72a8c24f2bf4155a40f4478fb74871ebee8b9` | `52dd0dddcc373f57d4e5fb65b3e72a8c24f2bf4155a40f4478fb74871ebee8b9` |
| `obs3_lower.bin` | 0 | 124.35 | 120.49 | 625,980 | 413,129,728 | `cd2a23693205ef26e460b9a9be7d91cfc4c70eaf443b6f261db5a80076636e31` | `cd2a23693205ef26e460b9a9be7d91cfc4c70eaf443b6f261db5a80076636e31` |
| `obs3_upper.bin` | 0 | 119.75 | 118.25 | 626,108 | 413,129,728 | `0e2ecae38744c5d3cdd9b7ef3854ae2788ec061cbcc4e9c165956bc29c8e7b1e` | `0e2ecae38744c5d3cdd9b7ef3854ae2788ec061cbcc4e9c165956bc29c8e7b1e` |
| `obs4_lower.bin` | 0 | 111.85 | 110.83 | 625,980 | 413,129,728 | `201c494cba42437c2694988e4d429faf208d06557678f2d522bad92dfc01254a` | `201c494cba42437c2694988e4d429faf208d06557678f2d522bad92dfc01254a` |
| `obs4_upper.bin` | 0 | 111.82 | 111.06 | 625,976 | 413,129,728 | `bab3a00aba721a49e2ce0fd5e1e316328cfcbf6042420ad6526daf03c91d3b61` | `bab3a00aba721a49e2ce0fd5e1e316328cfcbf6042420ad6526daf03c91d3b61` |
| `obs5_lower.bin` | 0 | 114.90 | 114.15 | 626,108 | 413,129,728 | `898409d6d371275959f50bdbea3192b0225af61bebbece70253ab55ddeaa384d` | `898409d6d371275959f50bdbea3192b0225af61bebbece70253ab55ddeaa384d` |
| `obs5_upper.bin` | 0 | 114.94 | 114.19 | 626,108 | 413,129,728 | `a08e8e962d4c1902e1e635005794049b7f977e708ae7d0b9b87869127f3158cb` | `a08e8e962d4c1902e1e635005794049b7f977e708ae7d0b9b87869127f3158cb` |
| `obs6_lower.bin` | 0 | 107.68 | 106.86 | 626,108 | 413,129,728 | `30b888d29e78cce80bf002a478a0132a66e3907801c3861bbc4be22e025c3074` | `30b888d29e78cce80bf002a478a0132a66e3907801c3861bbc4be22e025c3074` |
| `obs6_upper.bin` | 0 | 90.92 | 90.20 | 625,976 | 413,129,728 | `78b2fe0869ec6cc92f4ee2dd93b7657bfb2195183d92616edb6dc8b7a2bcd5f1` | `78b2fe0869ec6cc92f4ee2dd93b7657bfb2195183d92616edb6dc8b7a2bcd5f1` |
| `obs7_lower.bin` | 0 | 88.43 | 87.69 | 625,980 | 413,129,728 | `03eb2f380a0be47cb2864acfeaa9463fc38255e261cfe6043948463c7ac02a28` | `03eb2f380a0be47cb2864acfeaa9463fc38255e261cfe6043948463c7ac02a28` |
| `obs7_upper.bin` | 0 | 82.01 | 81.34 | 626,108 | 413,129,728 | `7d15c18d869bc6aac50c9741beb3d273f0e0e9ca03086c8f68f74c0670b4dc32` | `7d15c18d869bc6aac50c9741beb3d273f0e0e9ca03086c8f68f74c0670b4dc32` |
| `obs8_lower.bin` | 0 | 78.65 | 77.97 | 625,976 | 413,129,728 | `dacdf9c26e1935445f1b80be6c8ce4fc54c30203d5d37bc4182183bd19db3afe` | `dacdf9c26e1935445f1b80be6c8ce4fc54c30203d5d37bc4182183bd19db3afe` |
| `obs8_upper.bin` | 0 | 79.86 | 79.18 | 625,976 | 413,129,728 | `f2ac7370f2c74179ed92c5354331e276aac4ee4b44c0f020c2eab1c615002ab2` | `f2ac7370f2c74179ed92c5354331e276aac4ee4b44c0f020c2eab1c615002ab2` |

**Build totals:** 1,701.01 s summed wall time (28 min 21 s), 1,676.06 s user CPU, and 626,108 KiB highest per-file peak RSS. The sixteen binaries total 6,610,075,648 logical bytes (6.61 GB, about 6.16 GiB) and occupied 6,610,092,032 allocated bytes on the measured filesystem. The build directory, including the executable and logs, occupied 6,610,399,232 allocated bytes. A partial snapshot is renamed in place without making another copy.

**Full strict replay:** exit 0, 149.77 s wall, 142.39 s user CPU, 1,658,036 KiB peak RSS reported by GNU `time`. This RSS is not the sum of the simultaneously live verifier processes. Each snapshot contributed 49,204,104 finite motif rows checked at both endpoints and 8,035,467 cross-level envelope checks shared by the two endpoints.
