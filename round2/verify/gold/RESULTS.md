# Timed clean builds

Each side was built sequentially in its own fresh directory. The measured build used the two files extracted by `git show` from source commit `edd297185538b5ec399c04cb7fc80ebe2c07310a`. `sha256sum -c inputs.sha256` authenticated both extractions. The compile command, shown with local directory names, was:

```sh
g++ -std=c++17 -O3 -fopenmp -Wall -Wextra -Werror -pedantic \
  pinned/build_quintet.cpp -o pinned/build_quintet
```

From each clean side directory, the timed command was:

```sh
export LC_ALL=C OMP_WAIT_POLICY=PASSIVE OMP_DYNAMIC=FALSE
export OMP_NUM_THREADS=10 OMP_THREAD_LIMIT=10 OMP_PROC_BIND=FALSE
/usr/bin/time -v -o time.txt ../pinned/build_quintet \
  "$preset" "$side" "quintet_$side.bin" > build.log 2>&1
sha256sum "quintet_$side.bin"
```

| Preset | Side | Exit | Wall (s) | User (s) | System (s) | Peak RSS (KiB) | Bytes | Historical match |
|---|---|---:|---:|---:|---:|---:|---:|---|
| Gold1 | Lower | 0 | 178.81 | 591.95 | 6.74 | 2,455,920 | 1,588,521,616 | Yes |
| Gold1 | Upper | 0 | 188.46 | 607.24 | 7.10 | 2,456,132 | 1,588,521,616 | Yes |
| Gold2 | Lower | 0 | 533.21 | 2,634.71 | 43.95 | 8,791,816 | 6,344,902,168 | Yes |
| Gold2 | Upper | 0 | 505.02 | 2,496.13 | 42.84 | 8,791,748 | 6,344,902,168 | Yes |

All four builds matched their pinned SHA-256 digests, and `cmp` returned 0 against all four historical snapshots:

```text
81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986  gold1-lower.bin
8a36c3afdd2a8f4bbc03461a306af04594c46ada6ef037646261da4d364d9383  gold1-upper.bin
41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413  gold2-lower.bin
f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876  gold2-upper.bin
```

## Released source identity

`git show` of release commit `a2c768885cb0af90cbd3364da13bbf56d45c24da` yielded the two builder inputs in [inputs.sha256](inputs.sha256). A second comparison streamed the pinned and released Git blobs through `cmp`; all eight comparisons succeeded:

```text
build_quintet.cpp: identical
quintet_format.hpp: identical
verify_quintet.cpp: identical
verify_tail_bernstein.py: identical
cert_check.py: identical
audit_independent.py: identical
certificate_gold1.json: identical
certificate_gold2.json: identical
```

Verifier source was extracted from that release commit. The strict public wrapper was extracted from commit `cd2f3a83156ae5fac47a6321b241a41e48643a98`; its SHA-256 is `067ad87698d13cd2becc8c1eb98cdc2ba026611b4cdbeb2d3225c399d1597d46`. The source, receipts and Gold manifests were unchanged between those commits (`git diff --exit-code` on the pinned paths).

## Strict historical verification of the rebuilt pairs

From clean release Git sources, with each pair directory containing links to its two freshly rebuilt binaries:

```sh
export LC_ALL=C OMP_WAIT_POLICY=PASSIVE OMP_DYNAMIC=FALSE OMP_PROC_BIND=FALSE
export OMP_NUM_THREADS=4 OMP_THREAD_LIMIT=4 NEON2_THREADS=4
/usr/bin/time -v -o "$pair/verify-time.txt" \
  bash round2/verify/neon2/run.sh "$preset-check" "$pair" \
  > "$pair/replay.json" 2> "$pair/verify.log"
cmp "$pair/replay.json" "round2/proofs/neon2/quintet-cert/certificate_$preset.json"
sha256sum "$pair/replay.json"
```

Wall time includes compilation, the original verifier and exact combiner, the independent sampler, and the wrapper receipt comparison. The separate `cmp` confirms the emitted JSON. GNU `time` reports the largest individual-process RSS; the tree column sums concurrent RSS across the verifier process tree, sampled every 0.2 seconds via Linux `/proc`. The sampled value can count shared pages more than once and is not a strict upper bound.

| Pair | Exit | Wall (s) | Max RSS (KiB) | Sampled tree peak (KiB) | Exact receipt |
|---|---:|---:|---:|---:|---|
| Gold1 | 0 | 83.48 | 3,222,740 | 4,743,664 | Yes |
| Gold2 | 0 | 123.67 | 12,745,804 | 18,793,872 | Yes |

Gold1 receipt SHA-256: `d26d22a31d8f18790f517c451708a35b2721cbb4809613d1316c05df5071faeb`.

Gold1 exact endpoints: `2825181729593748861/4994114995300000000` and `439124781075521909/776183857175000000`, outward interval `[0.565702177913, 0.565748407438]`.

Gold2 receipt SHA-256: `5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284`.

Gold2 exact endpoints: `35083855339718511/62017003637500000` and `127993535256830043/226242120050000000`, outward interval `[0.565713486333, 0.565736986679]`. A separate `Fraction` check confirmed both intervals, both lower endpoints above `0.5657`, and the exact gap `142/251 − upper = 3697635659207/56786772132550000000 > 0` for Gold2.

The original verifier checks the finite rows and integer envelopes, signed drift and comparability. The Python stages verify the exact Bernstein tail, combine and round the bounds outward, and run the original 47-tree oracle per side. The separately executed independent sampler checks projections, selected rows and envelopes, the tail, and 142 adversarial trees per side. Gold1 verifies 191,375,747 finite rows and 29,879,919 envelope inequalities per side; Gold2 verifies 767,588,844 finite rows and 106,953,150 envelope inequalities per side. The exponent interval uses the separate existence theorem; this replay does not reprove it.

The four binaries occupy 15,866,847,568 bytes. The full run directory, including sources, logs, executables, linked pairs and receipts, used 15,867,756,544 allocated bytes after both verifications (`du -B1 -s`).
