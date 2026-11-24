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
| Gold2 | Upper | — | — | — | — | — | — | In progress |

Gold1 matched both pinned digests and `cmp` returned 0 against both historical snapshots:

```text
81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986  gold1-lower.bin
8a36c3afdd2a8f4bbc03461a306af04594c46ada6ef037646261da4d364d9383  gold1-upper.bin
41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413  gold2-lower.bin
```

Gold2 lower also matched its pinned digest and passed `cmp` against the historical snapshot.

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
