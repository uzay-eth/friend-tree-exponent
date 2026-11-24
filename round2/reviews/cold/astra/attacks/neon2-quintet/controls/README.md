# Claim 6: numerical negative controls and arithmetic checks

On the clean detached source at `b04cc50c1e391e0420eb7b006ac781d484c68543`, the unmodified production executable rejected two private raw mutations at the **first numerical envelope family**. Independently reconstructed row/envelope counts and signed endpoint arithmetic match the owner's `gold.json` on both sides. No arithmetic or coverage discrepancy was found in these targeted checks.

The executable used was `round2/.work/neon2-claim6-7M6Ys3rJ/verify_quintet`, SHA-256 `dc5312577c2d96673ae24edac7231b8fdd43d392bdb881fa43640d0fe63a2ed9`. Its two clean input copies are siblings of that executable. No changes were made to those inputs or the detached checkout.

## Raw rebound controls

`offsets.cpp` uses only the pinned layout header to locate two `g` coefficients and the corresponding env23 midpoint/radius. Independently, the first internal row is `(D,c)=(2,2)`; sender `(2,1)` is exact in both spaces, so the envelope test reduces to `|g((3,2),(2,1)) − g((2,2),(2,1)) − midpoint| ≤ radius`. The original inputs `[old, gain, midpoint, radius]` were `[121598525227, 25376042183, -96222483044, 0]`, giving equality `0 ≤ 0`.

| Raw field in private lower copy | Byte offset | Replacement | Exact `lhs > rhs` | Actual verifier |
|---|---:|---:|---|---|
| `g((2,2),(2,1))` | 16700976 | `-17592186044416` | `17713784569643 > 0` | exit **2**, `ERROR: env23 verification failed` (23.98 s) |
| env23 midpoint for the same row/sender | 1531698976 | `35184372088832` | `35280594571876 > 0` | exit **2**, `ERROR: env23 verification failed` (20.87 s) |

Both replacements are at inclusive allowed magnitude bounds, so they pass the structural checks. The new expected bindings were `7b946eda8bfd7786a6d89297433e57cb94e6d95db954ab0bf9779aec845ce008` and `e5aa36ae0b24cdefe90f6eb9b249ed4e71f7555f401f443d97f4bec1244cb819`, respectively. The stock C++ executable internally computes the new input hash but has **no expected-hash argument**; the surrounding Python wrapper normally compares the emitted digest on successful runs. These tests passed the freshly hashed bytes directly to the executable, with no stale-hash gate or `--perturb` argument. On numerical failure it deliberately emits no stdout JSON. Each private copy was deleted after its run. The verifier reads and checks basic magnitudes over the snapshot, but exits at env23 before the remaining envelope families or any drift table scan.

```sh
R=round2/.work/neon2-claim6-7M6Ys3rJ
C=round2/attacks/neon2-quintet/controls
T=$(mktemp -d neon2-controls.XXXXXXXX)
g++ -std=c++17 -O2 -I "$R/source/round2/neon2/quintet-cert" "$C/offsets.cpp" -o "$T/offsets"
taskset -c 5 nice -n 19 timeout 95s /usr/bin/time -f 'wall=%e rss=%M exit=%x' \
  python3 -B "$C/mutations.py" "$R/quintet_lower.bin" "$R/verify_quintet" "$T/offsets" "$T"
```

Outer control exit **0**, `wall=66.02 rss=1604444`, after asserting both exact verifier failures. The offsets helper emitted `1588521616 16700976 16701016 1531698976 1533839968`.

## Independent coverage and signed global arithmetic

From the 512-byte headers alone, `coverage.py` enumerates the admissible internal descriptors, their degree-sum masks, and canonical symmetry multiplicities, independently of the C++ layout routines. It compares those counts to the owner's already-produced receipt. Both modes match all nine finite row counters and all three aggregate envelope counters:

```
f_in=4328 f_tab=1996669 g=10167794 h3=86843988
p4=29674890 s3=17410176 p5=15108750 b5=26374140 s4=3795012
envelope 1d=5747301 2d=20036106 3d=4096512
```

Cross-level projection checks use ordered blocks even where storage identifies reflection, and unordered tuples only within the declared symmetric slots. C++ checks all finite `f` rows including the analytic seam; `(INF,INF)` is omitted from the drift scan for `g` and separately forced to zero. The Python tail bound covers both gain-rate extremes with a signed polynomial on the full outer rectangle. Its concrete checks include the seam and leaf-stock extremes. This coverage check confirms counts, not a second independent scan of every row value.

Lower signed drift uses `−(ℒ−λ)`, upper uses `+(ℒ−λ)`, and ceilings are taken *after* introducing that sign. Positive comparability is required before dividing the nonnegative residual-per-core ceiling by `m_min`; the resulting error is subtracted on the lower side and added on the upper. Allocating `g` to the child gives `C−1` edges. The lower comparison recharge uses `Σ(c−2)=−2` for a nontrivial core tree, so positive recharge yields a conservative lower bound; `C=1` uses its own unshifted case. The independent integer/rational composition from the owner's receipt produced:

| Side | `m_min` | residual ticks/core | Exact endpoint |
|---|---:|---:|---|
| Lower | 99882299906 | 19186 | `2825181729593748861/4994114995300000000` |
| Upper | 155236771435 | 19783 | `439124781075521909/776183857175000000` |

`verify_quintet` emits valid residual caps, which may be positive; exit 0 from it alone is not an assertion that a witness satisfies the nominal trial. The wrapper's corrected endpoints are essential.

## Signed 128-bit corner

`overflow.py` uses Python unbounded integers and checks the generic maxima `DF=4096`, `|weight|≤2^44`, `|env|≤2^45`, `Q≤2^40` and sub-tick scaling 1024. At the lower `f` seam, bound the three quadratic terms and the linear term separately, then multiply by the maximum twice-rate and sub-tick scale, including the ceiling addition. This gives a **93-bit** absolute intermediate bound `7428244698622088964015652863`. In the planted extreme algebraic corner at `D=4095,c=0`, oriented gain has numerator after sub-tick scaling `7426431162318714430543626240`, denominator `33554432`, directed sub-ticks `221324895689449144320` and final ticks `216137593446727680`, within `int64_t`. Its opposite sign confirms correct negative ceiling. Additional checked routes (killing term, centered gain, 3D envelope difference and envelope source) use at most 95 bits with rounding allowance.

The independent global motif multiplicity bound is 1047 for the production headers and 62888 at the generic declared arity limits. Even multiplying generic multiplicities by `INT64_MAX` residual caps, higher aggregation fits signed 128 at `580039420653723140950616` (79 bits) before the actual verifier's checked narrowing. Pure potential-weight comparison plus the largest allowed recharge stays within `1110873980937437184` signed 64-bit ticks. These are targeted arithmetic guards, not an assertion that every arbitrary allowed input fits the checked 64-bit residual output.

```sh
for side in lower upper; do
  python3 -B "$C/coverage.py" "$R/quintet_$side.bin" "$R/gold.json"
  python3 -B "$C/overflow.py" "$R/quintet_$side.bin"
done
```

All four commands exited 0. They read just the headers and the existing receipt, not the giant tables. No production replay or board edit was made in this lane.
