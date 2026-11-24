#!/usr/bin/env bash
set -euo pipefail

# One-command reproducible build, pure-GMP verification, independent Python
# exact-Fraction tree-drift check, and negative-control suite for Neon's
# B3-FullH + Endpoint-Only Coarse Level-5 (version=2) certified enclosure of mu.

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
L5_DIR="$ROOT/round2/cert_l5"
OUT_DIR="${B3F_L5_OUT_DIR:-${TMPDIR:-.}/b3f_l5_run}"
THREADS="${B3F_L5_THREADS:-4}"
PROFILE="${B3F_L5_PROFILE:-20}"
mkdir -p "$OUT_DIR"

echo "=== 1. Compiling version=2 builder and pure-GMP verifier ==="
gcc -std=c11 -O3 -march=native -ffp-contract=off -fopenmp -Wall -Wextra -Werror \
  "$L5_DIR/cert_b3f_l5_build.c" -lm -o "$OUT_DIR/cert_b3f_l5_build"
gcc -std=c11 -O3 -march=native -fopenmp -Wall -Wextra -Werror \
  "$L5_DIR/cert_b3f_l5_verify_gmp.c" -lgmp -o "$OUT_DIR/cert_b3f_l5_verify_gmp"

case "$PROFILE" in
  10)
    ARGS=(--D0 10 --KT 9 --LF 4000 --KF 80 --DP4 8 --KP4 6 --DYT 7 --KYT 5 --DYB 8 --KYB 6 --DS3 6 --KS3 5 --Dc 6 --buckets 2 --l5mask 3 --Q 10000000 --cushion 5e-7 --nit 22)
    ;;
  12)
    ARGS=(--D0 12 --KT 11 --LF 6000 --KF 100 --DP4 10 --KP4 8 --DYT 9 --KYT 7 --DYB 10 --KYB 8 --DS3 8 --KS3 6 --Dc 7 --buckets 2 --l5mask 3 --Q 10000000 --cushion 5e-7 --nit 18)
    ;;
  20)
    ARGS=(--D0 20 --KT 19 --LF 16000 --KF 200 --DP4 12 --KP4 9 --DYT 11 --KYT 8 --DYB 12 --KYB 9 --DS3 10 --KS3 8 --Dc 7 --buckets 2 --l5mask 3 --Q 10000000 --cushion 6e-7 --nit 22)
    ;;
  24)
    ARGS=(--D0 24 --KT 23 --LF 18000 --KF 220 --DP4 12 --KP4 10 --DYT 11 --KYT 9 --DYB 12 --KYB 10 --DS3 10 --KS3 8 --Dc 7 --buckets 2 --l5mask 3 --Q 10000000 --cushion 6e-7 --nit 22)
    ;;
  *)
    echo "Unsupported B3F_L5_PROFILE=$PROFILE (expected 10, 12, 20, or 24)" >&2
    exit 2
    ;;
esac

echo "=== 2. Building & Pure-GMP Verifying B3-FullH-$PROFILE + Coarse-L5 (Lower & Upper) ==="
"$OUT_DIR/cert_b3f_l5_build" --mode 2 "${ARGS[@]}" --threads "$THREADS" --out "$OUT_DIR/b3f${PROFILE}_l5_lo.bin"
sha256sum "$OUT_DIR/b3f${PROFILE}_l5_lo.bin"
"$OUT_DIR/cert_b3f_l5_verify_gmp" "$OUT_DIR/b3f${PROFILE}_l5_lo.bin" --threads "$THREADS" --tree-check 30
python3 -B "$L5_DIR/tree_identity_l5_check.py" "$OUT_DIR/b3f${PROFILE}_l5_lo.bin" 30

"$OUT_DIR/cert_b3f_l5_build" --mode 1 "${ARGS[@]}" --threads "$THREADS" --out "$OUT_DIR/b3f${PROFILE}_l5_hi.bin"
sha256sum "$OUT_DIR/b3f${PROFILE}_l5_hi.bin"
"$OUT_DIR/cert_b3f_l5_verify_gmp" "$OUT_DIR/b3f${PROFILE}_l5_hi.bin" --threads "$THREADS" --tree-check 30
python3 -B "$L5_DIR/tree_identity_l5_check.py" "$OUT_DIR/b3f${PROFILE}_l5_hi.bin" 30

echo "=== 3. Impossible-Exponent Override Controls ==="
check_reject() {
  local log="$1"; shift
  local rc=0
  "$@" >"$log" 2>&1 || rc=$?
  cat "$log"
  if [ "$rc" -ne 1 ] || ! grep -q '=== VERDICT: FAIL ' "$log"; then
    echo "ERROR: expected clean mathematical rejection (rc=1, VERDICT: FAIL), got rc=$rc" >&2
    exit 1
  fi
}
check_reject "$OUT_DIR/ov_lo.log" "$OUT_DIR/cert_b3f_l5_verify_gmp" "$OUT_DIR/b3f${PROFILE}_l5_lo.bin" --threads "$THREADS" --tree-check 0 --override-pq 566000 1000000
check_reject "$OUT_DIR/ov_hi.log" "$OUT_DIR/cert_b3f_l5_verify_gmp" "$OUT_DIR/b3f${PROFILE}_l5_hi.bin" --threads "$THREADS" --tree-check 0 --override-pq 565000 1000000

echo "=== ALL B3-FullH-$PROFILE + COARSE-L5 CHECKS PASSED ==="
