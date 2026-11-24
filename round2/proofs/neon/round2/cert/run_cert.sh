#!/usr/bin/env bash
set -euo pipefail

# One-command reproducible build, pure-GMP verification, UBSan check,
# independent Python exact-Fraction tree-drift check, and negative-control suite
# for Neon's B3-FullH + Level-4 certified enclosure of mu.

ROOT="$(cd "$(dirname "$0")/../.." && pwd)"
CERT_DIR="$ROOT/round2/cert"
OUT_DIR="${B3F_OUT_DIR:-${TMPDIR:-.}/b3f_cert_run}"
THREADS="${B3F_THREADS:-4}"
PROFILE="${B3F_PROFILE:-18}"
mkdir -p "$OUT_DIR"

echo "=== 1. Compiling builder, pure-GMP verifier, and UBSan verifier ==="
gcc -O3 -march=native -fopenmp -Wall -Wextra "$CERT_DIR/cert_b3f_build.c" -lm -o "$OUT_DIR/cert_b3f_build"
gcc -O3 -march=native -fopenmp -Wall -Wextra "$CERT_DIR/cert_b3f_verify_gmp.c" -lgmp -o "$OUT_DIR/cert_b3f_verify_gmp"
gcc -O2 -fsanitize=undefined -fno-sanitize-recover=all -fopenmp -Wall -Wextra "$CERT_DIR/cert_b3f_verify_gmp.c" -lgmp -o "$OUT_DIR/cert_b3f_verify_gmp_ubsan"

echo "=== 2. Fast D0=10 UBSan + Python exact-Fraction + Negative Controls ==="
"$OUT_DIR/cert_b3f_build" --mode 2 --D0 10 --KT 9 --LF 4000 --KF 80 \
  --DP4 6 --KP4 4 --DYT 6 --KYT 4 --DYB 6 --KYB 4 --DS3 5 --KS3 3 \
  --nit 24 --threads "$THREADS" --out "$OUT_DIR/d10_lo.bin"
"$OUT_DIR/cert_b3f_verify_gmp_ubsan" "$OUT_DIR/d10_lo.bin" --threads "$THREADS" --tree-check 15
python3 -B "$CERT_DIR/tree_identity_check.py" "$OUT_DIR/d10_lo.bin" 10
python3 -B "$CERT_DIR/negative_controls.py" "$OUT_DIR/cert_b3f_verify_gmp_ubsan" "$OUT_DIR/d10_lo.bin"

"$OUT_DIR/cert_b3f_build" --mode 1 --D0 10 --KT 9 --LF 4000 --KF 80 \
  --DP4 6 --KP4 4 --DYT 6 --KYT 4 --DYB 6 --KYB 4 --DS3 5 --KS3 3 \
  --nit 24 --threads "$THREADS" --out "$OUT_DIR/d10_hi.bin"
"$OUT_DIR/cert_b3f_verify_gmp_ubsan" "$OUT_DIR/d10_hi.bin" --threads "$THREADS" --tree-check 15
python3 -B "$CERT_DIR/tree_identity_check.py" "$OUT_DIR/d10_hi.bin" 10
python3 -B "$CERT_DIR/negative_controls.py" "$OUT_DIR/cert_b3f_verify_gmp_ubsan" "$OUT_DIR/d10_hi.bin"

if [ "$PROFILE" = "20" ]; then
  echo "=== 3. Production B3-FullH-20-L4 Lower & Upper Certificates (Pure GMP) ==="
  "$OUT_DIR/cert_b3f_build" --mode 2 --D0 20 --KT 19 --LF 16000 --KF 200 \
    --DP4 12 --KP4 9 --DYT 11 --KYT 8 --DYB 12 --KYB 9 --DS3 10 --KS3 7 \
    --nit 24 --threads "$THREADS" --out "$OUT_DIR/b3f20_lo.bin"
  sha256sum "$OUT_DIR/b3f20_lo.bin"
  "$OUT_DIR/cert_b3f_verify_gmp" "$OUT_DIR/b3f20_lo.bin" --threads "$THREADS" --tree-check 30
  python3 -B "$CERT_DIR/tree_identity_check.py" "$OUT_DIR/b3f20_lo.bin" 10

  "$OUT_DIR/cert_b3f_build" --mode 1 --D0 20 --KT 19 --LF 16000 --KF 200 \
    --DP4 12 --KP4 9 --DYT 11 --KYT 8 --DYB 12 --KYB 9 --DS3 10 --KS3 7 \
    --nit 24 --threads "$THREADS" --out "$OUT_DIR/b3f20_hi.bin"
  sha256sum "$OUT_DIR/b3f20_hi.bin"
  "$OUT_DIR/cert_b3f_verify_gmp" "$OUT_DIR/b3f20_hi.bin" --threads "$THREADS" --tree-check 30
  python3 -B "$CERT_DIR/tree_identity_check.py" "$OUT_DIR/b3f20_hi.bin" 10

  echo "=== 4. Production Impossible-Exponent Negative Controls ==="
  if "$OUT_DIR/cert_b3f_verify_gmp" "$OUT_DIR/b3f20_lo.bin" --threads "$THREADS" --tree-check 0 --override-pq 566000 1000000; then
    echo "ERROR: lower cert unexpectedly passed impossible exponent 566000/1000000!" >&2; exit 1
  fi
  if "$OUT_DIR/cert_b3f_verify_gmp" "$OUT_DIR/b3f20_hi.bin" --threads "$THREADS" --tree-check 0 --override-pq 565000 1000000; then
    echo "ERROR: upper cert unexpectedly passed impossible exponent 565000/1000000!" >&2; exit 1
  fi
  echo "=== ALL B3-FullH-20-L4 CHECKS PASSED ==="
  exit 0
fi

echo "=== 3. Production B3-FullH-18 Lower & Upper Certificates (Pure GMP) ==="
"$OUT_DIR/cert_b3f_build" --mode 2 --D0 18 --KT 17 --LF 14000 --KF 180 \
  --DP4 10 --KP4 7 --DYT 9 --KYT 6 --DYB 10 --KYB 7 --DS3 8 --KS3 6 \
  --nit 24 --threads "$THREADS" --out "$OUT_DIR/b3f18_lo.bin"
sha256sum "$OUT_DIR/b3f18_lo.bin"
"$OUT_DIR/cert_b3f_verify_gmp" "$OUT_DIR/b3f18_lo.bin" --threads "$THREADS" --tree-check 30
python3 -B "$CERT_DIR/tree_identity_check.py" "$OUT_DIR/b3f18_lo.bin" 15

"$OUT_DIR/cert_b3f_build" --mode 1 --D0 18 --KT 17 --LF 14000 --KF 180 \
  --DP4 10 --KP4 7 --DYT 9 --KYT 6 --DYB 10 --KYB 7 --DS3 8 --KS3 6 \
  --nit 24 --threads "$THREADS" --out "$OUT_DIR/b3f18_hi.bin"
sha256sum "$OUT_DIR/b3f18_hi.bin"
"$OUT_DIR/cert_b3f_verify_gmp" "$OUT_DIR/b3f18_hi.bin" --threads "$THREADS" --tree-check 30
python3 -B "$CERT_DIR/tree_identity_check.py" "$OUT_DIR/b3f18_hi.bin" 15

echo "=== 4. Production Impossible-Exponent Negative Controls ==="
if "$OUT_DIR/cert_b3f_verify_gmp" "$OUT_DIR/b3f18_lo.bin" --threads "$THREADS" --tree-check 0 --override-pq 566000 1000000; then
  echo "ERROR: lower cert unexpectedly passed impossible exponent 566000/1000000!" >&2; exit 1
fi
if "$OUT_DIR/cert_b3f_verify_gmp" "$OUT_DIR/b3f18_hi.bin" --threads "$THREADS" --tree-check 0 --override-pq 565000 1000000; then
  echo "ERROR: upper cert unexpectedly passed impossible exponent 565000/1000000!" >&2; exit 1
fi

echo "=== ALL B3-FullH-18 CHECKS PASSED ==="
