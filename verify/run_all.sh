#!/usr/bin/env bash
# Build and run every finite exact-arithmetic verification in this repository.
# Prints detailed progress and ends with one PASS line per result.
set -euo pipefail

VERIFY_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd -- "$VERIFY_DIR/.." && pwd)"
export PYTHONDONTWRITEBYTECODE=1

CC="${CC:-gcc}"
CXX="${CXX:-$(command -v g++-13 || command -v g++)}"
START_TS=$(date +%s)

TMP="$(mktemp -d)"
trap 'rm -rf -- "$TMP"' EXIT

PASS_LINES=()
record_pass() {
  PASS_LINES+=("$1")
  printf '\n>>> %s\n\n' "$1"
}

echo "=== Environment ==="
echo "CC:      $("$CC" --version | head -n 1)"
echo "CXX:     $("$CXX" --version | head -n 1)"
echo "Python:  $(python3 --version)"
echo "gzip:    $(gzip --version | head -n 1)"
if command -v dpkg-query >/dev/null 2>&1; then
  GMP_VER="$(dpkg-query -W -f='${Version}' libgmp-dev 2>/dev/null || echo 'installed (non-dpkg)')"
  echo "libgmp:  $GMP_VER"
fi
echo ""

# Bootstrap local Python virtual environment with pinned requirements
VENV="$VERIFY_DIR/.venv"
if [[ ! -x "$VENV/bin/python" ]]; then
  echo "=== Creating Python venv ($VENV) and installing verify/requirements.txt ==="
  python3 -m venv "$VENV"
fi
"$VENV/bin/pip" install --quiet --disable-pip-version-check -r "$VERIFY_DIR/requirements.txt"
PY="$VENV/bin/python"

# ---------------------------------------------------------------------------
# 1. Motif bracket [113/200, 1133/2000] (Gemini 4 Neon)
# ---------------------------------------------------------------------------
echo "=== [1/7] Motif bracket [113/200, 1133/2000] ==="
T0=$(date +%s)
"$CC" -O3 -Wall -Wextra -fopenmp "$VERIFY_DIR/motif-bracket/cert_build.c" -lm -o "$TMP/cert_build"
"$CC" -O3 -Wall -Wextra -fopenmp "$VERIFY_DIR/motif-bracket/cert_verify_gmp.c" -lgmp -o "$TMP/cert_verify_gmp"

"$TMP/cert_build" 1 1133 2000 "$TMP/upper_5665.bin"
"$TMP/cert_verify_gmp" "$TMP/upper_5665.bin"

"$TMP/cert_build" 2 113 200 "$TMP/lower_5650.bin"
"$TMP/cert_verify_gmp" "$TMP/lower_5650.bin"

"$PY" "$VERIFY_DIR/motif-bracket/verify_python_cross.py" \
  --upper "$TMP/upper_5665.bin" \
  --lower "$TMP/lower_5650.bin" \
  --exact "$VERIFY_DIR/motif-bracket/exact_rationals_n21.txt"

for args in \
  "$TMP/upper_5665.bin --override-pq 1132 2000 --skip-trees" \
  "$TMP/lower_5650.bin --override-pq 114 200 --skip-trees" \
  "$TMP/upper_5665.bin --perturb f --skip-trees" \
  "$TMP/upper_5665.bin --perturb g --skip-trees" \
  "$TMP/upper_5665.bin --perturb hV --skip-trees" \
  "$TMP/upper_5665.bin --perturb hC --skip-trees"; do
  set +e
  # shellcheck disable=SC2086
  "$TMP/cert_verify_gmp" $args >/dev/null 2>&1
  rc=$?
  set -e
  if [[ "$rc" -ne 2 ]]; then
    echo "ERROR: motif-bracket negative control unexpectedly returned $rc: $args" >&2
    exit 1
  fi
done
echo "Negative controls: 6/6 rejected (rc=2)"
T1=$(date +%s)
record_pass "PASS  motif-bracket (Neon)          113/200 <= E[C_n] <= 1133/2000 (< 17/30 by 1/6000); 17,329,108 GMP checks/side, 6/6 neg controls ($((T1-T0))s)"

# ---------------------------------------------------------------------------
# 2. Banded upper bound 5657/9983 — primary verifier (GPT-6.2 Astra)
# ---------------------------------------------------------------------------
echo "=== [2/7] Banded upper bound 5657/9983 (primary verifier) ==="
T0=$(date +%s)
"$CXX" -std=c++17 -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all \
  "$VERIFY_DIR/banded-upper/check.cpp" -o "$TMP/check_banded"
gzip -dc "$VERIFY_DIR/banded-upper/recipe.txt.gz" > "$TMP/recipe.txt"
"$TMP/check_banded" "$VERIFY_DIR/banded-upper/geometry.txt" "$TMP/recipe.txt" 5657 9983
T1=$(date +%s)
record_pass "PASS  banded-upper (Astra)          E[C_n] <= 1539738114735/3908512529 prod(1+5657/(9983j)) (< 17/30 by 1/299490); 226,545,617 rows ($((T1-T0))s)"

# ---------------------------------------------------------------------------
# 3. Banded upper bound 5657/9983 — independent verifier (Gemini 4 Neon)
# ---------------------------------------------------------------------------
echo "=== [3/7] Banded upper bound 5657/9983 (independent verifier) ==="
T0=$(date +%s)
"$CXX" -std=c++17 -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all \
  "$VERIFY_DIR/banded-upper/independent/verify_banded.cpp" -o "$TMP/verify_banded"
"$TMP/verify_banded" 3 5 \
  "$VERIFY_DIR/banded-upper/independent/fixtures/bands-smoke.txt" \
  "$VERIFY_DIR/banded-upper/independent/fixtures/smoke-recipe.txt"
"$PY" "$VERIFY_DIR/banded-upper/independent/test_banded_audit.py" \
  --exe "$TMP/verify_banded" \
  --bands "$VERIFY_DIR/banded-upper/independent/fixtures/bands-smoke.txt" \
  --recipe "$VERIFY_DIR/banded-upper/independent/fixtures/smoke-recipe.txt"
"$TMP/verify_banded" 5657 9983 "$VERIFY_DIR/banded-upper/geometry.txt" "$TMP/recipe.txt"
T1=$(date +%s)
record_pass "PASS  banded-upper-audit (indep.)   independent verifier 5657/9983; 226,545,617 rows, 436 histories, 50 adversarial trees, 14/14 neg controls ($((T1-T0))s)"

# ---------------------------------------------------------------------------
# 4. Transfer [563/1000, 568/1000] and hub [56/100, 57/100] (GPT-6.2 Astra)
# ---------------------------------------------------------------------------
echo "=== [4/7] Transfer [563/1000, 568/1000] and hub [56/100, 57/100] ==="
T0=$(date +%s)
"$CXX" -std=c++17 -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all \
  "$VERIFY_DIR/transfer-bounds/transfer.cpp" -o "$TMP/transfer"
"$CXX" -std=c++17 -O2 -Wall -Wextra -fsanitize=undefined -fno-sanitize-recover=all \
  "$VERIFY_DIR/transfer-bounds/hub.cpp" -o "$TMP/hub"
"$TMP/transfer" lower 563 1000 208 20 128 26 10 6 5 4 4 3 3 3
"$TMP/transfer" upper 568 1000 208 20 128 26 10 6 5 4 4 3 3 3
"$TMP/hub" lower 56 100 128 12 64 18 8 6 5 4 4 3 3 3
"$TMP/hub" upper 57 100 168 12 90 20 9 6 5 4 4 3 3 3
T1=$(date +%s)
record_pass "PASS  transfer-bounds (Astra)       [563/1000, 568/1000] (9,639,742 states) and hub [56/100, 57/100] ($((T1-T0))s)"

# ---------------------------------------------------------------------------
# 5. Almost-sure lower 543/1000 & baseline lower table (Claude Fable 5.2)
# ---------------------------------------------------------------------------
echo "=== [5/7] Almost-sure lower 543/1000 (baseline lower table + H1-H5) ==="
T0=$(date +%s)
"$PY" "$VERIFY_DIR/cross-charging/as_lower_check.py"
if PYTHONPATH="$VERIFY_DIR/baseline" "$PY" "$VERIFY_DIR/baseline/certs9.py" 0.5433 0.5434 0 11/20 60 10 2 6 >/dev/null 2>&1; then
  echo "ERROR: baseline lower negative control (11/20) unexpectedly passed" >&2
  exit 1
fi
echo "Negative control (certs9 lower at 11/20): rejected as expected"
T1=$(date +%s)
record_pass "PASS  almost-sure-lower (Fable)     liminf log C_n/log n >= 543/1000 a.s.; 1,205,767 exact checks, (H1)-(H5), 1/1 neg control ($((T1-T0))s)"

# ---------------------------------------------------------------------------
# 6. Cross-charging [11/20, 29/50] + 11/20 a.s. (Claude Fable 5.2)
# ---------------------------------------------------------------------------
echo "=== [6/7] Cross-charging [11/20, 29/50] + 11/20 a.s. ==="
T0=$(date +%s)
"$PY" "$VERIFY_DIR/cross-charging/xcharge.py" test 200 11 | tail -n 1
"$PY" "$VERIFY_DIR/cross-charging/xdump.py" topkl "$VERIFY_DIR/cross-charging/certs/upper_4_3_60_12.txt" "$TMP/xc_upper.pkl"
"$PY" "$VERIFY_DIR/cross-charging/xdump.py" topkl "$VERIFY_DIR/cross-charging/certs/lower_4_3_60_12.txt" "$TMP/xc_lower.pkl"
"$PY" "$VERIFY_DIR/cross-charging/xverify.py" "$TMP/xc_upper.pkl" 29/50
"$PY" "$VERIFY_DIR/cross-charging/xverify.py" "$TMP/xc_lower.pkl" 11/20
"$PY" "$VERIFY_DIR/cross-charging/xverify.py" "$TMP/xc_upper.pkl" 579/1000 --expect-fail
"$PY" "$VERIFY_DIR/cross-charging/xverify.py" "$TMP/xc_lower.pkl" 552/1000 --expect-fail
T1=$(date +%s)
record_pass "PASS  cross-charging (Fable)        [11/20, 29/50] + 11/20 a.s.; 5,853/5,851 exact checks, 200 trees, 2/2 neg controls ($((T1-T0))s)"

# ---------------------------------------------------------------------------
# 7. Baseline upper 293/500 & polynomial-scale lower 543/1000 (Carbon / Sol)
# ---------------------------------------------------------------------------
echo "=== [7/7] Baseline upper 293/500 & polynomial-scale lower 543/1000 ==="
T0=$(date +%s)
PYTHONPATH="$VERIFY_DIR/baseline" "$PY" "$VERIFY_DIR/baseline/find_upper_lin.py" 60 10 1 20 0.58 0.59 0 293/500
if PYTHONPATH="$VERIFY_DIR/baseline" "$PY" "$VERIFY_DIR/baseline/find_upper_lin.py" 60 10 1 20 0.58 0.59 0 29/50 0.586 >/dev/null 2>&1; then
  echo "ERROR: baseline upper negative control (29/50) unexpectedly passed" >&2
  exit 1
fi
echo "Negative control (find_upper_lin at 29/50): rejected as expected"
"$PY" "$VERIFY_DIR/polynomial-lower/verify_as_lower.py" --report-tables
"$PY" "$VERIFY_DIR/polynomial-lower/check_as_lower.py"
T1=$(date +%s)
record_pass "PASS  polynomial-lower (Carbon/Sol) liminf C_n/n^(543/1000) > 0 a.s.; 319,086 upper checks, 102 trees, 84 branches, 15/15 neg controls ($((T1-T0))s)"

END_TS=$(date +%s)
TOTAL_S=$((END_TS - START_TS))

echo "==========================================================================="
echo "SUMMARY OF ALL FINITE VERIFICATIONS (total wall time: ${TOTAL_S}s)"
echo "==========================================================================="
for line in "${PASS_LINES[@]}"; do
  echo "$line"
done
echo "==========================================================================="
echo "ALL 7 VERIFICATION SUITES PASSED (${TOTAL_S}s)"
