#!/usr/bin/env bash
# Compact, fail-closed replay of the light checks behind the paper.
# Run from anywhere: bash verify/compact_replay.sh
# Needs: Python 3.10+ (standard library only), g++ with C++17 and __int128, sha256sum.
# It does not replay the large drift certificates; see verify/README.md and verify/refined-lower/README.md.
set -euo pipefail
ROOT="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/.." && pwd)"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT
cd "$ROOT"
die(){ echo "COMPACT REPLAY FAILED: $*" >&2; exit 1; }

echo "=== [1/5] Theorem PREC: exact finite checks (verify/mean-arithmetic) ==="
python3 -I verify/mean-arithmetic/check_window.py 17 1200 | tee "$WORK/window.out"
[ "$(tail -n1 "$WORK/window.out")" = "RESULT: PASS" ] || die "check_window.py"

echo "=== [2/5] Theorem ENC lower endpoint: receipt arithmetic and planted defects (verify/refined-lower) ==="
(cd verify/refined-lower && sha256sum --quiet -c MANIFEST.sha256) || die "refined-lower manifest"
python3 -I verify/refined-lower/check_ratio.py --self-test | tee "$WORK/ratio.out"
[ "$(tail -n1 "$WORK/ratio.out")" = "RESULT: PASS" ] || die "check_ratio.py"

echo "=== [3/5] Theorem ORD: certificate digests ==="
Q=round3/q3-birth-order
sha256sum --quiet -c <<EOF || die "stochastic-order certificate digests"
90689e44d08fdd94bb9edc76cae6d14669a5ac39e84d4da35d19a82d8a517c80  $Q/cert_all_m.cpp
a10dce484242acb727973baf48c1bdac448dce9e761c506f55326a66b58c3210  $Q/cert_small_m.cpp
EOF

echo "=== [4/5] Theorem ORD: exact integer certificates, primitive audit, planted defects ==="
g++ -O3 -std=c++17 -Wall -Wextra "$Q/cert_all_m.cpp" -o "$WORK/cert_all_m"
g++ -O3 -std=c++17 -Wall -Wextra "$Q/cert_small_m.cpp" -o "$WORK/cert_small_m"
"$WORK/cert_all_m" | tee "$WORK/all.out"
grep -q "Boxes PASSED: 52, FAILED: 0" "$WORK/all.out" || die "cert_all_m"
"$WORK/cert_small_m" | tee "$WORK/small.out"
grep -q "Boxes PASSED: 84, FAILED: 0" "$WORK/small.out" || die "cert_small_m"
python3 -IB "$Q/audit_primitives.py"
python3 -IB "$Q/planted_defects.py"

echo "=== [5/5] Printed enclosure matches the receipts ==="
python3 -B paper/checks/check_literals.py --enclosure | tee "$WORK/enclosure.out"
tail -n1 "$WORK/enclosure.out" | grep -q '^PASS:' || die "check_literals.py --enclosure"

echo "=== ALL COMPACT CHECKS PASSED ==="
