#!/usr/bin/env bash
set -euo pipefail

# Portable verification runner for the Round 3 Section 7 public checkpoint.
# Requires: g++ (C++17, __int128 support), Python 3.10+ (standard library only).
# Run from the repository root: bash round3/verify_checkpoint.sh

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
R3="$ROOT/round3"
TMPDIR_WORK="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_WORK"' EXIT

echo "=== [1/7] Verify byte-identical Q3 C++ SHA-256 digests ==="
sha256sum -c <<EOF
90689e44d08fdd94bb9edc76cae6d14669a5ac39e84d4da35d19a82d8a517c80  $R3/q3-birth-order/cert_all_m.cpp
a10dce484242acb727973baf48c1bdac448dce9e761c506f55326a66b58c3210  $R3/q3-birth-order/cert_small_m.cpp
EOF

echo "=== [2/7] Q3: Build & run exact __int128 certificates + primitive/defect audits ==="
g++ -O3 -std=c++17 -Wall -Wextra "$R3/q3-birth-order/cert_all_m.cpp" -o "$TMPDIR_WORK/cert_all_m"
g++ -O3 -std=c++17 -Wall -Wextra "$R3/q3-birth-order/cert_small_m.cpp" -o "$TMPDIR_WORK/cert_small_m"
"$TMPDIR_WORK/cert_all_m"
"$TMPDIR_WORK/cert_small_m"
python3 -IB "$R3/q3-birth-order/audit_primitives.py"
python3 -IB "$R3/q3-birth-order/planted_defects.py"

echo "=== [3/7] Q2: Run exact rational urn-tail counterexample/repair & kernel audit ==="
python3 -IB "$R3/q2-growing-degrees/audit_exact.py"

echo "=== [4/7] Q5: Run worst single-parent KS finite/variational checker ==="
python3 -IB "$R3/q5-low-degree-forest/worst_hub_upper_check.py"

echo "=== [5/7] Q6: Run two-attachment six-rule & repaired-constant exact checker ==="
python3 -IB "$R3/q6-multiple-attachments/check_delta2.py"

echo "=== [6/7] Q7: Run beta(p)=0 on [927/1000, 1) outward rational certificate ==="
(cd "$R3/q7-mixed-attachments" && python3 -B check_elasticity_match.py)

echo "=== [7/7] Q8: Run pinned exact-Fraction k=9 verifier ==="
python3 -IB "$R3/q8-longer-walks/verify_k9.py"

echo "=== ALL ROUND 3 PACKAGED CHECKS PASSED ==="
