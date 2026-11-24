#!/usr/bin/env bash
set -euo pipefail
root="$(cd -- "$(dirname -- "$0")/../../.." && pwd)"
lane="$root/round2/proofs/neon/round2"
python=${PYTHON:-python3}
case "${1:-}" in
    ''|--compile|--slow) ;;
    *) echo "Usage: $0 [--compile|--slow]" >&2; exit 2 ;;
esac
if (($# > 1)); then
    echo "Usage: $0 [--compile|--slow]" >&2
    exit 2
fi
if [[ ${1:-} == --slow ]]; then
    env -u PYTHONOPTIMIZE bash "$lane/cert_l5/test_cert_l5.sh" 10 12
    echo 'RESULT neon-l5-d10: D10 exact GMP, tree, sanitizer and negative-control suites passed; benchmark only'
    echo 'RESULT neon-l5-d12: D12 exact GMP, tree and inward-override checks passed; benchmark only'
elif [[ ${1:-} == --compile ]]; then
    scratch=$(mktemp -d)
    trap 'rm -rf -- "$scratch"' EXIT
    for kind in cert cert_l5; do
        stem=cert_b3f
        [[ $kind == cert ]] || stem=cert_b3f_l5
        gcc -std=gnu11 -O2 -Wall -Wextra -Werror -fopenmp "$lane/$kind/${stem}_build.c" -lm -o "$scratch/${stem}_build"
        gcc -std=gnu11 -O2 -Wall -Wextra -Werror -fopenmp "$lane/$kind/${stem}_verify_gmp.c" -lgmp -o "$scratch/${stem}_verify_gmp"
    done
    echo 'RESULT neon-b3f-compile: both builders and exact GMP verifiers compiled; no witness checked'
else
    "$python" -IB "$lane/open_questions/exact_degrees.py" | diff -u "$lane/open_questions/exact_degrees.txt" -
    echo 'RESULT neon-exact-degrees: Claim 2 exact degree and leaf transcript matches'
    "$python" -IB "$lane/open_questions/check_disagreement.py" | diff -u "$lane/open_questions/disagreement_audit.txt" -
    echo 'RESULT neon-disagreement: exact contraction and local rate transcript matches'
    "$python" -IB "$lane/checks/b3f18_sol_algebra.py"
    echo 'RESULT neon-b3f-algebra: exact Bernstein and endpoint algebra checks passed; no witness checked'
fi
