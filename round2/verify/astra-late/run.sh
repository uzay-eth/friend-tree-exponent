#!/usr/bin/env bash
set -euo pipefail
cd -- "$(dirname -- "$0")/../.."
verify=$PWD/verify/astra-late
bundle=$PWD/proofs/astra-late
cd "$bundle/source"
python=${PYTHON:-python3}
"$python" -IB -c 'import sys; sys.exit(sys.flags.optimize != 0 or sys.version_info < (3, 10))'
independent() {
  for audit in audit_claim17 audit_claims_19_20 audit_claims_22_23 audit_claim26 audit_claim27; do
    printf '\n== %s ==\n' "$audit"
    "$python" -IB "$bundle/independent/$audit.py"
  done
  printf '%s\n' 'RESULT astra-late-gemini: independent finite audits for Claims 17, 19, 20, 22, 23, 26 and 27'
}
case "${1:-all}" in
all|quick)
  for check in \
    checks/limit_law_sol.py checks/left_tail_sol.py checks/limit_density_sol.py \
    checks/left_prefactor_sol.py open/limit-law/right/check.py \
    checks/age_order_sol.py checks/low_forest_sol.py open/forest-age/independent.py \
    checks/forest-marked-repair-sol.py checks/marked-spectrum-sol.py \
    open/age-order/spectrum/quenched/check.py open/age-order/moments/check.py \
    open/forest-age/fixed-shapes/diagonal/check.py open/forest-age/macroscopic/check.py \
    checks/eventual-domination-sol.py open/age-order/hub-ranking/check.py; do
    printf '\n== %s ==\n' "$check"
    "$python" -IB "round2/$check"
  done
  "$python" -IB "$verify/with-originals.py" round2/checks/effective-uniform-sol.py
  "$python" -IB "$verify/with-originals.py" round2/open/age-order/hub-ranking/all-powers/positive/check.py
  printf '%s\n' 'RESULT astra-late-fixtures: finite amplitude, forest, label-order and moment controls'
  "$python" -IB "$verify/d18-archived.py"
  printf '%s\n' 'RESULT astra-late-d18: archived high101 endpoint arithmetic and finite registration controls'
  "$python" -IB "$verify/gold2-archived.py"
  "$python" -IB "$bundle/independent/astra28-lifetime_recompute.py" "$PWD"
  printf '%s\n' 'RESULT astra-late-gold2: archived-input rational exclusion through 425 and frontier receipt hashes'
  "$python" -IB -c 'from runpy import run_path; a=run_path("round2/checks/limit-cusp-sol.py"); a["budget"](); a["kernel_control"](); print("PASS scoped cusp budgets and model fixtures; last width", a["critical_model"]())'
  printf '%s\n' 'RESULT astra-late-cusp: model scaling and translation fixtures, excluding historical source pins'
  if [[ ${1:-all} == all ]]; then independent; fi
  ;;
independent)
  independent
  ;;
threshold)
  "$python" -IB round2/open/age-order/hub-ranking/poisson-tail/check.py
  "$python" -IB round2/open/age-order/global-transport/quarter/b20.py
  "$python" -IB round2/open/age-order/global-transport/quarter/check.py
  ;;
*) printf '%s\n' 'usage: run.sh [all|quick|independent|threshold]' >&2; exit 2 ;;
esac
