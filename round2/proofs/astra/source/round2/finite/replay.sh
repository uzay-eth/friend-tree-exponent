#!/usr/bin/env bash
# Rebuild and independently verify the six finite rational matrices (Python 3.12).
set -euo pipefail
cd -- "$(dirname -- "$0")/../.."
replay_matrix=$(mktemp)
trap 'rm -f -- "$replay_matrix"' EXIT
python3 -B round2/finite/check_blocks.py
python3 -B round2/finite/finite_matrix.py --case 3,1 --case 3,3 --case 3,9 --case 4,1 --case 4,2 --case 4,3 --out "$replay_matrix"
cmp -- "$replay_matrix" round2/outputs/finite-matrices.json
python3 -B round2/finite/verify_matrix.py "$replay_matrix"
python3 -B round2/checks/one_port_neon.py
python3 -B round2/checks/finite_operator_neon.py
