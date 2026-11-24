#!/usr/bin/env bash
set -euo pipefail

# One-command reproduction script for the Hartle-Krapivsky-Redner-Zhang (arXiv:2604.01540)
# Fig. 6 fitting-protocol experiment on IR, DAN, and PAN growing trees.
# Uses at most 2 CPU cores at a time.

DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR"

PYTHON="../.venv/bin/python"          # the round-2 environment (round2/requirements.txt)
if [ ! -x "$PYTHON" ]; then
    PYTHON="python3"
fi

echo "=== Step 1: Exact analytical roots and small-N tree enumeration ==="
"$PYTHON" validate_exact.py

echo "=== Step 2: Compile fast C++20 simulator ==="
g++ -O3 -march=native -std=c++20 -pthread sim.cpp -o sim

echo "=== Step 3: Validate fast simulator against exact enumeration (N=4..8) and naive simulator (N=100) ==="
./sim validate 10000000 200000

echo "=== Step 4: Run 80 protocol replicates (2400 trees/rep) - COUPLED (same trees across all grid sizes) ==="
./sim run --model IR  --mode coupled --replicates 80 --trees 2400 --threads 2 --seed 1001 --out ir_coupled
./sim run --model DAN --mode coupled --replicates 80 --trees 2400 --threads 2 --seed 2001 --out dan_coupled
./sim run --model PAN --mode coupled --replicates 80 --trees 2400 --threads 2 --seed 3001 --out pan_coupled

echo "=== Step 5: Run 60 protocol replicates (2400 trees/rep) - INDEPENDENT (fresh trees per grid size) ==="
./sim run --model IR  --mode indep   --replicates 60 --trees 2400 --threads 2 --seed 4001 --out ir_indep
./sim run --model DAN --mode indep   --replicates 60 --trees 2400 --threads 2 --seed 5001 --out dan_indep
./sim run --model PAN --mode indep   --replicates 60 --trees 2400 --threads 2 --seed 6001 --out pan_indep

echo "=== Step 6: Extract exact Fig. 6 data from arXiv:2604.01540 PDF vector stream ==="
# The extracted stream is checked in (cvsn_stream.txt).  It is regenerated only if the PDF of the paper has been
# placed at ../lit/raw/2604.01540.pdf (the PDF itself is not part of the repository) and mutool is installed.
if [ -f ../lit/raw/2604.01540.pdf ] && command -v mutool >/dev/null; then
    mutool show -b ../lit/raw/2604.01540.pdf 553 > cvsn_stream.tmp && mv cvsn_stream.tmp cvsn_stream.txt
fi
"$PYTHON" parse_fig6.py

echo "=== Step 7: Analyze all replicates, fit variants, local slopes, and CV(C) ==="
"$PYTHON" analyze.py "$DIR" | tee analysis_output.txt

echo "=== Step 8: Verify independent-mode SD(b) against exact delta-method error propagation ==="
"$PYTHON" verify_sampling_theory.py
