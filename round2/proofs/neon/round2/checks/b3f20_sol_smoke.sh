#!/usr/bin/env bash
set -euo pipefail
export OMP_NUM_THREADS=2 UBSAN_OPTIONS=halt_on_error=1 TMPDIR="$PWD"
for mode in ${*:-2 1}; do
  cert=./b3f20_sol_d10_m${mode}.bin
  ./b3f20_sol_build_ubsan --mode "$mode" --D0 10 --KT 9 --LF 4000 --KF 80 \
    --DP4 6 --KP4 4 --DYT 6 --KYT 4 --DYB 6 --KYB 4 --DS3 5 --KS3 3 \
    --nit 24 --cushion 0.000002 --threads 2 --out "$cert"
  ./b3f20_sol_verify_ubsan "$cert" --threads 2 --tree-check 15
  python3 -B - "$cert" <<'PY'
# Use the exact shipped suite, replacing its hard-coded four-thread switch.
from pathlib import Path
import sys
p = Path('round2/cert/negative_controls.py')
src = p.read_text()
assert src.count('"--threads", "4"') == 1
sys.argv = [str(p), './b3f20_sol_verify_ubsan', sys.argv[1]]
src = src.replace('"--threads", "4"', '"--threads", "2"')
src = src.replace('    return r.returncode, r.stdout',
                  '    assert "runtime error:" not in r.stdout, r.stdout\n'
                  '    return r.returncode, r.stdout')
exec(compile(src, str(p), 'exec'), {'__name__': '__main__'})
PY
done
