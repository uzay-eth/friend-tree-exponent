#!/usr/bin/env bash
set -euo pipefail
here=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
cd -- "$here/../../proofs/astra/source"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1
unset PYTHONOPTIMIZE PYTHONPATH
python=${PYTHON:-python3}
export PATH
# The original finite replay invokes python3 directly. Use the chosen
# interpreter's directory for that script too.
python=$(command -v -- "$python")
PATH="$(dirname -- "$python"):$PATH"
"$python" -B -c 'import sys; sys.exit(0 if __debug__ and sys.version_info >= (3, 11) else "Python 3.11+ with assertions enabled required")'
check() {
    "$python" -B "round2/checks/$1.py"
}
result() {
    printf 'RESULT astra-%s: %s\n' "$1" "$2"
}

case ${1:-modest} in
    modest)
        work=$(mktemp -d -t astra-modest.XXXXXXXX)
        trap 'rm -rf -- "$work"' EXIT
        bash round2/finite/replay.sh                     # 1
        result 1 'six rational finite matrices and one-port controls'
        check fringe_neon                               # 2
        result 2 'finite marked-fringe generator and weight controls'
        check asymptotics_carbon                        # 3
        check standalone_3_sol
        result 3 'mean transform, fixed-size transfer and degree-law fixtures'
        "$python" -B - <<'PY'                         # 4
from runpy import run_path
f = run_path('round2/checks/operator_carbon.py')
for name in ('test_3_generator_dynkin_and_degree_law_identities',
             'test_4_both_edge_transfer_cancellations', 'test_5_planted_errors'):
    f[name]()
PY
        result 4 'finite generator, edge cancellations and negative controls'
        check finite_kernel_neon                        # 5
        result 5 'small rational cuts and path sums for capped kernel'
        check strong_limit_carbon                       # 6
        check standalone_6_sol
        result 6 'strong-limit and toy contour controls'
        check hub_mass_carbon                           # 7
        check standalone_7_sol
        result 7 'finite hub mass and stopping controls'
        check local_completeness_neon                    # 8
        result 8 'finite phase-averaged drift and grammar controls'
        check spatial_limit_hy4                         # 10
        check fringe_marginals_hy4
        result 10 'local-limit and marginal finite fixtures'
        check hub_spectrum_hy4                          # 11
        check hub_children_sol
        result 11 'rank, intermediate-degree and hub-child fixtures'
        (
            cd round2/finite/gauge
            "$python" -B excursion_check.py
            excursion() {
                printf '$ excursion_rows.py %s\n' "$*"
                "$python" -B excursion_rows.py "$@"
            }
            excursion --size 9 --s .5665 --test A
            excursion --size 9 --s .5665 --test B
            excursion --size 9 --s .54001 --test B --virtual
            excursion --size 9 --s .565615 --test half
            excursion --size 9 --s .565809 --test half
            excursion --size 7 --s 1 --test scalar
        ) > "$work/excursion.log"
        "$python" -B round2/checks/excursions_hy4.py "$work/excursion.log"
        check scalar_hy4
        check physical_trigger_hy4
        check response_transcendence_hy4
        result 12 'scalar-return fixtures and six regenerated excursion cases'
        check uncapped_gauge_sol                        # 13
        result 13 'current-parent uncapped rational weight controls'
        check zero_contraction_sol                      # 14
        "$python" -B round2/influence/zero/critical_controls.py
        g++ -O3 -std=c++20 round2/influence/zero/pair/pair.cpp -o "$work/pair"
        "$work/pair" 48 > "$work/pair.rows" 2> "$work/pair.log"
        grep -Ex 'slice DA=[0-9]+ states=[0-9]+ cumulative=22956300' "$work/pair.log"
        grep -Fx 'max initial finite 8600511260 at 3 0; tail 8484503829' "$work/pair.log"
        result 14 'signed scalar controls and 22,956,300-row native cap-48 scan'
        "$python" -B - <<'PY'                           # 9 and 15 (frozen arithmetic)
from fractions import Fraction
from pathlib import Path
from runpy import run_path

root = Path('round2')
a = run_path(str(root / 'joint-enclosure/arithmetic.py'))
read = lambda p: (root / p).read_text()
rows = a['read_ledgers'](read('outputs/occupation-affine-verified.txt'),
                         read('five-next/split-reference.tsv'), a['CEILING'])
old = tuple(a['LEGACY'][s] for s in a['SIDES'])
lo, hi = old
assert Fraction('0.5656581227') <= lo < Fraction('0.5656581228')
assert Fraction('0.5657979940') < hi <= Fraction('0.5657979941')
assert Fraction('0.565615') < lo <= hi < Fraction('0.565809')
ends = a['enclosure'](rows)
ticks = a['ticks'](ends)
a['crosscheck'](rows, ends, ticks)
a['plants'](rows, ends, ticks)
print('Astra 9 frozen ledger:', *old)
print('Astra 15 frozen ledger:', *ends, 'outward ticks:', ticks)
PY
        result 9 'original physical enclosure recomposed from the frozen ledgers'
        result 15 'joint enclosure recomposed from the frozen ledgers'
        check generation_sol                           # 16
        result 16 'physical depth-transfer fixtures'
        ;;
    physical)
        "$python" -B - <<'PY'
import os
from pathlib import Path
from runpy import run_path
from tempfile import TemporaryDirectory
from types import SimpleNamespace

if hasattr(os, 'sched_getaffinity'):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
with TemporaryDirectory(prefix='astra-physical-') as work:
    run_path('round2/occupation-tail/replay.py')['run'](
        SimpleNamespace(snapshots=None), Path(work))
PY
        ;;
    physical-audits)
        "$python" -B - <<'PY'
import os
from pathlib import Path
from runpy import run_path
import subprocess
import sys
from tempfile import TemporaryDirectory

if hasattr(os, 'sched_getaffinity'):
    os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
with TemporaryDirectory(prefix='astra-physical-audits-') as work:
    work = Path(work)
    subprocess.run([sys.executable, '-B', 'round2/motif-cert/replay.py',
                    '--work', work / 'base'], check=True)
    run_path('round2/checks/replay-physical-reviews.py')['audit'](work / 'base', work / 'audit')
PY
        ;;
    spatial-full)
        check operator_carbon
        ;;
    joint-preflight|joint-cold)
        mode=${1#joint-}
        "$python" -B round2/joint-enclosure/replay.py "$mode"
        ;;
    *)
        printf '%s\n' 'usage: run.sh [modest|physical|physical-audits|spatial-full|joint-preflight|joint-cold]' >&2
        exit 2
        ;;
esac
