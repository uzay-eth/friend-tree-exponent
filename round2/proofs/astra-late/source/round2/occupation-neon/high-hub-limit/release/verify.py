"""Verify the high101 ledger and optionally re-run one independent solve."""

from contextlib import redirect_stdout
from io import StringIO
from pathlib import Path
import runpy
import sys


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[3]
COMPANION = ROOT / 'round2/checks/high-hub-sol.py'
EXTERNAL = ROOT / 'round2/occupation-neon/interpolate/edge/higher/transport/replay.py'
EXTERNAL_SHA = 'a92a5e830ac6ed092eeb2be1c73709244f66aef9df1a227cf969785467a7e870'


def main():
    args = sys.argv[1:]
    if len(args) not in (2, 4) or (len(args) == 4 and tuple(args[2:]) not in
                                   (('--solve', 'lo'), ('--solve', 'hi'))):
        raise SystemExit('usage: verify.py LOWER UPPER [--solve {lo,hi}]')
    if not __debug__:
        raise SystemExit('run with assertions enabled')
    cold = runpy.run_path(str(COMPANION))
    assert cold['PINS'].pop(EXTERNAL) == EXTERNAL_SHA
    if EXTERNAL.exists():
        assert cold['digest'](EXTERNAL) == EXTERNAL_SHA
    sys.argv[1:] = [*args[:2], str(HERE), *args[2:]]
    text = StringIO()
    with redirect_stdout(text):
        cold['main']()
    result = text.getvalue()
    reference = HERE / ('cold-hi.log' if args[2:] == ['--solve', 'hi'] else 'cold-lo.log')
    expected = reference.read_text()
    if len(args) == 2:
        expected = ''.join(expected.splitlines(keepends=True)[:3])
    assert result == expected, result
    print(result, end='')
    print('PASS high101: original D18 payloads and exact aa511af + high101 ledger')


if __name__ == '__main__':
    main()
