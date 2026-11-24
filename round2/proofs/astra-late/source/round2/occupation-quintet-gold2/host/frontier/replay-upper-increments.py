#!/usr/bin/env python3
"""Replay the two disjoint post-FY Gold2 upper increments with corrected units."""
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

GOLD = Path(__file__).resolve().parents[2]
SOURCE = GOLD/'post-fy/upper-frontier/frontier.py'
PIN = '1e2071a42f1b39c52c647f6a3b0f0205b4f12309d107aa8874e4b22afd7a9592'
BASE_PIN = 'd8a9a8a65d6b652151c41de18fa7a0f4e8aec4a6ed8e89395412e4ddd5bc3c23'
OLD = b'self.round(w/width)'
NEW = b'self.round(d*w/width)'


def main():
    assert __debug__
    raw = SOURCE.read_bytes()
    assert sha256(raw).hexdigest() == PIN and raw.count(OLD) == 1
    assert sha256((SOURCE.parents[1]/'replay.py').read_bytes()).hexdigest() == BASE_PIN
    data = raw.replace(OLD, NEW)
    # Full exterior receive R changes the payoff by -slope*R. A gain
    # numerator width of d*R therefore needs a per-rate shift -slope.
    for d in (2, 6, 70):
        R, slope = Q(3, 2), Q(7, 3)
        w, width = -slope*R, d*R
        assert d*w/width == -slope and w/width > -slope
    m = module_from_spec(spec_from_file_location('gold2_corrected_increments', SOURCE))
    exec(compile(data, str(SOURCE), 'exec'), m.__dict__)
    print('corrected increment source', PIN, 'executed', sha256(data).hexdigest(), flush=True)
    original = sys.argv
    try:
        for name, args in (
            ('middle3', '--receiver 120 --rstock 10 --donor 7 --middle-min 3 --fixtures'),
            ('stock', '--old 19 --receiver 79 --rstock 12 --mother-min-stock 9 --donor 7'),
        ):
            print('corrected upper increment', name, flush=True)
            sys.argv = [str(SOURCE), *args.split()]
            m.main()
    finally:
        sys.argv = original
    print('PASS corrected original Gold2 upper H3 increments')


if __name__ == '__main__':
    main()
