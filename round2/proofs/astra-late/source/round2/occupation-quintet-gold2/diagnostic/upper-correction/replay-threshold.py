#!/usr/bin/env python3
"""Replay the frozen 79-cap H3 mask with correctly scaled reward slopes."""
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
GOLD = HERE.parents[1]
SOURCE = GOLD/'post-fy/replay.py'
PIN = 'd8a9a8a65d6b652151c41de18fa7a0f4e8aec4a6ed8e89395412e4ddd5bc3c23'
OLD = b'shift=self.round(w/width) if width else 0'
NEW = b'shift=self.round(d*w/width) if width else 0'


def load(*, corrected=False):
    raw = SOURCE.read_bytes()
    assert sha256(raw).hexdigest() == PIN and raw.count(OLD) == 1
    if corrected:
        raw = raw.replace(OLD, NEW)
        # Remove the frozen expectation of the old, incorrectly scaled value.
        old = b"if a.kind=='threshold' and a.side=='lower':\n                    assert total==Q(1084423030952180046412413909,1704710415481443450880000000000000)"
        assert raw.count(old) == 1
        raw = raw.replace(old, b'if True:\n                    print("corrected exact ticks/core", total*s.scale, flush=True)')
        print('threshold source', PIN, 'executed', sha256(raw).hexdigest(), flush=True)
    spec = spec_from_file_location('gold2_threshold', SOURCE)
    m = module_from_spec(spec)
    exec(compile(raw, str(SOURCE), 'exec'), m.__dict__)
    return m


if __name__ == '__main__':
    assert __debug__ and len(sys.argv) == 2 and sys.argv[1] in ('lower', 'upper')
    sys.argv[1:] = ['--kind', 'threshold', '--side', sys.argv[1]]
    load(corrected=True).main()
