#!/usr/bin/env python3
"""A nonzero size-seven boundary continuation, via the next exact shape shell."""
import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bulk import Rewards, arms, downward, run, stage_bounds, upward, v, ztuple  # noqa: E402


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('blob', type=Path)
    cli.add_argument('--cap', type=int, default=8, choices=(8, 9, 10))
    cli.add_argument('--side', choices=('both', 'floor', 'ceiling'), default='both')
    a = cli.parse_args()
    states, arcs = arms(a.cap)
    reward = Rewards(a.blob, states)
    # At size seven, stop paying current reward but keep its attached shape
    # and all cap-respecting successors; only the new reward is valued here.
    polys = {x: reward.polynomials(x) if v.size(x)>7 else None for x in states}
    inf = {x: p[0] if p else ztuple() for x, p in polys.items()}
    pres = {x: stage_bounds(p, inf[x], 101) for x, p in polys.items()}
    print('witness',reward.c.sign,'D>=101 k<=12 shell', (8,a.cap),'states',len(states),flush=True)
    for upper in (False, True):
        if a.side != 'both' and upper != (a.side == 'ceiling'):
            continue
        fmt = upward if upper else downward
        for role, result in run(states, arcs, polys, inf, pres, 101, 12, 101, upper).items():
            print('ceiling' if upper else 'floor',role,*(fmt(x) for x in result),fmt(sum(result)),
                  'exact',*(str(x) for x in result),flush=True)


if __name__ == '__main__':
    main()
