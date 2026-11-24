#!/usr/bin/env python3
"""Print the exact rational outputs of the frozen high101 bulk producer."""
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from bulk import HI, LO, Rewards, arms, run, stage_bounds, ztuple  # noqa: E402


def main(blob):
    states, arcs = arms(7)
    reward = Rewards(Path(blob), states)
    polys = {x: reward.polynomials(x) for x in states}
    inf = {x: p[0] if p else ztuple() for x, p in polys.items()}
    pres = {x: stage_bounds(p, inf[x], 101) for x, p in polys.items()}
    print(f'witness {reward.c.sign} M=101 T=101 K=12 cap=7 LO={LO} HI={HI}', flush=True)
    for upper in (False, True):
        result = run(states, arcs, polys, inf, pres, 101, 12, 101, upper)
        for role, values in result.items():
            print('ceiling' if upper else 'floor', role, *(str(x) for x in values), str(sum(values)), flush=True)


if __name__ == '__main__':
    main(sys.argv[1])
