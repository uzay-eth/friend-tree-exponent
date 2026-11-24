#!/usr/bin/env python3
"""Cross-check sparse rewards and compressed hazards on complete physical hosts."""
from collections import defaultdict
from pathlib import Path
import sys

from check import FAMILIES, J, OCC, Q, Rewards, arms, plain, shape, v

sys.path.insert(0, str(OCC / 'higher' / 'older'))
from replay import Local  # noqa: E402


def tree(leaves, *children):
    return (2, leaves, tuple(sorted(children)))


def fixtures(blob):
    states, arcs = arms(7)
    states = [x for x in states if x[2]]
    reward = Rewards(blob, states)
    independent = Local(reward.c.c)
    counts = defaultdict(int)
    for d, k in ((32, 1), (32, 17), (49, 5), (65, 5), (100, 10), (101, 18), (128, 34)):
        for mother_degree, pool in ((2, (J,)), (11, (tree(5), J, tree(0, J)))):
            other = tuple(pool[j % len(pool)] for j in range(k-1))
            for x in states:
                # Keep the tracked arm at physical index 2 even for equal siblings.
                host = (2, d-k-1, (x, *other))
                whole = tree(mother_degree-2, host)
                p, tags = shape(whole, 1)
                ch = [[] for _ in tags]
                for i, j in enumerate(p[1:], 1):
                    ch[j].append(i)
                expected = dict.fromkeys(FAMILIES, Q(0))
                for f, ids in independent.occurrences(ch, tags):
                    if f not in expected or ids[:2] != (1, 2):
                        continue
                    ext = [(Q(1, sum(tags[p[i]])) if p[i] not in ids else Q(0),
                            sum(Q(1, sum(tags[j])) for j in ch[i] if j not in ids))
                           for i in ids]
                    value = independent.charge(f, tuple(tags[j] for j in ids), ext)
                    assert value >= 0
                    expected[f] += value
                    counts[f] += bool(value)
                u = Q(1, mother_degree) + sum(Q(1, v.degree(y)) for y in other)
                a, b = reward.row(x, d, k)
                assert tuple(a + b*u for a, b in zip(a, b)) == tuple(expected.values()), (d, k, x, u)

                tagged = (3, x[1], x[2])
                full = tree(mother_degree-2, tree(d-k-1, tagged, *other))
                projected = defaultdict(Q)
                for y, r in v.row(full)[3].items():
                    h, = [child for child in y[2] if v.degree(child) >= d]
                    child, = [child for child in h[2] if child[0] == 3]
                    state = (v.degree(h), len(h[2]), plain(child))
                    if state != (d, k, x):
                        projected[state] += r
                v.row.cache_clear()
                want = defaultdict(Q, {(d+1, k, x): d-k-1+u+Q(1, v.degree(x)),
                                       (d, k+1, x): Q(d-k-1, d)})
                want[d, k, v.plus(x)] += Q(1, d)
                for y, r in arcs[x].items():
                    want[d, k, y] += r
                # The whole tree also emits intrinsic transitions out of the arm cap.
                assert {y: r for y, r in projected.items() if v.size(y[2]) <= 7} == {
                    y: r for y, r in want.items() if r and v.size(y[2]) <= 7}, (d, k, x)
                total = d-k-1+u+Q(d-k, d)+v.size(x)
                assert sum(projected.values()) == total, (d, k, x)
                counts['hosts'] += 1
    assert all(counts[f] for f in FAMILIES)
    print(f'witness {reward.c.sign}: {counts["hosts"]} exact complete-host reward and projected-arrow matches; '
          f'nonzero occurrences {dict((f, counts[f]) for f in FAMILIES)}', flush=True)


if __name__ == '__main__':
    for blob in map(Path, sys.argv[1:]):
        fixtures(blob)
