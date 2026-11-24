#!/usr/bin/env python3
"""Check compiled occupation scores by individual physical rings on small trees."""
import argparse
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cache
from itertools import combinations
import os
from pathlib import Path
import subprocess
import sys
import tempfile

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "motif-high"))
from oracle import Certificate

DISCOUNT = F(1133, 2000)
LEAF = (2, 1)


@cache
def geometry(state):
    parents, leaves = state
    near = [[] for _ in parents]
    for y, x in enumerate(parents[1:], 1):
        near[x].append(y)
        near[y].append(x)
    cd = tuple(map(len, near))
    degree = tuple(a + b for a, b in zip(leaves, cd))
    assert all(d >= 2 for d in degree)
    return near, degree, cd


def transitions(state):
    """Sum physical initiator rings, including every send involving the founder."""
    parents, leaves = state
    near, degree, _ = geometry(state)
    rates = defaultdict(F)
    for v, neighbours in enumerate(near):
        # Every adjacent leaf ringing at v produces a new leaf at v.
        rates[v, False] += leaves[v]
        for w in neighbours:
            rates[w, False] += F(1, degree[v])
        rates[v, True] += F(leaves[v], degree[v])
    for (recipient, promotion), rate in rates.items():
        if not rate:
            continue
        new = list(leaves)
        if promotion:
            new[recipient] -= 1
            yield (parents + (recipient,), tuple(new) + (1,)), rate, promotion
        else:
            new[recipient] += 1
            yield (parents, tuple(new)), rate, promotion


@cache
def subsets(parents, size):
    result = []
    for members in combinations(range(1, len(parents)), size):
        if all(parents[j] in members for j in members[1:]):
            shape = "".join(str(members.index(parents[j])) for j in members[1:]) or "-"
            result.append((members, shape))
    return result


def quartet(parents, members):
    root = members[0]
    kids = {v: [w for w in members if parents[w] == v] for v in members}
    a = kids[root]
    if len(a) == 3:
        return "T", (root, *a)
    if len(a) == 2:
        mid = next(x for x in a if kids[x])
        other = next(x for x in a if x != mid)
        return "F", (root, other, mid, kids[mid][0])
    mid = a[0]
    if len(kids[mid]) == 2:
        return "Y", (root, mid, *kids[mid])
    v = kids[mid][0]
    return "P", (root, mid, v, kids[v][0])


def read_dump(filename):
    rows = {}
    with filename.open() as f:
        for line in f:
            shape, *numbers = line.split()
            *flat, score, source = map(int, numbers)
            key = shape, tuple(flat)
            assert key not in rows and 0 <= score <= 2**40 and source >= 0
            rows[key] = score, source
    return rows


def physical_test(cert, rows, caps, purpose, max_size):
    rmax = 5 if purpose == "quartet" else 2
    v = rows["-", LEAF][0]
    counts = Counter()
    mutations = set()

    @cache
    def score(state):
        parents, _ = state
        _, degree, core = geometry(state)
        total = 0
        selected = {}
        for r in range(1, rmax+1):
            for members, shape in subsets(parents, r):
                flat = []
                for j, x in enumerate(members):
                    lim, clim = caps[r-1][j]
                    flat.extend((min(lim, clim)+1, 0) if degree[x] > lim or core[x] > clim
                                else (degree[x], core[x]))
                z, f = rows.get((shape, tuple(flat)), (0, 0))
                total += z
                selected[members] = (z, f)
        return total, selected

    @cache
    def delta(family, states, at):
        old = cert.weight(family, states)
        nxt = list(states)
        l, k = nxt[at]
        nxt[at] = l+1, k
        return cert.sign * (cert.weight(family, tuple(nxt))-old)

    def check(state):
        parents, leaves = state
        near, degree, core = geometry(state)
        t0, positions = score(state)
        actual_source = defaultdict(F)
        false_source = defaultdict(F)
        families = set()
        if purpose == "quartet":
            lim, clim = caps[-1][0]
            for members, _ in subsets(parents, 4):
                if any(degree[x] > lim or core[x] > clim for x in members):
                    continue
                family, vertices = quartet(parents, members)
                states = tuple((leaves[x], core[x]-1) for x in vertices)
                for i, x in enumerate(vertices):
                    d = delta(family, states, i)
                    for z in near[x]:
                        if z == 0 or z in vertices:
                            continue
                        ident = tuple(sorted((*vertices, z)))
                        reciprocal = F(1, degree[z])
                        actual_source[ident] += max(0, d)/F(2)-reciprocal*d
                        false_source[ident] += max(0, -d)/F(2)+reciprocal*d
                        if d:
                            families.add(family)
        else:
            start, _ = cert.sections["f"]
            for members, _ in subsets(parents, 2):
                x, z = members if purpose == "parent" else members[::-1]
                d, k = degree[x], core[x]-1
                if not cert.par["D2"] <= d < cert.par["DF"]-1:
                    continue
                change = cert.sign*(cert.data[start+d*(d+1)//2-1+k]
                                    -cert.data[start+d*(d-1)//2-1+k])
                reciprocal = F(1, degree[z])
                actual_source[members] = max(0, change)/F(2)-reciprocal*change
                false_source[members] = max(0, -change)/F(2)+reciprocal*change
        for members, (_, src) in positions.items():
            assert actual_source[members] >= src, (state, members, src, actual_source[members])
            if src > false_source[members]:
                mutations.add("reverse donor incidence")
        for f in families:
            counts[f] += 1
        total_credit = sum(actual_source.values(), F())
        gt, b = F(), F()
        for nxt, rate, promotion in transitions(state):
            t1, after = score(nxt)
            gt += rate*(t1-t0)
            counts["transitions"] += 1
            if promotion:
                b += rate
                fresh = len(parents)
                old_change = sum(after[group][0]-z for group, (z, _) in positions.items())
                extensions = 0
                for group in positions:
                    if len(group) < rmax and nxt[0][-1] in group:
                        extensions += after[group+(fresh,)][0]
                assert t1-t0 == old_change+extensions+v
                if t1-t0 != old_change+v:
                    counts["positive extension"] += 1
                    mutations.add("drop new connected sets")
        requirement = DISCOUNT*t0-gt+v*b
        assert total_credit >= requirement, (state, total_credit, requirement)
        # This mutation actually violates the assembled test.
        if total_credit < requirement + 2**40*b:
            mutations.add("oversized newborn credit")
        counts["states"] += 1
        counts["positive physical sources"] += bool(total_credit)

    layer = {((-1,), (2,))}
    for n in range(3, max_size+1):
        next_layer = set()
        for state in layer:
            check(state)
            if n < max_size:
                next_layer.update(nxt for nxt, _, _ in transitions(state))
        layer = next_layer
    # Decorations of all four shapes ensure independent exercise of the
    # parent, child, sibling and forward donor directions, including pooling.
    for parents in [(-1,0,1,2,3,4), (-1,0,1,1,2,3), (-1,0,1,2,2,3),
                    (-1,0,1,1,1,2), (-1,0,1,1,2,2,3,3)]:
        cd = [parents.count(i)+(i!=0) for i in range(len(parents))]
        for donor in range(1, len(parents)):
            for d_high in (3, 5, 72, 74, 85, 97, 1401):
                leaves = [max(1 if c == 1 else 0, 3-c) for c in cd]
                leaves[donor] = max(0, d_high-cd[donor])
                check((parents, tuple(leaves)))
    for high in (71, 72, 73, 78, 84, 85, 86, 1400):
        for neighbour in (2, 3, 4, 5, 6, 8, 11, 32, 193):
            check(((-1,0,1), (2,high-2,neighbour-1)))
            check(((-1,0,1), (2,neighbour-2,high-1)))
    assert {"reverse donor incidence", "oversized newborn credit", "drop new connected sets"} <= mutations, (purpose, counts, mutations)
    if purpose == "quartet":
        assert all(counts[f] for f in "PFYT"), counts
    print(purpose, "states", counts["states"], "transitions", counts["transitions"],
          "positive physical sources", counts["positive physical sources"],
          "positive newborn extensions", counts["positive extension"],
          "planted errors", len(mutations))
    score.cache_clear()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("binary", type=Path)
    parser.add_argument("snapshot", type=Path)
    parser.add_argument("--size", type=int, default=10)
    args = parser.parse_args()
    cert = Certificate(args.snapshot)
    with tempfile.TemporaryDirectory() as tmp:
        filename = Path(tmp, "rows")
        for purpose in ("quartet", "parent", "child"):
            if purpose == "quartet":
                d, c = (8,6,5,4,4), (6,4,4,4,4)
                caps = [[(d[r-1],c[r-1])]*r for r in range(1,6)]
                tail = []
            else:
                d, c = (96,85,6,4,4), (8,7,5,4,4)
                pair = [(d[1],c[1]),(d[2],c[2])]
                caps = [[(d[0],c[0])], pair if purpose == "parent" else pair[::-1]]
                tail = [purpose]
            command = [str(args.binary.resolve()), str(args.snapshot),
                       ",".join(map(str,d)), ",".join(map(str,c)), "1133", "2000", *tail]
            subprocess.run(command, env={**os.environ, "FT2_RESERVE_DUMP":str(filename)},
                           check=True, stdout=subprocess.DEVNULL)
            physical_test(cert, read_dump(filename), caps, purpose, args.size)


if __name__ == "__main__":
    main()
