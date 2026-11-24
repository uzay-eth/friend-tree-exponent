"""All-age original D18 rewards of a band of oriented four-core paths."""

import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import file_digest
from math import ceil, floor, lcm
from pathlib import Path
import sys

from path import HIGH, LOW

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "higher"))
from extract import Cert, Extract, PAYLOAD_SHA


GRID = 10**18
VIRGIN = (2, 0)
OVERFLOW = (0, 0)
ZEROS = (0, 0)


def original(path, side):
    with path.open("rb") as stream:
        digest = file_digest(stream, "sha256").hexdigest()
    if digest != PAYLOAD_SHA[side]:
        raise ValueError(f"{side} D18 payload SHA-256 mismatch")
    return Extract(Cert(path))


def physical_row(table, path):
    states = tuple((d - k - 1, k, 1) for d, k in path)
    if not table.supported("P4", states):
        return None
    base, free = table.row("P4", states)
    coefficients = {(i, parent): Q(value, table.c.S) for i, parent, _, value in free}
    if coefficients.keys() - {(0, True), (1, False)}:
        raise AssertionError("unexpected exterior port on the declared mask")
    return base / table.c.S, coefficients.get((0, True), 0), coefficients.get((1, False), 0)


def clean_fringe(path):
    """Sufficient physical support for a P4 reward disjoint from the frozen payments."""
    side_children = path[1][1] - 1
    min_size = sum(d - k for d, k in path) + 2 * side_children
    return min_size > 20 or (min_size > 12 and side_children > 0)


def solve(tables, d_caps, k_middle, donor_cap, mother_cap, theta, single):
    a_caps = (1, k_middle, 1, 0)
    # The six charges share one state graph, with independent outward bounds per objective.
    names = ("mass", "lower", "upper", "contrast", "unpaid-lower", "unpaid-upper")
    zeros = (ZEROS,) * len(names)
    unit = lcm(1_000_000, *range(2, max(*d_caps, donor_cap + 1, mother_cap + 1) + 1))
    reciprocal = [0] + [unit // d for d in range(1, max(*d_caps, donor_cap + 1, mother_cap + 1) + 1)]
    death = int(unit * LOW)
    death_range = int(unit * (HIGH - LOW))
    side_tail = reciprocal[max(2, donor_cap + 1)]
    mother_tail = reciprocal[max(2, mother_cap + 1)]
    seen = set()
    unpaid = set()
    states_visited = 0
    new_side = VIRGIN if donor_cap else OVERFLOW

    @cache
    def charge(path):
        if len(path) < 4 or (single and path != ((2, 1), (4, 3), (2, 1), (4, 0))):
            return zeros, zeros, zeros
        lo, hi = (physical_row(t, path) for t in tables)
        if lo is None or hi is None:
            assert lo is None and hi is None
            return zeros, zeros, zeros
        seen.add(path)
        clean = clean_fringe(path)
        if clean:
            unpaid.add(path)
        diff = tuple(u - theta * v for u, v in zip(lo, hi))
        empty = (0, 0, 0)
        scores = ((1, 0, 0), lo, hi, diff, lo if clean else empty, hi if clean else empty)
        return tuple(tuple((floor(v * GRID), ceil(v * GRID)) for v in src) for src in zip(*scores))

    def bellman(toll, edges):
        total = death + sum(a for a, _, _, _ in edges)
        out = []
        for j in range(len(names)):
            pair = []
            for upper in (0, 1):
                den = total
                num = toll[j][upper] * unit
                choices = [(0, death_range)]
                for low_rate, high_rate, successor, slope in edges:
                    value = successor[j][upper] + slope[j][upper]
                    num += low_rate * value
                    if high_rate > low_rate:
                        choices.append((value, high_rate - low_rate))
                for value, width in sorted(choices, reverse=bool(upper)):
                    if (value * den <= num) if upper else (value * den >= num):
                        break
                    den += width
                    num += width * value
                bound = -((-num) // den) if upper else num // den
                drift = toll[j][upper] * unit - (death + (death_range if ((bound < 0) == bool(upper)) else 0)) * bound
                for low_rate, high_rate, successor, slope in edges:
                    delta = successor[j][upper] + slope[j][upper] - bound
                    drift += (high_rate if ((delta > 0) == bool(upper)) else low_rate) * delta
                if (drift > 0 if upper else drift < 0):
                    raise AssertionError("discounted drift check failed")
                pair.append(bound)
            out.append(tuple(pair))
        return tuple(out)

    @cache
    def future(path, sides, mother):
        nonlocal states_visited
        n = len(path)
        if any(d > cap or k > kmax for (d, k), cap, kmax in zip(path, d_caps, a_caps)):
            return zeros
        # A short prefix registers future children of its current last vertex.
        # Its already-registered extended paths continue on their own.
        if n < 4 and path[-1][1] == a_caps[n - 1]:
            return zeros
        if n >= 2:
            assert len(sides) == path[1][1] - int(n >= 3)
        states_visited += 1
        intercept, parent_slope, child_slope = charge(path)
        edges = []

        def edge(a, b, value, slope=zeros):
            assert 0 <= a <= b
            if b:
                edges.append((a, b, value, slope))

        for i, (d, k) in enumerate(path):
            leaves = d - 1 - k
            gain = future(path[:i] + ((d + 1, k),) + path[i+1:], sides, mother)
            fixed = leaves * unit
            if i:
                fixed += reciprocal[path[i-1][0]]
            if i + 1 < n:
                fixed += reciprocal[path[i+1][0]]
            edge(fixed, fixed, gain)
            if i == 0:
                edge(reciprocal[mother], reciprocal[mother] if mother else mother_tail, gain, parent_slope)
            if i == 1 and sides:
                known = sum(reciprocal[dd] for dd, _ in sides)
                edge(known, known + sides.count(OVERFLOW) * side_tail, gain, child_slope)
            if leaves:
                promoted = path[:i] + ((d, k + 1),) + path[i+1:]
                plain_sides = tuple(sorted(sides + (new_side,))) if i == 1 else sides
                after = future(promoted, plain_sides, mother)
                if i == n - 1 and n < 4:
                    extended = future(promoted + (VIRGIN,), sides, mother)
                    after = tuple((a+c, b+e) for (a,b),(c,e) in zip(after, extended))
                edge(leaves * reciprocal[d], leaves * reciprocal[d], after)
        for (d, k), count in Counter(sides).items():
            if not d:
                continue
            j = sides.index((d, k))
            def replace(new):
                return future(path, tuple(sorted(sides[:j] + (new,) + sides[j+1:])), mother)
            leaves = d - 1 - k
            gain = leaves * unit + reciprocal[path[1][0]]
            edge(count*gain, count*(gain + k*reciprocal[2]),
                 replace((d+1, k) if d < donor_cap else OVERFLOW))
            if leaves:
                edge(count*leaves*reciprocal[d], count*leaves*reciprocal[d], replace((d, k+1)))
        if mother:
            own_child = reciprocal[path[0][0]]
            edge(own_child, own_child + mother*unit - 3*reciprocal[2],
                 future(path, sides, mother+1 if mother < mother_cap else 0))
        return bellman(intercept, edges)

    entries = {m: future((VIRGIN,), (), m) for m in (0, *range(2, mother_cap + 1))}
    bounds = {}
    for i, name in enumerate(names):
        vals = [s*Q(v[i][endpoint], GRID) for v in entries.values()
                for endpoint in (0, 1) for s in (LOW, HIGH)]
        bounds[name] = min(vals), max(vals)
    for label in ("mass", "lower", "upper", "unpaid-lower", "unpaid-upper"):
        assert bounds[label][1] >= 0
    print("physical D18 rows", len(seen), "unpaid rows", len(unpaid), "live states", states_visited,
          "checked signed drift inequalities", 2*len(names)*states_visited, flush=True)
    for name, (lo, hi) in bounds.items():
        print(name, lo, hi, f"[{float(lo):.9g},{float(hi):.9g}]")
    future.cache_clear()
    return bounds


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lower", type=Path, required=True)
    parser.add_argument("--upper", type=Path, required=True)
    parser.add_argument("--degrees", type=int, nargs=4, default=(4,5,3,5), metavar=("A","B","C","D"))
    parser.add_argument("--middle-arity", type=int, default=3)
    parser.add_argument("--donor", type=int, default=2)
    parser.add_argument("--mother", type=int, default=3)
    parser.add_argument("--theta", type=Q, default=Q(1))
    parser.add_argument("--single", action="store_true")
    args = parser.parse_args()
    if not all(2 <= d <= 16 for d in args.degrees) or not 1 <= args.middle_arity < args.degrees[1]:
        parser.error("tagged degrees must be from 2 to 16; invalid middle arity")
    if any(v == 1 or v < 0 for v in (args.donor, args.mother)):
        parser.error("memory caps must be zero or at least two")
    ts = original(args.lower, "lower"), original(args.upper, "upper")
    print("original payloads authenticated; mask", args.degrees, args.middle_arity,
          "donor/mother", args.donor, args.mother, "theta", args.theta, flush=True)
    solve(ts, args.degrees, args.middle_arity, args.donor, args.mother, args.theta, args.single)


if __name__ == "__main__":
    main()
