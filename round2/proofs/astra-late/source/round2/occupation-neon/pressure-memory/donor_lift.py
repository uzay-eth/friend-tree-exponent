"""Price each exterior donor of an oriented D18 P4 starting from its actual birth."""

import argparse
from fractions import Fraction as Q
from functools import cache
from math import ceil, floor, lcm

from path_band import GRID, HIGH, LOW, OVERFLOW, VIRGIN, original
from pathlib import Path


NO_VALUE = (0, 0, 0, 0)
NO_DONOR = (-1, -1)
OLD_DEGREES = (12, 8, 6, 16)


def eligible(path, exclude_previous):
    if len(path) != 4:
        return False
    outside = sum(k for d, k in path) - 3
    min_size = sum(d - k for d, k in path) + 2 * outside
    if not (min_size > 20 or (min_size > 12 and outside > 0)):
        return False
    return not (exclude_previous and all(d <= limit for (d, k), limit in zip(path, OLD_DEGREES))
                and path[0][1] == 1 and path[1][1] <= 3 and path[2][1] == 1 and path[3][1] == 0)


def interval(v):
    return floor(v * GRID), ceil(v * GRID)


class Reward:
    def __init__(self, tables, exclude_previous):
        self.tables = tables
        self.exclude_previous = exclude_previous
        self.row = cache(self.row)
        self.toll = cache(self.toll)

    def row(self, path):
        """Return the pointwise minimum and the five physical pressure coefficients."""
        if not eligible(path, self.exclude_previous):
            return None
        x = tuple((d-k-1, k, 1) for d, k in path)
        if not self.tables[0].supported("P4", x):
            assert not self.tables[1].supported("P4", x)
            return None
        results = []
        for table in self.tables:
            base, extra = table.row("P4", x)
            min_charge = base / table.c.S
            slopes = [Q(0)] * 5
            for i, parent, count, weight in extra:
                slot = 0 if parent else i + 1
                assert not parent or i == 0
                slopes[slot] = Q(weight, table.c.S)
                min_charge += count * min(slopes[slot], Q(0)) / 2
            assert min_charge >= 0
            results.append((min_charge, tuple(slopes)))
        return tuple(results)

    def toll(self, path, part):
        row = self.row(path)
        if row is None:
            return NO_VALUE, NO_VALUE
        if part == "minimum":
            return (*interval(row[0][0]), *interval(row[1][0])), NO_VALUE
        idx = 0 if part == "mother" else 1 + int(part)
        slopes = [r[1][idx] for r in row]
        # For a negative slope the unspent term is |slope|*(1/2-p).
        return ((*interval(max(0, -slopes[0])/Q(2)), *interval(max(0, -slopes[1])/Q(2))),
                (*interval(slopes[0]), *interval(slopes[1])))


def extrema(toll, edges, unit, s0, ds):
    minimum_total = s0 + sum(edge[0] for edge in edges)
    out = []
    for j, source in enumerate(toll):
        upper = bool(j % 2)
        total = minimum_total
        numerator = source * unit
        choices = [(0, ds)]
        for lo, hi, future, slope in edges:
            next_value = future[j] + slope[j]
            numerator += lo * next_value
            if hi > lo:
                choices.append((next_value, hi - lo))
        for next_value, increment in sorted(choices, reverse=upper):
            if upper and next_value * total <= numerator:
                break
            if not upper and next_value * total >= numerator:
                break
            numerator += increment * next_value
            total += increment
        value = -((-numerator) // total) if upper else numerator // total
        # Verify the drift sign directly against all adapted-rate choices.
        residual = source * unit - (s0 + (ds if ((value < 0) == upper) else 0)) * value
        for lo, hi, future, slope in edges:
            delta = future[j] + slope[j] - value
            residual += (hi if ((delta > 0) == upper) else lo) * delta
        assert residual <= 0 if upper else residual >= 0
        out.append(value)
    return tuple(out)


def lifetime(scores, part, cap, arity, side_cap, mother_cap, budget):
    """One distinguished donor, summed over every possible choice and chronological order."""
    connect = int(part) if part.isdigit() else -1
    donor_needed = connect >= 0
    max_d = max(*cap, side_cap+1, mother_cap+1, 2)
    unit = lcm(1_000_000, *range(2, max_d+1))
    reciprocal = [0] + [unit//d for d in range(1, max_d+1)]
    side_tail = reciprocal[max(2, side_cap+1)]
    mother_tail = reciprocal[max(2, mother_cap+1)]
    s0, ds = int(unit*LOW), int(unit*(HIGH-LOW))
    calls = 0
    paid = set()

    @cache
    def future(path, y, mother):
        nonlocal calls
        n = len(path)
        for i, (d, k) in enumerate(path):
            demand = int(n < 4 and i == n-1) + int(i == connect and y == NO_DONOR)
            if d > cap[i] or k + demand > arity[i]:
                return NO_VALUE
        calls += 1
        if calls > budget:
            raise RuntimeError("live-state budget exceeded")
        tagged = y != NO_DONOR
        sources = end_slope = NO_VALUE
        if n == 4 and (not donor_needed or tagged) and scores.row(path) is not None:
            paid.add(path)
            sources, end_slope = scores.toll(path, part)
        edges = []

        def edge(lo, hi, dest, slope=NO_VALUE):
            assert 0 <= lo <= hi
            if hi:
                edges.append((lo, hi, dest, slope))

        for i, (d, k) in enumerate(path):
            leaves = d - k - 1
            grown = future(path[:i] + ((d+1,k),) + path[i+1:], y, mother)
            internal = leaves*unit
            if i:
                internal += reciprocal[path[i-1][0]]
            if i+1 < n:
                internal += reciprocal[path[i+1][0]]
            edge(internal, internal, grown)
            if i == 0:
                edge(reciprocal[mother], reciprocal[mother] if mother else mother_tail, grown,
                     end_slope if part == "mother" else NO_VALUE)
            known_donor = tagged and i == connect
            unseen = k - int(i+1 < n) - int(known_donor)
            assert unseen >= 0
            edge(0, unseen*reciprocal[2], grown)
            if known_donor:
                yd = y[0]
                edge(reciprocal[yd], reciprocal[yd] if yd else side_tail, grown, end_slope)
            if leaves:
                promoted = path[:i] + ((d,k+1),) + path[i+1:]
                successors = [future(promoted, y, mother)]
                if i == n-1 and n < 4:
                    successors.append(future(promoted + (VIRGIN,), y, mother))
                if i == connect and not tagged:
                    successors.append(future(promoted, VIRGIN if side_cap else OVERFLOW, mother))
                combined = tuple(map(sum, zip(*successors)))
                edge(leaves*reciprocal[d], leaves*reciprocal[d], combined)
        if tagged and y != OVERFLOW:
            d, k = y
            leaves = d-k-1
            grown = (d+1,k) if d < side_cap else OVERFLOW
            fixed = leaves*unit + reciprocal[path[connect][0]]
            edge(fixed, fixed+k*reciprocal[2], future(path, grown, mother))
            if leaves:
                edge(leaves*reciprocal[d], leaves*reciprocal[d], future(path, (d,k+1), mother))
        if mother:
            own_child = reciprocal[path[0][0]]
            edge(own_child, own_child+mother*unit-3*reciprocal[2],
                 future(path,y,mother+1 if mother < mother_cap else 0))
        return extrema(sources, edges, unit, s0, ds)

    entries = [future((VIRGIN,), NO_DONOR, m) for m in (0,*range(2,mother_cap+1))]
    bounds = []
    for j in (0, 2):
        lo = min(v[j] for v in entries)
        hi = max(v[j+1] for v in entries)
        # These individual incidence measures and rewards are positive on every tree.
        bounds.append((max(Q(0), min(LOW*lo,HIGH*lo)/GRID), max(LOW*hi,HIGH*hi)/GRID))
    print("component",part,"censuses",len(paid),"live states",calls,"residuals",4*calls,flush=True)
    future.cache_clear()
    return bounds


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--lower", type=Path, required=True)
    ap.add_argument("--upper", type=Path, required=True)
    ap.add_argument("--degrees", type=int, nargs=4, default=[6,5,4,6])
    ap.add_argument("--arities", type=int, nargs=4, default=[2,3,2,1])
    ap.add_argument("--side-cap", type=int, default=3)
    ap.add_argument("--mother-cap", type=int, default=3)
    ap.add_argument("--parts", nargs="+", choices=["minimum","mother","0","1","2","3"],
                    default=["minimum","mother","0","1","2","3"])
    ap.add_argument("--include-old-band", action="store_true")
    ap.add_argument("--budget", type=int, default=1500000)
    args = ap.parse_args()
    if any(not (2 <= d <= 16 and 0 <= k < d) for d,k in zip(args.degrees,args.arities)):
        ap.error("invalid degree/arity caps")
    if any(k < 1 for k in args.arities[:3]) or len(set(args.parts)) != len(args.parts):
        ap.error("path arities or set of objectives are invalid")
    if any(v == 1 or v < 0 for v in (args.side_cap, args.mother_cap)):
        ap.error("memory degrees must be zero or at least two")
    scores = Reward((original(args.lower, "lower"), original(args.upper, "upper")),
                    not args.include_old_band)
    print("authenticated both D18 originals", "degrees",args.degrees,"arities",args.arities,
          "side",args.side_cap,"mother",args.mother_cap,"old included",args.include_old_band,flush=True)
    totals = [[Q(0),Q(0)] for _ in range(2)]
    for part in args.parts:
        prices = lifetime(scores,part,args.degrees,args.arities,args.side_cap,args.mother_cap,args.budget)
        for j, (lo, hi) in enumerate(prices):
            totals[j][0] += lo
            totals[j][1] += hi
            print("lower" if j == 0 else "upper",lo,hi,f"[{float(lo):.9g}, {float(hi):.9g}]",flush=True)
    for name,(lo,hi) in zip(("lower","upper"),totals):
        print("sum of listed parts",name,lo,hi,f"[{float(lo):.9g}, {float(hi):.9g}]")


if __name__ == "__main__":
    main()
