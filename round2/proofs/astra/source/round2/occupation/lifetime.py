#!/usr/bin/env python3
"""Lifetime ceiling for an upper comparison profile in common integer units."""
import argparse
import json
import sys
from array import array
from fractions import Fraction
from pathlib import Path


def ceildiv(x, y):
    return -(-x // y)


def ceiling(values, maximum, degree, cores, l):
    boundary = ceildiv(maximum * l.denominator, l.numerator)
    next_row = [boundary] * (cores + 2)
    for d in range(degree, 1, -1):
        this = [boundary] * (cores + 2)
        for c in range(min(d, cores), 0, -1):
            ceiling_source = max(0, values[d*(d-1)//2-1+c-1])
            leaves = d-c
            bud = 2*leaves
            ratios = []
            for rate in (2*d*leaves, 2*d*leaves + c*d):
                den = l.numerator*2*d + l.denominator*(rate + bud)
                num = l.denominator*(2*d*ceiling_source + rate*next_row[c] + bud*this[c+1])
                ratios.append(ceildiv(num, den))
            this[c] = min(boundary, max(ratios))
        next_row = this
    return next_row[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("profile", type=Path)
    p.add_argument("comparison", type=Path)
    p.add_argument("--degree", type=int, default=1399)
    p.add_argument("--cores", type=int, default=60)
    p.add_argument("--discount", type=Fraction, default=Fraction(113,200))
    args = p.parse_args()
    meta = json.loads(args.comparison.read_text())
    values = array("q")
    with args.profile.open("rb") as stream:
        values.fromfile(stream, meta["finite_rows"])
        assert not stream.read(1)
    if sys.byteorder != "little":
        values.byteswap()
    assert min(values) >= 0 and max(values) <= meta["upper"]
    assert 2 <= args.cores <= args.degree and args.degree*(args.degree+1)//2-1 <= len(values)
    value = ceiling(values, meta["upper"], args.degree, args.cores, args.discount)
    print(f"D <= {args.degree}, c <= {args.cores}, discount {args.discount}, "
          f"newborn lifetime ceiling {value}/{meta['unit']}")


if __name__ == "__main__":
    main()
