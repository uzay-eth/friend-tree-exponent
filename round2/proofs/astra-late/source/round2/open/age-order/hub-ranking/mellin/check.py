"""Directed bounds for the degree-mass moment after a serial or parallel placement."""

import argparse
from fractions import Fraction as Q
from math import ceil, factorial
from pathlib import Path

GRID = 10**15
HERE = Path(__file__).resolve().parent


def middle(lo, hi, degree=30, return_cap=40):
    """Downward bound, normalized by Gamma(1+p), for every lo <= p <= hi."""
    assert 1 <= lo <= hi
    q, z = hi.numerator, hi.denominator
    boundary = Q(1)
    after = [0, 0, GRID]
    for leaves in range(3, degree + 1):
        boundary *= (lo + leaves - 2) / (leaves - 2)
        after.append((GRID * boundary).__floor__())
    for d in range(degree-1, 1, -1):
        rescue = 0
        for k in range(return_cap-1, 1, -1):
            growth = k * (d*(k-1) + 1)
            rescue = z*(d*GRID + growth*rescue) // (z*(d+growth) + d*k*q)
        current = [rescue * after[1] // GRID]
        for leaves in range(1, d+1):
            numerator = leaves*z*(d*after[leaves+1] + current[leaves-1])
            current.append(numerator // (d*q + leaves*(d+1)*z))
        after = current
    return Q(after[1], GRID)


def beta_moment(alpha, order):
    value = Q(1)
    for j in range(order):
        value *= (alpha+j) / (1+j)
    return value


def positive_tip(lo, hi, cutoff=48):
    """Two parent-class ceilings for every exponent in an interval."""
    n = lo.__floor__()
    assert lo >= 1 and hi <= n+1

    def promoted(k):
        r = Q(1, k)
        # A_p(r) is a decreasing convex beta moment when 0<r<1.
        upper = beta_moment(r, n) * (1 - (lo-n)*(1-r)/(n+1))
        return (1+hi/r) * upper

    tip = {}
    tip[cutoff] = ceil(GRID * promoted(cutoff) / (1+cutoff*lo))
    for k in range(cutoff-1, 1, -1):
        bound = (GRID * promoted(k) + k*tip[k+1]) / (1+k*(1+lo))
        tip[k] = max(tip[k+1], ceil(bound))
    return Q(tip[2], GRID), Q(tip[3], GRID)


def checks(export):
    witnesses = [(p, middle(Q(p), Q(p))) for p in (1,2,3)]
    for (p, v), threshold in zip(witnesses, (Q(".56319"), Q(".50833"), Q(".48054")), strict=True):
        assert v > threshold
        print("serial middle", p, v, flush=True)

    # The low-power ledger reads the same 256 directed boxes from a plain TSV.
    low_rows = []
    for j in range(128, 384):
        lo, hi = Q(j,128), Q(j+1,128)
        low_rows.append(f"{lo} {hi} {middle(lo,hi)}\n")
    data = "".join(low_rows).encode()
    path = HERE / "lower-1-3.tsv"
    if export:
        path.write_bytes(data)
    else:
        assert data == path.read_bytes()

    margins = []
    for j in range(3*512, 7*512):
        lo, hi = Q(j,512), Q(j+1,512)
        m = middle(lo,hi)
        tip2, tip3 = positive_tip(lo,hi)
        margins.append(m - (4*tip2+2*tip3)/3)
    worst = min(margins)
    assert worst > 0, (worst, margins.index(worst))

    high = sum((Q((-1)**j, factorial(j)) for j in range(8)), Q())
    debit = (4*beta_moment(Q(1,2), 7) + 2*beta_moment(Q(1,3), 7))/3
    assert high - debit > Q(21,1000)
    print("PASS 256 lower-power boxes and 2048 signed boxes on [3,7], minimum", worst)
    print("PASS p>=7 normalized serial lower", high, "signed margin at least", high-debit)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--export", action="store_true")
    checks(parser.parse_args().export)
