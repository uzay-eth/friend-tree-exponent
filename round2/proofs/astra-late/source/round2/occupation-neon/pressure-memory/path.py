"""All-age occupation of one physical P4 census, with every side child born live."""

import argparse
from fractions import Fraction
from functools import cache
from math import ceil, floor

Q = Fraction
ZERO = (0, 0)
GRID = 10**16
LOW, HIGH = Q(565673, 10**6), Q(565775, 10**6)
TARGET = ((2, 1), (4, 3), (2, 1), (4, 0))
VIRGIN = (2, 0)
TAIL = (0, 0)
# Constant, incoming parent at A, other children at B, normalized original D18 units.
LOWER = (Q(7699396257, 12500000000000000000), Q(374338137, 10**11), Q(447940521, 10**11))
UPPER = (Q(360422966727456859, 10**20), -Q(8568217, 4000000000), -Q(253320197, 10**11))
SCORES = {
    "contrast": tuple(a - b for a, b in zip(LOWER, UPPER)),
    "mass": (Q(1), Q(0), Q(0)),
    "lower": LOWER,
    "upper": UPPER,
}


def signed_ratio(constant, known_total, terms, want_max):
    """Extremize a ratio by selecting optional rates in continuation-value order."""
    value = constant / known_total
    for successor, increment in sorted(terms, reverse=want_max):
        if (successor <= value if want_max else successor >= value):
            break
        known_total += increment
        constant += successor * increment
        value = constant / known_total
    return value


def control_row(reward, edges):
    total = LOW
    numerators = [reward * GRID, reward * GRID]
    terms = [[(Q(0), HIGH - LOW)], [(Q(0), HIGH - LOW)]]
    for lower_rate, upper_rate, successors, credit in edges:
        if not lower_rate <= upper_rate:
            raise ArithmeticError("inverted physical rate range")
        total += lower_rate
        for side, endpoint in enumerate(successors):
            action = endpoint + credit * GRID
            numerators[side] += lower_rate * action
            if upper_rate > lower_rate:
                terms[side].append((action, upper_rate - lower_rate))
    return (floor(signed_ratio(numerators[0], total, terms[0], False)),
            ceil(signed_ratio(numerators[1], total, terms[1], True)))


def stationary_enclosure(donor_cap, mother_cap, score):
    if donor_cap == 1 or mother_cap == 1 or min(mother_cap, donor_cap) < 0:
        raise ValueError("degree caps must be zero or at least two")
    new_donor = VIRGIN if donor_cap else TAIL
    tail_p = Q(1, max(2, donor_cap + 1))
    tail_mother_p = Q(1, max(2, mother_cap + 1))

    @cache
    def future(path, sides, mother):
        n = len(path)
        if any(d > goal_d or k > goal_k for (d, k), (goal_d, goal_k) in zip(path, TARGET)):
            return ZERO
        if n < len(TARGET) and path[-1][1] == TARGET[n-1][1]:
            return ZERO  # this prefix will never emit another extension
        final = path == TARGET
        intercept, parent_slope, child_slope = score if final else (0, 0, 0)
        edges = []
        for i, (d, k) in enumerate(path):
            l = d - k - 1
            gain_path = path[:i] + ((d + 1, k),) + path[i+1:]
            gain_value = future(gain_path, sides, mother)
            known = Q(l)
            if i:
                known += Q(1, path[i-1][0])
            if i + 1 < n:
                known += Q(1, path[i+1][0])
            if known:
                edges.append((known, known, gain_value, 0))
            if i == 0:
                low_p = Q(1, mother) if mother else Q(0)
                high_p = low_p if mother else tail_mother_p
                edges.append((low_p, high_p, gain_value, parent_slope))
            if i == 1 and sides:
                low_p = sum((Q(1, dd) for dd, kk in sides if dd), Q(0))
                high_p = low_p + sides.count(TAIL) * tail_p
                edges.append((low_p, high_p, gain_value, child_slope))
            if not l:
                continue
            promote_path = path[:i] + ((d, k + 1),) + path[i+1:]
            plain_sides = tuple(sorted(sides + (new_donor,))) if i == 1 else sides
            value = future(promote_path, plain_sides, mother)
            if i + 1 == n and n < len(TARGET):
                extended = future(promote_path + (VIRGIN,), sides, mother)
                value = (value[0] + extended[0], value[1] + extended[1])
            birth = Q(l, d)
            edges.append((birth, birth, value, 0))

        for i, (d, k) in enumerate(sides):
            if not d:
                continue
            l = d - k - 1
            next_donor = (d + 1, k) if d < donor_cap else TAIL
            gain_sides = tuple(sorted(sides[:i] + (next_donor,) + sides[i+1:]))
            rate = l + Q(1, path[1][0])
            edges.append((rate, rate + Q(k, 2), future(path, gain_sides, mother), 0))
            if l:
                promote_sides = tuple(sorted(sides[:i] + ((d, k + 1),) + sides[i+1:]))
                birth = Q(l, d)
                edges.append((birth, birth, future(path, promote_sides, mother), 0))
        if mother:
            next_mother = mother + 1 if mother < mother_cap else 0
            own_child = Q(1, path[0][0])
            # The mother itself has one core parent, in addition to the tracked A.
            edges.append((own_child, own_child + mother - Q(3, 2),
                          future(path, sides, next_mother), 0))
        return control_row(intercept, edges)

    initial_degrees = (0,) + tuple(range(2, mother_cap + 1))
    entries = {d: future((VIRGIN,), (), d) for d in initial_degrees}
    interval = (min(v[0] for v in entries.values()), max(v[1] for v in entries.values()))
    endpoints = [s * Q(v, GRID) for s in (LOW, HIGH) for v in interval]
    count = future.cache_info().currsize
    future.cache_clear()
    return min(endpoints), max(endpoints), count, entries


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--donor", type=int, default=2)
    parser.add_argument("--mother", type=int, default=3)
    parser.add_argument("--score", choices=SCORES, default="contrast")
    parser.add_argument("--theta", type=Q, default=Q(1),
                        help="upper-witness multiplier in the contrast")
    args = parser.parse_args()
    score = SCORES[args.score]
    if args.score == "contrast":
        score = tuple(a - args.theta * b for a, b in zip(LOWER, UPPER))
    low, high, states, entries = stationary_enclosure(args.donor, args.mother, score)
    print(args.score, "donor", args.donor, "mother", args.mother, "states", states)
    print("Palm birth degree: discounted endpoint ticks", entries)
    print("stationary lower and upper (exact)", low, high)
    print("outward at 1e-15", floor(low * 10**15) / 10**15, ceil(high * 10**15) / 10**15)
