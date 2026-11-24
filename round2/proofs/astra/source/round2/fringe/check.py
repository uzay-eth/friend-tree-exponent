#!/usr/bin/env python3
"""Replay the rational local inequalities and symbolic fringe balances."""

from collections import Counter
from fractions import Fraction as F
from functools import cache

DAMPING = F(3, 20)
PROMOTION_BOUND = F(9, 25)
RESIDUAL_BOUND = F(51, 100)


def h(core_degree):
    return (2, 1, F(17, 20), F(39, 50), F(7, 10))[min(core_degree, 4)]


def frac_mark(denominator):
    return F(1, denominator) if denominator else F(0)


# A state is (leaves, birth-mark denominator, sorted tuple of core children).
# Zero is the special mark of the physical centre and of an auxiliary J_0.
def node(leaves, mark, children=()):
    return leaves, mark, tuple(sorted(children))


def degree(state, planted):
    return state[0] + len(state[2]) + planted


def gain(state):
    return node(state[0] + 1, *state[1:])


def newborn(mark):
    return node(1, mark)


def promote(state, planted):
    return node(state[0] - 1, state[1], state[2] + (newborn(degree(state, planted)),))


def replace_at(state, path, replacement):
    if not path:
        return replacement
    children = list(state[2])
    i = path[0]
    children[i] = replace_at(children[i], path[1:], replacement)
    return node(*state[:2], children)


def cores(state, path=(), planted=0, parent_degree=None):
    yield path, state, planted, parent_degree
    for i, child in enumerate(state[2]):
        yield from cores(child, path + (i,), 1, degree(state, planted))


@cache
def physical_transitions(state):
    """Accumulate all choices of a clock and neighbour in the physical tree."""
    rates = Counter()
    for path, x, planted, _ in cores(state):
        leaves, _, children = x
        inv = F(1, degree(x, planted))
        if leaves:
            rates[replace_at(state, path, gain(x))] += leaves
            rates[replace_at(state, path, promote(x, planted))] += leaves * inv
        for i, y in enumerate(children):
            rates[replace_at(state, path + (i,), gain(y))] += inv
            rates[replace_at(state, path, gain(x))] += F(1, degree(y, 1))
    return rates


@cache
def intrinsic_transitions(state):
    """Intrinsic changes from the target's gain rate and the promoter's rate."""
    rates = Counter()
    for path, x, planted, parent in cores(state, planted=1):
        leaves, _, children = x
        inv = F(1, degree(x, planted))
        r = leaves + sum(F(1, degree(y, 1)) for y in children)
        if parent is not None:
            r += F(1, parent)
        if r:
            rates[replace_at(state, path, gain(x))] += r
        if leaves:
            rates[replace_at(state, path, promote(x, planted))] += leaves * inv
    return rates


@cache
def occupation(state):
    return Counter(x for path, x, _, _ in cores(state) if path)


def subset_row(row, x, planted, omit_old_credit=False):
    inv = F(1, degree(x, planted))
    row[newborn(degree(x, planted))] += x[0] * inv
    for y in x[2]:
        mark = frac_mark(y[1])
        assert mark >= inv
        row[gain(y)] += inv
        if not omit_old_credit:
            row[y] += mark - inv


@cache
def mass(state):
    """The physical weighted descendant sum before division at the root."""
    c = 1 + len(state[2])
    return h(c) + sum(frac_mark(y[1]) * weight(y) for y in state[2])


@cache
def weight(state):
    return mass(state) / (frac_mark(state[1]) + DAMPING)


def check_physical_generator(state):
    left = Counter()
    initial = occupation(state)
    for successor, rate in physical_transitions(state).items():
        for y, count in occupation(successor).items():
            left[y] += rate * count
        for y, count in initial.items():
            left[y] -= rate * count

    internal = Counter()
    for x, count in initial.items():
        change = Counter()
        for y, rate in intrinsic_transitions(x).items():
            change[y] += rate
            change[x] -= rate
        change[x] -= frac_mark(x[1])
        subset_row(change, x, 1)
        for y, value in change.items():
            internal[y] += count * value

    right = internal.copy()
    subset_row(right, state, 0)
    assert left == right
    wrong = internal.copy()
    subset_row(wrong, state, 0, omit_old_credit=True)
    return left != internal, left != wrong


@cache
def check_residual(state):
    a = frac_mark(state[1])
    w = weight(state)
    intrinsic = sum(rate * (weight(y) - w) for y, rate in intrinsic_transitions(state).items())
    copied = sum(frac_mark(y[1]) * weight(y) for y in state[2])
    assert all(weight(gain(y)) == weight(y) for y in state[2])
    assert copied == (a + DAMPING) * w - h(1 + len(state[2]))
    assert intrinsic <= PROMOTION_BOUND * w
    assert intrinsic - a * w + copied <= RESIDUAL_BOUND * w - h(1 + len(state[2]))
    return True


def main():
    small = 0
    largest = F(0)
    argmax = None
    for d in range(2, 100):
        for c in range(d + 1):
            lhs = (1 - F(c, d)) * (h(c + 1) - h(c) + 1 / (1 + DAMPING * d))
            ratio = lhs / h(c)
            if ratio > largest:
                largest, argmax = ratio, (d, c)
            assert lhs <= PROMOTION_BOUND * h(c)
            small += 1
    assert F(1, 16) < PROMOTION_BOUND * h(4)
    print(f"Local promotion inequality: {small} finite rows. Maximum ratio {largest} at (D,c)={argmax}.")
    print("All D >= 100: 1/16 < (9/25)(7/10).")

    law = {node(2, 0): F(1)}
    tested = 0
    mutations = [0, 0]
    target = [F(1), F(4, 3), F(19, 12), F(647, 360), F(17189, 8640)]
    mean_list = []
    for n in range(3, 10):
        next_law = Counter()
        mean = F(0)
        for state, prob in law.items():
            tested += 1
            rejected = check_physical_generator(state)
            mutations = [total + new for total, new in zip(mutations, rejected)]
            for y in occupation(state):
                assert check_residual(y)
            mean += prob * (1 + occupation(state).total())
            events = physical_transitions(state)
            assert sum(events.values()) == n
            if n < 9:
                for nxt, rate in events.items():
                    next_law[nxt] += prob * rate / n
        assert sum(law.values()) == 1
        if n < 8:
            assert mean == target[n - 3]
        mean_list.append(str(mean))
        law = next_law
    assert all(mutations)
    print(f"Physical generators: {tested} sorted marked states, n=3..9. All symbolic occupation differences vanish.")
    print(f"Residual drift and copy-mass identity: {check_residual.cache_info().currsize} distinct planted fringes.")
    print("Mean core counts, n=3..9: " + ", ".join(mean_list))
    print(f"Negative controls: missing centre immigration rejected on {mutations[0]} rows; omitted old-copy credit on {mutations[1]}.")


if __name__ == "__main__":
    main()
