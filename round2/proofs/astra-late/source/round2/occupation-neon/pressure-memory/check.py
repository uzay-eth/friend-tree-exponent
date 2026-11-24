"""Rational cutoff budgets for unknown donor reciprocals in the stationary core law."""

from fractions import Fraction
from itertools import product

Q = Fraction
# The earlier independently checked prior suffices for all prices here.
RATE = Q(565673, 10**6)
GRID = 10**24
PARENT_SLOPE = Q("0.00374338137") + Q("0.00214205425")
CHILD_SLOPE = Q("0.00447940521") + Q("0.00253320197")


def ceil_div(numerator, denominator):
    return -(-numerator // denominator)


def tick_up(value):
    return Q(ceil_div(value.numerator * GRID, value.denominator), GRID)


def donor_tail_budgets(cuts, rate=RATE):
    """Upper prices for E[w_J(D)] and E[k w_J(D)], with w_J=(1/(J+1)-1/D)_+."""
    if not cuts or any(j < 2 for j in cuts) or rate <= 0:
        raise ValueError("positive rate and donor cutoffs of at least two required")
    p_bound = GRID
    mean_lower = GRID * rate.denominator // rate.numerator
    mean_upper = mean_lower + (Q(mean_lower, GRID) < 1 / rate)
    prices = {}
    for degree in range(2, max(cuts) + 2):
        gain_numerator = (2 * degree - 1) * rate.denominator
        discounted = gain_numerator + 2 * rate.numerator
        p_bound = ceil_div(p_bound * gain_numerator, discounted)
        holding = 2 * rate.denominator * GRID
        mean_lower += holding // discounted
        mean_upper += ceil_div(holding, discounted)
        cutoff = degree - 1
        if cutoff not in cuts:
            continue
        price = Q(p_bound, GRID)
        factor = 1 + rate * degree / (degree + Q(1, 2) + rate)
        correction = factor - Q(GRID, mean_lower)
        if correction <= 0:
            raise ValueError("weighted telescoping factor is nonpositive")
        prices[cutoff] = (
            tick_up(price / (degree * factor)),
            tick_up(price * Q(mean_upper, GRID) / (degree * correction)),
        )
    return prices


def capped_p4_error(prices):
    child, parent = prices
    # Fixing a hidden child leaves at most two choices for the visible continuation.
    # A hidden parent can contribute to at most three occurrences per core child.
    return 2 * CHILD_SLOPE * child + 3 * PARENT_SLOPE * parent


def chebyshev(order, z):
    old, current = 1, z
    if order == 0:
        return old
    for _ in range(1, order):
        old, current = current, 2 * z * current - old
    return current


def slab_reciprocal(degree, left, order):
    """The polynomial (1 - T_n(3-2d/a)/T_n(3))/d on [a,2a]."""
    z = 3 - Q(2 * degree, left)
    return (1 - Q(chebyshev(order, z), chebyshev(order, 3))) / degree


def controls():
    for rate in (Q(1, 2), RATE, Q(3, 4)):
        tail = Q(1)
        tilted_mean = 1 / rate
        prev = {}
        for degree in range(2, 85):
            q = degree - Q(1, 2)
            tail *= q / (q + rate)
            tilted_mean += 1 / (q + rate)
            prev[degree] = (tail, tilted_mean)
        for j in range(2, 80):
            p, u = prev[j]
            p_next, u_next = prev[j + 1]
            a = 1 + rate * j / (j + Q(1, 2) + rate)
            assert (p / j - p_next / (j + 1)) * j * (j + 1) == p * a
            weighted = (p * u / j - p_next * u_next / (j + 1)) * j * (j + 1)
            assert p * (a * u - 1) <= weighted < p * a * u
        cuts = {4, 8, 16, 32}
        prices = donor_tail_budgets(cuts, rate)
        for j in cuts:
            part = sum(prev[d][0] / (d * (d + 1)) for d in range(j + 1, 85))
            part_with_k = sum(prev[d][0] * prev[d][1] / (d * (d + 1))
                              for d in range(j + 1, 85))
            assert part < prices[j][0] and part_with_k < prices[j][1]
    # Small chronological core trees test both incidence bounds. A larger
    # labelled fixture below attains the child and parent multiplicities.
    def path_counts(parents, leaves):
        n = len(parents)
        children = [[] for _ in range(n)]
        for y in range(1, n):
            children[parents[y]].append(y)
        states = [(len(children[y]), leaves[y]) for y in range(n)]
        child_multiplicity = [0] * n
        parent_multiplicity = [0] * n
        for a in range(1, n):
            if states[a] != (1, 0):
                continue
            b = children[a][0]
            if states[b] != (3, 0):
                continue
            for c in children[b]:
                if states[c] != (1, 0):
                    continue
                d = children[c][0]
                if states[d] != (0, 3):
                    continue
                parent_multiplicity[parents[a]] += 1
                for y in children[b]:
                    if y != c:
                        child_multiplicity[y] += 1
        assert all(m <= 2 for m in child_multiplicity)
        assert all(m <= 3 * len(children[y]) for y, m in enumerate(parent_multiplicity))
        return child_multiplicity, parent_multiplicity

    labelled = 0
    for size in range(2, 7):
        for prefix in product(*(range(i) for i in range(1, size))):
            par = (-1,) + prefix
            kids = [par.count(i) for i in range(size)]
            leaves = tuple(0 if v else 3 for v in kids)
            path_counts(par, leaves)
            labelled += 1
    par = (-1, 0, 1, 2, 2, 2, 3, 4, 5)
    child, parent = path_counts(par, (1, 0, 0, 0, 0, 0, 3, 3, 3))
    assert max(child) == 2 and parent[0] == 3
    # Decreasing a reciprocal at a virgin's registration would introduce false deficit.
    for degree in range(2, 50):
        deficit = Q(1, 2) - Q(1, degree)
        assert deficit == sum((Q(1, j * (j + 1)) for j in range(2, degree)), Q(0))
    # n-th differences vanish for the degree-(n-1) polynomial on every slab.
    for order in range(1, 12):
        for left in (16, 32, 64):
            values = [slab_reciprocal(d, left, order) for d in range(left, 2 * left + 1)]
            assert all(abs(1 - d * q) <= Q(1, chebyshev(order, 3))
                       for d, q in zip(range(left, 2 * left + 1), values))
            for _ in range(order):
                values = [b - a for a, b in zip(values, values[1:])]
            assert all(v == 0 for v in values)
    print(f"PASS: 234 telescopes, 12 partial-tail comparisons, {labelled} core shapes,"
          " both sharp P4 incidence factors and 33 reciprocal polynomial controls")


if __name__ == "__main__":
    controls()
    cuts = {18, 64, 128, 512, 1024, 4096, 8192, 16384, 32768, 65536, 131072, 262144}
    print("cutoff   child tail        arity-weighted tail   P4 contrast cost")
    budgets = donor_tail_budgets(cuts)
    for cutoff, prices in budgets.items():
        child, parent = prices
        value = capped_p4_error(prices)
        print(f"{cutoff:6d}   {float(child):.8e}    {float(parent):.8e}       {float(value):.8e}")
    cutoff, order = 65536, 10
    mass_ceiling = CHILD_SLOPE + 3 * PARENT_SLOPE * (1 - RATE)
    error = capped_p4_error(budgets[cutoff]) + mass_ceiling / chebyshev(order, 3)
    assert error < Q(6, 10**9)
    print(f"PASS: full two-role P4 donor approximation < 6e-9 per core"
          f" ({float(error):.10e}; T_{order}(3)={chebyshev(order, 3)})")
