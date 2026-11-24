"""Exact physical jump-chain fixture for the insertion forecast identity."""

from fractions import Fraction as Q
from functools import cache


@cache
def rows(tree):
    n = len(tree)
    out = []
    for p, ns in enumerate(tree):
        prob = sum((Q(1, n * len(tree[w])) for w in ns), Q(0))
        nxt = tuple(xs + (n,) if w == p else xs for w, xs in enumerate(tree)) + ((p,),)
        out.append((prob, nxt))
    assert sum(p for p, _ in out) == 1
    return tuple(out)


def fixture(first, last, terminal, threshold):
    @cache
    def endpoint(tree, v):
        if len(tree) == terminal:
            return Q(len(tree[v]) >= threshold)
        return sum((p * endpoint(nxt, v) for p, nxt in rows(tree)), Q(0))

    @cache
    def forecast(tree, left):
        return sum((p * endpoint(nxt, len(tree)) for p, nxt in rows(tree)), Q(0)) if left <= len(tree) + 1 <= last else Q(0)

    @cache
    def remaining(tree, left):
        if len(tree) >= last:
            return Q(0)
        return forecast(tree, left) + sum((p * remaining(nxt, left) for p, nxt in rows(tree)), Q(0))

    def future_direct(tree, left):
        if len(tree) == terminal:
            return Q(sum(len(tree[v - 1]) >= threshold for v in range(left, last + 1)))
        return sum((p * future_direct(nxt, left) for p, nxt in rows(tree)), Q(0))

    checked = caught = 0

    def visit(tree):
        nonlocal checked, caught
        n = len(tree)
        left = max(first, n + 1)
        assert future_direct(tree, left) == remaining(tree, left)
        checked += 1
        if n < last:
            for _, nxt in rows(tree):
                # The new vertex is in the terminal count's jump, but it must
                # not be included in the future forecast integral after birth.
                true = remaining(nxt, max(first, n + 2))
                wrong = true + (endpoint(nxt, n) if first <= n + 1 <= last else 0)
                caught += wrong != true
                visit(nxt)

    visit(((1, 2), (0,), (0,)))
    assert caught > 0
    return checked, caught


if __name__ == '__main__':
    cases = [(4, 5, 7, 2), (5, 6, 8, 3), (4, 5, 8, 4)]
    rows_checked, wrong = map(sum, zip(*(fixture(*case) for case in cases)))
    print(f'{rows_checked} exact conditional forecast rows; {wrong} already-inserted-tip controls caught')
