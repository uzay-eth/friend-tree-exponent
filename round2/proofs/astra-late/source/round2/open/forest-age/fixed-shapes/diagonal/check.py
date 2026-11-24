"""Rational controls for the residence sandwich and rare-component sampling."""

from fractions import Fraction as Q
from math import comb


def exposure(left, right, entrance, exit, eps):
    """Length of allowed birth times for an age-residence interval."""
    keep = max(Q(0), min(right, -entrance) - max(left, eps - exit))
    seen = max(Q(0), min(right, eps - entrance) - max(left, -exit))
    now = max(Q(0), min(right, -entrance) - max(left, -exit))
    return keep, now, seen


def main():
    windows = loss_mutants = hits = 0
    for left in (-10, -3, -1):
        for right in (Q(-1, 2), 0, Q(1, 2)):
            for start in (0, Q(1, 3), 2):
                for span in (Q(1, 7), Q(1, 2), 3):
                    for eps in (Q(1, 13), Q(2, 3)):
                        keep, now, seen = exposure(left, right, start, start + span, eps)
                        assert 0 <= keep <= now <= seen
                        assert now - keep <= eps and seen - now <= eps
                        windows += 1
                        no_late = exposure(left, min(right, 0), start, start + span, eps)[2]
                        loss_mutants += no_late != seen
    for total in (12, 50, 150):
        for red in (1, 2, 5):
            for draws in (1, 4, 8):
                exact = Q(comb(total - red, draws), comb(total, draws))
                prod = Q(1)
                for j in range(draws):
                    prod *= 1 - Q(red, total - j)
                iid = (1 - Q(red, total)) ** draws
                assert exact == prod <= iid
                assert iid - exact <= Q(red * draws * (draws - 1), 2 * total * (total - draws + 1))
                hits += 1
    assert loss_mutants > 0
    print(f'{windows} rational residence sandwiches; {loss_mutants} omitted-late-birth controls')
    print(f'{hits} exact hypergeometric checks with finite-population error')


if __name__ == '__main__':
    main()
