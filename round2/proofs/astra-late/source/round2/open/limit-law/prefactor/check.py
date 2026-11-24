#!/usr/bin/env python3
"""Exact finite controls for the departure law, coupling, and tail scaling."""
from fractions import Fraction as Q


def main():
    survive = Q(1)
    for k in range(4, 257):
        depart = survive / (k - 1)
        survive -= depart
        assert depart == Q(2, (k - 2) * (k - 1))
        assert survive == Q(2, k - 1)
    print('PASS: exact departure and survival weights through size 256')

    # At a hub ring, c core neighbours each incur precisely 1/d failure
    # rate. A rejection-free promotion increments c but leaves d intact;
    # an inert root sends to its slot without adding internally.
    cases = 0
    for d in range(2, 80):
        for c in range(1, d + 1):
            slots = [Q(1, d)] * d
            failure = sum(slots[:c])
            success = sum(slots[c:])
            assert failure == Q(c, d) and success == 1 - failure
            if c < d:
                assert sum(slots[:c + 1]) - failure == Q(1, d)
            cases += 1
    print(f'PASS: {cases} exact hub degree/core-neighbour states')

    # Toy poststar law at mu=1/2: the scaled limit D takes 1,2 with
    # probabilities 1/3,2/3, so the predicted prefactor is exactly 1.
    # For h=x^2=1/m the two departure thresholds are m and 4m.
    for m in (8, 32, 128, 512):
        tail = lambda j: Q(2, j - 2)
        end = 8*m
        truncated = sum(
            Q(2, (k-2)*(k-1)) * (Q(1, 3)*(k >= m) + Q(2, 3)*(k >= 4*m))
            for k in range(4, end + 1)
        )
        scaled = m * (truncated + tail(end+1))
        assert scaled == m * (tail(m)/3 + 2*tail(4*m)/3)
        assert 0 < scaled - 1 < Q(4, m)
    print('PASS: exact fixed-size mixture with a two-atom toy amplitude')


if __name__ == '__main__':
    main()
