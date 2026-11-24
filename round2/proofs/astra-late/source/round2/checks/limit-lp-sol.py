"""Exact beta normalization, crossing drift, and mesh exponent controls."""

from fractions import Fraction as Q
from math import comb, factorial


def main():
    rows = 0
    for a in range(2, 13):
        for b in range(2, 13):
            n = a + b - 2
            norm = Q(factorial(a + b - 1), factorial(a - 1) * factorial(b - 1))
            assert norm == (n + 1) * comb(n, a - 1)
            assert norm * Q(factorial(a - 1) * factorial(b - 1), factorial(n + 1)) == 1
            for x in (Q(0), Q(1, 10), Q(1, 2), Q(9, 10), Q(1)):
                atoms = [comb(n, k) * x**k * (1 - x)**(n - k) for k in range(n + 1)]
                assert sum(atoms) == 1
                assert norm * x**(a - 1) * (1 - x)**(b - 1) == (n + 1) * atoms[a - 1]
                rows += 1
    print(f"{rows} exact beta/binomial density rows, including endpoints")

    def f(j):
        return min(1, Q(2, j + 1))

    for j in range(1, 201):
        drift = j * (f(j + 1) - f(j)) + Q(j, j + 1) * (f(j - 1) - f(j)) + Q(1, j + 1) * (1 - f(j))
        assert drift == (Q(-1, 3) if j == 1 else Q(2 - j, (j + 1) * (j + 2))) <= 0
    print("200 exact first-send drift rows with zero immigration")

    for mu in (Q(51, 100), Q(3, 5), Q(9, 10)):
        for c in (Q(1, 1000), Q(1, 20), Q(1, 8)):
            assert (1 - Q(1, 4)) / 2 - 2 * Q(1, 16) == Q(1, 4)
            xi = 1 / (2 + mu + 2/c)
            a, b = xi/c, mu*xi
            assert c*a == (1 - b) / 2 - a == b/mu == xi
    print("9 exact optimized mesh balances; original inverse contribution exponent 1/4")


if __name__ == "__main__":
    main()
