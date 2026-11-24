"""Exact arithmetic for the conditional Yule transition and oscillation control."""

from fractions import Fraction as Q
from math import comb, factorial


def main():
    # At times log(a), log(b), the Yule state j becomes k with the
    # negative-binomial kernel. Bayes using the respective gamma densities
    # must give a Poisson increment of mean (b-a)*u. Exponents cancel.
    rows = 0
    for a, b in ((1, 2), (2, 3), (2, 7), (4, 9)):
        for j in range(1, 10):
            for k in range(j, j + 10):
                for u in (Q(1, 3), Q(1), Q(5, 2)):
                    transition = comb(k - 1, j - 1) * Q(a, b)**j * (1 - Q(a, b))**(k-j)
                    posterior_ratio = (b**k * u**(k-1) / factorial(k-1)) / (
                        a**j * u**(j-1) / factorial(j-1))
                    assert transition * posterior_ratio == ((b-a)*u)**(k-j) / factorial(k-j)
                    rows += 1
    print(f"{rows:,} exact Yule/Gamma Bayes rows match Poisson increments")
    for x in (Q(1, 5), Q(1, 2), Q(4, 5)):
        q = 1 - x
        values = [(q if v % 2 == 0 else 1/q) * q**v for v in range(2, 16)]
        assert all(a < b for a, b in zip(values[::2], values[1::2]))
    print("3 exact oscillation controls reverse adjacent signs without changing the rate")


if __name__ == "__main__":
    main()
