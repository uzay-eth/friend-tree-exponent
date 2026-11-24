"""Finite arithmetic for the degree-two urn cap and a toy marked Mellin law."""

from fractions import Fraction as Q
from math import comb


def main():
    rows = 0
    for m in range(3, 28):
        for i in range(1, 20):
            x = Q(i, 20)
            # Integrate (m-1)(m-2)t(1-t)^(m-3) from x to 1
            # by polynomial expansion, independent of the proposed formula.
            tail = (m - 1) * (m - 2) * sum(
                Q((-1)**j * comb(m - 3, j), j + 2) * (1 - x**(j + 2))
                for j in range(m - 2)
            )
            cap = (1 - x)**(m - 2) * (1 + (m - 2) * x)
            next_cap = (1 - x)**(m - 1) * (1 + (m - 1) * x)
            assert tail == cap and cap - next_cap == (m - 1) * x*x * (1 - x)**(m - 2)
            rows += 1
    print(f"{rows} exact beta-tail and monotonicity rows")

    # For mu=1/2 and atoms at squares all integrals here are rational.
    # Integrate in s=sqrt(y): y^(mu-1)dy=2ds.
    atoms = [(1, Q(2)), (2, Q(3)), (3, Q(1))]  # (sqrt(u), intensity)
    mellin = 2 * sum(s * w for s, w in atoms)
    assert mellin == 22
    for k in range(1, 21):  # R > k^2
        masses = [2 * (Q(s, k)) * w / mellin for s, w in atoms]
        assert sum(masses) == Q(1, k)
        assert [m / sum(masses) for m in masses] == [2 * s * w / mellin for s, w in atoms]
    print("20 exact marked-Mellin rows: Pareto tail and unchanged biased atom weights")


if __name__ == "__main__":
    main()
