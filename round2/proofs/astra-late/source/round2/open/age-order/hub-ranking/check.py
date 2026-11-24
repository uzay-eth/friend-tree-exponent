"""Exact boundary cases and consequences of the all-tree hub gap."""

from fractions import Fraction
from pathlib import Path
from runpy import run_path


def main():
    hostile = Path(__file__).resolve().parents[1] / "global-transport" / "hostile"
    scalar = run_path(str(hostile / "discharge" / "check.py"))
    tree = scalar["H"]
    floor = Fraction(371, 12500)
    assert floor == Fraction(17, 1000) + 2 * Fraction(634, 100000)

    for degree, boundary_ticks in enumerate((18940948, 20181908, 22577249), 1):
        star = tree["spider"](degree, 1)
        total_weight = sum(tree["profile"](star)[1])
        score = tree["score"](star)[0] * total_weight
        assert score > Fraction(boundary_ticks, 10**8) > floor

    lower_leaf = scalar["L"](2)
    assert lower_leaf > Fraction(1, 4)
    assert scalar["G"](2) / 2 == Fraction(189409485129, 2 * 10**12) > floor
    barrier_checks = 0
    for degree in range(3, 100):
        for leaves in range(1, degree - 1):
            value = Fraction(1, leaves + 1)
            gain = degree * (leaves + 1)
            numerator = 1 + gain * Fraction(1, leaves + 2) + leaves * Fraction(1, leaves)
            defect = numerator - (1 + gain + leaves) * value
            assert defect == 1 - Fraction(degree, leaves + 2) <= 0
            barrier_checks += 1
    assert (Fraction(3, 4) - (1 - Fraction(108, 305))) / 2 == Fraction(127, 2440)

    telescopes = 0
    for earlier in range(3, 65):
        summed_bound = Fraction(0)
        for later in range(earlier + 1, earlier + 81):
            old_size = later - 2
            summed_bound += floor / (old_size * (old_size + 1))
            assert summed_bound == floor * (Fraction(1, earlier - 1) - Fraction(1, later - 1))
            telescopes += 1
    print(f"PASS seed, stars, founder gap, {barrier_checks} first-send rows and {telescopes} label intervals")


if __name__ == "__main__":
    main()
