"""Rational checks of the Palm transform and ordered sibling sampling."""

from fractions import Fraction as F
from itertools import permutations, product
from math import prod


def check_palm():
    # Use synthetic mu=1/2 and perfect-square positive amplitudes so the
    # primitive z=exp(-mu*y) has a rational cutoff sqrt(w). Include W=0.
    amplitudes = [(0, F(1,6)), (1, F(1,3)), (4, F(1,4)), (9, F(1,4))]
    roots = {0:0, 1:1, 4:2, 9:3}
    normalizer = sum(p*roots[w] for w,p in amplitudes)
    assert normalizer == F(19,12)
    tests = 0
    for w,p in amplitudes:
        for c,d in ((F(0),F(1)), (F(1,3),F(2,3)), (F(2,3),F(1))):
            for a,b in ((F(0),F(1)), (F(1,7),F(3,7))):
                # z uniform in the unnormalized parent measure. In the
                # product law c<exp(-mu*E)<d and a<V<b.
                actual = p*(min(roots[w], d*roots[w])-c*roots[w])*(b-a)
                target = normalizer*(p*roots[w]/normalizer)*(d-c)*(b-a)
                assert actual == target
                if w > 1 and actual:
                    assert actual != p*(d-c)*(b-a)  # rejects missing tilt
                tests += 1
    # Different W values have different centered-age support even though
    # their overshoot laws agree; positivity cutoff cannot be ignored.
    assert sum(p*max(0,min(1,roots[w])-F(1,2)) for w,p in amplitudes) == F(5,12)
    print(f'Palm transform: {tests} rational boxes; normalizer {normalizer}; missing-tilt control PASS')


def check_siblings():
    # Independent cells of marked Poisson proposals. Enumerate proposal
    # counts up to 8 and compare three members to shifted survival sums.
    from math import factorial
    intensities = (F(1,3),F(1,2),F(1,4))
    tests = 0
    for m in (1,2,3):
        for choices in product(range(len(intensities)),repeat=m):
            used = tuple(choices.count(i) for i in range(len(intensities)))
            # For each bin independently enumerate the number of ordered
            # assignments of distinct proposal IDs, and multiply the
            # unnormalized Poisson weights (exp(-sum intensities) cancels).
            actual = F(1)
            for x,k in zip(intensities,used):
                actual *= sum(F(len(tuple(permutations(range(n),k))),factorial(n))*x**n
                              for n in range(9))
            target = prod(x**k*sum(x**j/factorial(j) for j in range(9-k))
                          for x,k in zip(intensities,used))
            assert actual == target
            tests += 1
    assert len(tuple(permutations(range(3),2))) == 6 != 3**2
    print(f'Ordered siblings: {tests} rational factorial boxes; repeated-ID control PASS')


if __name__ == '__main__':
    check_palm()
    check_siblings()
