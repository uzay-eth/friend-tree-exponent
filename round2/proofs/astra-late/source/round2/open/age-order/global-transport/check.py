"""Rational controls for the block estimate and the first-promotion charge."""
from fractions import Fraction as Q
from math import comb


def compositions(total, parts):
    if parts == 1:
        yield (total,)
    else:
        for first in range(total + 1):
            for tail in compositions(total - first, parts - 1):
                yield (first,) + tail


def blocks():
    count = 0
    for f in range(1, 6):
        for s in range(11):
            words = list(compositions(s, f + 1))
            assert len(words) == comb(s + f, f)
            for j in range(s):
                assert sum(sum(word[:-1]) <= j for word in words) == comb(j + f, f)
            h = sum((sum((Q(1, 1 + j) for j in range(s - w[-1], s)), Q(0))
                     for w in words), Q(0)) / len(words)
            assert h == (1 - Q(1, len(words))) / f <= Q(1, f)
            count += 1
    return count


def charges():
    count = 0
    for m in (4, 9, 40):
        for d in (2, 3):
            for cap in (Q(1, 3), Q(1)):
                for k in range(1, 7):
                    survival = Q(1)
                    cumulative = [Q(0)]
                    for i in range(k):
                        survival *= 1 - Q(1, (m + i) * d)
                        cumulative.append(min(cap, 1 - survival))
                    q = [Q(3, 4)**(i + 1) for i in range(k + 1)]
                    charge = sum((cumulative[i+1] - cumulative[i]) * q[i]
                                 for i in range(k)) + (cap - cumulative[-1]) * q[-1]
                    assert charge == cap*q[-1] + sum(cumulative[i+1] * (q[i] - q[i+1])
                                                      for i in range(k))
                    assert charge < cap*q[0]
                    count += 1
    return count


if __name__ == '__main__':
    print(f'PASS {blocks()} exact conditional red-block harmonic bounds')
    print(f'PASS {charges()} exact early promotion charge bounds')
