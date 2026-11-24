"""Exact physical colour-kernel and cutoff controls."""
from fractions import Fraction as Q
from math import comb, factorial


def edges(pairs):
    g = [[] for _ in range(max(max(a, b) for a, b in pairs) + 1)]
    for a, b in pairs:
        g[a].append(b)
        g[b].append(a)
    return tuple(map(tuple, g))


def grow(g, target):
    return tuple(a + (len(g),) if j == target else a for j, a in enumerate(g)) + ((target,),)


def cut(g, tag, red):
    return Q(sum(g[j] != (tag,) for j in red), len(red))


def check_tree(old, parent, theta):
    tag = len(old)
    g = grow(grow(old, parent), tag)
    red = frozenset({tag + 1})
    total = 0

    def walk(g, red, depth, p, q, r0, b0):
        nonlocal total
        # Marginal tilt likelihood conditional on the colour word, with physical
        # uniform sender and neighbour at the selected colour.
        r, b = len(red), len(g) - len(red)
        assert (r, b) == (1 + r0, len(old) + 1 + b0)
        current = cut(g, tag, red)
        if not depth:
            total += 1
            return p, q
        pp = qq = Q(0)
        for colour, senders in enumerate([set(range(len(g))) - red, red]):
            csum = Q(0)
            for sender in senders:
                for target in g[sender]:
                    new = grow(g, target)
                    newred = red | {len(g)} if colour else red
                    w = Q(1, len(senders) * len(g[sender]))
                    csum += w * cut(new, tag, newred)
                    p1, q1 = walk(new, newred, depth - 1,
                                  p * w * len(senders) / len(g),
                                  q * w * (theta if colour else 1 - theta),
                                  r0 + colour, b0 + 1 - colour)
                    pp += p1
                    qq += q1
            if colour:
                assert csum <= current
            else:
                # Condition on selecting a blue sender other than the tag.
                silent = sum(cut(grow(g, target), tag, red) / len(g[sender])
                             for sender in senders if sender != tag for target in g[sender])
                assert silent == (b - 1) * current
                assert csum <= current + Q(1, b * r)
        assert pp == p and qq == q
        return pp, qq

    walk(g, red, 3, Q(1), Q(1), 0, 0)
    return total


def budgets():
    # Integrate the reciprocal using the exact finite truncation remainder
    # of the negative-binomial sum, bounded by the omitted probability.
    for theta in [Q(1, 4), Q(1, 2), Q(3, 4), Q(999, 1000)]:
        u = 1 - theta
        for f in range(1, 9):
            p = inv = Q(0)
            for j in range(20):
                w = comb(j + f, f) * u**(f + 1) * theta**j
                p += w
                inv += w / (j + 1)
            assert inv <= u * (1 - u**f) / (theta * f) <= inv + (1 - p) / 21
    log_upper = sum((-1)**(j+1) * Q(1, j * 2**j) for j in range(1, 8))
    assert 1 - log_upper >= Q(297, 500)
    e_upper = sum(Q(1, factorial(j)) for j in range(5)) + Q(1, 100)
    assert e_upper < Q(11, 4)
    assert Q(4, 11) * Q(254, 255) >= Q(9, 25)
    assert Q(255, 768) >= Q(83, 250)
    c = Q(83, 250) * Q(297, 500) * Q(9, 25)
    # pi > 25/8 follows also by a 128-pair lower arctan(1) sum.
    pi_lower = 4 * sum((Q(1, 4*j+1) - Q(1, 4*j+3) for j in range(128)), Q(0))
    assert pi_lower > Q(25, 8) and 800 * c*c > 4


if __name__ == '__main__':
    star = edges([(0, i) for i in range(1, 8)])
    bare = edges([(0, 1), (0, 2), (0, 3), (1, 4), (4, 5), (2, 6), (3, 7)])
    s = [sum(Q(1, len(bare[j])) for j in row) / (len(row)*(len(row)+1))
         for row in bare]
    adverse = grow(grow(bare, 0), 1)
    assert s[1] > s[0] and [j for j in adverse[1] if len(adverse[j]) == 1] == [9]
    paths = sum(check_tree(old, parent, theta) for old, parent in [(star, 0), (star, 1), (bare, 1)]
                for theta in (Q(1, 9), Q(999, 1000)))
    print(f'PASS {paths} exact physical arrow histories: stars and nonzero zero-stock adverse fixture')
    budgets()
    print('PASS 32 reciprocal controls and exact cutoff budget Y=256')
