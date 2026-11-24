"""Exact threshold harmonic, physical one-birth checks, and the six-unit cutoff."""
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import comb, factorial

LEVELS = (Q(1, 20), Q(2, 9), Q(3, 13), Q(1, 3), Q(4, 5), Q(97, 100))


@cache
def beta_tail(size, red, x):
    if red <= 0:
        return Q(0)
    if red == size:
        return Q(1)
    return sum(comb(size - 1, i) * x**i * (1 - x)**(size - 1 - i)
               for i in range(red))


@cache
def stock_score(size, stock, x):
    # e times the probability. The common factor cancels in one-birth comparisons.
    return sum((beta_tail(size, stock - k, x) / factorial(k)
                for k in range(stock)), Q(0))


def check_harmonics():
    count = 0
    for n in range(2, 30):
        for j in range(1, n):
            for x in LEVELS:
                gain = stock_score(n + 1, j + 1, x)
                hold = stock_score(n + 1, j, x)
                loss = stock_score(n + 1, j - 1, x)
                assert gain >= hold >= loss
                assert n*stock_score(n, j, x) == j*gain + (n-j-1)*hold + loss
                count += 1
        for x in LEVELS:
            assert stock_score(n, 1, x) == (1 - x)**(n-1)
    # Omitting the k=1 Poisson term breaks the j=1 harmonic equation.
    n, x = 7, Q(1, 3)
    wrong = beta_tail(n+1, 2, x) + (n-2)*beta_tail(n+1, 1, x)
    assert wrong != n * beta_tail(n, 1, x)
    return count


def check_physical():
    cases, events, tight = 0, 0, 0
    for n in range(2, 8):
        for parents in product(*(range(v) for v in range(1, n))):
            adj = [set() for _ in range(n)]
            for v, parent in enumerate(parents, start=1):
                adj[v].add(parent)
                adj[parent].add(v)
            arrows = [(v, w, Q(1, len(adj[v]))) for v in range(n) for w in adj[v]]
            for tag, neighbours in enumerate(adj):
                leaves = sorted(v for v in neighbours if len(adj[v]) == 1)
                for j in range(len(leaves) + 1):
                    tracked = set(leaves[:j])
                    weights = {a: Q(0) for a in (-1, 0, 1)}
                    for sender, receiver, weight in arrows:
                        change = int(receiver == tag) - int(sender == tag and receiver in tracked)
                        weights[change] += weight
                    assert sum(weights.values()) == n
                    assert weights[1] >= j and weights[-1] == Q(j, len(neighbours)) <= 1
                    for x in (Q(3, 13), Q(4, 5)):
                        actual = sum(v * stock_score(n+1, j+d, x)
                                     for d, v in weights.items() if v)
                        assert actual >= n * stock_score(n, j, x)
                        if j and weights[1] == j and weights[-1] == 1:
                            assert actual == n * stock_score(n, j, x)
                            tight += 1
                    cases += 1
                    events += len(arrows)
    assert tight > 0
    return cases, events, tight


def cutoff():
    # The lower odd Taylor polynomial of exp(-u^3), integrated with t=u^3.
    a = Q(3, 2)
    gamma_third_lower = 3 * sum(Q((-1)**k, factorial(k)*(3*k+1)) * a**(3*k+1)
                                for k in range(8))
    assert gamma_third_lower > Q(13, 5)
    # Machin's formula: arctan(1/5)>1/5-1/(3*125), arctan(1/239)<1/239.
    pi_lower = 16 * (Q(1, 5) - Q(1, 375)) - Q(4, 239)
    assert pi_lower > Q(25, 8)
    assert 6 * pi_lower * Q(231, 1000)**2 > 1
    assert 36 * Q(3031, 10000)**3 > 1
    serial = sum(Q((-1)**k, factorial(k)) for k in range(10))
    assert serial > Q(36787, 100000)

    caps = []
    for s, prefactor in ((Q(1, 2), Q(231, 1000)),
                         (Q(2, 3), Q(3031, 10000) * Q(5, 13))):
        d = s * (1 - (1+s)/6) * Q(7, 9)
        assert s > (1+s)*d/6
        caps.append(prefactor * (1 - d/6))
    debit = Q(4, 3)*caps[0] + Q(2, 3)*caps[1]
    assert debit == Q(8322209, 22744800)
    assert Q(36787, 100000) - debit > Q(19, 10000)
    assert 20 * Q(3, 13) / (1-Q(3, 13)) == 6
    assert 21 * Q(2, 9) / (1-Q(2, 9)) == 6
    return debit


if __name__ == '__main__':
    rows = check_harmonics()
    cases, events, tight = check_physical()
    debit = cutoff()
    print(f'PASS {rows} exact harmonic rows, {cases} tracked physical states, '
          f'{events} directed event checks, {tight} saturated-star equalities')
    print(f'PASS normalized debit <= {debit}; normalized surplus > 19/10000 '
          'whenever (n+1)x/(1-x) >= 6')
