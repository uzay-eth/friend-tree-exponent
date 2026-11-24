"""Independent bounded replay of the 141 band and all-label splice."""

from fractions import Fraction as Q
from itertools import product
from math import factorial


def weights(parents):
    nbr = [set() for _ in range(len(parents)+1)]
    for v, u in enumerate(parents, 1):
        assert 0 <= u < v
        nbr[v].add(u)
        nbr[u].add(v)
    deg = list(map(len, nbr))
    s = [sum(Q(1, deg[v]) for v in row)/Q(deg[u]*(deg[u]+1))
         for u, row in enumerate(nbr)]
    A = sum(d*q for d, q in zip(deg, s))
    paired = hi = Q(0)
    for v, u in enumerate(parents, 1):
        paired += abs(s[u]-s[v])
        if s[u] != s[v] and deg[u if s[u] > s[v] else v] == 1:
            hi += abs(s[u]-s[v])
    assert 3*hi <= A and paired <= A
    return sorted(deg), A, hi


def reachable():
    count = 0
    for n in range(2, 8):
        for parents in product(*(range(v) for v in range(1, n))):
            weights(parents)
            count += 1
    da, aa, xa = weights((0, 0, 0, 1, 2, 3))
    db, ab, xb = weights((0, 0, 0, 1, 4, 5))
    assert da == db == [1, 1, 1, 2, 2, 2, 3]
    assert (aa, xa, ab, xb) == (Q(59,24), Q(1,12), Q(167,72), 0)
    print(f'PASS {count} reachable parent arrays; equal-degree fixtures have X/A=2/59 and 0')


def budgets():
    t = Q(5,6)
    # Reconstruct the coefficient from the protected-prefix bounds.
    log_hi = sum(((-1)**(j+1)*Q(5,11)**j/j for j in range(1, 12)), Q(0))
    exp_lo = sum(((-t)**j/factorial(j) for j in range(12)), Q(0))

    def c(b):
        return (b-1)/(Q(16,5)*b+1)*(1-log_hi/t)*exp_lo*(1-t*t/(2*(b-t)))

    gamma_lo = sum(((-1)**j*Q(3,3*j+1)/factorial(j) for j in range(6)), Q(0))
    assert gamma_lo == Q(141077,58240) > Q(24223,10000)
    assert Q(9,5)*Q(242,241) < Q(3,4)*gamma_lo
    assert Q(9,5)*Q(42,41) < Q(4,5)*gamma_lo
    for b, x, root, k in ((141, Q(1), Q(57,25), Q(3,4)),
                          (20, Q(10,11), Q(3,2), Q(4,5))):
        y = b*x
        assert root**6 < y
        cost = Q(4,3)+Q(2,3)*k/root
        # In the scaled test drop sqrt(1-x); x=1 is only an algebraic corner.
        attenuation = 1 if b == 141 else 1-x
        assert c(b)**2*Q(157,50)*y > cost**2*attenuation
    assert Q(4,3) < 2
    print('PASS gamma constants and strict rational corner margins (141, and b=20 at 10/11)')


def splice():
    founder = {1}
    early = set(range(2, 20))
    bridge = set(range(20, 501))
    assert not (founder & early or early & bridge or founder & bridge)
    assert founder | early | bridge == set(range(1, 501))
    assert (min(early)+1, max(early)+1, min(bridge)+1) == (3, 20, 21)
    print('PASS boundary ledger: 1=2; 2→3 through 19→20; 20→21 onward')


if __name__ == '__main__':
    reachable()
    budgets()
    splice()
