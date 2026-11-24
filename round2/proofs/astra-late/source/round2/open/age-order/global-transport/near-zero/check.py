"""Compact exact checks of the loss and indexing budgets; no tower expansion."""
from fractions import Fraction as Q


def check():
    # Above the floor, d times the error decreases termwise.
    d = 2**10
    assert Q(1, d**5) + Q(2, d) + 1 + Q(2, d**6) + 8 < 10
    assert Q(31, 32)**32 < Q(1, 2)
    assert Q(1, 32) < Q(7, 8)*Q(5, 8)*Q(1, 16)
    assert 8*10 > 4*10+1  # family mean > 2M using e > 2
    assert Q(1, 2) <= 1-Q(1, d)  # fresh family's finite-age rate

    loss = Q(20, 2**10)  # coefficient of N^-2
    hub_num = Q('0.02968')
    assert loss < hub_num
    assert loss < min(Q(127, 2440), Q('0.022645238637'))
    # With n<N, n(n+1)<=N^2; the largest needed birth size is n+2.
    for N in range(2, 20):
        for n in range(2, N):
            assert n*(n+1) <= N*N and n+2 <= N+2
        for m in range(1,N+3):
            assert 1-Q(m,2*(N+2)) >= Q(1,2)
    print('PASS finite-late-hit and fresh-stock budgets: mass loss < 10/d')
    print(f'PASS true-amplitude transfer and base margins: share loss < ({loss})/N^2')


if __name__ == '__main__':
    check()
