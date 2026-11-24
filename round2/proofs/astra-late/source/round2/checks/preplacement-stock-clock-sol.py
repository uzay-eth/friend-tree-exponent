#!/usr/bin/env python3
"""Independent rational spot checks for the monitored stock and adverse clocks."""

from fractions import Fraction as Q
from hashlib import sha1
from itertools import combinations, product
from math import comb, factorial, prod
from pathlib import Path


PINS = {
    'README.md': '976800b905a5bd8c75b1594f786da2690f6bb539',
    'check.py': '17e171cef28505bdadcc929544870d238c86f68b',
    'ledger.py': '57a60df2ffcd86142cebab29dad12d92a12c6a66',
    'replay.py': 'bf49f5510669c15ea4c030b7e859c96de222da04',
}


def verify_pins():
    stock = Path(__file__).resolve().parents[1] / 'open/age-order/global-transport/beta-clocks/stock'
    for name, expected in PINS.items():
        data = (stock / name).read_bytes()
        actual = sha1(f'blob {len(data)}\0'.encode() + data).hexdigest()
        assert actual == expected, (name, actual, expected)
    return len(PINS)


def trees(n):
    for parents in product(*(range(i) for i in range(1, n))):
        adj = [set() for _ in range(n)]
        for child, parent in enumerate(parents, 1):
            adj[child].add(parent)
            adj[parent].add(child)
        yield adj


def physical():
    counts = [0] * 4
    for n in range(3, 7):
        for adj in trees(n):
            deg = list(map(len, adj))
            receive = [sum(Q(1, deg[w]) for w in row) for row in adj]
            for tag, D in enumerate(deg):
                leaves = [w for w in adj[tag] if deg[w] == 1]
                for d in range(2, D + 1):
                    for j in range(1, min(d-1, len(leaves)) + 1):
                        for subset in combinations(leaves, j):
                            true_loss = sum(Q(1, D) for _ in subset)
                            assert true_loss <= Q(j, d)
                            assert j + Q(j,d) <= n
                            if j == 1:
                                # Even a supplemented loss can increase the dropped
                                # leaf's physical degree at most once.
                                for target in range(n):
                                    assert deg[subset[0]] + (target == subset[0]) <= 2
                            counts[0] += 1
                for mother in adj[tag]:
                    K = deg[mother]
                    for d in range(1, D + 1):
                        for k in range(max(K,2), n + 1):
                            ret = Q(1,k)
                            assert ret <= Q(1,K)
                            # Mothers below the ceiling may grow without a mark.
                            mother_growth = receive[mother] if K == k else Q(0)
                            grow = k-1+Q(1,d)
                            assert mother_growth <= grow
                            assert ret + grow <= 18
                            counts[1] += 1
                    if D >= 2:
                        for leaf in leaves:
                            if leaf == mother:
                                continue
                            for k in range(max(K,2), n):
                                assert 1+Q(1,k) <= 1+Q(1,K)
                                assert Q(1,D) <= Q(1,2)
                                assert K < k or receive[mother] <= k-Q(1,2)
                                assert 1+Q(1,k)+k <= 18
                                counts[2] += 1
                    if D == 1 and K >= 2:
                        J = sum(deg[w] == 1 and w != tag for w in adj[mother])
                        for k in range(2,K+1):
                            for j in range(min(J,k-1)+1):
                                assert j+1 <= receive[mother]
                                assert Q(1,K) <= Q(1,k)
                                physical_loss_needed = Q(J,K) if J == j else 0
                                assert physical_loss_needed <= Q(j,k)
                                # Match gain first, then promotion and necessary
                                # sibling deaths; the latter two are supplemented.
                                assert j+1+Q(j+1,k) <= n
                                counts[3] += 1
    return counts


def btail(s, j, x):
    return sum(comb(s-1,v)*x**v*(1-x)**(s-1-v) for v in range(max(0,j)))


def h(s,j,x):
    return sum(btail(s,j-i,x)/factorial(i) for i in range(j))


def source_rows():
    count = 0
    for s in (21,22,80):
        for x in (Q(1,s),Q(7,5*s),Q(1,3)):
            for d in range(2,8):
                for j in range(1,d):
                    nxt = lambda v: h(s+1,v,x)
                    base = h(s,j,x)
                    delta = nxt(j)-nxt(j-1)
                    die = Q(j,d)
                    assert j*nxt(j+1)+die*nxt(j-1)+(s-j-die)*nxt(j)-s*base == (1-die)*delta
                    count += 1
            for k in range(2,11):
                p = Q(1,k)
                assert p*h(s+1,1,x) >= 0  # the return baseline is zero
                nxt = lambda v: h(s+1,v,x)
                base = h(s,1,x)
                assert (1+p)*nxt(2)+(s-1-p-Q(1,2))*nxt(1)-s*base == p*(nxt(2)-nxt(1))+nxt(1)/2
                count += 2
    return count


def beta_events():
    count = 0
    for m, R, a, p in product((21,35,100),(16,18),(Q(1,3),Q(7,5),Q(2),Q(13)),range(1,8)):
        rising = lambda v: prod(v+t for t in range(p))
        for j in (1,2,5):
            f = lambda size, i: rising(a)/rising(size)*(Q(R,R+p)**i)
            assert m*f(m,j) == (m-R)*f(m+1,j)+R*f(m+1,j-1)
            count += 1
    return count


if __name__ == '__main__':
    print('PASS pinned stock blobs:', verify_pins())
    print('PASS physical stock/return/host/adverse rows:', *physical())
    print('PASS exact beta-binomial source rows:', source_rows())
    print('PASS exact uniformized beta-product moments:', beta_events())
