"""Exact hub certificates for selected old trees and their local scores."""
from fractions import Fraction as Q
from functools import cache
from pathlib import Path
from runpy import run_path


ROOT = Path(__file__).resolve().parents[2]
BRANCH = run_path(str(ROOT / 'preplacement-transfer' / 'check.py'))
EARLY = run_path(str(ROOT / 'early' / 'check.py'))
STOP = 112
send = BRANCH['p']


@cache
def return_lower(d, k):
    if k == 32:
        return Q(0)
    # K is the newly promoted neighbour's degree before any return.
    return (1 + k*(k-1+Q(1,d))*return_lower(d,k+1))/(1+k*(k-1+Q(1,d)))


@cache
def rescue(d, l):
    """Track stock; on last loss, try the new core's first return arrow."""
    if d >= 20:
        return 1 - Q(1, d)**l
    if not l:
        return rescue(d + 1, 1) * return_lower(d,2)
    return (d * rescue(d + 1, l + 1) + rescue(d, l - 1)) / (d + 1)


@cache
def upper(d, j):
    return send(d, j, STOP) if d <= STOP else Q(1, j + 1)


def grow(a, u):
    return a[:u]+(a[u]+(len(a),),)+a[u+1:]+((u,),)


def profile(a):
    d = tuple(map(len, a))
    r = tuple(sum(Q(1, d[v]) for v in nbrs) for nbrs in a)
    w = tuple(x/(k+1) for x, k in zip(r, d))
    s = tuple(x/k for x, k in zip(w, d))
    return d, w, s


def leaf_upper(a, v):
    (p,) = a[v]
    j = sum(len(a[u]) == 1 for u in a[p])-1
    if len(a[p]) >= j+2:
        return upper(len(a[p]), j)
    return 1-EARLY['silent_bound'](j, max(32,j), 1, 1)


def leaf_lower(a, v):
    (p,) = a[v]
    return rescue(2,1)*return_lower(1,len(a[p])) if len(a[p]) <= 32 else Q(0)


def serial_lower(parent_degree):
    a = Q(1,parent_degree)
    extra = a*(rescue(3,2)-rescue(2,1))/(parent_degree+Q(3,2)+a)
    return rescue(2,1) + extra - upper(2,0)


def score(a):
    d,w,s = profile(a)
    m = len(a)
    serial = sum(wi*serial_lower(di+1) for wi,di in zip(w,d))
    loss = Q(0)
    for u, nbrs in enumerate(a):
        for v in nbrs:
            if s[u] <= s[v]:
                continue
            p = grow(grow(a,u),v)
            loss += (s[u]-s[v])*max(Q(0),leaf_upper(p,m)-leaf_lower(p,m+1))
    return (serial-loss)/sum(w), loss/sum(w)


def spider_certificate():
    # For k>=32, write h=1/k. Bounds are affine in h. Drop centre credit.
    other_tip = rescue(2,1)*return_lower(1,3)
    for h in (Q(0), Q(1,32)):
        w_arm = (1+h)/3
        s_arm = w_arm/2
        serial = w_arm*serial_lower(3)+serial_lower(2)/4
        debit = s_arm*upper(3,1)+(Q(1,4)-s_arm)*(upper(2,0)-other_tip)
        assert serial-debit > Q(6,1000)
    assert all(score(spider(k,2))[0] > 0 for k in range(2,32))


def path(n):
    return tuple(tuple(v for v in (u-1,u+1) if 0 <= v < n) for u in range(n))


def spider(k, length):
    a = ((),)
    for _ in range(k):
        p = 0
        for _ in range(length):
            a = grow(a,p)
            p = len(a)-1
    return a


def examples():
    for n in (4, 5, 8, 16, 32):
        yield f'path-{n}', path(n)
    for n in (4, 8, 16, 32):
        a = path(n)
        for u in range(n):
            a = grow(a,u)
        yield f'comb-{n}', a
    for k in (3, 4, 8, 16, 32, 64, 128, 256):
        for l in (2, 3):
            yield f'spider-{k}x{l}', spider(k,l)


def crowded(spine, teeth, length):
    a = path(spine)
    for u in range(spine):
        for _ in range(teeth):
            p = u
            for _ in range(length):
                a = grow(a,p)
                p = len(a)-1
    return a


if __name__ == '__main__':
    for d in range(1,20):
        for k in range(2,32):
            assert return_lower(d,k) >= return_lower(d,k+1)
    for d in range(2,20):
        for l in range(1,d):
            assert rescue(d,l) >= rescue(d,l-1)
            assert rescue(d+1,l+1) >= rescue(d,l)
    lo = rescue(2,1)-upper(2,0)
    print(f'PASS monotone rescue table; middle lower {float(rescue(2,1)):.9f}, serial gap lower {float(lo):.9f}')
    spider_certificate()
    print('PASS all length-two spiders for k>=2; k>=32 numerator > .006 per arm')
    for name,a in examples():
        g,b = score(a)
        assert g > 0, name
        print(f'{name:16s} signed lower / A {float(g): .9f}, edge debit / A {float(b): .9f}')
    runs = [(length, spine, t, score(crowded(spine,t,length))[0])
            for length in (1,2,3) for spine in (2,4,8) for t in (2,4,8,16)]
    worst = min(runs,key=lambda row: row[3])
    assert worst[3] > Q(42428843,10**9)
    print(f'crowded combs: {sum(x[3]>0 for x in runs)}/{len(runs)} positive; '
          f'worst (length, spine, teeth)={worst[:3]}: {float(worst[3]):.9f}')
