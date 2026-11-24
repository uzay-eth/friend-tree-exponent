#!/usr/bin/env python3
"""Exact flat physical arrows versus permanent marked-depth rows."""
from collections import Counter, defaultdict
from fractions import Fraction as F
from functools import cache


Z = F(7, 5)


@cache
def children(p):
    return tuple(tuple(v for v, u in enumerate(p) if u == x) for x in range(len(p)))


def degrees(p, slot):
    return tuple(len(ch) + (slot if not x else 1) for x, ch in enumerate(children(p)))


def physical(p, a, slot):
    """Each individual source with one channel per real adjacent vertex."""
    D = degrees(p, slot)
    for x, ch in enumerate(children(p)):
        for y in (*ch, *((p[x],) if x else ())):
            b = a
            if D[y] == 1:
                # Promote a labelled leaf; its source is necessarily its parent.
                assert p[y] == x
                b = (*a[:y], F(1, min(D[x], 20)), *a[y+1:])
            yield p+(y,), b+(None,), F(1, D[x])


@cache
def fringe(p, a, x):
    ps, marks = [-1], [a[x]]

    def visit(v, i):
        for y in children(p)[v]:
            ps.append(i)
            marks.append(a[y])
            visit(y, len(ps)-1)

    visit(x, 0)
    return tuple(ps), tuple(marks)


def score(p, a):
    # Isomorphism-invariant, with shape and all marks, including the top.
    return F(len(p)**2 + sum(len(ch)**2 for ch in children(p)), 7) + sum(
        (depth(p,x)+1)**2 * m for x, m in enumerate(a) if m is not None)


def depth(p, x):
    return 0 if not x else 1+depth(p,p[x])


def copies_and_births(p, a, x, slot, mutant):
    D = degrees(p,slot)
    d = D[x]+(mutant == 'free-extra-slot' and x == 0 and slot == 0)
    p0 = F(1,d)
    for y in children(p)[x]:
        if D[y] > 1:
            q0 = a[y] if mutant != 'current-child-kill-cancel' else F(1,min(d,20))
            py, ay = fringe(p,a,y)
            yield (Z if mutant != 'no-copy-z' else 1)*(q0-p0), py, ay
            yield (Z if mutant != 'no-copy-z' else 1)*p0, py+(0,), ay+(None,)
        else:
            yield (Z if mutant != 'no-birth-z' else 1)*p0, (-1,0), (F(1,min(d,20)),None)


def compare(p,a,slot,mutant=None):
    def T(q,m):
        return sum(Z**depth(q,x)*score(*fringe(q,m,x)) for x,d in enumerate(degrees(q,slot))
                   if d>1 and (x or slot))

    # For a free centre only the nonfounder cores appear in T.
    base = T(p,a)
    lhs = sum(r*(T(q,b)-base) for q,b,r in physical(p,a,slot))
    rhs = F(0)
    for x,d in enumerate(degrees(p,slot)):
        if d<=1:
            continue
        val = F(0)
        if x or slot:
            q,b = fringe(p,a,x)
            val += sum(r*(score(n,c)-score(q,b))*(Z if mutant=='intrinsic-z' else 1)
                       for n,c,r in physical(q,b,1))
            val -= (F(1,2) if mutant=='kill-inert' and not x else b[0])*score(q,b)
        val += sum(r*score(q,b) for r,q,b in copies_and_births(p,a,x,slot,mutant))
        rhs += Z**depth(p,x)*val
    return lhs==rhs


def main():
    mutations=('no-copy-z','no-birth-z','intrinsic-z','kill-inert',
               'current-child-kill-cancel','free-extra-slot')
    misses=Counter()
    for slot, seed in ((1,(-1,0)),(0,(-1,0,0))):
        reached={(seed,(F(0),)+(None,)*(len(seed)-1))}
        layer=reached.copy()
        for _ in range(4):
            layer={(n,c) for p,a in layer for n,c,r in physical(p,a,slot)}
            reached |= layer
        # Cap seam and children whose parent's degree changed since they were born.
        p=(-1,0,1,0)+(0,)*20
        a=(F(0),F(1,2))+(None,)*(len(p)-2)
        reached.add((p,a))
        for p,a in reached:
            assert compare(p,a,slot)
            for m in mutations:
                misses[m] += not compare(p,a,slot,m)
        print(f'PASS {len(reached)} labelled states, slot={slot}')
    assert all(misses[m]>0 for m in mutations)
    print('Rejected mutations:', dict(misses))


if __name__=='__main__':
    main()
