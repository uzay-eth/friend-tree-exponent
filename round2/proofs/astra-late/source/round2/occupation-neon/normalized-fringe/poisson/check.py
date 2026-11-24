#!/usr/bin/env python3
"""Exact first-event Poisson fit and a null direction in the fitted flow tests."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
BASE = HERE.parent.parent
sys.path[:0] = [str(BASE/'ratio-transfer'), str(BASE/'check')]
from crowns import J, degree, gain, moves, n
from replay import Cert, Pair, embedding, shape, tree

S = Q('.56572365')  # reference trial, not an estimate or an assumption on mu
T = Q(1)
A, B = Q('.565673'), Q('.565775')
# Standalone GMP potential ceilings; no copied all-stock floor literals.
C0P, C0Q = Q('1.84124847315'), Q('1.95820675182')


@cache
def bracket(z):
    return tree(z[0], *(bracket(y) for y in z[1]))


def birth(z):
    return Q(z[0], degree(z))


def physical_moves(z):
    """Same planted hazards as crowns.moves, aggregated without listing leaves."""
    p,x = shape(bracket(z),1)
    kids = [tuple(i for i in range(1,len(p)) if p[i]==j) for j in range(len(p))]

    def after(j,promotion):
        def go(i):
            children = [go(k) for k in kids[i]]
            if i==j and promotion:
                children.append(J)
            return x[i][0]+((-1 if promotion else 1) if i==j else 0),tuple(sorted(children))
        return go(0)

    out = defaultdict(Q)
    for i,(l,k,h) in enumerate(x):
        donors = kids[i]+((p[i],) if i else ())
        r = l+sum(Q(1,sum(x[j])) for j in donors)
        if r:
            out[after(i,False)] += r
        if l:
            out[after(i,True)] += Q(l,l+k+h)
    return dict(out)


@cache
def row(z, s):
    """A-s: intrinsic arrows, exterior child transfer and newborn score."""
    out = defaultdict(Q, {z: -s})
    for y, r in physical_moves(z).items():
        out[z] -= r
        out[y] += r
    for y in z[1]:
        out[gain(y)] += Q(1, degree(z))
        out[y] -= Q(1, degree(z))
    out[J] += birth(z)
    assert sum(out.values()) == birth(z)-s
    return {y: r for y, r in out.items() if r}


def dot(c, h):
    return sum(r*h(y) for y,r in c.items())


def whole_row_fixture(z):
    """Compare T(A-s) to (L0-s)T as coefficient vectors, not just on the fit."""
    def nodes(y):
        yield y
        for w in y[1]:
            yield from nodes(w)

    projected, direct = defaultdict(Q), defaultdict(Q)
    for x in nodes(z):
        for y, r in row(x,S).items():
            projected[y] += r
        direct[x] -= S
    for nxt, r in moves(z):
        for y in nodes(nxt):
            direct[y] += r
        for y in nodes(z):
            direct[y] -= r
    assert {y:r for y,r in projected.items() if r} == {y:r for y,r in direct.items() if r}


def rref(m, cols):
    m = [list(map(Q, v)) for v in m]
    pivots = []
    for col in range(cols):
        j = len(pivots)
        pivot = next((i for i in range(j,len(m)) if m[i][col]), None)
        if pivot is None:
            continue
        m[j], m[pivot] = m[pivot], m[j]
        m[j] = [x/m[j][col] for x in m[j]]
        for i in range(len(m)):
            if i != j and m[i][col]:
                v = m[i][col]
                m[i] = [x-v*y for x,y in zip(m[i],m[j])]
        pivots.append(col)
    return m, pivots


def pilot(pair, cap):
    @cache
    def apex(z):
        p, x = shape(bracket(z), 1)
        v = [Q(0), Q(0)]
        for kind, ids in embedding(p):
            if ids[0] == 0:
                for i,w in enumerate(pair.weight(kind,tuple(x[j] for j in ids))):
                    v[i] += w
        return tuple(v)

    def h(z):
        p,q = apex(z)
        return p+T*q

    def source(z,s=S):
        return dot(row(z,s),h)

    layer = {J}
    states = []
    for k in range(2,cap+1):
        assert all(n(z)==k for z in layer)
        states.extend(sorted(layer))
        layer = {y for z in layer for y,r in moves(z) if r}
    # Zero-virgin fit: include one birth-flow score to make the square system.
    basis = states[1:]
    assert states[0] == J
    m, pivots = rref([[row(z,S).get(y,0) for y in basis]+[birth(z)-S,source(z)]
                     for z in states],len(states))
    assert pivots == list(range(len(states)))
    u = dict(zip(basis,[v[-1] for v in m[:-1]]))
    beta = m[-1][-1]
    assert u.get(J,0) == 0

    def correction(z,s=S):
        return dot(row(z,s),lambda y:u.get(y,0))+beta*(birth(z)-s)

    assert all(source(z)==correction(z) for z in states)
    for z in states+sorted(layer):
        assert physical_moves(z)==dict(moves(z))
        whole_row_fixture(z)

    # An a priori stationary bound on the fitted correction's complement.
    # Outside the cap only C*Delta can see nonzero finite fringe values.
    # The sum of child rates is <=1. All arithmetic is absolute and outward.
    jumps = [u.get(gain(z),0)-u.get(z,0) for z in states]
    jump_bound = max(map(abs,jumps))
    boundary = max(abs(source(z)-correction(z)) for z in layer)
    u_cap = max(map(abs,u.values()))
    # Use the same score at the true discount; its stationary mean is zero.
    # On the matched set, the residual then consists only of the mismatch
    # of the trial and true discounts. C0 bounds the stationary mean of H,
    # not its apex-wise pointwise maximum.
    converted = max(S-A,B-S)*(C0P+T*C0Q+u_cap+abs(beta))
    print('finite matched / next layer',len(states),len(layer),'birth coefficient',float(beta),
          'max next-layer defect',float(boundary),'certified outside transfer cap',float(jump_bound))
    print('stationary discount-conversion cap',float(converted))
    # Uniform in the incoming [A,B], the old nonnegative full charge totals
    # bound the charge on any unpaid subset (for the production allocation).
    # The fixed-trial row minus the actual-discount row is (mu-S)*(H-U-beta).
    print('all-slack limits from old ceilings over [a,b]',float((B-A)*C0P),float((B-A)*C0Q),
          'D-tE interval at conditional mu=s',float(-T*(B-S)*C0Q),float((S-A)*C0P))

    # Same birth rate (zero); no outgoing terms see *any* fitted test, at any s.
    # Hold both potential averages exactly fixed, not just their chosen mixture.
    stocks = (9,12,16,17,18,19,23,31,63,13999,14000,14001,14002,20000,10**6)
    def null(k, nested):
        witnesses = [(0,((stock,((18,()),) if nested else ()),)*k) for stock in stocks]
        columns = [[Q(1), Q(1,degree(z[1][0])), *apex(z)] for z in witnesses]
        m, pivots = rref(list(zip(*columns)),len(witnesses))
        directions = []
        for j in range(len(witnesses)):
            if j in pivots:
                continue
            v = [Q(0)]*len(witnesses)
            v[j] = 1
            for i,l in enumerate(pivots):
                v[l] = -m[i][j]
            d = sum(w for w in v if w>0)
            directions.append([w/d for w in v])
        return max((abs(sum(w*source(z) for w,z in zip(v,witnesses))),k,nested,witnesses,v)
                   for v in directions)

    _, k, nested, witnesses, v = max(null(k,nested) for k in (1,2,3,4) for nested in (False,True))
    assert sum(v)==0 and sum(w for w in v if w>0)==1
    for z in witnesses:
        c = row(z,S)
        assert birth(z)==0 and n(z)>cap and c.keys().isdisjoint(states)
        assert degree(z)==k+1 and all(degree(child)>cap for child in z[1])
        assert {degree(y) for y in physical_moves(z)}=={k+1,k+2}
        # The root shell generator agrees in mass after averaging: only
        # an intrinsic gain is possible, at k / D_child.
        shell = lambda y: Q(degree(y)==k+1)
        assert dot(c,shell)==-S-Q(k,degree(z[1][0]))
    for i in (0,1):
        assert sum(w*apex(z)[i] for w,z in zip(v,witnesses))==0
    assert sum(w/degree(z[1][0]) for w,z in zip(v,witnesses))==0
    diff = sum(w*(source(z)-correction(z)) for w,z in zip(v,witnesses))
    for s in (A,B):
        assert sum(w*(source(z,s)-correction(z,s)) for w,z in zip(v,witnesses))==diff
    mean_h = sum(w*h(z) for w,z in zip(v,witnesses) if w>0)
    assert mean_h == -sum(w*h(z) for w,z in zip(v,witnesses) if w<0)
    assert mean_h>0 and abs(diff)/mean_h > Q('1e-8')
    print('flow-null normalized gap',float(abs(diff)/mean_h),'signed per-core gap',float(diff))
    print('null direction children / nested',k,nested,[(z[1][0][0],str(w)) for z,w in zip(witnesses,v) if w])


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('payloads',nargs=2,type=Path)
    parser.add_argument('--cap',type=int,default=6,choices=range(4,9))
    args = parser.parse_args()
    pilot(Pair(*(Cert(f) for f in args.payloads)),args.cap)
