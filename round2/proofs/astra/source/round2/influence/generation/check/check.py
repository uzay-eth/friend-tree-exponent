#!/usr/bin/env python3
"""Exact physical-generation jets: labelled physical tree vs positive fringe."""

import argparse
from collections import defaultdict
from fractions import Fraction as F
from functools import cache

J = (F(0), 1, ())
MAX_DEPTH = 3


def need(ok, why):
    if not ok:
        raise ValueError(why)


def fringe_degree(x):
    return 1+x[1]+len(x[2])


def with_leaf(x):
    return x[0],x[1]+1,x[2]


@cache
def intrinsic(x,cap):
    """Fringe incidence recursion, independent of the labelled physical engine."""
    a,l,children = x
    p = F(1,fringe_degree(x))
    out = defaultdict(F)
    out[with_leaf(x)] += l+sum(F(1,fringe_degree(y)) for y in children)
    if l:
        newborn = (max(p,F(1,cap)) if cap else p,1,())
        out[a,l-1,tuple(sorted(children+(newborn,)))] += l*p
    for i,y in enumerate(children):
        for nxt,r in (*intrinsic(y,cap),(with_leaf(y),p)):
            out[a,l,tuple(sorted(children[:i]+(nxt,)+children[i+1:]))] += r
    return tuple((y,r) for y,r in out.items() if r)


def positive_jet(seed, steps, mode, cap=None):
    """Coefficients of s times occupation in x=1/s, i.e. generator powers."""
    flow = {(0,seed): F(1)}
    result = []
    for n in range(steps+1):
        col = [F(0)]*(MAX_DEPTH+1)
        birth = [F(0)]*(MAX_DEPTH+1)
        nxt = defaultdict(F)
        for (k,x), mass in flow.items():
            a,l,children = x
            p = F(1,fringe_degree(x))
            col[k] += mass
            old = intrinsic(x,cap)
            nxt[k,x] -= mass*(sum(r for _,r in old)+a)
            for y,r in old:
                nxt[k,y] += mass*r
            if k < MAX_DEPTH:
                birth[k+1] += mass*l*p
                mark = max(p,F(1,cap)) if cap else p
                nxt[k+1,(mark,1,())] += mass*l*p
            for y in children:
                need(y[0] >= p,'invalid birth mark')
                d = k if mode == 'old-depth' else k+1
                if d <= MAX_DEPTH:
                    if mode != 'omit-plain':
                        nxt[d,y] += mass*(y[0]-p)
                    nxt[d,with_leaf(y)] += mass*p
        result.append((tuple(col),tuple(birth)))
        flow = {x:v for x,v in nxt.items() if v}
    return result


def literal_seed():
    # Undirected vertex IDs: root 0, one real leaf 1, silent slot is external.
    # Tag -1 means a leaf; the ordinary ghost-root core has generation 0.
    return ((1,), (0,)), (0,-1)


def physical_step(tree):
    adj,tags = tree
    n = len(adj)
    for sender,ns in enumerate(adj):
        r = F(1,len(ns)+(sender == 0))
        for target in ns:
            child = list(adj)
            child[target] = child[target]+(n,)
            child.append((target,))
            new = list(tags)+[-1]
            generation = None
            if tags[target] < 0:
                need(tags[sender] >= 0,'leaf selected by noncore')
                new[target] = tags[sender]+1
                generation = new[target]
            yield (tuple(child),tuple(new)),r,generation


def physical_jet(steps):
    flow = {literal_seed():F(1)}
    result = []
    for n in range(steps+1):
        col,birth = ([F(0)]*(MAX_DEPTH+1) for _ in range(2))
        nxt = defaultdict(F)
        for state,mass in flow.items():
            adj,tags = state
            for k in tags:
                if 0 <= k <= MAX_DEPTH:
                    col[k] += mass
            diagonal = F(len(adj))-F(1,1+len(adj[0]))
            nxt[state] -= mass*diagonal
            total = F(0)
            for y,r,g in physical_step(state):
                total += r
                if g is not None and 0 < g <= MAX_DEPTH:
                    birth[g] += mass*r
                if n < steps:
                    nxt[y] += mass*r
            need(total == diagonal,'physical rate mass')
        result.append((tuple(col),tuple(birth)))
        flow = {x:v for x,v in nxt.items() if v}
    return result


def min_path(steps,h=3):
    """Positive (unlimited-in-size) lower bounds via first n real jumps."""
    s = F('0.56572365')
    flow = {literal_seed():F(1)}
    ans = [F(1)]+[F(0)]*h
    for i in range(steps):
        nxt = defaultdict(F)
        for state,mass in flow.items():
            adj,_ = state
            weight = mass/(s+len(adj)-F(1,len(adj[0])+1))
            for child,r,g in physical_step(state):
                if g is not None and g<=h:
                    ans[g] += weight*r
                if i+1<steps:
                    nxt[child] += weight*r
        flow = nxt
    return ans


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--steps',type=int,default=5)
    p.add_argument('--plant',choices=('old-depth','omit-plain','root-mark'))
    args = p.parse_args()
    need(4 <= args.steps <= 6,'fixture scope is 4–6 jumps')
    physical = physical_jet(args.steps)
    positive = positive_jet((F(1,2),1,()) if args.plant=='root-mark' else J,args.steps,args.plant)
    for n,((x,b),(y,c)) in enumerate(zip(physical,positive)):
        need(x==y,f'occupation mismatch derivative {n}: {x} != {y}')
        need(b==c,f'direct-birth flux mismatch derivative {n}: {b} != {c}')
        if n:
            old = physical[n-1][1]
            need(x==old,'occupation derivative is birth flux')
        print('order',n,'occupancy',*map(str,x),'flux',*map(str,b))
    for cap in (2,20):
        need(positive_jet(J,args.steps,None,cap) == physical,f'cap-{cap} mark/true degree separation')
    print('Capped birth marks 20 and stress cap 2: exact match; actual physical rates unchanged')
    lower = min_path(args.steps)
    print('discounted lower prefixes (rounded down to 10^-10)',
          *[f'{(x*10**10//1)//10**10}.{(x*10**10//1)%10**10:010d}' for x in lower])
    print('Ratios of stopped lower prefixes (diagnostic, not lower bounds):',
          *[f'{float(lower[j+1]/lower[j]):.8f}' for j in range(3)])
    print('PASS: literal tagged physical jets and positive marked-fringe transfer jets, depths 0..3')


if __name__ == '__main__':
    try:
        main()
    except ValueError as exc:
        raise SystemExit(f'REJECT: {exc}') from exc
