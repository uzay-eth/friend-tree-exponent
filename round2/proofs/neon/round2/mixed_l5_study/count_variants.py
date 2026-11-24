#!/usr/bin/env python3
"""Exact L5 row counts; h=0/1 exact apexes, one common BIG apex per bucket."""
from collections import Counter
from functools import cache
from itertools import combinations_with_replacement, product
from math import comb, prod


@cache
def size(s):
    return 1+sum(map(size,s))


shapes = [()]
for n in range(2,6):
    shapes += [s for a in range(1,n) for s in combinations_with_replacement(shapes,a)
               if sum(map(size,s)) == n-1]


@cache
def cnt(s, role, d, big):
    a = len(s)
    def exact(h):
        return sum(max(0,k-h-a+1) for k in range(2,d))
    n = exact(1) if role == 'child' else exact(0) + exact(1)
    if a == (0 if role == 'child' else 1):
        n += big
    return n*prod(comb(cnt(c,'child',d,big)+multiplicity-1,multiplicity)
                  for c,multiplicity in Counter(s).items())


def count_l4_coverage():
    """D0=12 table row/role incidences with at least one external credit.

    Counts are storage-coordinate diagnostics, not physical tree frequencies.
    In symmetric child positions, equal slots still count twice as receivers.
    """
    def endpoint(role):
        hs = (0,1) if role == 'apex' else (1,)
        exact = [(d,k,h) for h in hs for d in range(2,12)
                 for k in range(int(role=='apex'),d-h+1)]
        nt = 193 if role=='apex' else 105
        return exact + [(12,0,1)]*(nt-len(exact))
    def interior(min_k,max_k,cap,h):
        return [(d,k,h) for d in range(2,cap) for k in range(min_k,min(max_k,d-h)+1)]
    p,w = endpoint('apex'),endpoint('child')
    p4=interior(1,8,10,1)
    yt=interior(2,7,9,1)
    ybu=interior(2,8,10,0)+interior(2,8,10,1)
    s3=interior(3,6,8,0)+interior(3,6,8,1)
    # Compress each slot to flags that affect counts; canonicalize symmetric
    # groups by combinations with repetition of *actual control indices*.
    def flag(x, a, deg, h_inside):
        d,k,h = x
        return d<6, d<12 and k+h>a+h_inside, 6<=d<12, deg>1
    def group(table, count, a, deg, hi):
        if count==1:
            return Counter((flag(x,a,deg,hi),) for x in table)
        return Counter(tuple(sorted(flag(table[i],a,deg,hi) for i in ids))
                       for ids in combinations_with_replacement(range(len(table)),count))
    groups = {
        'P4': [group(p,1,1,1,0),group(p4,1,1,2,1),group(p4,1,1,2,1),group(w,1,0,1,1)],
        'Yt': [group(p,1,1,1,0),group(yt,1,2,3,1),group(w,2,0,1,1)],
        'Yb': [group(ybu,1,2,2,0),group(p4,1,1,2,1),group(w,1,0,1,1),group(w,1,0,1,1)],
        'S3': [group(s3,1,3,3,0),group(w,3,0,1,1)],
    }
    totals = {False: Counter(), True: Counter()}
    assert [prod(sum(x.values()) for x in parts) for parts in groups.values()] == [
        26263440,22554945,25004700,4763640]
    for parts in groups.values():
        for terms in product(*(x.items() for x in parts)):
            n = prod(mult for _,mult in terms)
            slots = sum((flags for flags,_ in terms), ())
            if all(f[0] for f in slots):
                continue
            big_interior=any(it and not small for small,_,_,it in slots)
            for small, ext, midrange, it in slots:
                if not ext:
                    continue
                counts=totals[it]
                counts['total'] += n
                if midrange:
                    counts['large_receiver'] += n
                elif big_interior:
                    counts['large_spectator_interior'] += n
                else:
                    counts['supported'] += n
    for kind, c in [('A',totals[True]),('B',totals[False])]:
        assert c['total']==c['large_receiver']+c['large_spectator_interior']+c['supported']
        print(f'{kind} mixed row/receiver incidences: {dict(c)}')


if __name__ == '__main__':
    for d,big in ((6,1),(7,1),(6,2),(7,2)):
        rows = sum(cnt(s,'apex',d,big) for s in shapes if size(s) == 5)
        qbytes = 8*rows
        # Same dense indexing metadata, cache and backing fiber storage as
        # coarse=6 build. Even D=7 with two buckets fits the existing ns<64.
        print(f'Dc=D5={d} BIG buckets={big} entries={rows} qbytes={qbytes} '
              f'projected_payload_D0=12:{qbytes+8958228} D0=20:{qbytes+8961876}')
    count_l4_coverage()
