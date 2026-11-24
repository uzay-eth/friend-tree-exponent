#!/usr/bin/env python3
"""Reconstruct interpolated occupation scores and test them against individual rings."""
import argparse
from fractions import Fraction as Q
from functools import cache
from itertools import product
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"occupation"))
from check_five import DISCOUNT, geometry, subsets, transitions

D = [[36],[16,15],[14,10,8],[12,8,5,5]]
C = [[10],[7,6],[5,4,4],[4,4,3,3]]


def argument(xs):
    return ",".join(str(v) for x in xs for v in x)


def audit(filename):
    rows = {}
    with filename.open() as stream:
        for line in stream:
            shape, *data = line.split()
            *states, score, source = map(int,data)
            key = shape,tuple(zip(states[::3],states[1::3],states[2::3]))
            assert key not in rows and 0 <= score <= 2**40 and source >= 0
            rows[key] = score,source
    newborn = rows["-",((2,1,0),)][0]

    @cache
    def at(shape, physical, corrupt):
        n = len(physical)
        if any(c > cap for (_,c),cap in zip(physical,C[n-1])):
            return 0,0
        coords = []
        for (d,c),knot in zip(physical,D[n-1]):
            if d < knot:
                coords.append([((d,c,0),Q(1))])
            else:
                w = Q(knot,d)
                if corrupt == "flat":
                    w = Q(1)
                elif corrupt == "shift":
                    w = Q(knot-1,d)
                coords.append([((knot,c,1),w),((knot,c,2),1-w)])
        total = [Q(0),Q(0)]
        for terms in product(*coords):
            weight = Q(1)
            for _,w in terms: weight *= w
            if weight:
                vals = rows[shape,tuple(s for s,_ in terms)]
                for j in range(2): total[j] += weight*vals[j]
        return tuple(total)

    @cache
    def terms(state,corrupt=""):
        p,_ = state
        _,ds,cs = geometry(state)
        score = source = 0
        for n in range(1,5):
            for nodes,shape in subsets(p,n):
                v,f = at(shape,tuple((ds[x],cs[x]) for x in nodes),corrupt)
                score += v
                source += f
        return score,source

    counts = [0,0,0]
    caught = set()
    gain_identities = 0
    # Compare the actual one-site jump against the expression used in a
    # tail corner, without the other lower-bound gaps masking a bad weight.
    exercised = set()
    for (shape,labels), (near, _) in rows.items():
        tail = [i for i,member in enumerate(labels) if member[2]]
        if len(tail)!=1 or labels[tail[0]][2]!=1:
            continue
        pos = tail[0]
        kind = shape,pos
        if kind in exercised: continue
        knot,c,_ = labels[pos]
        infinity = list(labels)
        infinity[pos] = knot,c,2
        far = rows[shape,tuple(infinity)][0]
        if near==far: continue
        exercised.add(kind)
        for d in (knot,knot+1,2*knot,1201):
            rate = d-c+Q(1,3)
            reference = Q(knot,d)*rate/(d+1)*(far-near)
            def measured(mutation):
                members = [(di,ci) for di,ci,_ in labels]
                members[pos] = d,c
                base = at(shape,tuple(members),mutation)[0]
                members[pos] = d+1,c
                new = at(shape,tuple(members),mutation)[0]
                return rate*(new-base)
            assert measured("") == reference
            for mutation in ("flat","shift"):
                if measured(mutation) != reference: caught.add(mutation)
            gain_identities += 1
    def test(state, mutants):
        initial,source = terms(state)
        moves = list(transitions(state))
        birth = sum(rate for _,rate,p in moves if p)
        gt = sum(rate*(terms(nxt)[0]-initial) for nxt,rate,_ in moves)
        required = DISCOUNT*initial-gt+newborn*birth
        assert source >= required, (state,source,required)
        if source < required+birth*2**40: caught.add("newborn")
        counts[0] += 1
        counts[1] += len(moves)
        counts[2] += bool(source)
        for corruption in mutants:
            if corruption in caught: continue
            wrong,_ = terms(state,corruption)
            gwrong = sum(rate*(terms(nxt,corruption)[0]-wrong) for nxt,rate,_ in moves)
            if source < DISCOUNT*wrong-gwrong+newborn*birth:
                caught.add(corruption)

    layer = {((-1,),(2,))}
    for size in range(3,10):
        for state in layer: test(state,())
        if size < 9:
            layer = {nxt for state in layer for nxt,_,_ in transitions(state)}
    for p in [(-1,0,1,2,3,4),(-1,0,1,1,2,3),(-1,0,1,2,2,3),(-1,0,1,1,1,2),
              (-1,0,1,1,2,2,3,3),(-1,0,1,1,1,1,1,1)]:
        core = [p.count(i)+(i!=0) for i in range(len(p))]
        for tag in range(80):
            degrees = (2,3,4,5,6,7,8,9,11,15,16,17,35,36,37,77,1009,1000000003)
            leaves = tuple(max(c,2,degrees[(11*i + tag + tag//(1+i))%len(degrees)])-c
                           for i,c in enumerate(core))
            test((p,leaves),("flat","shift"))
    print("physical states",counts[0],"one-ring transitions",counts[1],"tail gain identities",gain_identities,
          "positive sources",counts[2],"caught",','.join(sorted(caught)))
    assert gain_identities>=4 and {"flat","shift","newborn"} <= caught


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver",type=Path)
    parser.add_argument("snapshot",type=Path)
    parser.add_argument("side",type=int,choices=(1,2))
    args = parser.parse_args()
    with TemporaryDirectory() as work:
        rows = Path(work)/"rows"
        subprocess.run([str(args.solver.resolve()),str(args.snapshot),str(args.side),"63","4",
                        argument(D),argument(C)],env={**os.environ,"FT2_TAIL_DUMP":str(rows)},check=True)
        audit(rows)


if __name__ == "__main__":
    main()
