#!/usr/bin/env python3
"""Cold exact Neon D18 ledger and high-parent Bellman/physical fixtures."""
from fractions import Fraction as Q
from hashlib import file_digest
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import comb
from pathlib import Path
import argparse
import sys

ROOT = Path(__file__).resolve().parents[2]
BASE = ROOT/'round2/occupation-neon'
S = 10**24
LO, HI = Q('.565673570322'), Q('.565775')
PINS = {
    'lower': '9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
    'upper': '1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41',
}
FIRST = (Q('2.512846393e-6'), Q('2.025752646e-6'))
HIGHER = (Q('1.276902849680e-6'), Q('1.501029458674e-6'))
SUFFIX = (Q(58961819710386489,10**23), Q(46635843154908967,10**23))
DEN = (Q('1.356833021311'), Q('1.367169995210'))
TARGET = (565676634491,565771874464)
TRIAL = (Q('.565673'), HI)


def load(name,path):
    spec=spec_from_file_location(name,path)
    mod=module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def physical_census():
    # These are intensities on a concrete tree, independent of the table
    # readers. Recompute the projected state after each labelled arrow.
    parents=(-1,0,1,1,1,1,4)
    leaves=(6,13,1,1,0,2,1)
    old=[[] for _ in parents]
    def link(a,i,g): g[a].append(i);g[i].append(a)
    for i,a in enumerate(parents[1:],1): link(a,i,old)
    for i,n in enumerate(leaves):
        for _ in range(n):
            v=len(old);old.append([]);link(i,v,old)
    def state(g):
        dd=list(map(len,g)); hub=1
        kids=[w for w in g[hub] if w!=0 and dd[w]>=2]
        return dd[hub],len(kids),sum(dd[w]==2 and any(dd[x]==1 for x in g[w]) for w in kids)
    d,k,j=state(old)
    assert (d,k,j)==(18,4,2)
    rate={}
    for a,neighbors in enumerate(old):
        for w in neighbors:
            new=[u[:] for u in old]+[[]];link(w,len(old),new)
            q=state(new)
            rate[q]=rate.get(q,Q(0))+Q(1,len(neighbors))
    assert rate[d+1,k,j]==13+Q(2,2)+Q(1,2)+Q(1,3)+Q(1,7)
    assert rate[d,k+1,j+1]==Q(13,18)
    assert rate[d,k,j-1]==2*(Q(3,2)+Q(1,18))
    assert set(rate)=={(d,k,j),(d+1,k,j),(d,k+1,j+1),(d,k,j-1)}
    print('PASS: physical asymmetric D=18 census: two virgins, bare child, degree-3 young child, degree-7 parent')


def overlap_controls():
    # Predicates are on a charged row's actual owner, not on the enclosing
    # fringe in which it might also appear as a descendant.
    def labels(f,n,c,k,d,w):
        return {
            'f': f=='f',
            'size12': f in ('hV','hC','P4','Yt','Yb','S3') and n<=12 and d<=16,
            'suffix-H': f in ('hV','hC','P4','Yt','Yb','S3') and 13<=n<=20 and 3<=c<=4 and d<=16,
            'suffix-g': f=='g' and n<=20 and c<=4 and d<18 and w<18,
            'HP-g': f=='g' and 18<=d<=191 and 1<=k<=10 and w==2,
            'HP-pair': f=='hC' and 18<=d<=191 and 4<=k<=10 and w==2,
        }
    # For overlap only this is a superset: omit some arm and support tests.
    for f,n,c,k,d,w in product(('f','g','hC','P4'),(12,13,18,20,21),(3,4,5),(3,4,10,11),(16,17,18,191,192),(2,3,18)):
        assert sum(labels(f,n,c,k,d,w).values())<=1
    assert labels('g',20,3,2,18,2)['HP-g']  # overlaps observed-virgin young objective
    assert labels('hC',22,5,4,18,2)['HP-pair']
    assert not any(labels('g',22,5,4,18,3).values())  # nonvirgin mixed edge
    print('PASS: charged-row seam predicates; overlapping young/mixed strata distinguished')


def mixed_floors(C,d,k):
    # Separately maintained general physical all-eight charge extractor.
    # Independently optimize over a per-donor [0,1/2] box, no imported
    # high-parent formulas, star checker, or Bellman solver.
    host=(d-k-1,k,1); virgin=(1,0,1)
    def unit(f,count):
        a,terms=C.terms(f,(host,)+(virgin,)*count)
        assert all(i==0 for i,_,_,_ in terms)
        val=a+sum(n*min(Q(0),dq/2) for _,_,n,dq in terms)
        assert val>=0
        return val
    return unit('g',1), unit('hC',2) if k>=4 else Q(0)


def backwards(rows,mode):
    # Store one complete degree layer. Unlike the output replay, explicitly
    # evaluate both ratios and test both cross-multiplied inequalities.
    next_degree={}
    for d in range(191,1,-1):
        level={}
        for k in reversed(range(min(d-1,10)+1)):
            l=d-1-k
            b=Q(l,d)
            e,p=rows.get((d,k),(0,0))
            for j in range(k+1):
                a=j*(Q(3,2)+Q(1,d))
                r=(e*j if mode!=2 else 0)+(p*comb(j,2) if mode!=1 else 0)
                reward=r*S+b*level.get((k+1,j+1),0)+a*level.get((k,j-1),0)
                old=next_degree.get((k,j),0)
                corners=[(reward+g*old,HI+g+a+b) for g in (l+Q(j,2),l+Q(k+1,2))]
                v=min(n/q for n,q in corners)//1
                assert all(n>=q*v for n,q in corners)
                level[k,j]=v
        next_degree=level
    return LO*Q(next_degree[0,0],S)


def support_overlap_witness(charges):
    # Degree-3 young phases overlap the nonvirgin mixed pilot. Check
    # representative instantiated rows, allowing a zero floor.
    host=(13,4,1); v=(1,0,1); young=(2,0,1)
    for f,arms in (('g',(v,)),('hC',(v,v)),('g',(young,))):
        val,free=charges.terms(f,(host,)+arms)
        assert val+sum(num*min(0,w/2) for _,_,num,w in free)>=0


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('lower',type=Path);p.add_argument('upper',type=Path)
    a=p.parse_args()
    physical_census();overlap_controls()
    sys.path.insert(0,str(BASE/'fringe'))
    general=load('cold_general_charges',BASE/'fringe/charges.py')
    published=load('comparison_high_parent',BASE/'high-parent/replay.py')
    for i,(side,path) in enumerate(zip(PINS,(a.lower,a.upper))):
        before=path.stat()
        assert before.st_size==3495091968
        with path.open('rb') as f: assert file_digest(f,'sha256').hexdigest()==PINS[side]
        cert=general.Cert(path)
        sign=(1,-1)[i]
        assert (cert.mode,cert.lam,cert.D,cert.K,cert.S)==((2,1)[i],TRIAL[i],18,17,10**11)
        cold=general.Charges(cert); comparison=published.Physical(cert)
        support_overlap_witness(cold)
        rows={}
        for d,k in product(range(18,192),range(1,11)):
            r=mixed_floors(cold,d,k)
            assert r==published.charges(comparison,d,k),(d,k,r)
            rows[d,k]=r
        joint,edge,pair=(backwards(rows,j) for j in (0,1,2))
        assert joint>edge+pair>0
        # Build the common numerator from its independent ledger terms
        # (suffix fraction comes from the cap-20 replay, not PAID).
        cumulative=FIRST[i]+HIGHER[i]+SUFFIX[i]
        total=cumulative+joint
        x=TRIAL[i]+sign*total/DEN[i]
        ticks=sign*((sign*x*10**12)//1)
        assert ticks==TARGET[i]
        print(side,'1740 exact g/hC comparisons; joint',joint,'edge',edge,'pair',pair,flush=True)
        print(side,'first+size12+size20',cumulative,'total',total,'exact endpoint',x,
              f'outward 0.{ticks:012d}',flush=True)
        after=path.stat()
        assert (before.st_ino,before.st_size,before.st_mtime_ns)==(after.st_ino,after.st_size,after.st_mtime_ns)


if __name__=='__main__': main()
