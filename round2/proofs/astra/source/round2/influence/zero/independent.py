#!/usr/bin/env python3
"""Cold finite replay: vertex-ID arrows and full signed generator, no Green solve."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from pathlib import Path
import sys

from solver import occupations, enclosure, upper, lower

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'current-mark'/'eval'))
from check import directed


@cache
def n(S): return 1+S[0]+sum(n(Y) for Y in S[1])


def replay(N,s,plant=False):
    todo=[(1,())]; rows={}
    while todo:
        S=todo.pop()
        if S in rows: continue
        A,_=directed(S)
        d=1+S[0]+len(S[1]); p=Q(1,d)
        es=defaultdict(Q,{Z:r for Z,r in A.items() if n(Z)<=N})
        for Y in S[1]:
            es[Y[0]+1,Y[1]]+=p
            if not plant: es[Y]-=p
        todo.extend(es)
        rows[S]=(s+n(S)-p,S[0]*p,dict(es))
    domain,V,T,f=occupations(N,s)
    x,err=enclosure(T,f,[Q(1)]*len(T))
    h=dict(zip(domain,x))
    candidate={S:b+sum(t*h[Y] for Y,t in a.items()) for S,(b,a) in V.items()}
    if rows.keys()!=candidate.keys(): raise ValueError('domain mismatch')
    err=max(abs(d*candidate[S]-b-sum(t*candidate[Z] for Z,t in a.items())) for S,(d,b,a) in rows.items())
    if err>Q(1,10**12): raise ValueError(f'full signed residual {upper(err)}')
    # At N=3, J gain or promotion are the only successors. The bare pair
    # pays zero direct reward and p[U(plus(J))-U(J)] at p=1/2.
    if N==3:
        x=Q(2,3)/(s+Q(8,3)); z=(Q(1,2)+x+Q(1,4)*x/(s+Q(5,2)))/(s+Q(3,2)+Q(1,4)/(s+Q(5,2)))
        if abs(z-candidate[1,()])>Q(1,10**15): raise ValueError('hand N=3 differs')
    print(f'PASS N={N}: {len(rows)} explicit signed rows, full residual <= {upper(err)}')


if __name__=='__main__':
    try: replay(int(sys.argv[1]),Q('0.565728'),'--plant' in sys.argv)
    except ValueError as exc: sys.exit(f'REJECT: {exc}')
