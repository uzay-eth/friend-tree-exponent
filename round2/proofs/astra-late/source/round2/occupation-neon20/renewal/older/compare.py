"""Matched inner extension and frozen D18 pointwise diagnostics; read-only imports."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec
from functools import cache
import sys

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=spec_from_file_location(name,path);m=module_from_spec(spec);spec.loader.exec_module(m);return m


def diagnostic(side):
    A=load('d20renew',HERE.parent/'replay.py');E=A.Rows(A.Cert(Path('b3f20_'+('lo' if side=='lower' else 'hi')+'.bin'),side,verify=False))
    F=A.Fringe(E,18,4)
    @cache
    def v(t):
        if len(A.shape(t)[0])>4:return Q(0)
        n,d,rates=F.rates(t)
        if n>18:return Q(0)
        a,e=(sum(z[i] for z in F.charge(t).values()) if n>=17 else Q(0) for i in range(2))
        b=a+sum(r*v(u) for u,r in rates);u=v(A.P.plus(t)) if n<18 else 0
        val=A.down(min(b/d,(b+(u+e)/2)/(d+Q(1,2))))
        assert val>=0 and d*val<=b and (d+Q(1,2))*val<=b+(u+e)/2
        return val
    print(side,'complete size17..18 <=4 cores lower',float(A.LO*v(A.P.newborn(2))),flush=True)
    # Resolve D18's extract as its own module: never use its price in D20.
    sys.modules.pop('extract',None)
    B=load('d18rows',HERE.parents[2]/'occupation-neon'/'renewal'/'old-flow'/'rows.py')
    T=B.Rows(B.Cert(Path('b3f18_'+('lo' if side=='lower' else 'hi')+'.bin')))
    for x in (((8,3,1),(0,1,1),(3,0,1),(3,0,1)),((9,2,1),(0,1,1),(3,0,1),(1,0,1)),
              ((11,3,1),(0,1,1),(5,0,1),(5,0,1)),((11,3,1),(0,1,1),(5,0,1),(1,0,1))):
        s=x[:-1];r=Q(1,sum(x[-1]));a,ports=E.affine('hV',s)
        low20=a+sum(n*min(0,c/2) for _,_,n,c in ports)
        c=next(c for i,h,n,c in ports if i==0 and not h)
        low18=T.floor('hV',s)
        excess=c*r-min(0,c/2)
        assert min(low18,low20,excess)>=0
        print(side,'fixture',x,'D18/D20 floor and D20 edge excess',*[float(z) for z in (low18,low20,excess)],flush=True)


if __name__=='__main__':
    p=ArgumentParser();p.add_argument('side',choices=('lower','upper'));diagnostic(p.parse_args().side)
