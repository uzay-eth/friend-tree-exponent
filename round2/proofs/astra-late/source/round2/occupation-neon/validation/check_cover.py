"""Independent sparse operand check and fixed-size residual forecast check."""
from fractions import Fraction as Q
from functools import cache
from pathlib import Path
import sys

from forecast import LO,HI,A,B,C,tail

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'higher/older'))
from replay import Cert, Local


def check_row():
    rows=[Local(Cert(Path(f'b3f18_{side}.bin'))) for side in ('lo','hi')]
    x=((0,1,1),)*4
    for p,d in ((Q(0),2),(Q(1,2),8),(Q(1,7),123),(Q(1,77),2**21)):
        ext=((p,0),(0,0),(0,0),(0,Q(1,d)))
        a,b=[row.charge('P4',x,ext) for row in rows]
        assert a-Q('1.2415')*b==A+B*p+C/d
    print('four paired sparse full-tail coefficient checks passed')


def fixed_size():
    # Freeze one admissible policy only as a comparison baseline, never as
    # a model for external inputs. Follow it until d>K or face death.
    K=12; s=(LO+HI)/2; p=Q(1,10)
    def operands(d,k):
        l=d-k-1
        g=l+Q(k+2,4); b=Q(l,d); q=A+B*p+C/d
        return g,b,q,s+3+Q(1,d)+p+g+b
    u,v=tail(K,True,True)
    @cache
    def seed(d,k):
        if d>K:return u-v/d
        g,b,q,z=operands(d,k)
        return (q+g*seed(d+1,k)+(b*seed(d,k+1) if b else 0))/z
    # Deliberately approximate with an independent fixed-size cut at d=8.
    @cache
    def approximate(d,k):
        if d>8:return Q()
        g,b,q,z=operands(d,k)
        return (q+g*approximate(d+1,k)+(b*approximate(d,k+1) if b else 0))/z
    # Chronological defect transport independently forward on (d,k).
    flow={(2,0):Q(1)}; add=Q()
    for d in range(2,K+1):
        for k in range(d):
            w=flow.pop((d,k),Q());g,b,q,z=operands(d,k)
            defect=(q+g*approximate(d+1,k)+(b*approximate(d,k+1) if b else 0))/z-approximate(d,k)
            add+=w*defect
            for y,rate in (((d+1,k),g),((d,k+1),b)):
                if rate:flow[y]=flow.get(y,Q())+w*rate/z
    seam=sum(w*(u-v/d) for (d,k),w in flow.items())
    assert approximate(2,0)+add+seam==seed(2,0)
    print('fixed-policy exact forward/backward telescoping:',float(approximate(2,0)),float(add),float(seam),float(seed(2,0)))


if __name__=='__main__':
    check_row()
    fixed_size()
