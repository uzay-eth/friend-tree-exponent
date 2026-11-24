"""Source-normalized bounds for old-copy exit coefficients."""
from fractions import Fraction as Q
from pathlib import Path
import re

LO,HI=Q('.565676636772'),Q('.565771871520')
A=Q('.56564');R=Q('1.006');THETA=Q('.976433602362')


def majorant(s):
    # (C+k e_2^T)g <= R*g, Cg<=theta*g on the same Perron vector.
    # Multiply by the positive untriggered inverse and telescope.
    gamma=Q('.565775')  # a rational common discount above both trials
    # The independent direct trigger is 1/[2(s+2)] at the source.
    # Obtain an extra discount gain on the whole source by factoring its
    # first virgin holding time, which has intrinsic + death rate 2.
    assert A<=s<=gamma
    return (1+(R-1)/(1-THETA))*(A+2)/(s+2)


def bounds(s):
    # Outward trial bounds from the reviewed full/nonvirgin pilot and exact
    # finite 2-copy propagation (D<=5, k<=3, 96 uniformized steps).
    data={LO:('0.3388205426','0.4620492561','0.0008341635188'),
          HI:('0.3387544144','0.4618377419','0.0008339351325')}
    l01,u01,l2=map(Q,data[s]);u=majorant(s)
    # Same full-domain source; disjoint chronological classes.
    u2=u-l01
    assert 0<l2<u2<1
    return l01,u01,l2,u2,u


def fixture():
    # Exact transform: U=floor + P U in a stopped marked first-passage
    # fixture.  The retained harmonic mass is unchanged across an exit.
    # State 0 goes to 1 or exits directly; state 1 exits or triggers.
    p01,p03,p12,h=Q(1,4),Q(1,5),Q(1,2),Q(1,7)
    terminal=(Q(3,8),Q(6,7))
    u1=h+p12*terminal[0];u0=p01*u1+p03*terminal[1]
    mass=(p01*p12*terminal[0]+p03*terminal[1])/u0
    assert mass+p01*h/u0==1
    # Copy-weighted Cauchy, without moving the boundary.
    m=(p01*p12,p03)
    r=tuple(v/u0 for v in terminal)
    assert (sum(x*y for x,y in zip(m,r)))**2<sum(m)*sum(x*y*y for x,y in zip(m,r))
    # Difference martingale in source-normalized units: coefficient 2
    # can still get almost all the slack after knowing f0+f1 and F.
    for s in (LO,HI):
        l,u,l2,u2,F=bounds(s)
        eps=Q(1,10**8)
        x=(l,l2,F-l-l2)
        y=(l,u2-eps,eps)
        assert sum(x)==sum(y)==F and min(*x,*y)>0 and l<=u
    print('exact harmonic mass, Cauchy, and scalar-complement fixtures passed')


def main():
    root=Path(__file__).resolve().parents[7]
    log=(root/'outputs'/'finite-physical-trigger-nineteen-fine.log').read_text()
    theta=Q(int(re.search(r'contraction <= (\d+)',log)[1]),10**12)
    assert theta==THETA and 'upper PF multiplier 503/500' in log
    g2=Q(int(re.search(r'blocks 2=(\d+)',log)[1]),10**12)
    assert R-g2/(2*(A+Q(1,20)+Q(3,2)))<=theta
    p=(Q('.5658')-A)/(Q('.5658')-Q('.51'))
    assert Q(50,7)**p.numerator<R**p.denominator
    fixture()
    for s in (LO,HI):
        l,u,l2,u2,F=bounds(s)
        print('trial',s,'F upper',float(F),'0+1 interval',float(l),float(u),
              '2 interval',float(l2),float(u2),'0+1+2 lower',float(l+l2),flush=True)
        print('unresolved two-copy budget',float(u2-l2),
              'max normalized unresolved budget',float(1-(l+l2)/F),flush=True)


if __name__=='__main__':main()
