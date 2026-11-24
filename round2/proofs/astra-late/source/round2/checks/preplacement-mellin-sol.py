"""Exact rational enclosures of the serial-middle discounted stock test."""

from fractions import Fraction as Q
from math import factorial


def box(left, right, cutoff, stop, scale=10**15):
    """Lower and upper integer ticks at every state for p in [left,right]."""
    def outward(lo, hi):
        return lo.numerator//lo.denominator, -((-hi.numerator)//hi.denominator)

    V = {}
    for j in range(cutoff):
        if j < 2:
            V[cutoff,j]=(0,0)
            continue
        # A positive polynomial in p even when p is noninteger.
        lo=hi=Q(scale,factorial(j-2))
        for k in range(1,j-1):
            lo*=left+k
            hi*=right+k
        V[cutoff,j]=outward(lo,hi)
    min_return=None
    for d in range(cutoff-1,1,-1):
        R=(0,0)
        for k in range(stop-1,1,-1):
            send=Q(1,k)
            gain=k-1+Q(1,d)
            prev=R
            R=outward((send*scale+gain*prev[0])/(right+send+gain),
                      (send*scale+gain*prev[1])/(left+send+gain))
            # Prove the needed sign without differencing two close bounds:
            # R(k)-R(k+1) has sign send-(p+send)R(k+1).
            slack=send*scale-(right+send)*prev[1]
            min_return=slack if min_return is None else min(min_return,slack)
        hi,lo=V[d+1,1][1],V[d+1,1][0]
        V[d,0]=outward(Q(lo*R[0],scale),Q(hi*R[1],scale))
        for j in range(1,d):
            g=V[d+1,j+1]; loss=V[d,j-1]
            V[d,j]=outward((d*g[0]+loss[0])/(d+1+right*d/j),
                           (d*g[1]+loss[1])/(d+1+left*d/j))
    seams=[V[d+1,j+1][0]-V[d,j][1] for d in range(2,cutoff) for j in range(1,d)]
    return min(seams), min_return, V[2,1], V[3,2]


def drift():
    rows=0
    for p in (Q(1),Q(5,4),Q(2),Q(19,7),Q(7)):
        for j in range(2,24):
            gain=1+p/(j-1)
            loss=1-p/(j+p-2)
            assert j*(gain-1)+(loss-1)-p==p*(1/Q(j-1)-1/(j+p-2))>=0
            rows+=1
    print(f'PASS {rows} exact terminal drift rows, including equality at p=1')


def large_power():
    def rising(a,p):
        out=Q(1)
        for i in range(p): out*=a+i
        return out/factorial(p)
    lower=sum(Q((-1)**k,factorial(k)) for k in range(8))
    debit=(4*rising(Q(1,2),7)+2*rising(Q(1,3),7))/3
    assert lower==Q(103,280)>debit
    print('PASS p>=7 unit-loss Taylor floor strictly beats the descending Gamma debit at p=7')


if __name__=='__main__':
    drift()
    large_power()
    for p in (1,2,3):
        seam,ret,mid,_=box(Q(p),Q(p),30,40)
        assert seam>=0 and ret>=0
        print(f'p={p}: baseline lower tick {mid[0]}/10^15 (D0=30,K=40)')
    for D,mesh in ((8,32),(12,32),(20,64)):
        minimum=10**15
        for i in range(mesh,7*mesh):
            seam,ret,mid,_=box(Q(i,mesh),Q(i+1,mesh),D,40)
            assert seam>0 and ret>0, (D,i,seam,ret)
            minimum=min(minimum,mid[0])
        print(f'PASS {6*mesh} rational boxes on [1,7] at D0={D}, K=40; uniform lower tick {minimum}/10^15',flush=True)
