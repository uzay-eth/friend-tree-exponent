"""Bare-degree-two suppression with actual low endpoint exposures."""
from itertools import combinations_with_replacement, product
from pathlib import Path
from runpy import run_path

E=run_path(str(Path(__file__).resolve().parents[1]/'edge-dual/check.py'))
Q,G,U,L,P,B,F,D2,H=(E[n] for n in ('Q','G','U','L','P','B','F','D2','H'))

u=U(3,0)
f=max(0,u-L(3))
g=G(3)


def low(d,ds,r):
    s,C,v=P(d,tuple(x or r for x in ds))
    return Q(1,d),L(d+1),s,v,(d,ds,r)


def port(x,t):
    _, l,s,v,_=t
    return max(0,(x-s)*max(0,u-l),(s-x)*max(0,v-L(3)))


def highport(x,t):
    return max(H,x*max(0,u-t[1]))


def debit(A,B,k,pa=port,pb=port):
    x,y=A[0],B[0]
    ss=((x+y)/6,) if k==1 else ((Q(1,2)+x)/6,)+(Q(1,6),)*(k-2)+((Q(1,2)+y)/6,)
    return pa(ss[0],A)+pb(ss[-1],B)+f*sum(abs(x-y) for x,y in zip(ss,ss[1:]))-2*g*sum(ss)


def verify():
    levels=[]
    for d in range(2,10):
        rows=[]
        for c in range(1,d+1):
            if d==c==2: continue
            for r in (0,2): rows.append(low(d,(2,)+(r,)*(c-1),0))
        # For fixed degree and stock, only coefficient extrema matter.
        group={}
        for p in rows: group.setdefault(len(p[-1][1]),[]).append(p)
        levels.append([p for rows in group.values() for p in (min(rows,key=lambda p:p[2]),max(rows,key=lambda p:p[2]))])
    highs=[(Q(1,d),L(d+1)) for d in range(10,31)]+[(Q(0),Q(0)),(Q(1,31),Q(0))]
    def bill(t):
        d,ds,_=t[-1]
        return Q('.021') if d==2 else Q('.0037')
    lows=sum(levels,[])
    for k in (1,2,3):
        p=max(debit(A,Z,k)-bill(A)-bill(Z) for A,Z in product(lows,repeat=2))
        q=max(debit(A,Z,k,pb=highport)-bill(A)-D2 for A,Z in product(lows,highs))
        print(f'PASS chain length {k}: residual low-low {float(p):.12f}, low-high {float(q):.12f}')
        assert max(p,q)<B
    for d in range(2,10):
        worst=Q(1)
        for c in range(1,d+1):
            if d==c==2:continue
            for ds in combinations_with_replacement((0,*range(2,10)),c):
                # Relax several actual degree-two slots as chain ports. A direct
                # degree-two terminal is kept separately by maximizing the bill.
                for r in (0,10):
                    s,C,v=P(d,tuple(e or r for e in ds))
                    costs=[max(0,E['upper_out'](d,e,s,v)-B) for e in ds]
                    if c==1:
                        val=C-E['A']-max(B+costs[0],B+bill(low(d,ds,r)) if ds[0]==2 else 0)
                    else:
                        costs=[max(cost,bill(low(d,ds,r))) if e==2 else cost for cost,e in zip(costs,ds)]
                        val=C+E['A']*(c-2)-B-sum(costs)+(F*ds.count(0) if (d,c) in ((3,2),(3,3),(4,4)) else 0)
                    worst=min(worst,val)
        print(f'PASS endpoint degree {d}: minimum residual {float(worst):.12f}',flush=True)
        assert worst>0


if __name__=='__main__': verify()
