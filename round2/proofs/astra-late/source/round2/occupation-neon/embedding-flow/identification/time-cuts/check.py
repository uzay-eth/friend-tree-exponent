"""Birth-conditioned edge holding cuts; use the existing stationary LP."""
import argparse
from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path[:0]=[str(HERE.parent),str(HERE.parents[1])]
from solve import IncidenceModel, Pattern, optimize
from exact import install_cuts

GRID=10**15


def stock(u):
    return u[0]-u[1]-1


def rates(p,x,mother=Q(1,2)):
    """Exact included-edge rates, adaptive nonnegative outside input intervals."""
    d,k=p
    e,j=x
    a,b=Q(stock(p))+Q(1,e),Q(stock(x))+Q(1,d)
    return (a,a+Q(k-1,2)+mother),(b,b+Q(j,2)),Q(stock(p),d),Q(stock(x),e)


def price(model, reward, maximize, mother=Q(1,2)):
    """Discounted reward from an actual promotion of the distinguished child."""
    D=max(p[0] for p,x in reward)
    E=max(x[0] for p,x in reward)
    ext=max if maximize else min
    s=model.high if not maximize else model.low

    @cache
    def V(p,x):
        d,k=p
        e,j=x
        if d>D or e>E:
            return Q(0)
        a,b,bp,bx=rates(p,x,mother)
        const=reward.get((p,x),0)
        if bp:
            const+=bp*V((d,k+1),x)
        if bx:
            const+=bx*V(p,(e,j+1))
        v,w=V((d+1,k),x),V(p,(e+1,j))
        vals=[(const+g*v+h*w)/(s+bp+bx+g+h) for g in a for h in b]
        bound=ext(vals)
        out=Q(-(-bound*GRID//1) if maximize else bound*GRID//1,GRID)
        assert all((out-z)*(1 if maximize else -1)>=0 for z in vals)
        return out

    return {t:Q(stock(t),t.degree)*V((t.degree,t.arity+1),(2,0))
            for t in model.levels[1] if t.degree<=D and stock(t)}


def objective(model,reward):
    row=defaultdict(Q)
    for ((d,k),(e,j)), w in reward.items():
        model.mass(row,Pattern(d,k,(Pattern(e,j),)),w)
    return dict(row)


def install(model,rewards, stratify=False):
    for reward in rewards:
        assert all(v>=0 for v in reward.values())
        cost=objective(model,reward)
        for maximize in (False,True):
            sg=1 if maximize else -1
            row=defaultdict(Q,{i:sg*v for i,v in cost.items()})
            p0=price(model,reward,maximize)
            for t,v in p0.items():
                model.mass(row,t,-sg*v)
            if stratify:
                for d in range(3,model.cap+1):
                    pd=price(model,reward,maximize,Q(1,d))
                    for t,v in pd.items():
                        for k in range(1,d):
                            model.mass(row,Pattern(d,k,(t,)),-sg*(v-p0[t]))
                # Split the unobserved birth-grandmother tail too. Use a
                # coarsened upper degree-independent reciprocal on the tail.
                pd=price(model,reward,maximize,Q(1,model.cap+1))
                for t,v in pd.items():
                    model.mass(row,t,-sg*(v-p0[t]))
                    for d in range(2,model.cap+1):
                        for k in range(1,d):
                            model.mass(row,Pattern(d,k,(t,)),sg*(v-p0[t]))
            model.record(model.upper_rows,row,0)


def compare(matched, graded, test):
    out=[]
    for m in (matched,graded):
        c=objective(m,test)
        out.append(tuple(optimize(m,c,up)[0] for up in (False,True)))
    widths=[b-a for a,b in out]
    return out, widths


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--cap',type=int,default=4)
    p.add_argument('--order',type=int,default=2)
    p.add_argument('--share-through',type=int,default=2)
    p.add_argument('--held-out',action='store_true')
    a=p.parse_args()
    m,n=(IncidenceModel(a.cap,a.order,a.share_through) for _ in range(2))
    for v in (m,n):
        install_cuts(v)
    assert a.cap>=4 and a.order>=2
    # All single edge-census payoffs, optionally whole virgin-child shell sums.
    unit=[{((t.degree,t.arity),(t.children[0].degree,t.children[0].arity)):Q(1)}
          for t in n.levels[2]]
    shell=[{key:Q(1) for r in unit for key in r
            if key[0][0]==d and key[1]==(2,0)} for d in range(2,a.cap+1)]
    install(n,unit if a.held_out else unit+shell,True)
    print('columns/base upper/added',len(m.variables),len(m.upper_rows),len(n.upper_rows)-len(m.upper_rows),flush=True)
    # Closed-form child/parent residence from degree-two promotion into
    # the edge. At (2,1)->(2,0), fixed kill rate is 5/2, outside adds [0,1/2].
    anchor={((2,1),(2,0)):Q(1)}
    t=Pattern(2,0)
    assert abs(price(n,anchor,False)[t]-Q(1,2)/(n.high+3))<=Q(1,GRID)
    assert abs(price(n,anchor,True)[t]-Q(1,2)/(n.low+Q(5,2)))<=Q(1,GRID)
    for name,target in [('virgin pair',anchor),('D3 to virgin',shell[1]),('D4 to virgin',shell[2]),
                        ('D4,k2 to virgin',{((4,2),(2,0)):Q(1)}),
                        ('D4,k2 to D3,k1',{((4,2),(3,1)):Q(1)})]:
        intervals,widths=compare(m,n,target)
        print(name,[[float(x) for x in v] for v in intervals],
              'width reduction',float(1-widths[1]/widths[0]) if widths[0]>0 else None,flush=True)
