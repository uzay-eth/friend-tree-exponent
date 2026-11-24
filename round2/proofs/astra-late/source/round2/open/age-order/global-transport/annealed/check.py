"""Exact checks on actual successive jumps; no enumeration of old histories."""
from fractions import Fraction as Q


def tree(n, edges):
    a=[[] for _ in range(n)]
    for u,v in edges:
        a[u].append(v)
        a[v].append(u)
    return tuple(tuple(row) for row in a)


def rates(a):
    r=[Q(0) for _ in a]
    for adj in a:
        for v in adj:
            r[v]+=Q(1,len(adj))
    assert sum(r)==len(a)
    return r


def jump(a,u):
    n=len(a)
    return a[:u]+(a[u]+(n,),)+a[u+1:]+((u,),)


def pending(a,j):
    return Q(1,len(a[a[j][0]])) if len(a[j])==1 else Q(0)


def inspect(a):
    n=len(a); r=rates(a)
    d=list(map(len,a)); qs=[Q(1,k+1) for k in d]
    w=[r[i]*qs[i] for i in range(n)]; s=[w[i]/d[i] for i in range(n)]
    total=sum(w)
    E=sum((s[v]-s[u])*(qs[u]-qs[v]) for u in range(n) for v in a[u] if u<v)
    first=shared=old_step=younger_hazard=Q(0)
    for u,R in enumerate(r):
        b=jump(a,u); r2=rates(b)
        total2=sum(R2/(len(adj)+1) for R2,adj in zip(r2,b))
        delta=(3*d[u]+4-2*R)/(2*(d[u]+1)*(d[u]+2))-sum(qs[v] for v in a[u])/ (d[u]*(d[u]+1))
        assert total2-total==delta
        old_step+=R/n*total2
        for v,R2 in enumerate(r2):
            c=jump(b,v); p=R/n*R2/(n+1)
            first+=p*(v==n)
            shared+=p*(pending(c,n)-pending(c,n+1))/(n+2)
            younger_hazard+=p*pending(c,n+1)
    assert first==total/(n*(n+1))
    assert shared==(E-total/2)/(n*(n+1)*(n+2))
    assert shared<=-sum(w[i]*qs[i] for i in range(n))/(n*(n+1)*(n+2))<0
    assert first+shared>0
    assert younger_hazard==old_step/(n+1)
    assert old_step<(1+Q(3,2*n))*total
    # Geometric age weights. The infinite upper charge is bounded by its
    # geometric sum, without any approximation of the actual pending law.
    for ratio,floor in ((Q(1,2),Q(1)),(Q(1,3),Q(1,2)),(Q(3,7),Q(3,4))):
        assert younger_hazard*ratio/(1-ratio)/(n+2)<first*floor
    return shared


if __name__=='__main__':
    fixtures=[tree(n,edges) for n,edges in (
        (2,[(0,1)]), (3,[(0,1),(1,2)]),
        (6,[(j,j+1) for j in range(5)]),
        (12,[(0,j) for j in range(1,12)]),
        (13,[(j,j+1) for j in range(6)]+[(j,j+7) for j in range(6)]),
        (15,[(0,j) for j in (1,5,10)]+[(j,j+1) for j in (1,2,3,5,6,7,8,10,11,12,13)]),
    )]
    one_jump_count=0
    for a in fixtures:
        assert inspect(a)<0
        mean=sum(p/len(a)*inspect(jump(a,u)) for u,p in enumerate(rates(a)))
        assert mean<0
        one_jump_count+=len(a)
    print(f'PASS {len(fixtures)} physical seeds, {one_jump_count} true old-jump successors, both birth orders and discounted budgets')
