"""Exact finite checks of the mean-share port algebra and terminal bounds."""

from functools import cache
from fractions import Fraction as Q
from itertools import product


@cache
def graph(par):
    a=[[] for _ in range(len(par)+1)]
    for i,v in enumerate(par,1): a[v].append(i); a[i].append(v)
    d=[len(x) for x in a]
    r=[sum(Q(1,d[j]) for j in x) for x in a]
    return a,d,r


def terminal(par,v):
    assert len(par)>=2
    a,d,_=graph(par)
    l=sum(d[w]==1 for w in a[v])
    if d[v]==1:
        p=d[a[v][0]]
        return Q(2,15*p), Q(1,p)-Q(1,(p+1)*(2*p+1))
    rho=sum(Q(1,d[w]) for w in a[v] if d[w]>1)
    lo=max(Q(0),Q(l-1),Q(2,5) if l else Q(0))
    hi=min(Q(d[v]),Q(l)+rho)
    return lo,hi


@cache
def endpoint(par,u,v,h):
    a,d,r=graph(par); n=len(a)
    if sorted(v if w==u else u if w==v else w for w in a[u])==sorted(a[v]): return Q(0),Q(0)
    x,y=terminal(par,u); s,t=terminal(par,v)
    lo,hi=Q(x-t,n),Q(y-s,n)
    if h:
        rows=[endpoint(par+(w,),u,v,h-1) for w in range(n)]
        lo=max(lo,sum(q*l for q,(l,_) in zip(r,rows))/n)
        hi=min(hi,sum(q*b for q,(_,b) in zip(r,rows))/n)
    assert lo<=hi
    return lo,hi


def apply(interval,w):
    lo,hi=interval
    return (lo*w,hi*w) if w>=0 else (hi*w,lo*w)


def port(par,h):
    a,d,r=graph(par); n=len(a)
    weight=[r[u]/(d[u]*(d[u]+1)) for u in range(n)]
    terms=[]
    for u in range(n):
        terms.append(apply(endpoint(par+(u,n),n,n+1,h),r[u]/(d[u]+1)))
        for v in a[u]:
            if u<v: terms.append(apply(endpoint(par+(u,v),n,n+1,h),weight[v]-weight[u]))
    return tuple(sum(x)/ (n*(n+1)) for x in zip(*terms))


def effect(par,u,v,h):
    """Dirichlet mass for a finite degree-share boundary, exactly."""
    d=graph(par)[1]
    if not h: return Q(d[u]-d[v],len(d))
    return sum(q*effect(par+(w,),u,v,h-1) for w,q in enumerate(graph(par)[2]))/len(d)


def check_ports():
    checked=0
    for n in range(2,6):
        for par in product(*(range(i) for i in range(1,n))):
            a,d,r=graph(par)
            for h in (0,1,2):
                def z(u,v): return effect(par+(u,v),n,n+1,h)
                direct=sum(r[u]*rr*z(u,v) for u in range(n) for v,rr in enumerate(graph(par+(u,))[2]))
                parts=sum(r[u]/(d[u]+1)*z(u,n) for u in range(n))
                for u,ns in enumerate(a):
                    for v in ns:
                        if u<v: parts+=(r[v]/(d[v]*(d[v]+1))-r[u]/(d[u]*(d[u]+1)))*z(u,v)
                assert direct==parts
                checked+=1
    return checked


def unrooted(par):
    a,_,_=graph(par)

    def word(v,p):
        return '('+''.join(sorted(word(w,v) for w in a[v] if w!=p))+')'

    return min(word(v,-1) for v in range(len(a)))


def advance(law):
    merged={}
    for par,mass in law:
        for w,r in enumerate(graph(par)[2]):
            child=par+(w,)
            key=unrooted(child)
            rep,old=merged.get(key,(child,Q(0)))
            merged[key]=(rep,old+mass*r/(len(par)+1))
    ans=list(merged.values())
    assert sum(m for _,m in ans)==1
    return ans


if __name__=='__main__':
    print('exact finite-share port rows:',check_ports())
    h=3
    l,u=endpoint((0,0),0,1,h)
    assert 0<l<=u
    print('mean gap 2-3:',l/2,u/2)
    law=[((0,),Q(1))]
    for n,shapes in enumerate((1,1,2,3,6),2):
        assert len(law)==shapes
        lo,hi=Q(0),Q(0)
        for par,m in law:
            l,u=port(par,h)
            lo+=m*l
            hi+=m*u
        print(f'mean gap {n+1}-{n+2}:',lo,hi)
        assert lo<=hi
        assert lo>0
        law=advance(law)
