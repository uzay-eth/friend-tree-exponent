"""Bounded independent stochastic and mixed-ledger controls."""
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
from runpy import run_path


def arrows(t):
    return [sum(Q(1,len(t[v])) for v in row) for row in t]


def birth(t,u):
    return t[:u]+(t[u]+(len(t),),)+t[u+1:]+((u,),)


def cancellation():
    cases=0
    for n in range(2,7):
        for ps in product(*(range(i) for i in range(1,n))):
            t=((1,),(0,))
            for u in ps[1:]:
                t=birth(t,u)
            rates=arrows(t)
            for u in range(n):
                fresh=arrows(birth(t,u))
                assert rates[u]*fresh[n]==rates[u]/(len(t[u])+1)
                assert fresh[u]==rates[u]+1
                for v in range(u):
                    delta=lambda q: rates[q]/(len(t[q])*(len(t[q])+1))
                    assert rates[u]*fresh[v]-rates[v]*arrows(birth(t,v))[u] == (
                        delta(v)-delta(u) if v in t[u] else 0)
                    cases+=1
    print(f'PASS {cases} literal two-order physical coefficient comparisons')


def stochastic():
    # Independent backward tables, no producer imports.
    R={}
    for d in range(1,20):
        R[d,32]=Q(0)
        for k in range(31,1,-1):
            hop=k*(k-1+Q(1,d))
            R[d,k]=(1+hop*R[d,k+1])/(1+hop)
            assert R[d,k]>=R[d,k+1]
    W={(20,l):1-Q(1,20)**l for l in range(20)}
    for d in range(19,1,-1):
        W[d,0]=R[d,2]*W[d+1,1]
        for l in range(1,d):
            W[d,l]=(d*W[d+1,l+1]+W[d,l-1])/(d+1)
            assert W[d+1,l+1]>=W[d,l]>=W[d,l-1]
    P={(112,j):Q(1,j+1) for j in range(111)}
    for d in range(111,1,-1):
        for j in range(d-1):
            P[d,j]=(1+d*(j+1)*P[d+1,j+1]+j*P.get((d,j-1),0))/(1+d*(j+1)+j)
            assert P[d,j]>=P[d+1,j+1]
            assert P[d,j]<=Q(1,j+1)
            if j:
                assert P[d,j]<=P[d,j-1]
    assert W[2,1]-P[2,0]>Q(16694235,10**8)
    print('PASS 570 return stages, 171 rescue stages, 6105 upper stages and unrounded stock barrier')


def ledger():
    root=Path(__file__).resolve().parents[1]/'open/age-order/global-transport/hostile'
    ref=run_path(str(root/'edge-dual/check.py'))
    G,U,L=map(ref.get,('G','U','L'))
    a,b=Q('.00634'),Q('.017')
    # Core tree with all three chain types, direct high payments, high tail,
    # and a real direct degree-two terminal beside suppressed degree-two chains.
    vertices=[(0,1,2,3,4,5,6),(0,1,3,2,4,5,6),(0,1,1,2,3,5,7)]
    for lens in vertices:
        ds=[113,10,3,4,2,7,2,16]
        pairs=[(0,1),(0,2),(0,3),(1,4),(2,5),(3,6),(1,7)]
        adj=[[] for _ in ds]
        for (u,v),length in zip(pairs,lens):
            for _ in range(length):
                z=len(ds);ds.append(2);adj.append([])
                adj[u].append(z);adj[z].append(u);u=z
            adj[u].append(v);adj[v].append(u)
        s=[]; C=[]
        for d,ns in zip(ds,adj):
            l=d-len(ns); r=l+sum(Q(1,ds[z]) for z in ns); q=r/(d*(d+1))
            serial=d*q*G(d+1)+Q(l,2*d)*G(2)
            if l:
                sl=Q(1,2*d)
                serial-=l*max(0,(sl-q)*max(0,U(2,0)-L(d+1)),
                               (q-sl)*max(0,U(d+1,l-1)-L(2)))
            s.append(q); C.append(serial)
        debit=0
        for u,ns in enumerate(adj):
            for v in ns:
                if u<v:
                    debit+=max(0,(s[u]-s[v])*max(0,U(ds[u]+1,ds[u]-len(ns))-L(ds[v]+1)),
                                 (s[v]-s[u])*max(0,U(ds[v]+1,ds[v]-len(adj[v]))-L(ds[u]+1)))
        assert sum(C)-debit>b+2*a==Q(371,12500)
    print('PASS 3 mixed original-tree ledgers with direct and suppressed ports; floor 371/12500')


if __name__=='__main__':
    cancellation()
    stochastic()
    ledger()
