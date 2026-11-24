"""Low component message test, retaining terminal templates and high ports."""
from fractions import Fraction as Q
from functools import cache
from itertools import combinations_with_replacement
from pathlib import Path
from runpy import run_path

M = run_path(str(Path(__file__).resolve().parents[1]/'mixed-stock'/'check.py'))
G,U,L,profile,a = (M[k] for k in ('G','U','L','profile','a'))
T = 10**8
HIGH = 10


def ceil(x):
    return -((-x)//1)


@cache
def point(d, ds):
    s,C,_ = profile(d,ds)
    return (T*s)//1,(T*C)//1


@cache
def setup(d,ds):
    lo,C = point(d,tuple(x if x!=HIGH else 0 for x in ds))
    hi,E = point(d,ds)
    return lo,hi,min(C,E),d-len(ds)


@cache
def slopes(d,l,e,j):
    return tuple(ceil(T*max(0,x)) for x in (U(d+1,l)-L(e+1),U(e+1,j)-L(d+1)))


def raw(d,l,x,e,j,y,z):
    f,g=slopes(d,l,e,j)
    return max((x-y+1)*f,(z-x+1)*g,0)


def debit(d,l,x,e,ds):
    y,z,_,j=setup(e,ds)
    return (raw(d,l,x,e,j,y,z)+T-1)//T


@cache
def high_debit(d, l, x, y):
    # The high endpoint has s<=1/11. Keep its outgoing send
    # exposure bounded jointly with its s, instead of multiplying maxima.
    return max(ceil(Q((y+1)*ceil(T*U(d+1,l)),T)),
               ceil(T*max(0,min(Q(12,1000), Q(1,11)-Q(x,T)))))


def run(cap):
    degrees=tuple(range(2,cap+1))+(HIGH,)
    rows=[]
    incoming={}
    worst_high=(0,None)
    worst_low=(0,None)
    worst_root=(0,None)
    A=int(a*T)
    for d in degrees[:-1]:
        for c in range(1,d+1):
            for ds in combinations_with_replacement(degrees,c):
                x,y,C,l=setup(d,ds)
                cost=high_debit(d,l,x,y) if HIGH in ds else 0
                msgs={e:(int((a+Q(1,22))*T)-cost if e==HIGH else setup(e,(d,))[2]
                         -max(debit(d,l,x,e,(d,)),debit(d,l,y,e,(d,)))) for e in set(ds)}
                total=C+sum(msgs[e] for e in ds)
                worst_root=min(worst_root,(total,(d,ds))) if total<0 else worst_root
                if HIGH in ds:
                    val=total-msgs[HIGH]-cost-A
                    if val<worst_high[0]: worst_high=val,(d,ds)
                for e in set(ds)-{HIGH}:
                    b=(d,ds,total-msgs[e]); rows.append((e,b))
                    k=d,e,l
                    u,v=incoming.get(k,(x,y))
                    incoming[k]=min(x,u),max(y,v)
    for e,(d,ds,branch) in rows:
        if ds==(e,): continue
        lo,hi,term,k=setup(d,(e,))
        x,y,_,j=setup(d,ds)
        for l in range(e):
            span=incoming.get((e,d,l))
            if span is None: continue
            u,v=span
            # Difference between convex, piecewise affine rounded debit bounds.
            # Use all integer neighbours of each line crossing.
            pts={u,v}
            for low,high,stock in ((lo,hi,k),(x,y,j)):
                f,g=slopes(e,l,d,stock)
                if f+g:
                    p=((high+1)*g+(low-1)*f)//(f+g)
                    pts.update(i for i in (p,p+1) if u<=i<=v)
            val=branch-term-max((raw(e,l,t,d,j,x,y)-raw(e,l,t,d,k,lo,hi)+T-1)//T for t in pts)
            if val<worst_low[0]: worst_low=val,(e,l,d,ds)
    print(cap,len(rows),'high',worst_high,'low',worst_low,'root',worst_root,flush=True)


@cache
def endpoint(d, ds):
    s,C,u=profile(d,ds)
    return d,d-len(ds),s,C,u,L(d+1)


def pair(A,B):
    _,_,s,_,u,l=A
    _,_,t,_,v,j=B
    return max(0,(s-t)*max(0,u-j),(t-s)*max(0,v-l))


def port(A):
    # Keep outgoing exposure joint in high s and high U.
    _,_,s,_,u,l=A
    return max(s*u, Q(12,1000)*(1-l))


def variants(d,ds):
    for r in (0,10):
        yield endpoint(d,tuple(r if x==HIGH else x for x in ds)), Q(1,r) if r else Q(0)


def chains():
    # For d>=41 use 1/(2(d+1)); for other degrees use actual stock upper.
    for d in range(10,42):
        for l in range(d):
            assert Q(d+l,2*d*(d+1))*U(d+1,l)<Q(12,1000)
    assert Q(1,84)<Q(12,1000)
    ends=list(variants(2,(2,HIGH)))
    mid=endpoint(2,(2,2))
    mins=[]
    for k in (1,2,3):
        vals=[]
        for A,r in (variants(2,(HIGH,HIGH)) if k==1 else ends):
            if k==1:
                # two independent reciprocals; include the off-diagonal corner.
                for B,h,j in ((A,r,r),(endpoint(2,(0,HIGH)),Q(0),Q(1,10))):
                    vals.append(B[3]-2*port(B))
                continue
            for B,h in ends:
                vs=[A]+[mid]*(k-2)+[B]
                vals.append(sum(v[3] for v in vs)-sum(pair(u,v) for u,v in zip(vs,vs[1:]))-port(A)-port(B))
        mins.append(min(vals))
    print('connector endpoint lower bounds, lengths 1/2/3:',*[float(x) for x in mins],flush=True)
    assert all(v>-Q(1,22) for v in mins)
    assert mid[3]>0 and pair(mid,mid)==0
    # Boundary chains ending in a low terminal.
    vals=[]
    for e in range(2,10):
        term=endpoint(e,(2,))
        for A,r in variants(2,(e,HIGH)):
            vals.append((A[3]+term[3]-port(A)-pair(A,term),e,1))
        for A,r in ends:
            B=endpoint(2,(2,e))
            for vs in ([A,B,term],[A,mid,B,term]):
                vals.append((sum(v[3] for v in vs)-port(A)-sum(pair(u,v) for u,v in zip(vs,vs[1:])),e,len(vs)-1))
    worst=min(vals)
    print('pendant endpoint lower bound:', float(worst[0]), 'terminal, length',worst[1:],flush=True)
    assert worst[0]>a


def optimize_paths():
    """All paths of low cores with exactly two core slots; arbitrary stocks."""
    ds = range(2,10)
    scale=10**10
    def down(x): return (x*scale)//1
    def up(x): return ceil(x*scale)
    states=[(i,d,j) for i in (*ds,0,10) for d in ds for j in (*ds,0,10)]
    pts={s:endpoint(s[1],(s[0],s[2])) for s in states}
    weights={s:down(v[3]) for s,v in pts.items()}
    edges={s: [(t,weights[t]-up(pair(v,pts[t])))
               for j in (*ds,0,10) if (t:=(s[1],s[2],j)) in pts] for s,v in pts.items()}
    # A prescribed high boundary has no owned core credit in these messages.
    dist={s:weights[s]-up(port(v)) for s,v in pts.items() if s[0] not in ds}
    for step in range(1,len(states)+1):
        nxt=dict(dist)
        for s,x in dist.items():
            for t,w in edges[s]:
                z=x+w
                if z<nxt.get(t,float('inf')): nxt[t]=z
        if nxt==dist: break
        dist=nxt
    else:
        raise AssertionError('No stabilized path lower potential')
    conn=min((x-up(port(pts[s])),s) for s,x in dist.items() if s[2] not in ds)
    pendant=min((x+down(endpoint(s[2],(s[1],))[3])-up(pair(pts[s],endpoint(s[2],(s[1],)))),s)
                for s,x in dist.items() if s[2] in ds)
    print(f'PASS stocked paths: {len(states)} states, closure step {step}; '
          f'connector {conn}, pendant {pendant}, scale {scale}',flush=True)
    assert Q(conn[0],scale)>-Q(1,22) and Q(pendant[0],scale)>a


def singleton():
    rows=[]
    for d in range(2,10):
        for c in range(2,d+1):
            for r in (0,10):
                A=endpoint(d,(r,)*c)
                rows.append((A[3]-c*port(A)+Q(c-1,22),d,c,r))
    worst=min(rows)
    assert worst[0]>Q(22,1000)
    print(f'PASS {len(rows)} singleton junction endpoints; worst excess {float(worst[0]):.9f}',flush=True)


def fixture():
    # Two high internal hubs, a stocked connector, a high-port-only junction,
    # high terminals, and a stocked pendant terminating at a low core.
    tree=M['make'](([(0,2),(2,3),(3,1),(1,4),(4,5),(4,6),(0,7),(7,8),(8,9)],
                    (10,18,3,6,5,14,24,7,4,2)),[(0,2),(1,3)])
    H=M['D']['H']; deg,_,s=H['profile'](tree)
    cores=[i for i,d in enumerate(deg) if d>1]
    A={i:endpoint(deg[i],tuple(deg[j] for j in tree[i] if deg[j]>1)) for i in cores}
    assert all(s[i]==A[i][2] for i in cores)
    lower=sum(v[3] for v in A.values())-sum(pair(A[i],A[j]) for i in cores for j in tree[i] if j in A and j>i)
    assert lower>Q(1,22)+2*a
    print(f'PASS mixed physical fixture: {len(tree)} vertices, numerator lower {float(lower):.9f}',flush=True)


if __name__=='__main__':
    import sys
    if len(sys.argv)>1: run(int(sys.argv[1]))
    else:
        chains()
        optimize_paths()
        singleton()
        fixture()
