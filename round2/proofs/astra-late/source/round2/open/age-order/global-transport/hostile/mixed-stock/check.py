"""Mixed-stock hubhood discharge, with rational endpoint checks."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

D = run_path(str(Path(__file__).resolve().parents[1] / 'discharge' / 'check.py'))
G, U, L, profile = (D[k] for k in ('G','U','L','profile'))
a = Q(634,100000)
b = Q(1,22)
g = Q(16694235,10**8)
p = Q(611,1000)


def scalar():
    assert G(2)>Q(2,11) and U(2,0)<p and L(2)>Q(1,4)
    assert D['H']['rescue'](2,1)-D['H']['upper'](2,0)>g
    assert 8*a>b
    for e in range(2,10):
        for d in (0,10):
            s,C,u = profile(e,(d,))
            assert C-s*u>a and C-b>a
    assert g*Q(9,11)+G(2)*Q(9,20)-b>a
    # Analytic large-degree endpoints, including l=0 (no leaf toll).
    assert g*Q(113,114)+G(2)/2-p/2 < 0
    assert (g*Q(113,114)+G(2)/2-p/2)/2+a*Q(109,2)>b
    assert G(2)>Q(2,113) and g*Q(113,114)/2>b
    rows = 0
    for d in range(10,113):
        for c in range(2,d+1):
            assert profile(d,(0,)*c)[1]+a*(c-2)>b, (d,c)
            assert profile(d,(2,)*c)[1]+a*(c-2)>b, (d,c)
            rows += 1
    print(f'PASS 16 terminal endpoints, {rows:,} mixed-stock profiles at both endpoints, and large-degree inequalities')


def make(skeleton, terminal_degrees):
    old,deg = skeleton
    edges=list(old)
    degrees=list(deg)
    for u,e in terminal_degrees:
        edges.append((u,len(degrees)))
        degrees.append(e)
    tree=[[] for d in degrees]
    for u,v in edges:
        tree[u].append(v); tree[v].append(u)
    for u,d in enumerate(degrees):
        assert 0<len(tree[u])<=d
        for _ in range(d-len(tree[u])):
            tree[u].append(len(tree)); tree.append([u])
    return tuple(map(tuple,tree))


def control(tree):
    H=D['H']; old_deg,_,ss=H['profile'](tree)
    cores={u for u,d in enumerate(old_deg) if d>1}
    adj={u:[v for v in tree[u] if v in cores] for u in cores}
    C={}
    for u in cores:
        s,C[u],_=profile(old_deg[u],tuple(sorted(old_deg[v] for v in adj[u])))
        assert s==ss[u]
    internal=[u for u in cores if len(adj[u])>1]
    B=len(internal)
    assert B and min(old_deg[u] for u in internal)>=10
    csum=sum(len(adj[u])-2 for u in internal)
    assert len(cores)-B==2+csum
    for u in internal:
        assert C[u]+a*(len(adj[u])-2)>b
    core_debit=Q(0)
    for u in cores:
        for v in adj[u]:
            if v<=u: continue
            d,e=old_deg[u],old_deg[v]
            lu,lv=d-len(adj[u]),e-len(adj[v])
            loss=max(Q(0),(ss[u]-ss[v])*max(0,U(d+1,lu)-L(e+1)),
                     (ss[v]-ss[u])*max(0,U(e+1,lv)-L(d+1)))
            if len(adj[u])==1:
                assert C[u]-loss>a
            elif len(adj[v])==1:
                assert C[v]-loss>a
            else:
                assert loss<=b
            core_debit+=loss
    lower=sum(C.values())-core_debit
    assert lower>b+2*a
    return float(lower)


if __name__=='__main__':
    scalar()
    fixtures=[make(([(0,1),(1,2)],(10,11,113)),
                   [(0,2)]*9+[(1,3)]*7+[(2,17)]*3),
              make(([(0,1),(0,2),(0,3)],(10,10,12,40)),
                   [(0,2)]*6+[(1,2)]*2+[(2,4)]*8+[(3,10)]*35)]
    for t in fixtures:
        print('PASS mixed physical fixture:',len(t),'vertices; signed lower',control(t))
