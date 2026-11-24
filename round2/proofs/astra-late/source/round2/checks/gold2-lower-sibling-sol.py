#!/usr/bin/env python3
"""Independent replay of the original Gold2 compact sibling-only floor."""
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from math import lcm
from pathlib import Path

GOLD = Path(__file__).resolve().parents[1]/'occupation-quintet-gold2'
GRID = 1 << 20
J = (2,1)
OVER = (0,0)
ORDERS = ('hmys','hmsy','hsmy')
PARENT = dict(m='h',y='m',s='h')
D = dict(zip('hmys',(79,5,5,3)))
C = dict(zip('hmys',(8,4,3,3)))
EXPECTED = (5254160983682774032410464907,3989030756978086473725022297,
            7120570409236883412243761907)


def clip(x): return max(0,-x)


def put(p,i,z): return p[:i]+(z,)+p[i+1:]


def chronology():
    # Exhaust all increasing-parent trees on seven vertices, rejecting the
    # founder as a mother. Direct enumeration uses no prefix templates.
    hits=Counter()
    for tail in product(*(range(i) for i in range(1,7))):
        pa=(-1,*tail); expected=Counter()
        for y in range(1,len(pa)):
            m=pa[y]; h=pa[m]
            if h<=0: continue
            for s in range(h+1,len(pa)):
                if s!=m and pa[s]==h: expected[h,m,y,s]+=1
        seen=Counter()
        for order in ORDERS:
            groups=[[()]]+[[] for _ in order]
            for i in range(1,len(pa)):
                for k in range(3,-1,-1):
                    parent=order.index(PARENT[order[k]]) if k else None
                    for head in groups[k]:
                        if k==0 or head[parent]==pa[i]: groups[k+1].append((*head,i))
            for p in groups[4]: seen[tuple(p[order.index(x)] for x in 'hmys')]+=1
            hits[order]+=len(groups[4])
        assert seen==expected and all(x==1 for x in seen.values())
    assert all(hits[o]>0 for o in ORDERS)
    print('chronology: 720 labelled trees, incidence totals',dict(hits),flush=True)


class Physical:
    def __init__(self,h):
        self.h=h; self.p=h.p
        assert h.sign==1 and h.scale==2**40

    @cache
    def v(self,k,inner,tips):
        return self.p.row(k,range(len(inner)),range(len(inner),len(inner)+len(tips)),inner+tips)

    def h3(self,p): return self.v('h3',(p[1],),tuple(sorted(p[::2]))) or 0

    def ext(self,p,i):
        if i==1: return self.v('s3',(p[1],),tuple(sorted((*p[::2],J))))
        return self.v('p4',(p[i],p[1]),(J,p[2-i]))

    @cache
    def row(self,p):
        # Outward lower sign = -1. The coefficient is in payload ticks
        # per reciprocal rate unit, with no Bellman grid or LCM here.
        v=self.h3(p); a=self.h.lam*v; slopes=[]
        for i,(d,c) in enumerate(p):
            ns=p[max(0,i-1):i]+p[i+1:i+2]; n=c-len(ns)
            post=put(p,i,(d+1,c)); gain=v-self.h3(post)
            covered=self.ext(p,i) is not None
            a+=(d-c+sum(Q(1,x[0]) for x in ns))*gain
            if n and covered: a+=n*self.h.center(d,c)*gain
            slopes.append(gain if n and not covered else 0)
            if d>c:
                born=put(p,i,(d,c+1)); a+=Q(d-c,d)*(v-self.h3(born)-(self.ext(born,i) or 0))
        d,c=p[1]
        for i in (0,2):
            a-=(Q(1,p[i][0])-self.h.center(d,c))*(self.h.g((d+1,c),p[2-i])-self.h.g((d,c),p[2-i]))
        return a,tuple(slopes)

    @cache
    def box(self,p):
        a,b=self.row(p)
        v=a+sum(max(0,z)*Q(n,2) for z,n in zip(b,(p[0][1]-1,p[1][1]-2,p[2][1]-1)))
        assert v<=self.h.golden['verify']['caps']['h3']
        return v,b[0]

    @cache
    def reward(self,p,s):
        if p[0][0]<19: return 0
        box,b=self.box(p)
        r=Q(1,s[0]) if s!=OVER else (Q(1,4) if b>=0 else Q(0))
        improved=box-max(0,b/2)+b*r
        return GRID*(clip(improved)-clip(box))//1


def two_children(r):
    phys=r.h.physical
    edges=((0,1),(1,2),(1,3),(1,4),(2,5),(2,6),(5,7))
    adj={i:{b if a==i else a for a,b in edges if i in (a,b)} for i in range(8)}
    path=(1,2,5); identities=positive=nonadd=0
    for degrees in product((2,4),(19,38,65,79),(3,5),(2,3,4,9),(2,3,4,11),(2,5),(2,9),(2,17)):
        p={i:(d,len(adj[i])) for i,d in enumerate(degrees)}
        if any(d<c for d,c in p.values()): continue
        table=phys.subsets(r.p,adj,p)
        def get(tab,k,ids): return tab.get((k,tuple(sorted(ids))), (0,))[0]
        base=get(table,'h3',path); signed=r.h.lam*base
        for i in path:
            d,c=p[i]; out=adj[i]-set(path)
            gain=phys.subsets(r.p,*phys.grow(adj,p,i,None))
            rate=d-c+sum(Q(1,degrees[j]) for j in adj[i])
            for j in out:
                if any((k,tuple(sorted((*path,j)))) in table for k in ('p4','s3')):
                    rate-=Q(1,degrees[j])-r.h.center(d,c)
            signed+=rate*(base-get(gain,'h3',path))
            if d>c:
                birth=phys.subsets(r.p,*phys.grow(adj,p,i,8))
                signed+=Q(d-c,d)*(base-get(birth,'h3',path)-sum(get(birth,k,(*path,8)) for k in ('p4','s3')))
            if i==2:
                for j,k in ((1,5),(5,1)):
                    signed-=(Q(1,degrees[j])-r.h.center(d,c))*(get(gain,'g',(i,k))-get(table,'g',(i,k)))
        marks=tuple(p[i] for i in path); a,slopes=r.row(marks)
        assert signed==a+sum(b*sum(Q(1,degrees[j]) for j in adj[i]-set(path)) for i,b in zip(path,slopes))
        box,b=r.box(marks)
        z=[max(0,b/2)-b*Q(1,degrees[j]) for j in (3,4)]
        per=[clip(box-u)-clip(box) for u in z]
        assert 0<=sum(per)<=clip(box-sum(z))-clip(box)<=clip(signed)-clip(box)
        for k,j in enumerate((3,4)):
            mark=p[j] if degrees[j]<=D['s'] else OVER
            assert 0<=Q(r.reward(marks,mark),GRID)<=per[k]
        identities+=1; positive+=all(per)
        nonadd+=sum(per)<clip(box-sum(z))-clip(box)
    assert positive and nonadd
    print('two real additional children:',identities,'labelled row identities;',positive,
          'both charges positive;',nonadd,'strict superadditivity cases',flush=True)


def replay(r):
    lo,hi=r.h.prior; den=lcm(hi.denominator,*range(2,80)); inv={n:den//n for n in range(2,80)}
    discount=int(den*hi); totals=[]; terminal={}
    for order in ORDERS:
        nxt={}
        for n in range(4,0,-1):
            roles=order[:n]
            if n==4 and terminal:
                perm=tuple('hmys'.index(x) for x in roles)
                nxt={tuple(p[i] for i in perm):z for p,z in terminal.items()}; continue
            adj=tuple(tuple(j for j,y in enumerate(roles) if PARENT.get(x)==y or PARENT.get(y)==x) for x in roles)
            parent=roles.index(PARENT[order[n]]) if n<4 else -1
            axes=[((OVER,) if x=='s' else ())+tuple((d,c) for d in range(D[x],1,-1)
                  for c in range(min(d,C[x]),max(1,len(adj[i])+(i==0))-1,-1)) for i,x in enumerate(roles)]
            cur={}; positive=0
            for p in product(*axes):
                reward=r.reward(p[:3],p[3]) if n==4 else 0
                positive+=reward>0
                num=reward*den; mass=discount; ports=[]
                for i,(d,c) in enumerate(p):
                    if not d:
                        assert roles[i]=='s' and i!=parent
                        continue
                    g=OVER if roles[i]=='s' and d==3 else (d+1,c)
                    gval=cur.get(put(p,i,g),0); born=put(p,i,(d,c+1))
                    bval=cur.get(born,0)+(nxt.get((*born,J),0) if i==parent else 0)
                    fixed=(d-c)*den+sum(inv[p[j][0]] for j in adj[i] if p[j]!=OVER)
                    width=(c-len(adj[i]))*inv[2]+sum(inv[4] for j in adj[i] if p[j]==OVER)
                    birth=(d-c)*inv[d]
                    num+=fixed*gval+birth*bval; mass+=fixed+birth
                    if width: ports.append((gval,width))
                a,b=num,mass
                for v,w in sorted(ports):
                    if v*b>=a: break
                    a+=v*w; b+=w
                z=a//b
                # Exact worst Bellman residual on the box. No reward slope
                # on any rate: only the unshifted successor v appears here.
                assert num-mass*z+sum(w*min(0,v-z) for v,w in ports)>=0
                cur[p]=z
            nxt=cur
            if n==4: terminal=cur
            print(order,n,len(cur),'positive',positive,flush=True)
        totals.append(lo*Q(nxt[(J,)],GRID*r.h.scale))
    assert totals==[Q(n,71500737144954801917997875200000000) for n in EXPECTED]
    ticks=sum(totals)*r.h.scale
    assert ticks==Q(16363762149897743918379249111,65029541606195200000000)
    print('HOLDS sibling only',ticks,'ticks/core',float(ticks),flush=True)


def main():
    path=GOLD/'host/reader.py'; data=path.read_bytes()
    assert sha256(data).hexdigest()=='6be2e8aa5835ad0f1087179af5d84664d1176ed3ab071ca64f58af8d628fc910'
    host=module_from_spec(spec_from_file_location('gold2_sibling_host',path))
    exec(compile(data,str(path),'exec'),host.__dict__)
    chronology()
    with host.open_host('lower') as h:
        assert h.docs['lower_sha256']=='41565d49c10809eb1f732794822359410af3ed804577476563bef17f27946413'
        r=Physical(h); two_children(r); replay(r)


if __name__=='__main__': main()
