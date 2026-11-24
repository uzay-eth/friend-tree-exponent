#!/usr/bin/env python3
"""GOLD2 physical H3 chain and individual additional-child occupation floors."""
import argparse
from hashlib import sha256
from fractions import Fraction as Q
from functools import lru_cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations, permutations, product
from math import lcm
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = spec_from_file_location('gold2_host', HERE.parent/'host/reader.py')
host = module_from_spec(spec)
reader_bytes=Path(spec.origin).read_bytes()
assert sha256(reader_bytes).hexdigest()=='6be2e8aa5835ad0f1087179af5d84664d1176ed3ab071ca64f58af8d628fc910'
exec(compile(reader_bytes,spec.origin,'exec'),host.__dict__)
J = (2,1)
GRID = 2**20
PARENT = {'m':'h','y':'m','s':'h','h':'u'}
REG = ('hmys','hmsy','hsmy')
OVER = (0,0)


def gain(p): return p[0]+1,p[1]
def bud(p): return p[0],p[1]+1
def swap(p,i,q): return p[:i]+(q,)+p[i+1:]


class Rows:
    def __init__(self,h):
        self.h,self.p,self.sig = h,h.p,-h.sign
        self.tags = {'h3':('3m',),'p4':('4i','4i'),'s3':('S3',)}
        self.seams = {k:h.h[k] for k in ('d2','K2','d3m','K3m','d3e','K3e','d4i','K4i','SP4','dS3','KS3','d4e','K4e')}
        assert self.seams['d4e'] <= self.seams['d3e'] <= self.seams['d2']
        assert self.seams['K4e'] <= self.seams['K3e'] <= self.seams['K2']
        assert self.seams['d4i'] <= self.seams['d3e'] and self.seams['K4i'] <= self.seams['K3e']
        assert self.seams['d3m'] <= self.seams['d2'] and self.seams['K3m'] <= self.seams['K2']
        self.affine = lru_cache(maxsize=220000)(self.affine)
        self.h3_value = lru_cache(maxsize=250000)(self.h3_value)

    def support(self,k,inner):
        ids = tuple(self.p.ins[t].get(z,-1) for t,z in zip(self.tags[k],inner,strict=True))
        return min(ids)>=0 and (k!='p4' or tuple(sorted(ids)) in self.p.block['p4'])

    def value(self,k,inner,ends):
        return self.p.row(k,range(len(inner)),range(len(inner),len(inner)+len(ends)),inner+ends) or Q(0)

    def h3(self,p): return self.h3_value(p[1],*sorted((p[0],p[2])))

    def h3_value(self,m,x,y): return self.value('h3',(m,),(x,y))

    @staticmethod
    def ext(p,i):
        if i==1: return 's3',(p[1],),(p[0],p[2],J)
        return 'p4',(p[i],p[1]),(J,p[2-i])

    def affine(self,p):
        """Exact signed original physical row; coordinates are uncovered physical rates."""
        assert self.support('h3',(p[1],))
        b = []
        old = self.h3(p)
        A = -self.h.lam*old
        for i,(d,c) in enumerate(p):
            neighbors = p[max(0,i-1):i]+p[i+1:i+2]
            n = c-len(neighbors)
            assert n>=0 and d>=c
            jumped = self.h3(swap(p,i,gain(p[i])))-old
            rate = d-c+sum(Q(1,q[0]) for q in neighbors)
            covered = self.support(*self.ext(p,i)[:2])
            if covered and n:
                rate += n*self.h.center(d,c)
            b.append(Q(0) if covered or not n else self.sig*jumped)
            A += rate*jumped
            if d>c:
                promoted = swap(p,i,bud(p[i]))
                A += Q(d-c,d)*(self.h3(promoted)-old+self.value(*self.ext(promoted,i)))
        m = p[1]
        for donor,tip in ((p[0],p[2]),(p[2],p[0])):
            A += (Q(1,donor[0])-self.h.center(*m))*(self.h.g(gain(m),tip)-self.h.g(m,tip))
        return self.sig*A,tuple(b)


def registration():
    parents = (None,0,1,1,2,1,3,2,1,5,8,5,6,3,2,9)
    for canonical,orders in (('hmy',('hmy',)),('hmys',REG),('uhmy',('uhmy',))):
        got = []
        for order in orders:
            prefixes = [set() for _ in order]
            for t in range(1,len(parents)):
                for j in range(len(order)-1,0,-1):
                    k = order.index(PARENT[order[j]])
                    prefixes[j].update(p+(t,) for p in prefixes[j-1] if p[k]==parents[t])
                prefixes[0].add((t,))
            got.extend(tuple(p[order.index(i)] for i in canonical) for p in prefixes[-1])
        expected = []
        for group in combinations(range(1,len(parents)),len(canonical)):
            for p in permutations(group):
                tags = dict(zip(canonical,p,strict=True))
                if all(parents[tags[r]]==tags[PARENT[r]] for r in canonical[1:]):
                    expected.append(p)
        assert len(got)==len(set(got)) and set(got)==set(expected) and got
        print('registrations',canonical,len(got),flush=True)


def physical_checks(r):
    """Labelled birth/deletion oracle at the production degree and arity seams."""
    edges=((0,1),(1,2),(0,3),(0,4),(1,5),(2,6),(2,7))
    adj={i:set() for i in range(8)}
    for u,v in edges: adj[u].add(v);adj[v].add(u)
    phys=r.h.physical
    count=0
    for degrees in ((19,5,5,9,3,2,4,5),(23,8,4,87,3,2,4,5),
                    (24,7,9,88,3,2,4,5),(87,19,4,251,3,2,4,5),
                    (88,20,24,252,3,2,4,5),(260,79,105,4000,3,2,4,5)):
        labels={i:(d,len(adj[i])) for i,d in enumerate(degrees)}
        old=phys.subsets(r.p,adj,labels)
        up={i:phys.subsets(r.p,*phys.grow(adj,labels,i,None)) for i in adj}
        promo={i:phys.subsets(r.p,*phys.grow(adj,labels,i,8)) for i in adj if degrees[i]>len(adj[i])}
        def w(rows,k,ids): return rows.get((k,tuple(sorted(ids))),(Q(0),))[0]
        for (k,ids),(val,mid) in old.items():
            if k!='h3': continue
            m,=mid; x,y=sorted(set(ids)-{m}); seq=(x,m,y); subset=set(ids)
            a=-r.h.lam*val
            for i in seq:
                d,c=labels[i]; ext=adj[i]-subset
                jump=w(up[i],k,ids)-val
                a+=(d-c+sum(Q(1,degrees[j]) for j in adj[i]))*jump
                if i in promo:
                    enlarged=(*ids,8)
                    a+=Q(d-c,d)*(w(promo[i],k,ids)-val+sum(w(promo[i],t,enlarged) for t in ('p4','s3')))
                for j in ext:
                    augmented=tuple(sorted((*ids,j)))
                    if any((t,augmented) in old for t in ('p4','s3')):
                        a-=(Q(1,degrees[j])-r.h.center(d,c))*jump
            for donor,keep in ((x,y),(y,x)):
                a+=(Q(1,degrees[donor])-r.h.center(*labels[m]))*(w(up[m],'g',(m,keep))-w(old,'g',(m,keep)))
            A,b=r.affine(tuple(labels[i] for i in seq))
            true=A+sum(z*sum(Q(1,degrees[j]) for j in adj[i]-subset) for i,z in zip(seq,b,strict=True))
            assert true==r.sig*a
            assert true<=r.h.golden['verify']['caps']['h3']
            count+=1
    assert count
    print('literal labelled H3 rows',count,flush=True)


def bellman(total,mass,flex,allow_negative=False):
    """Floor of the minimum quotient, then verify the exact worst signed drift."""
    x,y = total,mass
    for v,w in sorted(flex):
        if v*y>=x: break
        x+=v*w; y+=w
    z=x//y
    assert (allow_negative or z>=0) and total-mass*z+sum(w*min(0,v-z) for v,w in flex)>=0
    return z


class Cohort:
    def __init__(self,rows,a):
        self.r,self.a = rows,a
        self.D = dict(zip('hmysu',(a.mother,a.middle,a.young,a.sibling,a.upstream),strict=True))
        self.C = dict(zip('hmysu',(a.mother_stock,a.middle_stock,a.young_stock,a.sibling,a.upstream),strict=True))
        self.den = lcm(rows.h.prior[1].denominator,*range(2,max(self.D.values())+2))
        self.disc = rows.h.prior[1].numerator*(self.den//rows.h.prior[1].denominator)
        self.inv = {i:self.den//i for i in range(2,max(self.D.values())+2)}
        self.mask = lru_cache(maxsize=220000)(self.mask)

    def mask(self,h,m,y):
        dh,ch=h; dm,cm=m; dy,cy=y
        a=self.a
        if not (a.min_mother <= dh and dh>dy and self.r.support('h3',(m,))): return None
        if a.exclude_strict and (ch,cm,cy)==(2,2,1) and dm<=5 and max(dh,dy)<=65: return None
        if a.outside_threshold and 19<=dh<=79 and ch<=8 and dm<=5 and cm<=4 and dy<=18 and cy<=3:
            return None
        A,b=self.r.affine((h,m,y))
        # Mother's single upstream slot and all her other children have separate ownership.
        B=A+sum(max(0,z)*Q(n,2) for z,n in zip(b,(ch-1,cm-2,cy-1),strict=True))
        assert B<=self.r.h.golden['verify']['caps']['h3'],(h,m,y,B)
        return B,b[0],A,b

    def reward(self,tags,kind):
        x=self.mask(*(tags[r] for r in 'hmy'))
        if x is None: return 0
        B,b,*_=x
        if kind=='box': return GRID*max(0,-B)//1
        who='s' if kind=='sibling' else 'u'
        q=tags[who]
        # The sentinel retains the same actual outside vertex and its ancestry.
        reciprocal=(Q(1,self.D[who]+1) if b>=0 else Q(0)) if q==OVER else Q(1,q[0])
        z=max(0,b/2)-b*reciprocal
        assert z>=0
        return GRID*(max(0,z-B)-max(0,-B))//1

    def sweep(self,orders,kind):
        bank={}; result=[]
        for order in orders:
            nxt={}
            for n in range(len(order),0,-1):
                roles=order[:n]
                adj=tuple(tuple(j for j,t in enumerate(roles) if PARENT.get(t)==r or PARENT.get(r)==t) for r in roles)
                nxt_m=roles.index(PARENT[order[n]]) if n<len(order) else -1
                # Reuse only after reindexing the terminal dictionary into this order.
                if nxt_m<0 and kind=='sibling' and bank:
                    perm=tuple('hmys'.index(r) for r in roles)
                    nxt={tuple(p[j] for j in perm):v for p,v in bank.items()}; continue
                axes=[]
                for i,r in enumerate(roles):
                    mn=len(adj[i])+(i==0)
                    ax=tuple((d,c) for d in range(self.D[r],1,-1)
                             for c in range(min(d,self.C[r]),max(1,mn)-1,-1))
                    axes.append((OVER,)+ax if r in 'su' else ax)
                cur={}; positive=0
                for p in product(*axes):
                    tags=dict(zip(roles,p)) if nxt_m<0 else None
                    reward=self.reward(tags,kind) if tags is not None and kind!='physical' else 0
                    affine=self.mask(*p) if tags is not None and kind=='physical' else None
                    if affine and affine[2]+sum(min(0,b)*Q(c-len(ns),2) for b,(_,c),ns in zip(affine[3],p,adj))>=0:
                        affine=None
                    positive+=reward>0
                    positive+=affine is not None
                    total,mass=reward*self.den,self.disc
                    flex=[]
                    shifted=[]
                    for i,(d,c) in enumerate(p):
                        r=roles[i]
                        visible=adj[i]
                        if p[i]==OVER:
                            if i==nxt_m:
                                # Absorbing upstream still promotes. D>=cut+1 and
                                # c/D<=1; zero is a valid lower birth-rate bound.
                                assert r=='u'
                            continue
                        after_gain=OVER if r in 'su' and d==self.D[r] else (d+1,c)
                        vg=cur.get(swap(p,i,after_gain),0)
                        lo=(d-c)*self.den+sum(self.inv[p[j][0]] for j in visible if p[j]!=OVER)
                        width=(c-len(visible))*self.inv[2]+sum(self.inv[self.D[roles[j]]+1] for j in visible if p[j]==OVER)
                        post=swap(p,i,bud(p[i])); birth=(d-c)*self.inv[d]
                        vb=cur.get(post,0)
                        if i==nxt_m:
                            vb+=nxt.get(post+(J,),0)
                        total+=lo*vg+birth*vb; mass+=lo+birth
                        if width:
                            flex.append((vg,width))
                            if affine: shifted.append((vg-(-GRID*affine[3][i]//1),width))
                    z=bellman(total,mass,flex)
                    if affine:
                        # max(0,-P) bounds both zero and the affine branch.
                        # Work in exact integer rate units, rounding branch
                        # reward contributions downward, including the slope.
                        v=total+(-GRID*affine[2]*self.den//1)
                        z=max(z,bellman(v,mass,shifted,allow_negative=True))
                    cur[p]=z
                print(kind,order,n,'states',len(cur),'positive rewards',positive,flush=True)
                nxt=cur
                if nxt_m<0 and kind=='sibling':
                    assert roles=='hmys'
                    bank=cur
            v=self.r.h.prior[0]*Q(nxt[(J,)],GRID*self.r.h.scale)
            result.append(v)
            print('order floor',kind,order,v,float(v),flush=True)
        return sum(result,Q(0))


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--side',choices=('lower','upper'),default='lower')
    for k,v in (('mother',37),('middle',5),('young',5),('sibling',2),('upstream',2),
                ('mother-stock',7),('middle-stock',4),('young-stock',3),('min-mother',19)):
        cli.add_argument('--'+k,type=int,default=v)
    cli.add_argument('--kind',choices=('box','physical','sibling','upstream','joint'),default='joint')
    cli.add_argument('--exclude-strict',action='store_true')
    cli.add_argument('--outside-threshold',action='store_true',help='exclude the companion D19–79/c8, D5/c4, D18/c3 H3 floor')
    a=cli.parse_args()
    assert __debug__ and min(a.mother,a.middle,a.young,a.sibling,a.upstream)>=2
    registration()
    with host.open_host(a.side) as h:
        r=Rows(h)
        print('GOLD2',a,h.lam,h.prior,h.golden['verify']['sha256'],r.seams,flush=True)
        physical_checks(r)
        c=Cohort(r,a); out=[]
        for k,orders in (('box',('hmy',)),('physical',('hmy',)),('sibling',REG),('upstream',('uhmy',))):
            if a.kind not in (k,'joint'): continue
            if a.kind=='joint' and k=='physical': continue
            v=c.sweep(orders,k); out.append(v)
            print('floor',a.side,k,v,float(v),flush=True)
        v=sum(out,Q(0)); print('TOTAL',a.side,a.kind,v,float(v),flush=True)


if __name__=='__main__': main()
