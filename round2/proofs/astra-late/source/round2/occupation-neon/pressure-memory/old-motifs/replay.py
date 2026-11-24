"""Exact signed discounted tagged flow with one persistent side donor."""
from argparse import ArgumentParser
from functools import cache
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'renewal'/'old-flow'))
from rows import Rows, Cert, Q, SPEC, label, up, bud

THETA=Q('1.2415')
NB=(1,0,1)


def extrema(a,d,knobs,sgn):
    """Optimize (a+sum(r*m))/(d+sum(r)) over independent closed boxes."""
    v=a/d
    for m,width in sorted(knobs,key=lambda q:sgn*q[0],reverse=True):
        if sgn*(m-v)<=0:break
        a+=width*m;d+=width;v=a/d
    return v


class Pilot:
    # None: absent; 0: stopped tracking above H (still p <= 1/(H+1));
    # -1: no-memory baseline. A positive state stores the actual donor degree.
    def __init__(self,L,H,f,j,caps,s,donor_cap,ds=Q(0)):
        self.L,self.H,self.f,self.j,self.caps,self.s,self.K,self.ds=L,H,f,j,caps,s,donor_cap,ds
        self.p=SPEC[f]

    @cache
    def affine(self,x):
        if len(x)<len(self.p):return 0,(0,)*len(x),0
        a,c,p=self.L.affine(self.f,x)
        b,d,q=self.H.affine(self.f,x)
        return a-THETA*b,tuple(u-THETA*v for u,v in zip(c,d)),p-THETA*q

    @cache
    def V(self,x,m,sgn):
        n=len(x);j=self.j;p=self.p[:n];s=self.s
        if any(sum(y)>d or y[1]>k for y,(d,k) in zip(x,self.caps)):return Q(0)
        has=j<n and x[j][1]>p.count(j)
        assert m==-1 or has==(m is not None)
        a,c,e=self.affine(x);d=s;knobs=[]

        def action(lo,hi,v):
            nonlocal a,d
            a+=lo*v;d+=lo
            if hi>lo:knobs.append((v,hi-lo))

        for i,y in enumerate(x):
            adj=[k for k in range(n) if p[k]==i]+([p[i]] if i else [])
            g=y[0]+sum(Q(1,sum(x[k])) for k in adj)
            u=self.V(label(x,i,up(y)),m,sgn)
            action(g,g,u)
            k=y[1]-p.count(i)
            if i==j and m not in (-1,None):
                k-=1
                if m:action(Q(1,m),Q(1,m),u+c[i])
                else:action(0,Q(1,self.K+1),u+c[i])
            action(0,Q(k,2),u+c[i])
            if not i:action(0,Q(1,2),u+e)
            if y[0]:
                b=Q(y[0],sum(y));z=label(x,i,bud(y))
                # On our masks the selected slot never creates an extension.
                # Its first extra child stays attached and is degree two at birth.
                born=2 if i==j and m is None else m
                v=self.V(z,born,sgn)
                if n<len(self.p) and self.p[n]==i:
                    assert i!=j
                    v+=self.V(z+(NB,),m,sgn)
                action(b,b,v)
        if m not in (-1,0,None):
            # Donor's parent is the monitored receiver. Its other D−1
            # neighbours contribute between zero and one each. Promotions of
            # the donor do not change its degree or the tracked descriptor.
            lo=Q(1,sum(x[j]))
            action(lo,lo+m-1,self.V(x,m+1 if m<self.K else 0,sgn))
        action(0,self.ds,0)
        return extrema(a,d,knobs,sgn)


def at(d,k):return (d-k-1,k,1)


def cases(f):
    # Couple two successive extension generations (length 2→3→4) on P4/Yt.
    # hV starts at two tags so the apex's existing side-donor identity persists.
    for d in (18,24,64):
        x=(at(d,5),at(3,0))
        if f=='hV':
            cap=((d+1,5),(4,1),(4,1))
            yield 'apex',x,0,cap
            yield 'terminal',(x[0],at(3,1),at(3,0)),2,cap
            yield 'through-terminal',x,2,cap
        elif f=='P4':
            cap=((d+1,5),(4,1),(4,1),(4,1))
            yield 'apex',x,0,cap
            # For the terminal donor, require a full-tag entry with no side
            # child yet; the first promotion creates a real virgin donor.
            yield 'terminal',(x[0],at(3,1),at(3,1),at(3,0)),3,cap
            yield 'through-terminal',x,3,cap
        elif f=='Yt':
            cap=((d+1,5),(4,2),(4,1),(4,1))
            yield 'apex',x,0,cap
            yield 'terminal',(x[0],at(3,2),at(3,0),at(3,0)),3,cap
            yield 'through-terminal',x,3,cap


def check_extrema():
    from itertools import product
    a=Q(-3,7); d=Q(9,5);k=[(Q(1,5),Q(2,3)),(Q(-3,7),Q(1,9)),(Q(-4,3),Q(2,7))]
    v=[(a+sum(r*m for r,(m,_) in zip(rs,k)))/(d+sum(rs))
       for rs in product(*[(0,b) for _,b in k])]
    assert extrema(a,d,k,-1)==min(v) and extrema(a,d,k,1)==max(v)


if __name__=='__main__':
    q=ArgumentParser(description=__doc__)
    q.add_argument('lo',type=Path);q.add_argument('hi',type=Path)
    q.add_argument('--family',action='append',choices=('P4','hV','Yt'))
    q.add_argument('--cap',type=int,default=5)
    q.add_argument('--discount',type=Q,help='fixed discount; default relaxes over the pinned bracket')
    a=q.parse_args();assert a.cap>=2
    L,H=map(lambda x:Rows(Cert(x)),(a.lo,a.hi))
    assert (L.sign,H.sign,L.c.D,H.c.D,L.c.S,H.c.S)==(1,-1,18,18,10**11,10**11)
    s=Q('.565673570322') if a.discount is None else a.discount
    ds=Q('.565775')-s if a.discount is None else Q(0)
    print('discount range',s,s+ds,flush=True)
    check_extrema(); count=0
    for f in a.family or ('P4','Yt','hV'):
        for tag,x,j,caps in cases(f):
            P=Pilot(L,H,f,j,caps,s,a.cap,ds)
            ranges=[]
            for sgn,pick in ((-1,min),(1,max)):
                base=P.V(x,-1,sgn)
                opts=(0,*range(2,a.cap+1)) if tag=='apex' else (None,)
                mem=pick(P.V(x,m,sgn) for m in opts)
                assert sgn*(base-mem)>=0
                if f=='P4' and tag=='through-terminal' and sgn==-1 and a.discount is None:
                    assert base<0<mem
                ranges.append((base,mem));count+=base!=mem
            print(f,tag,'D',sum(x[0]),'raw/memory',[(float(x),float(y)) for x,y in ranges],
                  'inward',[float(abs(x-y)) for x,y in ranges], 'states',P.V.cache_info().currsize,flush=True)
            P.V.cache_clear();P.affine.cache_clear()
    print('strict inward endpoints',count,flush=True)
    if not a.family and a.cap==5 and a.discount is None:
        assert count==51
