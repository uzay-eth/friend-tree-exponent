"""Exact charged-slot curl, gauge boundary and stopped Bellman comparison."""
from argparse import ArgumentParser
from functools import cache
from itertools import product
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[2]/'renewal'/'old-flow'))
from rows import Rows, Cert, Q, SPEC, label, up, bud

HI=Q('.565775')
THETA=Q('1.2415')


def gain(x,i):
    return label(x,i,up(x[i]))


class Face:
    def __init__(self,lo,hi,f,t,cap,full=False):
        self.lo,self.hi,self.f,self.t=lo,hi,f,t
        self.p=SPEC[f]
        self.j=len(self.p)-1
        self.cap=cap
        self.full=full

    @cache
    def source(self,x):
        a,c,e=self.lo.affine(self.f,x)
        if not self.t:return a,c,e
        b,d,f=self.hi.affine(self.f,x)
        return a-self.t*b,tuple(u-self.t*v for u,v in zip(c,d)),e-self.t*f

    @cache
    def z(self,x):
        """Anchor at the lowest legal stock with the same terminal arity."""
        if self.full:return -self.weight(x)
        j=self.j
        l,k,h=x[j]
        if l==max(0,2-k-h):return Q(0)
        y=label(x,j,(l-1,k,h))
        return self.z(y)-self.source(y)[1][j]

    def live(self,x):
        return all(sum(y)<=d and y[1]<=k for y,(d,k) in zip(x,self.cap))

    def row(self,x,gauged):
        a,c,e=self.source(x)
        x0=self.z(x) if gauged else Q(0)
        independent={}
        control=[]
        for i,y in enumerate(x):
            adj=[j for j in range(len(x)) if self.p[j]==i]+([self.p[i]] if i else [])
            g=gain(x,i)
            lam=y[0]+sum(Q(1,sum(x[j])) for j in adj)
            dz=self.z(g)-x0 if gauged else 0
            a+=lam*dz
            independent[g]=independent.get(g,0)+lam
            if y[0]:
                b=Q(y[0],sum(y));u=label(x,i,bud(y))
                independent[u]=independent.get(u,0)+b
                if gauged:a+=b*(self.z(u)-x0)
            m=y[1]-self.p.count(i)
            if m:control.append((g,Q(m,2),c[i]+dz))
            if i==0:control.append((g,Q(1,2),e+dz))
            if gauged and i==self.j and m:assert c[i]+dz==0
        if gauged and self.full and self.f!='hV':assert all(e==0 for _,_,e in control)
        return a-HI*x0,independent,control

    @cache
    def V(self,x,gauged,sign,boundary=True):
        if not self.live(x):
            return -self.z(x) if gauged and boundary else Q(0)
        a,known,controls=self.row(x,gauged)
        a+=sum(r*self.V(y,gauged,sign,boundary) for y,r in known.items())
        d=HI+sum(known.values())
        vals=[(a+sum(r*(e+self.V(y,gauged,sign,boundary)) for r,(y,_,e) in zip(cs,controls)))/
              (d+sum(cs)) for cs in product(*((0,b) for _,b,_ in controls))]
        return (min if sign<0 else max)(vals)

    def width(self,x,gauged):
        return sum(b*abs(e) for _,b,e in self.row(x,gauged)[2])

    def omega(self,x,i):
        _,c,e=self.source(x)
        return c[i]+(e if i==0 and not x[i][1]-self.p.count(i) else 0)

    def curl(self,x,i,j):
        return self.omega(x,i)+self.omega(gain(x,i),j)-self.omega(x,j)-self.omega(gain(x,j),i)

    def shared_curl(self,x,i,j):
        a,b=gain(x,i),gain(x,j)
        s=self.source
        return s(x)[1][i]+s(a)[1][j]-s(x)[1][j]-s(b)[1][i]

    def weight(self,x):
        l=self.lo.sign*self.lo.value(self.f,x)/self.lo.c.S
        h=self.hi.sign*self.hi.value(self.f,x)/self.hi.c.S
        return l-self.t*h

    def full_gradient(self,x,i):
        return self.weight(gain(x,i))-self.weight(x)

    def full_running(self,x):
        """Independent reconstruction after the full weight cancellation."""
        def one(R):
            # Import physical subembedding enumeration as used by the extractor.
            from rows import embedding
            n=len(x);a=(HI-R.c.lam)*R.value(self.f,x)
            for d in range(n):
                for g,keep in embedding(self.p):
                    if len(keep)!=n-1 or d in keep:continue
                    j,=[u for u in keep if self.p[u]==d or self.p[d]==u]
                    y=tuple(x[u] for u in keep);parent=self.p[j]==d
                    if R.live(g,y):
                        a+=(Q(1,sum(x[d]))-R.center(x[j],parent))*R.transfer(g,y,keep.index(j),parent)
            return a*R.sign/R.c.S
        return one(self.lo)-self.t*one(self.hi)


def at(d,k):return (d-k-1,k,1)


def seeds(f):
    for d in (18,24,32,64,65):
        if f=='hV':
            for j in (2,3,4,6,8):
                for k in (1,2):
                    yield (at(d,5),at(j+1,1),at(max(j,k+1),k))
        elif f=='P4':
            for j in (3,4,6,8):yield (at(d,5),at(4,1),at(4,1),at(j,1))
        else:
            for j in (3,4,6,8):yield (at(d,5),at(4,2),at(4,1),at(j,1))


def rank(lo,hi,f,t):
    xs=tuple(seeds(f));n=len(xs[0]);F=Face(lo,hi,f,t,())
    widths=[(F.width(x,False),F.width(x,True),x) for x in xs]
    curls=[(abs(F.curl(x,i,n-1)),x,i) for x in xs for i in range(n-1)]
    count=sum(bool(z) for z,_,_ in curls)
    active=[(x,i) for x in xs for i in range(n-1)
            if x[i][1]>F.p.count(i) and x[n-1][1]>F.p.count(n-1)]
    cc=[F.shared_curl(x,i,n-1) for x,i in active]
    nz=sum(bool(v) for v in cc)
    # A scalar gauge cannot change the circulation. The best correction to
    # all four oriented edges leaves one residual at least |curl|/4.
    if nz:
        q,(x,i)=max(zip(map(abs,cc),active))
        natural = [F.source(y)[1][j]+F.z(gain(y,j))-F.z(y)
                   for y,j in ((x,i),(gain(x,i),n-1),(x,n-1),(gain(x,n-1),i))]
        assert abs(natural[0]+natural[1]-natural[2]-natural[3])==q
        assert max(map(abs,natural))>=q/4
        print(f,'protected cycle coefficient floor',str(q/4),float(q/4),flush=True)
    for x in xs:
        assert F.source(x)[1][-1]==(0 if f=='hV' and F.lo.moved(f,x,n-1,False)
                                    else F.full_gradient(x,n-1))
        if f!='hV':
            a,c,e=F.source(x)
            assert e==F.full_gradient(x,0)
            assert all(c[i]==F.full_gradient(x,i) for i in range(n) if x[i][1]>F.p.count(i))
    if f!='hV':
        G=Face(lo,hi,f,t,(),full=True)
        running=[G.full_running(x) for x in xs]
        assert running==[G.row(x,True)[0] for x in xs]
        print(f,'full gauge running range',float(min(running)),float(max(running)),flush=True)
    high=max(curls)
    before=sum(a for a,_,_ in widths);after=sum(b for _,b,_ in widths)
    print(f,'theta',str(t),'faces',len(xs),'static mean widths',float(before/len(xs)),float(after/len(xs)),
          'all-slot curls',count,'/',len(curls),'max',high[0],high[1:],
          'shared active child curls',nz,'/',len(cc),'max active',str(max(map(abs,cc),default=0)),flush=True)
    # Choose the biggest reduction in pointwise box width; no occupation is inferred from this ranking.
    _,x=max((a-b,x) for a,b,x in widths)
    cap=tuple((sum(y)+1,y[1]+1) for y in x)
    B=Face(lo,hi,f,t,cap)
    old=[B.V(x,False,s) for s in (-1,1)]
    exact=[B.V(x,True,s)+B.z(x) for s in (-1,1)]
    bad=[B.V(x,True,s,False)+B.z(x) for s in (-1,1)]
    assert old==exact
    if f!='hV':
        A=Face(lo,hi,f,t,cap,full=True)
        assert old==[A.V(x,True,s)+A.z(x) for s in (-1,1)]
    print(f,'ranked matched/stopped',x,'raw and gauge',*[float(q) for q in old],
          'dropped exit falsely shifts',*[float(u-v) for u,v in zip(bad,old)],flush=True)
    return nz, old!=bad


def verify_birth(lo,hi):
    # In a positive hierarchy the newly produced full tag carries -z at its
    # post-promotion state; putting +b z in the preceding reward cancels it.
    for f in ('hV','P4','Yt'):
        F=Face(lo,hi,f,THETA,())
        A=Face(lo,hi,f,THETA,(),full=True)
        for x in seeds(f):
            y=label(x,F.j,(1,0,1))
            assert F.z(y)==0
            i=F.p[F.j]
            pre=label(y[:-1],i,(y[i][0]+1,y[i][1]-1,y[i][2]))
            assert label(pre,i,bud(pre[i]))+((1,0,1),)==y
            b=Q(pre[i][0],sum(pre[i]))
            assert b*A.z(y)!=0 or not A.weight(y)
            # The prefix action sees both the new tagged value and reward.
            raw=Q(17,10**6);a=Q(3,10**6)
            assert a+b*raw==a+b*A.z(y)+b*(raw-A.z(y))


if __name__=='__main__':
    p=ArgumentParser(description=__doc__)
    p.add_argument('lo',type=Path);p.add_argument('hi',type=Path)
    a=p.parse_args()
    L,H=map(lambda x:Rows(Cert(x)),(a.lo,a.hi))
    assert (L.c.mode,H.c.mode,L.c.D,H.c.D,L.c.S,H.c.S)==(2,1,18,18,10**11,10**11)
    for lo,hi,t in ((L,H,Q(0)),(H,L,Q(0)),(L,H,THETA)):
        print('orientation',lo.sign,hi.sign,'theta',t,flush=True)
        results=[rank(lo,hi,f,t) for f in ('hV','P4','Yt')]
        assert [x[0] for x in results[1:]]==[0,0]
        assert any(b for _,b in results)
    verify_birth(L,H)
