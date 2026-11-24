#!/usr/bin/env python3
"""Sparse Claim 6 singleton census in inward signed-drift convention."""
from fractions import Fraction as Q
import hashlib
import mmap
from math import comb
from pathlib import Path
import struct

SOURCE = 'b04cc50c1e391e0420eb7b006ac781d484c68543'
SIZE = 1588521616
PINS = {'lower': (2,56570237,'81dc1f250d49083716e0998dcb63b615bd55401c4e35ea388ce098d29b4a3986'),
        'upper': (1,56574828,'8a36c3afdd2a8f4bbc03461a306af04594c46ada6ef037646261da4d364d9383')}
FIELDS = ('magic mode P Q S S_r alpha_lo_num f_p0 f_p1 f_p2 DF d2 K2 d3m K3m d3e K3e '
          'd4i K4i SP4 dS3 KS3 d4e K4e d5i K5i SP5 d5o K5o d5c K5c SB5 dS4 KS4 d5e K5e Q_K nknots').split()


class Snapshot:
    def __init__(self,path,side,*,verify=True):
        self.path=Path(path);self.side=side
        mode,p,digest=PINS[side]
        with self.path.open('rb') as f:
            import os
            self.stat=os.fstat(f.fileno())
            if self.stat.st_size!=SIZE: raise ValueError('wrong production payload size')
            raw=struct.unpack('<Q63q',f.read(512));h=dict(zip(FIELDS,raw))
            if (h['magic'],h['mode'],h['P'],h['Q'],h['S'],h['S_r'],h['DF'])!=(0x4654325F55355631,mode,p,10**8,1<<40,1<<30,2000):
                raise ValueError('wrong pinned production header')
            if verify:
                f.seek(0)
                if hashlib.file_digest(f,'sha256').hexdigest()!=digest: raise ValueError('wrong payload SHA-256')
            self.data=mmap.mmap(f.fileno(),0,access=mmap.ACCESS_READ)
        self.authenticated=verify; self.header=h; self.lam=Q(p,h['Q']);self.scale=h['S'];self.sign=1 if mode==2 else -1
        d2,k=h['d2'],h['K2'];self.DF=h['DF'];self.d2=d2;self.K2=k
        self.inner={ (d,c):i for i,(d,c) in enumerate((d,c) for d in range(2,d2) for c in range(1,min(d,k)+1)) }
        count=len(self.inner);self.tail_base=count+max(0,d2-1-k);self.inf=self.tail_base+k+1
        self.off_f=count*(1+h['nknots']);self.off_g=self.off_f+self.fid(self.DF,0)
        self.newborn_f=self.f(2,1);self.ginf=self.g(self.inf,0)
        self.tail_edges=tuple(self.g(i,0) for i in range(self.tail_base,self.inf))
        self.check_unchanged()

    def check_unchanged(self):
        s=self.path.stat();fields=('st_dev','st_ino','st_size','st_mtime_ns','st_ctime_ns')
        if any(getattr(s,f)!=getattr(self.stat,f) for f in fields): raise ValueError('payload modified/replaced')

    def close(self): self.data.close()
    def __enter__(self): return self
    def __exit__(self,*_): self.close()

    @staticmethod
    def fid(d,c): return (d-2)*(d+3)//2+c
    def item(self,i): return struct.unpack_from('<q',self.data,512+8*i)[0]
    def g(self,i,j):
        a,b=sorted((i,j)); return self.item(self.off_g+b*(b+1)//2+a)
    def f(self,d,c):
        if not 2<=d or not 0<=c<=d: raise ValueError('nonphysical census')
        if d<self.DF:return self.item(self.off_f+self.fid(d,c))
        h=self.header; q=Q(d-c,d)
        return h['f_p0']+h['f_p1']*q if self.sign<0 else q*(h['f_p0']-h['f_p1']*q-Q(h['f_p2'],d))

    def newborn_edge(self,d,c_after):
        if not 1<=c_after<=d: raise ValueError('nonphysical promotion')
        if d<self.d2:
            i=self.inner[d,c_after] if c_after<=self.K2 else len(self.inner)+d-self.K2-1
            return self.g(i,0)
        tail=self.tail_edges[min(c_after,self.K2+1)-1]
        return self.ginf+Q(self.d2,d)*(tail-self.ginf)

    def center(self,d,c):
        i=self.inner.get((d,c)); return None if i is None else Q(self.item(i),self.header['S_r'])

    def components(self,d,c):
        l=d-c; old=self.f(d,c); df=self.f(d+1,c)-old; center=self.center(d,c)
        birth=Q(0) if not l else Q(l,d)*(self.f(d,c+1)-old+self.newborn_f+self.newborn_edge(d,c+1))
        return old,df,birth,center

    def row_ticks(self,d,c):
        old,df,birth,center=self.components(d,c)
        return (self.sign*((d-c+(0 if center is None else c*center))*df+birth-self.lam*old),
                self.sign*df if center is None else Q(0))

    def row(self,d,c):
        """(A,B) s.t. original inward signed singleton residual = A+B*rho."""
        a,b=self.row_ticks(d,c);return Q(a,self.scale),Q(b,self.scale)

    def corners(self,d,c):
        a,b=self.row(d,c);return a,a+b*Q(c,2)

    def tail_floor(self,c,degree=None):
        """All D>=degree at fixed core degree, exact one-variable Bernstein."""
        cutoff=max(self.DF,c,degree or self.DF)
        X=Poly([0,1]); one=Poly([1]); plus=one+X; q=one-c*X
        h=self.header
        a,b,m,e=(0,h['f_p0'],-h['f_p1'],-h['f_p2']) if self.sign>0 else (h['f_p0'],h['f_p1'],0,0)
        w=a+b*q+m*q*q+e*q*X
        jump=b*(one-q)*plus+m*(one-q)*(2*q+(one+q)*X)+e*((one-2*q)*X-q*X*X)
        bud=-b*X+m*(-2*q*X+X*X)-e*X*X
        edge=self.ginf+X*(self.d2*(self.tail_edges[min(c+1,self.K2+1)-1]-self.ginf))
        xmin=Q(1,cutoff)
        vals=[]
        for j,r in enumerate((q,q+Q(c,2)*X)):
            residual=self.sign*(r*jump+(q*(bud+self.newborn_f+edge)-self.lam*w)*plus*plus)
            for d in (cutoff,cutoff+1,3*cutoff+17):
                x=Q(1,d)
                assert residual(x)==self.corners(d,c)[j]*self.scale*(1+x)**2
            floor=min(residual.bernstein(xmin)); vals.append(floor/((1+xmin)**2 if floor>=0 else 1)/self.scale)
        return min(vals)

    def stock_floor(self,c):
        """Return (finite minimum, argmin degree, analytic tail floor, union floor)."""
        if not 0<=c<self.DF: raise ValueError('stock-floor core degree beyond finite scan')
        f,d=min((min(self.corners(d,c)),d) for d in range(max(2,c),self.DF))
        t=self.tail_floor(c)
        return f,d,t,min(f,t)


class Poly(tuple):
    def __new__(cls,cs): return tuple.__new__(cls,tuple(map(Q,cs)))
    @staticmethod
    def lift(x): return x if isinstance(x,Poly) else Poly([x])
    def __add__(self,x):
        b=self.lift(x); n=max(len(self),len(b))
        return Poly((self[i] if i<len(self) else 0)+(b[i] if i<len(b) else 0) for i in range(n))
    __radd__=__add__
    def __neg__(self): return Poly(-v for v in self)
    def __sub__(self,x): return self+-self.lift(x)
    def __rsub__(self,x): return self.lift(x)+-self
    def __mul__(self,x):
        b=self.lift(x); out=[Q(0)]*(len(self)+len(b)-1)
        for i,a in enumerate(self):
            for j,v in enumerate(b): out[i+j]+=a*v
        return Poly(out)
    __rmul__=__mul__
    def __call__(self,x):
        y=Q(0)
        for a in reversed(self): y=y*x+a
        return y
    def bernstein(self,xmax):
        n=len(self)-1
        return tuple(sum(self[j]*xmax**j*Q(comb(i,j),comb(n,j)) for j in range(i+1)) for i in range(n+1))
