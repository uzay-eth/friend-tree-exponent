"""Physical charged motifs of the pinned Neon D18 potentials."""
from collections import defaultdict
from fractions import Fraction as Q
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'higher'))
from extract import Cert, Extract, NB, SPEC, embedding, shape, up, bud, label


class Charges(Extract):
    """The same transfers as GMP, evaluated at actual integer degrees."""
    def supported(self,f,x):
        if f=='hV':return 1<=x[1][1]<=self.c.K
        if f=='hC':return 2<=x[0][1]<=self.c.K
        return super().supported(f,x)

    @lru_cache(maxsize=100000)
    def weight(self,f,x):
        if not self.supported(f,x):return Q(0)
        if f=='f':return Q(int(self.c.f[self.c.vi(*x[0])]))
        c=self.c; a,*b=x
        def raw(ss):
            p,*q=ss
            if f=='g':return c.g[c.state[p],c.state[q[0]]]
            if f=='hV':return c.hvtab[c.par[p],c.hv[q[0]],c.child[q[1]]] if q[0] in c.hv else 0
            if f=='hC':return c.hctab[c.hc[p],c.tp(*q)] if p in c.hc else 0
            if f=='P4':return c.p4tab[c.par[p],c.p4[q[0]],c.p4[q[1]],c.child[q[2]]]
            if f=='Yt':return c.yttab[c.par[p],c.yt[q[0]],c.tp(*q[1:])]
            if f=='Yb':return c.ybtab[c.ybu[p],c.ybw[q[1]],c.child[q[0]],c.child[q[2]]]
            return c.s3tab[c.s3[p],c.tt(*q)]
        slots=[c.embed(y) if not ((f=='P4' and i in (1,2)) or (f=='Yt' and i==1) or
                                  (f=='Yb' and i in (0,2)) or (f=='S3' and i==0)) else [(y,Q(1))]
               for i,y in enumerate(x)]
        return sum(int(raw(ss))*prod(ww) for rr in product(*slots) for ss,ww in (zip(*rr),))

    def transferred(self,f,x,i,upward):
        c=self.c;y=x[i];d=sum(y); _,k,h=y
        if d<c.D:
            if not self.route(f,x,i,upward):return Q(0),Q(0)
            r=c.rpar if upward else c.rch
            return Q(int(r[c.state[y]]),c.Sr),self.weight(f,label(x,i,up(y)))-self.weight(f,x)
        if k>c.K or f not in ('f','g'):return Q(0),Q(0)
        j=(1-h)*(c.K+2)+k
        r=Q(int((c.trpar if upward else c.trch)[j]),c.Sr)
        if f=='f':return r,Q(int(c.f2[j])*d+int(c.f3[j]),d**3)
        others=c.embed(x[1-i]);D=c.D
        def controls(w):
            inf=c.state[('inf',)]; t1=c.state[('T1',k,h)]; u=c.state[w]
            if i==0:
                return int(c.g[inf,u])-int(c.g[t1,u]),int(c.yp[j,c.child[w]])
            return int(c.g[u,inf])-int(c.g[u,t1]),int(c.yw[c.par[w],k])
        a,b=(sum(w*q[z] for t,w in others for q in (controls(t),)) for z in (0,1))
        return r,(3*Q(D,d)**2*(1-Q(D,d))*a/D+Q(D,d)**3*b)

    @lru_cache(maxsize=80000)
    def terms(self,f,x):
        """Uncentered signed row plus exact outside donor coefficients."""
        if not self.supported(f,x):return Q(0),()
        p=SPEC[f]; rows=embedding(p); n=len(p)
        old=self.weight(f,x); val=-self.c.lam*old; free=[]
        for i,(L,k,h) in enumerate(x):
            ds=[j for j in range(n) if p[j]==i]+([p[i]] if i else [])
            dp=self.weight(f,label(x,i,up(x[i])))-old
            val+=(L+sum(Q(1,sum(x[j])) for j in ds))*dp
            if L:
                xx=label(x,i,bud(x[i]))+(NB,)
                extra=sum(self.weight(t,tuple(xx[j] for j in ids)) for t,ids in embedding(p+(i,)) if len(ids)==n+1)
                if f=='f':extra+=self.weight('f',(NB,))
                val+=Q(L,sum(x[i]))*(self.weight(f,xx[:-1])-old+extra)
            for par,num in ((True,h-int(i!=0)),(False,k-p.count(i))):
                if num:
                    center,coeff=self.transferred(f,x,i,par)
                    val+=num*center*coeff
                    free.append((i,par,num,self.sign*(dp-coeff)/self.c.S))
        for g,ids in rows:
            if len(ids)!=n-1:continue
            donor=next(j for j in range(n) if j not in ids)
            i,=tuple(j for j in ids if p[j]==donor or p[donor]==j)
            sub=tuple(x[j] for j in ids)
            if not self.supported(g,sub):continue
            center,coeff=self.transferred(g,sub,ids.index(i),p[i]==donor)
            val+=(Q(1,sum(x[donor]))-center)*coeff
        return self.sign*val/self.c.S,tuple(free)

    @lru_cache(maxsize=80000)
    def rows(self,z,h=1,apex=True):
        p,x=shape(z,h); children=[[] for _ in p]
        for i,j in enumerate(p[1:],1):children[j].append(i)
        out=defaultdict(lambda:[Q(0),Q(0)])
        for f,ids in embedding(p):
            if apex and ids[0]:continue
            v,terms=self.terms(f,tuple(x[j] for j in ids)); slope=Q(0)
            for j,par,num,k in terms:
                i=ids[j]
                if par and p[i]<0:
                    assert num==1
                    slope+=k
                else:
                    ds=[p[i]] if par else [j for j in children[i] if j not in ids]
                    assert len(ds)==num
                    v+=sum(Q(1,sum(x[j])) for j in ds)*k
            assert v+min(Q(0),slope/2)>=0,(f,tuple(x[j] for j in ids),v,slope)
            out[f][0]+=v;out[f][1]+=slope
        return dict(out)
