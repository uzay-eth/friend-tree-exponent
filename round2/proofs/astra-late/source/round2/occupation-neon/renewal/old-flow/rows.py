"""Single charged D18 higher motif, retaining all unobserved adjacent stubs."""
from functools import cache
from itertools import product
from math import prod
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[2]/'higher'))
from extract import Extract, Cert, Q, SPEC, NB, embedding, label, up, bud


class Rows(Extract):
    @cache
    def value(self,f,x):
        if all(sum(y)<self.c.D for y in x):return self.weight(f,x)
        return sum(self.weight(f,tuple(v for v,w in a))*prod(w for v,w in a)
                   for a in product(*(self.c.embed(y) for y in x)))

    def live(self,f,x):
        if f=='hV':return x[1][2]==1 and 1<=x[1][1]<=self.c.K
        if f=='hC':return 2<=x[0][1]<=self.c.K
        return self.supported(f,x)

    def center(self,y,par):
        c=self.c
        if sum(y)<c.D:
            return Q(int((c.rpar if par else c.rch)[c.state[y]]),c.Sr)
        k,h=y[1:];i=(1-h)*(c.K+2)+min(k,c.K+1)
        return Q(int((c.trpar if par else c.trch)[i]),c.Sr)

    def moved(self,f,x,i,par):
        k,h=x[i][1:]
        if f=='g' and sum(x[i])>=self.c.D:
            return k<=self.c.K and (i==0 and k>=(1 if par else 2) or i==1 and k>=1)
        return self.route(f,x,i,par)

    def transfer(self,f,x,i,par):
        c=self.c;y=x[i];d=sum(y)
        if f!='g' or d<c.D:return self.value(f,label(x,i,up(y)))-self.value(f,x)
        t=Q(c.D,d);k,h=y[1:];a=(1-h)*(c.K+2)+k
        amp=Q(0)
        for v,w in c.embed(x[1-i]):
            z=[v,v];z[i]=('inf',)
            f0=int(c.g[c.state[z[0]],c.state[z[1]]]);z[i]=('T1',k,h)
            f1=int(c.g[c.state[z[0]],c.state[z[1]]])
            ref=int(c.yp[a,c.child[v]] if i==0 else c.yw[c.par[v],k])
            amp+=w*(3*t*t*(1-t)*Q(f0-f1,c.D)+t**3*ref)
        return amp

    @cache
    def affine(self, f, x):
        """Constant and signed unresolved parent/child receive coefficients."""
        c = self.c
        p = SPEC[f]
        n = len(p)
        if not self.live(f,x):
            return Q(0), (Q(0),)*n, Q(0)
        own = self.value(f,x)
        a = -c.lam*own
        child = [Q(0)]*n
        external = Q(0)
        for i,z in enumerate(x):
            L,k,h = z
            neigh = [j for j in range(n) if p[j]==i] + ([p[i]] if i else [])
            delta = self.value(f,label(x,i,up(z)))-own
            a += (L+sum(Q(1,sum(x[j])) for j in neigh))*delta
            if L:
                born = label(x,i,bud(z))+(NB,)
                new = sum(self.value(g,tuple(born[j] for j in ids))
                          for g,ids in embedding(p+(i,)) if len(ids)==n+1)
                a += Q(L,sum(z))*(self.value(f,born[:-1])-own+new)
            for parent,count in ((True,h-int(i!=0)),(False,k-p.count(i))):
                assert count>=0
                if count and self.moved(f,x,i,parent):
                    a += count*self.center(z,parent)*delta
                elif count and parent:
                    external = delta
                elif count:
                    child[i] = delta
        for donor in range(n):
            for g,keep in embedding(p):
                if len(keep)!=n-1 or donor in keep:
                    continue
                j, = [u for u in keep if p[u]==donor or p[donor]==u]
                small = tuple(x[u] for u in keep)
                i = keep.index(j)
                par = p[j]==donor
                if self.live(g,small):
                    assert self.moved(g,small,i,par)
                    a += (Q(1,sum(x[donor]))-self.center(x[j],par))*self.transfer(g,small,i,par)
        factor = Q(self.sign,self.c.S)
        return a*factor, tuple(z*factor for z in child), external*factor

    def controls(self, f, x):
        """Bound the total unknown gain at each receiver, shared with drift."""
        a, coeff, exterior = self.affine(f,x)
        p = SPEC[f]
        r = []
        for i,z in enumerate(x):
            k = z[1] - p.count(i)
            r.append([(Q(0),Q(0)),(Q(k,2),Q(k,2)*coeff[i])] if k else [(Q(0),Q(0))])
        return a,r,exterior

    def floor(self,f,x):
        a,cols,ext = self.controls(f,x)
        return max(0,a+sum(min(z for _,z in v) for v in cols)+min(0,ext/2))
