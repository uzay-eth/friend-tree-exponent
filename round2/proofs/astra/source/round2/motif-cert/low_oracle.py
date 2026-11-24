"""Independent Fraction formulas for the source, knot transfers, rates and newborns."""
from fractions import Fraction as Q
from functools import lru_cache
from pathlib import Path
from random import Random
import argparse
import subprocess
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"motif-high"))
from oracle import Certificate, HEADER, INF, NB, bud, degree, gain, inv


class Local(Certificate):
    def __init__(self,path):
        super().__init__(path)
        hh=HEADER.unpack_from(self.mapping)
        self.coef=hh[7:10]
        self.grid=tuple(Q(n,self.par["QK"]) for n in hh[38:38+self.par["nk"]])
        self.l2=self.levels[2]
        self.ft=self.sections["f"][0]
        self.fty=self.sections["fy"][0]

    @lru_cache(maxsize=32768)
    def f(self,s):
        D=degree(s)
        if D<self.par["DF"]:
            return self.data[self.ft+D*(D-1)//2-1+s[1]]
        q=Q(s[0],D)
        a,b,c=self.coef
        return a+b*q if self.sign==1 else a*q-b*q*q-c*q/D

    @lru_cache(maxsize=32768)
    def knots(self,fam,u):
        n=len(self.grid)
        off=self.fty+(fam*self.l2.inf+u)*n
        return tuple(self.data[off:off+n])

    def fy(self,fam,u,x):
        y=self.knots(fam,u)
        for i,(a,b) in enumerate(zip(self.grid,self.grid[1:])):
            if a<=x<=b: return ((b-x)*y[i]+(x-a)*y[i+1])/(b-a)
        raise ValueError("reciprocal outside the grid")

    @lru_cache(maxsize=32768)
    def corner_state(self,i):
        if i==self.l2.inf: return INF
        if i<self.l2.inf: return self.l2.internal[i]
        k=i-self.l2.inf-1
        return self.l2.d0-k-1,k

    def vertex(self,s,target):
        w=self.f(s)
        D=degree(s)
        diff=self.f(gain(s))-w
        if D<self.l2.d0:
            r=s[0]+self.center(0,s)+s[1]*self.center(1,s)
        else:
            r=s[0]+(Q(s[1]+1,2) if self.sign*diff>=0 else 0)
        born=0 if s[0]==0 else self.f(bud(s))-w+self.f(NB)+self.weight("E",(bud(s),NB))
        return self.sign*(r*diff+Q(s[0],D)*born-target*w)

    def edge(self,labels,target):
        states=tuple(self.corner_state(i) for i in labels)
        W=lambda a,b:self.sign*self.tables["E"][self.offset("E",0,(a,b))]
        total=-target*W(*labels)
        for j,(v,s) in enumerate(zip(labels,states)):
            i=1-j;other=labels[i]; neighbor=states[i]
            cur=W(*labels)
            def entry(end):
                return W(end,other) if j==0 else W(other,end)
            if s is INF: continue
            L,k=s; D=degree(s)
            if v>self.l2.inf:
                jump=entry(self.l2.inf)-cur
                known=inv(neighbor) if other<self.l2.inf else 0
                lower=max(0,Q(self.l2.d0-1-k+known,self.l2.d0+1)) if k<=self.l2.kt else 0
                total+=max(jump,lower*jump)
                if k<=self.l2.kt:
                    d=entry(v+1)-cur
                    total+=max(d,Q(max(0,self.l2.d0-1-k),self.l2.d0)*d)
                continue
            nxt=self.l2.mu(gain(s))[0][0]; delta=entry(nxt)-cur
            total+=(inv(neighbor)-self.center(1-j,s))*self.sign*(self.f(gain(s))-self.f(s))
            total+=L*delta
            if other<self.l2.inf: total+=inv(neighbor)*delta
            elif other>self.l2.inf:
                d0=(W(nxt,self.l2.inf)-W(v,self.l2.inf) if j==0 else
                    W(self.l2.inf,nxt)-W(self.l2.inf,v))
                total+=Q(max(d0,delta),self.l2.d0)
            if L:
                post=bud(s)
                change=entry(self.l2.id[post])-cur
                newborn=(self.weight("C",(post,neighbor,NB)) if j==0
                         else self.weight("V",(neighbor,post,NB)))
                total+=Q(L,D)*(change+self.sign*newborn)
            # FY source on the edge providing a donor, remainder on the edge receiving it.
            if s not in self.act["V"] and k:
                count=k if j else 1
                inc,con=(0,1) if j else (1,0)
                total+=count*self.sign*self.fy(inc,v,inv(neighbor))
                total+=count*max(x*delta-self.sign*y for x,y in zip(self.grid,self.knots(con,v)))
            elif s in self.act["V"]:
                total+=(k*self.center(3,s) if j else self.center(2,s))*delta
            if j==0 and k>=2:
                if s in self.act["C"]: total+=(k-1)*self.center(4,s)*delta
                else:
                    total+=(k-1)*(self.sign*self.fy(2,v,inv(neighbor))+
                           max(x*delta-self.sign*y for x,y in zip(self.grid,self.knots(2,v))))
        return total


def probes(c,exe,path,target,queries):
    inp="".join(f"{name} {a} {b}\n" for name,a,b in queries)
    res=subprocess.run([exe,"probe",path,str(target.numerator),str(target.denominator)],
                       input=inp,text=True,capture_output=True,check=True)
    vals=[int(x) for x in res.stdout.splitlines()]
    if len(vals)!=len(queries): raise ArithmeticError("incomplete probe response")
    return dict(zip(queries,vals))


def audit(path,exe,target,count=1200):
    c=Local(path);rnd=Random(1789)
    rows=set(); l=c.l2; cases=(2,3,c.par["DC"]-1,c.par["DH"]-1,l.d0-1,l.d0,l.d0+1,c.par["DF"]-1)
    for D in cases:
        for k in range(D): rows.add(("v",D-k-1,k))
    for _ in range(count):
        D=rnd.randrange(2,c.par["DF"]);k=rnd.randrange(D)
        rows.add(("v",D-1-k,k))
    representatives={l.inf,l.inf+1,l.n-1}
    for D in (2,3,5,10,14,22,23,24,28,32,38,l.d0-2,l.d0-1):
        if 2<=D<l.d0:
            for k in {0,1,2,min(9,D-1),D-1}:
                if 0<=k<D: representatives.add(l.id[D-1-k,k])
    for p in sorted(representatives & l.parent.keys()):
        for w in representatives: rows.add(("g",p,w))
    parents=list(l.parent)
    for _ in range(count): rows.add(("g",rnd.choice(parents),rnd.randrange(l.n)))
    rows=sorted(rows);lookup=probes(c,exe,path,target,rows)
    largest=Q(0)
    for key in rows:
        tag,a,b=key
        exact=c.vertex((a,b),target) if tag=="v" else c.edge((a,b),target)
        excess=lookup[key]-exact
        if not 0<=excess<32: raise ArithmeticError(f"{key} outward rounding error {excess}")
        largest=max(largest,excess)
    print(f"Fraction low corner oracle: vertex={sum(k[0]=='v' for k in rows)} "
          f"edge={sum(k[0]=='g' for k in rows)} largest rounding excess={largest} ticks PASS")

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("snapshot");p.add_argument("--target",required=True,type=Q)
    p.add_argument("--low",default="round2/outputs/mc-low")
    args=p.parse_args();audit(args.snapshot,args.low,args.target)
