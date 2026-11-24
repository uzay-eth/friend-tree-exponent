"""Full-degree coefficient-filtered copy envelopes and exact first-exit transfer."""
from fractions import Fraction as Q
from pathlib import Path
from functools import cache
from argparse import ArgumentParser
from collections import defaultdict
from math import lcm
import sys

HERE=Path(__file__).resolve()
sys.path.insert(0,str(HERE.parents[2]))
from replay import root, arm_step, J, rnd, LO, HI
sys.path.insert(0,str(HERE.parents[7]/'finite'/'gauge'))
from physical_trigger_sector import geometry, physical, virgin

S=LO
REPS=(2,3,4,5,6,8,10,12,15,20)


class Bound:
    def __init__(self,H,n):
        self.H,self.n=H,n
        self.alpha=root(S,H,20,-1,lambda _:1)[2][0]
        h={i:Q(1,2)/(S+Q(1,i)+Q(3,2)) for i in J}
        self.v0,self.e0=self.reset(h)
        assert max(self.v0.values())==self.v0[20] and self.v0[20]>=Q(1,2)/(S+Q(1,20)+Q(3,2))
        self.R0={j:root(S,H,j,-1,self.v0.__getitem__,Q(1,2)) for j in J}
        self.C0={j:max(self.v0[20]/(S+Q(1,j)),*(v for row in self.R0[j][2:] for v in row)) for j in J}
        assert all(self.C0[j]==self.v0[20]/(S+Q(1,j)) for j in J)
        self.A1={i:{j:self.current(i,j,2,0,0) for j in J} for i in REPS}
        self.B1={i:self.A1[min(r for r in REPS if r>=i)] for i in J}
        b={i:self.table(i,self.B1[i],{} )[2][0] for i in J}
        self.v1,self.e1=self.reset(b)
        self.R1={i:self.table(i,self.B1[i],self.v1) for i in J}
        self.A2={i:{j:self.current(i,j,2,0,1) for j in J} for i in REPS}
        self.B2={i:self.A2[min(r for r in REPS if r>=i)] for i in J}
        b={i:self.table(i,self.B2[i],{})[2][0] for i in J}
        self.v2,self.e2=self.reset(b)
        self.R2={i:self.table(i,self.B2[i],self.v2) for i in J}

    def reset(self,h):
        H=self.H;w=h
        for _ in range(self.n):
            w={i:h[i]+root(S,H,i,-1,w.__getitem__)[2][0] for i in J}
        e=max(Q(0),*(h[i]+root(S,H,i,-1,w.__getitem__)[2][0]-w[i] for i in J))/(1-self.alpha)
        return {i:rnd(w[i]+e,-1) for i in J},e

    def table(self,i,branch,reset):
        lam=S+Q(1,i)
        old=root(S,self.H,i,-1,lambda j:branch[j]+reset.get(j,0))
        # Exactly one virgin-root promotion fee b*branch[2]. Only the virgin
        # state can reach that action, and it is never reached from an old root.
        old[2][0]+=rnd(branch[2]/(2*(lam+Q(3,2))),-1)
        return old

    def R(self,T,i,d,k,fee):
        if d<=self.H:return T[i][d][k]
        lam=S+Q(1,i)
        return fee*(1/lam-Q(k+1,d)/(lam+1))

    def cap0(self,j):return self.C0[j]

    def fee1(self,j,m):return Q(1,m)*self.cap0(m)/(S+Q(1,j))

    @cache
    def current(self,i,j,d,k,n,port=None):
        """First and second copy fee for an already-grown attached child."""
        H=self.H;lam=S+Q(1,i);g=S+Q(1,j);a=Q(1,j)
        if d>=H:
            factor=1/(g*lam)-Q(k+1,d+1)/((g+1)*(lam+1))
            if not n:return rnd(a*self.v0[20]*factor,-1)
            return rnd(a*(self.v1[20]+self.B1[j][20])*factor+
                       a/lam*self.fee1(j,20)*(1/lam-Q(k+1,d)/(lam+1)),-1)
        l=d-k-1;b=Q(l,d);R=self.R1 if n else self.R0
        post=self.current(i,j,d,k+1,n,port) if l else 0
        if n:post+=a/lam*self.fee1(j,min(d,20))
        r,r1=R[j][d][k],R[j][d+1][k]
        return rnd(arm_step(a*r+b*post,lam+b,(l,l+Q(k,2)),(0,min(a,port) if port else a),
                            self.current(i,j,d+1,k,n,port),r1-r,-1,True),-1)

    @cache
    def local(self,t):
        i,L,ch=t;d=1+L+len(ch);k=len(ch);lam=S+Q(1,i)
        p0=self.R(self.R0,i,d,k,self.v0[20])
        p1=self.R(self.R1,i,d,k,self.v1[20]+self.B1[i][20])
        p2=self.R(self.R2,i,d,k,self.v2[20]+self.B2[i][20])
        for y in ch:
            j=y[0];a=Q(1,j)
            c=sum(self.fee1(j,u[0]) for u in y[2])
            # Bound a detached old y's future reset and birth fee, plus
            # retained grandchildren and all later grandchildren of y.
            dy,ky=1+y[1]+len(y[2]),len(y[2])
            x,z=(self.current(i,j,dy,ky,n,Q(1,d)) for n in (0,1))
            p1+=x;p2+=z+a/lam*c
        return tuple(rnd(x,-1) for x in (p0,p1,p2))


def trans(t):
    n,d,_=physical.statistics(t);p=Q(1,d);b=Q(t[1],d)
    a=list(physical.intrinsic_moves(t));c=[]
    if d>2 and b:a.append((physical.newborn(min(d,20)),b))
    for y in t[2]:c.extend(((y,physical.mark(y[0])-p),(physical.plus(y),p)))
    return physical.mark(t[0])+n-p,a,c


def fixtures():
    # One exact marked seam with unequal shifted/unshifted exits. An exit at
    # tally one reads coefficient one, not two, to produce source coefficient two.
    y=(3,1,((2,1,()),));t=(2,5,(y,));p=Q(1,7);a=Q(1,3)
    _,_,c=trans(t)
    assert c==[(y,a-p),(physical.plus(y),p)]
    x0=(Q(1,7),Q(2,9),Q(3,11));x1=(Q(1,5),Q(4,9),Q(7,10))
    coeff=[Q(0)]+[(a-p)*x0[k]+p*x1[k] for k in range(3)]
    z=Q(3,4)
    poly=lambda x:sum(v*z**k for k,v in enumerate(x))
    moment=lambda x:sum(k*k*v*z**k for k,v in enumerate(x))
    m=lambda x:sum(k*v*z**k for k,v in enumerate(x))
    assert poly(coeff)==z*((a-p)*poly(x0)+p*poly(x1))
    assert moment(coeff)==z*sum(r*(moment(x)+2*m(x)+poly(x)) for r,x in ((a-p,x0),(p,x1)))
    assert coeff[2]==(a-p)*x0[1]+p*x1[1]
    print('exact shifted-seam tally and log-z moment fixtures passed',flush=True)


def transfer(s,H,K,N,B):
    ini=physical.newborn(2);keys=[ini];idx={ini:0};edges=[];seams=[]
    rate=K*(H-1)+2;unit=lcm(s.denominator,*range(2,H+1));grid=10**20
    for t in keys:
        d,a,c=trans(t);a.append((t,rate-d));assert rate>=d
        row=[];flow=[[Q(0)]*3 for _ in (0,1)]
        for ex,arrows in enumerate((a,c)):
            inner=defaultdict(int)
            for u,r in arrows:
                if not r:continue
                h,k=geometry(u)
                if h>H or k>K:
                    for n,x in enumerate(B.local(u)):flow[ex][n]+=r*x
                    continue
                if u not in idx:idx[u]=len(keys);keys.append(u)
                assert (r*unit).denominator==1
                inner[idx[u]]+=int(r*unit)
            row.append(tuple(inner.items()))
        edges.append(row)
        seams.append([int(rnd((flow[0][n]+(flow[1][n-1] if n else 0))*unit,-1)*grid) for n in range(3)])
    P=[B.local(t) for t in keys]
    up=[[int(v[k]*grid) for v in P] for k in range(3)]
    low=[[0]*len(keys) for _ in range(3)]
    z=[[0]*len(keys) for _ in range(3)]
    term=[int(virgin(t))*unit//2*grid for t in keys]
    den=int((rate+s)*unit);power=1
    def move(x,i,k,a,c):
        return sum(r*x[k][j] for j,r in a)+(sum(r*x[k-1][j] for j,r in c) if k else 0)
    for _ in range(N):
        L=[[0]*len(keys) for _ in range(3)];U=[[0]*len(keys) for _ in range(3)];Z=[[0]*len(keys) for _ in range(3)]
        for i,(a,c) in enumerate(edges):
            for k in range(3):
                t=term[i]*power if k==0 else 0
                L[k][i]=t+move(low,i,k,a,c)
                U[k][i]=t+seams[i][k]*power+move(up,i,k,a,c)
                Z[k][i]=move(z,i,k,a,c)+seams[i][k]*power
        low,up,z=L,U,Z;power*=den
    l,u,b=([Q(v[k][0],power*grid) for k in range(3)] for v in (low,up,z))
    assert all(0<l[k]<u[k] for k in range(3))
    print('trial',s,'states',len(keys),'seam envelope',*[float(x) for x in b],
          'step envelope',*[float(u[k]-l[k]-b[k]) for k in range(3)],flush=True)
    print('filtered 0,1,2 intervals',[(float(l[k]),float(u[k])) for k in range(3)],flush=True)
    old=Q('.9157613370') if s==LO else Q('.9157808984')
    assert u[2]<old
    print('f2 interval width reduction vs source-only diagnostic',float((old-l[2])/(u[2]-l[2])),flush=True)
    return l,u


if __name__=='__main__':
    p=ArgumentParser(description=__doc__);p.add_argument('--rounds',type=int,default=7)
    p.add_argument('--degree',type=int,default=5);p.add_argument('--steps',type=int,default=96)
    p.add_argument('--both',action='store_true');a=p.parse_args()
    fixtures();B=Bound(21,a.rounds)
    print('newborn envelope',*[float(v[2]) for v in (B.v0,B.v1,B.v2)],'reset remainders',
          *map(float,(B.e0,B.e1,B.e2)),flush=True)
    for s in (LO,HI) if a.both else (LO,):transfer(s,a.degree,3,a.steps,B)
