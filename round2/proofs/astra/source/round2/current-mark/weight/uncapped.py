#!/usr/bin/env python3
"""Exact endpoint proof, with uncapped or coarsened current reciprocal marks."""
from fractions import Fraction as F
from math import ceil
import argparse

GRID=10**9
COARSE=(2,3,4,5,6,7,8,9,10,11,12,13,14,15,20,25,30,35,45,60,80)


def require(ok,msg):
    if not ok: raise ValueError(msg)


def model(M,d,labels):
    j=2
    m=[]
    pool=set(labels or ())
    for D in range(M+2):
        if D in pool: j=D
        m.append(j if labels else max(2,D))
    z=[1/(1+d*i) for i in m]
    k=[(a-b)/d for a,b in zip(z,z[1:])]
    if labels:
        tz,tk=z[M],0*d
    else:
        # D>=M: z(D)<=z(M), r*k(D)<=D*k(D)<=1/[d*(1+d*(M+1))],
        # and k(D)<=k(M). Bounds have no assumption on the actual arity.
        tz=z[M]+1/(d*(1+d*(M+1)))
        tk=k[M]
    return z,k,tz,tk


def analytic_tail(d):
    # P_i = C/i, C=(d+1/2)/d^2. Its difference pays root-gain decay.
    return (d+(F(1,2) if isinstance(d,F) else .5))/d**2


def propose(M,d,labels=None,smooth=False):
    z,k,tz,tk=model(M,d,labels)
    u=.5+d
    P=z[:]
    # Uniform boost guarantees room for rational rounding; the checker below
    # consumes none of the floating point calculations as evidence.
    for _ in range(100):
        old=P[:]
        K=max(k[j]/old[j] for j in range(2,M))
        K=max(K,tk/old[M])
        if smooth:
            C=analytic_tail(d)
            P[M],P[M+1]=C/M,C/(M+1)
        else:
            P[M]=P[M+1]=tz+K/M
        for i in range(M-1,1,-1):
            c=max(z[i],P[i+1],z[i]+max(0,i-1.5)*(u*k[i]+P[i+1]-z[i])/(u+max(0,i-1.5)))
            for j in range(2,M+1):
                pj=1/j
                diff=old[j]-old[j+1] if j<M else 0.
                drop=k[j] if j<M else tk
                ratio=drop/old[j]
                for r in ((pj,i-2+pj) if j<M else (0,i-2+pj)):
                    c=max(c,(u*(z[i]+r*k[i]+ratio/i)+r*P[i+1])/(u+r+diff/(i*old[j])))
            P[i]=c
        if max(abs(x-y) for x,y in zip(old,P))<1e-12: break
    return P


def profile(P,d,rho,extra=None, exact=False):
    M=len(P)-2
    g=rho-d
    H=M*[None]+[[P[M]/g]*(M+1)]
    # A plus copy changes h only at the first child. For the uncapped gauge
    # the extra gain rate there is at most a+d (also true for a coarsening).
    for i in range(M-1,1,-1):
        H[i]=[0]*(i+1)
        for c in range(i,0,-1):
            l=i-c
            b=(F(l,i) if exact else l/i)
            B=b*(H[i][c+1]+P[i]) if l else 0
            upper=l+c*(F(1,2) if exact else .5)+(d+(F(1,2) if exact else .5) if extra is None else extra)
            v=max((r*H[i+1][c]+B)/(g+r+b) for r in (l,upper))
            v=max(v,F(1,100) if exact else .01)
            H[i][c]=F(ceil(GRID*v),GRID) if exact else v
    return H


def check(M,d,rho,P,labels=None,plant=None,smooth=False):
    z,k,tz,tk=model(M,d,labels)
    A=F(1,2)
    require(rho>d and len(P)==M+2,'table parameters')
    require(all(1>=P[i]>=max(P[i+1],z[i]) for i in range(2,M+1)),'edge monotonicity')
    K=max([k[j]/P[j] for j in range(2,M)]+[tk/P[M]])
    if smooth:
        C=analytic_tail(d)
        require(P[M]==C/M and P[M+1]==C/(M+1),'tail edge mismatch')
        # r<=D, u in [d,d+1/2]. For P=C/D, the root gain factor
        # r*(P(D)-P(D+1)-u*k(D)) is nonnegative if
        # C >= u D(D+1)/[(1+dD)(1+d(D+1))]. Otherwise use r<=D.
        # A universal, simpler sufficient tail is checked below.
        require(M*d>=1 and d+A<1 and C>=1/d+K,'infinite source degree corner')
        require(C >= (d+A)/(d*d),'source tail difference bound')
    else:
        require(P[M]>=tz+K/M and P[M]==P[M+1],'infinite source degree corner')
    slack=None
    corners=0
    def record(v,where):
        nonlocal slack,corners
        require(v>=0,f'edge corner {where}, margin={v}')
        slack=v if slack is None else min(slack,v)
        corners+=1
    for i in range(2,M):
        for a in (F(0),A):
            u=a+d
            for r in (F(0),F(i)-F(3,2)):
                record(u*(P[i]-z[i]-r*k[i])+r*(P[i]-P[i+1]),(i,r,a))
            for j in range(2,M+1):
                ps=(F(1,j),) if j<M else (F(0),F(1,M))
                slip=P[j]-P[j+1] if j<M else F(0)
                drop=k[j] if j<M else tk
                for p in ps:
                    for r in (p,i-2+p):
                        send=F(1,i) if plant!='wrong-sender' else A
                        v=(u*(P[i]-z[i]-r*k[i])+r*(P[i]-P[i+1]))*P[j]
                        v+=send*(P[i]*slip-u*drop)
                        record(v,(i,j,r,a))
    H=profile(P,d,rho,exact=True)
    virgin=H[2][1]
    require(virgin<=1,f'virgin normalization: {float(virgin):.10f}>1')
    H[2][1]=F(1)
    local=0
    for i in range(2,M):
        for c in range(1,i+1):
            l=i-c
            gain=H[i+1][c]-H[i][c]
            b=F(l,i)*(H[i][c+1]-H[i][c]+P[i]) if l else 0
            for r in (l,l+F(c,2)+A+d):
                require(r*gain+b<=(rho-d)*H[i][c],f'local {(i,c,r)}')
                local+=1
    return corners,local,slack,virgin,H


def run(M,d,rho,labels=None,plant=None,smooth=False):
    guess=propose(M,float(d),labels,smooth)
    P=[F(ceil(x*1.00002*GRID),GRID) for x in guess]
    if smooth: P[M],P[M+1]=analytic_tail(d)/M,analytic_tail(d)/(M+1)
    out=check(M,d,rho,P,labels,plant,smooth)
    n,b,s,v,H=out
    print(f'PASS cutoff={M} labels={len(labels) if labels else "uncapped"} rho={rho} d={d} edge={n} local={b}')
    print(f'edge min={s} virgin={v} tail edge={P[M]} h_min={min(x for row in H[2:] for x in row[1:])}')
    return P,H


def exact_inputs(M,d,labels=None):
    P=[F(ceil(x*1.00002*GRID),GRID) for x in propose(M,float(d),labels,not labels)]
    if not labels: P[M],P[M+1]=analytic_tail(d)/M,analytic_tail(d)/(M+1)
    return P


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--M',type=int,default=160)
    p.add_argument('--d',type=F,default=F('.34'))
    p.add_argument('--rho',type=F,default=F('.553'))
    p.add_argument('--labels')
    p.add_argument('--coarse',action='store_true')
    p.add_argument('--constant-tail',action='store_true')
    p.add_argument('--plant',choices=('wrong-sender',))
    args=p.parse_args()
    if args.coarse:
        require(args.labels is None, 'choose --coarse or --labels')
        run(80,F('.34'),F('.56495'),COARSE,args.plant)
        return
    labels=[int(x) for x in args.labels.split(',')] if args.labels else None
    require(labels is None or sorted(set(labels))==labels and labels[0]==2 and labels[-1]<=args.M,'invalid labels')
    run(args.M,args.d,args.rho,labels,args.plant,not labels and not args.constant_tail)


if __name__=='__main__': main()
