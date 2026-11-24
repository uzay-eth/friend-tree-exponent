"""Copy-count tilt certificate and exact positive coefficient propagation."""
from argparse import ArgumentParser
from collections import defaultdict
from fractions import Fraction as Q
from math import lcm
from pathlib import Path
import sys

HERE=Path(__file__).resolve()
sys.path.insert(0,str(HERE.parents[2]))
from first_copy import LO, HI
sys.path.insert(0,str(HERE.parents[6]/'finite'/'gauge'))
from physical_trigger_sector import geometry, physical, virgin

A=Q('.565'); C=Q('.988137'); VMIN=Q('.126364920098')
D=Q(3,20); RHO=Q(51,100)
CERTS=((400,Q('1.5'),Q('1.00101'),Q('.994')),
       (250,Q('1.344'),Q('1.00118'),Q('.997')),
       (180,Q('1.27'),Q('1.00132'),Q('.9991')))
CAP2=1/(2*(D+Q(1,20))*VMIN*(1-CERTS[0][3]))


def bounds(low,seam,step):
    # Global (unrestricted geometry) continuations with the direct-only
    # resolvent are used after each seam or step cutoff.
    alpha=Q('.659')
    sig=LO+Q(1,20)-Q('.36');gamma=LO+Q(1,20)
    # Zero-copy bound: any arbitrary nonvirgin state has b v_max as
    # source, and even a virgin satisfies h<= b, as v_max<1.
    v=Q(1,4)/(1-alpha);assert v<1
    c=Q(13,14)/gamma
    # Apply the intrinsic inverse to copy W; then propagate the resulting
    # newborn source through all direct-reset generations.
    M=(Q(13,20)/sig)*(1+5/(gamma*(1-alpha))*Q(13,14))
    vals=[c,c*M,c*M*M]
    return [v+sum(seam[j]*vals[k-j] for j in range(k+1))+step*max(vals[:k+1])
            for k,v in enumerate(low)]


def verify_tilt(n,z0,z,theta):
    s0=RHO+(z0-1)*(Q(1,2)+D)
    trial=A+(s0-A)/n
    cap=1/(2*(D+Q(1,20))*VMIN*(1-theta))
    assert trial<=LO and z**n<z0
    assert C**(n-1)*Q(50,7)<theta**n
    assert LO>=RHO+(z-1)*(Q(1,2)+D)
    p=(Q('.5665')-A)/(Q('.5665')-Q('.462'))
    a,b=p.numerator,p.denominator
    assert Q('3.489')**a*Q('.970148')**(b-a)<C**b
    milestones=[]
    for eps in (Q(1),Q(1,100),Q(1,10**5),Q(1,10**8)):
        k=1
        while cap/z**k>eps:k*=2
        lo,hi=k//2,k
        while lo+1<hi:
            m=(lo+hi)//2
            if cap/z**m<=eps:hi=m
            else:lo=m
        assert cap/z**hi<=eps and cap/z**(hi-1)>eps
        milestones.append(hi)
    print('tilt radius',z,'anchor discount',trial,'unused discount',LO-trial,
          'prefactor',float(cap),'decay',float(1/z),'copy cutoffs',milestones,flush=True)


def transitions(t):
    n,d,_=physical.statistics(t);p=Q(1,d);b=Q(t[1],d)
    base=list(physical.intrinsic_moves(t));copies=[]
    if b and d>2:base.append((physical.newborn(min(d,20)),b))
    for y in t[2]:
        a=physical.mark(y[0]);assert a>=p
        copies.extend(((y,a-p),(physical.plus(y),p)))
    return physical.mark(t[0])+n-p,base,copies


def solve(s,H,K,N):
    keys=[physical.newborn(2)];seen={keys[0]};rows=[]
    rate=K*(H-1)+2;unit=lcm(*range(2,H+1),s.denominator)
    for t in keys:
        d,base,copies=transitions(t)
        base.append((t,rate-d));assert rate>=d
        row=[]
        for arrows in (base,copies):
            inside=defaultdict(int);outside=Q(0)
            for u,r in arrows:
                if not r:continue
                h,k=geometry(u)
                if h>H or k>K:
                    # Independent first-exit accounting in weighted units.
                    outside+=r*weight(u)
                    continue
                if u not in seen:keys.append(u);seen.add(u)
                assert (r*unit).denominator==1
                inside[u]+=int(r*unit)
            row.append((inside,outside))
        rows.append(row)
    idx={t:i for i,t in enumerate(keys)}
    arcs=[[[ (idx[u],r) for u,r in a.items()] for a,b in row] for row in rows]
    E=lcm(*(b.denominator for row in rows for a,b in row))*unit
    exits=[[int(b*E) for a,b in row] for row in rows]
    weights=[weight(t) for t in keys];L=lcm(*(w.denominator for w in weights))
    dwell=int((rate+s)*unit);power=1
    sources=[int(virgin(t))*unit//2 for t in keys]
    v=[[0]*len(keys) for _ in range(4)];v[3]=[int(w*L) for w in weights]
    seams=[[0]*len(keys) for _ in range(3)]
    for _ in range(N):
        y=[[0]*len(keys) for _ in range(4)];e=[[0]*len(keys) for _ in range(3)]
        for i,(a,c) in enumerate(arcs):
            y[0][i]=sum(r*v[0][j] for j,r in a)+sources[i]*power
            y[3][i]=sum(r*v[3][j] for aa in (a,c) for j,r in aa)
            for k in (1,2):y[k][i]=sum(r*v[k][j] for j,r in a)+sum(r*v[k-1][j] for j,r in c)
            for k in (0,1,2):
                e[k][i]=sum(r*seams[k][j] for j,r in a)
                if k==0:e[k][i]+=unit*exits[i][0]*power
                if k==1:e[k][i]+=unit*exits[i][1]*power
                if k:e[k][i]+=sum(r*seams[k-1][j] for j,r in c)
        v,seams=y,e;power*=dwell
    lower=[Q(v[k][0],power) for k in (0,1,2)]
    seam=[Q(seams[k][0],E*power) for k in (0,1,2)]
    future=Q(v[3][0],L*power)
    assert all(x>0 for x in lower+seam) and future>0
    # Full un-tilted residual kernel, including any number of direct resets,
    # has the same W factor on arbitrary initial states as the tilted one.
    print('trial',s,'states',len(keys),'steps',N,'lower',*[float(x) for x in lower],
          'first exit W',*[float(x) for x in seam],'step remainder W',float(future),flush=True)
    return lower,seam,future


def weight(t):
    c=1+len(t[2]);h=(1,Q(17,20),Q(39,50),Q(7,10))[min(c,4)-1]
    return (h+sum(physical.mark(y[0])*weight(y) for y in t[2]))/(physical.mark(t[0])+D)


def fixtures():
    y=(3,3,((2,1,()),));t=(2,6,(y,));a,p=Q(1,3),Q(1,8)
    assert transitions(t)[2]==[(y,a-p),(physical.plus(y),p)]
    assert physical.plus(y)[2]==y[2]
    for u in (t,y,physical.plus(y),physical.newborn(20)):
        assert weight(u)==weight(physical.plus(u))
        assert sum(r*weight(v) for v,r in transitions(u)[2])<=(physical.mark(u[0])+D)*weight(u)
    # Nested surviving geometry: exact marker polynomial through two copies,
    # with intrinsic steps represented separately and no compressed arm reset.
    t=(2,3,((3,3,((3,1,()),)),))
    def walk(u,k):
        if not k:return [Q(1),Q(0),Q(0)]
        den,a,b=transitions(u);ans=[Q(0)]*3
        for ex,arrows in enumerate((a,b)):
            for v,r in arrows:
                if r:
                    for i,z in enumerate(walk(v,k-1)[:3-ex]):ans[i+ex]+=r*z/(LO+den)
        return ans
    c=walk(t,2);den,_,b=transitions(t)
    direct=sum(r*rr/((LO+den)*(LO+transitions(u)[0])) for u,r in b for _,rr in transitions(u)[2]);
    assert c[2]==direct>0
    print('marked shared-subtree and exact two-copy polynomial fixtures passed',flush=True)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__);p.add_argument('--degree',type=int,default=4)
    p.add_argument('--cores',type=int,default=3);p.add_argument('--steps',type=int,default=64)
    a=p.parse_args()
    for row in CERTS:verify_tilt(*row)
    fixtures()
    for s in (LO,HI):
        l,b,r=solve(s,a.degree,a.cores,a.steps)
        # W-resolvent bound; never use F(mu)=1 at a trial.
        print('coefficient upper from seam and step budget',*[float(x) for x in bounds(l,b,r)],flush=True)
        # Import only outward, independently budgeted all-degree results from
        # the previous pilot; neither relies on subtracting from F(mu)=1.
        p0,p01=(Q('.2842876910'),Q('.3388205426')) if s==LO else (Q('.2842536518'),Q('.3387544144'))
        assert p0>l[0] and p01>l[0]+l[1]
        print('combined trial lower for 0/1/2',float(p01+l[2]),flush=True)
        p01u=Q('.4620492561') if s==LO else Q('.4618377419')
        u2=min(bounds(l,b,r)[2],(CAP2-p01)/CERTS[0][2]**2)
        print('combined trial upper for 0/1/2',float(min(CAP2,p01u+u2)),flush=True)
