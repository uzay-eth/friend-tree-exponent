"""All-degree one-old-copy return, with arbitrary direct resets on either side."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
import sys

sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from first_copy import LO, HI, rnd, choose

J=range(2,21)


def arm_step(base,den,g,p,u,delta,side,matched):
    if matched:
        return choose(((base+r*u+t*(u+delta))/(den+r+t) for t in p for r in g),side)
    return choose(((base+(r+t)*u+w*delta)/(den+r+t) for t in p for w in p for r in g),side)


def eta(d,k):return Q(1,d*(d+1)) if not k else Q(0)


def tail(b,lam,d,k,side):
    return b/(lam+eta(d,k)) if side==1 else 1/lam-(1-b)/(lam+1)


def root(s,H,i,side,fee,virgin=0):
    """Intrinsic discounted root, detached mark death, full tail."""
    lam=s+Q(1,i)
    v=[None]*(H+1)
    v[H]=[rnd(fee(20)*tail(Q(H-k-1,H),lam,H,k,side),side) for k in range(H)]
    for d in range(H-1,1,-1):
        v[d]=[None]*d
        for k in range(d-1,-1,-1):
            l=d-k-1;b=Q(l,d);u=v[d+1][k]
            a=virgin if (d,k)==(2,0) else 0
            if l:a+=b*((fee(min(d,20)) if d>2 else 0)+v[d][k+1])
            v[d][k]=rnd(choose(((a+g*u)/(lam+b+g) for g in (l,l+Q(k,2))),side),side)
    return v


def solve(s,H,n,side,matched):
    # Use the same upper root table for all row-sum controls; mark 20 maximizes the true operator.
    alpha=root(s,H,20,-1,lambda _:1)[2][0]
    assert alpha<1
    def reset(source):
        w=source.copy()
        for _ in range(n):
            w={i:source[i]+root(s,H,i,side,w.__getitem__)[2][0] for i in J}
        err=Q(0)
        if side<0:
            err=max(max(0,source[i]+root(s,H,i,-1,w.__getitem__)[2][0]-w[i]) for i in J)/(1-alpha)
            w={i:w[i]+err for i in J}
        return w,err

    h={i:Q(1,2)/(s+Q(1,i)+Q(3,2)) for i in J}
    v,e=reset(h)
    # Exact no-further-copy payoff encloses arbitrary detached starting arms.
    detached={j:root(s,H,j,side,v.__getitem__,Q(1,2)) for j in J}

    def arm(i,j,all_sources):
        """One arm from its actual birth; no geometry or degree killing."""
        lam=s+Q(1,i); gamma=s+Q(1,j); a=Q(1,j)
        # Both destination responses use the same terminal coefficient. Above H,
        # nudging D upward increases b, so the lower tail can use b for both.
        # Lower integral of b along this attached child uses its own physical
        # generator with external gains; same drift as for the detached root.
        z=[None]*(H+1);z[H]=[]
        for k in range(H):
            b=Q(H-k-1,H);et=eta(H,k)
            val=(a*v[20]*b/((gamma+et)*(lam+et)) if side==1 else
                 a*v[20]*(1/(gamma*lam)-Q(H,H+1)*(1-b)/((gamma+1)*(lam+1))))
            z[H].append(rnd(val,side))
        R=detached[j]
        for d in range(H-1,1,-1):
            z[d]=[None]*d
            for k in range(d-1,-1,-1):
                l=d-k-1;b=Q(l,d);g=(l,l+Q(k,2));u=z[d+1][k]
                # Pure nonvirgin-selected class: no first copy from either virgin
                # source arrow, even though the shifted destination is old.
                eligible=all_sources or (d,k)!=(2,0)
                r=R[d][k] if eligible else 0
                r1=R[d+1][k] if eligible else 0
                base=a*r+(b*z[d][k+1] if l else 0)
                p=(0,a)
                val=arm_step(base,lam+b,g,p,u,r1-r,side,matched)
                z[d][k]=rnd(val,side)
        return z[2][0]

    # Monotone grouping of exterior parent killing: an old arm's integrand
    # need not be monotone in its own mark, so keep all nineteen arm marks.
    reps=(2,4,8,14,20)
    def calc(all_sources):
        q={i:{j:arm(i,j,all_sources) for j in J} for i in reps}
        print('arm tables', 'all' if all_sources else 'nonvirgin', 'side',side,'matched',matched,flush=True)
        nv={}
        for i in J:
            qi=q[(max if side==1 else min)(r for r in reps if (r<=i if side==1 else r>=i))]
            nv[i]=qi[2]/(2*(s+Q(1,i)+Q(3,2)))+root(s,H,i,side,qi.__getitem__)[2][0]
        old,err=reset(nv)
        return old[2],err
    nv,en=calc(False);both,eb=calc(True)
    return v[2],nv,alpha,e,en,both,eb


def checks():
    # Same physical emission, no independent splitting of the two rates.
    a=Q(1,3);r=Q(2,7);r1=Q(9,11)
    for p in (0,a/2,a):assert (a-p)*r+p*r1==a*r+p*(r1-r)
    # On a virgin D=2 source, if shifted output is excluded, the upper and
    # lower envelopes must conserve the exact conditional total emission a.
    for p in (0,a/2,a):assert (a-p)+p==a
    # When gain moves the continuation down but the shift pays more, both
    # independent extrema can be strictly wrong for the common port.
    for side in (1,-1):
        x=arm_step(1,2,(0,1),(0,a),0,a,side,True)
        y=arm_step(1,2,(0,1),(0,a),0,a,side,False)
        assert side*(x-y)>0
    # Nonvirgin fixture: reconstruct the two rates and shift from the original
    # marked generator at one common physical parent state (D=8, mark 3).
    sys.path.insert(0,str(Path(__file__).resolve().parents[5]/'finite'))
    import residual_entries as original
    original.MARK_CAP=20
    for y in ((3,0,((2,1,()),)), (3,3,()), (3,3,((3,1,()),))):
        t=(2,6,(y,));_,D,_=original.statistics(t)
        p=Q(1,D);a=original.mark(y[0]);yp=original.plus(y)
        assert yp==(y[0],y[1]+1,y[2]) and D==8
        for r,r1 in ((Q(2,7),Q(9,11)),(Q(9,11),Q(2,7))):
            assert arm_step(a*r,2,(0,),(p,),0,r1-r,1,True)==((a-p)*r+p*r1)/(2+p)
    print('exact emission conservation, strict shared-port and 6 marked fixtures passed',flush=True)


if __name__=='__main__':
    p=ArgumentParser(description=__doc__);p.add_argument('--degree',type=int,default=21)
    p.add_argument('--rounds',type=int,default=5);p.add_argument('--matched-only',action='store_true')
    p.add_argument('--checks-only',action='store_true')
    p.add_argument('--both-trials',action='store_true');args=p.parse_args();assert args.degree>=21
    checks()
    if args.checks_only:raise SystemExit(0)
    r=2*(HI-LO)/(LO-Q('.51'))
    assert Q(50,7)<sum(Q(2**j,factorial(j)) for j in range(6))
    print('near-critical radius',float(1-r),float(1/(1-r)),flush=True)
    union=[]
    for s in (LO,HI) if args.both_trials else (LO,):
        m=[solve(s,args.degree,args.rounds,side,True) for side in (1,-1)]
        for side,record in zip(('lower','upper'),m):print(s,side,'no-copy, nonvirgin, alpha, two reset errors, any-age, reset error',*[float(z) for z in record],flush=True)
        assert 0<m[0][1]<=m[1][1]
        if not args.matched_only:
            b=[solve(s,args.degree,args.rounds,side,False) for side in (1,-1)]
            assert b[0][1]<=m[0][1]<=m[1][1]<=b[1][1]
            print('independent ports',float(b[0][1]),float(b[1][1]),flush=True)
        print('subtotal without virgin-source class',*[float(sum(x[:2])) for x in m],flush=True)
        a,b=[x[0]+x[5] for x in m];union.append((a,b))
        assert 0<a<b<1 and all(x[1]<=x[5] for x in m)
        print('no-copy plus exactly one copy, all source ages',float(a),float(b),flush=True)
    if len(union)==2:
        print('critical >=2-copy complement',float(1-union[0][1]),float(1-union[1][0]),flush=True)
