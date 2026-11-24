#!/usr/bin/env python3
"""Exact apex-owned P4/Yt correction, with a live unbounded donor degree."""
from fractions import Fraction as Q
from functools import cache
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from pathlib import Path
import argparse
import gzip
import hashlib
import json
import re
import sys

HERE = Path(__file__).resolve().parent


def load(name,path):
    spec = spec_from_file_location(name,path)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


O = load('older',HERE.parents[1]/'higher'/'older'/'replay.py')
I = load('corrector',HERE.parent/'corrector.py')
sys.path.insert(0,str(HERE.parents[1]/'check'))
from crowns import n, degree, moves

J = (1,())
BARE = (0,1,1)
BRANCH = (0,2,1)
ZERO = Q(0)
GRID = 10**22
STOCK = HERE.parents[1]/'rows'


def all_floors():
    """Authenticate the integer operands, recompute, and compare all ten rows."""
    raw=(STOCK/'pinned-fg-slabs.json.gz').read_bytes()
    assert hashlib.sha256(raw).hexdigest()=='63aab8b10d07b43abd80e5852e60b72f5513bc23618d83ce6b0fd3ba3d5ab321'
    log=(STOCK/'bands.log').read_bytes()
    assert hashlib.sha256(log).hexdigest()=='61aa4d4b927f5a1416fb744651455fe0529b2b05df20f857fb63204e12f8bc48'
    vals=json.loads(gzip.decompress(raw))
    R=load('stock_check',STOCK/'check.py')
    checked={int(mode):{} for mode in (1,2)}
    for mode,k,a,b in re.findall(rb'^band ([12]) ([0-4]) finite minimum \(Fraction\((\d+), (\d+)\)',log,re.M):
        checked[int(mode)][int(k)]=Q(int(a),int(b))
    result={}
    for mode,name in ((1,'upper'),(2,'lower')):
        key=vals[name]['f']
        assert (key['sign'],key['LF'],key['S'],key['Sr'],key['D0'])==((1 if mode==2 else -1),14000,10**11,10**9,18)
        floors=tuple(R.envelope(row) for _,row in sorted(R.table(key).items()))
        assert len(floors)==5 and tuple(checked[mode][k] for k in range(5))==floors
        result[mode]=floors
    print('authenticated exact all-stock floors: 5/5 each sign',flush=True)
    return result


def down(x):
    return Q(x*GRID//1,GRID)


def band(E,kind,x,receiver,has_donor):
    """Exact affine row in shared apex send and, when present, 1/d_donor."""
    ext = [(ZERO,ZERO)]*4
    a = E.charge(kind,x,ext)
    ext[0] = (Q(1),ZERO)
    p = E.charge(kind,x,ext)-a
    ext[0] = (ZERO,ZERO)
    if has_donor:
        ext[receiver] = (ZERO,Q(1))
        w = E.charge(kind,x,ext)-a
    else:
        w = ZERO
    assert min(a+i*p+j*w for i in (0,Q(1,2)) for j in (0,Q(1,2)))>=0
    return a,p,w


def profiles(E,s):
    # Prefix graphs are full crowns on these faces. Every trigger is a
    # protected-slot promotion, hence the only transition making a copy.
    p4 = (BARE,BARE,BARE,BARE)
    yt = (BARE,BRANCH,BARE,O.NB)
    out = {}
    for name,x,receiver,fixed in (('P4',p4,3,Q(3)),('Yt',yt,2,Q(4))):
        # Sum the six directed inside sends plus one leaf -> virgin gain
        # on the Yt face. Its virgin promotion (1/2) is also a kill.
        a,p,w = band(E,name,x,receiver,True)
        # Allow all donor growth forever. For z=A+B/d, use
        # A_track(1/d)=−g/[d(d+1)] and 0<=g<=d.
        # The numerator is affine in the common apex input r. Allow
        # arbitrary k and descendant geometry at the live donor.
        # Constant + linear matching with conservative 1/d^2 cross-term.
        inequalities = [(-1,0,ZERO),(-1,Q(-1,2),ZERO)] # z>=0
        for r in (ZERO,Q(1,2)):
            for e in (0,1):
                # g=d*e at the two endpoints; exact gain loss is e*u/(1+u).
                # Multiply residual by 1+u and use degree-3 Bernstein.
                H=s+fixed+r
                for j in range(4):
                    u1=Q(j,6)
                    u2=Q(j*(j-1),24)
                    u3=Q(j*(j-1)*(j-2),48)
                    inequalities.append((H+(H+1)*u1+u2,
                                         (H+e)*u1+(H+1)*u2+u3,
                                         (a+p*r)*(1+u1)+w*(u1+u2)))
        vals=[]
        for (ai,bi,ci),(aj,bj,cj) in combinations(inequalities,2):
            det=ai*bj-aj*bi
            if not det: continue
            A,B=(ci*bj-cj*bi)/det,(ai*cj-aj*ci)/det
            if all(i*A+j*B<=v for i,j,v in inequalities):
                vals.append((A+B/2,A,B))
        v,A,B=max(vals)
        assert all(i*A+j*B<=q for i,j,q in inequalities)
        out[name] = (v,a,p,w,A,B)
    return out


def solve(E,s,upper=False):
    pv = profiles(E,s)
    # P4 needs four tip promotions to hit the size-six face. Yt needs
    # a second promotion at the middle and a promotion at its first tip.
    # On the all-degree continuation the donor's entire subtree can grow.
    p4c = J
    for _ in range(4):
        p4c = (0,(p4c,))
    ytc = (0,((0,tuple(sorted((J,(0,(J,)))))),))
    assert n(p4c)==6 and n(ytc)==7

    @cache
    def V(z):
        size = n(z)
        for f,face in (('P4',p4c),('Yt',ytc)):
            if z==face:
                if upper:
                    _,a,p,w,_,_=pv[f]
                    return (a+max(0,p/2)+max(0,w/2))/(s+(3 if f=='P4' else 4))
                return pv[f][0]
        if size>=7:
            return ZERO
        intrinsic = moves(z,ZERO)
        end = sum(q*V(t) for t,q in intrinsic)
        base = s+size-Q(1,degree(z))
        # An exterior gain changes only the root leaf stock, and may
        # enter or leave the active face. Keep the same exterior rate.
        leaf,children = z
        plus = V((leaf+1,children))
        rates=(ZERO,Q(1,2))
        if upper:
            v=max((end+r*plus)/(base+r) for r in rates)
            assert all((base+r)*v>=end+r*plus for r in rates)
        else:
            v=down(min((end+r*plus)/(base+r) for r in rates))
            assert all((base+r)*v<=end+r*plus for r in rates)
        return v

    return V(J),pv,V


def verify_lift(E,s,trial):
    """Sum corrected whole-tree potential changes from individual arrows."""
    v0,pv,V = trial
    from random import Random

    def score(z):
        l,kids = z
        if l==0 and len(kids)==1:
            l1,b = kids[0]
            if l1==0 and len(b)==1:
                l2,c = b[0]
                if l2==0 and len(c)==1 and c[0][0]==0 and len(c[0][1])==1:
                    A,B=pv['P4'][-2:]
                    return A+B/degree(c[0][1][0])
            if l1==0 and len(b)==2 and J in b:
                other=b[1-b.index(J)]
                if other[0]==0 and len(other[1])==1:
                    A,B=pv['Yt'][-2:]
                    return A+B/degree(other[1][0])
        return V(z)

    def W(p,x):
        kids=[[] for _ in x]
        for j,i in enumerate(p[1:],1):kids[i].append(j)
        @cache
        def z(i):return x[i][0],tuple(sorted(z(j) for j in kids[i]))
        return sum(score(z(i))-v0 for i in range(1,len(x)))

    rng=Random(194)
    cases=[O.sample(i,rng) for i in range(3,40)]
    for donor in (1,2,19,130):
        # Put both faces under a founder, and vary the actual donor stock.
        cases.extend(O.fixture(p+[0,0]+[5]*donor) for p in
                     ([-1,0,1,2,3,4],[-1,0,1,2,2,3,4]))
    hit=set()
    for p,x in cases:
        paid,_=E.tree(p,x)
        old=W(p,x)
        gain=[Q(t[0])+sum(Q(1,sum(x[j])) for j in range(len(x)) if p[j]==i or p[i]==j)
              for i,t in enumerate(x)]
        drift=-s*old
        for i,(l,k,h) in enumerate(x):
            drift+=gain[i]*(W(p,O.label(x,i,O.up(x[i])))-old)
            if l:
                drift+=Q(l,l+k+h)*(W(p+(i,),O.label(x,i,O.bud(x[i]))+(O.NB,))-old)
        # Only the protected top-owned P4/Yt faces are booked. On this face
        # the states already imply the stated external-child counts.
        reward=ZERO
        for f,ids,value in paid:
            xx=tuple(x[i] for i in ids)
            if ids[0]!=0 and ((f=='P4' and xx==(BARE,)*4) or (f=='Yt' and xx==(BARE,BRANCH,BARE,O.NB))):
                reward+=value
                hit.add(f)
        assert reward+drift>=s*v0*(len(x)-1),(p,x,reward,drift)
    assert hit=={'P4','Yt'}
    return len(cases)


def historic_control(c, ch, rows, floors):
    """Rebuild only the withdrawn coefficient; build drift by whole scores."""
    s=Q('.565770220167')
    wrong=list(floors[1]); wrong[1]*=1000
    original=I.FLOOR
    try:
        I.FLOOR={**floors,1:wrong}
        _,v0,M,_,tail=I.solve(c,ch,rows,64,s)
    finally:
        I.FLOOR=original
    d=(c.lam-s)*M/s
    for L in (14001,10**6):
        for m in (7,8,12):
            kids=[[i+1] for i in range(m-1)]+[[]]
            x=tuple((L,len(kids[i]),int(i!=0)) for i in range(m))

            def potential(state, children):
                assert min(map(sum,state))>64
                return c.potential(children,state)/c.S-sum(tail[min(y[1],c.KF+1)]-v0 for y in state[1:])+d

            # All old vertices stay in the tail, and the new virgin has
            # exactly zero correction. Do not consult the claimed margin.
            old=potential(x,kids)
            drift=s*old
            for i,(l,k,h) in enumerate(x):
                rate=l+sum(Q(1,sum(x[j])) for j in (i-1,i+1) if 0<=j<m)
                drift-=rate*(potential(O.label(x,i,O.up(x[i])),kids)-old)
                newkids=[row[:] for row in kids]+[[]]; newkids[i].append(m)
                child=O.label(x,i,O.bud(x[i]))+(O.NB,)
                promoted=c.potential(newkids,child)/c.S-sum(tail[min(y[1],c.KF+1)]-v0 for y in child[1:-1])+d
                drift-=Q(l,sum(x[i]))*(promoted-old)
            if drift<0:
                print('withdrawn upper whole-tree counterexample: cores',m,'leaves per core',L,
                      'signed drift',float(drift),flush=True)
                return
    raise AssertionError('no historic counterexample')


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('payload',type=Path,nargs=2)
    a = p.parse_args()
    floors=all_floors()
    assert {k:tuple(map(Q,v)) for k,v in I.FLOOR.items()}==floors
    for path,mode in zip(a.payload,(2,1)):
        c=I.D.Cert(path)
        assert (c.mode,c.D,c.K,c.S)==(mode,18,17,10**11)
        E=O.Local(c)
        s=c.lam+(1 if mode==2 else -1)*Q('0.00001')
        trial=solve(E,s)
        val,prof=trial[:2]
        print('whole-tree lift controls',verify_lift(E,s,trial),flush=True)
        print('mode',mode,'target',float(s),'P4/Yt base credit',float(s*val),flush=True)
        for f,(v,A,p,w,z0,z1) in prof.items():
            print(f,'at virgin donor',float(v),'coefficients',*(float(x) for x in (A,p,w)),flush=True)
        assert val>0
        ch=I.D.Charge(c)
        rows=I.operands(c,64)
        if mode==1:
            historic_control(c,ch,rows,floors)
        I.FLOOR=floors
        margin,_,M,norm,_=I.solve(c,ch,rows,64,s)
        ceiling=s*solve(E,s,upper=True)[0]
        assert margin+s*val<=margin+ceiling<0
        print('first-body margin',float(margin),'combined margin',float(margin+s*val),
              'selected-entry credit ceiling',float(ceiling),flush=True)
        lo,hi=0,10_000_001
        direction=1 if mode==2 else -1
        while hi-lo>1:
            mid=(hi+lo)//2
            ss=c.lam+direction*Q(mid,10**12)
            m=I.solve(c,ch,rows,64,ss)[0]
            if m>=0:lo=mid
            else:hi=mid
        base=c.lam+direction*Q(lo,10**12)
        s=base+direction*(Q(0) if mode==2 else Q('8e-8'))
        v,pv,_=solve(E,s)
        margin,_,M,norm,_=I.solve(c,ch,rows,64,s)
        assert margin+s*v>0 and (mode==2 or margin<0)
        assert norm+2*max(max(A,A+B/2) for *_,A,B in pv.values())<Q('.01')
        print('corrected pointwise base',float(base),'lift endpoint',float(s),
              'base-only margin at lift',float(margin),'margin with lift',float(margin+s*v),flush=True)


if __name__=='__main__':
    main()
