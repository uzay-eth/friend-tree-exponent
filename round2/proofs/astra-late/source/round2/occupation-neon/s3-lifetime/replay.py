#!/usr/bin/env python3
"""Chronological S3 slots and a real fourth-child donor, on the true D18 host."""
import argparse
from fractions import Fraction as Q
from hashlib import file_digest
from itertools import combinations_with_replacement, product
from math import comb
from pathlib import Path
import sys

OCC=Path(__file__).resolve().parents[1]
sys.path[:0]=[str(OCC/'fringe'),str(OCC.parent/'attacks'/'neon1')]
from charges import Charges, Cert, NB
from extract import PAYLOAD_SHA
from charges import up, bud, label

LO=Q('.565673570322');HI=Q('.565775');GRID=10**23
Z=Q(0)
BARE=(0,1,1)
# Each observed tip is virgin, bare degree two, or has degree >=3. The
# last class is only used as an exterior donor; its S3 payoff is zero.
V,B,A=0,1,2
CAP5_FLOOR={
    'lower':Q(900106797407393717984871399,25000000000000000000000000000000000),
    'upper':Q(454237538259898302104217603,25000000000000000000000000000000000),
}


def floor(q):return Q(q*GRID//1,GRID)


def lo_fmt(q):return f'{float(q):.14g}'


def size_floor(l,k,labels,ids,cap):
    assert len(ids)==k
    return 1+l+sum(sum(labels[i])+labels[i][1] if i<len(labels) else cap+1 for i in ids)


def unpaid(l,k,labels,ids,cap):return k>=4 and size_floor(l,k,labels,ids,cap)>12


def check_support():
    from audit import build, physical_view

    def focal_fringe(physical):
        fringe={1}
        for child,parent in enumerate(physical[2:],2):
            if parent in fringe:fringe.add(child)
        return fringe

    parents=(-1,0,1,1,1,1)
    cases=((0,(1,1,1,1),9),(1,(1,1,1,1),10),
           (1,(3,1,1,1),12),(1,(3,2,1,1),13),(1,(3,2,2,2),15))
    for stock,leaves,expected_size in cases:
        physical=build(parents,(1,stock)+leaves)
        children,states,cores=physical_view(physical)
        assert states[1]==(stock,4,1) and children[1]==[2,3,4,5]
        # Count from the ordinary tree, severing the edge to the parent.
        fringe=focal_fringe(physical)
        actual_size=len(fringe);actual_cores=len(fringe.intersection(cores))
        assert (actual_size,actual_cores)==(expected_size,5)
        tips=tuple(states[i] for i in children[1])
        size=size_floor(stock,4,tips,tuple(range(4)),5)
        assert size==actual_size
        paid_by_size12=actual_size<=12
        paid_by_size20=12<actual_size<=20 and actual_cores in (3,4)
        assert unpaid(stock,4,tips,tuple(range(4)),5)==(actual_size>12)
        assert not (unpaid(stock,4,tips,tuple(range(4)),5) and (paid_by_size12 or paid_by_size20))

    # Include a real core grandchild and a tip beyond the tracking cap.
    physical=build(parents+(2,),(1,1,0,6,1,1,1))
    children,states,_=physical_view(physical)
    labels=tuple(states[i] for i in children[1] if sum(states[i])<=5)
    ids=tuple(range(4)) # three exact labels and one overflow
    assert len(labels)==3 and states[2]==(0,1,1)
    assert size_floor(1,4,labels,ids,5)==15<len(focal_fringe(physical))==16
    assert unpaid(1,4,labels,ids,5)
    print('PASS ordinary-tree support: 9/10/12 rejected, 13/15 accepted; core-grandchild and overflow checked',flush=True)


def charges(c,d,k,j):
    return row(c,d,k,(NB,)*j+(BARE,)*(3-j))


def row(c,d,k,tips):
    xs=((d-k-1,k,1),)+tips
    a,free=c.terms('S3',xs)
    slope=Z
    for i,par,num,v in free:
        if i==0:
            assert num==(1 if par else k-3)
            slope=v
        else:
            assert not par and num==xs[i][1]
            a+=num*min(Z,v/2)
    assert a+min(Z,slope*Q(k-2,2))>=0,(xs,a,slope)
    return a,slope


def controls(c):
    # Independently compute old triple potential differences and the
    # three deleted-tip hC transfers, using nine physical donor fixtures.
    from math import prod
    from functools import cache
    C=c.c;tests=0;strict=0
    @cache
    def val(x):
        if x[0] not in C.s3:return Z
        return sum(int(C.s3tab[C.s3[x[0]],C.tt(*ss)])*prod(ww)
                   for items in product(*(C.embed(t) for t in x[1:])) for ss,ww in [zip(*items)])
    for host in ((0,4,1),(1,4,1),(0,6,1)):
        for tips in ((NB,)*3,(NB,BARE,(2,0,1)),((1,1,1),(0,2,1),(25,0,1))):
            x=(host,)+tips;rec=Q(1,sum(host));old=val(x)
            out_rate=Q(1,7)+Q(host[1]-3,5)
            rho=[sum(Q(1,sum(t)) for t in tips)+out_rate]+[rec+Q(t[1],4) for t in tips]
            direct=-C.lam*old
            for i,(l,k,h) in enumerate(x):
                direct+=(l+rho[i])*(val(label(x,i,up(x[i])))-old)
                if l:direct+=Q(l,sum(x[i]))*(val(label(x,i,bud(x[i])))-old)
            center=Q(int(C.rch[C.state[host]]),C.Sr)
            for j,t in enumerate(tips):
                rest=tips[:j]+tips[j+1:]
                def hc(top):
                    if top not in C.hc:return Z
                    return sum(int(C.hctab[C.hc[top],C.tp(a,b)])*u*v for (a,u),(b,v) in product(*map(C.embed,rest)))
                direct+=(Q(1,sum(t))-center)*(hc(up(host))-hc(host))
            direct*=Q(c.sign,C.S)
            a,free=c.terms('S3',x)
            replay=a+sum(s*(Q(1,7) if par else Q(n,5) if i==0 else Q(n,4)) for i,par,n,s in free)
            assert replay==direct>=0
            strict+=bool(direct and any(s for _,_,_,s in free));tests+=1
    assert strict
    for k in range(3,7):
        assert sum(comb(j,2) for j in range(2,k))==comb(k,3)
    print('PASS',tests,'independent physical S3 transfer/generator rows and chronological multiplicities',flush=True)


def census(c,cap=3,mask='unpaid',minimum=2):
    """All physical triples with tip degree <= cap; arbitrary older arms."""
    labels=tuple((d-k-1,k,1) for d in range(2,cap+1) for k in range(d));n=len(labels)
    host={s for s in c.c.s3 if s[-1]==1};dm=max(map(sum,host));km=max(s[1] for s in host)
    good={(sum(s),s[1]) for s in host if mask=='all' or s[1]>=4}
    trip=tuple(ids for ids in combinations_with_replacement(range(n),3) if sum(labels[ids[-1]])>=minimum)
    tab={(d,k):{ids:row(c,d,k,tuple(labels[i] for i in ids)) for ids in trip} for d,k in good}
    lists={k:tuple(reversed(tuple(combinations_with_replacement(range(n+1),k)))) for k in range(km+1)}
    val={}
    for d in range(dm,1,-1):
        for k in range(min(d-1,km),-1,-1):
            l=d-k-1;b=Q(l,d)
            table=tab.get((d,k))
            for ids in lists[k]:
                count=[ids.count(i) for i in range(n+1)]
                q=sum((Q(count[i],sum(x)) for i,x in enumerate(labels)),Z)
                u=val.get((d+1,k,ids),Z);v=val.get((d,k+1,(0,)+ids),Z)
                a,r=Z,Z
                if table and (mask=='all' or unpaid(l,k,labels,ids,cap)):
                    for ii in combinations_with_replacement((i for i in range(n) if count[i]),3):
                        if ii not in table:continue
                        p,slope=table[ii]
                        m=count[ii[0]]*count[ii[1]]*count[ii[2]]
                        if not m:continue
                        if ii[0]==ii[2]:m=comb(count[ii[0]],3)
                        elif ii[0]==ii[1]:m=comb(count[ii[0]],2)*count[ii[2]]
                        elif ii[1]==ii[2]:m=count[ii[0]]*comb(count[ii[1]],2)
                        if m:
                            a+=m*(p-slope*sum((Q(1,sum(labels[i])) for i in ii),Z));r+=m*slope
                transitions=[];num=b*v;den=HI+b
                for i,x in enumerate(labels):
                    count_i=count[i]
                    if not count_i:continue
                    li,ki,_=x;di=sum(x);pos=ids.index(i)
                    def successor(j):return val[d,k,tuple(sorted(ids[:pos]+(j,)+ids[pos+1:]))]
                    gain=successor(labels.index((li+1,ki,1)) if di<cap else n)
                    gr=count_i*(li+Q(1,d))
                    br=Q(count_i*li,di)
                    num+=gr*gain;den+=gr+br
                    if br:num+=br*successor(labels.index((li-1,ki+1,1)))
                    if ki:transitions.append((Q(count_i*ki,2),gain))
                # Fractional linear extreme for independent tip loss rates:
                # in an optimizer take an upper endpoint iff its successor
                # value is below the result. Check every sorted prefix.
                found=[]
                for rho in (q,q+Q(1,2)+Q(count[-1],cap+1)):
                    g=l+rho;xx=num+a+r*rho+g*u;yy=den+g; candidates=[xx/yy]
                    for rate,t in sorted(transitions,key=lambda it:it[1]):
                        xx+=rate*t;yy+=rate;candidates.append(xx/yy)
                    found.append(min(candidates))
                result=floor(min(found));assert result>=0
                val[d,k,ids]=result
    return LO*val[2,0,()],len(tab)*len(trip)


def broad(c,mask='unpaid',split=True):
    """Choose all triples of active degree-two tips, never sample a tip twice."""
    values={};host={s for s in c.c.s3 if s[-1]==1};dmax=max(map(sum,host));kmax=max(s[1] for s in host)
    rows={(sum(x),x[1],j):charges(c,sum(x),x[1],j) for x in host for j in range(4)}
    for d in range(dmax,1,-1):
        for k in range(min(d-1,kmax),-1,-1):
            l=d-k-1;b=Q(l,d)
            # a virgin becomes bare by promoting; either can age out by
            # a gain. A real bare tip receives from its own core child.
            for v in range(k+1):
                for w in range(k-v+1):
                    g=values.get((d+1,k,v,w),Z);u=values.get((d,k+1,v+1,w),Z)
                    start=values.get((d,k,v-1,w),Z)
                    # Same-level bare→aged decreases w, and virgin→bare
                    # decreases v. Both destinations already exist.
                    bl=values.get((d,k,v,w-1),Z);vb=values.get((d,k,v-1,w+1),Z)
                    vr=Q(v)*(1+Q(1,d)); vb_r=Q(v,2)
                    base=Z;coef=Z
                    ids=(0,)*v+(1,)*w+(2,)*(k-v-w)
                    if (d,k,0) in rows and (mask=='all' or unpaid(l,k,(NB,BARE),ids,2)):
                        for j in range(4):
                            count=comb(v,j)*comb(w,3-j)
                            if count:
                                a,r=rows[d,k,j]
                                aa=a-r*Q(3,2) # outside tip sum = total child sum minus 3/2
                                if split:base+=count*aa;coef+=count*r
                                else:base+=count*(a+min(Z,r*Q(k-2,2)))
                    def ratio(q,external,bare):
                        rate=l+q+external;br=w*(Q(1,d)+bare)
                        return (base+coef*(q+external)+rate*g+b*u+vr*start+vb_r*vb+br*bl)/(HI+rate+b+vr+vb_r+br)
                    vals=[ratio(q,p,r) for q,p,r in product((Q(v+w,2),Q(v+w,2)+Q(k-v-w,3)),
                                                           (Z,Q(1,2)),(Z,Q(1,2)))]
                    assert min(vals)>=0
                    values[d,k,v,w]=floor(min(vals))
    return LO*values[2,0,0,0],rows


def tagged_reward(a,r,k,z,aware):
    # k−2 exterior donors: the mother and k−3 unselected children.
    baseline=a+min(Z,r)*Q(k-2,2)
    return baseline+(max(Z,r)*z+max(Z,-r)*(Q(1,2)-z) if aware else Z)


def check_tagged_regression(c):
    if c.sign!= -1:return
    from audit import build, graph, physical_view

    physical=build((-1,0,1,2,2,2,2,3,4),(1,0,1,0,0,1,2,1,1))
    children,states,cores=physical_view(physical)
    host=2;selected=children[host][:3];donor=children[host][3]
    l,k,_=states[host]
    tips=tuple(states[i] for i in selected)
    assert states[host]==(1,4,1) and tips==(BARE,BARE,NB)
    crown={cores[host]}
    for child,parent in enumerate(physical):
        if parent in crown:crown.add(child)
    assert (len(crown),len(crown.intersection(cores)))==(13,7)
    assert unpaid(l,k,(NB,BARE),(B,B,V,A),2)
    degrees=tuple(map(len,graph(physical)))
    z=Q(1,degrees[cores[donor]]);w=Q(1,degrees[physical[cores[host]]])
    assert (z,w)==(Q(1,3),Q(1,2))
    a,r=row(c,sum(states[host]),k,tips)
    assert r==Q(-2031689,10**11)
    reduced=a+r*(z+w)
    legacy_flat=a+r*Q(k-3,2)
    legacy_aware=legacy_flat+(-r)*(Q(1,2)-z)
    assert legacy_flat-reduced==-r/3>0 and legacy_aware-reduced==-r/2>0
    flat=tagged_reward(a,r,k,z,False);aware=tagged_reward(a,r,k,z,True)
    assert 0<=flat<aware==reduced
    print('PASS 13-vertex upper negative regression: legacy flat/aware exceed reduced row by',
          legacy_flat-reduced,legacy_aware-reduced,flush=True)


def tagged(c,rows,aware=True):
    """One tuple: first three promotions, donor: the fourth physical child."""
    v={};kmax=max(k for _,k,_ in rows);dmax=max(d for d,_,_ in rows)
    for d in range(dmax,1,-1):
        for k in range(min(d-1,kmax),-1,-1):
            l=d-k-1;b=Q(l,d);n=min(k,4)
            for x in reversed(tuple(product((V,B,A),repeat=n))):
                # Fixed position labels are chronological, not sorted by
                # degree. Promotion of the host appends to an empty slot.
                dx=v.get((d+1,k,x),Z);bk=v.get((d,k+1,x+(V,) if n<4 else x),Z)
                active=x[:3];valid=unpaid(l,k,(NB,BARE),x+(V,)*(k-n),2) and all(t!=A for t in active) and (d,k,active.count(V)) in rows
                a,r=rows[d,k,active.count(V)] if valid else (Z,Z)
                # All older-state and unselected reciprocals use admissible
                # corners; the same tagged donor drives hub gain and reward.
                boxes=[(Q(1,2),) if t!=A else (Z,Q(1,3)) for t in x]
                def calc(p,other,*qs):
                    gain=l+p+other+sum(qs);base=gain*dx+b*bk
                    rate=HI+gain+b
                    if valid:
                        tagged_q=qs[3]
                        reward=tagged_reward(a,r,k,tagged_q,aware)
                        assert 0<=reward<=a+r*(tagged_q+p+other)
                        base+=reward
                    for i,t in enumerate(x):
                        if t==A:continue
                        y=x[:i]+(A,)+x[i+1:]
                        if t==V:
                            g0=1+Q(1,d);prom=Q(1,2)
                            base+=prom*v[d,k,x[:i]+(B,)+x[i+1:]];rate+=prom
                        else:g0=Q(1,d)
                        # The bare tip's child is a real untracked donor.
                        # Min over its endpoints for each labelled slot;
                        # root gain does not depend on that donor.
                        base+=g0*v[d,k,y];rate+=g0
                    return base,rate
                candidates=[]
                for p,other,*q in product((Z,Q(1,2)),(Z,Q(k-n,2)),*boxes):
                    base,den=calc(p,other,*q)
                    # All bare loss rates enter linearly; the minimizing
                    # fraction is at a rectangular vertex.
                    vs=[v[d,k,x[:i]+(A,)+x[i+1:]] for i,t in enumerate(x) if t==B]
                    candidates.extend((base+sum(r*y for r,y in zip(extra,vs)))/(den+sum(extra))
                                      for extra in product((Z,Q(1,2)),repeat=len(vs)))
                ans=floor(min(candidates));assert 0<=ans<=min(candidates)
                v[d,k,x]=ans
    return LO*v[2,0,()]


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('payloads',nargs='+',type=Path)
    p.add_argument('--tip-cap',type=int,choices=(2,3,4,5,6),default=5)
    p.add_argument('--large-only',type=int,help='directly price triples with maximum tip degree >= threshold')
    args=p.parse_args();assert __debug__
    check_support()
    for path in args.payloads:
        c=Charges(Cert(path));side='lower' if c.sign==1 else 'upper'
        with path.open('rb') as fp:assert file_digest(fp,'sha256').hexdigest()==PAYLOAD_SHA[side]
        controls(c)
        check_tagged_regression(c)
        all_,rows=broad(c,'all');masked,_=broad(c);frozen,_=broad(c,split=False)
        t0=tagged(c,rows,False);t1=tagged(c,rows,True)
        assert all_>=masked>=frozen and t1>=t0 and masked>=t1
        match,_=census(c,2);assert match==masked
        larger,n=census(c,args.tip_cap)
        assert larger>=masked
        if args.tip_cap==5:
            assert larger==CAP5_FLOOR[side]
            print('PASS',side,'broad cap-five floor unchanged from the independently audited result',flush=True)
        print(side,'true host labels',len(rows)//4,'row labels including tip age',len(rows),flush=True)
        for key,r in (('all-age all owners (overlapping diagnostic)',all_),('unpaid k>=4 and true size>12',masked),('unpaid split corners',frozen),
                      ('ordered triple conservative donor',t0),('ordered triple real fifth vertex',t1),
                      ('fifth donor credit',t1-t0),(f'broader tip cap {args.tip_cap}, {n} true row classes',larger),
                      ('difference of certificate floors (not additive)',larger-masked)):
            print(key,r,lo_fmt(r),flush=True)
        if args.large_only:
            assert 2<args.large_only<=args.tip_cap
            value,count=census(c,args.tip_cap,minimum=args.large_only)
            assert 0<=value<=larger
            print('direct large-only reward',args.large_only,args.tip_cap,'charge classes',count,value,lo_fmt(value),flush=True)


if __name__=='__main__':main()
