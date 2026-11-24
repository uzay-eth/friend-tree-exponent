#!/usr/bin/env python3
"""Physical ceilings and a positive occupation for the exact donor-flux apertures."""
import argparse
from functools import cache
import hashlib
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path
from fractions import Fraction as Q

spec=spec_from_file_location('donor_flux',Path(__file__).resolve().parents[1]/'replay.py')
F=module_from_spec(spec);spec.loader.exec_module(F)
HI=F._older.HI; LO=F._older.LO
GRID=10**24
C_MIN={'lower':Q(113465364604,4*10**11),'upper':Q(115462837236,4*10**11)}
ROWS=(('Yt',((4811,8,1),(0,2,1),(0,1,1),(1,1,1)),2),
      ('P4',((78,2,1),(0,2,1),(0,1,1),(7,0,1)),1),
      ('P4',((1260,9,1),(0,1,1),(0,2,1),(1,1,1)),3),
      ('P4',((18,2,1),(0,1,1),(0,1,1),(1,1,1)),3))


def down(x):return Q(x*GRID//1,GRID)
def ceil(x):return -down(-x)


@cache
def host_ceiling(D,k):
    """Upper stable-age mass for one degree/arity, no exterior law assumption."""
    tab={}
    for d in range(D,1,-1):
        for j in range(min(k,d-1),-1,-1):
            l=d-j-1; b=Q(l,d)
            r=int((d,j)==(D,k)); v=tab.get((d+1,j),0); w=tab.get((d,j+1),0)
            gains=(Q(l),l+Q(j+1,2))
            a=ceil(max((r+b*w+g*v)/(LO+b+g) for g in gains))
            assert all(a*(LO+b+g)>=r+b*w+g*v for g in gains)
            tab[d,j]=a
    return HI*tab[2,0]


def upper_total(side,kind,x,i,slope):
    D=sum(x[0]); k=x[0][1]; n=x[i][1]-F.SPEC[kind].count(i)
    assert n==1
    # Count choices of the middle and (possibly two) tips. At the fixed
    # descriptor their actual child counts bound embeddings pointwise.
    count=k
    if kind=='P4':count*=x[1][1]*x[2][1]
    if kind=='Yt':
        # Exactly two children at the fork; the two designated descriptors
        # are different, so a matching occurrence has one assignment.
        assert x[1][1]==2 and x[2]!=x[3]
    cap=host_ceiling(D,k)*count
    # The outside edge is unique in the chosen slot. This caps either signed
    # excess for every donor life (full send at D=2, or full loss at infinity).
    reward=abs(slope)*cap/2
    print('fixed',kind,'host',D,k,'receiving slot',i,'slope',float(slope),
          'host mass ceiling',float(cap/count),'max embeddings/host',count,
          'reciprocal or loss mass ceiling',float(cap/2),'credit ceiling at supplied denominator',float(reward/F.DEN[side]),
          'true normalized ceiling using GMP c0',float(reward/C_MIN[side]),flush=True)
    return cap,reward


def script(E,side):
    """Positive exact support: lower P4 at (21,k=2), terminal donor degree 3."""
    kind,x,i=ROWS[3]
    # Birth order: apex, middle, middle, terminal, extra apex child, donor.
    # No restriction on the eventual exterior parent's degree. Track chosen
    # physical births and prevent only projected jumps during this one script.
    parents=(-1,0,1,2,0,3)
    steps=[(0,'p'),(1,'p'),(2,'p'),(0,'g'),(0,'p'),(3,'g'),(3,'p')]
    leaves=[1]; lower=Q(1)
    for u,job in steps:
        k=[parents[1:len(leaves)].count(j) for j in range(len(leaves))]
        ds=[l+c+1 for l,c in zip(leaves,k)]; ext=Q(1,2)
        # exact total projected rates when r_parent=0; donors in projection
        # include every child created in the script so only one live port.
        gain=[Q(l)+sum((Q(1,ds[v]) for v in range(len(ds))
                       if parents[v]==j or parents[j]==v),Q(0)) for j,l in enumerate(leaves)]
        bud=[Q(l,d) for l,d in zip(leaves,ds)]
        rate=gain[u] if job=='g' else bud[u]
        assert rate>0
        lower=lower*rate/(HI+sum(gain)+sum(bud)+ext)
        if job=='p':
            assert parents[len(leaves)]==u
            leaves[u]-=1;leaves.append(1)
        else:leaves[u]+=1
    k=[parents[1:].count(j) for j in range(len(leaves))]
    no=tuple((l,c,1) for l,c in zip(leaves,k));assert no[1:4]==x[1:] and no[5]==F.NB
    t=(x[0],)+no[1:]; assert sum(t[0])==21
    # Track the same donor through degree 2,3,4; gains remain live while all
    # other projected jumps kill. A donor promotion is killed, not a reset.
    def episode(source,hi):
        discount=LO if hi else HI
        tab={}
        for host in range(21,sum(no[0])-1,-1):
            for d in (4,3,2):
                ds=[sum(y) for y in t];ds[5]=d;ds[0]=host
                ls=[y[0] for y in t];ls[5]=d-1;ls[0]=host-3
                gain=[l+sum((Q(1,ds[v]) for v in range(len(ds)) if parents[v]==j or parents[j]==v),Q(0))
                      for j,l in enumerate(ls)]
                b=sum(Q(l,q) for l,q in zip(ls,ds))
                base=discount+sum(gain)+b
                adv=tab.get((host+1,d),0)
                fixed=(source(d) if host==21 else 0)+gain[5]*tab.get((host,d+1),0)+gain[0]*adv
                ends=[(fixed+r*adv)/(base+r) for r in (0,Q(1,2))]
                tab[host,d]=(max if hi else min)(ends)
        return tab[sum(no[0]),2],tab[21,2]

    send=lambda d:Q(1,d)
    loss=lambda d:Q(1,2)-Q(1,d)
    sm=LO*lower*episode(send,False)[0]; lm=LO*lower*episode(loss,False)[0]
    assert sm>lm>0
    slope=next(c for j,parent,n,c in F.Flux(E).slopes(kind,x) if j==i and not parent)
    assert (slope>0)==(side=='lower')
    v=slope*sm if slope>0 else -slope*lm
    print('script exact support send lower',float(sm),'loss lower',float(lm),
          'signed excess lower (outward 1e-30 ticks)',v*10**30//1,float(v),flush=True)
    for name,r in [('send',send),('loss',loss)]:
        L,U=episode(r,False)[1],episode(r,True)[1]
        assert 0<L<=U
        print('conditional isolated donor episode',name,'integral bracket',float(L),float(U),flush=True)


def memory_control():
    # At any given degree an actual gain moves d to d+1, never resets.
    # Degree-monotonicity yields this lower rate uniformly over other children.
    for d in (2,3,4):
        assert 1/Q(d)-1/Q(d+1)==Q(1,d*(d+1))
    don=(2,3,4); start=sum(Q(1,2)-Q(1,d) for d in don)
    append=(*don,2)
    assert sum(Q(1,2)-Q(1,d) for d in append)==start
    print('exact add-zero/no-reset memory controls 4',flush=True)


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('blob',type=Path);a=p.parse_args()
    with a.blob.open('rb') as f:digest=hashlib.file_digest(f,'sha256').hexdigest()
    C=F.Cert(a.blob); side='lower' if C.mode==2 else 'upper';assert F.PAYLOAD_SHA[side]==digest
    print(side,'pinned',digest,flush=True);E=F.Local(C); FF=F.Flux(E);memory_control()
    for kind,x,i in ROWS:
        slope=next(c for j,p,n,c in FF.slopes(kind,x) if j==i and not p)
        upper_total(side,kind,x,i,slope)
    script(E,side)


if __name__=='__main__':main()
