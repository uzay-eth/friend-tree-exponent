#!/usr/bin/env python3
"""All-time host-degree payment for live selected young arms."""
from functools import cache
from fractions import Fraction as Q
from math import comb
from charges import Cert, Charges

LO,HI=Q('.565675597658'),Q('.565772677869')

GRID=10**22
ZERO=Q(0)


def floor(x):return Q(x*GRID//1,GRID)


class Young:
    def __init__(self,c,degree,arity,virgins,gate=12,band='high',degree3=False,
                 host_floor=0, family='all'):
        self.c=c;self.degree=degree;self.arity=arity;self.virgins=virgins
        self.gate=gate;self.band=band;self.degree3=degree3
        self.host_floor=host_floor; self.family=family

    @cache
    def affine_toll(self,d,k,m,n):
        """Nonnegative selected young edges, cherries and stars at d>=18.

        Tracked children have k=0 and degree 2 or 3. The unknown arms
        enter only the common apex gain/charge. No crown cap is used.
        """
        if m+n==0 or d<self.host_floor or (self.band=='high' and d<18) or (self.band=='low' and d>=18):return ZERO,ZERO
        host=(d-k-1,k,1); tip=(1,0,1); old=(2,0,1)
        assert 0<m+n<=k and min(host)>=0
        a=b=ZERO
        for name,need in (('g',1),('hC',2),('S3',3)):
            if self.family!='all' and name!=self.family:continue
            for j in range(max(0,need-n),min(m,need)+1):
                size=comb(m,j)*comb(n,need-j)
                q,free=self.c.terms(name,(host,)+(tip,)*j+(old,)*(need-j))
                a+=size*q
                for i,par,num,dp in free:
                    assert i==0
                    if par:
                        assert num==1
                        b+=size*dp
                    else:
                        assert num==k-need
                        a+=size*(Q(m-j,2)+Q(n-(need-j),3))*dp
                        b+=size*dp
        return a,b

    def solve(self):
        # v(d,k,m) with k capped and m the exact number of selected tips.
        # Newborns may overshoot m cap: track no more than M (all are
        # young at their actual births), replacing no existing tip.
        cap,K,M=self.degree,self.arity,self.virgins
        nxt={}
        for d in range(cap,1,-1):
            cur={}
            for k in range(min(d-1,K),-1,-1):
                L=d-k-1; p=Q(1,d); b=L*p
                for m in range(min(k,M)+1):
                    for n in range((min(k,M)-m if self.degree3 else 0)+1):
                        # The crown has at least d+m+2*n ordinary vertices.
                        r,s=self.affine_toll(d,k,m,n) if d+m+2*n>self.gate else (ZERO,ZERO)
                        y=nxt.get((k,m,n),ZERO)
                        bk=cur.get((k+1,m+int(m+n<M),n),ZERO)
                        bad=m*Q(1,2); grow=m*(1+p); leave=n*(2+Q(2,3)+p)
                        c=b*bk
                        if m:c+=bad*cur[k,m-1,n]+grow*cur[k,m-1,n+int(self.degree3)]
                        if n:c+=leave*cur[k,m,n-1]
                        g=L+Q(m,2)+Q(n,3); top=Q(k-m-n+1,2)
                        corners=[ZERO,top]
                        if s and 0<-r/s<top:corners.append(-r/s)
                        ans=min((max(ZERO,r+s*u)+c+(g+u)*y)/(HI+b+bad+grow+leave+g+u) for u in corners)
                        cur[k,m,n]=floor(ans)
            nxt=cur
        return LO*nxt[0,0,0]


if __name__=='__main__':
    import argparse
    from pathlib import Path
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('blobs',type=Path,nargs=2)
    p.add_argument('--degree',type=int,default=65)
    p.add_argument('--arity',type=int,default=8)
    p.add_argument('--virgins',type=int,default=4)
    p.add_argument('--band',choices=('low','high'),default='high')
    p.add_argument('--degree3',action='store_true')
    p.add_argument('--host-floor',type=int,default=0)
    p.add_argument('--family',choices=('all','g','hC','S3'),default='all')
    args=p.parse_args()
    for b in args.blobs:
        c=Charges(Cert(b)); y=Young(c,args.degree,args.arity,args.virgins,band=args.band,degree3=args.degree3,
                                  host_floor=args.host_floor,family=args.family)
        print(c.sign,'all-crown',args.band,'host/young-arm lower reward',float(y.solve()),flush=True)
