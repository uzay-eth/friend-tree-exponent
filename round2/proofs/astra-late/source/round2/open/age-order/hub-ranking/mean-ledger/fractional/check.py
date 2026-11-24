"""Interval stock-sensitive moment upper and signed ledger."""
from fractions import Fraction as Q
from functools import cache
from itertools import combinations_with_replacement, product
from math import factorial, lcm
from pathlib import Path
from runpy import run_path
import argparse

ROOT=Path(__file__).resolve().parents[2]
FM=run_path(str(ROOT/'fractional-middle/check.py'))
up,down=FM['up'],FM['down']


class Bounds:
    def __init__(self,lo,hi,stop=64):
        self.lo,self.hi,self.stop=lo,hi,stop
        self.m=FM['middle'](lo,hi)
        self.bonuses={}
        self.rows={}
        for k in range(stop,1,-1):
            a=1+Q(1,k)
            gm=FM['gamma_upper']
            cap=(1-1/(hi+a+1))*gm(a,lo,hi)+(gm(a-Q(1,2),lo,hi)+gm(1+Q(1,k+1),lo,hi))/(2*(lo+a+1))
            cap=min(1+Q(1,k),up(cap))
            for j in range(k-1):
                if k==stop:
                    self.rows[k,j]=self.tail(k,j)
                    continue
                gain=k*(j+1)
                numerator=cap+gain*self.U(k+1,j+1)+j*self.rows.get((k,j-1),0)
                bound=up(numerator/(1+gain+j+lo*k))
                self.rows[k,j]=max(self.U(k+1,j+1),bound)
        for k in range(2,stop):
            for j in range(k-1):
                assert self.U(k,j)>=self.U(k+1,j+1),(lo,hi,k,j,'gain')
                assert not j or self.U(k,j)<=self.U(k,j-1),(lo,hi,k,j,'stock')
        self.l={32:Q(0)}
        for k in range(31,1,-1):
            self.l[k]=down((self.m/k+k*self.l[k+1])/(hi+k+Q(1,k)))
        self.g=self.m-self.U(2,0)

    def tail(self,k,j):
        return (1+Q(1,k))*min(Q(1,j+1),1/(1+self.lo*k))

    def U(self,k,j):
        return self.rows.get((k,j)) if k<=self.stop else self.tail(k,j)

    def L(self,k):
        return self.l.get(k,Q(0))

    def boost(self):
        """Lower serial old-host bonus; use the same finite harmonic surrogate."""
        lo,hi=self.lo,self.hi
        S=FM['S']; D=40; K=96
        pol=[Q(1)]
        for j in range(1,D): pol.append(pol[-1]*(1+lo/j))
        e=sum(Q((-1)**j,factorial(j)) for j in range(18))
        nxt=[0]+[int(S*down(e*sum(pol[j-1-k]/factorial(k) for k in range(j)))) for j in range(1,D+1)]
        q,z=hi.numerator,hi.denominator
        for d in range(D-1,1,-1):
            if d==2: gain=Q(nxt[2],S)
            ret=0
            for k in range(K-1,1,-1):
                f=k*(d*(k-1)+1)
                ret=z*(d*S+f*ret)//(z*(d+f)+d*k*q)
            row=[ret*nxt[1]//S]
            for j in range(1,d+1): row.append(j*z*(d*nxt[j+1]+row[-1])//(d*q+j*(d+1)*z))
            nxt=row
        assert self.m==Q(nxt[1],S)
        loss=Q(nxt[0],S)
        H=self.m
        for k in range(63,1,-1):
            x=(1+Q(1,k))*gain+loss/2
            y=hi+Q(3,2)+Q(1,k)
            u=k-Q(1,2)
            H=max(self.m,down(min(x/y,(x+u*H)/(y+u))))
            self.bonuses[k]=H-self.m

    @cache
    def profile(self,d,c,r):
        l=d-c
        w=(l+Q(r))/(d+1);s=w/d;t=Q(1,2*d)
        v=(w+l*t)*self.g+w*self.bonuses.get(d+1,0)+l*t*self.bonuses.get(2,0)
        if l:
            v-=l*(max(0,t-s)*max(0,self.U(2,0)-self.L(d+1))
                  +max(0,s-t)*max(0,self.U(d+1,l-1)-self.L(2)))
        return s,v

    def ledger(self,a,b,fee,ft,D=32,transfer=Q(0)):
        # All chain-port bills have first claim; only actual direct ports
        # consume the baseline excess.
        z=Q(D,D+1)*self.g+(self.g-self.U(2,0))/2
        half=Q(D,D+1)*self.g/2+self.g/4-self.U(2,0)/(4*(D+1))
        tail_debit=(1+Q(1,D+1))/(D+1)
        if min(half,min(z,0)/2+(a-fee-transfer)*D/2)-tail_debit<2*a+b or tail_debit/2>b:
            return False,'tail'
        if self.g<=0 or Q(D,D+1)*self.g-tail_debit<=0: return False,'seed or large star'
        for d in range(2,D):
            if self.profile(d,0,0)[1]<=0:return False,('star',d)
            for c in range(1,d+1):
                if d==c==2: continue
                bill=ft if d==2 else fee
                U=self.U(d+1,d-c)
                failed=False
                for k in ((0,c) if Q(2*d-c,2*d*(d+1))*U<=b else range(c+1)):
                    for r in (Q(k,2),Q(c,2)):
                        s,C=self.profile(d,c,r)
                        x=C+a*(c-2)-b-k*bill-(c-k)*(max(0,s*U-b)+(transfer if d>=15 else 0))
                        if x<0:failed=True
                if failed:
                    if d>8:return False,('coarse core',d,c)
                    cut=14 if d<=4 else 6
                    for k in range(c+1):
                        for es in combinations_with_replacement((0,*range(2,cut+1)),c-k):
                            for t in (Q(0),Q(1,cut+1)):
                                r=Q(k,2)+sum(Q(1,e) if e else t for e in es)
                                s,C=self.profile(d,c,r)
                                excess=sum(max(0,max(0,s-(Q(1,d*e*(e+1)) if e else 0))*max(0,U-(self.L(e+1) if e else 0))-b) for e in es)
                                x=C+a*(c-2)-b-k*bill-excess+(transfer*es.count(0) if d<=3 else 0)
                                if x<0:return False,('core',d,c,k,es,t,x)
        end=[]
        for d in range(2,10):
            for c in range(1,d+1):
                if d==c==2: continue
                for r in (Q(1,2),Q(c,2)):
                    s,_=self.profile(d,c,r)
                    end.append((Q(1,d),self.L(d+1),s,self.U(d+1,d-c),ft if d==2 else fee))
        K=max(31,self.stop)
        for d in range(10,K):
            h=max(((d+l)*self.U(d+1,l)/(2*d*(d+1)) for l in range(d)))
            end.append((Q(1,d),self.L(d+1),None,h,fee))
        for r in (Q(0),Q(1,K)):
            end.append((r,0,None,(1+Q(1,K+1))*min(Q(1,2*(K+1)),1/((K+1)*(1+(K+1)*self.lo))),fee))
        u=self.U(3,0);v=max(0,u-self.L(3));g=self.g+self.bonuses.get(3,0)
        X=lcm(12,*(6*r.denominator for r,_,_,_,_ in end),*(t.denominator for _,_,t,_,_ in end if t is not None))
        Y=lcm(g.denominator,v.denominator,u.denominator,self.L(3).denominator,
              *(z.denominator for _,ell,_,h,f in end for z in (ell,h,f)))
        end=[(int(r*X)//6,int(max(0,u-ell)*Y),None if t is None else int(X*t),
              int(Y*h),int(Y*max(0,h-self.L(3))),int(Y*f)*X) for r,ell,t,h,f in end]
        slope=int(Y*v);credit=int(2*Y*g)
        @cache
        def port(s,x):
            _,fwd,t,h,rev,f=x
            return (max(s*fwd,h*X) if t is None else max(0,(s-t)*fwd,(t-s)*rev))-f
        for k in (1,2,3):
            worst=None
            for x,y in product(end,repeat=2):
                r,t=x[0],y[0]
                ss=(r+t,) if k==1 else (X//12+r,)+(X//6,)*(k-2)+(X//12+t,)
                cost=port(ss[0],x)+port(ss[-1],y)+slope*sum(abs(r-t) for r,t in zip(ss,ss[1:]))-credit*sum(ss)
                worst=cost if worst is None else max(cost,worst)
            if worst>b*X*Y:return False,('chain',k,b-Q(worst,X*Y))
        return True,None


if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('--start',type=int,default=0)
    parser.add_argument('--stop',type=int,default=473)
    args=parser.parse_args()
    assert 0<=args.start<args.stop<=473
    for j in range(args.start,args.stop):
        lo,hi=Q(j,512),Q(j+1,512)
        if j<36:
            box=Bounds(lo,hi,48)
            box.boost()
            ok,why=box.ledger(*map(Q,('.015','.015','.007','.022')),96,Q('.005'))
        else:
            box=Bounds(lo,hi,32)
            if j<128: box.boost()
            ok,why=box.ledger(*map(Q,('.011','.011','.007','.020')),64)
        assert ok,(lo,hi,why)
        if j in (35,127,472): print('PASS cells through',hi,flush=True)
    print(f'PASS {args.stop-args.start} cells [{Q(args.start,512)},{Q(args.stop,512)}]; normalized numerator >= .033 with at least two cores',flush=True)
