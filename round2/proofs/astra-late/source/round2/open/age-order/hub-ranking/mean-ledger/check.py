"""Rational physical moment ledger; rates use actual opposite parent degree."""
from fractions import Fraction as Q
from functools import cache
from itertools import product
from pathlib import Path
from runpy import run_path

TIP = run_path(str(Path(__file__).resolve().parents[1]/'tip-moments/check.py'))


class Ledger:
    def __init__(self,p,middle):
        self.p=p
        self.middle=middle
        cap=TIP['cap']
        self.tips={64:cap(64,p,True)/(1+64*p)}
        for k in range(63,1,-1):
            self.tips[k]=(cap(k,p,True)/k+self.tips[k+1])/(p+1+Q(1,k))
        self.g=middle-self.U(2)
        self.good={32:Q(0)}
        for k in range(31,1,-1):
            self.good[k]=(middle/k+k*self.good[k+1])/(p+k+Q(1,k))
            assert self.good[k]>=self.good[k+1]

    @cache
    def U(self,k):
        if k in self.tips:
            return self.tips[k]
        return Q(1,k) if self.p==1 else Q(k+1,2*k*k)

    def L(self,k):
        return self.good.get(k,Q(0))

    @cache
    def profile(self,d,c,receive):
        l=d-c
        w=(l+Q(receive))/(d+1); s=w/d; t=Q(1,2*d)
        C=(w+l*t)*self.g-l*(max(0,t-s)*max(0,self.U(2)-self.L(d+1))
                           +max(0,s-t)*max(0,self.U(d+1)-self.L(2)))
        return s,C

    def chains(self,b,fee,terminal_fee):
        H=10
        end=[]
        for d in range(2,H):
            for c in range(1,d+1):
                if d==c==2: continue
                for r in (Q(0),Q(c-1,2)):
                    s,_=self.profile(d,c,Q(1,2)+r)
                    end.append((Q(1,d),self.L(d+1),s,self.U(d+1),terminal_fee if d==2 else fee))
        # Enclose all large endpoint reciprocals above the lower-bound cutoff.
        for d in range(H,31):
            end.append((Q(1,d),self.L(d+1),None,self.U(d+1)/(d+1),fee))
        for r in (Q(0),Q(1,31)):
            end.append((r,0,None,Q(1,32**2),fee))
        u=self.U(3);v=max(0,u-self.L(3))
        def port(s,x):
            _,l,t,own,fee=x
            f=max(0,u-l)
            return (max(s*f,own) if t is None else max(0,(s-t)*f,(t-s)*max(0,own-self.L(3))))-fee
        for k in range(1,4):
            def cost(x,y):
                r,t=x[0],y[0]
                S=((r+t)/6,) if k==1 else ((Q(1,2)+r)/6,)+(Q(1,6),)*(k-2)+((Q(1,2)+t)/6,)
                return port(S[0],x)+port(S[-1],y)+v*sum(abs(j-i) for i,j in zip(S,S[1:]))-2*self.g*sum(S)
            worst=max(cost(x,y) for x,y in product(end,repeat=2))
            assert b-worst>0
            print(f'PASS p=1 chain length {k}: slack > {int(10**9*(b-worst))}/10^9')

    def stars(self,floor,cutoff=64):
        assert self.g>floor  # seed: total serial weight one, edge coefficient zero.
        for d in range(2,cutoff):
            assert self.profile(d,0,0)[1]>floor
        # Large stars: centre credit minus centre debit alone suffices.
        assert Q(cutoff,cutoff+1)*(self.g-Q(1,cutoff+1))>floor

    def first(self):
        a,b,fee,terminal_fee=map(Q,('.0085','.008','.004','.010'))
        worst=min((C+a*(c-2)-b-k*price-c*max(0,s*self.U(d+1)-b),d,c,k,r)
                  for d in range(2,64) for c in range(1,d+1) if (d,c)!=(2,2)
                  for k,r in ((0,0),(0,Q(c,2)),(c,Q(c,2)))
                  for price in (terminal_fee if d==2 else fee,)
                  for s,C in (self.profile(d,c,r),))
        assert worst==(Q('.0005'),3,3,0,0)
        assert Q(1,65**2)<b
        z=Q(64,65)*self.g+(self.g-self.U(2))/2
        assert z>0 and self.L(2)>Q(1,65)
        # l/d>=1/2 uses the stock credit; l/d<1/2 uses c>d/2 instead.
        assert z/2-2*a-b>0
        assert (a-fee)*32-2*a-b>0
        print('PASS p=1 retained cores: exact minimum 1/2000 and high-degree tail')
        self.chains(b,fee,terminal_fee)
        self.stars(Q(1,40))
        assert a*2+b==Q(1,40)

    def second(self):
        gain=self.g-Q(1,50)
        for d in range(2,64):
            for c in range(1,d+1):
                for r in (0,Q(c,2)):
                    s,C=self.profile(d,c,r)
                    E=C-c*s*self.U(d+1)-Q(1,50)*(d*s+Q(d-c,2*d))
                    assert E >= 0 and (E>0 or c==d and r==0)
        U=Q(66,2*65**2)
        assert gain>U and Q(64,65)*(gain-U)+(gain-self.U(2))/2>0
        self.stars(Q(3,100))
        print('PASS p=2 sitewise discharge with 1/50 of serial weight retained; stars and seed')


if __name__=='__main__':
    Ledger(1,Q('.56319')).first()
    Ledger(2,Q('.50833')).second()
