"""Exact interval enclosure and all-degree moment discharge on [1,3]."""
from fractions import Fraction as Q
from functools import cache
from hashlib import sha1
from itertools import product
from math import factorial, lcm, prod
from pathlib import Path
import sys


TABLE = Path(__file__).resolve().parents[2]/'mellin/lower-1-3.tsv'
TABLE_BLOB = 'defa8793b129e4b3c3f5ef6f23fd4033fc4ff8fc'
SCALE = 10**12
D = 32
A,B,F,FT = map(Q, ('.0085','.008','.004','.0131'))


def up(x):
    return Q(-((-x*SCALE)//1),SCALE)


def down(x):
    return Q((x*SCALE)//1,SCALE)


@cache
def integer_moment(r,j):
    return prod(r+i for i in range(j))/factorial(j)


def chord(r,p):
    assert 0<r<=1
    j=p//1
    return (1-p+j)*integer_moment(r,j)+(p-j)*integer_moment(r,j+1)


def ceiling(r,lo,hi):
    """A_p(r) bounded using the decreasing Beta moment or its shifted chord."""
    if r<=1:
        return chord(r,lo)
    s=r-1
    return (1+hi/s)*chord(s,lo)


class Interval:
    def __init__(self,lo,hi,M):
        self.M=M
        cap={}
        for k in range(64,1,-1):
            a=1+Q(1,k)
            p=1+lo+a; q=1+hi+a
            large=ceiling(a,lo,hi)
            small=(ceiling(a-Q(1,2),lo,hi)+ceiling(1+Q(1,k+1),lo,hi))/2
            cap[k]=up((1-1/q)*large + small/p)
        self.u={64:min(Q(1,64),up(cap[64]/(1+64*lo)))}
        for k in range(63,1,-1):
            self.u[k]=min(Q(1,k),up((cap[k]/k+self.u[k+1])/(1+lo+Q(1,k))))
        self.g=M-self.U(2)
        self.l={32:Q(0)}
        for k in range(31,1,-1):
            self.l[k]=down((M/k+k*self.l[k+1])/(hi+k+Q(1,k)))

    def U(self,k):
        return self.u.get(k,Q(1,k))

    def L(self,k):
        return self.l.get(k,Q(0))

    def profile(self,d,c,receive):
        l=d-c
        w=(l+receive)/(d+1); s=w/d; t=Q(1,2*d)
        C=(w+l*t)*self.g-l*(max(0,t-s)*max(0,self.U(2)-self.L(d+1))
                           +max(0,s-t)*max(0,self.U(d+1)-self.L(2)))
        return w,s,C

    def stars(self,floor):
        return (self.g>=floor and min(self.profile(d,0,Q(0))[2] for d in range(2,D))>=floor
                and Q(D,D+1)*(self.g-Q(1,D+1))>=floor)

    def sitewise(self):
        h=self.g-B-2*A
        u=Q(1,D+1)
        if not (h>u and Q(D,D+1)*(h-u)+(h-self.U(2))/2>=0 and self.stars(B+2*A)):
            return False
        # At receive zero the stock fraction can be cancelled once l>0;
        # at receive c/2 check the actual centre direction (piecewise-linear C in receive).
        for d in range(2,D):
            for c in range(1,d+1):
                for r in (Q(0),Q(c,2)):
                    w,s,C=self.profile(d,c,r)
                    if C < c*s*self.U(d+1)+(B+2*A)*(w+Q(d-c,2*d)):
                        return False
        return True

    def chains(self):
        end=[]
        for d in range(2,10):
            for c in range(1,d+1):
                if d==c==2: continue
                for r in (Q(1,2),Q(c,2)):
                    _,s,_=self.profile(d,c,r)
                    end.append((Q(1,d),self.L(d+1),s,self.U(d+1),FT if d==2 else F))
        for d in range(10,31):
            end.append((Q(1,d),self.L(d+1),None,self.U(d+1)/(d+1),F))
        for r in (Q(0),Q(1,31)):
            end.append((r,0,None,Q(1,32**2),F))
        u=self.U(3); v=max(0,u-self.L(3))
        # Clear denominators once for the Cartesian product. All the tests
        # below still compare the exact same rational expressions.
        xs=lcm(12,*(6*x.denominator for x,_,_,_,_ in end),
               *(t.denominator for _,_,t,_,_ in end if t is not None))
        ys=lcm(self.g.denominator,v.denominator,u.denominator,self.L(3).denominator,
               *(z.denominator for _,ell,_,own,fee in end for z in (ell,own,fee)))
        end=[(int(x*xs)//6,int(max(0,u-ell)*ys),None if t is None else int(t*xs),
              int(own*ys),int(max(0,own-self.L(3))*ys),int(fee*ys)*xs)
             for x,ell,t,own,fee in end]
        gamma=int(2*self.g*ys)
        slope=int(v*ys)

        @cache
        def port(s,x):
            _,f,t,own,reverse,fee=x
            return (max(s*f,xs*own) if t is None else max(0,(s-t)*f,(t-s)*reverse))-fee

        worst=None
        for k in range(1,4):
            for x,y in product(end,repeat=2):
                r,t=x[0],y[0]
                if k==1: S=(r+t,)
                else: S=(xs//12+r,)+(xs//6,)*(k-2)+(xs//12+t,)
                cost=port(S[0],x)+port(S[-1],y)+slope*sum(abs(j-i) for i,j in zip(S,S[1:]))-gamma*sum(S)
                if worst is None or cost>worst[0]: worst=cost,k,x,y
        return B-Q(worst[0],xs*ys),worst[1:]

    def suppressed(self):
        u=Q(1,D+1)
        if not (u*u<B and self.L(2)>u and (A-F)*D/2>2*A+B and self.stars(B+2*A)):
            return False,'tail precondition or star'
        z=Q(D,D+1)*self.g+(self.g-self.U(2))/2
        if z<2*A+B: return False,'tail stock floor'
        for d in range(2,D):
            fee=FT if d==2 else F
            for c in range(1,d+1):
                if d==c==2: continue
                ks=(0,c) if self.U(d+1)/(d+1)<=B else range(c+1)
                for k,r in ((k,r) for k in ks for r in (Q(k,2),Q(c,2))):
                    _,s,C=self.profile(d,c,r)
                    residual=C+A*(c-2)-B-k*fee-(c-k)*max(0,s*self.U(d+1)-B)
                    if residual<0: return False,('retained',d,c,k,r,float(residual))
        z,witness=self.chains()
        return z>=0, ('chain',float(z),witness[0])


def main(path):
    counts={'sitewise':0,'suppressed':0}
    bad=[]
    last=Q(1)
    data=path.read_bytes()
    blob=sha1(f'blob {len(data)}\0'.encode()+data).hexdigest()
    assert blob==TABLE_BLOB, (path,blob)
    print('PASS authenticated 1/128 middle table',blob,flush=True)
    for line in data.decode().splitlines():
        lo,hi,m=map(Q,line.split())
        assert lo==last and hi-lo==Q(1,128) and hi<=3 and 0<m<1
        last=hi
        I=Interval(lo,hi,m)
        if I.sitewise(): counts['sitewise']+=1; continue
        ok,info=I.suppressed()
        if ok: counts['suppressed']+=1
        else: bad.append((str(lo),str(hi),float(m),float(I.g),info))
    assert last==3
    print('PASS counts',counts,flush=True)
    for row in bad: print('UNCLOSED',row)
    assert not bad, f'{len(bad)} unclosed bins'
    print('PASS entire closed exponent interval [1,3], uniform normalized numerator 1/40')


if __name__=='__main__':
    main(Path(sys.argv[1]) if len(sys.argv)>1 else TABLE)
