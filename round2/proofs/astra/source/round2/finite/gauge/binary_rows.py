"""Certified two-mark row and row-difference bounds using scalar root tables."""
from fractions import Fraction as F
from functools import cache
from math import ceil, floor
from binary_weight import A,BETA,DAMP,mark,solve

Q=F(49,125)
S=F(113123,200000)
N=32
UNIT=10**12
H,KAPPA=solve(64,Q,False)
VT=None
def z(D):
    a=mark(D)
    if VT is None:
        return a/(a+DAMP)
    return virtual_weight.credit(D)
def eta(D,c):
    return F(1,D*(D+1)) if c==1 else F(0)
def g(col,D):
    if isinstance(col,tuple):
        return col[0] if mark(D)==A else col[1]
    return int(col==0 or col==mark(D))

@cache
def alpha(col):
    ratios=[]
    profiles=[(A,64,H)] if VT is None else [(mark(i),virtual_weight.CUTOFF,h) for i,h in VT.items()]
    for aa,cut,h in profiles:
        for D in range(2,cut):
            for c in range(1,D+1):
                v=h[D,c]
                dd=h[D+1,c]-v
                bud=F(D-c,D)
                dc=(h[D,c+1]-v if c<D else 0)
                worst=(D-c)*dd + F(c-1,2)*max(0,dd) + bud*(dc+z(D))
                slack=(S-DAMP)*v-worst
                assert slack>0
                ratios.append(bud*g(col,D)/(v+slack/(aa+DAMP)))
        v=h[cut,1]
        ratios.append(g(col,cut)/(v+((S-DAMP)*v-z(cut))/(aa+DAMP)))
    return max(ratios)

def tables(lam,outer,reward,terminal):
    """Each reward and terminal pair is in already outward-rounded integer units."""
    t={}
    for D in range(N,1,-1):
        for c in range(D,0,-1):
            if D==N:
                t[D,c]=terminal(c)
                continue
            L=D-c
            b=F(L,D)
            lpay,upay=reward(D,c)
            lg,ug=t[D+1,c]
            lb,ub=t.get((D,c+1),(0,0))
            los,his=[],[]
            for e in (F(0),F(c-1,2)+outer):
                gain=L+e
                div=lam+gain+b
                los.append((lpay+gain*lg+b*lb)/div)
                his.append((upay+gain*ug+b*ub)/div)
            t[D,c]=floor(min(los)),ceil(max(his))
    return t

def direct_pair(lam,tolls):
    """Toll bounds are functions of the root degree, paid at rate b."""
    def terminal(c):
        b=F(N-c,N)
        lo,hi=tolls(N)
        return floor(UNIT*lo*b/(lam+eta(N,c))),ceil(UNIT*hi*(1/lam-(1-b)/(lam+1)))
    def reward(D,c):
        b=F(D-c,D)
        lo,hi=tolls(D)
        return floor(UNIT*b*lo),ceil(UNIT*b*hi)
    return tables(lam,0,reward,terminal)

@cache
def first(tag,col):
    lam=S+tag
    coef=KAPPA*alpha(col)/(lam-Q)
    return direct_pair(lam,lambda D:(F(g(col,D)),g(col,D)+coef*z(D)))

@cache
def child(outer,tag,col,forced):
    """Before a child copy, old-grandchild cost is handled separately."""
    lam,v=S+outer,S+tag
    c1=KAPPA*alpha(col)/(v-Q)
    before=first(tag,col)
    def terminal(c):
        b=F(N-c,N)
        e=eta(N,c)
        lo=g(col,N)*b/((lam+e)*(v+e))
        hi=((g(col,N)+c1*z(N))/v+c1*z(N)/(lam-Q))/lam
        return floor(UNIT*lo),ceil(UNIT*hi)
    def reward(D,c):
        b=F(D-c,D)
        lo,hi=before[D,c]
        if forced:
            lop,hip=before[D+1,c]
            lo,hi=min(lo,lop),max(hi,hip)
        return lo,hi+ceil(UNIT*c1*b*z(D)/(lam-Q))
    return tables(lam,tag if forced else F(0),reward,terminal)

@cache
def second(outer,col):
    def toll(D):
        tag=mark(D)
        lo,hi=child(outer,tag,col,True)[2,1]
        return g(col,D)+tag*F(lo,UNIT),g(col,D)+tag*F(hi,UNIT)
    return direct_pair(S+outer,toll)

def show(lo,hi):
    return f'[{floor(lo*UNIT)}/{UNIT},{ceil(hi*UNIT)}/{UNIT}]'

if __name__=='__main__':
    import argparse
    p=argparse.ArgumentParser()
    p.add_argument('--s',type=F,default=S)
    p.add_argument('--degree',type=int,default=N)
    p.add_argument('--virtual',action='store_true')
    a=p.parse_args()
    S,N=a.s,a.degree
    if a.virtual:
        import virtual_weight
        VT,KAPPA=virtual_weight.verify(False)
        Q,DAMP=virtual_weight.Q,virtual_weight.DAMP
    assert S>Q+DAMP and N>=16
    print('s',S,'degree cutoff',N,'binary weight',Q+DAMP)
    for col in (A,BETA):
        for tag in (A,BETA):
            lo,hi=second(tag,col)[2,1]
            print('true K',tag,col,show(F(lo,UNIT),F(hi,UNIT)))
        lo,hi=child(BETA,A,col,False)[2,1]
        print('true row difference',col,show((A-BETA)*F(lo,UNIT),(A-BETA)*F(hi,UNIT)))
        for outer in (A,BETA):
            lo,hi=child(outer,BETA,col,False)[2,1]
            print('high-degree continuation P',outer,col,
                  show(BETA*F(lo,UNIT),BETA*F(hi,UNIT)))
    for x in (F(1,3),F(1,5)):
        col=(x,F(1))
        lo1,hi1=second(A,col)[2,1]
        lo2,hi2=second(BETA,col)[2,1]
        print('binary Collatz, first component',x,
              show(min(lo1/x,lo2)/UNIT,max(hi1/x,hi2)/UNIT))
