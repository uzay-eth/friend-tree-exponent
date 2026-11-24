"""Two actual marks, four birth-degree bounds used only for the Lyapunov score."""
from fractions import Fraction as F
from math import ceil

TAGS=(2,4,8,16)
FRONT={2:F(1),4:F('0.963'),8:F('0.931'),16:F('0.974')}
A=F(1,2)
B=F(1,16)
DAMP=F(7,50)
Q=F(2,5)
CUTOFF=40
UNIT=10**9

def tag(D):
    return max(i for i in TAGS if i<=D)

def actual(D):
    return A if D<16 else B

def credit(D):
    a=actual(D)
    return a*FRONT[tag(D)]/(a+DAMP)

def verify(verbose=True):
    tail=F(ceil(UNIT*credit(CUTOFF)/Q),UNIT)
    tables={}
    for i in TAGS:
        h={(CUTOFF,c):tail for c in range(1,CUTOFF+1)}
        aa=actual(i)
        pre=F(1,i)*(1+(A+DAMP)/aa)
        for D in range(CUTOFF-1,1,-1):
            for c in range(D,0,-1):
                ls=D-c;bud=F(ls,D)
                nxt=h.get((D,c+1),0)+credit(D)
                vals=[((ls+e)*h[D+1,c]+bud*nxt)/(Q+ls+e+bud)
                      for e in (F(0),F(c-1,2)+pre)]
                h[D,c]=F(ceil(UNIT*max(vals)),UNIT)
        initial=h[2,1]
        assert initial<=FRONT[i]
        h[2,1]=FRONT[i]
        for D in range(2,CUTOFF):
            for c in range(1,D+1):
                for e in (F(0),F(c-1,2)+pre):
                    v=h[D,c]
                    drift=(D-c+e)*(h[D+1,c]-v)+F(D-c,D)*(h.get((D,c+1),0)-v+credit(D))
                    assert drift<=Q*v
        assert Q*tail>=credit(CUTOFF)
        assert max(h[D,1]+credit(D) for D in range(2,CUTOFF+1))<=2
        tables[i]=h
        if verbose:
            print('virtual tag',i,'real mark',aa,'required',initial,'allowed',FRONT[i])
    lo=min(min(t.values()) for t in tables.values())
    hi=max(max(t.values()) for t in tables.values())
    delta=F(113,200)-Q-DAMP
    assert delta*Q*hi/(B+DAMP)<=lo
    if verbose:
        print('checked endpoint rows',len(TAGS)*sum(2*D for D in range(2,CUTOFF)),
              'min',lo,'max',hi,'rho',Q+DAMP,'size exponent',delta)
    ratio=max(F(1),max(t[D+1,c]/t[D,c] for t in tables.values()
                        for D in range(2,CUTOFF) for c in range(1,D+1)))
    return tables,ratio

if __name__=='__main__':
    verify()
