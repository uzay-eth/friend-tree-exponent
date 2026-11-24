from fractions import Fraction as F
from functools import cache
from math import prod

@cache
def hub_lo(l,c):
    a=[F(1)]
    for k in range(1,25): a.append(a[-1]/(c+k))
    return sum(a[:l])/(sum(a)+a[-1]/(c+24))

@cache
def promo_table(c,K=60):
    # Region with k leaves, one marked, >= c cores; first-step majorant
    # for promotion (boundary k=K+1 uses elementary 1/k bound).
    f={K+1:F(1,K+1)}
    A,B={},{}
    for k in range(1,K+1):
        r=k*(k+c)
        den=r+k-(k-1)*A.get(k-1,0)
        A[k]=r/den
        B[k]=(1+(k-1)*B.get(k-1,0))/den
        if k==1: # ensure Fraction when all input terms in row 1 were integers
            A[k]=F(r,r+k); B[k]=F(1,r+k)
    for k in range(K,0,-1): f[k]=A[k]*f[k+1]+B[k]
    assert all(f[k]>=f[k+1] for k in range(1,K+1))
    return f

def uleaf(l,c):
    if c: return min(F(1,l),promo_table(c)[l])
    # Star: first event at rate l+1; either marked leaf, another leaf,
    # or one birth; iterate using base F=1/(l-1) above cutoff 40.
    z=F(1,40)
    for k in range(40,l-1,-1):
        z=(F(1,k)+F(k-1,k)*uleaf(k-1,1)+k*z)/(k+1)
    return z
uleaf=cache(uleaf)

S=10**10
def lower(x): return (x*S).__floor__()
def upper(x): return (x*S).__ceil__()

def pref(N):
    inc=[[] for _ in range(N+1)]
    for i in range(N+1):
        for j in range(N+1):
            inc[i].append(lower(hub_lo(i,j)))
    return inc

if __name__=='__main__':
    for n in (3,4,5,9,10,11):
        with open(f'exact-{n}.in','w') as f:
            for row in pref(n): print(*row,file=f)
            print(*(upper(uleaf(l,0)) if l>=2 else S for l in range(n+1)),file=f)
            for c in range(n+1):
                print(*(upper(uleaf(l,c)) if l and c else S for l in range(n+1)),file=f)
            print(*(upper(F(1,max(1,l-1))) for l in range(n+1)),file=f)
