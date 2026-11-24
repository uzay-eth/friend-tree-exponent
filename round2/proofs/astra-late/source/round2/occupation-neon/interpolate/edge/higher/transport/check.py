"""Independent exact check of a synthetic, signed tagged transport graph."""
from fractions import Fraction as Q
from functools import cache
from itertools import product
import subprocess
import sys

A,K=20,2
CAP=(A,2,2,3,3)
AR=(K,1,1,1,1)
U=10**18
J=(2,0)
S=(Q('.565673'),Q('.565775'))


def toll(x,family):
    if len(x)!=5 or x[0][0]<(2 if family else 18) or x[0][1]<2:return Q(0),Q(0)
    d=x[0][0]
    return Q(19-d,1000),Q(d-19,100)


@cache
def ref(x,upper,family):
    n=len(x)
    if any(d>CAP[i] or k>AR[i] for i,(d,k) in enumerate(x)):return Q(0)
    a,b=toll(x,family);rates=[]
    p=(-1,0,1,0,3) if family else (-1,0,1,2,3)
    for i,(d,k) in enumerate(x):
        def y(t):return x[:i]+(t,)+x[i+1:]
        if d==0:
            rates.append(((1-Q(k+1,CAP[i]+1),Q(1)),ref(y((0,k+1)),upper,family)))
            continue
        L=d-k-1
        g=ref(y((0,k) if i==4 and d==CAP[i] else (d+1,k)),upper,family)
        if i==3 and n==5:
            dd=x[4][0]
            rates.append(((Q(1,dd),)*2 if dd else (0,Q(1,CAP[4]+1)),g+b))
        fixed=L+sum(Q(1,x[j][0]) for j in range(min(n,4)) if p[j]==i or p[i]==j)
        unseen=k-p[:n].count(i)+int(i==0)
        assert unseen>=0
        rates.append(((fixed,fixed+Q(unseen,2)),g))
        if L:
            post=y((d,k+1));z=ref(post,upper,family)
            if n<5 and i==p[n]:z+=ref(post+(J,),upper,family)
            rates.append(((Q(L,d),)*2,z))
    vals=tuple(q for r,q in rates)
    change=max if upper else min
    raw=min if upper else max
    assert raw(S)>0
    # Enumerate all physical corner choices independently of the native sort.
    return change((a+sum(r*v for r,v in zip(rs,vals)))/(s+sum(rs))
                  for s in S for rs in product(*(set(r) for r,_ in rates)))


def check(binary,family):
    lines=[f'20 2 2 1 2 1 3 1 3 1 20 {family}']
    dom=[[(d,k) for d in range(CAP[i],1,-1) for k in range(min(d-1,AR[i]),-1,-1)] for i in (1,2,3)]
    for d in range(20,1 if family else 17,-1):
        for k,*x in product(range(3),*dom):
            a,b=toll(((d,k),*x,J),family);lines.append(f'{int(a*U)} {int(a*U)} {int(b*U)} {int(b*U)}')
    run=subprocess.run([binary],input='\n'.join(lines)+'\n',text=True,check=True,capture_output=True)
    ll,hh=map(int,run.stdout.split()[:2]);lo,hi=ref((J,),False,family),ref((J,),True,family)
    assert Q(ll,U)<=lo<=hi<=Q(hh,U)
    assert max(lo-Q(ll,U),Q(hh,U)-hi)<Q('1e-12') and lo<hi and not (lo==hi==0)
    print('PASS exact independent signed corner sweep', ('P4','Yb')[family],ref.cache_info().currsize,'reference states;',run.stdout.strip())
    print('native excess errors',float(lo-Q(ll,U)),float(Q(hh,U)-hi))


if __name__=='__main__':
    for family in (0,1):
        ref.cache_clear()
        check(sys.argv[1],family)
