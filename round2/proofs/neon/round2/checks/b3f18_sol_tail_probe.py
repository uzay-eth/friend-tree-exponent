"""Compare GMP tail row results against independently recovered rational controls."""
from fractions import Fraction as Q
from itertools import product
from random import Random
import subprocess
from b3f18_sol_algebra import B, controls


def fit(xs,ys):
    # Invert the Vandermonde interpolation by Newton divided differences.
    a=list(ys)
    for j in range(1,len(xs)):
        for i in range(len(xs)-1,j-1,-1): a[i]=(a[i]-a[i-1])/(xs[i]-xs[i-j])
    p=[a[-1]]
    for x,v in reversed(list(zip(xs[:-1],a[:-1]))):
        p=[v-x*p[0]]+[b-x*c for b,c in zip(p,p[1:])]+[p[-1]]
    return p


rng=Random(3023)
stream=[]; wants=[]
for d,mode,n,h,overflow,xfer in product((10,18), (1,2),range(3),range(2),range(2),range(2)):
    if overflow and xfer: continue
    # A synthetic nonoverflow cutoff above D0 exercises the fractional tmax branch.
    kt=d+2; k=kt+1 if overflow else kt-rng.randrange(4)
    c=k+h; ext=c-n
    for trial in range(4):
        F=list(map(Q,(rng.randrange(-10**11,10**11) for _ in range(4))))
        Fb=[F[0]]+list(map(Q,(rng.randrange(-10**11,10**11) for _ in range(3))))
        Up=list(map(Q,(rng.randrange(-10**11,10**11) if n else 0 for _ in range(4))))
        d0,d1,d2,u3=Up; d3=u3-n*F[3]
        ds=[d0,d1,d2,d3]
        a=F[0]-F[1]; b=F[1]-F[2]; r0=Q(2,17); r1=Q(2,17)+Q(n,18)
        rb=104725116; Y=rng.randrange(-10**11,10**11); center=Q(rb,10**9)
        first=[]; second=[]; third=[]; opt=max if mode==1 else min
        if overflow:
            first=[Q(0),a]; second=[Q(0),2*b]; third=[Q(0)]
            # Fit the overflow numerator directly from gain at rational nodes.
            xs=[Q(j,4) for j in (1,2,3)]
            ys=[]
            for t in xs:
                s=1-t; z=d*t/(d+t)
                y=(d/t+1)*(B(F,z)-B(F,t))-3*t*s*s*a-3*t*t*s*2*b
                ys.append(y*(d+t)**2/t**3)
            third += [a/b for a,b in zip(controls(fit(xs,ys),2,Q(1)),controls([d*d,2*d,1],2,Q(1)))]
            choices=[opt(first)+d0/(3*d),opt(second)+(3*d1-d0)/(3*d),opt(third)+opt(((d0-3*d1+3*d2)/d,d3/d))]
        else:
            m=min(Q(1),Q(d,c))
            xs=[m*j/5 for j in (1,2,3,4,5)]
            for ri,re in product((r0,r1),(Q(0),Q(ext,2))):
                # Recover low-degree coefficients from evaluations with degree
                # 7 after rational denominator clearing, rather than use C's
                # expanded coefficient formulas.
                xx=[m*j/8 for j in range(1,9)]; vals=[]
                def direct(t):
                    s=1-t; D=d/t; z=d*t/(d+t)
                    value=(D-c+ri+re)*(B(F,z)-B(F,t))+(1-c/D)*(B(Fb,t)-B(F,t))+B(ds,t)/D
                    return value-(re-center)*(3*t*t*s*a/d+t**3*Y) if xfer else value
                for t in xx: vals.append(direct(t)*(d+t)**3)
                p=fit(xx,vals)
                one=p[1]/(3*d**3)
                two=(p[2]-one*(9*d*d-6*d**3))/(3*d**3)
                assert p[0]==0
                first.append(one); second.append(two)
                ys=[(direct(t)-3*t*(1-t)**2*one-3*t*t*(1-t)*two)*(d+t)**3/t**3 for t in xs]
                third += [a/b for a,b in zip(controls(fit(xs,ys),4,m),controls([d**3,3*d*d,3*d,1],4,m))]
            choices=[opt(z) for z in (first,second,third)]
        for kind in (2,3,4):
            args=[d,kt,kind,k,h,mode,n,ext,xfer,r0.numerator,r0.denominator,r1.numerator,r1.denominator,rb,Y,*F,*Fb[1:],*Up]
            assert all(int(i)==i for i in args)
            stream.append(' '.join(str(int(x)) for x in args)); wants.append(choices[kind-2])

res=subprocess.run(['./sol-b3f-tail-probe'],input='\n'.join(stream)+'\n',text=True,capture_output=True,check=True)
assert not res.stderr
got=list(map(Q,res.stdout.splitlines()))
assert len(got)==len(wants)
for i,(x,y) in enumerate(zip(got,wants)):
    assert x==y,(i,stream[i],x,y)
print(f'C/GMP tail probe: {len(got)} exact T1/T2/T3 extrema agree with independently recovered Fraction controls')
