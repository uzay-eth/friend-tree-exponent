"""Independent small-catalogue checks and positive controls for search.py."""
from fractions import Fraction
from itertools import product
from math import gcd

import mpmath as mp
import search as s

mp.mp.dps=95


def check_mode(mode, h, width):
    hh=(h,)*width
    rows=s.coefficient_rows(hh[1:],mode)
    passing=sum(gcd(*map(int,row),a0)==1 for row in rows for a0 in range(-h,h+1))
    assert passing==s.row_count(hh,mode)
    print(f'{mode}: exhaustive height {h} primitive count {passing}: OK')


for mode,width in [('poly',4),('field',5),('kernel',6),('sparse',4)]:
    check_mode(mode,2,width)

# Direct integration independent of the hypergeometric representation.
for z in (-1,.5,1,2,3):
    x=mp.mpf('.56572366')
    assert abs(s.fz(z).f(x)-mp.quad(lambda t:t**x*mp.exp(z*(t-1)),[0,1]))<mp.mpf('1e-90')
print('F[z] versus independent quadrature at all 5 z: OK')

t=s.kernel_templates([s.F],1)[0]
ws = [s.Window(mp.mpf(x),mp.mpf('2.1e-8')) for x in ('0.549734972504','0.771192302407','0.56572366')]
got=[[] for _ in ws]
assert s.scan_template(t,ws,got)==0
for w,hits in zip(ws,got):
    brute=[]
    lo,hi=w.x-w.radius,w.x+w.radius
    for q0,q1,q2,p1,p2,p0 in product((-1,0,1),repeat=6):
        q=[v for v in (q0,q1,q2) if v]
        if not q or q[0]<0:
            continue
        f=lambda x:(p0+p1*x+p2*x*x)+(q0+q1*x+q2*x*x)*s.F.at(x)
        if f(lo)*f(hi)<0:
            brute.append(mp.findroot(f,(lo,hi)))
    assert sorted(mp.nstr(r,70) for r in brute)==sorted(mp.nstr(h.root,70) for h in hits)
assert any(h.formula=='-1 +1*F[1] +1*x = 0' for h in got[0])
assert any(h.formula=='-1 +1*F[1] +1*x^2 = 0' for h in got[1])
print('Independent height-1 HKRZ enumeration; PAN, DAN and friend centers:',[len(h) for h in got],': OK')

h=100
w=s.Window(mp.mpf('0.5723'),mp.mpf('0.0011'))
n,hitlists=s.scan_rational([w],h)
ratios={Fraction(p,q) for q in range(2,h+1) for p in range(1,q)}
want={r for r in ratios if Fraction('0.5712')<=r<=Fraction('0.5734')}
assert n==len(ratios) and {Fraction(h.formula) for h in hitlists[0]}==want
print(f'Independent rational sieve: {n} reduced rational candidates, {len(want)} hits: OK')

# Exact rejection of the degree-6 near relation noted by PSLQ. The derivative
# is Lipschitz with the coefficient-weighted bound on |P''| for 0<x<1.
a=[34600289,-95405394,60930212,-26222804,-36140252,52814813,60233216]
lo,hi=map(Fraction,(s.IR_LO,s.IR_HI))
mid,r=(lo+hi)/2,(hi-lo)/2
p=sum(v*mid**i for i,v in enumerate(a))
d=sum(i*v*mid**(i-1) for i,v in enumerate(a) if i)
m2=sum(i*(i-1)*abs(v) for i,v in enumerate(a))
assert abs(p)>abs(d)*r+m2*r*r/2
print('Exact rational Taylor bound rules out the PSLQ degree-6 near relation throughout the certified IR-0 bracket: OK')
