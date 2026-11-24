#!/usr/bin/env python3
"""Independent integer bounds; no C++ arithmetic or witness table scans."""
from fractions import Fraction
from math import comb
from pathlib import Path
from struct import unpack
import sys


def ceildiv(n,d): return -((-n)//d)


with Path(sys.argv[1]).open('rb') as f:
    h=unpack('<Q63q',f.read(512))
assert h[0]==0x4654325F55355631 and h[4:6]==(1<<40,1<<30)

W=16<<40; E=32<<40; D=4096; SR=1<<30; SUB=1024
i128=(1<<127)-1; i64=(1<<63)-1
# Finite lower f at seam; each numerator term bounded separately. Positive
# denominator slack is added before integer ceiling for signed 128 checks.
df=3*W*D**2+W*D
seam=df*2*D*SUB+2*D**2-1
bud=4*W*D**2*SUB+D**2-1
kill=(1<<40)*W*SUB+(1<<40)-1
center=2*W*256*SR*SUB+SR-1  # centered rows are below d2, never at seam
env3=128**3*(2*W+E)
src=(256*SR)*(2*E)*SUB+256*SR-1
for key,v in dict(f_lower_seam=seam, f_bud=bud, f_kill=kill, f_center=center,
                  env3_difference=env3, env_source=src).items():
    assert 0<v<=i128
    print(key,v,'bits',v.bit_length())

# Legal extreme values simultaneously realize the sign of each lower-tail
# seam term at c=0, D_old=4095. Compute the directed cap independently.
df_n=-W*D**2-W*D**2-W*D-W*D**2
n=-2*(D-1)*df_n; d=2*D**2
subticks=ceildiv(n*SUB,d); ticks=ceildiv(subticks,SUB)
exact=Fraction(n,d)
assert ticks==ceildiv(exact.numerator,exact.denominator) and ticks<=i64
assert ceildiv(-n*SUB,d)==-(n*SUB//d)
print('lower seam corner df_n',df_n,'scaled signed gain numerator',n*SUB,
      'denominator',d,'subticks',subticks,'directed ticks',ticks)


def multiplicity(a,b,c,d,e,f,g):
    return comb(a,2)+comb(b,2)*(b-1)+comb(c,3)+comb(d,2)*(d-1)**2+e*comb(e-1,2)*(f-1)+comb(g,4)


for name,ks in [('production',(h[14],h[18],h[21],h[25],h[28],h[30],h[33])),
                ('universal',(32,24,24,16,16,16,16))]:
    m=multiplicity(*ks)
    w=(2+m)*W+(D-2)*(1<<40)
    hc=m*i64
    assert w<=i64 and hc<=i128
    print(name,'arity',ks,'multiplicity',m,'absolute shifted comparison bound',w,
          'higher signed128 aggregate before checked narrowing',hc)
