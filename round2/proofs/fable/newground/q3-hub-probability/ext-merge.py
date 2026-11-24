import collections
import pathlib
import re
import sys
from fractions import Fraction as Q

vals=collections.defaultdict(lambda:[Q(0),Q(0)])
masses=collections.defaultdict(int)
counts=collections.defaultdict(lambda:[0,0])
den=None
for stem in sys.argv[1:]:
    lines=pathlib.Path(stem+'.log').read_text().splitlines()
    for line in lines:
        t=line.split()
        if t[0]=='m':
            m,mode,nf,nr,d=map(int,t[1::2]); assert den is None or den==d; den=d
            counts[m,mode][0]+=nf; counts[m,mode][1]+=nr
        else:
            lo,hi=map(Q,t[-2:]); a=vals[m,mode,t[0],int(t[1])]; a[0]+=lo; a[1]+=hi
    run=pathlib.Path(stem+'.run').read_text()
    for i,u,z,d in re.findall(r'^mass (\d+) (\d+) (\d+) (\d+)$',run,re.M):
        assert int(d)==den; masses[int(i),int(u)]+=int(z)
assert masses and all(z==den for z in masses.values())
for (m,i),(nf,nr) in sorted(counts.items()):
    expected={10:(106,719),22:(5623756,97055181),23:(14828074,268282855),24:(39299897,743724984)}
    assert (nf,nr)==expected[m]
    print(f'm {m} improved {i} free {nf} rooted {nr} denominator {den}')
    for (mm,ii,t,u),(lo,hi) in sorted(vals.items()):
        if (m,i)==(mm,ii):
            assert (i,u) in masses
            print(t,u,f'{float(lo):.12g}',f'{float(hi):.12g}',lo,hi)
