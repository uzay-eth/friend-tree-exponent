"""Exact consequences of the two-moment Schur enclosures."""
from fractions import Fraction as F

from pathlib import Path
import re
import sys


def verify_inputs(path):
    """Check exactly the four rounded constants used below against the recorded true rows."""
    active=None
    lines={}
    for line in path.read_text().splitlines():
        if line.startswith('$ '):
            match=re.search(r'--s (\S+) --test (\S+)',line)
            active=(F(match[1]),match[2]) if match else None
        elif active:
            match=re.fullmatch(r'(true row 1/2|true positive dormant row) \[([^,]+),([^\]]+)\]',line)
            if match:
                lines[active+(match[1],)]=(F(match[2]),F(match[3]))
    assert lines[F('.5665'),'A','true row 1/2'][0]>=F('.5225363')
    assert lines[F('.5665'),'B','true row 1/2'][0]>=F('.0706297')
    assert lines[F('.5665'),'A','true positive dormant row'][0]>=F('.2995629')
    assert lines[F('.54001'),'B','true positive dormant row'][1]<F('4.002')
    print('verified all four inputs against',path.name)


if len(sys.argv)>1:
    verify_inputs(Path(sys.argv[1]))
u,anchor = F('.5665'),F('.54001')
kA,kB,dA = F('.5225363'),F('.0706297'),F('.2995629')
critical=1-(1+dA)*kB/(1-kA)
assert critical<F('.807760')
print('critical excess upper',critical, '< .807760')
for s,cap in ((F('.5630'),F('.9980')),(F('.5631'),F('.9920')),
              (F('.564'),F('.9395')),(F('.565'),F('.8844')),
              (F('.565615'),F('.8522')),(F('.565809'),F('.8423'))):
    p=(u-s)/(u-anchor)
    assert F('4.002')**p.numerator*F('.807760')**(p.denominator-p.numerator)<cap**p.denominator
    hmin,w=F('.501218471'),F('.974')/(F(1,16)+F('.14'))
    factor=1+w/((1-cap)*hmin)
    print('s',s,'positive dormant return <',cap,'time bound multiplier <',factor.__ceil__())
for s,cap in ((F('.565'),F('.8844')),(F('.565615'),F('.8522')),
              (F('.565809'),F('.8423'))):
    # F = k_A + k_beta (1+d_A) sum d_beta^j.
    # The deep anchor bounds apply at later s by chronological monotonicity.
    err=F('.836')*(1+F('4.685'))*cap**100/(1-cap)
    print('s',s,'100-step scalar tail <',(10**12*err).__ceil__(),'/ 1000000000000')

l=F('.565')
row_floor=kA+kB
for s,D,rlo,rhi,Flo,Fhi in (
        tuple(map(F,('.565615','.8522','.95919','1.05819','.8468','1.2414'))),
        tuple(map(F,('.565809','.8423','.94667','1.04515','.8090','1.1746')))):
    low_power=(s-l)/(l-anchor)
    high_power=(u-s)/(u-anchor)
    assert F('5.435')**low_power.numerator*rlo**low_power.denominator < 1
    assert F('5.435')**high_power.numerator < rhi**high_power.denominator
    assert Flo<1-(1-rlo)*(1+(1-row_floor)/(1-D))
    assert 1+(rhi-1)*(1+(rhi-row_floor)/(1-D))<Fhi
    print('s',s,'prior-assisted scalar between',Flo,Fhi)
