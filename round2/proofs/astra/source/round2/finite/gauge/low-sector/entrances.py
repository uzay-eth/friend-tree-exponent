"""Save outward discounted masses on the four-low model's three-core high exits."""
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path
import sys

from sector import Sector, UNIT, geometry


if __name__=='__main__':
    s=F('.565809')
    model=Sector(s,4,0,7)
    lo,hi,tail=model.occupation(16,27)
    # 1.351 is the independently verified beta-column model majorant.
    assert tail*F('1.351')<F(1,10**16)
    mass=defaultdict(lambda:[0,0])
    for l,h,edges in zip(lo,hi,model.exits):
        for tree,p in edges:
            if geometry(tree)[1:]!=(3,1):
                continue
            bound=mass[tree]
            bound[0]+=l*p.numerator//p.denominator
            bound[1]+=(h*p.numerator+p.denominator-1)//p.denominator
    dest=Path(sys.argv[1])
    with dest.open('x') as f:
        f.write(f'# mass units {UNIT}, 27 copy generations, s=.565809, beta seed\n')
        for tree,(l,h) in sorted(mass.items()):
            f.write(f'{l}\t{h}\t{tree!r}\n')
    print('three-core one-high labelled exits',len(mass),'written to',dest)
