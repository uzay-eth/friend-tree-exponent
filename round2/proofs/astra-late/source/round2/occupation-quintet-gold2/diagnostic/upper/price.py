#!/usr/bin/env python3
"""Upper-only, original GOLD2 H3 complement of the D24 matched mask."""
import argparse
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SOURCE=ROOT/'higher/cohort.py'
PIN='b82b6768a91ee1ae66f319b078b9076d04878546705a3538795fff73ec06057c'


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--min-mother',type=int,default=80)
    p.add_argument('--mother',type=int,default=480)
    p.add_argument('--mother-stock',type=int,default=12)
    p.add_argument('--middle',type=int,default=2)
    p.add_argument('--middle-stock',type=int,default=2)
    p.add_argument('--young',type=int,default=12)
    p.add_argument('--young-stock',type=int,default=3)
    a=p.parse_args()
    assert 19<=a.min_mother<=a.mother<=1998 and 2<=a.middle<=5 and 2<=a.middle_stock<=4
    assert 2<=a.young<=64 and 2<=a.mother_stock<=15 and 1<=a.young_stock<=8
    a.side='upper'; a.sibling=a.upstream=2; a.outside_threshold=True; a.exclude_strict=False
    data=SOURCE.read_bytes();assert sha256(data).hexdigest()==PIN
    spec=spec_from_file_location('upper_gold2_cohort',SOURCE)
    module=module_from_spec(spec);exec(compile(data,str(SOURCE),'exec'),module.__dict__)
    with module.host.open_host('upper') as h:
        assert h.docs['upper_sha256']=='f1dcd7ad3e6df82a985eb3dd6fd793b5a3e38c78368d94e4b135ab174ab92876'
        rows=module.Rows(h)
        module.registration();module.physical_checks(rows)
        c=module.Cohort(rows,a)
        assert c.mask((79,8),(5,4),(18,3)) is None
        if a.min_mother<=80:
            assert c.mask((80,8),(5,4),(18,3)) is not None
        print('upper original GOLD2',h.docs['upper_sha256'],a,flush=True)
        floor=c.sweep(('hmy',),'physical')
        assert floor>0
        print('EXACT upper physical H3 complement',floor,'ticks/core',floor*2**40,'units/core',float(floor),flush=True)


if __name__=='__main__':main()
