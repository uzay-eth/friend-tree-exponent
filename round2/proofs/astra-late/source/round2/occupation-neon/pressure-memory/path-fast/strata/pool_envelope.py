"""Apply the physical sided birth law to one joint P4 support-envelope receipt."""
import argparse
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('receipt',type=Path)
    p.add_argument('--limit',type=int,default=18)
    a=p.parse_args()
    data=json.loads(a.receipt.read_text())
    cap,grid=data['mother_cap'],data['grid']
    lo,hi=map(Q,data['discount'])
    classes=('overflow',*map(str,range(2,cap+1)))
    if (data['conditional'].keys()!={'envelope'} or a.limit<cap or data['include_old_band'] or
            grid!=10**18 or (lo,hi)!=(Q('.565673'),Q('.565775')) or
            data['channels']!=['lower_lo','lower_hi','upper_lo','upper_hi']):
        p.error('need one complementary P4 envelope with all its classes resolved')
    rows=data['conditional']['envelope']
    if set(rows)!=set(classes) or any(len(v)!=4 or any(type(x) is not int for x in v) or
                                     v[0]>v[1] or v[2]>v[3] for v in rows.values()):
        p.error('invalid conditional entries')
    spec=spec_from_file_location('physical_birth',Path(__file__).resolve().parents[2]/'birth-law'/'replay.py')
    b=module_from_spec(spec)
    spec.loader.exec_module(b)
    audit=b.load('audit',b.R2/'checks'/'stationary-p4-sol.py')
    lp=b.build(a.limit,True,audit.class_kernel())
    scale=10**6
    for j,side in [(0,'lower'),(2,'upper')]:
        def cost(kind,m,lane):
            if kind!='occ':
                return 0
            row=rows[str(m[0])] if 2<=m[0]<=cap else rows['overflow']
            return scale*m[2]*Q(row[j+lane],grid)
        L=lp.lower(lambda k,m:cost(k,m,0))/scale
        H=-lp.lower(lambda k,m:-cost(k,m,1))/scale
        lower=max(0,min(lo*L,hi*L),min(lo*Q(x[j],grid) for x in rows.values()))
        upper=min(max(lo*H,hi*H),max(hi*Q(x[j+1],grid) for x in rows.values()))
        if lower>upper:
            raise ValueError('crossed envelope')
        print(side,b.fmt((lower,upper)),lower,upper,flush=True)


if __name__=='__main__':
    main()
