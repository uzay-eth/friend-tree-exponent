"""Price a complete conditional donor_lift receipt under the shared physical birth law."""
import argparse
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
PARTS = {'minimum','mother','0','1','2','3'}


def read_entries(files, refine=None):
    """Join disjoint component shards from one cut; never sum across cuts."""
    header, rows = None, {}
    for path in files:
        data = json.loads(path.read_text())
        part = data.pop('conditional')
        if (header is not None and header != data) or rows.keys() & part.keys():
            raise ValueError('entry shards overlap or describe different cuts')
        header = data
        rows.update(part)
    if (rows.keys() != PARTS or header['grid'] != 10**18 or
            list(map(Q,header['discount'])) != [Q('.565673'),Q('.565775')] or
            header['channels'] != ['lower_lo','lower_hi','upper_lo','upper_hi']):
        raise ValueError('need one complete compatible six-component receipt')
    def labels(cap): return ('overflow', *map(str,range(2,cap+1)))
    base = labels(header['mother_cap'])
    fine = base
    if refine:
        data = json.loads(refine.read_text())
        extra = data.pop('conditional')
        cap = data.pop('mother_cap')
        if (data != {k:v for k,v in header.items() if k!='mother_cap'} or extra.keys() != {'mother'}
                or cap <= header['mother_cap']):
            raise ValueError('incompatible refined mother component')
        header['mother_cap'] = cap
        fine = labels(cap)
        rows['mother'] = extra['mother']
    for part,r in rows.items():
        if (set(r) != set(fine if part=='mother' else base) or
                any(len(v)!=4 or not all(type(x) is int for x in v) or v[0]>v[1] or v[2]>v[3]
                    for v in r.values())):
            raise ValueError('invalid class entries')
    return header,{d:tuple(sum(r.get(d,r['overflow'])[j] for r in rows.values()) for j in range(4))
                   for d in fine}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('files',type=Path,nargs='+')
    p.add_argument('--limit',type=int,default=18)
    p.add_argument('--refine-mother',type=Path,help='replace only the mother component with finer degree entries')
    a=p.parse_args()
    header,rows=read_entries(a.files,a.refine_mother)
    if a.limit < header['mother_cap']:
        p.error('LP cutoff must resolve each finite birth-mother class')
    spec=spec_from_file_location('physical_birth',HERE.parent/'birth-law'/'replay.py')
    bl=module_from_spec(spec)
    spec.loader.exec_module(bl)
    audit=bl.load('audit',bl.R2/'checks'/'stationary-p4-sol.py')
    lp=bl.build(a.limit,True,audit.class_kernel())
    grid=header['grid']
    scale=10**6
    lo,hi=map(Q,header['discount'])
    cap=header['mother_cap']
    print('degrees',header['degrees'],'arities',header['arities'],'mother',cap,
          'LP rows',len(lp.rows),'cols',len(lp.cols),flush=True)
    for j,name in [(0,'lower'),(2,'upper')]:
        def cell(c): return rows[str(c)] if 2<=c<=cap else rows['overflow']
        uniform=(max(0,lo*min(r[j] for r in rows.values())/grid),
                 hi*max(r[j+1] for r in rows.values())/grid)
        def cost(kind,m,lane):
            return scale*m[2]*Q(cell(m[0])[j+lane],grid) if kind=='occ' else 0
        l=lp.lower(lambda k,m:cost(k,m,0))/scale
        u=-lp.lower(lambda k,m:-cost(k,m,1))/scale
        ends=[s*v for s in (lo,hi) for v in (l,u)]
        box=max(uniform[0],min(ends)),min(uniform[1],max(ends))
        if box[0]>box[1]:
            raise ValueError('crossed certified pooled interval')
        print(name,'shared-class uniform',bl.fmt(uniform),'physical sided LP',bl.fmt(box),
              'exact',*map(str,box),flush=True)


if __name__=='__main__':
    main()
