"""Replay a streamed P4 stratum against the frozen exact donor_lift."""
import argparse
import ctypes as C
from fractions import Fraction as Q
from hashlib import sha256
import importlib.util
import json
from pathlib import Path
from time import process_time

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location('frozen_path_fast',HERE.parent/'replay.py')
base=importlib.util.module_from_spec(spec)
spec.loader.exec_module(base)
I,U=base.I,base.U
PARTS=base.PARTS


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--library',type=Path,required=True)
    ap.add_argument('--lower',type=Path,required=True)
    ap.add_argument('--upper',type=Path,required=True)
    ap.add_argument('--degrees',type=int,nargs=4,default=[6,5,4,6])
    ap.add_argument('--arities',type=int,nargs=4,default=[2,3,2,1])
    ap.add_argument('--mother-cap',type=int,default=3)
    ap.add_argument('--side-cap',type=int,default=3)
    ap.add_argument('--include-old-band',action='store_true')
    ap.add_argument('--budget',type=int,default=80000000)
    ap.add_argument('--parts',nargs='+',choices=(*PARTS,'envelope'),default=PARTS)
    ap.add_argument('--engines',nargs='+',choices=('frozen','plain','merged','strata'),default=['strata'])
    ap.add_argument('--reference',type=Path,help='exact conditional JSON from a matched frozen engine')
    ap.add_argument('--entries',type=Path)
    a=ap.parse_args()
    if len(set(a.parts))!=len(a.parts) or len(set(a.engines))!=len(a.engines) or not 0<a.budget<1<<63:
        ap.error('invalid parts, engines or budget')
    if any(not (0 <= k < d <= 16 and d >= 2 and k >= (i < 3))
           for i,(d,k) in enumerate(zip(a.degrees,a.arities))):
        ap.error('invalid path caps or D18 ownership boundary')
    if any(not (m==0 or 2<=m<=23) for m in (a.side_cap,a.mother_cap)):
        ap.error('invalid memory cap')
    if 'envelope' in a.parts and 'frozen' in a.engines:
        ap.error('frozen solver has no joint support envelope')
    if 'envelope' in a.parts and len(a.parts)!=1:
        ap.error('joint support envelope replaces the six-component sum')
    if base.GRID!=10**18 or (base.LOW,base.HIGH)!=(Q('.565673'),Q('.565775')):
        raise ValueError('frozen rates or grid changed')
    tag=dict(degrees=a.degrees,arities=a.arities,side_cap=a.side_cap,mother_cap=a.mother_cap,
             include_old_band=a.include_old_band,grid=base.GRID,discount=list(map(str,(base.LOW,base.HIGH))),
             channels=['lower_lo','lower_hi','upper_lo','upper_hi'])
    ref={}
    if a.reference:
        data=json.loads(a.reference.read_text())
        if any(data.get(k)!=v for k,v in tag.items()):
            raise ValueError('reference has different controls')
        ref=data['conditional']
        if not set(a.parts)<=ref.keys():
            raise ValueError('reference does not cover requested components')
    t=process_time()
    tables=tuple(base.original(p,side) for p,side in [(a.lower,'lower'),(a.upper,'upper')])
    tab,indices=base.native_tables(tables)
    dll=base.library(a.library)
    dll.strata.argtypes=(base.INDEX,base.INDEX,C.c_int,C.c_int,C.POINTER(base.Tables),
                         C.c_int,C.c_int,C.c_int,U,C.POINTER(I),C.POINTER(U),C.POINTER(I))
    dll.strata.restype=C.c_int
    arr=C.c_int*4
    args=(arr(*a.degrees),arr(*a.arities),a.side_cap,a.mother_cap,tab)
    print('authenticated originals',tag,f'auth_cpu={process_time()-t:.3f}',flush=True)
    newest={}
    for engine in a.engines:
        total=[[Q(0),Q(0)],[Q(0),Q(0)]]
        counts=[0,0,0]
        handle=None
        if engine=='frozen':
            handle=dll.prepare(*args,not a.include_old_band,base.LOAD())
            if not handle:
                raise RuntimeError('frozen prepare failed')
        try:
            for part in a.parts:
                t=process_time()
                ticks=(I*(4*max(1,a.mother_cap)))()
                span=(I*2)()
                component=-3 if part=='envelope' else PARTS.index(part)-2
                if handle:
                    stats=(U*9)()
                    if dll.price(handle,component,a.budget,ticks,stats,span):
                        raise RuntimeError(f'frozen solve failed: {part}')
                    paths,n,ne,nv,nr,slots,c,h,s=stats
                    fields=f'paid={paths} states={n} edges={ne} variable={nv} residuals={nr} slots={slots}'
                else:
                    stats=(U*11)()
                    mode={'plain':0,'merged':2,'strata':3}[engine]
                    if dll.strata(*args,a.include_old_band,component,mode,a.budget,ticks,stats,span):
                        raise RuntimeError(f'{engine} solve failed: {part}')
                    size,terminal,pay,absent,tagged,n,ne,nv,nr,slots,c=stats
                    fields=(f'skeleton={size} terminal={terminal} nonzero={pay} absent={absent} tagged={tagged} '
                            f'states={n} edges={ne} variable={nv} residuals={nr} peak_value_slots={slots}')
                rows=list(zip(*[iter(ticks)]*4))
                by_class=dict(zip(('overflow',*map(str,range(2,a.mother_cap+1))),map(list,rows)))
                if part in ref and ref[part]!=by_class or part in newest and newest[part]!=by_class:
                    raise ValueError(f'{engine}/{part}: exact conditional parity failed')
                newest[part]=by_class
                print(engine,part,fields,f'range={tuple(span)} checksum={c:016x} cpu={process_time()-t:.3f}',flush=True)
                for i,box in enumerate(base.endpoints(rows)):
                    for j,x in enumerate(box): total[i][j]+=x
                for i,x in enumerate((n,ne,nr)): counts[i]+=x
        finally:
            if handle: dll.release(handle)
        print(engine,'totals states/edges/residuals',*counts,flush=True)
        for name,(lo,hi) in zip(('lower','upper'),total):
            print(engine,name,lo,hi,f'[{float(lo):.12g}, {float(hi):.12g}]',flush=True)
    if a.entries:
        # Separate receipt: never overwrite even another experimental run.
        with a.entries.open('x') as f:
            json.dump(dict(tag,library_sha256=sha256(a.library.read_bytes()).hexdigest(),conditional=newest),f,indent=2)
            f.write('\n')


if __name__=='__main__':
    main()
