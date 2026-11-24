"""Compare identical ChainModel supports before/after adding exact slab rows."""
import argparse
from fractions import Fraction as Q
from pathlib import Path
from time import process_time

from bands import BandModel
from price_chains import high_charge
from solve import optimize


def main():
    a=argparse.ArgumentParser(description=__doc__)
    a.add_argument('--threshold',type=int,default=3)
    a.add_argument('--arity',type=int,default=1)
    a.add_argument('--order',type=int,default=1)
    a.add_argument('--inner',type=int,default=2)
    a.add_argument('--role',type=int,choices=range(3),default=0)
    a.add_argument('--chain-arity',type=int,nargs=3,default=(1,1,0))
    a.add_argument('--payloads',type=Path,nargs=2)
    a.add_argument('--mode',nargs='+',choices=('mass','reciprocal','difference'),default=['mass','reciprocal'])
    args=a.parse_args()
    if 'difference' in args.mode and (not args.payloads or args.order<3):
        a.error('difference needs both payloads and order >=3')
    m=BandModel(args.threshold,args.arity,args.order,inner=args.inner,wide_roles=(args.role,),
                chain_arity=args.chain_arity,local=False)
    costs={}
    for mode in args.mode:
        if mode=='difference':
            paid,cost=high_charge(m,args.payloads,args.role,mode,Q(1),18)
            print('authenticated original hV support cells',paid,flush=True)
        else:
            cost=m.chain_objective(high_role=args.role,high_slabs=(0,),reciprocal=mode=='reciprocal')
        costs[mode]=cost
    print('cut',m.boundaries,'order',m.m,'axis',args.role,'chain cells',len(m.chains),flush=True)
    boxes={}
    for name in ('global','slab'):
        t=process_time()
        if name=='slab': m.add_local()
        print(name,'columns',len(m.variables),'equalities',len(m.equalities),'inequalities',len(m.upper_rows),
              'seam occurrences',len(m.seams),'empty cells',len(m.empty),f'build_cpu={process_time()-t:.3f}',flush=True)
        for mode,cost in costs.items():
            out=[]
            for maximize in (False,True):
                t=process_time()
                v,err=optimize(m,cost,maximize,scale=10**7 if mode=='difference' else 10,method='highs-ipm')
                val=(max if maximize else min)(v*s for s in (m.low,m.high))
                print(name,mode,'upper' if maximize else 'lower',val,f'{float(val):.12g}',
                      f'repair={float(err):.4g} cpu={process_time()-t:.3f}',flush=True)
                out.append(val)
            if name=='slab':
                lo,hi=boxes[mode]
                out=max(lo,out[0]),min(hi,out[1])
                print('intersection',mode,*(f'{float(x):.12g}' for x in out),
                      'width reduction',f'{float(hi-lo-out[1]+out[0]):.12g}',flush=True)
            else: boxes[mode]=out


if __name__=='__main__':
    main()
