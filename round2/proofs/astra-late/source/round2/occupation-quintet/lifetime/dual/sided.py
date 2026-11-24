#!/usr/bin/env python3
"""Price the masked reduced reward in the existing physical birth-law LP."""
import argparse
import hashlib
from fractions import Fraction as Q
from json import loads

from export import HERE, COMMIT, RECEIPT_HASH, EXPECTED, load

EXPORT_HASH = {
    'lower':'d698df2bdd20fd6d5a802f84d18124a1fcd6f491bb4dde8ae5038ca5d27f6169',
    'upper':'a23fe15ac518c35f9e7dcc92bc1fa497253f16e7677a9d94e572815682253fa8',
}
OBJECTIVE_SCALE = 10**7


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--limit',type=int,default=24)
    a=p.parse_args()
    assert 19 <= a.limit <= 64
    r2=HERE.parents[2]
    lpmod=load('sided_birth',r2/'occupation-neon/pressure-memory/birth-law/replay.py')
    audit=load('birth_kernel',r2/'checks/stationary-p4-sol.py')
    pair=[]
    for side,h in EXPORT_HASH.items():
        raw=(HERE/f'{side}.json').read_bytes()
        assert hashlib.sha256(raw).hexdigest()==h
        pair.append(loads(raw))
    lo,hi=Q(pair[0]['LO']), Q(pair[0]['HI'])
    assert all((Q(z['LO']),Q(z['HI']))==(lo,hi) for z in pair)
    # Kernel coefficients are only physical birth-class lower probabilities.
    # They remain valid at any smaller upper discount, regardless of witness.
    assert hi <= audit.HI
    kernel=audit.class_kernel()
    lpmod.LO,lpmod.HI=lo,hi
    lp=lpmod.build(a.limit,True,kernel)
    print('limit',a.limit,'rows',len(lp.rows),'columns',len(lp.cols),flush=True)
    for z in pair:
        side=z['side']; grid=z['V_grid']
        assert z['source_commit']==COMMIT and z['receipt_sha256']==RECEIPT_HASH
        assert Q(z['V_J_ticks'])*grid==EXPECTED[side]
        data={(d,c):(Q(A),Q(B)) for d,c,A,B in z['rows']}

        def objective(kind,meta):
            if kind!='occ': return Q(0)
            d,k,b,p,u,v=meta
            if d>a.limit: return Q(0)
            A,B=data.get((d,k+1),(0,0))
            rho=u+v  # exactly the same two coordinates as age balance and edge marginals
            return A+B*rho

        assert all(objective(k,m)>=0 for k,_,m in lp.cols)
        lb=max(Q(0),lp.lower(lambda kind,m: objective(kind,m)/OBJECTIVE_SCALE)*OBJECTIVE_SCALE)
        stat=lo*lb
        revised=lo*Q(z['V_J_ticks'])+stat
        print(side,'exact reduced stationary floor',stat,'floor ticks',lpmod.dec(stat),
              'full f floor ticks',lpmod.dec(revised), 'exact full', revised,flush=True)


if __name__=='__main__':
    main()
