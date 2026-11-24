"""Freeze top observed paired physical rows, separately by family and cohort."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import sys

from probe import physical_view, snapshots

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / 'higher/older'))
from replay import Cert, Local, SPEC


def queries():
    lo,hi = [Local(Cert(Path(f'b3f18_{s}.bin'))) for s in ('lo','hi')]
    by_cohort=[]
    n=16384
    for first in (20261124,20261201):
        best={}
        for seed in range(first,first+8):
            par=next(snapshots(seed,(n,)))
            ch,st,core=physical_view(par)
            up=[-1]*len(st)
            for j,kids in enumerate(ch):
                for i in kids: up[i]=j
            rec=[Q(1,sum(x)) for x in st]
            fp=sha256(json.dumps(par,separators=(',',':')).encode()).hexdigest()
            for f,ids in lo.occurrences(ch,st):
                if f in ('f','g') or not ids[0]:continue
                x=tuple(st[j] for j in ids)
                t=SPEC[f]
                exterior=tuple((rec[up[j]] if i==0 else Q(),
                                sum((rec[a] for a in ch[j] if a not in ids),Q()))
                               for i,j in enumerate(ids))
                assert all(st[j][1]-t.count(i)==sum(a not in ids for a in ch[j]) for i,j in enumerate(ids))
                a,b=(e.charge(f,x,exterior) for e in (lo,hi))
                assert min(a,b)>=0
                mode=('inner' if max(map(sum,x))<18 else 'mixed-tail')
                for side,v in (('lo',a),('hi',b)):
                    key=(f,mode,side)
                    if not v or (key in best and best[key][0]>=v):continue
                    best[key]=(v, {'family':f,'states':x,'exterior':[list(map(str,r)) for r in exterior]},
                               {'seed':seed,'n':n,'cohort_first_seed':first,'host_labels':[core[j]+1 for j in ids],
                                'parents_sha256':fp,'lo':str(a),'hi':str(b),'degree_class':mode})
        by_cohort.append(best)
    seen=set()
    for best in by_cohort:
        for key,(_,q,p) in sorted(best.items()):
            token=json.dumps(q,separators=(',',':'))
            if token in seen:continue
            seen.add(token)
            yield token,p


if __name__=='__main__':
    for q,p in queries():
        print(q)
        print(json.dumps(p,separators=(',',':')),file=sys.stderr)
