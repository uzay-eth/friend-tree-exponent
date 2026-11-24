"""Project sampled physical trees onto every LP row, with cluster diagnostics."""
import argparse
from collections import defaultdict
from functools import lru_cache
from fractions import Fraction as F
from itertools import combinations, product

import numpy as np
from scipy.sparse import csr_matrix

from lp_level1 import Model
from lp_level2 import MotifModel


def audit(model, data, mu, batch):
    N, T = model.N, model.D3
    keys = {name:j for j,name in enumerate(model.names)}

    def bucket(d,c):
        return (d,c) if d<=N else (0,min(c,N+1))

    @lru_cache(maxsize=100000)
    def lift(p,c):
        def interp(s):
            d,cap = s
            if d:
                return [(F(1,d),0)]
            u = F(1,max(N+1,cap))
            return [(F(0),0),(u,1)]
        return [(keys['e',p,c,x,y],ix,iy,x,y)
                for (x,ix),(y,iy) in product(interp(p),interp(c))]

    @lru_cache(maxsize=100000)
    def hooks(p,c):
        return tuple(model.ext.get((mode,(p,c))) for mode in ('ancestor','sibling','child'))

    rows = list(model.exact_rows(mu,mu))
    ridx, cols, vals = [], [], []
    for i, (_, expr, _) in enumerate(rows):
        for j,value in expr.items():
            ridx.append(i); cols.append(j); vals.append(float(value))
    A = csr_matrix((vals,(ridx,cols)), shape=(len(rows),len(keys)))
    rhs = np.array([r for _,_,r in rows])
    samples, counts = [], []
    vec = np.zeros(len(keys))
    offsets, arr = data['offsets'], data['cores']
    size = 0
    for t,(start,end) in enumerate(zip(offsets[:-1],offsets[1:]),1):
        members = arr[start:end].tolist()
        sts = [bucket(d,c) for c,d,_ in members]
        xs = [1/d for _,d,_ in members]
        kids = [[] for _ in members]
        # The finite initial seed has no parent. Its absent parent credit is 0.
        ks = [0.]*len(members)
        for i,(_,_,p) in enumerate(members[1:],1):
            kids[p].append(i)
            ks[p] += xs[i]
        rho = [k+(xs[p] if p>=0 else 0.) for k,(_,_,p) in zip(ks,members)]
        for s in sts:
            vec[keys['p',s]] += 1
        for i,(_,d,p) in enumerate(members[1:],1):
            _,dp,g = members[p]
            s,sp = sts[i],sts[p]
            # Seed has one fewer parent; a seed with c=1 has no legal LP edge.
            if sp[1]==1:
                continue
            for j,ip,ic,x,y in lift(sp,s):
                pp = xs[p]*max(N+1,sp[1]) if not sp[0] else 0.
                pc = xs[i]*max(N+1,s[1]) if not s[0] else 0.
                vec[j] += (pp if ip else 1-pp)*(pc if ic else 1-pc)
            ja,js,jc = hooks(sp,s)
            if ja is not None and g>=0: vec[ja] += xs[g]
            if js is not None: vec[js] += ks[p]-xs[i]
            if jc is not None: vec[jc] += ks[i]
            if not T or d>T or dp>T:
                continue
            if g>=0 and members[g][1]<=T:
                st = sts[g],sp,s
                j = keys.get(('chain',*st))
                if j is not None:
                    vec[j] += 1
                    for k, x in enumerate((rho[g]-xs[p],rho[p]-xs[g]-xs[i],rho[i]-xs[p])):
                        j = keys.get(('motif-external','C',st,k))
                        if j is not None: vec[j] += x
            elif g>=0:
                vec[keys['omitted',('ancestor',(sp,s))]] += 1
            for z in kids[i]:
                if members[z][1]>T:
                    vec[keys['omitted',('child',(sp,s))]] += 1
            for z in kids[p]:
                if z!=i and members[z][1]>T and js is not None:
                    vec[keys['omitted',('sibling',(sp,s))]] += 1
        if T:
            for p, children in enumerate(kids):
                if members[p][1]>T or members[p][0]<3:
                    continue
                for i,z in combinations((i for i in children if members[i][1]<=T),2):
                    if sts[i]>sts[z]: i,z=z,i
                    st=sts[p],sts[i],sts[z]
                    vec[keys['wedge',*st]] += 1
                    for k,x in enumerate((rho[p]-xs[i]-xs[z],rho[i]-xs[p],rho[z]-xs[p])):
                        j=keys.get(('motif-external','V',st,k))
                        if j is not None: vec[j] += x
        size += len(members)
        if t%batch==0 or t==len(offsets)-1:
            samples.append(A@vec-size*rhs)
            counts.append(size)
            size = 0
            vec.fill(0)
    R, n = np.array(samples), np.array(counts)
    mean = R.sum(axis=0)/n.sum()
    se = (R-n[:,None]*mean).std(axis=0,ddof=1)/(len(n)**.5*n.mean())
    families = defaultdict(list)
    for i,(name,_,_) in enumerate(rows):
        family = name if isinstance(name,str) else name[0]
        if family in ('station-high','station-low','+','-'):
            family2 = name[1] if isinstance(name[1],str) else name[1][0]
            family += '/' + family2
        families[family].append(i)
    print(f'N={N} D3={T} rows={len(rows)} batches={len(n)} cores={n.sum()}')
    for family, indices in sorted(families.items()):
        i = max(indices,key=lambda j:mean[j]); excess = mean[indices]>5*se[indices]+1e-11
        print(f'{family:33} n={len(indices):6d} max={mean[i]:+.4g} SE={se[i]:.3g} '
              f'over5SE={sum(excess)}  {rows[i][0]}')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('sample')
    parser.add_argument('--configs', nargs='+', default=['8:0','8:6'])
    parser.add_argument('--batch',type=int,default=50)
    parser.add_argument('--mu',type=F,default=F('.56572366'))
    args=parser.parse_args()
    for cfg in args.configs:
        n,d = map(int,cfg.split(':'))
        audit((MotifModel if d else Model)(n,d),np.load(args.sample),args.mu,args.batch)
