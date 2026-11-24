#!/usr/bin/env python3
"""Exact finite raw residual controls for four signed D18/D20 potentials."""
import argparse
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import sys

import numpy as np
from scipy.optimize import linprog

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'rows'))
from extract import Cert, Rows, audit, children, embeddings, stamp  # noqa: E402

S = Q('.56572365')
NB = (1, 0, 1)
# Original global comparison bounds, rounded outward (not stationary fits).
CAP = tuple(map(Q, ('1.84124847315', '1.95820675182', '1.84219674916', '2.00854054323')))
FLOOR = (Q(0), Q(0), Q(113472636536, 4*10**11), Q(115454594604, 4*10**11))
FIXED = {('all', False): (533,535,537,538), ('all', True): (320,532,533,535,538),
         ('chain', True): (70,80,216,286,301), ('branch', True): (381,442,467,469,509)}


def add(x, i, y):
    return x[:i] + (y,) + x[i+1:]


def raw(c, rows, p, leaves):
    """Compressed physical generator, re-evaluating the whole potential."""
    ch = children(p)
    x = tuple((l,len(k),int(i>0)) for i,(l,k) in enumerate(zip(leaves,ch)))
    d = tuple(map(sum,x))
    assert min(d)>=2
    def value(t,y):
        return sum(rows.value(f,tuple(y[j] for j in ids)) for f,ids in embeddings(t))
    w = value(p,x)
    generator = Q(0)
    for i,(l,k,h) in enumerate(x):
        ns = ch[i] + ([p[i]] if h else [])
        g = l + sum(Q(1,d[j]) for j in ns)
        generator += g * (value(p,add(x,i,(l+1,k,h)))-w)
        if l:
            generator += Q(l,d[i])*(value(p+(i,),add(x,i,(l-1,k+1,h))+(NB,))-w)
    return w/c.S, (generator-S*w)/c.S


def independent(c,p,leaves):
    """Every ordinary directed arrow, and the separate whole-tree evaluator."""
    tree = audit.build(p,leaves)
    def val(par):
        return c.potential(*audit.physical_view(par)[:2])/c.S
    w = val(tree)
    rates = audit.arrows(tree)
    assert sum(rates.values()) == len(tree)
    return w, sum(r*(val(tree+(j,))-w) for j,r in rates.items())-S*w


def catalog():
    """Deterministic complete trees; zero stocks are legal at interior cores."""
    for n in (2,3,4,6,9):
        chain = (-1, *range(n-1))
        star = (-1, *(0,)*(n-1))
        for label,p in (('chain',chain),('star',star)):
            for l in (1,2,3,5,8,14,16,17,18,19,20,28,64,256,16000):
                yield label+'-hom',p,(l,)*n
            for a,b in ((0,1),(0,5),(1,2),(1,5),(1,18),(1,20),(1,64),(2,18),(17,19),(18,20)):
                degree = [p.count(i)+(i>0) for i in range(n)]
                for phase in range(2):
                    ls = tuple(max(2-d,(a,b)[(i+phase)%2]) for i,d in enumerate(degree))
                    yield label+'-het',p,ls
    for p in ((-1,0,0,1),(-1,0,1,1),(-1,0,0,1,1,2,2)):
        for seq in product((0,2,17,20), repeat=3):
            ls = tuple(max(2-p.count(i)-(i>0),seq[i%3]) for i in range(len(p)))
            yield 'branch',p,ls
    for n in (18,19,20,21,22):
        p = (-1,*(0,)*(n-1))
        yield 'star-arity-seam',p,(1,)*n
        yield 'star-arity-seam',p,(5,)+tuple(1+2*(i%2) for i in range(n-1))


def solve_exact(a,b):
    m = [list(row)+[y] for row,y in zip(a,b)]
    n = len(m)
    for i in range(n):
        j = next((j for j in range(i,n) if m[j][i]),None)
        if j is None:
            return None
        m[i],m[j]=m[j],m[i]
        k=m[i][i]
        m[i]=[v/k for v in m[i]]
        for j in range(n):
            if i!=j and m[j][i]:
                k=m[j][i]
                m[j]=[u-k*v for u,v in zip(m[j],m[i])]
    return [r[-1] for r in m]


def cert(rows, target, budget, contrast, fixed=None):
    # Express the target as a linear combination of actual residual rows.
    # Each equality is checked in exact Q after HiGHS chooses a basis.
    vectors = [r+(Q(1),) if contrast else r for r in rows]
    ids=fixed
    if ids is None:
        n = len(vectors)
        a = np.array(vectors,dtype=float).T
        z = np.max(np.abs(a),axis=1)
        out=linprog(np.ones(2*n),A_eq=np.hstack((a,-a))/z[:,None],
                    b_eq=np.array(target,dtype=float)/z,bounds=(0,None),method='highs')
        assert out.success, out.message
        ids=sorted(set(i%n for i in np.flatnonzero(out.x>1e-8)))
    assert len(ids)==len(target), ids
    coeff=solve_exact(list(zip(*(vectors[i] for i in ids))),target)
    assert coeff is not None
    assert tuple(sum(w*vectors[i][j] for i,w in zip(ids,coeff)) for j in range(len(target)))==target
    norm=sum(map(abs,coeff))
    assert norm*budget<1
    return norm,ids,coeff


def diagnostic(rows, columns):
    """Common-offset fit in a reference gauge, not used for the proof."""
    a=np.array(rows,dtype=float)[:,columns]
    m=len(columns)
    b=np.column_stack((a,np.ones(len(a))))
    aeq=[[1]*m+[0,0]]
    out=linprog([0]*(m+1)+[1],A_ub=np.column_stack((np.vstack((b,-b)),-np.ones(2*len(a)))),
                b_ub=np.zeros(2*len(a)),A_eq=aeq,b_eq=[1],bounds=[(None,None)]*(m+1)+[(0,None)],method='highs')
    assert out.success, out.message
    return out.fun, tuple(out.x[:m])


def main():
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--fixed',action='store_true',help='replay selected supports, no numerical optimization')
    cli.add_argument('--blobs',nargs=4,type=Path,default=list(map(Path,('b3f18_lo.bin','b3f18_hi.bin','b3f20_lo.bin','b3f20_hi.bin'))))
    a=cli.parse_args()
    pins=list(map(stamp,a.blobs))
    cs=[*(audit.Cert(p) for p in a.blobs[:2]),*(Cert(p,s,verify=False) for p,s in zip(a.blobs[2:],('lower','upper')))]
    assert [(c.D,c.K,c.mode,c.lam,c.S) for c in cs]==[(18,17,2,Q('.565673'),10**11),(18,17,1,Q('.565775'),10**11),(20,19,2,Q('.565680'),10**11),(20,19,1,Q('.565768'),10**11)]
    assert all(p[2]==3495091968 for p in pins[:2])
    es=list(map(Rows,cs))
    fixtures=list(dict.fromkeys(catalog()))
    computed=set().union(*FIXED.values()) if a.fixed else set(range(len(fixtures)))
    raw_rows={i:tuple(raw(c,e,fixtures[i][1],fixtures[i][2])[1]/len(fixtures[i][1]) for c,e in zip(cs,es))
              for i in sorted(computed)}
    print(f'trial {S}, {len(raw_rows)} physical trees, 4 raw columns',flush=True)
    chosen=set()
    for group,contrast in (('all',False),('all',True),('chain',True),('branch',True)):
        index=[i for i,(name,p,ls) in enumerate(fixtures) if group=='all' or
               (name.startswith(group) and max(ls)<=64)]
        if a.fixed:
            assert set(FIXED[group,contrast]) <= set(index)
            index=list(FIXED[group,contrast])
        result=[]
        for target in product(*zip(FLOOR,CAP)):
            t=target+(Q(0),) if contrast else target
            norm,ix,coeff=cert([raw_rows[i] for i in index],t,Q('1e-8'),contrast,
                               list(range(len(index))) if a.fixed else None)
            ids=[index[j] for j in ix]
            result.append((norm,t,ids,coeff)); chosen.update(ids)
        norm,corner,ids,coeff=max(result)
        print(group,'contrast' if contrast else 'absolute','exact box corners',len(result),
              'worst l1',float(norm),'uniform obstruction',float(1/norm),flush=True)
        print('worst corner',*[str(x) for x in corner],'rows',ids,flush=True)
        assert 1/norm>Q('1e-8')
        floor=(10**12/norm)//1
        assert Q(floor,10**12)*norm<=1
        print(' outward lower (units 10^-12)',floor,flush=True)
    for group in (() if a.fixed else ('all','chain','branch')):
        selected=[raw_rows[i] for i,(name,_,ls) in enumerate(fixtures) if group=='all' or
                  (name.startswith(group) and max(ls)<=64)]
        for columns in ((0,1),(2,3),(0,1,2,3)):
            obj,coef=diagnostic(selected,columns)
            print('sum-one diagnostic',group,columns,f'{obj:.10g}',tuple(round(float(x),8) for x in coef),flush=True)
    # Add small mixed-regime trees to check all motif families, not only LP rows.
    chosen.update(i for i,(t,p,ls) in enumerate(fixtures)
                  if t=='branch' and len(p)==7 and sum(ls)<26)
    for i in sorted(chosen):
        _,p,ls=fixtures[i]
        for j,(c,e) in enumerate(zip(cs,es)):
            out=independent(c,p,ls)
            assert raw(c,e,p,ls)==out and (i not in raw_rows or raw_rows[i][j]==out[1]/len(p)),(i,j)
    assert pins==list(map(stamp,a.blobs))
    print('independent ordinary-arrow whole-tree matches',len(chosen)*4,flush=True)


if __name__=='__main__':
    main()
