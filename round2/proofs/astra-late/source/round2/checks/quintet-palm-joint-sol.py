#!/usr/bin/env python3
"""Independent signed fixed-trial checker for the original quintet birth-Palm ratio."""
import argparse
from collections import defaultdict
from contextlib import closing
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import re
import sys

R = Path(__file__).resolve().parents[1]
CHECKS = R/'checks'
O = R/'occupation-quintet'
EXT = O/'edges/fy/extended/cut36-floor.log'
ROOTS = O/'conditioned-joint/roots.log'
EXTPIN = '56b31c406cfe1d794a85f6dea700211eaedbf297e5e1cd00470896e0bd28a824'
ROOTPIN = 'ef19be540999333210d6c503cbd9d54bbfd07a504c8654b59abe47fefc7b416e'
PREVIOUS = '23a20d846a206d16e451c871e17b0b1b26e595ecb4014a5921e0ec4b08dcffb3'
PRODUCER = '6d8899d1fcda135fc636c7723783e9dc9b5abf7015f6a1450099032c8f713f05'
STORED = {
    'lower': ('9223ca226bc472217b16b56318414bc94525603f50192230fe86aa4af0b043a6',Q('.565703539413800')),
    'upper': ('e93dcbb82544c69bc82c0e54722cd0c4641f64afc9d0aecba25c0b940e318eee',Q('.565746271349400')),
}


def read(path,pin):
    from hashlib import sha256
    raw = path.read_bytes()
    assert sha256(raw).hexdigest() == pin, path
    return raw


def load(path):
    spec = spec_from_file_location(path.stem,path)
    m = module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


def roots_and_extension():
    found = {}
    for line in read(ROOTS,ROOTPIN).decode().splitlines():
        side,bound,*_ = line.split()
        assert side in ('lower','upper') and bound == 'lower' and side not in found
        pairs = re.findall(r'Fraction\((\d+), (\d+)\)',line)
        v2,v3 = (Q(int(n),int(d)) for n,d in pairs)
        assert v2 > v3 > 0
        found[side] = (v2,v3)
    assert set(found) == {'lower','upper'}
    text = read(EXT,EXTPIN).decode()
    assert text.startswith('START lower bound ')
    parts = re.findall(r'^(?:mother pays|child pays with mother donor|child-child, both directed orders) (\S+)',text,re.M)
    [total] = re.findall(r'^TOTAL (\S+)',text,re.M)
    assert len(parts) == 3 and sum(map(Q,parts)) == Q(total) > 0
    return found, Q(total)*2**40


def objective(lp,cells,gold,roots,credits,trial,blend=Q(0),reduced=None):
    """In physical ticks, one cost per exact lifted column (including auxiliaries)."""
    assert 0 <= blend <= 1
    assert (reduced is not None) == (blend>0)
    v2,v3 = roots
    fixed = sum(credits)-gold['cap_f']-max(0,gold['cap_g'])
    out = []
    tail = 0
    for kind,coeff,m in lp.cols:
        if kind != 'occ':
            out.append(Q(0))
            continue
        d,k,b,p,child,mother = m
        assert dict(coeff).get(lp.rows['birth'],0) == b
        lo,hi,H = cells[d,k]
        source = b*(v2 if d == 2 else v3)
        if d == 25:
            tail += 1
            assert source == b*v3
        if blend:
            base,affine = reduced
            A,B = affine.get((d,k+1),(0,0)) if d <= 24 else (0,0)
            replacement = b*base+A+B*(child+mother)
            assert replacement >= 0
            source += blend*(replacement-source)
        out.append(source+fixed-H-trial*(hi if trial>=0 else lo))
    assert tail == 120
    return out


def repair(lp,cost,y,LO,uniform=False):
    """Exact dual bound for the signed objective in eta units; never clip at zero."""
    assert len(cost) == len(lp.cols) and set(y) == set(range(len(lp.rows)))
    assert lp.rhs == {'birth':1}
    defects = defaultdict(Q)
    residuals = []
    for (kind,entries,_), c in zip(lp.cols,cost,strict=True):
        r = sum((y[i]*co for i,co in entries),Q())-c
        residuals.append((kind,r))
        defects[kind] = max(defects[kind],r)
    assert set(defects) == {'occ','edge','stubs','tail-up','palm-slack'}
    assert all(r <= defects[kind] for kind,r in residuals)
    # eta masses: at most 1/LO for the first four, at most 1 for Palm slack.
    penalty = (sum(defects.values())/LO if uniform else
               defects['palm-slack']+sum(v for k,v in defects.items() if k!='palm-slack')/LO)
    return y[lp.rows['birth']]-penalty,dict(defects)


def proposal(lp,cost):
    from scipy.optimize import linprog
    scale,grid = 10**7,10**12
    A,rhs = lp.matrix()
    opts = dict(dual_feasibility_tolerance=1e-9,primal_feasibility_tolerance=1e-9)
    c = [float(z/scale) for z in cost]
    s = linprog(c,A_eq=A,b_eq=rhs,method='highs-ipm',options=opts)
    if not s.success:
        s = linprog(c,A_eq=A,b_eq=rhs,method='highs-ds')
    assert s.success,s.message
    return {i:Q(round(z*grid)*scale,grid) for i,z in enumerate(s.eqlin.marginals)}


def controls(lp,LO):
    # Real LP matrix: a signed constant along the birth row has exact eta value -3.
    birth = lp.rows['birth']
    cost = [-3*dict(col).get(birth,0) for _,col,_ in lp.cols]
    y = dict.fromkeys(range(len(lp.rows)),Q(0))
    y[birth] = Q(-3)
    assert repair(lp,cost,y,LO)[0] == -3
    palm = lp.rows['palm',0]
    y[palm] = Q(-2)
    _,defects = repair(lp,cost,y,LO)
    assert defects['palm-slack'] == 2


def perspective_control(lp,cells,gold,roots,credits,trial,prior):
    """Rebuild the homogeneous proposal in ratio units; certify separately."""
    from scipy.optimize import linprog
    from scipy.sparse import bmat, csc_matrix
    import numpy as np

    N = objective(lp,cells,gold,roots,credits,Q(0))
    x = objective(lp,cells,gold,roots,credits,trial)
    U = [Q(cells[m[:2]][1]) if kind=='occ' else Q(0) for kind,_,m in lp.cols]
    assert all(n-trial*u==cost for n,u,cost in zip(N,U,x,strict=True))
    assert all(u>0 if kind=='occ' else u==0 for u,(kind,_,_) in zip(U,lp.cols,strict=True))
    birth, age, mass = (lp.rows[r] for r in ('birth',('age',(2,0)),'mass'))
    virgins = 0
    for kind,terms,m in lp.cols:
        col = dict(terms)
        if kind!='occ':
            assert all(col.get(i,0)==0 for i in (birth,age,mass))
            continue
        d,k,b,p,u,v = m
        discount = b-col.get(mass,0)
        assert col.get(birth,0)==b and discount in prior
        if (d,k)==(2,0):
            assert b==Q(1,2) and u==0 and col[age]==Q(3,2)+v
            virgins += 1
        else:
            assert col[age]==-discount
    assert virgins==4
    matrix,rhs = lp.matrix()
    S = gold['verify']['S']
    eq = bmat([[matrix,csc_matrix(-np.asarray(rhs)[:,None])],
               [csc_matrix([[float(u/S) for u in U]]),csc_matrix((1,1))]],format='csc')
    # Use ratio units scaled by 1e6 directly, independently of producer scale.
    ans = linprog([float(n*10**6/S) for n in N]+[0],A_eq=eq,b_eq=[0]*len(rhs)+[1],method='highs-ipm')
    assert ans.success and ans.x[-1]>0,ans.message
    return ans.fun/10**6


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--side',choices=('lower','upper'),required=True)
    trials = parser.add_mutually_exclusive_group(required=True)
    trials.add_argument('--trial',type=Q,help='inward signed offset a from nominal')
    trials.add_argument('--endpoint',type=Q,help='translate public endpoint to a')
    trials.add_argument('--stored',action='store_true',help='verify the pinned producer vector with no replacement')
    parser.add_argument('--release-extension',action='store_true')
    parser.add_argument('--reduced-blend',type=Q,default=Q(0))
    parser.add_argument('--multipliers',type=Path,help='JSON exact rational strings in builder row insertion order, ticks')
    a = parser.parse_args()
    assert __debug__
    prev_path = CHECKS/'quintet-joint-sol.py'
    read(prev_path,PREVIOUS)
    prev = load(prev_path)
    for path in prev.PINS:
        prev.pinned(path)
    old = load(CHECKS/'quintet-stationary-sol.py')
    prior_gold = json.loads(prev.pinned(prev.WORK/'gold.json'))
    prior = tuple(Q(prior_gold[s]['certified_exact']) for s in prev.SIDES)
    K = prev.kernel(load(CHECKS/'stationary-p4-sol.py'),*prior)
    birth = load(R/'occupation-neon/pressure-memory/birth-law/replay.py')
    birth.LO,birth.HI,birth.SAT = *prior,7
    lp = birth.build(24,True,K)
    assert (len(lp.rows),len(lp.cols)) == (613,5554)
    controls(lp,prior[0])
    roots,ext = roots_and_extension()
    side = a.side
    i = prev.SIDES.index(side)
    direction = (1,-1)[i]
    gold = prior_gold[side]
    nominal = Q(gold['nominal'])
    if a.stored:
        assert not (a.release_extension or a.reduced_blend or a.multipliers)
        read(ROOTS.parent/'check.py',PRODUCER)
        pin,endpoint = STORED[side]
        receipt = json.loads(read(ROOTS.parent/f'{side}.json',pin))
        assert (receipt['side'],receipt['scale'],receipt['cost_scale'],receipt['dual_grid']) == (side,2**40,10**6,10**12)
        trial = Q(receipt['trial'])
        assert trial>0 and nominal+direction*trial==endpoint
    else:
        trial = a.trial if a.trial is not None else direction*(a.endpoint-nominal)
    assert trial is not None
    edge = prev.edge_inputs(i)
    extra = ext if side=='lower' and not a.release_extension else Q(0)
    sys.path.insert(0,str(prev.S))
    host = load(prev.S/'check.py')
    path = prev.WORK/f'quintet_{side}.bin'
    with closing(host.Comparison(path,side,gold)) as ref:
        env = old.Envelope(ref.snapshot,ref.raw,gold)
        cells = old.cells(env,ref)
        rr = None
        if a.reduced_blend:
            census = load(O/'census/census.py')
            with census.Snapshot(path,side,verify=False) as snap:
                rr = prev.reduced(side,gold,prior,snap)
                snap.check_unchanged()
        cost = objective(lp,cells,gold,roots[side],(*edge,extra),trial,a.reduced_blend,rr)
        if a.stored:
            integer = receipt['integer_dual']
            assert len(integer)==613 and all(type(z) is int for z in integer)
            y = {i:Q(z*receipt['cost_scale'],receipt['dual_grid']) for i,z in enumerate(integer)}
        elif a.multipliers:
            raw = json.loads(a.multipliers.read_bytes())
            assert isinstance(raw,list) and all(isinstance(z,(str,int)) for z in raw)
            y = {i:Q(z) for i,z in enumerate(raw)}
        else:
            y = proposal(lp,cost)
        result,defects = repair(lp,cost,y,prior[0])
        if a.stored:
            submitted, same_defects = repair(lp,cost,y,prior[0],uniform=True)
            assert result>=submitted>0 and same_defects==defects
            estimate = perspective_control(lp,cells,gold,roots[side],(*edge,extra),trial,prior)
            print(side,'same-law perspective proposal (diagnostic)',estimate,
                  'stored uniform-repair margin',submitted,flush=True)
        ref.check_identity()
    print(side,'trial',trial,'endpoint',nominal+direction*trial,'roots',roots[side],
          'edge floors',*edge,extra,'blend',a.reduced_blend,flush=True)
    print('613 rows / 5554 columns; signed eta lower',result,'approx',float(result),
          'defects',defects,flush=True)
    print('fixed-trial sign', 'certified' if result>=0 else 'not certified',flush=True)


if __name__ == '__main__':
    main()
