#!/usr/bin/env python3
"""Exact original-GOLD1 composition replay; no optimizer or producer compiler."""
import argparse
from contextlib import closing
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import re
import sys

R = Path(__file__).resolve().parents[1]
CHECKS = R/'checks'
O = R/'occupation-quintet'
BASE = O/'conditioned-joint'
COMPOSITION = BASE/'composition'
DONOR = O/'edges/old-donor'
S = 2**40
SIDES = ('lower','upper')
PINS = {
    CHECKS/'quintet-palm-joint-sol.py':'af2ebf9e53cf8bf50b57e7274a06c80487fe1354b1ef98bc158965007885046b',
    CHECKS/'quintet-h3-cohort-sol.py':'441a102eeef9ab9d67db40777ae047b42cc5a00486db80d8d0472035af7fdf4d',
    DONOR/'matched.log':'35b7573b85a005166032ff3596282d5ea3cbb07c4d77d1e38e2bbc8932c65e0d',
    DONOR/'lower-wide.log':'b81568c75a4e57b94db733acc4d9deaac2a4eb568b7c498091b6654316d24ee3',
}
OLD_MARGIN = (
    Q('246796378766836906157775417183400184888681/328929067814462075504886439151250000000000'),
    Q('337626417760847229096216385970226521841377/603831138983070476772255298925000000000000'),
)
RECEIPT_PINS = (
    '92d2e19fd91cd165f7ab8312406f3e3b6f00b061689bd81bb56579c59c78b5ed',
    'ca1c0184fb092767919b33332767ad0b47eb35181ca21616f441421334eb8bde',
)


def pinned(path, digest):
    raw = path.read_bytes()
    assert sha256(raw).hexdigest() == digest, path
    return raw


def load(path):
    spec = spec_from_file_location(path.stem, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def donor_floor(path, side):
    text = pinned(path, PINS[path]).decode()
    sections = re.split(r'^(lower|upper) (floor|ceiling) ', text, flags=re.M)
    [floor] = [body for s, k, body in zip(sections[1::3], sections[2::3], sections[3::3], strict=True)
               if (s,k) == (side,'floor')]
    values = re.findall(r'^(?:mother pays|mother donates|two children, both orders) (\S+)', floor, re.M)
    [total] = re.findall(r'^TOTAL (\S+)', floor, re.M)
    assert len(values) == 3 and sum(map(Q,values)) == Q(total) > 0
    return Q(total)


def stationary_upper(prev, side):
    text = prev.pinned(prev.S/'host24.log').decode()
    [(lo,hi)] = re.findall(rf'^{side} stationary denominator lower/upper in ticks (\S+) (\S+)', text, re.M)
    joint = prev.pinned(prev.S/'joint24.log').decode()
    [joint_hi] = re.findall(rf'^{side} joint H savings.* stationary denominator in ticks \S+ (\S+)$',joint,re.M)
    assert 0 < Q(lo) < Q(hi) <= Q(joint_hi)
    return Q(hi),Q(joint_hi)


def multipliers(receipt, side):
    assert receipt['side'] == side
    if set(receipt) != {'side','trial','integer_dual'}:
        assert (receipt['scale'],receipt['cost_scale'],receipt['dual_grid']) == (S,10**6,10**12)
    integer = receipt['integer_dual']
    assert len(integer) == 613 and all(type(n) is int for n in integer)
    trial = Q(receipt['trial'])
    assert trial > 0
    return trial, {i:Q(n*10**6,10**12) for i,n in enumerate(integer)}


def cost_at(N, U, trial, credit, occurrence):
    return [n+a*credit-trial*u for n,u,a in zip(N,U,occurrence,strict=True)]


def signed_check(lp, cost, dual, LO):
    assert (len(lp.rows), len(lp.cols), len(dual), len(cost)) == (613,5554,613,5554)
    assert lp.rhs == {'birth':1} and set(dual) == set(lp.rows.values())
    defect = dict.fromkeys(('occ','edge','stubs','tail-up','palm-slack'), Q(0))
    for (kind,coeff,_), c in zip(lp.cols,cost,strict=True):
        defect[kind] = max(defect[kind],sum(dual[i]*a for i,a in coeff)-c)
    uniform = dual[lp.rows['birth']]-sum(defect.values())/LO
    better = uniform+(1/LO-1)*defect['palm-slack']
    return uniform, better, defect


def replay(side,lp,N,U,occ,booked,roots,cells,data,palm,LO,nominal,direction,
           algebraic_trial,added,donor_ticks,h3_ticks,host,receipt):
    trial,dual = multipliers(receipt,side)
    assert trial > algebraic_trial
    c = cost_at(N,U,trial,added,occ)
    assert c == palm.objective(lp,cells,data,roots,(*booked,donor_ticks,h3_ticks),trial)
    base = cost_at(N,U,trial,0,occ)
    assert all(x-y == a*added for x,y,a in zip(c,base,occ,strict=True))
    uniform,tight,defects = signed_check(lp,c,dual,LO)
    independent,other_defects = palm.repair(lp,c,dual,LO)
    assert tight == independent and defects == other_defects
    assert tight >= uniform > 0,(side,uniform,defects)
    endpoint = nominal+direction*trial
    print(side,'613 multipliers / 5554 columns; trial',trial,'endpoint',endpoint,
          'outward',host.outward(endpoint,side),flush=True)
    print(side,'exact uniform eta margin',uniform,'approx',float(uniform),
          'Palm-specific adjustment',tight-uniform,'defects',defects,flush=True)


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--side', choices=SIDES, action='append')
    cli.add_argument('--controls-only', action='store_true')
    args = cli.parse_args()
    assert __debug__
    for path,digest in PINS.items():
        pinned(path,digest)
    palm = load(CHECKS/'quintet-palm-joint-sol.py')
    h3 = load(CHECKS/'quintet-h3-cohort-sol.py')
    pinned(CHECKS/'quintet-joint-sol.py',palm.PREVIOUS)
    prev = load(CHECKS/'quintet-joint-sol.py')
    for path in prev.PINS:
        prev.pinned(path)
    gold = json.loads(prev.pinned(prev.WORK/'gold.json'))
    LO,HI = (Q(gold[s]['certified_exact']) for s in SIDES)
    K = prev.kernel(load(CHECKS/'stationary-p4-sol.py'),LO,HI)
    birth = load(R/'occupation-neon/pressure-memory/birth-law/replay.py')
    birth.LO,birth.HI,birth.SAT = LO,HI,7
    lp = birth.build(24,True,K)
    assert (len(lp.rows),len(lp.cols)) == (613,5554)
    palm.controls(lp,LO)
    roots,ext = palm.roots_and_extension()
    sys.path.insert(0,str(prev.S))
    host = load(prev.S/'check.py')
    stationary = load(CHECKS/'quintet-stationary-sol.py')
    for side in args.side or SIDES:
        i = SIDES.index(side)
        direction = (1,-1)[i]
        data = gold[side]
        nominal = Q(data['nominal'])
        assert data['verify']['S'] == S and data['verify']['mode'] == (2,1)[i]
        donor = donor_floor(DONOR/('lower-wide.log' if i==0 else 'matched.log'),side)
        *components, h3_floor = map(Q,h3.EXPECTED[i])
        assert sum(components) == h3_floor > 0
        donor_ticks,h3_ticks = donor*S,h3_floor*S
        added = donor_ticks+h3_ticks
        old_pin, old_endpoint = palm.STORED[side]
        old_trial,old_y = multipliers(json.loads(pinned(BASE/f'{side}.json',old_pin)),side)
        assert nominal+direction*old_trial == old_endpoint
        Mhi,joint_hi = stationary_upper(prev,side)
        with closing(host.Comparison(prev.WORK/f'quintet_{side}.bin',side,data)) as ref:
            cells = stationary.cells(stationary.Envelope(ref.snapshot,ref.raw,data),ref)
            booked = (*prev.edge_inputs(i),ext if i==0 else Q(0))
            N = palm.objective(lp,cells,data,roots[side],booked,Q(0))
            occ = [int(kind=='occ') for kind,_,_ in lp.cols]
            U = [Q(cells[m[:2]][1]) if a else Q(0) for a,(_,_,m) in zip(occ,lp.cols,strict=True)]
            assert all(u>0 if a else u==0 for u,a in zip(U,occ,strict=True))
            old_cost = cost_at(N,U,old_trial,0,occ)
            assert old_cost == palm.objective(lp,cells,data,roots[side],booked,old_trial)
            old_uniform,old_exact,_ = signed_check(lp,old_cost,old_y,LO)
            assert old_uniform == old_exact == OLD_MARGIN[i] > 0
            print(side,'old eta margin',old_uniform,'stationary M upper ticks',Mhi,flush=True)
            print(side,'normalized donor / H3',donor,h3_floor,'ticks',donor_ticks,h3_ticks,flush=True)
            for label,credit in (('donor-only',donor_ticks),('H3-only',h3_ticks),('combined',added)):
                t = old_trial+credit/Mhi
                assert (t-old_trial)*Mhi == credit > 0
                print(side,'algebraic',label,'trial',t,'outward',host.outward(nominal+direction*t,side),flush=True)
                jt = old_trial+credit/joint_hi
                assert (jt-old_trial)*joint_hi == credit and jt<=t
                print(side,'algebraic joint24 ceiling',label,'outward',host.outward(nominal+direction*jt,side),flush=True)
            if not args.controls_only:
                receipt = json.loads(pinned(COMPOSITION/f'{side}.json',RECEIPT_PINS[i]))
                replay(side,lp,N,U,occ,booked,roots[side],cells,data,palm,LO,nominal,direction,
                       old_trial+added/Mhi,added,donor_ticks,h3_ticks,host,receipt)
            ref.check_identity()


if __name__ == '__main__':
    main()
