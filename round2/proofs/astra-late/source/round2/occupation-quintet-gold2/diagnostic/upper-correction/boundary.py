#!/usr/bin/env python3
"""Check physical signed slopes and the rational Bellman drift at mask boundaries."""
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
import json
from math import lcm
from pathlib import Path

HERE = Path(__file__).resolve().parent
GOLD = HERE.parents[1]


def load(path, pin=None):
    data = path.read_bytes()
    assert pin is None or sha256(data).hexdigest() == pin
    m = module_from_spec(spec_from_file_location(path.stem, path))
    exec(compile(data, str(path), 'exec'), m.__dict__)
    return m


def check(r, solver, disc, adj, p, cohort):
    A, slopes = r(p)
    caps = [Q(c-len(ns),2) for (_,c),ns in zip(p,adj,strict=True)]
    if A+sum(min(0,b*w) for b,w in zip(slopes,caps,strict=True)) >= 0:
        return False
    # Zero successor values make this an exit-row affine recurrence; its
    # quotient is checked independently at all exact physical corners.
    fixed = [n-c+sum(Q(1,p[j][0]) for j in ns) for (n,c),ns in zip(p,adj,strict=True)]
    births = [Q(n-c,n) for n,c in p]
    scale=solver.GRID
    common=lcm(2,*(n for n,c in p))
    widths=[int(common*w) for w in caps]
    bad=[]; good=[]
    for b,w in zip(slopes,widths,strict=True):
        if cohort:
            bad.append(-(-scale*b//1)); good.append(-scale*b//1)
        else:
            bad.append(Q(-scale*b,common)//1); good.append(-scale*b//1)
    for b,g,w in zip(slopes,good,widths,strict=True):
        assert not w or g<=-scale*b<g+1
    rates=list(product(*((Q(0),w) if w else (Q(0),) for w in caps)))
    def exact(z, clip):
        return min(scale*(max(0,v) if clip else v)-z*(disc+sum(fixed)+sum(births)+sum(xs))
                   for xs in rates for v in [-A-sum(b*x for b,x in zip(slopes,xs,strict=True))])
    def outward(offsets):
        terms=[]
        for f,w,b,z in zip(fixed,widths,births,offsets,strict=True):
            terms.extend(((int(common*f),w,z,0),(int(common*b),0,0,0)))
        if cohort:
            den=lcm(common,disc.denominator)
            return solver.bellman(-A*scale*den//1,int(den*(disc+sum(fixed)+sum(births))),
                                  [(z,w*(den//common)) for z,w in zip(offsets,widths,strict=True)],allow_negative=True)
        return solver.affine(-A,terms,common,disc,False)
    old,new=map(lambda a:max(0,outward(a)),(bad,good))
    assert exact(new,True)>=0
    gap=exact(old,True)
    if gap>=0:
        return False
    print('physical boundary',p,'d',common,'outward A',A,'outward slopes',slopes,'rate caps',caps,flush=True)
    print('old/corrected integer values',old,new,'true minimum residual old/corrected',gap,exact(new,True),flush=True)
    assert old>new
    return True


def main():
    b=load(HERE/'replay-threshold.py').load()
    cohort=load(GOLD/'higher/cohort.py','b82b6768a91ee1ae66f319b078b9076d04878546705a3538795fff73ec06057c')
    private=b.gold.PRIVATE
    wire=load(GOLD.parent/'attacks/neon2-quintet/gold2/physical/check.py')
    data=(private/'source/round2/neon2/quintet-cert/certificate_gold2.json').read_bytes()
    assert sha256(data).hexdigest()=='5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'
    docs=json.loads(data); discount=Q(docs['upper']['certified_exact'])
    for side in ('upper','lower'):
        print('79-cap',side,flush=True)
        p=wire.Gold(private/f'quintet_{side}.bin')
        try:
            with b.Snapshot(side) as s:
                r=b.Rows(s,p,docs[side]); assert s.lam==Q(docs[side]['nominal'])
                assert any(check(lambda p:r.triple(*p),b,discount,((1,2),(0,),(0,)),p,False)
                           for p in product(((2,2),(3,3),(5,4)),
                                            ((D,c) for D in (79,19,78,20) for c in range(8,1,-1)),
                                            ((D,c) for D in (18,2,3,17,4,5) for c in range(min(D,3),0,-1)))
                           if p[0] in r.mids and any(x*y>0 for x,y in zip(r.triple(*p)[1],(p[0][1]-2,p[1][1]-1,p[2][1]-1))))
                s.check_unchanged()
        finally:
            b.os.close(p.fd)
    print('80–480 upper',flush=True)
    with cohort.host.open_host('upper') as h:
        r=cohort.Rows(h)
        assert any(check(r.affine,cohort,discount,((1,),(0,2),(1,)),p,True)
                   for p in product(((D,c) for D in (480,80,479,81,88,252) for c in range(12,1,-1)),((2,2),),
                                    ((D,c) for D in (12,2,3,11,4,5) for c in range(min(D,3),0,-1)))
                   if any(x*y for x,y in zip(r.affine(p)[1],(p[0][1]-1,0,p[2][1]-1))))


if __name__ == '__main__':
    main()
