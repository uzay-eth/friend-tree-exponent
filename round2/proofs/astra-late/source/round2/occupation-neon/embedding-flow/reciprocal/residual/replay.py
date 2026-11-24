"""Matched all-degree reciprocal boxes/moments for the frozen post-EDGE residual."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import product
from pathlib import Path
import sys
from time import monotonic, process_time


HERE = Path(__file__).resolve().parent
OCC = HERE.parents[2]
PROJ = OCC/'interpolate/compose/residual/edge-projection'
sys.path.insert(0, str(HERE.parent))
from model import Model
from run import optimize

LO = Q('.565678056632067')
HI = Q('.565770469400619')
MASK = (18, 4)
DEN = {2: Q('1357493815758092773/1000000000000000000'),
       1: Q('273436004214812983/200000000000000000')}
BOOKED = {2: Q('18402664634233592780094779/178500000000000000000000000000000'),
          1: Q('101984612445685206436931/1700000000000000000000000000000')}
ZERO = Q(0)


def load(path, name):
    spec = spec_from_file_location(name, path)
    mod = module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


def sources():
    frozen = load(PROJ/'combined/replay.py', 'resolved_frozen_combined')
    edge_dir = OCC/'interpolate/edge'
    for manifest in (PROJ/'combined/release-inputs.sha256', PROJ/'combined/inputs.sha256',
                     edge_dir/'inputs.sha256', edge_dir.parent/'compose/inputs.txt'):
        frozen.pin_manifest(manifest)
    p = load(PROJ/'replay.py','resolved_frozen_projection')
    for rel, want in p.PINS.items():
        assert sha256((p.ROOT/rel).read_bytes()).hexdigest()==want, rel
    e = p.module(next(iter(p.PINS)), 'resolved_frozen_edge')
    return e,p,frozen


def rewards(e,p,frozen,path,side,fmax,gmax):
    mode = 2 if side=='lower' else 1
    before = path.stat()
    assert p.digest(path)==e.PAYLOAD_HASHES[mode]
    c = e.Charges(e.D.Cert(path))
    assert c.c.mode==mode and c.sign==2*mode-3
    trial,n,k = p.SPEC[mode]
    s = Q(trial)
    edge, census = e.Edge(c,n,k,'disjoint'), e.Census(c.c)
    v = edge.solve(s)
    margin,_,z,tail = census.solve(s,v)
    denominator = e.DEN[mode]+max(c.sign*(w-z[e.J]) for w in (*z.values(),*tail.values()))
    if c.sign==1:
        denominator += max(v.values())
    assert denominator==DEN[mode]
    net = margin+e.RETAINED[mode]
    assert margin<0<net and s+c.sign*net/denominator==frozen.BASE[mode]
    post = s+c.sign*(net+BOOKED[mode])/denominator
    assert (LO<=post if mode==2 else post<=HI)
    f = p.freward(e,c,census,s,v,z,tail,fmax)
    g = p.greward(e,c,edge,s,v,gmax)
    assert len(f)==sum(range(2,fmax+1))
    assert len(g)==sum(range(1,gmax))*sum(range(2,gmax+1))
    after = path.stat()
    assert (before.st_ino,before.st_size,before.st_mtime_ns,before.st_ctime_ns)==(
        after.st_ino,after.st_size,after.st_mtime_ns,after.st_ctime_ns)
    print(side,'authenticated',e.PAYLOAD_HASHES[mode],'EDGE margin',margin,'retained',e.RETAINED[mode],
          'denominator',denominator,'ead947b already booked',BOOKED[mode],flush=True)
    return mode,f,g


def separate_costs(model, f, g, fcut, gcut):
    fc,gc = defaultdict(Q),defaultdict(Q)
    count_f=count_g=0
    for cells in model.roots:
        x, = cells
        if not 0<x.degree<=fcut:
            continue
        a,b = f[x.label]
        model.add_observation(fc,cells,(0,),a)
        model.add_observation(fc,cells,(0,),b,kind=(0,'mother'))
        if x.arity:
            model.add_observation(fc,cells,(0,),b,kind=(0,'children'))
        assert min(a,a+b*Q(x.arity+1,2))>=0
        count_f+=1
    for cells in model.edges:
        p,x = cells
        if not (0<p.degree<=gcut and 0<x.degree<=gcut):
            continue
        a,b,d = g[p.label,x.label]
        model.add_observation(gc,cells,(0,0),a)
        # u = mother's pressure + sibling pressure at parent, v = child's own child pressure.
        model.add_observation(gc,cells,(0,0),b,kind=(0,'mother'))
        if p.arity>1:
            model.add_observation(gc,cells,(0,0),b,kind=(0,'children'))
        if x.arity:
            model.add_observation(gc,cells,(0,0),d,kind=(1,'children'))
        assert a+min(ZERO,b*Q(p.arity,2))+min(ZERO,d*Q(x.arity,2))>=0
        count_g+=1
    assert count_f==sum(range(2,fcut+1)) and count_g==sum(range(1,gcut))*sum(range(2,gcut+1))
    return dict(fc),dict(gc)


def pay_by_name(model,cost):
    return {model.variables[i][0]: x for i,x in cost.items() if x}


def check_objective(model, fc, gc, projector, f, g, fcut, gcut):
    """Compare each paid affine row with the frozen evaluator at physical box corners."""
    evaluate = projector.objective({k:v for k,v in f.items() if k[0]<=fcut},
                                   {k:v for k,v in g.items() if max(p[0] for p in k)<=gcut})
    observed = set()
    checked = 0
    for cells in model.roots+model.edges:
        if any(c.degree==0 for c in cells):
            continue
        idx = (0,)*len(cells)
        parent = cells[0]
        ports = {(0,'mother'):Q(1,2)}
        if model.spare(cells,0):
            ports[0,'children'] = Q(model.spare(cells,0),2)
        if len(cells)==2 and cells[1].arity:
            ports[1,'children'] = Q(cells[1].arity,2)
        cost = fc if len(cells)==1 else gc
        masses = model.blocks[cells,'mass'][idx]
        assert len(set(cost.get(i,ZERO) for i in masses))==1
        observed.update(masses)
        cols = {p:model.blocks[cells,p][idx][0] for p in ports}
        observed.update(cols.values())
        for corner in product(*((ZERO,max_) for max_ in ports.values())):
            values = dict(zip(ports,corner,strict=True))
            actual = cost.get(masses[0],ZERO) + sum(cost.get(cols[p],ZERO)*x for p,x in values.items())
            u = values[0,'mother'] + values.get((0,'children'),ZERO)
            if len(cells)==1:
                expected = evaluate('occ',(parent.label,None,u,None))
            else:
                child = cells[1]
                p,w = (parent.degree,parent.degree,parent.arity),(child.degree,child.degree,child.arity)
                expected = evaluate('edge',(p,w,Q(1,p[0]),Q(1,w[0]),u,values.get((1,'children'),ZERO)))
            assert actual==expected/10**6>=0
            checked+=1
    assert all(i in observed for i,x in add(fc,gc).items() if x)
    return checked


def add(f,g):
    return {i:f.get(i,ZERO)+g.get(i,ZERO) for i in f.keys()|g.keys()}


def price(model,cost,method,seconds):
    start = process_time(),monotonic()
    if seconds:
        # Use the common all-column exact repair for bounded LP proposals.
        pricing_source = HERE.parents[1]/'pricing/replay.py'
        pricing = sys.modules.get('resolved_residual_pricing') or load(pricing_source,'resolved_residual_pricing')
        value = pricing.Pricing(model).run(cost,'ipm' if method=='highs-ipm' else 'ds',seconds,1,1)
        return model.low*max(ZERO,value.value),value.repair,process_time()-start[0],monotonic()-start[1]
    last = None
    for proposal in dict.fromkeys((method,'highs-ds')):
        try:
            value,repair = optimize(model,cost,False,method=proposal,scale=10**8)
            break
        except ArithmeticError as fail:
            print('proposal failed',proposal,fail,flush=True)
            last=fail
    else:
        raise last
    # Physical supported rows are jointly nonnegative, so a failed-to-improve
    # dual lower bound can safely be replaced by the pointwise zero floor.
    return model.low*max(ZERO,value),repair,process_time()-start[0],monotonic()-start[1]


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--lower',type=Path,required=True)
    ap.add_argument('--upper',type=Path,required=True)
    ap.add_argument('--side',choices=('lower','upper'))
    ap.add_argument('--cuts',nargs='+',default=['5:4:4','7:6:4'],metavar='H:F:G')
    ap.add_argument('--arity',type=int,default=4)
    ap.add_argument('--order',type=int,default=1)
    ap.add_argument('--slabs',type=int,default=1)
    ap.add_argument('--method',choices=('highs-ipm','highs-ds'),default='highs-ipm')
    ap.add_argument('--seconds',type=float,help='time limit per LP proposal with checked fallback')
    ap.add_argument('--interval',choices=('edge-memory','original'),default='original')
    ap.add_argument('--variants',nargs='+',choices=('box','box-moments','share','moments'),
                    default=['box','share','moments'])
    ap.add_argument('--split',action='store_true',help='also solve individual f and g on the same models')
    ap.add_argument('--audit-only',action='store_true',help='verify extraction and all objective corners without LP search')
    a = ap.parse_args()
    assert __debug__
    cuts = [tuple(map(int,c.split(':'))) for c in a.cuts]
    if not all(2<=f<=18 and 2<=g<=17 and max(f,g)<h for h,f,g in cuts):
        ap.error('require 2<=F<=18, 2<=G<=17 and max(F,G)<H')
    if a.seconds is not None and not 0<a.seconds<=300:
        ap.error('seconds must be between 0 and 300')
    low,high = (LO,HI) if a.interval=='edge-memory' else (Q('.565673570322'),Q('.565775'))
    print('discount',a.interval,low,high,flush=True)
    print('degree and arity census hashes',*(sha256((HERE.parent/n).read_bytes()).hexdigest() for n in ('model.py','run.py')),flush=True)
    start = process_time(),monotonic()
    e,p,frozen = sources()
    tables = {}
    for side in ('lower','upper'):
        if a.side is None or a.side==side:
            tick = process_time(),monotonic()
            tables[side] = rewards(e,p,frozen,getattr(a,side),side,max(x[1] for x in cuts),max(x[2] for x in cuts))
            print('frozen extraction CPU/wall',*(round(v,3) for v in (process_time()-tick[0],monotonic()-tick[1])),flush=True)
    for h,fcut,gcut in cuts:
        by_side = {side:{} for side in tables}
        names = {}
        for label in a.variants:
            m = a.order if label in ('box-moments','moments') else 0
            slabs = a.slabs if label in ('box-moments','moments') else 1
            model = Model(h,a.arity,m,slabs,edge_balances=True,share=label in ('share','moments'),low=low,high=high)
            print('H:F:G',h,fcut,gcut,label,'arity/order/slabs',a.arity,m,slabs,'columns/eq/ineq',
                  len(model.variables),len(model.equalities),len(model.upper_rows),flush=True)
            for side,(mode,f,g) in tables.items():
                fc,gc = separate_costs(model,f,g,fcut,gcut)
                objective = (pay_by_name(model,fc),pay_by_name(model,gc))
                if side not in names:
                    names[side] = objective
                    print(side,'exact frozen objective corner checks',check_objective(model,fc,gc,p,f,g,fcut,gcut),flush=True)
                else:
                    assert names[side]==objective, 'objective differs between models'
                if a.audit_only:
                    continue
                vals = []
                for label2,cost in [('joint',add(fc,gc))]+([('f',fc),('g',gc)] if a.split else []):
                    v,r,cpu,wall = price(model,cost,a.method,a.seconds)
                    vals.append(v)
                    print(side,label,label2,'exact floor',v,'decimal',float(v),'repair before conversion',float(r),
                          'solver CPU/wall',round(cpu,3),round(wall,3),flush=True)
                by_side[side][label]=vals[0]
        for side,vals in by_side.items():
            if a.audit_only:
                continue
            mode = tables[side][0]
            best = max(vals.values())
            print(side,'best',float(best),'share of ead947b reward',float(best/BOOKED[mode]),
                  'exponent-equivalent before legacy debit',float(best/DEN[mode]),
                  'gain vs matched boxes',float(best-vals['box']) if 'box' in vals else 'n/a',
                  'positive difference vs ead947b reward (comparison only)',float(max(ZERO,best-BOOKED[mode])),
                  'frozen mask exact', (fcut,gcut)==MASK,flush=True)
    print('total CPU/wall',round(process_time()-start[0],3),round(monotonic()-start[1],3),flush=True)


if __name__=='__main__':
    main()
