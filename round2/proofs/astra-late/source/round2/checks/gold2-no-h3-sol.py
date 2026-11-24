#!/usr/bin/env python3
"""Verify the existing D1750 duals and compose 99/175 with no H3 credit or refit."""
import argparse
from collections import defaultdict
from decimal import Decimal, localcontext
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations_with_replacement, product
import json
from mmap import ACCESS_READ, mmap
from pathlib import Path
import re
from tarfile import open as open_tar

R = Path(__file__).resolve().parents[1]
G = R / 'occupation-quintet-gold2'
F = G / 'lifetime/pressure/frontier'
A = R / 'attacks/neon2-marginals'
S, GRID, SCALE = 2**40, 10**12, 10**6
CLAIM = 'cc8658e2e6354577716b44ff7f57bf6e3bb9831a1c7a7ac0ecf5c9bf717dbf9c'
GOLD = '5e11f9c0b6eb3bd43e880be608c37d6f0b1fdbcc632b91bb5e8037db52e13284'


def raw(path, pin):
    b = path.read_bytes()
    assert sha256(b).hexdigest() == pin, path
    return b


def doc(path, pin):
    return json.loads(raw(path, pin))


def load(path, pin):
    spec = spec_from_file_location(path.stem, path)
    mod = module_from_spec(spec)
    exec(compile(raw(path, pin), str(path), 'exec'), mod.__dict__)
    return mod


def claim_inputs(prior):
    archive = R / '.work/neon2-marginals-source/6c8afec.tar'
    pin = '2f2c2c1ee81577255aa4a8834b7c97850e7445f92008351bc5d2862e35cf7177'
    raw(archive, pin)
    with open_tar(archive) as tar:
        b = tar.extractfile('round2/neon2/stationary-marginals/certificates_all.json').read()
    assert sha256(b).hexdigest() == CLAIM
    evidence = doc(A/'replay.json', 'a3565adcb81ca54db62112c74a7c59b6056387e501691fc4c1badd14008891a6')
    assert evidence['source']['commit'] == '6c8afec5cd1b7f98bbcd7c6e7848c0d631a50f7d'
    assert evidence['source']['supplied_marginal_archive_sha256'] == pin
    assert evidence['replay']['byte_identical'] and evidence['replay']['output_sha256'] == CLAIM
    assert evidence['replay']['exit'] == evidence['signed_coboundary']['arithmetic_check_exit'] == 0
    j = json.loads(b)
    lo, hi = map(Q, j['mu_bracket_exact'])
    assert lo <= prior[0] <= prior[1] <= hi
    iv = []
    names = [f'p{d}_D{d}' for d in range(2, 7)] + ['pi_J_2_1', 'pi_bare_2_2', 'pi_bare_L0_fin']
    for i, name in enumerate(names, 1):
        o = j['observables'][str(i)]
        assert o['obs_name'] == name
        bounds = tuple(map(Q, o['certified_interval_exact']))
        for sign, side, bound in zip((-1, 1), ('lower', 'upper'), bounds, strict=True):
            t = o[side]
            v = t['verify']
            assert v['sha256'] == o[side+'_sha256'] == evidence['inputs']['snapshots'][f'obs{i}_{side}.bin']
            assert (v['obs_id'], v['obs_sign'], v['S'], v['DF']) == (i, sign, S, 1600)
            errs = []
            for k in ('lo', 'hi'):
                c = v['caps_'+k]
                f = max(c['f_in'], c['f_tab'], t['tail']['cap_ticks_'+k])
                errs.append(max(0, f + max(0, c['g']) + c['higher_per_core']))
            assert errs == evidence['signed_coboundary']['endpoint_errors_ticks'][f'obs{i}_{side}.bin']
            extra = Q(40, 2**1600)/lo if i == 8 and sign == 1 else 0
            assert Q(v['theta_P'], v['theta_Q']) + sign*Q(max(errs), S) + extra == bound == Q(t['certified_exact'])
        iv.append(bounds)
    (j0,j1),(p0,p1) = iv[5:7]
    cond = [(lo*(1/j1-1)-Q(3,2), hi*(1/j0-1)-Q(3,2)), (j0/(2*p1)-hi, j1/(2*p0)-lo)]
    for k, bounds in zip(('E_invD_parent_given_J','E_sum_invD_nbrs_given_2_2'), cond, strict=True):
        assert bounds == tuple(map(Q, j['corollaries'][k]['exact']))
    assert all(0 <= a <= b <= 1 for a,b in iv+cond)
    return iv, cond


def overlay(lp, iv, cond, only_lower_masses=False):
    names = [f'D={d}' for d in range(2,7)] + ['J','22','bare','rho|J','rho|22']
    degrees = dict(zip(names[:5], range(2,7), strict=True))
    bounds = dict(zip(names, iv+cond, strict=True))
    # The bare observable includes the physical degree tail. In the lift its
    # unresolved mass is a subset of the saturated overflow cell.
    def coefficients(x):
        d,k,_,_,u,v = x
        for name, (l,h) in bounds.items():
            den = 1
            if name in degrees: val = int(d == degrees[name])
            elif name == 'J': val = int((d,k) == (2,0))
            elif name == '22': val = int((d,k) == (2,1))
            elif name == 'bare': val = int(d <= 24 and k+1 == d)
            elif name == 'rho|J':
                den = int((d,k) == (2,0)); val = den*v
            else:
                den = int((d,k) == (2,1)); val = den*(u+v)
            yield name+' lo', val-l*den
            yield name+' hi', h*den-val
        yield 'bare tail outer', int((d,k) == (25,7))
    def tag(name): return ('marginal', name)
    lower = [f'D={d} lo' for d in range(2,7)]+['22 lo']
    keys = [tag(k) for k,_ in coefficients((2,0,0,0,0,0)) if not only_lower_masses or k in lower]
    for k in keys: lp.rows[k] = len(lp.rows)
    for i,(kind,col,x) in enumerate(lp.cols):
        if kind == 'occ':
            lp.cols[i] = kind, col+tuple((lp.rows[tag(k)],v) for k,v in coefficients(x) if v and tag(k) in lp.rows), x
    if not only_lower_masses:
        lp.col('bare-hidden', {tag('bare lo'):1, tag('bare hi'):-1, tag('bare tail outer'):-1})
    for k in keys: lp.col(('claim10-slack',k[1]), {k:-1})


def edge_ledger():
    def parse(path,pin,pat):
        [v] = re.findall(pat,raw(path,pin).decode(),re.M)
        return Q(v)
    old = dict(
        edge18=parse(G/'edges/edge18.log','1061d512219da5a38f264c07b4082f1db391b213a885f4cd955ea691946e2948',r'credit ticks/core (\S+)'),
        absorb12=S*parse(G/'edges/absorb12.log','eded191d098bb2c9296c5d7d2861a051bc303bce122953c518ebf1981d22d521',r'^TOTAL (\S+)'),
        old_payer=S*parse(G/'payer/README.md','0020b89ff3d30dcf654225d7f37626d62c8f21b60312b7b2b342036faac6904f',r'Exact normalized value: `(\d+/\d+)`'))
    edge = load(G/'true-drift/lower/check.py','6e08c9f7a31c1645f036e92449ebd914a8003f1927587d99101eeffacb2b322b')
    new = dict((k[0],v) for k,v in edge.RECEIPTS.items())
    assert list(edge.RECEIPTS) == [('wide-payer',251,12,125,3,5),('wide-host',251,18,34,8,5),('old-donor',251,11,40,5,9)]
    readme = raw(G/'true-drift/lower/README.md','9624ff55ce1efbb6a088b13f2c05da11b1190d872ff603c52d83aeab05e55e19').decode()
    assert list(map(Q,re.findall(r'^\| .*? \| `(\d+/\d+)` \|',readme,re.M))) == [*new.values(),edge.POST]
    assert sum(new.values())+edge.POST == Q(119326162107114155945908957363,6201700363750000000000000)
    return old, new, edge.POST, (edge.LO, edge.HI)


def masks():
    # Representatives of the joint threshold partition. Test even arity/degree
    # combinations that cannot occur; the masks only become smaller physically.
    counts = [0]*5
    for d,c,p,k,z in product((2,18,19,79,80,150,151,251,252), (1,2,10,11,12,13,17,18,19,22,23),
                             (2,5,6,18,19,34,35,40,41,125,126,252), (1,3,4,5,6,7,8,9), (2,5,6,9)):
        if max(d,p)<=18 or d>=252 or not 2<=c<=22 or (d<80 and c<=10): continue
        active = (d>18 and c<=17 and p<=5,
                  d<=150 and c<=10 and 6<=p<=34 and k<=7 and z<=5,
                  c<=12 and 35<=p<=125 and k<=3 and z<=5,
                  c<=18 and 6<=p<=34 and k<=8 and z<=5 and (d>150 or c>10 or k>7),
                  c<=11 and 6<=p<=40 and k<=5 and z>=6)
        assert sum(active)<=1
        counts = [n+a for n,a in zip(counts,active,strict=True)]
    assert all(counts)
    for d,p in combinations_with_replacement((2,3,4,10,11,18,19,25,26,251,252,253),2):
        a,b = sorted((d,p))
        assert sum((b<=18, 18<b<252 and (b<=25 or a<=3), a<=10 and b>=252, a>=252))<=1


def price(lp, witness, score, lo):
    assert (len(lp.rows),len(lp.cols),lp.rhs,witness['grid'],witness['scale']) == (
        witness['rows'],witness['checked_columns'],{'birth':1},GRID,SCALE)
    idx = {str(k):i for k,i in lp.rows.items()}
    assert len(idx) == len(lp.rows)
    sparse = dict(witness['dual'])
    assert len(sparse) == len(witness['dual']) and sparse.keys() <= idx.keys()
    assert all(type(x) is int and x != 0 for x in sparse.values())
    y = [0]*len(idx)
    for k,v in sparse.items(): y[idx[k]] = v
    defects = defaultdict(Q)
    for kind, column, x in lp.cols:
        c = score(x) if kind == 'occ' else 0
        assert c >= 0
        residual = sum(y[i]*a for i,a in column) - Q(GRID,SCALE)*c
        defects[str(kind)] = max(defects[str(kind)], residual)
    assert {k:v for k,v in defects.items() if v} == {k:Q(v) for k,v in witness['defects'].items()}
    floor = Q(SCALE,GRID)*(lo*y[idx['birth']]-sum(defects.values()))
    assert floor == Q(witness['price_ticks']) > 0
    return floor


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--check-host', action='store_true', help='recompute 307 cells on the authenticated original lower inode')
    parser.add_argument('--check-roots', action='store_true', help='rerun the original exact Bellman sweeps, not the LP fit')
    args = parser.parse_args()
    assert __debug__
    trial = doc(F/'repair/D1750-M5-T3-C12.json','5ffd94a84339edcaf48984aa3de19e5bbd79decda38bfacfb93db0a7871c0a78')
    root = doc(F/'marginals/lower-floor.json','95b05a81de95b5ef559414d006a23acc82f4790126db69d8509fe488277f979a')
    profile = doc(G/'host/conditioned24.json','3f98385e0662b4a5aea71dcfef03cdbd66012de01f0083456e243592007de790')
    gold = doc(R/'attacks/neon2-quintet/gold2/production.json', GOLD)
    g = gold['lower']; p = profile['results']['lower']
    prior = tuple(Q(gold[s]['certified_exact']) for s in ('lower','upper'))
    lo,hi = prior
    assert trial['prior'] == p['prior'] == root['prior'] == list(map(str,prior))
    assert root['gold2_receipt_sha256'] == profile['original_receipt_sha256'] == GOLD
    assert gold['lower_sha256'] == g['verify']['sha256'] == p['payload_sha256'] == root['gold2_payload_sha256']
    assert trial['claim10_sha256'] == CLAIM and trial['input_root_sha256'] == '95b05a81de95b5ef559414d006a23acc82f4790126db69d8509fe488277f979a'
    assert Q(root['discount']) == hi and Q(root['overflow_ticks']) == 0
    assert root['options'] == dict(mother=5,live_cut=700,live_arity=12,tip=3,donors=2,tag_limit=12)
    for path, digest in root['sources'].items():
        raw(R/path,digest)
    raw(F/'repair/replay.py','982ae0c93068356769f599c56409dcdbc7ffe853e267bcbd8e49ebb29d0bf5ab')
    assert (trial['states'],trial['corners']) == (7038933,17262)
    assert all(Q(v) >= Q(root['roots'][k]) >= 0 and (GRID*Q(v)).denominator == 1 for k,v in trial['roots'].items())
    roots = {tuple(map(int,k.split(','))):Q(v) for k,v in trial['roots'].items()}
    assert set(roots) == {(0,0)} | {(d,c) for d in range(2,6) for c in range(2,d+1)}
    old,new,post,edge_prior = edge_ledger()
    assert prior == edge_prior
    assert old == {k:Q(v) for k,v in p['edge_components_ticks'].items()}
    base = sum(old.values()); add = sum(new.values())
    assert base == Q(p['edge_floor']) == Q(trial['edge_ticks'])
    masks()
    lam = Q(g['nominal']); target = Q(trial['target']); delta = target-lam
    assert target == Q(99,175) and lo <= target <= hi and delta == Q(p['signed_shift']) > 0
    M,H,E = (g[k] for k in ('m_max','higher_per_core','err_per_core'))
    assert (g['cap_f'],g['cap_g'],H,E) == (55,792,18370,19217)
    assert E == g['cap_f']+max(0,g['cap_g'])+H == p['cap_error']
    assert (g['verify']['S'],g['verify']['caps']['g'],Q(p['nominal'])) == (S,g['cap_g'],lam)
    cells = {tuple(map(int,k.split(','))):tuple(map(int,v)) for k,v in p['cells'].items()}
    assert (p['limit'],p['saturation'],len(cells)) == (24,7,307) and not p['reduced_f']
    # The first comparison has an Euler correction absent from the upper one.
    assert all(g['m_min'] <= l and 0 < u <= M and 0 <= h <= H for l,u,h in cells.values())
    audit = load(R/'checks/gold2-host-h3-sol.py','ef62e72358be9f868e3a39ab996292f861b93ee264b1c7756c8364019a3697b7')
    if args.check_host or args.check_roots:
        host = audit.load(G/'host/reader.py')
        assert host.receipt() == gold
        with host.open_host('lower') as w:
            if args.check_host:
                assert cells == w.cells(24,7)
                caps = audit.cap_cells(w)
                assert all(h == caps[key] for key,(_,_,h) in cells.items())
                print('PASS all 307 original-payload comparison cells and independent higher cap counts',flush=True)
            if args.check_roots:
                path = G/'lifetime/pressure/replay.py'
                pressure = load(path,root['sources'][str(path.relative_to(R))])
                life = pressure.life
                with mmap(w.p.fd,0,access=ACCESS_READ) as mm:
                    t = life.tail.Snapshot(mm,w.golden)
                    W = {}
                    life.price(t,7999,30,hi,w.cap_f,False,(1751,W))
                    for cut in (700,1750):
                        values,n,k = pressure.live(t,w.cap_f,W,cut,12,5,3,2,hi,False)
                        baseline = {key:Q(root['roots'][f'{key[0]},{key[1]}']) for key in roots}
                        if cut == 700:
                            assert values == baseline and (n,k) == (root['states'],root['corners'])
                        else:
                            assert {key:max(v,baseline[key]) for key,v in values.items()} == roots
                            assert (n,k) == (trial['states'],trial['corners'])
                        print('PASS exact singleton roots at live cut',cut,'states/corners',n,k,flush=True)
    tags = load(R/'checks/gold2-frontier-sol.py','19aa4bb826f494fb08581fcd37ffcfbe26191aa89114eedee4bf51cf575dd06a').tags
    model = audit.load(R/'occupation-neon/pressure-memory/birth-law/replay.py')
    model.LO,model.HI,model.SAT = lo,hi,7
    iv, cond = claim_inputs(prior)
    def score(x):
        d,k,b,*_ = x
        assert b == 0 or d > 5 or (d,k+2) in roots
        _,U,h = cells[d,k]
        return b*roots.get((d,k+2),roots[0,0])+H-h+delta*(M-U)
    with localcontext() as ctx:
        ctx.prec = 35
        def show(label, n): print(label,n,f'({Decimal(n.numerator)/Decimal(n.denominator):.15f})',flush=True)
        show('original D18 + FY:',base)
        show('new FY only:',add)
        show('post-FY residual:',post)
        show('entire H3 credit excluded:',Q(trial['h3_ticks']))
        print('trial delta, M, H, E:',delta,M,H,E,flush=True)
        for label,shape,sign in [('Gold2 only',(623,5564),-1),('Claim10 conditional',(644,5586),1)]:
            lp = model.build(24,sided=True)
            tags(lp,12,prior)
            if sign == 1: overlay(lp,iv,cond)
            assert (len(lp.rows),len(lp.cols)) == shape
            w = trial['results'][label]
            if sign == 1:
                active = {k:v for k,v in w['dual'] if k.startswith("('marginal', ")}
                needed = [f'D={d} lo' for d in range(2,7)]+['22 lo']
                assert active.keys() == {str(('marginal',k)) for k in needed}
                assert all(v > 0 for v in active.values())
                assert set(w['defects']) == {'edge','occ'}
                print('Claim 10 dual uses only lower masses: D=2,3,4,5,6 and (D,c)=(2,2)',flush=True)
            P = price(lp,w,score,lo)
            if sign == 1:
                reduced = model.build(24,sided=True)
                tags(reduced,12,prior)
                overlay(reduced,iv,cond,only_lower_masses=True)
                assert (len(reduced.rows),len(reduced.cols)) == (629,5570)
                assert P == price(reduced,dict(w,rows=629,checked_columns=5570),score,lo)
                print('PASS same exact dual price with only the six Claim 10 lower mass rows',flush=True)
            show(label+' exact repaired joint floor:',P)
            no_post = base+add-E-delta*M+P
            margin = no_post+post
            # Old composition is only an independent accounting cross-check.
            assert base+Q(trial['h3_ticks'])-E-delta*M+P == Q(w['margin_ticks']) == -Q(w['deficit_ticks'])
            assert sign*margin > 0 and sign*no_post > 0
            show(label+' no-H3 margin:',margin)
            show(label+' no-H3/no-post-FY margin:',no_post)
    print('PASS two saved D1750 duals, Claim 10 endpoints, disjoint edge ledger and exact signed margins')


if __name__ == '__main__':
    main()
