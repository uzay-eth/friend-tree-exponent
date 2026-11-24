#!/usr/bin/env python3
"""Cold replay of original GOLD2 D24 plus matched 2417 H3."""
import argparse
import ast
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
import json
from math import ceil, comb, floor
from pathlib import Path
import re

R = Path(__file__).resolve().parents[1]
G = R/'occupation-quintet-gold2'
SIDES = ('lower','upper')
TARGET = Q(99,175)
S = 2**40
PINS = {
    G/'host/reader.py':'6be2e8aa5835ad0f1087179af5d84664d1176ed3ab071ca64f58af8d628fc910',
    G/'host/conditioned24.json':'3f98385e0662b4a5aea71dcfef03cdbd66012de01f0083456e243592007de790',
    G/'host/h3-joint/receipt.json':'ae6f5cde2ef710b4196c3c445a3786580e5ccc9dfd9206b2962f60628866c6ed',
    G/'lifetime/README.md':'17ad57b53161e61517a8ea776fa536abec39a76f967a1c57fe044b14fb418a5c',
    G/'edges/edge18.log':'1061d512219da5a38f264c07b4082f1db391b213a885f4cd955ea691946e2948',
    G/'edges/absorb12.log':'eded191d098bb2c9296c5d7d2861a051bc303bce122953c518ebf1981d22d521',
    G/'payer/README.md':'0020b89ff3d30dcf654225d7f37626d62c8f21b60312b7b2b342036faac6904f',
    G/'post-fy/threshold-lower.historical.log':'39d0d02575cb05df74776466d787cf8f8ef90399003ee0709e2e7c0e9d5e9191',
    G/'post-fy/threshold-upper.log':'c8609d448d63a9e425279b6adc157b83bee935fe8ca3510e9b3b6ac5cfb277f2',
    R/'occupation-neon/pressure-memory/birth-law/replay.py':'92a99ef12d8567798c4ed4dc3beb6535d3f18ef3ad5028aff0ec84e3684a5ba8',
}


def raw(p,pin=None):
    data=p.read_bytes()
    assert sha256(data).hexdigest()==(pin or PINS[p]),p
    return data


def load(p):
    spec=spec_from_file_location(p.stem,p)
    m=module_from_spec(spec)
    exec(compile(raw(p),str(p),'exec'),m.__dict__)
    return m


def outward(x,side):
    n=(floor if side=='lower' else ceil)(x*10**15)
    return f'{n//10**15}.{n%10**15:015d}'


def credits(side):
    def one(path,pat):
        [v]=re.findall(pat,raw(G/path).decode(),re.M)
        return Q(v)
    f=re.search(r'Conditioned floor \| \*\*([\d.]+)\*\* \| \*\*([\d.]+)\*\*',raw(G/'lifetime/README.md').decode()).groups()
    scan='post-fy/threshold-'+('lower.historical' if side=='lower' else 'upper')+'.log'
    options=ast.literal_eval(raw(G/scan).decode().splitlines()[0].split('original Gold2 ',1)[1])
    assert options['side']==side and not options['reserve'] and not options['ceiling']
    assert tuple(options[k] for k in ('host','stock','receiver','rstock','donor','tstock','old','young'))==(5,4,79,8,18,3,19,18)
    h3=S*one(scan,r'^H3 normalized (?:orientations .* total )?(\d+/\d+) ')
    edge={}
    if side=='lower':
        edge={
            'edge18':one('edges/edge18.log',r'credit ticks/core (\S+)'),
            'absorb12':S*one('edges/absorb12.log',r'^TOTAL (\S+)'),
            'old_payer':S*one('payer/README.md',r'Exact normalized value: `(\d+/\d+)`'),
        }
    return Q(f[SIDES.index(side)]),edge,h3


def cap_cells(w):
    """Count valid internal roles independently, using original family maxima."""
    p=w.p
    fam={}
    for name,tag in (('h3','3m'),('s3','S3'),('s4','S4')):
        mult={'h3':2,'s3':3,'s4':4}[name]
        fam[name]={mark:comb(mark[1],mult) for mark in p.ins[tag]}
    roles={t:list(d) for t,d in p.ins.items()}
    for name,tags in (('p4',('4i','4i')),('p5',('5i','5i','5i')),('b5',('5o','5c'))):
        best=defaultdict(int)
        for block in p.block[name]:
            marks=[roles[tag][i] for tag,i in zip(tags,block,strict=True)]
            if name=='p4':
                a,b=marks
                for x,y in ((a,b),(b,a)): best[x]=max(best[x],comb(x[1],2)*(y[1]-1))
            elif name=='p5':
                a,b,c=marks
                best[b]=max(best[b],comb(b[1],2)*(a[1]-1)*(c[1]-1))
            else:
                a,b=marks
                best[a]=max(best[a],a[1]*comb(a[1]-1,2)*(b[1]-1))
        fam[name]=best
    physical=defaultdict(int)
    caps=w.golden['verify']['caps']
    for name,tab in fam.items():
        for mark,n in tab.items(): physical[mark]+=n*max(0,caps[name])
    out=defaultdict(int)
    for (d,c),v in physical.items():
        key=(d,c-1) if d<=24 else (25,min(c-1,7))
        out[key]=max(out[key],min(w.higher_cap,v))
    return out


def price(lp,lo,score,scale):
    from scipy.optimize import linprog
    c=[Q(score(m),scale) if kind=='occ' else Q(0) for kind,_,m in lp.cols]
    assert min(c)>=0 and lp.rhs=={'birth':1}
    A,rhs=lp.matrix()
    y=linprog(list(map(float,c)),A_eq=A,b_eq=rhs,method='highs-ipm',
              options=dict(dual_feasibility_tolerance=1e-9,primal_feasibility_tolerance=1e-9))
    if not y.success: y=linprog(list(map(float,c)),A_eq=A,b_eq=rhs,method='highs-ds')
    assert y.success,y.message
    grid=10**12
    nums=[round(v*grid) for v in y.eqlin.marginals]
    defects=dict.fromkeys(('occ','edge','stubs','tail-up'),Q(0))
    for (kind,terms,_),v in zip(lp.cols,c,strict=True):
        defects[kind]=max(defects[kind],sum(nums[i]*a for i,a in terms)-grid*v)
    return scale*max(Q(0),(lo*nums[lp.rows['birth']]-sum(defects.values()))/grid)


def interval(side,old,f,edge,h3):
    direction=1 if side=='lower' else -1
    lam=Q(old['nominal'])
    Mlo,Mhi=map(Q,old['denominator'])
    Hsave=Q(old['saved_plus_reduced'])
    assert 0<Mlo<Mhi and (old['rows'],old['columns'],len(old['cells']))==(601,5542,307)
    assert not old['reduced_f'] and Q(old['f_floor'])==f
    assert {k:Q(v) for k,v in old.get('edge_components_ticks',{}).items()}==edge
    assert Q(old['edge_floor'])==sum(edge.values())
    net0=f+sum(edge.values())-old['cap_error']+Hsave
    assert net0>0 and h3>0
    assert lam+direction*net0/Mhi==Q(old['separate_endpoint'])
    net=net0+h3
    a=direction*(TARGET-lam)
    assert a==Q(old['signed_shift'])>0
    sep=net-a*Mhi
    margin=Q(old['margin'])+h3
    assert margin>=sep
    exact=lam+direction*net/Mhi
    assert exact==Q(old['separate_endpoint'])+direction*h3/Mhi
    print(side,'matched H3 ticks',h3,'net ticks',net,'M upper ticks',Mhi,flush=True)
    print(side,'endpoint',exact,'outward',outward(exact,side),
          'target separate/joint margin',sep,margin,flush=True)
    return exact,margin


def replay_final(side,w,model,cells,old,record,f,edge,h3,separate):
    d=w.sign; lam=w.lam; L,M,H=w.global_low,w.global_high,w.higher_cap
    lo,hi=w.prior
    net=f+sum(edge.values())+h3-old['cap_error']
    dhi=Q(old['denominator'][1]); saved=Q(old['saved_plus_reduced'])
    assert net>0 and 0<=saved<=H and 0<L<=Q(old['denominator'][0])<=dhi<=M
    exact_fields=dict(nominal=lam,net=net,saved=saved,h3_floor=h3,f_floor=f,edge_floor=sum(edge.values()),
                      scalar_endpoint=lam+d*net/M,scalar_target_margin=net-d*(TARGET-lam)*M,
                      d24_separate_endpoint=separate,error=old['cap_error'],higher_cap=H,
                      h3_family_cap=w.golden['verify']['caps']['h3'])
    assert all(Q(record[k])==v for k,v in exact_fields.items())
    assert record['d24_denominator']==old['denominator'] and record['original_prior']==old['prior']
    assert record['original_denominator']==list(map(str,(L,M)))
    assert record['payload_sha256']==old['payload_sha256']
    log='threshold-'+('lower.historical' if d==1 else 'upper')+'.log'
    assert record['h3_log']==log and record['h3_log_sha256']==PINS[G/'post-fy'/log]

    def trial(point,restricted,key):
        r=record[key]
        assert (Q(r['point']),r['restricted'])==(point,restricted)
        prior=(lo,hi) if not restricted else (lo,point) if d==1 else (point,hi)
        assert lo<=point<=hi and tuple(map(Q,r['prior']))==prior
        model.LO,model.HI,model.SAT=*prior,7
        lp=model.build(24,sided=True)
        assert (len(lp.rows),len(lp.cols))==(601,5542)
        a=d*(point-lam); assert a>0
        def score(t):
            _,u,h=cells[t[:2]]
            return H-h+a*(M-u)
        rebate=price(lp,prior[0],score,10**7)
        margin=net-a*M+rebate
        assert (rebate,margin,margin>0)==(Q(r['rebate']),Q(r['margin']),r['proved'])
        print(side,key,'exact margin',margin,'approx',float(margin),flush=True)
        return margin

    uniform=separate+d*max(0,trial(separate,False,'uniform_trial'))/dhi
    assert uniform==Q(record['uniform_endpoint'])
    fixed=d*Q(ceil(d*uniform*10**12)+1,10**12)
    assert d*(fixed-uniform)>0
    margin=trial(fixed,True,'stronger_fixed_trial')
    final=fixed if margin>0 else uniform
    assert Q(record['final_endpoint'])==final and record['outward']==outward(final,side)
    print(side,'final',final,'outward',record['outward'],flush=True)
    return final


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--side',choices=SIDES,action='append')
    p.add_argument('--arithmetic-only',action='store_true')
    a=p.parse_args(); assert __debug__
    host=load(G/'host/reader.py')
    base=json.loads(raw(G/'host/conditioned24.json'))
    assert base['original_receipt_sha256']==host.RECEIPT
    docs=host.receipt()
    final=json.loads(raw(G/'host/h3-joint/receipt.json'))
    assert final['host_sha256']==PINS[G/'host/conditioned24.json'] and final['original_receipt_sha256']==host.RECEIPT
    assert Q(final['target'])==TARGET
    got={}
    for side in a.side or SIDES:
        old=base['results'][side]; doc=docs[side]
        assert old['payload_sha256']==docs[side+'_sha256']
        assert old['cap_error']==doc['err_per_core'] and Q(old['nominal'])==Q(doc['nominal'])
        f,edge,h3=credits(side)
        separate,_=interval(side,old,f,edge,h3)
        if a.arithmetic_only: continue
        with host.open_host(side) as w:
            cells=w.cells(24,7)
            assert cells=={tuple(map(int,k.split(','))):tuple(map(int,v)) for k,v in old['cells'].items()}
            caps=cap_cells(w)
            assert all(h==caps[k] for k,(_,_,h) in cells.items())
            m=load(R/'occupation-neon/pressure-memory/birth-law/replay.py')
            lo,hi=w.prior
            assert (lo,hi)==tuple(map(Q,old['prior']))
            m.LO,m.HI,m.SAT=lo,hi,7
            lp=m.build(24,sided=True)
            assert (len(lp.rows),len(lp.cols))==(601,5542)
            L,M,H=w.global_low,w.global_high,w.higher_cap
            def bound(fn,scale=S): return price(lp,lo,lambda t:fn(*cells[t[:2]]),scale)
            dlo=L+bound(lambda l,u,h:l-L)
            dhi=M-bound(lambda l,u,h:M-u)
            save=bound(lambda l,u,h:H-h,10**7)
            shift=Q(old['signed_shift'])
            j=bound(lambda l,u,h:H-h+shift*(M-u),10**7)
            assert (dlo,dhi)==tuple(map(Q,old['denominator'])) and save==Q(old['saved_plus_reduced'])
            assert j==Q(old['joint_rebate'])
            assert Q(old['margin'])==max(f+sum(edge.values())-old['cap_error']+j-shift*M,
                                        f+sum(edge.values())-old['cap_error']+save-shift*dhi)
            print(side,'307 census cells; exact paired denominator, cap savings, and joint objective',
                  '601 rows / 5542 columns (four exact repairs)',flush=True)
            got[side]=replay_final(side,w,m,cells,old,final['results'][side],f,edge,h3,separate)
    if set(got)==set(SIDES):
        assert TARGET<got['lower']<got['upper']
        print('PASS independent GOLD2 H3+D24',*(outward(got[s],s) for s in SIDES))


if __name__=='__main__': main()
