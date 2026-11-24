"""Compose authenticated Gold2 host cells and matched original H3 floors."""
import argparse
import ast
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path
import re
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from reader import R2, RECEIPT, TARGET, pinned, receipt

S=2**40
HOST='3f98385e0662b4a5aea71dcfef03cdbd66012de01f0083456e243592007de790'
LOGS={
    'lower':('threshold-lower.historical.log','39d0d02575cb05df74776466d787cf8f8ef90399003ee0709e2e7c0e9d5e9191',Q(1084423030952180046412413909,1550425090937500000000)),
    'upper':('threshold-upper.log','c8609d448d63a9e425279b6adc157b83bee935fe8ca3510e9b3b6ac5cfb277f2',Q(3304204122287698944634328103,6201700363750000000000)),
}
LP='7dd40ac41e5d95ae41843da51de436fd45b781f1'


def bind(path,key):
    raw=path.read_bytes()
    if sha256(raw).hexdigest()!=key: raise ValueError(f'changed input: {path}')
    return raw


def inputs():
    original=receipt()
    profile=json.loads(bind(HERE.parent/'conditioned24.json',HOST))
    assert profile['original_receipt_sha256']==RECEIPT
    result={}
    for side,(name,digest,ticks) in LOGS.items():
        x=profile['results'][side]; gold=original[side]; v=gold['verify']
        assert x['limit']==24 and x['saturation']==7 and not x['reduced_f']
        assert x['payload_sha256']==v['sha256'] and x['nominal']==gold['nominal']
        assert tuple(map(Q,x['prior']))==tuple(Q(original[s]['certified_exact']) for s in LOGS)
        assert x['cap_error']==gold['err_per_core'] and Q(x['target'])==TARGET
        raw=bind(HERE.parents[1]/'post-fy'/name,digest).decode()
        options=ast.literal_eval(raw.splitlines()[0].split('original Gold2 ',1)[1])
        assert options['side']==side and options['kind'] in ('h3','threshold')
        assert not options['ceiling'] and not options['reserve']
        assert tuple(options[k] for k in ('host','stock','receiver','rstock','donor','tstock','old','young'))==(5,4,79,8,18,3,19,18)
        norm=re.search(r'^H3 normalized(?: orientations .*? total)? (\d+/\d+) ',raw,re.M)[1]
        assert Q(norm)*S==ticks and 'states 271350' in raw
        assert v['caps']['h3']==(408 if side=='lower' else 233)
        if side=='lower':
            components=x['edge_components_ticks']
            assert Q(components['old_payer'])==Q(1046834122639903380009,20000000000000000)
            assert Q(components['absorb12'])==Q(39526146756920663734311,10**17)
            assert sum(map(Q,components.values()))==Q(x['edge_floor'])
        else: assert Q(x['edge_floor'])==0
        result[side]=(x,gold,ticks)
    return result


def outward(x,d,places=15):
    z=d*((d*x*10**places)//1)
    return f'{z//10**places}.{z%10**places:0{places}d}'


def solve(side,x,gold,h3):
    d=1 if side=='lower' else -1
    lo,hi=map(Q,x['prior']); nominal=Q(x['nominal'])
    L,M,H=gold['m_min'],gold['m_max'],gold['higher_per_core']
    refined_lo,refined_hi=map(Q,x['denominator'])
    E=x['cap_error']; f,g=Q(x['f_floor']),Q(x['edge_floor']); b=f+g+h3-E
    assert 0<L<=refined_lo<=refined_hi<=M and b>0
    saved=Q(x['saved_plus_reduced']); assert 0<=saved<=H
    cells={tuple(map(int,k.split(','))):tuple(map(Q,v)) for k,v in x['cells'].items()}
    assert len(cells)==307 and all(L<=l and u<=M and 0<=h<=H for l,u,h in cells.values())
    scalar=nominal+d*b/M
    denom=nominal+d*(b+saved)/refined_hi
    scalar_target=b-d*(TARGET-nominal)*M
    if side=='lower':
        assert Q(x['margin'])+h3>0 and scalar_target>0

    def trial(s,restrict):
        # If the target conclusion already holds there is nothing to show;
        # in the contrary case the real measure has this narrower discount.
        low,high=(lo,min(hi,s)) if restrict and d==1 else (max(lo,s),hi) if restrict else (lo,hi)
        assert lo<=low<=high<=hi and lo<=s<=hi
        m=pinned('gold2_h3_trial_lp',R2/'occupation-neon/pressure-memory/birth-law/replay.py',LP)
        m.LO,m.HI,m.SAT=low,high,7
        lp=m.build(24,sided=True)
        assert (len(lp.rows),len(lp.cols))==(x['rows'],x['columns'])
        shift=d*(s-nominal); assert shift>=0
        def cost(kind,meta):
            if kind!='occ':return 0
            _,u,h=cells[meta[:2]]
            val=H-h+shift*(M-u)
            assert val>=0
            return val/10**7
        joint=max(0,low*10**7*lp.lower(cost))
        margin=b-shift*M+joint
        # A contradiction follows for mu at or beyond the trial in the adverse
        # direction, including equality. For a uniform-prior margin only,
        # one can also convert the margin into an explicit inward offset.
        print(side,'restricted' if restrict else 'uniform','trial',s,'margin',float(margin),flush=True)
        return dict(point=str(s),restricted=restrict,prior=list(map(str,(low,high))),
                    rebate=str(joint),margin=str(margin),proved=margin>0)

    p=trial(denom,False)
    uniform=denom+d*max(0,Q(p['margin']))/refined_hi
    # Deliberately test one stronger fixed point, just beyond the uniform
    # ratio. The outward grid point and margin are separately auditable.
    fixed=Q(d*((d*uniform*10**12).__ceil__()+1),10**12)
    z=trial(fixed,True)
    final=fixed if z['proved'] else uniform
    return dict(nominal=str(nominal),original_prior=x['prior'],payload_sha256=x['payload_sha256'],
                h3_log=LOGS[side][0],h3_log_sha256=LOGS[side][1],h3_family_cap=gold['verify']['caps']['h3'],
                f_floor=str(f),edge_floor=str(g),h3_floor=str(h3),error=E,higher_cap=H,net=str(b),
                original_denominator=list(map(str,(L,M))),d24_denominator=x['denominator'],saved=str(saved),
                scalar_endpoint=str(scalar),scalar_target_margin=str(scalar_target),
                d24_separate_endpoint=str(denom),uniform_trial=p,uniform_endpoint=str(uniform),
                stronger_fixed_trial=z,final_endpoint=str(final),outward=outward(final,d))


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--write',action='store_true');a=p.parse_args()
    out={side:solve(side,*v) for side,v in inputs().items()}
    assert Q(out['lower']['final_endpoint'])>TARGET and Q(out['lower']['final_endpoint'])<Q(out['upper']['final_endpoint'])
    doc=dict(profile='Gold2 D24 with conditioned f, lower edge rebates and matched original H3',
             host_sha256=HOST,original_receipt_sha256=RECEIPT,lp_git_blob=LP,target=str(TARGET),results=out)
    raw=(json.dumps(doc,indent=2,sort_keys=True)+'\n').encode()
    if a.write: (HERE/'receipt.json').write_bytes(raw)
    else: assert (HERE/'receipt.json').read_bytes()==raw
    print('outward',*(v['outward'] for v in out.values()),'receipt SHA-256',sha256(raw).hexdigest())
    print('PASS exact both-sign Gold2 H3+D24 composition')


if __name__=='__main__':main()
