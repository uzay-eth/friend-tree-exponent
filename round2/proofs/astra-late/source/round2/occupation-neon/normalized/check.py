"""Rational dual bounds for normalized physical stable-age occupation."""

import argparse
import json
from fractions import Fraction as Q
from pathlib import Path


LO, HI = Q('.565673570322'), Q('.565775')
MINUS, PLUS = Q('.565673'), Q('.565775')
UNIT = 10**6  # report rewards in micro-units
SCALE = 10**16
PATH = Path(__file__).resolve().parents[1] / 'rows' / 'charged-f.json'
BANDS = (
    ('256918177771749/100000000000000000000', '4272535757583/4000000000000000000',
     '6713770119873/5000000000000000000', '932799090941/625000000000000000',
     '78795139786299/50000000000000000000'),
    ('20697884432099/10000000000000000000', '861071198001/1000000000000000000',
     '5412445357533/5000000000000000000', '6016226911043/5000000000000000000',
     '31764287097597/25000000000000000000'),
)


def down(q):
    return Q(q*SCALE//1, SCALE)


def inner_rows(cap):
    data = json.loads(PATH.read_text())
    assert data['payload_sha256'] == {
        'lower':'9921eee052c489f07e8eba63d65c1dac5e85e4e495d8427e3d67a1a0642025ec',
        'upper':'1c4dec226208b02e2126dac5dcd911cfd714635bdaf322eac01e5ef0735fda41'}
    sides = []
    for name, sign, lam in (('lower',1,MINUS),('upper',-1,PLUS)):
        item = data[name]
        assert (item['D0'],item['S'],item['sign'],Q(item['lambda']))==(18,10**11,sign,lam)
        rows = {}
        for row in item['rows']:
            l,k,h = row['state']
            d = l+k+h
            if d>cap:
                continue
            assert h==1 and d<18
            p,c = map(Q,row['center'])
            f = row['f']
            val = (l+p+k*c)*(row['f_gain']-f)-lam*f
            if l:
                g = sum(Q(v)*Q(w) for v,w in row['g_bud_newborn'])
                val += Q(l,d)*(row['f_bud']+row['f_newborn']+g-f)
            value = sign*val/item['S']
            assert value==Q(row['signed_corners'][0])>=0
            rows[d,k] = UNIT*value
        assert len(rows)==sum(range(2,cap+1))
        sides.append(rows)
    return sides


class Age:
    def __init__(self, cap, overflow, parent=2):
        self.start = (2,0)
        self.rows = {}
        tail = cap+1
        for d in reversed(range(2,tail+1)):
            for k in reversed(range((overflow+1) if d==tail else d)):
                moves = []
                if d==tail:
                    bs = (Q(0) if k==overflow else max(Q(0),1-Q(k+1,tail)),Q(1))
                    for b in bs:
                        moves.append((b, ((b,(tail,k+1)),) if k<overflow else ()))
                else:
                    l = d-k-1
                    b = Q(l,d)
                    for g in (l,l+Q(1,parent)+Q(k,2)):
                        moves.append((b, tuple((r,y) for r,y in ((b,(d,k+1)),(g,(d+1,k))) if r)))
                self.rows[d,k] = [(s,b,jumps,s+sum(r for r,_ in jumps))
                                  for s in (LO,HI) for b,jumps in moves]

    def dual(self, payoff, t):
        v = {}
        for x,actions in self.rows.items():
            f = payoff.get(x,0)
            v[x] = down(min((s*(f+t*b)+sum(r*v[y] for r,y in jumps))/den
                            for s,b,jumps,den in actions))
            assert all(den*v[x] <= s*(f+t*b)+sum(r*v[y] for r,y in jumps)
                       for s,b,jumps,den in actions)
        return v[self.start]-max(t*LO,t*HI)

    def lower(self, payoff):
        base = self.dual(payoff,Q(0))
        best = base,Q(0)
        scale = 8*max((abs(x) for x in payoff.values()),default=Q(0))
        a,b = -scale,scale
        # This proposes dual multipliers; the exact Bellman inequalities
        # above certify every proposal independently of the search.
        for _ in range(14):
            x,y = (2*a+b)/3,(a+2*b)/3
            u,v = self.dual(payoff,x),self.dual(payoff,y)
            best = max(best,(u,x),(v,y))
            if u<v:
                a = x
            else:
                b = y
        return base,best[0],best[1]

    def interval(self, f):
        return self.lower(f),self.lower({x:-v for x,v in f.items()})


def ends(lo,hi):
    unit = 10**10
    def decimal(tick):
        a,b = divmod(abs(tick),unit)
        return f'{"-" if tick<0 else ""}{a}.{b:010d}'
    return f'[{decimal(lo*unit//1)}, {decimal(-(-hi*unit//1))}]'


def show(model, name, f):
    low, high = model.interval(f)
    print(name, 'flow only',ends(low[0],-high[0]),'with pi(b)',ends(low[1],-high[1]),
          'dual t', *(f'{float(v[2]):.6g}' for v in (low,high)),flush=True)
    return low[1],-high[1],low[2],high[2]


def birth_floor(d):
    prob = Q(1)
    for l in range(1,d):
        wanted = Q(l,l+1) if l==d-1 else l
        prob *= wanted/(HI+l+Q(l,l+1)+Q(1,2))
    return prob


def birth_mix(models, f, t):
    base = models[0][1].dual(f,t)
    assert sum(w for w,_ in models[1:])<=1
    # Class lower bounds, with unassigned birth mass kept at the baseline.
    return base+sum(w*max(Q(0),m.dual(f,t)-base) for w,m in models[1:])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cap',type=int,choices=range(3,10),default=5)
    args = parser.parse_args()
    cap = args.cap
    m = Age(cap,cap)
    minus,plus = inner_rows(cap)
    gamma = minus[2,0]/plus[2,0]
    total = {x:minus[x]+plus[x] for x in minus}
    diff = {x:minus[x]-gamma*plus[x] for x in minus}
    assert diff[2,0]==0
    ones = {x:Q(1) for x in m.rows}
    assert m.dual(ones,Q(0))==1 and m.dual({x:-1 for x in ones},Q(0))==-1
    print('cap',cap,'physical rows',len(minus),'normalized LP states',len(m.rows),'gamma',float(gamma))
    missing = {x:Q(1) for x in m.rows if x[0]==cap+1}
    show(m,'unresolved degree-tail mass',missing)
    young = {x:Q(x==(2,0)) for x in m.rows}
    show(m,'virgin mass',young)
    show(m,'retained lower f (micro)',minus)
    show(m,'retained upper f (micro)',plus)
    total_lo,_,t_sum,_ = show(m,'retained sum (micro)',total)
    e_lo,e_hi,t_lo,t_hi = show(m,'retained lower - gamma*upper (micro)',diff)
    # Analytic initial residence conditional on a physical promoter bin:
    # true parent degree can only increase thereafter. No frozen hazard.
    for birth_d in (2,4,8,18):
        conditional = LO/(LO+Q(3,2)+Q(1,birth_d))
        prob = birth_floor(birth_d)
        print('birth degree',birth_d,'Palm mass floor',f'{float(prob):.9g}',
              'conditional virgin mass floor',f'{float(conditional):.9g}',flush=True)
    models = [(Q(0),m)]+[(birth_floor(d),Age(cap,cap,d)) for d in range(3,9)]
    total_lo = max(total_lo,birth_mix(models,total,t_sum))
    e_lo = max(e_lo,birth_mix(models,diff,t_lo))
    e_hi = min(e_hi,-birth_mix(models,{x:-v for x,v in diff.items()},t_hi))
    print('after exact-degree birth strata: retained sum lower',ends(total_lo,total_lo),
          'retained contrast',ends(e_lo,e_hi),flush=True)
    # Use pre-existing denominator ceilings only; no denominator search.
    c0minus, c0plus = Q('1.458909469275165'), Q('1.518769126087492')
    width = max((s-MINUS)*c0minus+(PLUS-s)*c0plus for s in (LO,HI))*UNIT
    rem = width-total_lo
    assert rem>=0
    print('total all-class budget (micro)',ends(0,width),'remaining all-class cap (micro)',ends(0,rem))
    print('full paired difference cone-only enclosure (micro)',ends(e_lo-gamma*rem,e_hi+rem))
    band = tuple(sum(Q(side[k]) for side in BANDS)*UNIT for k in range(5))
    missing_floor = {x:band[x[1]] for x in missing if x[1] < min(cap,5)}
    tail_lo,_,t_tail,_ = show(m,'missing f-sum floor payoff (micro)',missing_floor)
    tail_lo = max(tail_lo,birth_mix(models,missing_floor,t_tail))
    print('remaining all-class enclosure using missing f floors (micro)',ends(tail_lo,rem))
    # Buy only a previously omitted degree class. The control generator is
    # unchanged; extend the descriptor one degree to read its exact rows.
    nxt = cap+1
    rm,rp = inner_rows(nxt)
    extra = {x:rm[x]-gamma*rp[x] for x in rm if x[0]==nxt}
    young_plus_one = Age(nxt,nxt)
    v,w,_,_ = show(young_plus_one,f'missing D={nxt} paired contrast only (micro)',extra)
    extra_sum = {x:rm[x]+rp[x] for x in extra}
    cost,*_ = show(young_plus_one,f'missing D={nxt} f-sum (micro)',extra_sum)
    joint_lo,joint_hi,*_ = show(young_plus_one,f'retained through D={nxt} joint contrast (micro)',diff|extra)
    joint_lo,joint_hi = max(joint_lo,e_lo+v),min(joint_hi,e_hi+w)
    print('combined retained contrast (micro)',ends(joint_lo,joint_hi),
          'remaining all-class cap (micro)',ends(0,rem-cost))
    print(f'D={nxt} contrast row hull (micro)',ends(min(extra.values()),max(extra.values())))
    # Above D=nxt charge only the two all-L floors, leaving their excess
    # (along with the other motifs) in the nonnegative residual budget.
    suffix = {x:(UNIT*Q(BANDS[0][x[1]]),UNIT*Q(BANDS[1][x[1]]))
              for x in young_plus_one.rows if x[0]==nxt+1 and x[1]<5}
    suffix_sum = {x:a+b for x,(a,b) in suffix.items()}
    suffix_diff = {x:a-gamma*b for x,(a,b) in suffix.items()}
    paid = total|extra_sum|suffix_sum
    contrast = diff|extra|suffix_diff
    paid_lo,*_ = show(young_plus_one,'disjoint combined f-sum (micro)',paid)
    con_lo,con_hi,*_ = show(young_plus_one,'disjoint combined f-contrast (micro)',contrast)
    residual = width-max(total_lo+cost,paid_lo)
    assert residual>=0
    print('after disjoint floor pairing, residual cap (micro)',ends(0,residual))
    print('after disjoint floor pairing, full contrast (micro)',
          ends(con_lo-gamma*residual,con_hi+residual))
    print('later k >= 5 absolute mass cap from birth clock',ends(Q(0),(1+LO)**-5))


if __name__=='__main__':
    main()
