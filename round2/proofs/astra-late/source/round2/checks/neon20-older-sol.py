#!/usr/bin/env python3
"""Independent exact D20 hV older-arm replay and same-witness ledger."""
import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
from itertools import permutations, product
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1] / 'occupation-neon20'
PINS = {
    'renewal/older/README.md': 'd3e1f6cc48d9828127254fee307fce3f40915849f3cb011c39be4b0f1d61f527',
    'renewal/older/replay.py': 'a0a2bdac60436a12156cc3b35c8a01526587b28c217affc839897c47ab15dc25',
    'rows/extract.py': 'f9352f511c7eb38f94aee4cf6bbce4a3286ae0a4dcb13f577578773ab4dbf1fe',
    'rows/lower.json': 'b59a2c20b31e7b09f88df436e94e75921d22ac424a759a58c3caf118f7692d4d',
    'rows/upper.json': 'eeb67cf94792f0a6372f72b02214fa07656a166a65935997d06703a8d9ed85ad',
    'renewal/replay.py': 'c72fc793e1aaff19d86d22fc72694e301dfa163b0d311d484217e97cb8bf4257',
}
sys.path.insert(0,str(ROOT/'rows'))
from extract import Cert, Rows

S, H = 10**24, Q('0.565768')
LO = Q('0.565680')
J = (2,0)   # true total degree, core children; parent slot always present
EDGES = ((0,1),(1,2),(0,3))  # A-B, B-C, A-U
ORDERS = ((1,2,3),(1,3,2),(3,1,2))
BASE = {'lower': (227220866275491127308960,5718952004892475098000),
        'upper': (222182631325898847818400,2470851317190905657040)}
DECIMALS = {'lower': ('4.03428774415e-8','2.76609076108e-8','4.29685004568e-8'),
            'upper': ('1.75177464202e-7','5.90407013113e-8','9.10074725771e-8')}


def census(x):
    d,k=x
    return d-k-1,k,1


def change(x,i,y):
    return x[:i]+(y,)+x[i+1:]


def bound(b,d,options):
    # Exact linear residual in rate coordinates; locate its last nonnegative
    # integer tick directly, without the production fractional greedy routine.
    def resid(w):
        return b-d*w+sum(cap*min(0,nxt-w) for cap,nxt in options)
    w=None
    for breakpoint in sorted({z for _,z in options}):
        if resid(breakpoint)<0:
            break
    else:
        breakpoint=None
    test = [(cap,nxt) for cap,nxt in options if breakpoint is None or nxt<breakpoint]
    w=(b+sum(cap*nxt for cap,nxt in test))/(d+sum(cap for cap,_ in test))
    assert resid(w)==0
    tick=w//1
    assert tick>=0 and resid(tick)>=0>resid(tick+1)
    return tick


def fixture():
    parent = {b:a for a,b in EDGES}
    assert tuple(order for order in permutations((1,2,3)) if order.index(parent[2])<order.index(2))==ORDERS
    # Distinguishable actual vertices: two apex siblings may exchange roles
    # only if the chosen B has its own designated C. Birth chronology does
    # not erase a pre-existing sibling when introducing a new registration.
    births=((0,1),(0,2),(1,3),(2,4),(0,5))
    for history in permutations(births):
        age={child:i for i,(_,child) in enumerate(history)}
        if any(p and age[p]>age[ch] for p,ch in history):
            continue
        physical={(b,c,u) for a,b in births if a==0 for p,c in births if p==b
                  for a,u in births if a==0 and u!=b}
        by_order=[]
        for order in ORDERS:
            by_order += [tags for tags in physical if tuple(sorted((1,2,3),key=lambda j:age[tags[j-1]]))==order]
        assert Counter(by_order)==Counter({x:1 for x in physical})


def solve(c,rows):
    tests=Counter()
    mutation=Counter()

    @cache
    def reward(x):
        r=tuple(map(census,x[:3]))
        a,ports=rows.affine('hV',r)
        floor=a+sum(n*min(0,s/2) for _,_,n,s in ports)
        assert floor>=0
        # Independently derive the apex child-port slope from the physical
        # hV gain difference and the eligibility for hV -> Yb transfer.
        ports_a=[(n,s) for j,pa,n,s in ports if j==0 and not pa]
        assert ports_a and ports_a[0][0]==x[0][1]-1
        beta=ports_a[0][1]
        d,k=x[0]
        raw=rows.value('hV',change(r,0,census((d+1,k))))-rows.value('hV',r)
        routed=r[0] in c.ybu and r[1] in c.ybw
        assert beta==(0 if routed else c.sign*raw/c.S)
        send=Q(1,x[3][0])
        actual=beta*send-min(0,beta/2)
        assert actual>=0
        if actual>0:
            tests['positive']+=1
            tests['positive slope' if beta>0 else 'negative slope']+=1
            mutation['recipient']+=actual != beta/d-min(0,beta/2)
            mutation['virgin']+=actual != beta/2-min(0,beta/2)
            mutation['flipped sign']+=actual != -beta*send-min(0,-beta/2)
        return S*actual

    pieces=[]
    for order in ORDERS:
        # Store coordinates in role order throughout, not birth order.
        @cache
        def v(x):
            m=sum(z is not None for z in x)
            for i,z in enumerate(x):
                if z is not None and (z[0]>(16 if i==0 else 6) or z[1]>(3 if i==0 else 2)):
                    return 0
            size_lb=m+sum(census(z)[0] for z in x if z is not None)
            # Threshold is on payment only: smaller reachable tags still advance.
            b=reward(x) if m==4 and size_lb>16 else Q(0)
            if m==4 and size_lb==16:
                tests['seam states']+=1
                mutation['size>=16']+=reward(x)>0
            d=H
            options=[]
            for i,z in enumerate(x):
                if z is None:
                    continue
                degree,k=z
                l=degree-k-1
                assert l>=0
                peers=tuple(j if i==p else p for p,j in EDGES if i in (p,j) and x[j if i==p else p] is not None)
                outside=k+1-len(peers)
                assert outside>=0
                plus=v(change(x,i,(degree+1,k)))
                inf=Q(l)+sum(Q(1,x[j][0]) for j in peers)
                options.append((Q(outside,2),plus))
                b+=inf*plus
                d+=inf
                if l:
                    after=change(x,i,(degree,k+1))
                    rate=Q(l,degree)
                    b+=rate*v(after)
                    d+=rate
                    if m<4 and (i,order[m-1]) in EDGES:
                        b+=rate*v(change(after,order[m-1],J))
                        tests['extensions']+=1
            w=bound(b,d,options)
            if m==4 and size_lb==17:
                tests['explicit corners']+=1
                assert all((d+sum(cs))*w<=b+sum(a*z for a,(_,z) in zip(cs,options))
                           for cs in product(*((0,a) for a,_ in options)))
            return w

        ticks=v((J,None,None,None))
        assert ticks>0
        print('order',''.join('ABCU'[i] for i in order),'newborn ticks',ticks,
              'states incl. exits',v.cache_info().currsize,flush=True)
        pieces.append(LO*Q(ticks,S))
        v.cache_clear()
    assert all(mutation[k]>0 for k in ('recipient','virgin','flipped sign','size>=16')),mutation
    print('support checks',dict(tests),'planted mismatches',dict(mutation),flush=True)
    return pieces


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('side',choices=BASE)
    args=p.parse_args()
    for f,h in PINS.items():
        assert sha256((ROOT/f).read_bytes()).hexdigest()==h,f
    fixture()
    s=args.side
    path=Path('b3f20_'+('lo' if s=='lower' else 'hi')+'.bin')
    initial=path.stat()
    assert initial.st_ino==(3690464 if s=='lower' else 3856843)
    c=Cert(path,s,verify=False)  # the unchanged baseline payloads already hashed independently
    pieces=solve(c,Rows(c))
    for x,h in zip(pieces,DECIMALS[s]):
        assert abs(x-Q(h))<Q(1,10**18),(x,h)
    extra=sum(pieces)
    base=Q(sum(BASE[s]),10**29)
    den=Q(json.loads((ROOT/'rows'/f'{s}.json').read_bytes())['recharge']['ceiling'])
    total=base+extra
    z=c.lam+c.sign*total/den
    outward=c.sign*(c.sign*z*10**15//1)
    print(s,'base',base,'extra',extra,'total',total,'denominator',den,'outward ticks',outward,flush=True)
    c.check_unchanged()


if __name__=='__main__':
    if not __debug__:
        raise SystemExit('assertions required')
    main()
