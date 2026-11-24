"""Signed comparison with the physical return and serial host still marked."""

from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

HERE = Path(__file__).resolve().parent
OLD = run_path(str(HERE.parent / 'ledger.py'))
CUT = OLD['CUT']


def prepare(middle, fav, host, cap, y):
    assert Q(1,2) <= y <= Q(17,16)**CUT
    u2 = cap(2,y)
    gain = middle-u2
    high = Q(17,16)/(CUT*y)
    if gain <= high:
        return None
    leaf_gain = host[2]-u2

    def profile(d,c,receive,own):
        l = d-c
        s = (l+receive)/(d*(d+1))
        t = Q(1,2*d)
        reward = d*s*(host.get(d+1,middle)-u2) + l*t*leaf_gain
        leaf_cost = l*max(0,t-s)*max(0,u2-fav.get(d+1,0))
        core_cost = l*max(0,s-t)*max(0,own-fav[2])
        return s,reward-leaf_cost-core_cost

    types = set()
    sites = []
    for d in range(2,CUT):
        own = cap(d+1,y)
        if profile(d,0,Q(0),own)[1] <= 0:
            return None
        for c in range(1,d+1):
            if d==c==2:
                continue
            for k in range(c+1):
                for receive in {Q(k,2),Q(c,2)}:
                    s,reward = profile(d,c,receive,own)
                    sites.append((c,k,d==2,reward,s*own))
            for receive in {Q(1,2),Q(c,2)}:
                s,_ = profile(d,c,receive,own)
                types.add((Q(1,d),s,own,d==2))
    for reciprocal in (Q(0),Q(1,CUT)):
        types.add((reciprocal,None,high/(CUT+1),False))
    bills = OLD['chain_bills'](types,host[3]-u2,cap(3,y))
    sparse = min(0,(Q(CUT,CUT+1)*(gain-high)+(leaf_gain-u2)/2)/2)
    rich = Q(CUT,CUT+1)*(gain-high)/2 + leaf_gain/4-u2/(2*(CUT+1))
    return sites,bills,sparse,rich,high
