"""Signed core-tree budget, with exact allowances for chains and all larger degrees."""

from fractions import Fraction as Q
from functools import cache
from itertools import combinations_with_replacement

CUT = 32
COEFF = 600_000_000
PRICE = 10**9


def ceil(x):
    return -((-x)//1)


def profile(d, c, receive, gain, u2, own):
    leaves = d - c
    s = (leaves + receive)/(d*(d+1))
    leaf = Q(1,2*d)
    reward = gain*(d*s+leaves*leaf)
    cost = leaves*(max(0,leaf-s)*u2 + max(0,s-leaf)*own)
    return s, reward-cost


def chain_bills(types, gain, u3):
    """Upper gross bill by number of degree-two terminal endpoints."""
    # Use rectangular integer enclosures; the unit is COEFF*PRICE.
    marked = []
    for reciprocal, s, own, terminal in types:
        x = COEFF * reciprocal/6
        bounds = None if s is None else (COEFF*s//1, ceil(COEFF*s))
        marked.append((x//1, ceil(x), bounds, ceil(PRICE*own), terminal))
    bare_price = ceil(PRICE*u3)
    bare_credit = PRICE*2*gain//1

    @cache
    def at_boundary(lo, hi, end):
        _, _, retained, cap, _ = end
        if retained is None:
            # 'cap' is the entire own outward cost s*U on this high port.
            return max(hi*bare_price, COEFF*cap)
        sl, sh = retained
        return max(0,(hi-sl)*bare_price,(sh-lo)*cap)

    answer = [0,0,0]
    for x,y in combinations_with_replacement(marked,2):
        terminal_count = x[-1]+y[-1]
        for length in (1,2,3):
            if length == 1:
                chain = [(x[0]+y[0],x[1]+y[1])]
            else:
                chain = [(COEFF//12+x[0],COEFF//12+x[1])]
                if length == 3:
                    chain.append((COEFF//6,COEFF//6))
                chain.append((COEFF//12+y[0],COEFF//12+y[1]))
            bill = at_boundary(*chain[0], x)+at_boundary(*chain[-1],y)
            bill += bare_price*sum(max(rhi-llo,lhi-rlo)
                                   for (llo,lhi),(rlo,rhi) in zip(chain,chain[1:]))
            bill -= bare_credit*sum(lo for lo,_ in chain)
            answer[terminal_count] = max(answer[terminal_count],bill)
    return [Q(v,COEFF*PRICE) for v in answer]


def prepare(middle, cap, y):
    assert Q(17,16)**CUT >= y
    high = Q(17,16)/(CUT*y)
    gain = middle-cap(2,y)
    if gain <= 0:
        return None
    types = set()
    sites = []
    for d in range(2,CUT):
        own = cap(d+1,y)
        if profile(d,0,Q(0),gain,cap(2,y),own)[1] <= 0:
            return None
        for c in range(1,d+1):
            if d==c==2:
                continue
            for k in range(c+1):
                for receive in {Q(k,2),Q(c,2)}:
                    s,reward=profile(d,c,receive,gain,cap(2,y),own)
                    sites.append((c,k,d==2,reward,s*own))
            for receive in {Q(1,2),Q(c,2)}:
                s,_=profile(d,c,receive,gain,cap(2,y),own)
                types.add((Q(1,d),s,own,d==2))
    for reciprocal in (Q(0),Q(1,CUT)):
        types.add((reciprocal,None,high/(CUT+1),False))
    bills = chain_bills(types,gain,cap(3,y))
    sparse = min(0,(Q(CUT,CUT+1)*(gain-high)+(gain-cap(2,y))/2)/2)
    rich = Q(CUT,CUT+1)*(gain-high)/2+gain/4-cap(2,y)/(2*(CUT+1))
    if gain<=high:
        return None
    return sites,bills,sparse,rich,high


def certify(prepared, tariffs):
    if prepared is None:
        return False
    sites,bills,sparse,rich,high=prepared
    a,base,fee,terminal=tariffs
    assert 0<fee<=a and base>0 and terminal>=fee
    if high/(CUT+1)>base or min(rich,(a-fee)*CUT/2+sparse)<=2*a+base:
        return False
    if any(bill>base+(2-i)*fee+i*terminal for i,bill in enumerate(bills)):
        return False
    for c,k,is_terminal,reward,own in sites:
        if reward+a*(c-2)-base < k*(terminal if is_terminal else fee)+(c-k)*max(0,own-base):
            return False
    return True
