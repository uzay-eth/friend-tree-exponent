"""Pair the opposite sharp-tail prices on each port of the old core tree."""

from fractions import Fraction as Q
from functools import cache
from itertools import combinations_with_replacement
from pathlib import Path
from runpy import run_path

HERE=Path(__file__).resolve().parent
LOCAL=run_path(str(HERE.parent/'ledger.py'))
X=LOCAL['COEFF']
Y=LOCAL['PRICE']
CUT=32
ceil=LOCAL['ceil']


def chains(groups,gain,u3,fav3,high):
    bare=ceil(Y*max(0,u3-fav3))
    paid=2*gain*Y//1
    def bill(data,l,h):
        if data is None:
            return max(h*ceil(Y*u3),X*ceil(Y*high/(CUT+1)))
        w,sl,lines=data
        return max(0,(h-sl)*w,*((sh-l)*cap for sh,cap in lines))
    ends=[]
    for d,(fav,typed) in groups.items():
        ends.append((Q(1,d),(ceil(Y*max(0,u3-fav)),min(lo for lo,_,_ in typed),
                            {(hi,cap) for _,hi,cap in typed}),int(d==2)))
    ends.extend([(Q(0),None,0),(Q(1,CUT),None,0)])
    result=[0,0,0]
    for (u,one,t),(v,two,z) in combinations_with_replacement(ends,2):
        p,q=X*u/6,X*v/6
        for seq in ((p+q,), (X//12+p,X//12+q), (X//12+p,Q(X,6),X//12+q)):
            boxes=[(a//1,ceil(a)) for a in seq]
            cost=bill(one,*boxes[0])+bill(two,*boxes[-1])
            cost+=bare*sum(max(h-a,b-l) for (l,h),(a,b) in zip(boxes,boxes[1:]))
            cost-=paid*sum(l for l,h in boxes)
            result[t+z]=max(result[t+z],cost)
    return [Q(v,X*Y) for v in result]


def conditions(middle,fav,host,U,y,high_constant=None):
    u2=U(2,0)
    gain=middle-u2
    leaf_gain=host[2]-u2
    high=(Q(17,16)/y if high_constant is None else high_constant)/CUT
    if gain<=high:
        return None
    groups={}
    cons={}
    def add(row,rhs):
        cons[row]=max(cons.get(row,-999),rhs)
    def profile(d,c,r):
        l=d-c
        s=(l+r)/(d*(d+1))
        t=Q(1,2*d)
        reward=d*s*(host.get(d+1,middle)-u2)+l*t*leaf_gain
        cost=l*max(0,t-s)*max(0,u2-fav.get(d+1,0))
        if l:
            cost+=l*max(0,s-t)*max(0,U(d+1,l-1)-fav[2])
        return s,reward-cost
    for d in range(2,CUT):
        if profile(d,0,Q(0))[1]<=0:
            return None
        typed=set()
        for c in range(1,d+1):
            if c==d==2:
                continue
            own=U(d+1,d-c)
            for r in (Q(1,2),Q(c,2)):
                s=(d-c+r)/(d*(d+1))
                typed.add((X*s//1,ceil(X*s),ceil(Y*max(0,own-fav[3]))))
            for k in range(c+1):
                direct=c-k
                last=6 if d<=4 else 1
                for es in combinations_with_replacement(range(2,last+2),direct):
                    known=[e for e in es if e<=last]
                    r=Q(k,2)+sum(Q(1,e) for e in known)
                    unknown=direct-len(known)
                    for rec in (r,r+Q(unknown,last+1)):
                        s,credit=profile(d,c,rec)
                        exposures=sorted(
                            [max(0,s-Q(1,d*e*(e+1)))*max(0,own-fav.get(e+1,0)) for e in known]
                            +[s*own]*unknown,reverse=True)
                        debit=0
                        for j in range(direct+1):
                            add((c-2,j-1,0 if d==2 else -k,-k if d==2 else 0),debit-credit)
                            if j<direct:
                                debit+=exposures[j]
        groups[d]=(fav.get(d+1,0),typed)
    bills=chains(groups,host[3]-u2,U(3,0),fav[3],high)
    for i,bill in enumerate(bills):
        add((0,1,2-i,i),bill)
    sparse=min(0,(Q(CUT,CUT+1)*(gain-high)+(leaf_gain-u2)/2)/2)
    rich=Q(CUT,CUT+1)*(gain-high)/2+leaf_gain/4-u2/(2*(CUT+1))
    add((0,1,0,0),high/(CUT+1))
    add((-2,-1,0,0),-rich)
    add((14,-1,-16,0),-sparse)
    add((1,0,-1,0),0)
    add((0,0,-1,1),0)
    return cons


def certify(cons,tariff):
    if cons is None or not all(x>0 for x in tariff):
        return False
    return all(sum(i*w for i,w in zip(row,tariff))>rhs for row,rhs in cons.items())
