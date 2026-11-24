"""Cold replay of the birth-scale arithmetic and uniform hub bases."""
from fractions import Fraction as Q
from hashlib import sha1
from itertools import combinations_with_replacement, product
from pathlib import Path
from runpy import run_path


AGE = Path(__file__).resolve().parents[1] / 'open/age-order'
HOSTILE = AGE / 'global-transport/hostile'


def authenticate():
    for folder, readme, checker in (
        ('global-transport/birth-window', '2f5e400501497715b464d2505c4c7c75cae2ef7d',
         'ef023018466ad515c16d415caeaeafdc170f8cf2'),
        ('hub-ranking', '5086e7ee4315e616d2b55f1113d2715ecc24be90',
         'e940014900eac1b8cd777a7598ffa31d3d4658bc'),
    ):
        for name, expected in (('README.md', readme), ('check.py', checker)):
            raw = (AGE / folder / name).read_bytes()
            assert sha1(f'blob {len(raw)}\0'.encode() + raw).hexdigest() == expected
    print('PASS four reviewed writer blobs')


def main():
    authenticate()
    endpoint = run_path(str(HOSTILE / 'mixed-stock/check.py'))
    G, U, L, P = (endpoint[k] for k in ('G', 'U', 'L', 'profile'))
    a, b, f, q = map(Q, ('.00634', '.017', '.00356', '.00625'))
    eta = Q(1, 2**27)
    delta = Q(1, 200000)
    D = 2**32
    stages = 19 + 20 + 20*32 + 1
    assert stages == 680 < 700 and Q(64, 64-31) < 2
    size = Q(1, 2**45)  # (2e^-16)^3 < 2^-45, e>2
    late = Q(1, 2**(1024-700))
    assert Q(19, 20)*128 > 121
    assert size+late+Q(40*D, 2**121)+Q(8, D) < eta
    assert size+Q(4*D, 2**64)+Q(8, D) < eta
    # e<4, and the true retained stock gives the extra D/4 multiplier.
    assert 5-131400 < 30-2*(2**16+128)
    assert 5-131400 < 30-2*128
    print('PASS 680-decision cap, Gamma size tail and both age-to-mass conversions')

    high = min(C+a*(c-2)-b-c*(max(0,s*u-b)+fee)
               for d in range(10,113) for c in range(2,d+1)
               for fee,e in ((0,0),(0,2),(f,4),(f,2),(q,2))
               for s,C,u in (P(d,(e,)*c),))
    g, p = (endpoint[k] for k in ('g', 'p'))
    bracket = g*Q(113,114)+G(2)/2-p/2
    assert bracket < 0 and Q(1,228) < b and G(2) > Q(2,114)
    high_tail = min(bracket/2+(a-f)*Q(113,2)-2*a-b,
                    g*Q(113,228)-2*q-b,
                    g*Q(9,11)+G(2)*Q(9,20)-Q(1,22)-a-q)

    def bills(d, c, ds, r):
        s,C,u=P(d,tuple(e or r for e in ds))
        charge = Q('.021') if d == 2 else Q('.0037')
        for e in ds:
            t = Q(1,d*e*(e+1)) if e else 0
            h = L(e+1) if e else 0
            direct=max(0,s-t)*max(0,u-h)-b
            C -= max(0, direct, charge if e == 2 else 0)
        return C+a*(c-2)-b+(f*ds.count(0) if (d,c) in ((3,2),(3,3),(4,4)) else 0)

    low = min((bills(d,c,ds,r),d,c,ds,r) for d in range(2,10)
              for c in range(1,d+1) if (d,c)!=(2,2)
              for ds in combinations_with_replacement((0,*range(2,10)),c) for r in (0,10))
    assert min(high,high_tail,low[0]) > delta > 4*eta
    print(f'PASS simultaneous retained budgets: minimum {low[0]}; after 4 eta > {delta-4*eta}')

    h = Q('.012')
    for d in range(10,42):
        for l in range(d):
            assert Q(d+l,2*d*(d+1))*U(d+1,l)<h
    assert Q(1,86)<h
    u=U(3,0)
    assert u>L(3) and G(3)>eta
    highs=[(Q(1,d),L(d+1),None,None,q) for d in range(10,31)]
    highs += [(x,0,None,None,q) for x in (Q(0),Q(1,31))]
    lows=[]
    for d in range(2,10):
        for c in range(1,d+1):
            if (d,c)!=(2,2):
                for r in (0,2):
                    s,_,v=P(d,(2,)+(r,)*(c-1))
                    lows.append((Q(1,d),L(d+1),s,v,Q('.021') if d==2 else Q('.0037')))

    def port(s, endpoint):
        _,ell,t,v,fee=endpoint
        slope=max(0,u-ell)
        if t is None:
            return max(h,s*slope)-fee
        return max(0,(s-t)*slope,(t-s)*max(0,v-L(3)))-fee

    margins=[]
    for left,right in ((highs,highs),(lows,highs),(lows,lows)):
        for length in (1,2,3):
            cost=[]
            for X,Y in product(left,right):
                x,y=X[0],Y[0]
                ss= ((x+y)/6,) if length==1 else ((Q(1,2)+x)/6,)+(Q(1,6),)*(length-2)+((Q(1,2)+y)/6,)
                cost.append(port(ss[0],X)+port(ss[-1],Y)
                            +(u-L(3))*sum(abs(i-j) for i,j in zip(ss,ss[1:]))-2*G(3)*sum(ss))
            margins.append(b-max(cost))
    assert min(margins) == Q(172282757,2400000000000) > delta > 6*eta
    print(f'PASS nine chain classes: minimum margin {min(margins)}; after 6 eta > {delta-6*eta}')

    old = endpoint['D']['H']
    floor = a*2+b
    assert floor == Q(371,12500)
    for d, tick in enumerate((18940948,20181908,22577249),1):
        t=old['spider'](d,1)
        score=old['score'](t)[0]*sum(old['profile'](t)[1])
        assert score>Q(tick,10**8)>floor
    assert L(2)>Q(1,4) and G(2)/2==Q(189409485129,2*10**12)>floor
    assert g*Q(4,5)>4*eta and Q(1,200)>4*eta and G(2)>eta
    early=old['EARLY']
    founder=Q(3,4)-1+early['silent_bound'](1,4)
    assert founder==Q(127,1220)>eta
    # Different earlier stopping sizes give different offsets but the same
    # conditional kernel at each size, and hence the same reciprocal identity.
    for m,i,j in ((2,3,4),(2,3,19),(4,7,12),(32,33,88)):
        assert 2<=m<=i-1<j-1
        assert sum(floor/(n*(n+1)) for n in range(i-1,j-1))==floor*(Q(1,i-1)-Q(1,j-1))
    print('PASS exact seed/star floor, founder seam and earlier-size reciprocal sum')


if __name__ == '__main__':
    main()
