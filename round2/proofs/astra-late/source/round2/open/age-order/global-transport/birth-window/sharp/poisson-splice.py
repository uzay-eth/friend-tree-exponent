"""Replay the large birth-scale and finite-label splice with the Poisson tail."""
from argparse import ArgumentParser
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
from runpy import run_path


ROOT = Path(__file__).resolve().parents[2]
QUARTER = run_path(str(ROOT / 'quarter/check.py'))
THIRD = ROOT / 'annealed/extension/third-physical'
BAND = run_path(str(THIRD / 'band/check.py'))
SCALED = run_path(str(THIRD / 'scaled/check.py'))
M = Q(36787,100000)


def finite(b, lo, hi):
    # Only the initial stock-one middle gets this extra first-step credit.
    # Keep the same improvement API with the common lower M: raising M also
    # raises the favourable payoff after promotion, which starts at stock one
    # at a *later* size. Compute that payoff before changing the serial constant.
    if lo>=Q(3,20) and BAND['ledger'](b,M,lo,hi)[0]:
        return True
    depth = 32 if b<24 else 14
    U = QUARTER['tips'](b,lo,hi,depth)
    L = QUARTER['favourable'](b,hi,depth)
    helper = QUARTER['improvement']
    helper.__globals__['M'] = M*(1+(1-hi)/(2*(b+1)))
    try:
        site, chain = helper(b,U,L)
        tariff = tuple(Q(v,1000) for v in (4,4,5,9,2))
        return not (site(*tariff) or chain(*tariff[1:4]))
    finally:
        helper.__globals__['M'] = M


def main():
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--start',type=int,default=20)
    parser.add_argument('--stop',type=int,default=49)
    args=parser.parse_args()
    assert 20<=args.start<args.stop<=49
    assert sum(Q((-1)**j,factorial(j)) for j in range(10)) > M
    assert 5 < Q(3,2)**32 and Q(33,32)*Q(3,2)<2
    assert SCALED['chains'](5,M) == (True,735)
    print('PASS all degrees at (b+1)x >= 5',flush=True)
    assert BAND['ledger'](20,M,Q(3,20),Q(1,2))[0]
    print('PASS b=20 on [3/20,1/2] with the unboosted lower',flush=True)
    def cover(b,lo,hi,depth):
        if (b+1)*lo>=5 or finite(b,lo,min(hi,Q(5,b+1))):
            return
        assert depth,(b,lo,hi)
        mid=(lo+hi)/2
        cover(b,lo,mid,depth-1)
        cover(b,mid,hi,depth-1)

    for b in range(args.start,args.stop):
        cuts = (Q(1,10),Q(11,100),Q(1,8),Q(3,20),Q(1,2)) if b<24 else tuple(
            [Q(1,10)] + [Q(1,10)+Q(i,100) for i in range(1,6)] + [Q(1,2)])
        for lo,hi in zip(cuts,cuts[1:]):
            cover(b,lo,hi,5)
        print(f'PASS b={b} on [1/10,1/2]',flush=True)


if __name__ == '__main__':
    main()
