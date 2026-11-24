"""Exact ledger replays using the Poisson-stock serial lower bound."""
from fractions import Fraction as Q
from hashlib import sha1
from math import factorial
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
band = run_path(str(HERE.parent / 'band/check.py'))
M = Q(36787, 100000)


if __name__ == '__main__':
    from argparse import ArgumentParser
    p = ArgumentParser(description=__doc__)
    p.add_argument('--finite', action='store_true')
    p.add_argument('--obstruction', action='store_true')
    p.add_argument('--start', type=int, default=20)
    p.add_argument('--stop', type=int, default=39)
    args = p.parse_args()
    assert sum((-1)**i*Q(1,factorial(i)) for i in range(10)) > M
    for file, expected in (
        ('band/check.py', 'e2aa75a98bf62981a91e33c3b3664e2f55d248c6'),
        ('check.py', 'e39f8bfdf66f19e34a59247fa3649ccaf884358b'),
    ):
        data = (HERE.parent / file).read_bytes()
        assert sha1(f'blob {len(data)}\0'.encode()+data).hexdigest() == expected
    if args.obstruction:
        x = Q(1,10)
        ok,where = band['ledger'](20,M,x,x)
        assert not ok and where == ('site',2,1,0,Q(0))
        caps = band['tips'](20,x,x)
        U = caps[2,0]
        assert (M-U)*Q(7,12)-U/12 < Q(8,1000)+max(0,caps[3,1]/6-Q(4,1000))
        print(f'LIMIT b=20, x=1/10: degree-two terminal misses reserve with U(2,0)={U}')
    elif args.finite:
        assert 20 <= args.start < args.stop <= 39
        for b in range(args.start,args.stop):
            seam = Q(5,b+1)
            cuts = [Q(1,8), *(Q(i,200) for i in range(26,48) if Q(i,200)<seam), seam]
            for x,y in zip(cuts,cuts[1:]):
                ok,where = band['ledger'](b,M,x,y)
                assert ok,(b,x,y,where)
            print(f'PASS b={b} over [1/8,{cuts[-1]}]',flush=True)
    else:
        scaled = run_path(str(HERE / 'check.py'))
        assert 2*3**32>5*2**33
        assert scaled['chains'](5,M) == (True,735)
        print('PASS symbolic ledger at y>=5')
