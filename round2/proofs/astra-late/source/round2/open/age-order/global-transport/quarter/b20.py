"""Separate conditional b=20 extension; does not alter the frozen ledger."""
from fractions import Fraction as Q
from hashlib import sha1
from pathlib import Path
from runpy import run_path


path = Path(__file__).with_name('check.py')
data = path.read_bytes()
assert sha1(f'blob {len(data)}\0'.encode()+data).hexdigest() == '5dab19ad0224b96b63bcaf0a8915c623609609c7'
A = run_path(str(path))


def cover(lo, hi, remaining=10):
    U = A['tips'](20,lo,hi)
    if A['ledger'](20,U)[0]:
        return 1
    site,chain = A['improvement'](20,U,A['favourable'](20,hi,24))
    for f,ft in ((4,8),(5,9),(5,10),(4,11)):
        fee,terminal = Q(f,1000),Q(ft,1000)
        failure = site(Q('.004'),Q('.004'),fee,terminal,Q('.002')) or chain(Q('.004'),fee,terminal)
        if failure is None:
            return 1
    assert remaining > 0,(lo,hi,failure)
    mid = (lo+hi)/2
    return cover(lo,mid,remaining-1)+cover(mid,hi,remaining-1)


if __name__ == '__main__':
    lo,hi = Q(1,10),Q(1,3)
    assert 20*hi/(1-hi) >= 6
    n = cover(lo,hi)
    print(f'PASS b=20 on [{lo},{hi}]: {n} interval boxes; Poisson cutoff above {hi}')
