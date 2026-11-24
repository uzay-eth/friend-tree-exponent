"""Exact transfer budgets for the small-power and below-one windows."""
from fractions import Fraction as Q
from math import factorial
from pathlib import Path
from runpy import run_path


HERE = Path(__file__).resolve().parent
BIRTH = HERE.parents[1] / 'global-transport' / 'birth-window' / 'check.py'
MEAN = HERE.parent / 'mean-ledger' / 'check.py'
FRACTIONAL = HERE.parent / 'fractional-middle' / 'check.py'


def main():
    birth = run_path(str(BIRTH))
    birth['check']()
    mean = run_path(str(MEAN))['Ledger'](1, Q('.56319'))
    mean.first()

    # A whole capped return, including all its sender-degree changes, ends
    # at intensity >= 1/31. Do not charge 31 again at each sender degree.
    retry_time = 40 + 20*31 + 1
    loss = Q(40*2**32, 2**121) + Q(8, 2**32)
    base_loss = Q(4*2**32, 2**64) + Q(8, 2**32)
    assert max(loss, base_loss) < Q(1,2**27)
    small_p = Q(1,10**9)
    eta = Q(1,2**27) + (retry_time+128+1)*small_p
    assert 4*eta < birth['delta'] and 6*eta < birth['delta']
    assert eta < min(Q(1,400), Q(127,1220), birth['G'](2), birth['G'](3))
    print(f'PASS small powers through {small_p}: endpoint difference loss < {float(eta):.12g}')

    # e > 8/3, so 1/e < 3/8. In the upper bound below use p >= 1-q.
    q = Q(1,15000)
    assert (1-q)*Q(65,24) > Q(8,3)
    # Sharpen the upper slope using e >= 65/24 and p >= 1-q.
    upper = Q(3,8)
    serial = 6*Q(3,8)+upper
    edge = 2*Q(3,8)+upper
    retained = Q(3,2)*(serial+edge)+edge
    chain = serial + Q(4,3)*edge
    assert retained*q < Q('.0005')
    assert chain*q < Q('.000291477')
    assert retained*q < Q(1,40)
    assert (Q(64,65)*mean.g+(mean.g-mean.U(2))/2)/2-Q(1,40) > retained*q
    assert Q(32)*Q('.0045')-Q(1,40) > retained*q
    assert mean.g > serial*q
    assert 1-serial*q > mean.U(2)
    print(f'PASS near-one powers from {1-q}: retained loss <= {retained*q}, chain loss <= {chain*q}')

    # Founder: endpoint comparisons for the log-concave lower/upper ratio.
    e_upper = Q(11,4)
    assert sum(Q(1,factorial(j)) for j in range(5)) + Q(1,100) < e_upper
    assert 2/e_upper > 1-Q(108,305)
    assert 3/e_upper**2 > mean.U(2)
    print('PASS founder comparison for all 0<p<=1')

    fractional = run_path(str(FRACTIONAL))
    p = Q(1,2)
    m = fractional['middle'](p,p)
    u,_ = fractional['side_bounds'](p,p,m)
    assert Q('.6240') < m < Q('.6241') and Q('.4878') < u(2) < Q('.4879')
    assert fractional['ledger'](p,p,m) == (False,'tail or seed')
    print('OPEN p=1/2: degree-only fractional ledger does not close')


if __name__ == '__main__':
    main()
