"""Cold replay: authenticate writer inputs and audit the near-zero margins."""
from fractions import Fraction as Q
from hashlib import sha1
from pathlib import Path
from runpy import run_path


AGE = Path(__file__).resolve().parents[1] / 'open/age-order'
WRITER = AGE / 'global-transport/near-zero'


def git_blob(data):
    return sha1(f'blob {len(data)}\0'.encode() + data).hexdigest()


def main():
    for name, expected in (
        ('README.md', 'b18258099359cc8bb7a5fefaa8f29ad1df554795'),
        ('check.py', '4452e5abae9b17816c1b7abac926baf1cd34b4ad'),
    ):
        assert git_blob((WRITER / name).read_bytes()) == expected, name
    print('PASS both frozen writer blobs')

    # The exp series gives e>1+1+1/2, hence e^d>8 for d>=3.
    assert Q(5, 2)**3 > 8
    # Stronger e>8/3 gives target-channel success >5/8.
    assert sum(Q(1, j) for j in (1, 1, 2, 6, 24)) > Q(8, 3)
    assert Q(7, 8)*Q(5, 8)*Q(1, 16) > Q(1, 32)
    assert 2*31**32 < 32**32
    # d=2^k, h=16k: h/2=8k > log_2(2d^4).
    assert 8*10 > 4*10+1
    d = 2**10
    assert 16*(1-Q(1, d)) > 8

    # Multiply each error by d. All nonconstant terms decrease in d.
    hit = Q(1, d**5) + Q(2, d)
    stock = 1 + Q(2, d**6)
    # Chebyshev: variance d at deviation d/2; Markov: mean one at d/4.
    reserve = d*(Q(d, (d//2)**2) + Q(1, d//4))
    assert reserve == 8 and hit + stock + reserve < 10
    # The loss-budget endpoint is closed: successful remaining mass is
    # > d/4 times the clock discount, strictly above the threshold.
    assert d//4 > 1
    print('PASS race, finite-hit, stock and Gamma-minus-Exp error budgets')

    half_amplitude_loss = Q(2*10, 1024)
    hub = Q(2968, 100000)
    assert half_amplitude_loss == Q(5, 256) < hub
    assert hub-half_amplitude_loss == Q(8119, 800000)
    early = run_path(str(AGE / 'early/check.py'))
    base23 = (early['silent_bound'](1, 4)-Q(1, 4))/2
    base34, _ = early['adjacent']((0,), 4)
    assert base23 == Q(127, 2440)
    assert base34 > Q('0.022645238637')
    assert half_amplitude_loss < min(base23, base34)
    print('PASS exact early seam: p2-p3 >= 127/2440; p3-p4 > .022645238637')
    print(f'PASS loss coefficient {half_amplitude_loss}; hub margin coefficient {hub-half_amplitude_loss}')


if __name__ == '__main__':
    main()
