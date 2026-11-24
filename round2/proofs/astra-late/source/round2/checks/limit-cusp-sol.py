"""Pin the cusp inputs and check the scale budget and model endpoint controls."""
from collections import defaultdict
from fractions import Fraction as Q
from hashlib import sha1
from math import comb, log
from pathlib import Path


def pins():
    root = Path(__file__).resolve().parents[2]
    for path, blob in (
        ('round2/open/limit-law/smoothing/cusp/README.md', 'c087feab5688dc3390ebbe14983f784771339b54'),
        ('round2/open/limit-law/smoothing/planted/README.md', '494a6357aadafe3cb33802d77935d3ceff172cbf'),
        ('round2/attacks/fable6-density/README.md', 'b2dabc52e9e2bebc293c2708b29ae9d166867550'),
    ):
        data = (root / path).read_bytes()
        assert sha1(f'blob {len(data)}\0'.encode()+data).hexdigest() == blob, path


def budget():
    mu_lo, mu_hi = Q(113, 200), Q(1133, 2000)
    assert Q(1, 2) < mu_lo <= mu_hi < 1
    assert sum(map(Q, ('8/9', '1/16', '7/144'))) == 1
    assert Q(9, 8)*2*mu_hi < Q(4, 3)
    assert 16*2 < 40
    assert Q(9, 8)*Q(41, 20)*mu_hi < Q(4, 3)
    assert 16*Q(41, 20) < 40
    # The exit size budget cannot be applied to the original one-core seed.
    assert Q(9, 8)*2*mu_lo > 1
    for mu in (mu_lo, mu_hi, Q(4, 7)):
        alpha = 1/mu-1
        critical = 1/(1-alpha)
        assert 0 < alpha < 1 and 1 < critical
        assert critical*(alpha-1) == -1 and alpha+1/critical == 1
        assert 2-(2*mu-1)/mu == alpha+1


def convolution(a, b):
    out = defaultdict(Q)
    for x, w in a.items():
        for y, v in b.items():
            out[x+y] += w*v
    return dict(out)


def difference(a, order):
    out = defaultdict(Q)
    for x, w in a.items():
        for k in range(order+1):
            out[x-k] += (-1)**(order-k)*comb(order, k)*w
    return {x: w for x, w in out.items() if w}


def kernel_control():
    product, increments = {0: Q(1)}, {0: Q(1)}
    tv_product = Q(1)
    for k in range(1, 5):
        kernel = {0: Q(1, k+1), k: Q(k, k+1)}
        delta = difference(kernel, 1)
        product = convolution(product, kernel)
        increments = convolution(increments, delta)
        tv_product *= sum(map(abs, delta.values()))
    assert difference(product, 4) == {x: w for x, w in increments.items() if w}
    assert sum(map(abs, increments.values())) <= tv_product


def critical_model():
    # On x∈[-h,1-h] the bare model is (x_+)^(3/4), p*=4.
    # Scale by h=1/N and integrate in t: the fourth power of the
    # normalized first-difference norm has a log coefficient (3/4)^4.
    # The monotone-integrand unit rectangles bound the model integral;
    # the bounded seam and cutoff only give O(1) in this normalization.
    alpha, p = 0.75, 4
    target = alpha**p
    last = None
    for n in (256, 4096, 65536):
        def inc(t):
            return ((t+1)**alpha-t**alpha)**p
        low = (0.25+sum(inc(t) for t in range(1, n)))/log(n)
        high = (1.25+sum(inc(t) for t in range(1, n-1)))/log(n)
        # The limit need not lie between the finite-n rectangles.
        # On [1,n-1], concavity bounds inc(t) between target/(t+1)
        # and target/t; [0,1] contributes between zero and one.
        analytic_low = (0.25+target*log(n/2))/log(n)
        analytic_high = (1.25+target*log(n-1))/log(n)
        assert max(low, analytic_low) < min(high, analytic_high)
        if last is not None:
            assert high-low < last
        last = high-low
    # On the left seam of the fourth difference only the rightmost
    # sample survives: for the same model its integral is exactly h^4/4.
    assert Q(3, 4)*4+1 == 4
    return last


if __name__ == '__main__':
    pins()
    budget()
    kernel_control()
    print('PASS cusp/planted/audit pins, exact repaired star-exit scale budgets and convolution control')
    print(f'PASS critical p=4 model log-coefficient brackets; last width {critical_model():.6f}')
