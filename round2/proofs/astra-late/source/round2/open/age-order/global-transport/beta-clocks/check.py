"""Exact finite-size race laws and sharp-tail comparisons for the random friend tree."""

from argparse import ArgumentParser
from fractions import Fraction as Q
from functools import cache
from hashlib import sha1
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
from runpy import run_path

HERE = Path(__file__).resolve().parent
SCALED = HERE.parent / "annealed/extension/third-physical/scaled/check.py"
LEDGER_PIN = "85e9ce571c1b6a842e6d3aa8faa71dbda147b2df"
GRID = 10**9
MIN_SIZE = 21
POISSON_AT_ZERO = Q(367879441, GRID)


def ceiling(x):
    return -((-x.numerator) // x.denominator)


def blob(path):
    data = path.read_bytes()
    return sha1(f"blob {len(data)}\0".encode() + data).hexdigest()


def rising(x, n):
    return prod(x + i for i in range(n))


def root_ticks(a, k, unit=GRID):
    """Floor of unit times the positive kth root, using integer comparisons."""
    assert a >= 0 and k >= 1
    num = a.numerator * unit**k
    den = a.denominator
    lo, hi = 0, unit
    while hi**k * den <= num:
        hi *= 2
    while hi - lo > 1:
        mid = (lo + hi) // 2
        if mid**k * den <= num:
            lo = mid
        else:
            hi = mid
    return lo


def race_checks():
    """First-step equations checked against the product law and both seam cases."""
    count = 0
    for m in range(4, 20):
        for p in range(1, 8):
            for h in (Q(1, 2), Q(3, 2), Q(8, 3), Q(m)):
                for r in (Q(1, 3), Q(1), Q(5, 2), Q(m)):
                    terminal = rising(r, p)
                    now = terminal * h / ((p + h) * rising(m, p))
                    following = terminal * h / ((p + h) * rising(m + 1, p))
                    assert m * now == h * terminal / rising(m + 1, p) + (m - h) * following
                    count += 1
    return count


def physical_checks():
    """Read the two monitored rates directly from labelled physical arrows."""
    count = 0
    for n in range(2,8):
        for parents in product(*(range(i) for i in range(1,n))):
            adj = [set() for _ in range(n)]
            for child,parent in enumerate(parents,1):
                adj[parent].add(child)
                adj[child].add(parent)
            for tag,neighbours in enumerate(adj):
                leaves = sorted(v for v in neighbours if len(adj[v])==1)
                for stock in range(1,min(len(leaves),len(neighbours)-1)+1):
                    chosen = set(leaves[:stock])
                    up = down = Q(0)
                    for sender, destinations in enumerate(adj):
                        for dest in destinations:
                            rate = Q(1,len(destinations))
                            if sender in chosen and dest==tag:
                                up += rate
                            if sender==tag and dest in chosen:
                                down += rate
                    assert up==stock and down==Q(stock,len(neighbours))
                    virtual = Q(stock,stock+1)
                    assert down<=virtual and up+virtual<=n
                    count += 1
    return count


def product_floor(left, right, stages=8, mesh=128, minimum=MIN_SIZE):
    """Integrate the positive source at each visited stock using beta clocks."""
    assert 0 < left <= right < minimum and 1 <= stages < minimum
    unit = 10**20
    t_max = right / (1 - right/minimum)
    ratio = Q(mesh+1, mesh)
    rates = [j+Q(j,j+1) for j in range(1,stages+1)]
    factors = [root_ticks(1/ratio**h.numerator, h.denominator, unit) for h in rates]
    coeffs = [
        [prod(b/(b-a) for b in rates[:j] if b != a) for a in rates[:j]]
        for j in range(1,stages+1)
    ]
    # Each endpoint e^(-h_i*t) is bounded separately. Only nonnegative
    # differences of survival probabilities enter the occupation lower bound.
    lo_surv = hi_surv = [unit]*stages
    previous = [unit]*stages
    lo_exp = hi_exp = unit
    source = [0]*stages
    while t_max * (Q(hi_exp,unit)-1) < minimum:
        next_lo = lo_exp*(mesh+1)//mesh
        next_hi = (hi_exp*(mesh+1)+mesh-1)//mesh
        lo_surv = [v*c//unit for v,c in zip(lo_surv,factors)]
        hi_surv = [(v*(c+1)+unit-1)//unit for v,c in zip(hi_surv,factors)]
        z = t_max*(Q(next_hi,unit)-1)
        base = max(Q(0), 1 - z/minimum)
        gamma_limit = exponential_lower(z, unit)
        increment = [
            Q(comb(minimum-1,v))*(left*lo_exp/(unit*minimum))**v
            * min(base**(minimum-1-v), gamma_limit) for v in range(stages)
        ]
        for i, weights in enumerate(coeffs):
            lower = sum(c*(lo_surv[l] if c>=0 else hi_surv[l]) for l,c in enumerate(weights))
            upper = sum(c*(hi_surv[l] if c>=0 else lo_surv[l]) for l,c in enumerate(weights))
            mass = max(0, previous[i]-ceiling(upper))
            previous[i] = max(0, lower//1)
            value = sum(increment[v]/factorial(i-v) for v in range(i+1))
            source[i] += (mass*value)//1
        lo_exp,hi_exp = next_lo,next_hi
    event_prob = Q(1)
    gain = Q(0)
    for j,h in enumerate(rates,1):
        gain += event_prob*source[j-1]/(unit*(j+1)*h)
        event_prob *= j/h
    return Q((POISSON_AT_ZERO*(1+gain)*GRID)//1,GRID)


def exponential_lower(z, unit):
    """Round exp(-z) down via a signed Taylor bound and repeated squaring."""
    squares = 0
    while z > Q(1,4):
        z /= 2
        squares += 1
    low = (unit*sum((-z)**i/factorial(i) for i in range(12)))//1
    for _ in range(squares):
        low = low*low//unit
    return Q(low,unit)


@cache
def wait_cells(panels, b0):
    # e^(-E) is uniform when E has standard exponential law. W is bounded
    # below on each equal-probability cell by taking its right endpoint.
    return tuple(
        b0 * (1 - Q(root_ticks(Q(i, panels), b0) + 1, GRID))
        if i < panels else Q(0)
        for i in range(1, panels + 1)
    )


@cache
def gamma_lower(k):
    """Lower bound for Gamma(1/k) using the signed Taylor series and rectangles."""
    small = sum(Q((-1)**j * k, (1 + j*k)*factorial(j)) for j in range(12))
    slices = 96
    step = Q(1, slices)
    tick_ratio = (GRID * (1 - step + step**2/2 - step**3/6)) // 1
    exp_ticks = GRID
    integral_ticks = 0
    for i in range(1, 9*slices + 1):
        exp_ticks = exp_ticks*tick_ratio // GRID
        if i > slices:
            root = root_ticks(Q(i, slices), k)
            integral_ticks += exp_ticks*root // (GRID*i)
    return small + Q(integral_ticks, GRID)


def price_caps(y, panels=128, minimum=MIN_SIZE):

    @cache
    def cap(k, cutoff):
        assert cutoff == y
        if k >= 32:
            return Q(2) / ((k - 1) * y)
        total = Q(0)
        for w in wait_cells(panels, minimum - 1):
            v = y + max(0, w)
            root = root_ticks(1 / v**(k - 1), k) + 1
            total += Q(root, GRID)
        cap_value = total / (panels * gamma_lower(k))
        return Q(ceiling(GRID * cap_value), GRID)
    return cap


BOXES = (
    ("23/10", "5/2", (".00689",".00620",".00516",".01524")),
    ("5/2", "11/4", (".009",".006",".006",".016")),
    ("11/4", "3", (".007",".005",".004",".015")),
    ("3", "7/2", (".006",".004",".004",".011")),
    ("7/2", "9/2", (".006",".004",".004",".011")),
)


def ledger_box(left, right, tariff, stages=8, panels=128, minimum=MIN_SIZE):
    local = run_path(str(HERE/"ledger.py"))
    cap = price_caps(left, panels, minimum)
    middle = product_floor(left, right, stages, minimum=minimum)
    prepared = local["prepare"](middle, cap, left)
    assert local["certify"](prepared, tuple(map(Q,tariff))), (left,right,tariff)
    print(f"PASS y in [{left},{right}], m >= {minimum}: "
          f"serial >= {middle}, U2 <= {cap(2,left)}, U3 <= {cap(3,left)}",
          flush=True)


if __name__ == "__main__":
    parser = ArgumentParser(description=__doc__)
    parser.add_argument("--box", nargs=2, type=Q)
    parser.add_argument("--tariff", nargs=4, type=Q)
    args = parser.parse_args()
    assert bool(args.box) == bool(args.tariff)
    assert blob(SCALED) == LEDGER_PIN
    assert sum(Q((-1)**j, factorial(j)) for j in range(16)) > POISSON_AT_ZERO
    print(f"PASS {race_checks()} exact stopped-beta moment rows", flush=True)
    print(f"PASS {physical_checks()} labelled stock subsets and both monitored rates", flush=True)
    boxes = [(*args.box,args.tariff)] if args.box else BOXES
    for left,right,tariff in boxes:
        ledger_box(Q(left),Q(right),tariff)
    if not args.box:
        assert all(Q(a[1])==Q(b[0]) for a,b in zip(BOXES,BOXES[1:]))
        print("PASS all finite sizes b>=20, (b+1)x>=23/10 with the reviewed 9/2 splice")
