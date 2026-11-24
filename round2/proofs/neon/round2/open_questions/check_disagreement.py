"""Independent exact fixtures for Long's aggregate coupling; no imports from Long."""
from fractions import Fraction as F
from itertools import product

S0, Q = F(11, 20), F(19, 20)


def certificate(terms=16, shift=F(5, 2), cap=Q):
    entrance, occupation = F(1), []
    for k in range(1, terms + 1):
        rate = k + shift
        occupation.append(entrance / (rate + S0))
        entrance *= rate / (rate + S0)
    assert sum(occupation) == (1 - entrance) / S0
    rows = [((d-1)*(1/S0 - d*sum(w/(d+k) for k,w in
                                  enumerate(occupation, 1)))/d, d)
            for d in range(2, 64)]
    assert all(0 < mass <= cap for mass, _ in rows)
    assert F(7, 8) < F(15, 16)**2 and F(15, 16) < cap
    assert F(40)/(S0**2*133) == F(16000, 16093) < 1
    return max(rows)


def local_row(d, lo, hi, pred, old, children, stub=0,
              count_children=True, count_returns=True, denom_shift=0):
    """Compute rates directly from two current lists of physical neighbours.

    An absent slot (None) is an unpaired child; old entries are (a,b,infected).
    pre-infection leaf stocks need not be equal on the two sides.
    """
    left = [pred[0]] if pred else []
    right = [pred[1]] if pred else []
    left += [a for a, _, _ in old] + [a for a, _ in children if a]
    right += [b for _, b, _ in old] + [b for _, b in children if b]
    ll, lr = lo-stub-len(left), hi-stub-len(right)
    if min(ll, lr) < 0 or len(old) + bool(pred) + stub > d:
        return False
    z, r, b = abs(hi-lo), len(children), sum(c for _, _, c in old)
    al = ll + sum(F(1, a) for a in left)
    ar = lr + sum(F(1, a) for a in right)
    assert al == lo-stub-sum(1-F(1, a) for a in left)
    assert ar == hi-stub-sum(1-F(1, a) for a in right)
    mismatch = abs(al-ar)
    p = max(F(ll, lo), F(lr, hi))
    assert p <= 1
    assert mismatch <= z + r*count_children + F(1+b*count_returns, 2)
    edge = abs(F(1, lo)-F(1, hi))
    assert edge <= F(z, d*(d+z+denom_shift))
    outward = sum(not c for _, _, c in old)*edge
    assert outward <= 1
    assert mismatch + p + outward <= z+r+b+F(5, 2)
    assert (len(old)*edge <= 1)  # arbitrary seed, no d-1 requirement
    return True


def audit():
    count = 0
    for d, dl, dr, pred, old, children in product(
        range(2, 9), range(4), range(4), [(2, 2), (2, 41), (41, 2), None],
        [[], [(2, 2, False)], [(2, 55, True)],
         [(60, 2, True)]*4, [(4, 4, False)]*5],
        [[], [(2, 2)], [(2, 19)], [(None, 17)], [(25, None)]*3],
    ):
        count += local_row(d, d+dl, d+dr, pred, old, children, int(pred is None))
    big = 10**60
    for a,b in ((big,big+1), (big+1,big), (big,big+10**12), (big,big)):
        count += local_row(big, a, b, (2, 19), [(2, 2, False)], [(None, 70)])
    assert count > 5000
    return count


def rejected(name, f):
    try:
        f()
    except AssertionError:
        print('REJECTED:', name)
        return
    raise AssertionError(f'negative control accepted: {name}')


def unit_jump_fixture(use_shared_source=False):
    # A degree-2 core has two common neighbours; the degree-3 copy
    # also has a leaf. Two common choices match at rate 1/3 each.
    # The two excess sends (1/6 each) can coincide with promotion
    # at the degree-3 copy on shared source clocks.
    events = [(F(1, 3), (0, 1, 0)), (F(1, 3), (0, 0, 1))]
    if use_shared_source:
        events = [(F(1, 3), (0, 1, 1))]
    assert sum(rate*jump[1] for rate,jump in events) == F(1, 3)
    assert sum(rate*jump[2] for rate,jump in events) == F(1, 3)
    assert all(sum(jump) == 1 for _,jump in events)


if __name__ == '__main__':
    m = certificate()
    assert m == (F(1669075573443141901662377825, 1772128093509182310390306279), 6)
    print(f'PASS: 62 exact reproduction rows; maximum {m[0]} at d={m[1]}')
    print(f'PASS: {audit()} exact coupled rate rows, including degrees 10^60')
    unit_jump_fixture()
    rejected('truncate the occupation sum at 8 without weakening the bound',
             lambda: certificate(terms=8))
    rejected('increase reservoir birth rate by one without recertifying',
             lambda: certificate(shift=F(7, 2)))
    rejected('omit unpaired newborn feedback',
             lambda: local_row(3, 3, 3, (2,2), [], [(None,17)], count_children=False))
    rejected('omit repeated outward feedback',
             lambda: local_row(5, 5, 5, (2,2), [(2,100,True)]*3, [], count_returns=False))
    rejected('over-strengthen reciprocal bound denominator',
             lambda: local_row(3, 3, 4, (2,2), [(2,2,False)], [], denom_shift=1))
    rejected('share a source ring across promotion and outward infection',
             lambda: unit_jump_fixture(use_shared_source=True))
