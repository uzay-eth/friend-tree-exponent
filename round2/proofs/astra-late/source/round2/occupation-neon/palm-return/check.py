#!/usr/bin/env python3
"""Exact physical birth-Palm minorant, signed occupation and adverse controls."""
from fractions import Fraction as Q
from hashlib import sha256
import json
from pathlib import Path

LO, HI = Q('.565673'), Q('.565775')
TYPES = tuple(range(2, 21))
J = (2, 0)
GRID = 10**16
CAP = 34
OCC = Path(__file__).resolve().parents[1]


def down(x):
    return Q(x * GRID // 1, GRID)


def states(cap):
    return [(d, k) for d in range(cap, 1, -1) for k in range(d-1, -1, -1)]


def rates(i, d, k):
    l = d-k-1
    return Q(l, d), Q(l), Q(l)+Q(k, 2)+Q(1, i)


def coherent(i):
    """One common positive suboccupation, by thinning and extra killing."""
    mass = {J: Q(1)}
    occupation = {}
    births = dict.fromkeys(TYPES, Q(0))
    exit_mass = Q(0)
    kill_mass = Q(0)
    for d, k in reversed(states(CAP)):
        b, g, h = rates(i, d, k)
        x = (d, k)
        w = mass.get(x, 0)/(HI+b+h)
        occupation[x] = w
        births[min(d, 20)] += w*b
        kill_mass += w*(h-g)
        if b:
            y = (d, k+1)
            mass[y] = mass.get(y, 0)+w*b
        if d < CAP:
            y = (d+1, k)
            mass[y] = mass.get(y, 0)+w*g
        else:
            exit_mass += w*g
    assert HI*sum(occupation.values()) + kill_mass + exit_mass == 1
    assert all(births.values()) and sum(births.values()) < 1
    degree = tuple(sum(LO*w for (d, _), w in occupation.items() if min(d, 20)==j) for j in TYPES)
    assert all(degree)
    return tuple(map(down, births.values())), degree


def lower(i, cap, payoff, *, normalized):
    """All-age lower price of a signed score that is zero beyond cap."""
    v = {}
    for d, k in states(cap):
        b, g, h = rates(i, d, k)
        a, r = payoff(d, k)
        next_bud = v.get((d, k+1), 0)
        next_gain = v.get((d+1, k), 0)
        actions = [(s+b+z, (s if normalized else 1)*(a+r*(z-g))
                    +b*next_bud+z*next_gain)
                   for s in ((LO, HI) if normalized else (HI,)) for z in (g, h)]
        v[d, k] = down(min(num/den for den, num in actions))
        assert all(den*v[d, k] <= num for den, num in actions)
    return v[J]


def completion(A, f):
    """Propose v by monotone robust iteration, verify using the exact dual."""
    v = [Q(0)]*len(f)
    for _ in range(48):
        w = [f_i+sum(a*x for a, x in zip(row, v)) for f_i, row in zip(f, A)]
        base = min(w)
        v = [down(x-base) for x in w]
    assert min(v) == 0
    bound = min(f_i-v_i+sum(a*x for a, x in zip(row, v))
                for f_i, v_i, row in zip(f, v, A))
    return max(min(f), bound)


def frows():
    data = (OCC/'rows'/'charged-f.json').read_bytes()
    assert sha256(data).hexdigest() == '415a712878e03b71c431a36bf4cf8eb11a047c2a47cafb5d9298fe23bc1561da'
    data = json.loads(data)
    tables = []
    for side, sign, lam in [('lower', 1, LO), ('upper', -1, HI)]:
        item = data[side]
        assert (item['D0'], item['sign'], item['S'], Q(item['lambda'])) == (18, sign, 10**11, lam)
        table = {}
        for row in item['rows']:
            l, k, h = row['state']
            assert h == 1
            d = l+k+h
            old = row['f']
            gain = row['f_gain']-old
            a = l*gain-lam*old
            if l:
                edge = sum(val*Q(w) for val, w in row['g_bud_newborn'])
                a += Q(l, d)*(row['f_bud']+row['f_newborn']+edge-old)
            r_p, r_c = map(Q, row['center'])
            psi = gain if d < 18 else Q(row['c2']*d+row['c3'], d**3)
            a += (r_p+k*r_c)*psi
            a, r = Q(sign, item['S'])*a, Q(sign, item['S'])*(gain-psi)
            corners = [a] if d < 18 else [a, a+r*Q(k+1, 2)]
            assert corners == list(map(Q, row['signed_corners'])) and min(corners) >= 0
            table[d, k] = a, r
        assert set(table) == set(states(64))
        tables.append(table)
    return tables


def family(mother_degree, child_degree, *, older_core=False):
    """Replay labelled real arrows; founder 0, nonfounder mother 1, child 4."""
    neighbors = [{1}, {0}]

    def add(donor, receiver):
        assert receiver in neighbors[donor]
        rate = Q(1, len(neighbors[donor]))
        neighbors[receiver].add(len(neighbors))
        neighbors.append({receiver})
        return rate

    add(1, 0)   # 0 becomes founder
    add(0, 1)   # 1 becomes a core; leaf 3 at 1
    add(3, 1)   # leaf 4 at 1
    for _ in range(mother_degree-3):
        add(3, 1)
    if older_core:
        add(1, 3)
    assert len(neighbors[4]) == 1
    assert add(1, 4) == Q(1, mother_degree)  # actual core birth
    assert len(neighbors[1]) == mother_degree
    birth_type = min(len(neighbors[1]), 20)
    first_leaf = len(neighbors)-1
    for _ in range(child_degree-2):
        add(first_leaf, 4)
    assert len(neighbors[4]) == child_degree
    assert all(len(neighbors[x]) == 1 for x in neighbors[4]-{1})
    mother_gain = sum(Q(1, len(neighbors[x])) for x in neighbors[1])
    return birth_type, Q(1, len(neighbors[1])), mother_gain


def obstruction(tables):
    assert family(20, 2)[:2] == (20, Q(1, 20))
    assert family(1000, 2)[:2] == (20, Q(1, 1000))
    # In the degree-two virgin, b=1/2 and the constant loss excluding the
    # mother pressure is 3/2. Bound the extra hazard from the mother at 20
    # below until she first grows (growth rate at most 20). Bound it above
    # on the whole virgin residence for a mother starting at 1000.
    gap = Q(1, 2)/(HI+Q(3, 2)) * (
        Q(1, 20)/(HI+Q(3, 2)+Q(1, 20)+20) - Q(1, 1000)/(LO+Q(3, 2)))
    assert gap > Q('.000430038')
    print('same-type genuine return T[20,2] uniform row gap > 0.000430038')

    # At fixed birth degree 3, the virgin loss h=3/2+1/3 agrees.
    # The full labelled generator on h before virgin departure differs:
    # only a maternal gain changes h, by 1/4-1/3. Hence S''=h²−A h.
    i, reciprocal, leaf_gain = family(3, 2)
    j, same, core_gain = family(3, 2, older_core=True)
    assert (i, reciprocal) == (j, same) == (3, Q(1, 3))
    second = (leaf_gain-core_gain)*(reciprocal-Q(1, 4))
    assert second == Q(1, 24)
    print('same exact birth degree 3: virgin survival second-derivative gap', second)

    for name, table in zip(('lower', 'upper'), tables):
        d, k = next(x for x in reversed(states(64)) if x[1] == 0 and table[x][1])
        a, slope = table[d, k]
        i, r, _ = family(20, d)
        j, z, _ = family(1000, d)
        assert i == j == 20 and abs((a+r*slope)-(a+z*slope)) == abs((r-z)*slope) > 0
        print(name, 'same-type exact D18 first-body gap at', (d, k), '=', abs((r-z)*slope))


def show(a, b):
    scale = 10**12
    def ticks(n):
        whole, frac = divmod(abs(n), scale)
        return f'{"-" if n<0 else ""}{whole}.{frac:012d}'
    return f'[{ticks(a*scale//1)}, {ticks(-(-b*scale//1))}]'


def ledger(B, prices, p, q, t, masked_lo, masked_hi):
    edge_lo, edge_hi = Q('.565677980686114'), Q('.565770513279953')
    p_ceiling, q_ceiling = Q('1.356833021311'), Q('1.367169995210')
    paid_p, paid_q = (completion(B, prices(lambda d, k: table[d, k])) for table in (p, q))
    rem_p = (edge_hi-LO)*p_ceiling-paid_p
    rem_q = (HI-edge_lo)*q_ceiling-paid_q
    assert min(paid_p, paid_q, rem_p, rem_q)>0
    full_lo, full_hi = masked_lo-t*rem_q, masked_hi+rem_p
    assert full_lo < 0 < full_hi
    print('masked f lower payments (micro/core):', show(10**6*paid_p, 10**6*paid_p),
          show(10**6*paid_q, 10**6*paid_q))
    print('whole contrast after separate omitted budgets (micro/core):', show(10**6*full_lo, 10**6*full_hi))
    odds_lo, odds_hi = ((x-LO)/(HI-x) for x in (edge_lo, edge_hi))
    print('full D18 product odds inherited from frozen EDGE:', show(odds_lo, odds_hi))
    assert edge_hi-edge_lo == Q('.000092532593839')


def main():
    assert __debug__
    tables = frows()
    obstruction(tables)
    A, O = zip(*(coherent(i) for i in TYPES))
    occ = list(map(sum, O))
    B = tuple(tuple(lower(i, CAP, lambda d, k: (Q(d-k-1, d)*(min(d, 20)==j), Q(0)),
                          normalized=False) for j in TYPES) for i in TYPES)
    assert all(0<a<=b for ar, br in zip(A, B) for a, b in zip(ar, br))
    rows = list(map(sum, B))
    assert max(rows)<1
    print('coherent normalized suboccupation mass', show(min(occ), max(occ)))
    print('stronger physical birth lower row sums', show(min(rows), max(rows)))
    print('worst-case completion deficit range', show(1-max(rows), 1-min(rows)))

    # The same instantaneous donor pressure is used in paired original rows
    # and in the gain successor. Keep an explicit mask, zero past the cut.
    cap = 18
    p, q = tables
    t = Q(31, 25)
    contrast = lambda d, k: tuple(x-t*y for x, y in zip(p[d, k], q[d, k]))

    def prices(f):
        return tuple(lower(i, cap, f, normalized=True) for i in TYPES)

    low = prices(contrast)
    neg = prices(lambda d, k: tuple(-x for x in contrast(d, k)))
    mixed_lo, mixed_hi = completion(B, low), -completion(B, neg)
    assert min(low) <= mixed_lo <= mixed_hi <= -min(neg)
    print('D18 masked f contrast (micro/core), arbitrary types:', show(10**6*min(low), -10**6*min(neg)))
    print('D18 masked f contrast (micro/core), physical Palm:', show(10**6*mixed_lo, 10**6*mixed_hi))
    ledger(B, prices, p, q, t, mixed_lo, mixed_hi)
    print('PASS exact physical birth/occupation minorants, signed duals and reachable obstructions')


if __name__ == '__main__':
    main()
