#!/usr/bin/env python3
"""Price both pre- and post-entrance arms with exact shared apex pressure."""
import argparse
from pathlib import Path

from check import HI, J, LO, Q, Rewards, add, arms, at, band_price, downward, scale, stage_bounds, upward, v, ztuple

GRID = 10**23


def enclosure(corners, upper):
    if upper:
        return tuple(Q(-(-max(col)*GRID//1), GRID) for col in zip(*corners))
    return tuple(Q(min(col)*GRID//1, GRID) for col in zip(*corners))


def run(states, arcs, polys, inf, pres, M, K, top, upper):
    s = LO if upper else HI
    overflow = band_price(states, arcs, inf, top, s, K, pres, upper)
    next_arm = [None] + list(reversed(overflow)) + [{}]
    roles = ('all', 'old', 'young')
    next_host = [{role: ztuple() for role in roles} for _ in range(K+2)]
    for k in range(K-1, -1, -1):
        beta = 1-Q(k+1, top)
        for role in ('all', 'young'):
            a = add(next_arm[k+1][J], next_host[k+1][role])
            b = (beta if not upper else Q(1))
            next_host[k][role] = enclosure([scale(a, b/(s+b))], upper)

    states_info = [(x, v.degree(x), v.size(x)-Q(1, v.degree(x)),
                    v.plus(x), arcs[x]) for x in states]
    for d in range(top-1, 1, -1):
        arm = [{} for _ in range(K+2)]
        host = [{role: ztuple() for role in roles} for _ in range(K+2)]
        p = Q(1,d)
        for k in range(min(K,d-1), 0, -1):
            leaf = d-k-1
            birth = leaf*p
            for x, dx, intrinsic, y, edges in states_info:
                a, b = ztuple(), ztuple()
                if d >= M and polys[x] is not None:
                    pa, pb = polys[x][1][k]
                    a, b = (tuple(at(f,p)/(1+p)**3 for f in pol) for pol in (pa,pb))
                for z, weight in edges.items():
                    a = add(a, scale(arm[k][z], weight))
                a = add(a, scale(arm[k].get(y,ztuple()), p))
                a = add(a, scale(arm[k+1].get(x,ztuple()), birth))
                den = s + intrinsic + p + birth
                a1 = next_arm[k][x]
                # The very same sibling/mother pressure enters reward and gain.
                corners = [scale(add(add(a,scale(b,u)),scale(a1,leaf+Q(1,dx)+u)),
                                 1/(den+leaf+Q(1,dx)+u)) for u in (Q(0),Q(k,2))]
                arm[k][x] = enclosure(corners, upper)
        for k in range(min(K,d-1), -1, -1):
            leaf = d-k-1
            birth = leaf*p
            for role in roles:
                pay = role == 'all' or (role == 'old') == (d < M)
                a = scale(add(host[k+1][role], arm[k+1].get(J,ztuple()) if pay else ztuple()), birth)
                grow = next_host[k][role]
                host[k][role] = enclosure([scale(add(a,scale(grow,g)),1/(s+birth+g))
                                           for g in (Q(leaf), leaf+Q(k+1,2))], upper)
        next_arm, next_host = arm, host
    return {role: scale(val, HI if upper else LO) for role, val in next_host[0].items()}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument('blob', type=Path)
    cli.add_argument('--degree', type=int, default=101)
    cli.add_argument('--tail', type=int, default=101)
    cli.add_argument('--arity', type=int, default=12, choices=range(5,18))
    cli.add_argument('--cap', type=int, default=7, choices=range(6,9))
    a = cli.parse_args()
    assert 32 <= a.degree <= a.tail
    states, arcs = arms(a.cap)
    reward = Rewards(a.blob, states)
    polys = {x: reward.polynomials(x) for x in states}
    inf = {x: p[0] if p else ztuple() for x, p in polys.items()}
    pres = {x: stage_bounds(p,inf[x],a.tail) for x, p in polys.items()}
    print('witness',reward.c.sign,'arm states',len(states),'M',a.degree,'tail',a.tail,'arity',a.arity,flush=True)
    lo = run(states,arcs,polys,inf,pres,a.degree,a.arity,a.tail,False)
    hi = run(states,arcs,polys,inf,pres,a.degree,a.arity,a.tail,True)
    for role in lo:
        assert all(l<=h for l,h in zip(lo[role],hi[role]))
        print(role, tuple(map(downward,lo[role])), tuple(map(upward,hi[role])),
              'joint family sum',downward(sum(lo[role])),upward(sum(hi[role])),flush=True)


if __name__ == '__main__':
    main()
