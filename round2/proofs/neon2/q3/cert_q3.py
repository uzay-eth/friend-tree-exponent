#!/usr/bin/env python3
"""Exact rational certificates from rooted-tree finite-state propagation."""
import argparse
from collections import defaultdict
from functools import lru_cache
from math import lcm

from symmetries import F, advance, decode, labelled, rooted


@lru_cache(maxsize=None)
def stock_lower(k, d, c):
    # Stock birth j, death j/max(d,j+c), killed at zero.
    # 1 - sum_{i>=k} prod_{j<=i} (death_j/birth_j) / sum_{i>=0} prod(...).
    if not k:
        return F(0)
    J = 28
    terms = [F(1)]
    for j in range(1, J + 1):
        terms.append(terms[-1] / max(d, j+c))
    tail = terms[-1] / (max(d, J+1+c)-1)
    survival = sum(terms[:k]) / (sum(terms) + tail)
    return max(survival, 1-F(1,d)**k)


@lru_cache(maxsize=None)
def arrow_table(c):
    # Boundary values f(0)=1 and f(J)=2/(J+1). Solve the exact
    # first-step system for j -> j+1, j-1 or success, with rates
    # j, j/(j+1+c), 1/(j+1+c), respectively.
    J = 30
    slope, const = [F(0)], [F(1)]
    for j in range(1, J):
        up = j*(j+1+c)
        den = up+j+1-j*slope[-1]
        const.append((1+j*const[-1])/den)
        slope.append(up/den)
    f = [F(0)]*J+[F(2,J+1)]
    for j in reversed(range(J)):
        f[j] = slope[j]*f[j+1]+const[j]
    assert all(0 <= b <= a <= 1 for a,b in zip(f,f[1:]))
    for j in range(1,J):
        assert j*(j+1+c)*(f[j+1]-f[j])+j*(f[j-1]-f[j])+1-f[j] == 0
    return f


def observables(adj, n):
    d = len(adj[0]); nbr = adj[0]
    k = sum(len(adj[x]) == 1 for x in nbr)
    upper = F(0)
    for x in nbr:
        r = sum(len(adj[y]) == 1 for y in adj[x] if y != 0)
        c = len(adj[x])-r-1
        upper += arrow_table(c)[r] if r < 30 else F(2,r+1)
        if upper >= 1:
            break
    lower = stock_lower(k,d,d-k) if d > 1 else F(0)
    return F(max(k-1,0),n), F(d,n), lower, min(F(1),upper)


def audit(n=7):
    # Independent per-initiator rational checker and two planted false inequalities.
    bad_lower = bad_upper = 0
    for size, _, _, adj in labelled(n):
        if size < 3:
            continue
        local = adj[:size]
        for u in range(size):
            d = len(local[u]); k = sum(len(local[y]) == 1 for y in local[u]);
            a = sum(F(1,size*len(local[y])) for y in local[u])
            gain = a
            loss = F(k,size*d)
            # A degree-one root has k=0 and on promotion acquires one leaf,
            # while its old neighbour continues to be a core.
            assert (d+a)/(size+1) <= F(d,size)
            assert (k-1+gain-loss)/(size+1) >= F(k-1,size)
            clipped = (gain*k + loss*max(k-2,0) + (1-gain-loss)*max(k-1,0))/(size+1)
            assert clipped >= F(max(k-1,0),size)
            bad_lower += (k+gain-loss)/(size+1) < F(k,size)
            bad_upper += (d+a)/(size+1) < F(d,size)
    assert bad_lower and bad_upper
    print(f'local martingale audit through n={n}: invalid L/n and reverse D/n controls rejected ({bad_lower}, {bad_upper} failures)')


def verify_labelled(N=8):
    """Check aggregate endpoints against raw histories, without signed propagation."""
    raw = defaultdict(lambda: [F(0)]*4)
    for n,w,d,adj in labelled(N):
        if n != N:
            continue
        for u in range(2, min(8,n)+1):
            # Reindex the root to zero, without canonical strings.
            perm = list(range(n))
            perm[u-1], perm[0] = perm[0], perm[u-1]
            y = observables([[perm[z] for z in adj[x]] for x in perm],n)
            for j in range(4):
                raw[u][j] += w*y[j]
    seeds, den = seed_distributions()
    for u, outcome in raw.items():
        w, d = seeds[f'u{u}'], den[u]
        for n in range(u,N):
            w,factor = advance(w,n)
            d *= factor
        assert [sum(x*observables(decode(k),N)[j] for k,x in w.items()) for j in range(4)] == outcome
        assert sum(w.values()) == d
    print(f'four certificate observables independently match the full labelled history sum at n={N}')


def seed_distributions():
    # Use separate labelled enumeration to seed distributions at label birth,
    # re-rooted at that label. Also keep both sides of signed orbit reductions.
    seeds = defaultdict(lambda: defaultdict(int)); den = {}
    for n,w,d,adj in labelled(8):
        den[n] = d
        if n == 2:
            seeds['u2'][rooted(adj,1)] += w
        else:
            seeds[f'u{n}'][rooted(adj,n-1)] += w
    return seeds, den


def bounds(N, full=False):
    seeds, den = seed_distributions()
    active = {'u2': (dict(seeds['u2']),1)}
    diffs = {}
    for n in range(2,N+1):
        if 3 <= n <= min(8, N):
            active[f'u{n}'] = (dict(seeds[f'u{n}']),den[n])
            if n == 3:
                active.update(T3c=({'(()())':1},1), T3l=({'((()))':1},1))
            if n == 4:
                active.update(P4i=({'((())())':1},1), P4e=({'(((())))':1},1))
            p, a = active[f'u{n-1}']; q, b = active[f'u{n}']
            assert a == b
            diffs[n] = ({k:p.get(k,0)-q.get(k,0) for k in p.keys()|q.keys() if p.get(k,0)!=q.get(k,0)},a)
        if n in (6, 10, 12, N):
            obs = {key:observables(decode(key),n) for key in active['u2'][0]}
            print(f'n={n}, rooted states={len(obs)}')

            def evaluate(w,d):
                out = [F(0)]*4
                for k,x in w.items():
                    lo,hi,pl,ph = obs[k]
                    if x < 0:
                        lo,hi,pl,ph = hi,lo,ph,pl
                    for j,y in enumerate((lo,hi,pl,ph)):
                        out[j] += x*y
                return [x/d for x in out]

            def display(name, vals):
                # Keep full mean fractions; encode conservative hub fractions
                # with a common million denominator unless explicitly requested.
                scale = 10**6
                hlo, hhi = vals[2:]
                lo = (hlo*scale).__floor__()
                hi = (hhi*scale).__ceil__()
                print(name, *vals[:2], *(vals[2:] if full else (f'{lo}/{scale}', f'{hi}/{scale}')),
                      '|', *(f'{float(x):.9f}' for x in vals))

            for name,(w,d) in active.items():
                assert sum(w.values()) == d
                vals = evaluate(w,d)
                display(name, vals)
            for v,(w,d) in diffs.items():
                assert sum(w.values()) == 0
                vals = evaluate(w,d)
                rev = evaluate({k:-x for k,x in w.items()},d)
                assert rev == [-vals[1],-vals[0],-vals[3],-vals[2]]
                assert not (rev[0]>0 and vals[0]>0) and not (rev[2]>0 and vals[2]>0)
                display(f'diff{v-1},{v}', vals)
        if n == N:
            return
        cache = {}
        def evolve(w,d):
            out = defaultdict(int)
            for key,x in w.items():
                if key not in cache:
                    cache[key],_ = advance({key:1},n)
                for t,a in cache[key].items():
                    out[t] += x*a
            return {k:x for k,x in out.items() if x}, d*n*lcm(*range(1,n))
        active = {k:evolve(w,d) for k,(w,d) in active.items()}
        diffs = {k:evolve(w,d) for k,(w,d) in diffs.items()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--size', type=int, default=12, choices=range(6,17))
    parser.add_argument('--full-hub-fractions', action='store_true')
    args = parser.parse_args()
    audit()
    verify_labelled()
    bounds(args.size, args.full_hub_fractions)
