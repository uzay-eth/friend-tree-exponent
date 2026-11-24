#!/usr/bin/env python3
"""Exact checks: ordinary Dirichlet killing and independent resolvent jets in 1/s."""

from functools import cache
from fractions import Fraction as F
from math import log
from random import Random


def orient(g, v, p=-1):
    return tuple(sorted(orient(g, x, v) for x in g[v] if x != p))


@cache
def graph(t):
    g = []

    def visit(branch, parent):
        v = len(g)
        g.append([])
        if parent != -1:
            g[v].append(parent)
            g[parent].append(v)
        for child in branch:
            visit(child, v)

    visit(t, -1)
    return tuple(map(tuple, g))


def add(g, v):
    a = list(g)
    a[v] += (len(g),)
    return tuple(a) + ((v,),)


@cache
def size(t):
    return len(graph(t))


@cache
def degree(t, inert):
    return len(t) + inert


@cache
def children(t):
    return tuple(c for c in t if c)


@cache
def fringes(t, inert):
    out = [t] if inert else []
    for c in children(t):
        out.extend(fringes(c, True))
    return tuple(out)


@cache
def cores(t, inert):
    return len(fringes(t, inert)) + (not inert)


@cache
def b(t, inert):
    return F(t.count(()), degree(t, inert))


@cache
def transitions(t, inert):
    """Enumerate every real clock and its uniform neighbour choices; omit inert exits."""
    g = graph(t)
    rates = {}
    for v, adj in enumerate(g):
        rate = F(1, len(adj) + (inert and v == 0))
        for x in adj:
            a = orient(add(g, x), 0)
            rates[a] = rates.get(a, 0) + rate
    return tuple(rates.items())


def gen(f, t, inert):
    return sum(q * (f(a) - f(t)) for a, q in transitions(t, inert))


@cache
def plus(t):
    return ((),) + t


def arrows(f, t, inert):
    return sum(f(plus(c)) - f(c) for c in children(t)) / F(degree(t, inert))


J = ((),)
P3 = ((), ())
CAP = 7
S = F(2)
ORDERS = range(1, 8)  # Laplace coefficients s^-1 .. s^-7 encode derivatives 0 .. 6.


def score(t, inert):
    return F(degree(t, inert) >= 3)


@cache
def vs_check(t, birth):
    """Independently apply the triangular killed intrinsic resolvent V_s."""
    if size(t) > CAP:
        return F(0)
    rates = transitions(t, True)
    reward = b(t, True) if birth else score(t, True)
    return (reward + arrows(lambda a: planted(a, birth), t, True)
            + sum(q * vs_check(a, birth) for a, q in rates)) / (S + sum(q for _, q in rates))


@cache
def planted_solution():
    """Solve both (s-L0-C Delta)u = b or psi by rational elimination."""
    states = [t for level in levels(J, CAP, True) for t in level]
    idx = {t: i for i, t in enumerate(states)}
    n = len(states)
    rows = []
    for t in states:
        row = [F(0)] * n + [score(t, True), b(t, True)]
        rates = transitions(t, True)
        row[idx[t]] = S + sum(q for _, q in rates)
        for a, q in rates:
            if a in idx:
                row[idx[a]] -= q
        for c in children(t):
            q = F(1, degree(t, True))
            row[idx[c]] += q
            row[idx[plus(c)]] -= q
        rows.append(row)
    for i in range(n):
        j = next(j for j in range(i, n) if rows[j][i])
        rows[i], rows[j] = rows[j], rows[i]
        pivot = rows[i][i]
        rows[i][i:] = [v / pivot for v in rows[i][i:]]
        for j in range(i + 1, n):
            q = rows[j][i]
            if q:
                rows[j][i:] = [a - q*b for a, b in zip(rows[j][i:], rows[i][i:])]
    sol = [(F(0), F(0))] * n
    for i in range(n - 1, -1, -1):
        sol[i] = tuple(rows[i][n+k] - sum(rows[i][j]*sol[j][k] for j in range(i+1, n))
                       for k in range(2))
    return dict(zip(states, sol))


@cache
def planted(t, birth):
    return planted_solution().get(t, (F(0), F(0)))[birth]


@cache
def free(t, birth):
    if size(t) > CAP:
        return F(0)
    rates = transitions(t, False)
    rhs = b(t, False) if birth else score(t, False) - S * planted(J, False)
    return (rhs + sum(q * free(a, birth) for a, q in rates)
            + arrows(lambda a: planted(a, birth), t, False)) / (S + sum(q for _, q in rates))


@cache
def combined(t, birth):
    return free(t, birth) + sum(planted(a, birth) - (0 if birth else planted(J, False))
                                for a in fringes(t, False))


@cache
def pj(t, n, birth):
    """Coefficient of s^-n in the unrestricted planted resolvent."""
    if n == 0:
        return F(0)
    if n == 1:
        return b(t, True) if birth else score(t, True)
    f = lambda a: pj(a, n - 1, birth)
    return gen(f, t, True) + arrows(f, t, True)


@cache
def fj(t, n, birth):
    if n == 0:
        return F(0)
    local = (b(t, False) if birth else score(t, False)) if n == 1 else F(0)
    return (local - (0 if birth else pj(J, n, False))
            + gen(lambda a: fj(a, n - 1, birth), t, False)
            + arrows(lambda a: pj(a, n - 1, birth), t, False))


@cache
def cj(t, n, birth):
    return fj(t, n, birth) + sum(pj(a, n, birth) - (0 if birth else pj(J, n, False))
                                  for a in fringes(t, False))


@cache
def raw_jet(t, n, birth):
    """Direct free-tree powers, without any oriented renewal decomposition."""
    if n == 1:
        if birth:
            g = graph(t)
            return sum(F(sum(len(g[x]) == 1 for x in adj), len(adj))
                       for adj in g if len(adj) > 1)
        return F(sum(len(adj) >= 3 for adj in graph(t)))
    return gen(lambda a: raw_jet(a, n - 1, birth), t, False)


@cache
def core_jet(t, n):
    if n == 1:
        return F(sum(len(adj) > 1 for adj in graph(t)))
    return raw_jet(t, n - 1, True)


@cache
def killed_answer(t, birth):
    """Direct killed physical resolvent of the RHS of (8) or (10)."""
    if size(t) > CAP:
        return F(0)
    rates = transitions(t, False)
    rhs = ((1 - planted(J, True)) * raw_jet(t, 1, True) if birth else
           raw_jet(t, 1, False) - S * planted(J, False) * core_jet(t, 1))
    return (rhs + sum(q * killed_answer(a, birth) for a, q in rates)) / (S + sum(q for _, q in rates))


@cache
def cutoff_correction(t, birth):
    """Discounted boundary value lost by killing the entire aggregate."""
    if size(t) > CAP:
        return combined(t, birth)
    rates = transitions(t, False)
    return sum(q * cutoff_correction(a, birth) for a, q in rates) / (S + sum(q for _, q in rates))


def levels(start, limit, inert):
    ts = {start}
    while ts:
        yield sorted(ts)
        if size(next(iter(ts))) == limit:
            break
        ts = {a for t in ts for a, _ in transitions(t, inert)}
        if not inert:
            ts = {orient(graph(t), v) for t in ts for v, adj in enumerate(graph(t))
                  if len(adj) > 1}


def tree_checks():
    caught = {}
    pcount = fcount = paircount = boundary = cutoff_roots = core_jumps = leaf_jumps = 0
    for level in levels(J, CAP, True):
        for t in level:
            for birth in (True, False):
                assert planted(t, birth) == vs_check(t, birth)
    for level in levels(J, 6, True):
        for t in level:
            pcount += 1
            for seed in (3, 13, 37):
                def w(a):
                    key = sum((i + 1)**2 * (len(adj) + 3)**3
                              for i, adj in enumerate(graph(a)))
                    return F((key * seed + size(a)**seed) % 101 - 50, 17)

                lhs = gen(lambda a: sum(map(w, fringes(a, True))), t, True)
                internal = sum(gen(w, a, True) for a in fringes(t, True))
                external = sum(arrows(w, a, True) for a in fringes(t, True))
                new = w(J) * sum(b(a, True) for a in fringes(t, True))
                assert lhs == internal + external + new
                if lhs != internal - external + new:
                    caught['reverse C Delta'] = True
                no_slot = sum(gen(w, a, False) for a in fringes(t, True))
                if lhs != no_slot + external + new:
                    caught['forget inert parent in L0'] = True

    for level in levels(P3, CAP, False):
        for t in level:
            g = graph(t)
            roots = [orient(g, v) for v, adj in enumerate(g) if len(adj) > 1]
            if size(t) <= 6:
                fcount += 1
                for t2 in roots:
                    paircount += 1
                    for birth in (True, False):
                        assert killed_answer(t, birth) == killed_answer(t2, birth)
                    for n in ORDERS:
                        assert cj(t, n, True) == cj(t2, n, True)
                        assert cj(t, n, False) == cj(t2, n, False)
                    if combined(t, True) != combined(t2, True):
                        cutoff_roots += 1
                assert any(len(adj) >= 2 for adj in g)
                assert sum(score(a, True) for a in fringes(t, False)) + score(t, False) == raw_jet(t, 1, False)
                if sum(score(a, False) for a in fringes(t, False)) + score(t, False) != raw_jet(t, 1, False):
                    caught['forget parent in degree score'] = True
                for n in ORDERS:
                    assert cj(t, n, True) == raw_jet(t, n, True) - sum(
                        pj(J, i, True) * raw_jet(t, n - i, True) for i in range(1, n))
                    assert cj(t, n, False) == raw_jet(t, n, False) - sum(
                        pj(J, i, False) * core_jet(t, n + 1 - i) for i in range(1, n + 1))
                for x, adj in enumerate(g):
                    is_core = len(adj) > 1
                    p = x if is_core else adj[0]
                    t0, t1 = orient(g, p), orient(add(g, x), p)
                    if is_core:
                        core_jumps += 1
                        assert sorted(fringes(t0, False)) == sorted(fringes(t1, False))
                    else:
                        leaf_jumps += 1
                        assert sorted(fringes(t0, False) + (J,)) == sorted(fringes(t1, False))
                        assert len(t0) == len(t1) and len(children(t0)) + 1 == len(children(t1))
                    assert combined(t1, True) - combined(t0, True) == free(t1, True) - free(t0, True) + (0 if is_core else planted(J, True))
                    assert combined(t1, False) - combined(t0, False) == free(t1, False) - free(t0, False)
                    for n in ORDERS:
                        diff = cj(t1, n, True) - cj(t0, n, True)
                        fdiff = fj(t1, n, True) - fj(t0, n, True)
                        assert diff == fdiff + (0 if is_core else pj(J, n, True))
                        assert cj(t1, n, False) - cj(t0, n, False) == fj(t1, n, False) - fj(t0, n, False)
                        if not is_core and diff != fdiff:
                            caught['drop virgin in leaf jump'] = True

            for birth in (True, False):
                rhs = ((1 - planted(J, True)) * raw_jet(t, 1, True) if birth else
                       raw_jet(t, 1, False) - S * planted(J, False) * cores(t, False))
                # Full jump rule with local killed resolvents, retaining component values
                # of the enlarged tree. An additional whole-tree cemetery loses those.
                lhs = S * combined(t, birth) - gen(lambda a: combined(a, birth), t, False)
                assert lhs == rhs
                assert combined(t, birth) - cutoff_correction(t, birth) == killed_answer(t, birth)
                assert S * killed_answer(t, birth) - gen(lambda a: killed_answer(a, birth), t, False) == rhs
                cemetery = lhs + sum(q * combined(a, birth) for a, q in transitions(t, False)
                                     if size(a) > CAP)
                if cemetery != rhs:
                    boundary += 1
                    caught['silently kill whole aggregate'] = True
    assert cutoff_roots and boundary
    print(f'Trees: {pcount} planted, {fcount} core-rooted free, {paircount} root comparisons; '
          f'{core_jumps} core and {leaf_jumps} leaf attachments. Exact Taylor orders 0..6.')
    print(f'Cutoff: {cutoff_roots} rooted pairs differ with truncated inputs; '
          f'{boundary} nonzero whole-tree cemetery defects.')
    return caught


def stopped_checks(caught):
    rng = Random(505)
    end, count = 256, 0
    paths = [([0, 0, 0] + [n for n in range(3, end + 1)], F(1, 5)),
             ([0, 0, 0] + [max(1, int(n / log(n)**2)) for n in range(3, end + 1)], F(1)),
             ([0, 0, 0] + [n for n in range(3, end + 1)], F(1))]
    for i in range(520):
        c = [0, 0, 0, 1]
        for n in range(4, end + 1):
            c.append(min(n, c[-1] + rng.randrange(3 if i % 5 == 0 else 2)))
        paths.append((c, F(rng.randint(1, 18), 9)))
    seen = set()
    for c, k in paths:
        total = F(0)
        tau = None
        for r in range(3, end + 1):
            total += F(c[r], r * (r+1))
            if total > k:
                tau = r
                break
        if tau is None:
            # Exact constant continuation; for c[n]=n (divergent continuation)
            # stopping occurs before end, so stopping makes the tail constant anyway.
            k = max(k, total + F(c[-1], end+1))
        b = [c[min(r, tau)] if tau is not None else c[r] for r in range(end+1)]
        assert all(0 <= b[r] <= r for r in range(3, end+1))
        tail = [F(0)] * (end + 2)
        tail[end+1] = F(b[end], end+1)
        for r in range(end, 2, -1):
            tail[r] = F(b[r], r * (r+1)) + tail[r+1]
        assert 0 <= tail[3] <= k + 1
        ms = {3, 10, end - 8}
        if tau is not None:
            ms.update({max(3, tau-1), tau, min(end, tau+1)})
            seen.add('stopped')
            if tau == 3:
                seen.add('tau=3')
        else:
            seen.add('unstopped')
        for m in ms:
            assert 0 <= tail[m+1] <= k + 1
            if tau == m:
                seen.add('tau=m')
            if tau == m+1:
                seen.add('tau=m+1')
            for n in {m, min(end, m+1), end}:
                inc = sum(F(b[r]-b[r-1], r) for r in range(m+1, n+1))
                assert inc == F(b[n], n) - F(b[m], m+1) + tail[m+1] - tail[n]
                paper_rhs = F(b[n], n) - F(b[m], m+1) + sum(
                    F(b[r], r*(r+1)) for r in range(m+1, n))
                if n > m:
                    assert inc == paper_rhs
                elif inc != paper_rhs:
                    caught['apply printed Abel equality at n=m'] = True
                assert inc <= F(b[n], n) + tail[m+1]
                assert inc <= tail[m+1]  # monotone continuation gives tail[n] >= B_n/n
                if inc != -F(b[m], m+1) + tail[m+1] - tail[n]:
                    caught['omit B_n/n from exact parts'] = True
                if inc != F(b[n], n) + tail[m+1] - tail[n]:
                    caught['omit lower endpoint in exact parts'] = True
            if tau is not None and tau <= m:
                assert tail[m+1] == F(c[tau], m+1)
        if any(c[r] - c[r-1] > 1 for r in range(4, end+1)):
            caught['accept multi-promotion as physical'] = True
        count += 1
    assert seen == {'tau=3', 'tau=m', 'tau=m+1', 'stopped', 'unstopped'}
    print(f'Stopped prefixes: {count} exact paths; {", ".join(sorted(seen))}.')


if __name__ == '__main__':
    negatives = tree_checks()
    stopped_checks(negatives)
    assert len(negatives) == 9, negatives
    print(f'Caught all {len(negatives)} negatives: {", ".join(sorted(negatives))}.')
