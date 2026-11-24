#!/usr/bin/env python3
"""Generate a standalone Fraction checker from the prior checker and audit delta."""
from pathlib import Path

here = Path(__file__).resolve().parent
s = (here.parent / 'mixed_l5_study/check_coarse.py').read_text()


def sub(a, b):
    global s
    assert s.count(a) == 1, a
    s = s.replace(a, b)


sub('from collections import Counter', 'from argparse import ArgumentParser\nfrom collections import Counter')
sub('D5 = DC = 6', 'D5 = DC = 6\nBUCKETS = 2\nFAULT = None')
sub('BIG = (-1, 0, 0)', '''BIG = (-1, 0, 0)
TAIL = (-2, 0, 0)
LARGE = (BIG, TAIL)''')
sub('    return Q(crc32(repr((s,d,apex)).encode()) % 37 - 18, 17)',
    '''    # Separate the fixture's tail factors from the inner envelope so a
    # wrong-bucket lookup is observable (the other factors remain signed).
    bias = 4 if s[0] < 0 else 0
    return bias + Q(crc32(repr((s,d,apex)).encode()) % 37 - 18, 17)''')
sub('''    a = tuple(s if sum(s) < DC else BIG for s in st)
    if any(s == BIG and d != 1 for s,d in zip(a,degrees(p))):''', '''    a = tuple(s if sum(s) < DC else TAIL if BUCKETS == 2 and sum(s) >= D0 else BIG for s in st)
    if any(s in LARGE and d != 1 for s,d in zip(a,degrees(p))):''')
sub('''def fiber(is_apex):
    # All physical inner BIG controls plus all appropriate tail controls. For
    # tail controls include infeasible corners too (a conservative superset).
    states = {INF}''', '''def fiber(is_apex, large):
    # T2/T3 and the shared INF are separate controls even on overflow.
    inner = large == BIG
    tail = BUCKETS == 1 or large == TAIL
    if FAULT == 'wrong-fiber' and large == TAIL:
        inner, tail = True, False
    states = {INF} if tail else set()''')
sub('''            states.update((l,k,h) for d in range(DC,D0) if (l:=d-k-h)>=0)
            states.update((i,k,h) for i in (-2,-3,-4))''', '''            if inner:
                states.update((l,k,h) for d in range(DC,D0) if (l:=d-k-h)>=0)
            if tail:
                states.update((i,k,h) for i in (-2,-3,-4))''')
sub('if BIG not in a:', 'if not any(s in LARGE for s in a):')
sub('if s!=BIG:', 'if s not in LARGE:')
sub('factor_bounds(ds[i],i==0)', 'factor_bounds(ds[i],i==0,s)')
sub('''def factor_bounds(d, apex):
    vals = [factor(s,d,apex) for s in fiber(apex)]''', '''def factor_bounds(d, apex, large):
    vals = [factor(s,d,apex) for s in fiber(apex,large)]''')
sub('project(p,s).count(BIG)>0', 'any(v in LARGE for v in project(p,s))')
sub('old_by_bud[u,n] = Q(st[u][0],sum(st[u]))*q5(pa,a)', '''if FAULT != 'drop-tail-birth' or TAIL not in project(pa,a):
                        old_by_bud[u,n] = Q(st[u][0],sum(st[u]))*q5(pa,a)
                elif FAULT == 'allow-overflow-birth' and st[u][1] > KT:
                    n = m+(len(adj),)
                    old_by_bud[u,n] = Q(st[u][0],sum(st[u]))''')
sub('assert x*f0 <= exact <= x*f1', "assert x*f0 <= exact <= x*f1, 'fiber source'")
sub('''            if a[i] == BIG:
                assert not g and not b
                counts['big_stays_big'] += 1
                continue
            r = Q(st[u][0]) + sum(Q(1,sum(st[w])) for w in adj[u] if w in n and a[n.index(w)]!=BIG)
            slack = Q(sum(w in n and a[n.index(w)]==BIG for w in adj[u]),DC) + Q(len(adj[u].difference(n)),2)
            rr = Q(st[u][0]) + sum(Q(1,sum(st[w])) for w in adj[u])''', '''            rr = Q(st[u][0]) + sum(Q(1,sum(st[w])) for w in adj[u])
            if a[i] in LARGE:
                assert not b
                if a[i] == TAIL or BUCKETS == 1:
                    assert not g
                    counts['tail_stays_tail' if a[i] == TAIL else 'big_stays_big'] += 1
                    continue
                dest = a[:i]+(TAIL,)+a[i+1:]
                jump = row(p,dest)-v
                crosses = sum(st[u]) == D0-1
                assert g == (jump if crosses else ZERO)
                exact = rr*g
                corner = (D0-1 if FAULT != 'drop-crossing' else 0)*jump
                low, high = min(ZERO,corner), max(ZERO,corner)
                assert low <= exact <= high, 'boundary crossing'
                if crosses:
                    counts['crossings'] += 1
                    counts['crossing_negative' if g < 0 else 'crossing_positive' if g > 0 else 'crossing_zero'] += 1
                else:
                    counts['inner_stays_inner'] += 1
                l5actual += exact
                l5lo += low
                l5hi += high
                continue
            r = Q(st[u][0]) + sum(Q(1,sum(st[w])) for w in adj[u] if w in n and a[n.index(w)] not in LARGE)
            caps = {BIG: Q(1,DC), TAIL: ZERO if FAULT == 'zero-tail-credit' else Q(1,D0)}
            slack = sum(caps[a[n.index(w)]] for w in adj[u] if w in n and a[n.index(w)] in LARGE) + Q(len(adj[u].difference(n)),2)''')
sub('assert low<=exact<=high', "assert low<=exact<=high, 'gain rate interval'")
sub("                counts['birth_with_big'] += BIG in project(p,x)", "                counts['birth_with_big'] += BIG in project(p,x)\n                counts['birth_with_tail'] += TAIL in project(p,x)")
sub('    assert born == old_by_bud', "    assert born == old_by_bud, 'newborn ownership'")
sub("if __name__ == '__main__':\n    rng", "def audit():\n    rng")
sub("        if rep%80<2:\n            leaves[1 if rep%80==0 else 0] = 3", '''        if rep%80<2:
            leaves[1 if rep%80==0 else 0] = 3
        elif rep%10 == 2:
            u = rng.randrange(n)
            leaves[u] = D0-1-len(a[u])''')
sub("    for k in ('strict_source','multiple_deletions','mixed_credits','birth_with_big','overflow_birth_excluded','gain_crosses_dc'):\n        assert counts[k], k\n    print('Exact coarse tree inequalities PASS:', ', '.join(f'{k}={v}' for k,v in sorted(counts.items())))", '''    required = ['strict_source','multiple_deletions','mixed_credits','birth_with_big','overflow_birth_excluded','gain_crosses_dc']
    if BUCKETS == 2:
        required += ['crossing_negative','crossing_positive','inner_stays_inner','tail_stays_tail','birth_with_tail']
    for k in required:
        assert counts[k], k
    return counts


if __name__ == '__main__':
    parser = ArgumentParser(description=__doc__)
    parser.add_argument('--d', type=int, choices=(6,7), required=True)
    parser.add_argument('--buckets', type=int, choices=(1,2), required=True)
    args = parser.parse_args()
    D5 = DC = args.d
    BUCKETS = args.buckets
    counts = audit()
    print(f'D={DC}, buckets={BUCKETS}: Exact tree inequalities PASS:', ', '.join(f'{k}={v}' for k,v in sorted(counts.items())), flush=True)
    faults = {'allow-overflow-birth': 'newborn ownership'}
    if BUCKETS == 2:
        faults.update({'drop-crossing': 'boundary crossing', 'wrong-fiber': 'fiber source',
                       'zero-tail-credit': 'gain rate interval', 'drop-tail-birth': 'newborn ownership'})
    for FAULT, expected in faults.items():
        for fn in (factor,row,q4,fiber,extrema,factor_bounds):
            fn.cache_clear()
        try:
            audit()
        except AssertionError as e:
            assert str(e) == expected, (FAULT, expected, str(e))
            print(f'planted {FAULT}: rejected ({e})', flush=True)
        else:
            raise AssertionError(f'planted {FAULT} survived')''')
(here / 'check_trees.py').write_text(s)
