"""Small labelled oracle and provenance pins for the early extension."""

from collections import defaultdict
from fractions import Fraction as F
from functools import cache
from hashlib import sha256
from pathlib import Path
from runpy import run_path


ROUND = Path(__file__).resolve().parents[1]
EARLY = ROUND / 'open/age-order/early'
PINS = {
    'check.py': '6864f1f5fa2ffc667ce70f1b8decd779cd4964faaa146872f9398ffd5da6f1e9',
    'further/check.py': '0ba2167d6d6b8b0b3bd7e8cc02c93822ee117e6810b1a7e78fc887753e84ed06',
    'further/README.md': '4073f8857583d82333a752ac30d8854d7cf34653b8896dd495170a08df4e1886',
}
for path, digest in PINS.items():
    assert sha256((EARLY / path).read_bytes()).hexdigest() == digest, path
W = run_path(str(EARLY / 'further/check.py'))
S = run_path(str(Path(__file__).with_name('early-age-sol.py')))
EDGE, attach, arrows, hub, core = (S[key] for key in ('EDGE', 'add', 'arrows', 'hub', 'core'))
U = 10**12


def distance(g, a, b):
    todo = [(a, -1, 0)]
    while todo:
        x, p, d = todo.pop()
        if x == b:
            return d
        todo.extend((y, x, d + 1) for y in g[x] if y != p)
    raise AssertionError('disconnected')


def rooted(g, at, parent=-1, tags=()):
    return (tuple(at == t for t in tags),
            tuple(sorted(rooted(g, i, at, tags) for i in g[at] if i != parent)))


@cache
def key(g, *tags):
    # Center-rooted canon with explicit colors, rather than rooting at tag 1.
    eccentric = [max(distance(g, i, j) for j in range(len(g))) for i in range(len(g))]
    return min(rooted(g, i, tags=tags) for i, e in enumerate(eccentric) if e == min(eccentric))


def down(x):
    return x.numerator // x.denominator


@cache
def absolute(g, v):
    lo, hi = hub(g, v)
    if len(g[v]) == 1:
        d = len(g[g[v][0]])
        lo = core(2, 1) / (1 + d*d)
    return down(lo*U), -down(-hi*U)


@cache
def law(g, u, v, h):
    if key(g,u,v) == key(g,v,u):
        return 0, 0
    a,b = absolute(g,u)
    c,d = absolute(g,v)
    lo,hi = a-d,b-c
    if h:
        successors = [law(attach(g,t),u,v,h-1) for t in range(len(g))]
        x,y = (sum(p*z[i] for p,z in zip(arrows(g),successors)) for i in (0,1))
        lo,hi = max(lo,down(x)),min(hi,-down(-y))
    assert lo <= hi
    return lo,hi


def canonical_multiset(moves):
    out = defaultdict(F)
    for dest, p in moves:
        k = key(dest)
        out[k] += p
    return out


def placement_coeffs(g, ordered):
    n = len(g)
    out = defaultdict(F)
    representatives = {}
    for i,j,p in ordered:
        dest = attach(attach(g,i),j)
        a,b = key(dest,n,n+1),key(dest,n+1,n)
        if a == b:
            continue
        if a < b:
            out[a] += p
            representatives[a] = (dest,n,n+1)
        else:
            out[b] -= p
            representatives[b] = (dest,n+1,n)
    return {k:w for k,w in out.items() if w},representatives


def check_gap(g, depths):
    n = len(g)
    r = [n*p for p in arrows(g)]
    q = [a/(len(ns)*(len(ns)+1)) for a,ns in zip(r,g)]
    den = n*(n+1)
    paired = [(i,n,r[i]/((len(g[i])+1)*den)) for i in range(n)]
    paired += [(i,j,(q[j]-q[i])/den) for i,ns in enumerate(g) for j in ns if i<j]
    brute = [(i,j,p*s) for i,p in enumerate(arrows(g))
             for j,s in enumerate(arrows(attach(g,i)))]
    assert sum(p for _,_,p in brute) == 1
    assert placement_coeffs(g,brute)[0] == placement_coeffs(g,paired)[0]
    for h in depths:
        lo = hi = F(0)
        for i,j,p in paired:
            if not p:
                continue
            dest = attach(attach(g,i),j)
            got = law(dest,n,n+1,h)
            assert got == W['bound'](W['encode'](dest,n,n+1),h)
            a,b = sorted(p*F(x,U) for x in got)
            lo,hi = lo+a,hi+b
        assert (lo,hi) == W['compare']([list(x) for x in g],h)[0]


def main():
    # Positive shape masses and target multiplicities without orbit division.
    # Concrete star fixture: center / individual leaf weights 3/4 and 1/12;
    # the three leaves give one output orbit of total probability 1/4.
    star = attach(attach(EDGE,0),0)
    assert arrows(star) == (F(3,4),F(1,12),F(1,12),F(1,12))
    assert sorted(canonical_multiset((attach(star,i),p) for i,p in enumerate(arrows(star))).values()) == [F(1,4),F(3,4)]
    # Two T3 completions lead to the same *unmarked* fork, but their tagged
    # contrasts cannot cancel: one compares the long-arm middle with the tip,
    # the other compares a short-arm leaf with the long-arm tip.
    t3 = attach(EDGE,0)
    chain, edge = attach(attach(t3,0),3), attach(attach(t3,0),1)
    assert key(chain) == key(edge)
    assert {key(chain,3,4),key(chain,4,3)}.isdisjoint({key(edge,3,4),key(edge,4,3)})
    labelled = {EDGE:F(1)}
    total_shapes = recursions = endpoint_checks = 0
    for n,count in enumerate((1,1,2,3,6,11,23,47,106),2):
        independently = canonical_multiset(labelled.items())
        from_writer = {}
        for code,p in W['classes'](n).items():
            g,_ = W['decode'](code)
            k = key(tuple(map(tuple,g)))
            assert k not in from_writer
            from_writer[k] = p
        assert independently == from_writer and len(from_writer) == count
        total_shapes += count
        if n <= 7:
            tested = set()
            for g in labelled:
                if key(g) in tested:
                    continue
                tested.add(key(g))
                for i in range(n):
                    assert absolute(g,i) == W['endpoint'](g,i)
                    endpoint_checks += 1
                depths = (0,1,2) if n <= 5 else (0,)
                check_gap(g,depths)
                recursions += len(depths)
        if 8 <= n < 10:
            representatives = {key(g): g for g in labelled}
            labelled = {representatives[k]: p for k,p in independently.items()}
        if n < 10:
            transitions = defaultdict(F)
            for g,w in labelled.items():
                for i,p in enumerate(arrows(g)):
                    transitions[attach(g,i)] += w*p
            labelled = transitions
    print(f'{total_shapes} free-shape weights through size 10; labelled through size 8; star orbit 1/4,3/4')
    print(f'{endpoint_checks} physical endpoint rows and {recursions} exact whole-gap controls through size 7')
    print('fork negative control: same free shape, distinct tagged contrasts')


if __name__ == '__main__':
    main()
