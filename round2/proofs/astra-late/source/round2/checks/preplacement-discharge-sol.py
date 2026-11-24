"""Physical rooting fixtures and scalar controls for the hubhood discharge."""
from fractions import Fraction as Q
from pathlib import Path
from runpy import run_path

ROOT = Path(__file__).resolve().parents[1] / 'open' / 'age-order'
D = run_path(str(ROOT / 'global-transport/hostile/discharge/check.py'))
H = D['H']
G, U, L = (D[k] for k in ('G', 'U', 'L'))
T = D['TICKS']


def make(core, sizes):
    g = [[] for _ in sizes]
    for u, v in core:
        g[u].append(v)
        g[v].append(u)
    for u, d in enumerate(sizes):
        assert len(g[u]) <= d
        for _ in range(d - len(g[u])):
            g[u].append(len(g))
            g.append([u])
    return tuple(map(tuple, g))


def ceil(x):
    return -((-x)//1)


def control(g, stocks=False):
    degrees, w, s = H['profile'](g)
    cores = [i for i, d in enumerate(degrees) if d > 1]
    nbr = {i: [j for j in g[i] if degrees[j] > 1] for i in cores}
    a = [(d, tuple(sorted(degrees[j] for j in g[i] if degrees[j] > 1))) for i, d in enumerate(degrees)]
    credit = {}
    for u in cores:
        d, ds = a[u]
        # Actual endpoint scores computed from explicitly appended vertices.
        base = w[u]*G(d+1) + sum(w[j]*G(2) for j in g[u] if degrees[j] == 1)
        for j in g[u]:
            if degrees[j] != 1:
                continue
            tree = H['grow'](H['grow'](g, u), j)
            x, y = len(g), len(g)+1
            if s[u] < s[j]:
                x, y = y, x
            px, py = tree[x][0], tree[y][0]
            others = sum(len(tree[z]) == 1 for z in tree[px]) - 1
            base -= abs(s[u]-s[j])*max(0, U(len(tree[px]), others)-L(len(tree[py])))
        assert base == D['profile'](d, ds)[1]
        credit[u] = base

    def debit(a, b):
        d, ds = a
        e, es = b
        x, _, ux = D['profile'](d, ds)
        y, _, uy = D['profile'](e, es)
        if stocks:
            return max(0, (x-y)*(max(0, ux-L(e+1))), (y-x)*(max(0, uy-L(d+1))))
        # Independent reproduction of the pinned rounded debit on exact profiles.
        xt, yt = (T*x)//1, (T*y)//1
        aa, bb = ceil(T*max(0, ux-L(e+1))), ceil(T*max(0, uy-L(d+1)))
        return Q(ceil(Q(max(0, (xt-yt+1)*aa, (yt-xt+1)*bb), T)), T)

    def charge(a):
        d, ds = a
        c = D['profile'](d, ds)[1]
        return c if stocks else Q((T*c)//1, T)

    def msg(a, e):
        b0 = (e, (a[0],))
        return charge(b0)-debit(a, b0)

    negative = bare = count = adverse = 0
    for root in cores:
        def descend(u, p):
            nonlocal negative, bare, adverse
            val = charge(a[u]) + sum(descend(j, u) for j in nbr[u] if j != p)
            if p is not None:
                edge = debit(a[p], a[u])
                val -= edge
                adverse += edge > Q(1, T)
                m = msg(a[p], degrees[u])
                if not stocks:
                    assert val >= m
                negative += (m < 0)
                bare += (len(nbr[u]) == degrees[u])
            return val
        assert descend(root, None) > 0
        if not stocks:
            assert charge(a[root])+sum(msg(a[root], degrees[j]) for j in nbr[root]) >= Q(27404273, T)
        count += 1
    return count, negative, bare, adverse


def extension():
    for e in range(2, 10):
        for parent in (0, 10):
            s, credit, u = D['profile'](e, (parent,))
            assert credit-s*u > Q(3, 500) and credit-Q(1, 110) > Q(3, 500)
    assert Q(9, 11)*Q(16694235, 10**8)+Q(9, 20)*G(2)-Q(1, 22) > Q(3, 500)
    assert 8*Q(3, 500) > Q(1, 22)


if __name__ == '__main__':
    assert L(2) > Q(1, 4)
    assert 2*(Q(9, 11)*Q(16694235, 10**8)+Q(9, 20)*G(2))-Q(1, 110) > Q(22124577, T)
    assert all(H['score'](H['spider'](d, 1))[0] > 0 for d in (2, 3, 4, 9, 40))
    early = H['EARLY']
    assert early['adjacent']((0,), 4)[0] >= Q(22645238637, 10**12)
    # Bare degree-nine centre; then two adjacent bare high-arity centres with repeated slots.
    trees = [make([(0, j) for j in range(1, 10)], [9]+[2]*9),
             make([(0, 1)] + [(0, i) for i in range(2, 10)] + [(1, i) for i in range(10, 18)], [9, 9]+[2, 3, 9, 2]*4),
             make([(i, i+1) for i in range(8)], [2, 3, 4, 2, 7, 9, 2, 3, 2])]
    rows = [control(g) for g in trees]
    assert sum(r[3] for r in rows) > 0
    print(f'PASS {sum(r[0] for r in rows)} cap-nine rerootings; {sum(r[3] for r in rows)} adverse debits and {sum(r[2] for r in rows)} bare child checks')
    extension()
    old = make([(0, 1)] + [(0, i) for i in range(2, 11)] + [(1, i) for i in range(11, 23)], [10, 13]+[2]*19+[17, 42])
    assert control(old, True)[0] == 23
    print('PASS five stars, sixteen terminal endpoint bounds and 23 additional bare-skeleton rerootings')
