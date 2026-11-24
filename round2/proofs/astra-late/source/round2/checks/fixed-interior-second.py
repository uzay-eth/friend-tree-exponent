"""Independent exact adverse fixtures for the fixed-interior proof."""

from fractions import Fraction as Q
from itertools import product


def tree(edges):
    g = {}
    for x, y in edges:
        g.setdefault(x, set()).add(y)
        g.setdefault(y, set()).add(x)
    return g


def grow(g, target):
    return tree([(x, y) for x in g for y in g[x] if x < y] +
                [(target, max(g) + 1)])


def outside(g, tag):
    """Literal arrow list; root-to-tag and stock-to-tag are excluded."""
    return {(x, y): Q(1, len(g[x])) for x in g for y in g[x] if y != tag}


def inert_plus_tag(g, tag):
    # Compute the proposed active channels from nontrivial components,
    # retaining each removed edge as an exterior slot at its root.
    active = {}
    for root in g[tag]:
        if len(g[root]) == 1:
            continue
        todo = [(root, tag)]
        while todo:
            x, previous = todo.pop()
            interior = g[x] - {tag}
            degree = len(interior) + (x == root)
            for y in interior:
                active[x, y] = Q(1, degree)
                if y != previous:
                    todo.append((y, x))
    active.update({(tag, y): Q(1, len(g[tag])) for y in g[tag]})
    return active


def clock_controls():
    count = arrows = prefix = 0
    for width, depth, offset in product((1, 2, 7, 31, 89), (0, 1, 5), (0, 1, 4, 19)):
        edges = [(0, 1)] + [(1, 2 + i) for i in range(width)]
        end = 2
        for _ in range(depth):
            new = max(max(e) for e in edges) + 1
            edges.append((end, new))
            end = new
        base = tree(edges)
        tag = max(base) + 1
        # 0 was a leaf on the old side, with post-birth degree two.
        g0 = grow(base, 0)
        b = Q(len(g0)) - Q(1, 2)
        g = g0
        pre = outside(g0, tag)
        for _ in range(offset):
            g = grow(g, tag)
        post = outside(g, tag)
        assert {e: w for e, w in pre.items() if e[0] != tag} == {
            e: w for e, w in post.items() if e[0] != tag}
        assert sum(post.values()) == b
        prefix += 1
        defects = [(0, 1), (1, 0), (tag, 0)]
        if offset:
            defects += [(tag, tag + 1)]
        states = [(g, 0)] + [(grow(g, y), 1) for x, y in defects if (x, y) in post]
        for h, f in states:
            literal = outside(h, tag)
            assert literal == inert_plus_tag(h, tag)
            assert all((root, tag) not in literal for root in h[tag])
            gain = sum(Q(1, len(h[y])) for y in h[tag])
            r = len(h[tag]) - Q(1, 2)
            assert len(h) - r == b + f
            assert sum(literal.values()) == b + f + r - gain
            assert sum(w for (x, _), w in literal.items() if x == tag) == 1
            assert (r - gain) >= 0
            count += 1
            arrows += len(literal)
    print(f'PASS {count} auxiliary-clock states, {arrows} active arrows, {prefix} exact prefixes')


def beta_monomial(a, b, red, black):
    ans = Q(1)
    for i in range(red):
        ans *= (a + i) / (a + b + i)
    for i in range(black):
        ans *= (b + i) / (a + b + red + i)
    return ans


def tilt_controls():
    count = 0
    for a, v, k in product((Q(1, 2), Q(1, 3), Q(1, 31)), (5, 29, 97), (0, 1, 4, 12)):
        b = v - a
        for word in product((0, 1), repeat=5):
            red, black = sum(word), len(word) - sum(word)
            # Exact joint reference weight of prefix and suffix under the
            # endpoint mixture, and after posterior resampling.
            prefix = beta_monomial(a, b, k, 0)
            assert beta_monomial(a, b, k + red, black) == (
                prefix * beta_monomial(a + k, b, red, black))
            count += 1
    print(f'PASS {count} prefix/suffix posterior moment identities')


def kernel_controls():
    count = 0
    for theta, y, delta in product((Q(1, 101), Q(1, 3), Q(5, 7), Q(100, 101)),
                                   (Q(1, 120), Q(2, 3), Q(5), Q(121)),
                                   (Q(0), Q(1, 6), Q(7, 2))):
        s = y / (1 - theta)
        log_drift = delta * ((1 - theta) / (1 + (1 - theta) * s) - 1 / s)
        assert log_drift / (1 - theta) == -delta / (y * (1 + y))
        assert log_drift * (1 + y) / (1 - theta) == -delta / y
        if delta:
            wrong = -delta / (s * (1 + s))  # whole-tree clock kernel in s
            assert wrong != log_drift
        count += 1
    print(f'PASS {count} signed kernel identities, distinguishing the whole-tree clock')


if __name__ == '__main__':
    clock_controls()
    tilt_controls()
    kernel_controls()
