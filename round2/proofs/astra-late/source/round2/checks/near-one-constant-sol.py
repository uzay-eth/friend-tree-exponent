"""Independent exact first-defect checks and solvable urn normalization control."""
from fractions import Fraction as Q
from itertools import product
from math import comb, exp, factorial, gamma, isclose, lgamma, log


def tree(parents):
    adj = [[] for _ in range(len(parents) + 1)]
    for v, p in enumerate(parents, 1):
        adj[p].append(v)
        adj[v].append(p)
    return adj


def attach(adj, v):
    out = [a.copy() for a in adj]
    out[v].append(len(out))
    return out + [[v]]


def component(adj, tag, root):
    todo = [root]
    for x in todo:
        todo.extend(y for y in adj[x] if y != tag and y not in todo)
    return todo


def audit(adj, tag):
    n0, degree = len(adj), len(adj[tag])
    a = sum(Q(1, len(adj[w])) for w in adj[tag])
    beta = n0 - a
    old_cores = [w for w in adj[tag] if len(adj[w]) > 1]
    prior = Q(1)
    for k in range(7):
        n = len(adj)
        d = len(adj[tag])
        outside = [(x, y, Q(1, len(row))) for x, row in enumerate(adj)
                   for y in row if y != tag]
        tag_arrows = [(x, y, q) for x, y, q in outside if x == tag]
        intrinsic = [(x, y, q) for x, y, q in outside if x != tag]
        assert d == degree + k and n == n0 + k
        assert sum(q for _, _, q in outside) == beta
        assert sum(q for _, _, q in tag_arrows) == 1
        assert sum(q for _, _, q in intrinsic) == beta - 1
        assert all(y not in adj[tag] or y in old_cores for _, y, _ in intrinsic)
        assert len(tag_arrows) == d
        assert sum(q for _, y, q in tag_arrows if y in old_cores) == Q(len(old_cores), d)
        assert sum(q for _, y, q in tag_arrows if y not in old_cores) == 1 - Q(len(old_cores), d)
        for x, y, q in outside:
            # Every first defect is a genuine physical outside jump, even tag->core.
            defect = attach(adj, y)
            cores = [w for w in defect[tag] if len(defect[w]) > 1]
            parts = [(w, component(defect, tag, w)) for w in cores]
            assert 1 + sum(len(nodes) - 1 for _, nodes in parts) == n + 1 - d
            rate = sum(len(nodes) - Q(1, len(defect[w])) for w, nodes in parts)
            physical = sum(Q(1, len(row)) for row in defect for z in row if z != tag)
            assert 1 + rate == physical
            assert len(defect) - sum(Q(1, len(defect[w])) for w in defect[tag]) >= beta + 1
            if x != tag:
                assert set(cores) == set(old_cores)
            elif y not in old_cores:
                assert dict(parts)[y] == [y, n]
        next_prior = prior * (1 - beta / n)
        assert prior - next_prior == prior * beta / n
        prior = next_prior
        adj = attach(adj, tag)
    return 7


def gamma_control():
    # Solvable ordinary Pólya urn. Defect at m changes blue from b to b+1.
    # Postdefect m*blue_fraction -> Gamma(b+1), so inverse moment = 1/Gamma(b+1).
    # Check exact injection of conditional tails into the first-defect mixture.
    b, a, t = 3, 2, Q(1, 3)
    n0 = a + b

    def beta_cdf(p, q):
        degree = p + q - 1
        return sum(comb(degree, j) * t**j * (1 - t)**(degree - j)
                   for j in range(p, degree + 1))

    weight, mixture = Q(1), Q(0)
    for m in range(n0, 45):
        mixture += weight * Q(b, m) * beta_cdf(b + 1, m - b)
        weight *= 1 - Q(b, m)
        assert mixture + weight * beta_cdf(b, m + 1 - b) == beta_cdf(b, a)
    for a, n0 in [(1.0, 2), (0.5, 5), (1.0 / 3, 7), (1.75, 8)]:
        b = n0 - a
        m = 100_000
        scaled = exp(lgamma(m - b) - lgamma(m) + b * log(m))
        assert isclose(scaled, 1, rel_tol=4e-4)
    assert gamma(n0 := 5) / gamma(2) / gamma(4) == factorial(n0 - 1) / factorial(3)


if __name__ == '__main__':
    checks = 0
    for n in range(2, 7):
        for parents in product(*(range(i) for i in range(1, n))):
            for tag in range(n):
                checks += audit(tree(parents), tag)
    print(f'PASS {checks} prefix cases: directed Palm weights, all exit types, branch deficit and urn gap')
    gamma_control()
    print('PASS 40 exact urn-mixture truncations and four Gamma survival-ratio controls')
