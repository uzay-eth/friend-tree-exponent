"""Independent tiny physical arrow and original-degree mask fixtures."""

from fractions import Fraction as Q


def graph(parents):
    result = [[] for _ in range(len(parents) + 1)]
    for child, parent in enumerate(parents, 1):
        result[child].append(parent)
        result[parent].append(child)
    return result


def rates(adj, slot=None):
    out = [Q(0)] * len(adj)
    for source, row in enumerate(adj):
        for target in row:
            out[target] += Q(1, len(row) + (source == slot))
    return out


def component(adj, top, cutoff):
    assert len(adj[top]) <= cutoff
    visited = {top}
    todo = [top]
    while todo:
        for target in adj[todo.pop()]:
            if len(adj[target]) <= cutoff and target not in visited:
                visited.add(target)
                todo.append(target)
    return visited


def main():
    birth = rates(graph((0,)), 0)
    assert birth == [1, Q(1, 2)]
    pair = 1 / sum(birth)
    star = birth[0] * pair / sum(rates(graph((0, 0)), 0))
    path = birth[1] * pair / sum(rates(graph((0, 1)), 0))
    assert (pair, star, path) == (Q(2, 3), Q(1, 4), Q(2, 15))
    print("inert pair/star/path occupations: 2/3, 1/4, 2/15")

    # 1 is bare between sources 0 and 2; each has a private leaf.
    base = (0, 1, 0, 2)
    rows = 0
    for cutoff in range(2, 22):
        for left in (cutoff - 1, cutoff, cutoff + 1):
            for right in (cutoff - 1, cutoff, cutoff + 1):
                adj = graph(base + (0,) * (left - 1) + (2,) * (right - 1))
                expected = left + 1 > cutoff and right + 1 > cutoff
                assert (component(adj, 1, cutoff) == {1}) == expected
                assert rates(adj)[1] == Q(1, left + 1) + Q(1, right + 1)
                rows += 1
    assert 1 - Q(4, 6) - 2 * Q(1, 2)**5 > 0
    print(f"{rows} original-degree masks and exact incoming rates on singleton plants")


if __name__ == "__main__":
    main()
