"""Build reachable low-component matrices and exact Perron certificates (standard library)."""

import argparse
import json
from collections import defaultdict
from fractions import Fraction as F
from math import ceil, floor, lcm


HUB = 'hub'
FOUNDER = (2, 2, ())
NEWBORN = (2, 1, ())


def decode(code):
    degree, leaves, adjacent = [], [], []

    def visit(item, parent):
        total, leaf_count, children = item
        i = len(degree)
        degree.append(total)
        leaves.append(leaf_count)
        adjacent.append(set())
        if parent is not None:
            adjacent[i].add(parent)
            adjacent[parent].add(i)
        for child in children:
            visit(child, i)
    visit(code, None)
    return degree, leaves, adjacent


def encode(degree, leaves, adjacent, vertices):
    def rooted(i, parent):
        children = tuple(sorted(rooted(j, i) for j in adjacent[i]
                                if j != parent and j in vertices))
        return degree[i], leaves[i], children

    return min(rooted(i, None) for i in vertices)


def branch(adjacent, start, excluded):
    found = set()
    queue = [start]
    while queue:
        i = queue.pop()
        if i == excluded or i in found:
            continue
        found.add(i)
        queue.extend(adjacent[i])
    return found


def outcomes(code, cutoff, cap):
    if code == HUB:
        yield F(1), (HUB, NEWBORN)
        return
    degree, leaves, adjacent = decode(code)
    vertices = set(range(len(degree)))
    for i in sorted(vertices):
        gain = leaves[i] + sum(F(1, degree[j]) for j in adjacent[i])
        if gain:
            if degree[i] + 1 == cutoff:
                pieces = [branch(adjacent, j, i) for j in sorted(adjacent[i])]
                result = (HUB,) + tuple(encode(degree, leaves, adjacent, p) for p in pieces)
            else:
                new_degrees = degree.copy()
                new_leaves = leaves.copy()
                new_degrees[i] += 1
                new_leaves[i] += 1
                result = (encode(new_degrees, new_leaves, adjacent, vertices),)
            yield gain, result

        if not leaves[i]:
            continue
        bud = F(leaves[i], degree[i])
        new_degrees = degree + [2]
        new_leaves = leaves + [1]
        new_leaves[i] -= 1
        new_adjacent = [ns.copy() for ns in adjacent] + [{i}]
        new_adjacent[i].add(len(degree))
        new_vertices = vertices | {len(degree)}
        if len(new_vertices) <= cap:
            yield bud, (encode(new_degrees, new_leaves, new_adjacent, new_vertices),)
            continue
        minimum_side = (cap + cutoff - 2)//(cutoff - 1)
        cuts = []
        for u in sorted(new_vertices):
            for v in sorted(new_adjacent[u]):
                if u < v:
                    a = branch(new_adjacent, u, v)
                    b = new_vertices - a
                    if min(len(a), len(b)) >= minimum_side:
                        cuts.append((a, b))
        assert cuts
        for a, b in cuts:
            yield bud/len(cuts), (encode(new_degrees, new_leaves, new_adjacent, a),
                                  encode(new_degrees, new_leaves, new_adjacent, b))


def reachable_matrix(cutoff, cap):
    codes = [FOUNDER]
    index = {FOUNDER: 0}
    rows = []
    for i, code in enumerate(codes):
        row = defaultdict(F)
        for rate, result in outcomes(code, cutoff, cap):
            row[i] -= rate
            for child in result:
                if child not in index:
                    index[child] = len(codes)
                    codes.append(child)
                row[index[child]] += rate
        rows.append({j: q for j, q in row.items() if q})
    return codes, rows


def classes(rows):
    # Tarjan on the off-diagonal positive graph.
    stack, active, order, earliest, result = [], set(), {}, {}, []

    def search(i):
        order[i] = earliest[i] = len(order)
        stack.append(i)
        active.add(i)
        for j, entry in rows[i].items():
            if i == j or entry <= 0:
                continue
            if j not in order:
                search(j)
                earliest[i] = min(earliest[i], earliest[j])
            elif j in active:
                earliest[i] = min(earliest[i], order[j])
        if earliest[i] == order[i]:
            component = []
            while True:
                j = stack.pop()
                active.remove(j)
                component.append(j)
                if j == i:
                    break
            result.append(sorted(component))
    for i in range(len(rows)):
        if i not in order:
            search(i)
    return sorted(result)


def integer_block(rows, members):
    index = {old: i for i, old in enumerate(members)}
    local = [{index[j]: x for j, x in rows[i].items() if j in index} for i in members]
    denominator = lcm(*(x.denominator for row in local for x in row.values()))
    integer = [{j: int(x*denominator) for j, x in row.items()} for row in local]
    return denominator, integer


def certify_block(rows, members):
    denominator, integer = integer_block(rows, members)
    n = len(members)
    if n == 1:
        q = F(integer[0].get(0, 0), denominator)
        return q, q, [1]
    shift = 1 + max(0, *(-row.get(i, 0) for i, row in enumerate(integer)))
    shifted = []
    for i, row in enumerate(integer):
        row = row.copy()
        row[i] = row.get(i, 0) + shift
        assert min(row.values()) > 0
        shifted.append(row)
    scale = 10**30
    weights = [scale]*n
    for k in range(40000):
        nxt = [sum(x*weights[j] for j, x in row.items()) for row in shifted]
        if k % 16 == 0:
            ratios = [F(y-shift*w, denominator*w) for y, w in zip(nxt, weights)]
            lo, hi = min(ratios), max(ratios)
            if hi-lo < F(1, 10**16):
                return lo, hi, weights
        high = max(nxt)
        weights = [max(1, (y*scale + high//2)//high) for y in nxt]
    raise RuntimeError(f'power iteration failed on a class of size {n}')


def char_oracle(rows, members):
    # Newton's identities on the integer matrix L*G. Encode coefficients for
    # det(X I-L*G); the eigenvalue lambda corresponds to X=L*lambda.
    denominator, a = integer_block(rows, members)
    n = len(members)
    power = [[int(i == j) for j in range(n)] for i in range(n)]
    traces = []
    coefficients = [1]
    for k in range(1, n+1):
        nxt = [[0]*n for _ in range(n)]
        for i, row in enumerate(a):
            for j, value in row.items():
                for col, x in enumerate(power[j]):
                    nxt[i][col] += value*x
        power = nxt
        traces.append(sum(power[i][i] for i in range(n)))
        value = -sum(coefficients[k-i]*traces[i-1] for i in range(1, k+1))
        assert value % k == 0
        coefficients.append(value//k)
    return denominator, coefficients


def build_case(cutoff, cap):
    codes, rows = reachable_matrix(cutoff, cap)
    group_records = []
    bounds = []
    for group in classes(rows):
        lo, hi, weights = certify_block(rows, group)
        bounds.append((lo, hi))
        group_records.append({'indices': group, 'weights': [str(w) for w in weights]})
    lo = max(x for x, _ in bounds)
    hi = max(y for _, y in bounds)
    maximizing = [i for i, (_, y) in enumerate(bounds) if y >= lo]
    if len(maximizing) != 1:
        raise RuntimeError(f'Perron class unresolved for {cutoff},{cap}: {maximizing}')
    winner = maximizing[0]
    denominator, polynomial = char_oracle(rows, group_records[winner]['indices'])
    decimal_scale = 10**14
    lower, upper = floor(decimal_scale*lo), ceil(decimal_scale*hi)
    print(f'M={cutoff} B={cap} types={len(codes)} classes={len(bounds)} '
          f'Perron_class={len(group_records[winner]["indices"])} '
          f'lambda in [{lower}/{decimal_scale},{upper}/{decimal_scale}]')
    return {'M': cutoff, 'B': cap, 'types': codes,
            'rows': [[[j, x.numerator, x.denominator] for j, x in sorted(r.items())]
                     for r in rows],
            'classes': group_records, 'winner': winner,
            'char_scale': denominator, 'char_coefficients': [str(x) for x in polynomial],
            'lower_numerator': lower, 'upper_numerator': upper,
            'bracket_denominator': decimal_scale}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--case', action='append', required=True, metavar='M,B')
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    cases = []
    for setting in args.case:
        cutoff, cap = map(int, setting.split(','))
        if cutoff < 3 or cap < 1:
            parser.error('each case needs M >= 3 and B >= 1')
        cases.append(build_case(cutoff, cap))
    with open(args.out, 'w', encoding='utf-8') as destination:
        json.dump(cases, destination, separators=(',', ':'))
        destination.write('\n')


if __name__ == '__main__':
    main()
