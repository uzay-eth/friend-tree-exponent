"""Exact generator and nested incidence checks on finitely planted physical trees."""

from bisect import bisect_left
from collections import defaultdict
from fractions import Fraction as Q
from itertools import product
from math import comb, prod

from chains import ChainModel
from model import Cell


class Planted:
    def __init__(self, parents, leaves, external_degree=23):
        self.p = tuple(parents)
        self.leaves = tuple(leaves)
        self.external_degree = external_degree
        self.children = tuple(tuple(j for j, parent in enumerate(parents) if parent == i)
                              for i in range(len(parents)))
        self.d = tuple(1 + l + len(kids) for l, kids in zip(leaves, self.children))
        assert self.p[0] == -1 and min(self.d) >= 2
        self.edges = [(p, j) for j, p in enumerate(self.p) if p >= 0]
        self.chains = [(self.p[p], p, j) for p, j in self.edges if self.p[p] >= 0]

    def jumps(self):
        gains = [Q(l) for l in self.leaves]  # All leaf-to-core arrows.
        gains[0] += Q(1, self.external_degree)
        for p, j in self.edges:
            gains[p] += Q(1, self.d[j])
            gains[j] += Q(1, self.d[p])
        for i, gain in enumerate(gains):
            leaves = list(self.leaves)
            if gain:
                leaves[i] += 1
                yield gain, Planted(self.p, leaves, self.external_degree)
                leaves[i] -= 1
            if leaves[i]:
                rate = Q(leaves[i], self.d[i])
                leaves[i] -= 1
                yield rate, Planted(self.p + (i,), leaves + [1], self.external_degree)

    def roles(self, geom):
        for j in geom:
            yield self.d[j], len(self.children[j])

    def cell(self, model, degree, arity):
        if degree < model.threshold:
            return Cell(degree, arity)
        if arity <= model.arity:
            y = Q(model.threshold, degree)
            # The degree at the interface belongs to the lower reciprocal interval.
            return Cell(0, arity, max(0, bisect_left(model.boundaries, y)-1))
        return None

    def observations(self, model):
        result = defaultdict(Q)
        for geom in self.chains:
            roles = tuple(self.roles(geom))
            cells = tuple(self.cell(model, *r) for r in roles)
            if None in cells or (cells, "mass") not in model.blocks:
                continue
            labels = tuple(c.label for c in cells)
            for e in model.tests(labels):
                result[labels, e] += prod(Q(model.threshold, d)**power
                                           for (d, _), power in zip(roles, e))
        return result

    def vector(self, model, discount):
        result = defaultdict(Q)
        split = ((model.high-discount)/(model.high-model.low),
                 (discount-model.low)/(model.high-model.low))
        for geom in self.edges + self.chains:
            roles = tuple(self.roles(geom))
            cells = tuple(self.cell(model, *r) for r in roles)
            if None in cells or (cells, "mass") not in model.blocks:
                continue
            axes = []
            for (degree, _), c in zip(roles, cells):
                if c.degree:
                    axes.append((Q(1),))
                    continue
                lo, hi = model.boundaries[c.slab:c.slab+2]
                x = (Q(model.threshold, degree)-lo)/(hi-lo)
                W = (1 + Q(1, degree))**model.m
                axes.append(tuple(comb(model.n, j)*x**j*(1-x)**(model.n-j)/W
                                  for j in range(model.n+1)))
            rates = {}
            for i, j in enumerate(geom):
                if i == 0:
                    mother = self.p[j]
                    rates[i, "mother"] = Q(1, self.d[mother] if mother >= 0 else self.external_degree)
                others = [c for c in self.children[j] if c not in geom]
                if others:
                    rates[i, "children"] = sum(Q(1, self.d[c]) for c in others)
            for idx in product(*(range(len(axis)) for axis in axes)):
                weight = prod(axis[j] for axis, j in zip(axes, idx))
                if not weight:
                    continue
                for col, fraction in zip(model.blocks[cells, "mass"][idx], split):
                    result[col] += weight*fraction
                for kind, value in rates.items():
                    result[model.blocks[cells, kind][idx][0]] += weight*value
        return result


def dot(row, vector):
    return sum((v*vector[i] for i, v in row.items() if i in vector), Q(0))


def check():
    opts = dict(threshold=4, arity=2, order=2, inner=3, wide_roles=(0, 1, 2),
                chain_arity=(2, 2, 2), max_wide=1, breaks=(Q(2, 9),))
    base = ChainModel(**opts, project=False)
    model = ChainModel(**opts)
    assert model.variables == base.variables
    projections = model.upper_rows[len(base.upper_rows):]
    trees = [
        ([-1,0,1], [1,1,1]),                 # a first surviving two-generation chain
        ([-1,0,0,1,1,3,4], [0,1,1,0,1,2,1]),
        ([-1,0,0,1,2,2,3,5], [15,0,0,1,1,2,2,1]),  # exact D=18 at a reciprocal cut
        ([-1,0,0,1,1,2,3,5], [16,0,2,1,5,0,1,1]),
        ([-1,0,1,1,2,3,2,6], [0,14,1,0,1,1,1,3]),
        ([-1,0,0,0,1,2,3,4,7], [14,0,0,2,1,1,1,4,1]),  # unresolved high mother, k=3
        ([-1,0,1,2,2,2,3,4,4], [17,0,15,1,1,1,1,1,1]),
        ([-1,0,1,2,3,4,4,5], [0,1,0,14,0,16,3,1]),
    ]
    balances = [item for item in model.balance.items() if len(item[0][0]) == 3]
    projection_checks = balance_checks = nonzero = 0
    entered = set()
    for parents, leaves in trees:
        tree = Planted(parents, leaves)
        values = tree.observations(model)
        generator = defaultdict(Q)
        for rate, target in tree.jumps():
            after = target.observations(model)
            for key in values.keys() | after.keys():
                generator[key] += rate*(after[key]-values[key])
                if not values.get(key) and after.get(key):
                    entered.add(key)
        for discount in (model.low, (2*model.low+model.high)/3):
            v = tree.vector(model, discount)
            for key, row in balances:
                expected = discount*values[key] - generator[key]
                assert dot(row, v) == expected, (tree.p, leaves, key, expected, dot(row,v))
                nonzero += bool(expected)
                balance_checks += 1
            for row, rhs in projections:
                a = dot(row, v)
                assert a <= rhs, (tree.p, leaves, a, rhs)
            projection_checks += len(projections)
    assert nonzero and entered and len(projections) > 0
    print("PASS", balance_checks, "chain balances,", nonzero, "nonzero,",
          len(entered), "distinct entrance tests")
    print("PASS", projection_checks, "reciprocal and persistent receive projections",
          "including degree 18 and out-of-census mothers")


if __name__ == "__main__":
    check()
