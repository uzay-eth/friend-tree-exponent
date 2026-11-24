"""Stationary inverse-degree balances for oriented three-core paths."""

from collections import defaultdict
from fractions import Fraction as Q
from functools import cache
from itertools import product
from math import comb
from typing import NamedTuple

from model import Model


class Basis(NamedTuple):
    threshold: int
    m: int
    n: int
    slabs: int = 1


@cache
def restrict_controls(values, left, right):
    degree = len(values) - 1
    # The ith blossom has degree-i copies of left and i copies of right.
    return tuple(sum(values[a+b] * comb(degree-i, a) * left**a * (1-left)**(degree-i-a)
                     * comb(i, b) * right**b * (1-right)**(i-b)
                     for a in range(degree-i+1) for b in range(i+1))
                 for i in range(degree+1))


class ChainModel(Model):
    def __init__(self, threshold, arity, order, subdivisions=1, *, inner=3,
                 wide_roles=(0, 1, 2), max_wide=1, chain_arity=None, project=True,
                 port_transport=True, balances=True, breaks=(), **options):
        self.breaks = tuple(map(Q, breaks))
        if self.breaks:
            if not all(0 < x < 1 for x in self.breaks) or tuple(sorted(set(self.breaks))) != self.breaks:
                raise ValueError("reciprocal cut points must increase strictly inside (0,1)")
            if subdivisions != 1:
                raise ValueError("choose either uniform subdivisions or explicit cut points")
            subdivisions = len(self.breaks) + 1
        self.boundaries = ((Q(0),) + self.breaks + (Q(1),) if self.breaks
                           else tuple(Q(i, subdivisions) for i in range(subdivisions+1)))
        self.basis = Basis(threshold, order, 2*order+2)
        super().__init__(threshold, arity, order, subdivisions, **options)
        if not 2 <= inner < threshold or not 0 <= max_wide <= 3:
            raise ValueError("chain inner degree must be in [2,H), wide count in [0,3]")
        if not set(wide_roles) <= {0, 1, 2}:
            raise ValueError("roles are numbered from 0 to 2")
        if chain_arity is None:
            chain_arity = (max(arity, threshold - 2),) * 3
        if len(chain_arity) != 3 or min(chain_arity) < 0:
            raise ValueError("three nonnegative arity caps are required")
        allowed = []
        for i, cap in enumerate(chain_arity):
            allowed.append([c for c in self.cells if (i == 2 or c.arity) and c.arity <= cap
                            and (i in wide_roles or 0 < c.degree <= inner)])
        self.chains = [chain for chain in product(*allowed)
                       if sum(not c.degree or c.degree > inner for c in chain) <= max_wide]
        self.chain_labels = {tuple(c.label for c in chain) for chain in self.chains}
        for chain in self.chains:
            self.make_block(chain, "mass")
            for kind in self.ports(chain):
                self.make_block(chain, kind)
                self.box_port(chain, kind, self.slots(chain, kind))
        mass = defaultdict(Q, {self.total: -1})
        for chain in self.chains:
            self.add_observation(mass, chain, (0, 0, 0))
        self.add_ineq(mass)  # A terminal vertex has only one ordered ancestor pair.
        if balances:
            for chain in self.chains:
                self.chain_balance(chain)
            self.chain_entrances()
            for (labels, exponents), row in self.balance.items():
                if len(labels) == 3:
                    self.add_eq(row)
        if project:
            for pair in self.edges:
                self.connect(pair, prepend=True, transport=port_transport)
                if pair[1].arity:
                    self.connect(pair, prepend=False, transport=port_transport)

    @cache
    def axis(self, c, e, form="old", factor="one"):
        if c.degree or not self.breaks:
            return super().axis(c, e, form, factor)
        full = Model.axis(self.basis, c._replace(slab=0), e, form, factor)
        return restrict_controls(full, *self.boundaries[c.slab:c.slab+2])

    @staticmethod
    def slots(cells, kind):
        i, direction = kind
        if direction == "mother":
            return 1
        return cells[i].arity - int(i < len(cells) - 1)

    @classmethod
    def ports(cls, cells):
        yield 0, "mother"
        for i in range(len(cells)):
            kind = i, "children"
            if cls.slots(cells, kind):
                yield kind

    def chain_gain(self, row, cells, e, i, sign, *, difference=False):
        form = "difference" if difference else "old"
        axes = [self.axis(c, exp) for c, exp in zip(cells, e)]
        axes[i] = self.axis(cells[i], e[i], form, "leaves")
        self.add_product(row, cells, "mass", axes, sign)
        axes[i] = self.axis(cells[i], e[i], form)
        for j in (i-1, i+1):
            if 0 <= j < 3:
                internal = axes.copy()
                internal[j] = self.axis(cells[j], e[j], factor="p")
                self.add_product(row, cells, "mass", internal, sign)
        for kind in self.ports(cells):
            if kind[0] == i:
                self.add_product(row, cells, kind, axes, sign)

    def chain_balance(self, cells):
        labels = tuple(c.label for c in cells)
        for e in self.tests(labels):
            row = self.row(labels, e)
            self.add_observation(row, cells, e, discount=True)
            for i, c in enumerate(cells):
                self.chain_gain(row, cells, e, i, 1, difference=not c.degree)
                self.add_observation(row, cells, e, factor=(i, "birth"))
        for i, c in enumerate(cells):
            if c.degree:
                degree = c.degree + 1 if c.degree + 1 < self.threshold else 0
                after = labels[:i] + ((degree, c.arity),) + labels[i+1:]
                if after in self.chain_labels:
                    for e in self.tests(after):
                        self.chain_gain(self.row(after, e), cells, e, i, -1)
            after = labels[:i] + ((c.degree, c.arity + 1),) + labels[i+1:]
            if after in self.chain_labels:
                for e in self.tests(after):
                    self.add_observation(self.row(after, e), cells, e, -1,
                                         factor=(i, "birth"))

    def chain_entrances(self):
        for pair in self.edges:
            apex, middle = pair
            after = apex.label, (middle.degree, middle.arity + 1), self.virgin.label
            if after in self.chain_labels:
                for e in self.tests(after):
                    self.add_observation(self.row(after, e), pair, e[:2], -1,
                                         factor=(1, "birth"))

    def connect(self, pair, *, prepend, transport):
        """Disintegrate one edge over its outer mother or each child of its tip."""
        selector = (0, "mother") if prepend else (1, "children")
        count = self.slots(pair, selector)
        choices = []
        omitted_reciprocal = Q(1, self.threshold)
        for donor in self.cells:
            if prepend and not donor.arity:
                continue
            chain = (donor,) + pair if prepend else pair + (donor,)
            if (chain, "mass") in self.blocks:
                choices.append((donor, chain))
            elif donor.degree:
                omitted_reciprocal = max(omitted_reciprocal, Q(1, donor.degree))
        if not choices:
            return
        for idx, cols in self.blocks[pair, "mass"].items():
            port = self.blocks[pair, selector][idx][0]
            missing = defaultdict(Q, {i: Q(count) for i in cols})
            known = defaultdict(Q)
            # given[j] integrates the entire old receive field j over covered extensions.
            given = {kind: defaultdict(Q) for kind in self.ports(pair)} if transport else {}
            for donor, chain in choices:
                weights = self.axis(donor, 0)
                reciprocals = self.axis(donor, 0, factor="p")
                for j, (weight, reciprocal) in enumerate(zip(weights, reciprocals)):
                    extended_idx = (j,) + idx if prepend else idx + (j,)
                    masses = self.blocks[chain, "mass"][extended_idx]
                    for col in masses:
                        missing[col] -= weight
                        known[col] += reciprocal
                    for kind, row in given.items():
                        role, direction = kind
                        new_kind = role + int(prepend), direction
                        if kind == selector:
                            for col in masses:
                                row[col] += reciprocal
                        if (chain, new_kind) in self.blocks:
                            row[self.blocks[chain, new_kind][extended_idx][0]] += weight
            # The selector itself sees the true reciprocal of the identified donor.
            lower = dict(known)
            lower[port] = -1
            self.add_ineq(lower)
            upper = defaultdict(Q, {port: 1})
            for col, v in known.items():
                upper[col] -= v
            for col, v in missing.items():
                upper[col] -= omitted_reciprocal * v
            self.add_ineq(upper)
            for kind, seen in given.items():
                base_port = self.blocks[pair, kind][idx][0]
                seen[base_port] -= count
                self.add_ineq(seen)
                bound = Q(self.slots(pair, kind), 2)
                self.add_ineq({i: -seen.get(i, 0) - bound * missing.get(i, 0)
                               for i in seen.keys() | missing.keys()})

    def chain_objective(self, *, high_role, high_slabs, overlap="all", reciprocal=False):
        """Incidence or true reciprocal of a specified high role, per unit of eta."""
        if overlap not in ("all", "virgin-tip", "bare-tip"):
            raise ValueError(overlap)
        cost = defaultdict(Q)
        for chain in self.chains:
            high = chain[high_role]
            if high.degree or high.slab not in high_slabs:
                continue
            if overlap != "all" and chain[2].label != (2, int(overlap == "bare-tip")):
                continue
            factor = (high_role, "p") if reciprocal else None
            self.add_observation(cost, chain, (0, 0, 0), factor=factor)
        return cost
