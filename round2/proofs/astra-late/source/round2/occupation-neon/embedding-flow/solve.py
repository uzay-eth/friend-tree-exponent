"""Stationary balances and shared exterior-donor marginals of rooted core patterns."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import cache
from itertools import combinations_with_replacement, permutations, product
from math import factorial, prod
from pathlib import Path
import sys
from typing import NamedTuple

Q = Fraction
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "pressure-memory"))
from path_band import original
from extract import SPEC

LO = Q("0.565678493279215")
HI = Q("0.565770152011715")


class Pattern(NamedTuple):
    degree: int
    arity: int
    children: tuple = ()


def partitions(total, first=1):
    if total == 0:
        yield ()
    for n in range(first, total + 1):
        for rest in partitions(total - n, n):
            yield (n,) + rest


@cache
def largest_degree(tree):
    return max(tree.degree, max((largest_degree(c) for c in tree.children), default=0))


def generate_patterns(cap, max_size, caps=None):
    limits = (cap,)*max_size if caps is None else tuple(caps)
    if len(limits) != max_size or limits[0] != cap or min(limits) < 2:
        raise ValueError("one degree limit, starting with the root limit, is required at every order")
    if any(a < b for a, b in zip(limits, limits[1:])):
        raise ValueError("degree limits must be nonincreasing with pattern size")
    levels = {}
    for n in range(1, max_size + 1):
        limit = limits[n-1]
        eligible = {size: [t for t in level if largest_degree(t) <= limit]
                    for size, level in levels.items()}
        current = []
        for partition in partitions(n - 1):
            multiplicities = Counter(partition)
            groups = [combinations_with_replacement(eligible[size], multiplicity)
                      for size, multiplicity in multiplicities.items()]
            for forest in product(*groups):
                children = tuple(sorted(child for group in forest for child in group))
                for d in range(max(2, len(children) + 1), limit + 1):
                    for k in range(len(children), d):
                        current.append(Pattern(d, k, children))
        levels[n] = tuple(sorted(current))
    return levels


@cache
def roles(tree):
    return ((),) + tuple((i,) + tail for i, child in enumerate(tree.children)
                         for tail in roles(child))


def vertex(tree, role):
    for i in role:
        tree = tree.children[i]
    return tree


def first_equivalent_role(tree, role):
    if not role:
        return ()
    child = tree.children[role[0]]
    return (tree.children.index(child),) + first_equivalent_role(child, role[1:])


def replace(tree, role, changed):
    """Return the canonical pattern and the new address of the changed vertex."""
    if not role:
        return changed, ()
    i = role[0]
    child, address = replace(tree.children[i], role[1:], changed)
    children = tuple(sorted(tree.children[:i] + (child,) + tree.children[i+1:]))
    return Pattern(tree.degree, tree.arity, children), (children.index(child),) + address


@cache
def topology(tree):
    return tuple(sorted(topology(child) for child in tree.children))


@cache
def automorphisms(tree):
    return (prod(factorial(m) for m in Counter(tree.children).values())
            * prod(automorphisms(child) for child in tree.children))


def lookup_family(tree):
    addresses = roles(tree)
    for family, parents in SPEC.items():
        if len(parents) != len(addresses):
            continue
        for order in permutations(addresses[1:]):
            data = ((),) + order
            if all(data[i][:-1] == data[parents[i]] for i in range(1, len(data))):
                return family, data
    raise ValueError("no D18 family for the given tree")


def add(row, index, coefficient):
    if coefficient:
        row[index] += coefficient


class IncidenceModel:
    def __init__(self, degree_cap, order, shared_order, low=LO, high=HI, *, caps=None):
        if not 0 <= shared_order <= order or low <= 0 or high < low:
            raise ValueError("invalid lift or discount")
        self.cap = degree_cap
        self.order = order
        self.shared_order = shared_order
        self.low, self.high = low, high
        depth = max(order, shared_order+1)
        self.caps = (degree_cap,)*depth if caps is None else tuple(caps)
        self.levels = generate_patterns(degree_cap, depth, self.caps)
        self.equalities = []
        self.upper_rows = []
        self.variables = []
        self.group_bounds = {}
        self.masses = {}
        self.flux = {}
        self.balance = {}
        self.projections = []

        self.total = self.variable("total", "total", 1/low)
        for size, patterns in self.levels.items():
            for tree in patterns:
                # Each concrete mass counts all embeddings of the labelled roles.
                group = ("mass", topology(tree))
                bound = (self.caps[size-1]-1)**(size-1) / low
                indices = [self.variable((tree, "low"), group, bound)]
                if size <= order and low < high:
                    indices.append(self.variable((tree, "high"), group, bound))
                self.masses[tree] = tuple(indices)
        self.virgin = Pattern(2, 0)

        for size in range(1, order+1):
            for tree in self.levels[size]:
                self.make_port(tree, (), True, 1)
                for role in roles(tree):
                    x = vertex(tree, role)
                    remaining = x.arity - len(x.children)
                    if remaining:
                        self.make_port(tree, role, False, remaining)
        self.global_rows()
        for size in range(1, order+1):
            for tree in self.levels[size]:
                self.stationary_row(tree)

    def variable(self, name, group, bound):
        i = len(self.variables)
        self.variables.append((name, group))
        self.group_bounds[group] = bound
        return i

    def record(self, rows, row, rhs):
        rows.append(({i: v for i, v in row.items() if v}, rhs))
        return len(rows)-1

    def mass(self, row, tree, coefficient):
        for i in self.masses[tree]:
            add(row, i, coefficient)

    def extension_cap(self, tree):
        """Largest resolved donor degree for this existing census, or zero."""
        size = len(roles(tree))
        if size >= len(self.caps) or largest_degree(tree) > self.caps[size]:
            return 0
        return self.caps[size]

    def make_port(self, tree, role, parent, count):
        key = tree, first_equivalent_role(tree, role), parent
        if key in self.flux:
            return
        size = len(roles(tree))
        group = ("port", topology(tree), parent)
        # Repeated and nonequivalent child slots are bounded simultaneously.
        limit = self.caps[size-1]
        cap = (limit-1)**(size-1) / self.low * size * limit / 2
        index = self.variable(key, group, cap)
        self.flux[key] = index
        row = defaultdict(int, {index: 1})
        donor_limit = self.extension_cap(tree) if size <= self.shared_order else 0
        max_unseen = Q(1, donor_limit+1 if donor_limit else 2)
        self.mass(row, tree, -count * max_unseen)
        if donor_limit:
            minimum = defaultdict(int, {index: -1})
            for donor in self.levels[1]:
                if donor.degree > donor_limit:
                    continue
                if parent:
                    if donor.arity == 0:
                        continue
                    extended_tree = Pattern(donor.degree, donor.arity, (tree,))
                else:
                    x = vertex(tree, role)
                    children = tuple(sorted(x.children + (donor,)))
                    new = Pattern(x.degree, x.arity, children)
                    extended_tree, _ = replace(tree, role, new)
                reciprocal = Q(1, donor.degree)
                self.mass(minimum, extended_tree, reciprocal)
                self.mass(row, extended_tree, max_unseen-reciprocal)
            self.projections.append(self.record(self.upper_rows, minimum, 0))
        self.projections.append(self.record(self.upper_rows, row, 0))

    def port(self, tree, role, parent):
        return self.flux.get((tree, first_equivalent_role(tree, role), parent))

    def global_rows(self):
        self.record(self.upper_rows, {self.total: -1}, -1/self.high)
        self.record(self.upper_rows, {self.total: 1}, 1/self.low)
        rate_tail = prod((Q(2*d-1, 2)/(Q(2*d-1, 2)+self.low)
                          for d in range(2, self.cap+1)), start=Q(1))
        mass = defaultdict(int, {self.total: -1})
        arity = defaultdict(int, {self.total: -1})
        yule = defaultdict(int, {self.total: 1-rate_tail})
        birth = defaultdict(int)
        omitted_birth = defaultdict(int, {self.total: -1})
        omitted_reciprocal = defaultdict(int, {self.total: 1-Q(2,self.cap+1)})
        for tree in self.levels[1]:
            d, k, _ = tree
            b = Q(d-k-1, d)
            self.mass(mass, tree, 1)
            self.mass(arity, tree, k)
            self.mass(yule, tree, -1)
            self.mass(birth, tree, b)
            self.mass(omitted_birth, tree, 1-b)
            self.mass(omitted_reciprocal, tree, Q(1+k,self.cap+1) - Q(1+k,d))
        for row, rhs in [(mass,0), (arity,0), (yule,0), (birth,1),
                         (omitted_birth,-1), (omitted_reciprocal,1)]:
            self.record(self.upper_rows, row, rhs)

    def gain(self, row, tree, role, sign):
        x = vertex(tree, role)
        fixed = x.degree-x.arity-1 + sum((Q(1,c.degree) for c in x.children), Q(0))
        if role:
            fixed += Q(1, vertex(tree, role[:-1]).degree)
        else:
            add(row, self.port(tree, role, True), sign)
        side = self.port(tree, role, False)
        if side is not None:
            add(row, side, sign)
        self.mass(row, tree, sign * fixed)

    def stationary_row(self, tree):
        row = defaultdict(int)
        for i, s in zip(self.masses[tree], (self.low, self.high)):
            add(row, i, s)
        for role in roles(tree):
            x = vertex(tree, role)
            d, k, ch = x
            self.gain(row, tree, role, 1)
            self.mass(row, tree, Q(d-k-1,d))
            if d > 2 and d-1 > k:
                prev, new_role = replace(tree, role, Pattern(d-1,k,ch))
                self.gain(row, prev, new_role, -1)
            if k > len(ch):
                prev, _ = replace(tree, role, Pattern(d,k-1,ch))
                self.mass(row, prev, -Q(d-k,d))
            if role and x == self.virgin:
                mother = vertex(tree, role[:-1])
                others = mother.children[:role[-1]] + mother.children[role[-1]+1:]
                previous = Pattern(mother.degree, mother.arity-1, others)
                source, _ = replace(tree, role[:-1], previous)
                self.mass(row, source, -Q(mother.degree-mother.arity,mother.degree))
        self.balance[tree] = self.record(self.equalities, row, int(tree == self.virgin))

    def objective(self, families, tables, theta, unpaid=False, *, price_cap=16):
        if not 2 <= price_cap <= 16:
            raise ValueError("inner D18 charge extraction requires degree at most 16")
        cost = defaultdict(int)
        for size in range(1, self.order+1):
            for tree in self.levels[size]:
                if largest_degree(tree) > price_cap:
                    continue
                family, addresses = lookup_family(tree)
                if family not in families:
                    continue
                xs = [vertex(tree, role) for role in addresses]
                labels = tuple((x.degree-x.arity-1, x.arity, 1) for x in xs)
                if unpaid:
                    n_side = sum(x.arity for x in xs)-(size-1)
                    size_min = sum(x.degree-x.arity for x in xs) + 2*n_side
                    if size_min <= 12 or (size_min <= 20 and size+n_side <= 4):
                        continue
                for table, factor in zip(tables, (Q(1), -theta)):
                    if not factor or not table.supported(family, labels):
                        continue
                    base, ports = table.row(family, labels)
                    scale = factor / (table.c.S*automorphisms(tree))
                    self.mass(cost, tree, base * scale)
                    for i, is_parent, _, slope in ports:
                        column = self.port(tree, addresses[i], is_parent)
                        if column is None:
                            raise ArithmeticError("charged port missing")
                        add(cost, column, slope*scale)
        return {i: val for i, val in cost.items() if val}


def optimize(model, cost, maximize, *, method="highs", scale=10**7):
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import coo_matrix

    def matrix(rows):
        rr, cc, values = [], [], []
        for j, (row, _) in enumerate(rows):
            for i, value in row.items():
                rr.append(j)
                cc.append(i)
                values.append(float(value))
        shape = len(rows), len(model.variables)
        return coo_matrix((values,(rr,cc)),shape=shape).tocsr(), [float(b) for _,b in rows]

    objective = {i: (-x if maximize else x) for i,x in cost.items()}
    coefficients = np.zeros(len(model.variables))
    for i, value in objective.items():
        coefficients[i] = float(value * scale)
    eq, b_eq = matrix(model.equalities)
    upper, b_upper = matrix(model.upper_rows)
    solution = linprog(coefficients, A_eq=eq, b_eq=b_eq, A_ub=upper, b_ub=b_upper,
                       bounds=(0,None), method=method,
                       options={"threads": 1, "dual_feasibility_tolerance": 1e-9})
    if not solution.success:
        raise ArithmeticError(solution.message)

    reduced = defaultdict(int, objective)
    center = Q(0)
    denominator = 10**13
    for rows, multipliers, signed in [(model.equalities,solution.eqlin.marginals,True),
                                     (model.upper_rows,solution.ineqlin.marginals,False)]:
        for (row,rhs), weight in zip(rows, multipliers):
            y = Q(round(weight*denominator), denominator*scale)
            if not signed:
                y = min(0,y)
            if y:
                center += rhs*y
                for i,value in row.items():
                    reduced[i] -= value*y
    worst = defaultdict(Q)
    for i, residual in reduced.items():
        group = model.variables[i][1]
        worst[group] = max(worst[group], -residual)
    repair = sum((model.group_bounds[group]*error for group,error in worst.items()), Q(0))
    value = center-repair
    return (-value if maximize else value), repair


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--lower", type=Path, required=True)
    parser.add_argument("--upper", type=Path, required=True)
    parser.add_argument("--cap", type=int, default=3)
    parser.add_argument("--caps", type=int, nargs="+",
                        help="nonincreasing degree limits, one per retained pattern size; overrides --cap")
    parser.add_argument("--price-cap", type=int, default=16,
                        help="largest motif degree in the charged objective; balance rows may extend higher")
    parser.add_argument("--order", type=int, default=4)
    parser.add_argument("--share-through", type=int, default=3)
    parser.add_argument("--families", nargs="+", choices=tuple(SPEC), default=["hV","hC","P4","Yt","Yb","S3"])
    parser.add_argument("--theta", type=Q, default=Q(1))
    parser.add_argument("--outside-size-reserve", "--unpaid", dest="unpaid", action="store_true",
                        help="exclude the historical size-12 and size-20 supports")
    parser.add_argument("--method", choices=("highs","highs-ds","highs-ipm"), default="highs")
    parser.add_argument("--mass-cuts", action="store_true",
                        help="include child- and mother-count projections at every retained degree threshold")
    args = parser.parse_args()
    if args.caps:
        args.cap = args.caps[0]
    if not 2 <= args.cap or not 2 <= args.price_cap <= 16 or not 1 <= args.order <= 4 or not 0 <= args.share_through <= args.order:
        parser.error("invalid cap or order")
    if args.unpaid and any(f in ("f","g") for f in args.families):
        parser.error("the complete-fringe size gate applies only to original higher families")
    tables = original(args.lower, "lower"), original(args.upper, "upper")
    model = IncidenceModel(args.cap, args.order, args.share_through, caps=args.caps)
    if args.mass_cuts:
        from identification.exact import install_cuts
        install_cuts(model)
    print("authenticated both originals; concrete patterns by size",
          {n:len(states) for n,states in model.levels.items()}, flush=True)
    print("columns, balance rows, inequality rows",
          len(model.variables),len(model.equalities),len(model.upper_rows),flush=True)
    objective = model.objective(set(args.families),tables,args.theta,args.unpaid,price_cap=args.price_cap)
    lower, rlo = optimize(model,objective,False,method=args.method)
    upper, rhi = optimize(model,objective,True,method=args.method)
    candidates = [v*s for v in (lower,upper) for s in (LO,HI)]
    print("unnormalized contrast",lower,upper)
    print("exact dual repairs",rlo,rhi)
    print("stationary contrast per core",min(candidates),max(candidates))
    print("decimal for readability",float(min(candidates)),float(max(candidates)),flush=True)


if __name__ == "__main__":
    main()
