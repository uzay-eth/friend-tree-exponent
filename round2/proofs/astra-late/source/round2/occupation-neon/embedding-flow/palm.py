"""Couple a physical birth-law LP to stationary rooted incidence balances."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

from solve import IncidenceModel, Pattern, optimize
from identification.exact import install_cuts

HERE = Path(__file__).resolve().parent
MEMORY = HERE.parent / "pressure-memory"


def load_module(name, path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def coarse_cell(degree, arity, limit):
    if degree <= 12:
        band = degree
    elif degree <= 18:
        band = 13
    elif degree <= limit:
        band = 19
    else:
        band = limit+1
    return band, min(arity, 3)


def install_sided_law(model, limit):
    """Match true finite-census masses, both receive fluxes, and common edges."""
    if limit < 18 or len(model.caps) < 2:
        raise ValueError("the sided representation requires a degree limit of at least 18 and two orders")
    source = load_module("sided_birth", MEMORY / "birth-law" / "replay.py")
    audit = source.load("birth_minorant", source.R2 / "checks" / "stationary-p4-sol.py")
    sided = source.build(limit, True, audit.class_kernel())
    sides = [dict() for _ in sided.rows]
    occupied = defaultdict(list)
    broad_edges = defaultdict(list)
    total = {model.total: -1}
    for j, (kind, values, meta) in enumerate(sided.cols):
        i = model.variable(("sided", j), ("sided", kind), 1/source.LO)
        for row, v in values:
            sides[row][i] = v
        if kind == "occ":
            d,k,b,p,u,v = meta
            total[i] = 1
            if d <= min(limit,model.cap):
                rate = b - dict(values).get(sided.rows["mass"], 0)
                occupied[Pattern(d,k)].append((i, rate, u, v))
        elif kind == "edge":
            a,b,_,_ = meta
            broad_edges[a,b].append(i)
    for name, r in sided.rows.items():
        model.record(model.equalities, sides[r], sided.rhs.get(name, 0))
    model.record(model.equalities, total, 0)

    for tree, cols in occupied.items():
        count = defaultdict(int, {i: 1 for i,_,_,_ in cols})
        model.mass(count, tree, -1)
        model.record(model.equalities, count, 0)
        discounted = {i: s for i,s,_,_ in cols}
        for i,s in zip(model.masses[tree], (model.low,model.high)):
            discounted[i] = -s
        model.record(model.equalities, discounted, 0)
        for parent in (True,False):
            port = model.port(tree, (), parent)
            if port is not None:
                row = {i: v if parent else u for i,_,u,v in cols}
                row[port] = -1
                model.record(model.equalities, row, 0)

    # Each nonfounder parent-child edge is counted once from either endpoint.
    # A partially retained degree bin gives a lower bound on its coarse mass.
    fine = defaultdict(list)
    for tree in model.levels[2]:
        child, = tree.children
        pair = (coarse_cell(tree.degree,tree.arity,limit),
                coarse_cell(child.degree,child.arity,limit))
        fine[pair].append(tree)

    def fully_read(band):
        if band == limit+1:
            return False
        maximum = band if band <= 12 else 18 if band == 13 else limit
        return maximum <= model.caps[1]

    for pair, trees in fine.items():
        row = defaultdict(int, {i:-1 for i in broad_edges[pair]})
        if not row:
            raise ValueError("exact embedding has no physical sided edge class")
        for tree in trees:
            model.mass(row,tree,1)
        rows = model.equalities if all(fully_read(c[0]) for c in pair) else model.upper_rows
        model.record(rows,row,0)
    return len(sided.cols),len(occupied),len(fine)


def read_prices():
    directory = MEMORY / "path-fast"
    pool = load_module("physical_path_prices", directory / "pool.py")
    header, entries = pool.read_entries(
        [directory/"wide-head-entries.json",directory/"wide-tail-entries.json"],
        refine=directory/"wide-mother12-entries.json",
    )
    return header, entries


def payoff(model, entries, column, unit):
    """The constant uses the total birth mass one, including all omitted mothers."""
    constant = Q(entries["overflow"][column], unit)
    cost = defaultdict(int)
    for tree in model.levels[1]:
        v = entries.get(str(tree.degree))
        if v is None:
            continue
        birth = Q(tree.degree-tree.arity-1, tree.degree)
        model.mass(cost, tree, birth*(Q(v[column],unit)-constant))
    return constant, dict(cost)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--caps", type=int, nargs="+", default=[18,5,3])
    parser.add_argument("--order", type=int, default=2)
    parser.add_argument("--share-through", type=int, default=2)
    parser.add_argument("--limit", type=int, default=18)
    parser.add_argument("--sided", action="store_true")
    parser.add_argument("--mass-cuts", action="store_true")
    parser.add_argument("--upper-support", action="store_true",
                        help="also optimize the upper conditional endpoint in each price vector")
    parser.add_argument("--method", choices=["highs-ipm","highs-ds","highs"],default="highs-ipm")
    args = parser.parse_args()
    header, entries = read_prices()
    if args.caps[0] < header["mother_cap"]:
        parser.error("the root degree cap must cover the whole finite conditional partition")
    model = IncidenceModel(args.caps[0],args.order,args.share_through,caps=args.caps)
    if args.mass_cuts:
        install_cuts(model)
    if args.sided:
        print("sided columns, coupled root censuses, coupled edge bins",
              *install_sided_law(model,args.limit),flush=True)
    print("patterns", [len(x) for x in model.levels.values()],
          "columns",len(model.variables),"equalities",len(model.equalities),
          "upper rows",len(model.upper_rows),flush=True)
    for side, j in (("lower",0),("upper",2)):
        for lane in range(1+args.upper_support):
            constant, cost = payoff(model,entries,j+lane,header["grid"])
            value, repair = optimize(model,cost,bool(lane),method=args.method)
            prices = [s*(constant+value) for s in (model.low,model.high)]
            bound = max(prices) if lane else min(prices)
            print(side,"ceiling" if lane else "floor","per core",bound,float(bound),
                  "unnormalized dual repair",float(repair),flush=True)


if __name__ == "__main__":
    main()
