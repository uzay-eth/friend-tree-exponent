"""Price a physical high-degree hV band and its attached chain incidences."""

import argparse
from collections import defaultdict
from fractions import Fraction as Q
from pathlib import Path
from time import process_time

from chains import ChainModel
from solve import optimize


def high_charge(model, payloads, role, mode, theta, minimum=18):
    from higher.hv import H, Reader, Role, contrast, moment_controls, opener

    if model.m < 3 or model.threshold > H or minimum < H:
        raise ValueError("the original hV tail requires order >=3, threshold <=18 and D>=18")
    readers = [Reader(opener(p, name, authenticate=True))
               for p, name in zip(payloads, ("lower", "upper"), strict=True)]
    cost = defaultdict(Q)
    paid = 0
    for chain in model.chains:
        selected = chain[role]
        if (selected.degree or model.boundaries[selected.slab+1] > Q(model.threshold, minimum)
                or any(not c.degree for i, c in enumerate(chain) if i != role)):
            continue
        roles = tuple(Role(c.degree or None, c.arity) for c in chain)
        if mode == "lower":
            row = readers[0].row(roles)
        elif mode == "upper":
            row = readers[1].row(roles)
        else:
            row = contrast(*readers, roles, theta if mode == "difference" else -theta)
        sign = -1 if mode == "upper" else 1
        bounds = [(Q(0), Q(1))] * 3
        bounds[role] = tuple(Q(H, model.threshold) * b
                             for b in model.boundaries[selected.slab:selected.slab+2])
        base, fields = moment_controls(row, model.m, bounds)
        for kind, coefficients in [("mass", base)] + [(k, net) for k, (_, net) in fields.items()]:
            for idx, v in coefficients.items():
                for column in model.blocks[chain, kind][idx]:
                    cost[column] += sign * v
        paid += 1
    if not paid:
        raise ValueError("no high-degree chain cell meets the selected support")
    return paid, cost


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--threshold", type=int, default=4)
    p.add_argument("--arity", type=int, default=2)
    p.add_argument("--chain-arity", type=int, nargs=3, default=(2, 1, 1))
    p.add_argument("--inner", type=int, default=3)
    p.add_argument("--role", type=int, choices=range(3), default=0)
    p.add_argument("--minimum", type=int, default=18)
    p.add_argument("--order", type=int, default=3)
    p.add_argument("--mode", nargs="+", choices=["mass", "reciprocal", "difference", "sum", "lower", "upper"],
                   default=["mass", "difference"])
    p.add_argument("--theta", type=Q, default=Q(1))
    p.add_argument("--payloads", nargs=2, type=Path)
    p.add_argument("--no-transport", action="store_true")
    p.add_argument("--no-chain-projection", action="store_true")
    p.add_argument("--no-chain-balances", action="store_true")
    p.add_argument("--method", choices=["highs", "highs-ipm", "highs-ds"], default="highs-ipm")
    args = p.parse_args()
    if not (args.inner < args.threshold <= args.minimum):
        p.error("require inner < threshold <= minimum")
    cuts = (Q(args.threshold, args.minimum),) if args.threshold < args.minimum else ()
    start = process_time()
    m = ChainModel(args.threshold, args.arity, args.order, breaks=cuts,
                   inner=args.inner, wide_roles=(args.role,), chain_arity=args.chain_arity,
                   project=not args.no_chain_projection, port_transport=not args.no_transport,
                   balances=not args.no_chain_balances)
    print("threshold", args.threshold, "reciprocal cuts", m.boundaries, "axis", args.role,
          "arity caps", tuple(args.chain_arity), "degree cap at other roles", args.inner, flush=True)
    print("columns", len(m.variables), "equalities", len(m.equalities), "upper rows", len(m.upper_rows),
          "chain cells", len(m.chains), "CPU build", round(process_time()-start, 3), flush=True)
    for mode in args.mode:
        if mode in ("mass", "reciprocal"):
            c = m.chain_objective(high_role=args.role, high_slabs=(0,), reciprocal=(mode == "reciprocal"))
            scale = 10
        else:
            if not args.payloads:
                p.error("the physical charge requires both original D18 payloads")
            paid, c = high_charge(m, args.payloads, args.role, mode, args.theta, args.minimum)
            print("authenticated two D18 witnesses, paid cells", paid, flush=True)
            scale = 10**7
        bounds = []
        for maximize in (False, True):
            v, err = optimize(m, c, maximize, scale=scale, method=args.method)
            physical = (max if maximize else min)(v*s for s in (m.low, m.high))
            bounds.append(physical)
            print(mode, "upper" if maximize else "lower", "exact per core", physical,
                  "decimal", float(physical), "raw repair", float(err), flush=True)
        print(mode, "width", float(bounds[1]-bounds[0]), "CPU", round(process_time()-start, 3), flush=True)


if __name__ == "__main__":
    main()
