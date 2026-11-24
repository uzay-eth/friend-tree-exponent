"""Price the conditioned newborn and the existing edge floors on one stationary law."""

import argparse
from contextlib import closing
from fractions import Fraction as Q
import json
from pathlib import Path
import re
import sys


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "stationary"))
import joint


ROOT_LOG = "ef19be540999333210d6c503cbd9d54bbfd07a504c8654b59abe47fefc7b416e"
PAYING_FY = "56b31c406cfe1d794a85f6dea700211eaedbf297e5e1cd00470896e0bd28a824"
COST_SCALE = 10**6
DUAL_GRID = 10**12


def conditioned_roots():
    records = joint.read_pinned(HERE / "roots.log", ROOT_LOG).decode().splitlines()
    pattern = (r"(lower|upper) lower states 72322 corners 750 roots "
               r"\(Fraction\((\d+), (\d+)\), Fraction\((\d+), (\d+)\)\) ")
    roots = {}
    for line in records:
        match = re.match(pattern, line)
        assert match, line
        side, *integers = match.groups()
        a, b, c, d = map(int, integers)
        assert side not in roots
        roots[side] = Q(a, b), Q(c, d)
    assert set(roots) == {"lower", "upper"}
    return roots


def edge_floors(scale):
    result = joint.edge_prices()
    path = HERE.parent / "edges/fy/extended/cut36-floor.log"
    raw = joint.read_pinned(path, PAYING_FY).decode().splitlines()
    [extra] = [Q(line.split()[1]) for line in raw if line.startswith("TOTAL ")]
    assert extra > 0 and raw[0].startswith("START lower bound ")
    result["lower"].append(scale * extra)
    return result


def column_costs(lp, cells, gold, weight_cap, roots, edge_credit):
    # eta = pi/mu. A constant per physical core is integrated against every
    # occupation column; the living mother's birth rate multiplies only V.
    second, older = roots
    fixed = sum(edge_credit) - gold["err_per_core"]
    numerators, denominators = [], []
    for kind, entries, meta in lp.cols:
        value, upper = 0, 0
        if kind == "occ":
            d, k, birth, *_ = meta
            assert dict(entries).get(lp.rows["birth"], 0) == birth
            source = birth * (second if d == 2 else older)
            _, upper, cap = cells[d, k]
            value = fixed + weight_cap - cap + source
        numerators.append(Q(value))
        denominators.append(Q(upper))
    return numerators, denominators


def ratio_proposal(lp, numerators, denominators, scale):
    """Perspective LP. Its floating solution chooses a trial for exact verification."""
    import numpy as np
    from scipy.optimize import linprog
    from scipy.sparse import bmat, csc_matrix

    matrix, rhs = lp.matrix()
    density = [float(w / scale) for w in denominators]
    equalities = bmat([[matrix, csc_matrix(-np.array(rhs, dtype=float)[:, None])],
                       [csc_matrix([density]), csc_matrix((1, 1))]], format="csc")
    objective = [float(x / COST_SCALE) for x in numerators] + [0.]
    target = [0] * len(rhs) + [1]
    answer = linprog(objective, A_eq=equalities, b_eq=target, method="highs-ipm")
    if not answer.success:
        answer = linprog(objective, A_eq=equalities, b_eq=target, method="highs-ds")
    if not answer.success:
        raise ValueError(answer.message)
    return Q.from_float(answer.fun) * COST_SCALE / scale


def suggest_dual(lp, cost):
    from scipy.optimize import linprog

    matrix, rhs = lp.matrix()
    c = [float(x / COST_SCALE) for x in cost]
    options = {"dual_feasibility_tolerance": 1e-9, "primal_feasibility_tolerance": 1e-9}
    answer = linprog(c, A_eq=matrix, b_eq=rhs, method="highs-ipm", options=options)
    if not answer.success:
        answer = linprog(c, A_eq=matrix, b_eq=rhs, method="highs-ds")
    if not answer.success:
        raise ValueError(answer.message)
    return [round(z * DUAL_GRID) for z in answer.eqlin.marginals]


def signed_margin(lp, costs, integer_dual, lower_prior):
    if len(integer_dual) != len(lp.rows) or any(type(y) is not int for y in integer_dual):
        raise ValueError("invalid exact dual")
    y = [Q(z * COST_SCALE, DUAL_GRID) for z in integer_dual]
    defects = {}
    for cost, (kind, terms, _) in zip(costs, lp.cols, strict=True):
        error = sum((y[i] * a for i, a in terms), Q()) - cost
        defects[kind] = max(defects.get(kind, Q()), error)
    assert set(defects) == {"occ", "stubs", "tail-up", "palm-slack", "edge"}
    assert lp.rhs == {"birth": 1}
    return y[lp.rows["birth"]] - sum(defects.values()) / lower_prior


def evaluate(side, lp, numerators, denominators, gold, scale, lower_prior, solve):
    path = HERE / f"{side}.json"
    if solve:
        optimum = ratio_proposal(lp, numerators, denominators, scale)
        # Take a small inward rounding allowance for the LP and rational repair.
        trial = Q((optimum * 10**13).__floor__() - 2, 10**13)
    else:
        saved = json.loads(path.read_bytes())
        assert saved["side"] == side and saved["scale"] == scale
        assert saved["cost_scale"] == COST_SCALE and saved["dual_grid"] == DUAL_GRID
        trial = Q(saved["trial"])
    assert trial > 0
    cost = [n - trial * d for n, d in zip(numerators, denominators, strict=True)]
    y = suggest_dual(lp, cost) if solve else saved["integer_dual"]
    margin = signed_margin(lp, cost, y, lower_prior)
    if margin <= 0:
        raise ValueError(f"uncovered trial {trial}: exact signed margin {margin}")
    direction = 1 if side == "lower" else -1
    endpoint = Q(gold["nominal"]) + direction * trial
    if solve:
        record = dict(side=side, scale=scale, cost_scale=COST_SCALE, dual_grid=DUAL_GRID,
                      trial=str(trial), integer_dual=y)
        path.write_text(json.dumps(record, separators=(",", ":")) + "\n")
    print(side, "inward displacement from nominal", trial, "exact signed eta margin", margin,
          "endpoint", endpoint, "outward", joint.host.outward(endpoint, side), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--solve", action="store_true", help="generate exact dual receipts")
    parser.add_argument("--side", choices=("lower", "upper"), action="append")
    args = parser.parse_args()
    assert __debug__
    joint.read_pinned(Path(joint.__file__),
                      "a9f3677f0b00780ff284a4836042af3355134b2904a14483c7a2f16ce8f55bb8")
    gold = joint.authenticated_receipt(args.work / "gold.json")
    low, high = (Q(gold[s]["certified_exact"]) for s in ("lower", "upper"))
    audit = joint.pinned_module("quintet_source_kernel", joint.R2 / "checks/stationary-p4-sol.py",
                                "b50cc2d2e51a6e7ab754e2d851abdd690cf20922")
    assert high <= audit.HI
    model = joint.host.birth_model(low, high, 7)
    lp = model.build(24, sided=True, kernel=audit.class_kernel())
    assert (len(lp.rows), len(lp.cols)) == (613, 5554)
    roots = conditioned_roots()
    for side in args.side or ("lower", "upper"):
        data = gold[side]
        with closing(joint.Comparison(args.work / f"quintet_{side}.bin", side, data)) as weights:
            cells = joint.host.census_cells(weights, 24, 7, data)
            extras = edge_floors(weights.scale)[side]
            numerators, denominators = column_costs(lp, cells, data, weights.higher_cap, roots[side], extras)
            print(side, "birth class lifetime ticks", *roots[side], "edge floors", *extras, flush=True)
            evaluate(side, lp, numerators, denominators, data, weights.scale, low, args.solve)
            weights.check_identity()
    print("PASS exact signed checks on every original degree-24 column")


if __name__ == "__main__":
    main()
