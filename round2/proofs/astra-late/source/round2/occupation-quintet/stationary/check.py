"""Price the five-vertex host allocation on the physical stationary census."""

import argparse
from contextlib import closing
from fractions import Fraction as F
import hashlib
import importlib.util
import json
from math import ceil, floor
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from weights import Comparison, R2, authenticated_receipt, physical


def birth_model(lower, upper, saturation):
    path = R2 / "occupation-neon/pressure-memory/birth-law/replay.py"
    raw = path.read_bytes()
    blob = hashlib.sha1(f"blob {len(raw)}\0".encode() + raw).hexdigest()
    if blob != "7dd40ac41e5d95ae41843da51de436fd45b781f1":
        raise ValueError("birth and edge balance implementation changed")
    spec = importlib.util.spec_from_file_location("quintet_birth_balance", path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    module.LO, module.HI, module.SAT = lower, upper, saturation
    return module


def census_cells(weights, limit, saturation, receipt):
    """Envelopes in exact cells and in the joint degree/offspring overflow."""
    if not 18 <= limit < weights.h["DF"] - 1 or not 1 <= saturation < limit:
        raise ValueError("choose 18 <= limit < DF-1 and 1 <= saturation < limit")
    cells = {}
    global_low, global_high = None, None
    for d in range(2, weights.h["DF"]):
        lower, upper = weights.finite(d)
        row_low, row_high = int(lower.min()), int(upper.max())
        global_low = row_low if global_low is None else min(row_low, global_low)
        global_high = row_high if global_high is None else max(row_high, global_high)
        if d <= limit:
            for i, (lo, hi) in enumerate(zip(lower, upper)):
                cells[d, i] = int(lo), int(hi), weights.host_cap.get((d, i+1), 0)
            continue
        for k in range(saturation+1):
            if k == saturation:
                lo = int(lower[k:].min())
                hi = int(upper[k:].max())
                cap = max((weights.host_cap.get((d, c), 0)
                           for c in range(k+1, min(d, weights.h["K3m"])+1)), default=0)
            else:
                lo, hi = int(lower[k]), int(upper[k])
                cap = weights.host_cap.get((d, k+1), 0)
            cell = limit+1, k
            old = cells.get(cell, (lo, hi, cap))
            cells[cell] = min(old[0], lo), max(old[1], hi), max(old[2], cap)
    for k in range(saturation+1):
        lo, hi = weights.far(k+1, saturated=k == saturation)
        old = cells[limit+1, k]
        cells[limit+1, k] = min(old[0], lo), max(old[1], hi), old[2]
    expected = receipt["verify"]["comp_fin"]
    assert (global_low, global_high) == (expected["min_W"], expected["max_W"])
    assert all(lo >= weights.global_low and hi <= weights.global_high
               and 0 <= cap <= weights.higher_cap for lo, hi, cap in cells.values())
    return cells


def physical_comparison(weights):
    """Check individual allocated values on labelled trees in two rootings."""
    checked = 0
    for name, tree in physical.fixtures(weights.h):
        core, marks = physical.physical_core(tree)
        values = physical.subsets(weights.snapshot, core, marks)
        if len(core) == 1:
            assert len(values) == 1
            continue
        for root in sorted({min(core), max(core)}):
            parent, stack = {root: None}, [root]
            while stack:
                u = stack.pop()
                for w in core[u]:
                    if w not in parent:
                        parent[w] = u
                        stack.append(w)
            owned = {u: F(0) for u in core}
            for (family, vertices), (score, inside) in values.items():
                if family == "g":
                    a,b = vertices
                    hosts = (a if parent[a] == b else b,)
                elif family == "p4":
                    hosts = inside
                elif family == "p5":
                    hosts = (inside[1],)
                else:
                    hosts = (inside[0],)
                for u in hosts:
                    owned[u] += F(score, len(hosts))
            assert sum(owned.values()) == sum(score for score,_ in values.values())
            assert sum(c-2 for _, c in marks.values()) == -2
            for u, (d,c) in marks.items():
                if d < weights.h["DF"]:
                    lo, hi = weights.finite(d)
                    low, high = int(lo[c-1]), int(hi[c-1])
                else:
                    low, high = weights.far(c)
                assert low <= owned[u] + weights.alpha*(c-2), (name, u, "lower")
                assert owned[u] <= high, (name, u, "upper")
                checked += 1
    return checked


def outward(value, direction, digits=15):
    ticks = (floor if direction == "lower" else ceil)(value * 10**digits)
    sign = "-" if ticks < 0 else ""
    units, part = divmod(abs(ticks), 10**digits)
    return f"{sign}{units}.{part:0{digits}d}"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path,
                        help="directory containing the two private copies and reproduced gold.json")
    parser.add_argument("--limit", type=int, default=24)
    parser.add_argument("--saturation", type=int, default=7)
    parser.add_argument("--side", choices=("lower", "upper"), action="append")
    parser.add_argument("--extra-lower", type=F, default=F(0), metavar="TICKS")
    parser.add_argument("--extra-upper", type=F, default=F(0), metavar="TICKS")
    args = parser.parse_args()
    if min(args.extra_lower, args.extra_upper) < 0:
        raise ValueError("extra rewards must be positive certified ticks per core")
    receipt = authenticated_receipt(args.work / "gold.json")
    prior = tuple(F(receipt[side]["certified_exact"]) for side in ("lower", "upper"))
    birth = birth_model(*prior, args.saturation)
    lp = birth.build(args.limit, sided=True)
    print("census cutoff", args.limit, "overflow offspring", args.saturation,
          "conservation rows", len(lp.rows), "columns", len(lp.cols), flush=True)
    for side in args.side or ("lower", "upper"):
        gold = receipt[side]
        sign = 1 if side == "lower" else -1
        added = args.extra_lower if sign == 1 else args.extra_upper
        nominal = F(gold["nominal"])
        with closing(Comparison(args.work/f"quintet_{side}.bin", side, gold)) as weights:
            cells = census_cells(weights, args.limit, args.saturation, gold)
            print(side, "weight counts", weights.weight_counts,
                  "physical host comparisons", physical_comparison(weights), flush=True)
            def price(fn):
                def objective(kind, meta):
                    if kind == "occ":
                        d, k, *_ = meta
                        return fn(cells[d,k])
                    return 0
                lower = lp.lower(objective)
                return max(F(0), prior[0] * lower)

            scale = weights.scale
            den_lo = weights.global_low + scale*price(lambda r: F(r[0]-weights.global_low, scale))
            den_hi = weights.global_high - scale*price(lambda r: F(weights.global_high-r[1], scale))
            reserve = weights.higher_cap - 10000*price(lambda r: F(weights.higher_cap-r[2], 10000))
            net = added - gold["cap_f"] - gold["cap_g"] - reserve
            denominator = den_hi if net >= 0 else den_lo
            offset = net / denominator
            def additional_gain(row):
                den_gain = (weights.global_high-row[1] if offset >= 0
                            else row[0]-weights.global_low)
                return (weights.higher_cap - row[2] + abs(offset)*den_gain) / F(10000)
            start_den = weights.global_high if offset >= 0 else weights.global_low
            margin = added - gold["err_per_core"] - offset*start_den + 10000*price(additional_gain)
            # A positive stationary excess at the candidate trial also improves
            # the endpoint at least by its size divided by the same upper bound.
            joint_offset = offset + max(F(0), margin)/den_hi
            endpoint = nominal + sign*joint_offset
            weights.check_identity()
            print(side, "stationary denominator lower/upper in ticks", den_lo, den_hi, flush=True)
            print(side, "scaled denominator", outward(den_lo/scale,"lower",12),
                  outward(den_hi/scale,"upper",12), "residual higher cap", reserve,
                  "outward ticks", ceil(reserve), flush=True)
            print(side, "extra ticks", added, "separate offset", offset,
                  "joint margin ticks", margin, "exact endpoint", endpoint,
                  "outward", outward(endpoint,side), flush=True)


if __name__ == "__main__":
    main()

