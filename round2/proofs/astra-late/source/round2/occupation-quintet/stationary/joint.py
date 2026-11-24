"""Joint stationary bound for the original five-vertex witness, in integer ticks."""

import argparse
from contextlib import closing
from fractions import Fraction as Q
import hashlib
import importlib.util
import json
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
from weights import Comparison, R2, RECEIPT_HASH, authenticated_receipt


def pinned_module(name, path, git_blob):
    data = path.read_bytes()
    if hashlib.sha1(f"blob {len(data)}\0".encode() + data).hexdigest() != git_blob:
        raise ValueError(f"input changed: {path}")
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


host = pinned_module("quintet_host_allocation", HERE / "check.py",
                     "bb1987fb9854c1e6194e84fafef74b32c44ee2cb")
EXPORTS = {
    "lower": "d698df2bdd20fd6d5a802f84d18124a1fcd6f491bb4dde8ae5038ca5d27f6169",
    "upper": "a23fe15ac518c35f9e7dcc92bc1fa497253f16e7677a9d94e572815682253fa8",
}


def read_pinned(path, expected):
    raw = path.read_bytes()
    if hashlib.sha256(raw).hexdigest() != expected:
        raise ValueError(f"input changed: {path}")
    return raw


def reduced_rewards(side, gold, prior):
    path = HERE.parent / "lifetime" / "dual" / f"{side}.json"
    ledger = json.loads(read_pinned(path, EXPORTS[side]))
    assert ledger["schema"] == "quintet-f-dual-v1" and ledger["side"] == side
    assert ledger["source_commit"] == "b04cc50c1e391e0420eb7b006ac781d484c68543"
    assert ledger["receipt_sha256"] == RECEIPT_HASH
    assert ledger["payload_sha256"] == gold["verify"]["sha256"]
    assert ledger["scale"] == 2**40 and ledger["cap_f_ticks"] == gold["cap_f"]
    assert Q(ledger["lambda"]) == Q(gold["nominal"])
    assert (Q(ledger["LO"]), Q(ledger["HI"])) == prior
    grid, cut, arity = ledger["V_grid"], ledger["cut"], ledger["arity"]
    assert (ledger["degree"], arity, cut) == (12000, 30, 64)
    values = {(d,c): Q(v,grid) for d,c,v in ledger["values"]}
    rows = {(d,c): (Q(a),Q(b)) for d,c,a,b in ledger["rows"]}
    expected = {(d,c) for d in range(2,cut+1) for c in range(1,min(d,arity)+1)}
    assert set(rows) == expected and len(ledger["rows"]) == len(expected)
    assert set(values) == expected | {(cut+1,c) for c in range(1,arity+1)}
    assert values[2,1] == Q(ledger["V_J_ticks"])
    for (d,c),(a,b) in rows.items():
        assert min(a, a+b*Q(c,2)) >= 0, (d,c)
    return values[2,1], rows


def edge_prices():
    folder = HERE.parent / "edges"
    receipts = ((folder / "replay.log", "4d2cc123b2e9826db5cc8315e86d283b36b1e1231ed6f4d3adff0bfacf555c0d",
                 "mask degree 18 ", "credit ticks/core "),
                (folder / "fy" / "stock8.log", "deb154ef476a27a643b7c4ff5f1d12cca07f74a474cd989d193c84032632a236",
                 "TOTAL ", "TOTAL "))
    prices = {s: [] for s in EXPORTS}
    for path, digest, start, marker in receipts:
        side = None
        found = {}
        for line in read_pinned(path, digest).decode().splitlines():
            token = line.split(maxsplit=1)
            if token and token[0] in prices:
                side = token[0]
            if line.startswith(start):
                assert side is not None and side not in found
                found[side] = Q(line.split(marker,1)[1].split()[0])
        assert set(found) == set(prices) and all(v > 0 for v in found.values())
        for side, value in found.items():
            prices[side].append(value)
    return prices


def estimate(side, gold, weights, cells, lp, limit, prior, extra):
    newborn, reduced = reduced_rewards(side, gold, prior)
    receive = {}
    for kind, _, meta in lp.cols:
        if kind != "occ":
            continue
        d,k,_,_,child,mother = meta
        # The degree overflow cell represents *every* D>limit. Even if a
        # coefficient is stored at D=limit+1, this cell cannot use it.
        a,b = reduced.get((d,k+1), (0,0)) if d <= limit else (0,0)
        receive[meta] = a + b*(child+mother)
        assert receive[meta] >= 0

    def price(cost, normalization):
        # The group bounds certify the repairs in eta=pi/mu; every cell cost
        # in this calculation is nonnegative, so multiplication by LO is safe.
        def objective(kind, meta):
            if kind != "occ":
                return 0
            d,k = meta[:2]
            value = cost(cells[d,k], receive[meta])
            assert value >= 0
            return Q(value, normalization)
        return max(Q(0), prior[0]*normalization*lp.lower(objective))

    scale = weights.scale
    lo = weights.global_low + price(lambda r,f: r[0]-weights.global_low, scale)
    hi = weights.global_high - price(lambda r,f: weights.global_high-r[1], scale)
    assert 0 < lo <= hi
    combined = price(lambda r,f: weights.higher_cap-r[2]+f, 10**7)
    fixed = prior[0]*newborn + sum(extra) - gold["err_per_core"]
    net = fixed + combined
    offset = net / (hi if net >= 0 else lo)
    # Anchor the affine prediction at its current signed offset. This gives
    # a new lower bound on the same ratio even if the dual has a repair cost.
    base = weights.global_high if offset >= 0 else weights.global_low
    def joint(row, slack):
        improvement = base-row[1] if offset >= 0 else row[0]-base
        return weights.higher_cap-row[2]+slack+abs(offset)*improvement
    margin = fixed - offset*base + price(joint, 10**7)
    increment = max(0, margin) / hi
    signed_offset = offset + increment
    nominal = Q(gold["nominal"])
    direction = 1 if side == "lower" else -1
    endpoint = nominal + direction*signed_offset
    print(side, "f lifetime birth floor", prior[0]*newborn, "g18 and outside FY", *extra)
    print(side, "joint H savings plus reduced f per core", combined,
          "stationary denominator in ticks", lo, hi)
    print(side, "separate offset", offset, "joint trial margin", margin,
          "extra offset", increment, flush=True)
    print(side, "exact endpoint", endpoint, "outward", host.outward(endpoint,side),
          "inward from original", host.outward(direction*(endpoint-Q(gold['certified_exact'])),"lower"), flush=True)
    return endpoint


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path, help="private copies of both original payloads and gold.json")
    parser.add_argument("--limit", type=int, default=24)
    parser.add_argument("--side", choices=("lower", "upper"), action="append")
    parser.add_argument("--without-edges", action="store_true", help="matched f and host calculation only")
    args = parser.parse_args()
    if not 19 <= args.limit <= 64:
        raise ValueError("physical Palm classes and exported rows require 19 <= limit <= 64")
    gold = authenticated_receipt(args.work / "gold.json")
    prior = tuple(Q(gold[s]["certified_exact"]) for s in ("lower", "upper"))
    audit = pinned_module("quintet_physical_birth_kernel", R2 / "checks/stationary-p4-sol.py",
                          "b50cc2d2e51a6e7ab754e2d851abdd690cf20922")
    assert prior[1] <= audit.HI
    sided = host.birth_model(*prior, 7)
    lp = sided.build(args.limit, sided=True, kernel=audit.class_kernel())
    extras = edge_prices() if not args.without_edges else {s:(0,0) for s in EXPORTS}
    print("physical LP", args.limit, "conservation rows", len(lp.rows), "columns", len(lp.cols), flush=True)
    for side in args.side or ("lower", "upper"):
        with closing(Comparison(args.work/f"quintet_{side}.bin", side, gold[side])) as weights:
            cells = host.census_cells(weights, args.limit, 7, gold[side])
            estimate(side, gold[side], weights, cells, lp, args.limit, prior, extras[side])
            weights.check_identity()


if __name__ == "__main__":
    main()
