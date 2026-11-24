"""Certify the original gold1 trials with additional edge and three-core rewards."""
import argparse
from contextlib import closing
from fractions import Fraction as Q
import json
from pathlib import Path
import re
import sys


HERE = Path(__file__).resolve().parent
QUINTET = HERE.parents[1]
sys.path.insert(0, str(QUINTET / "stationary"))
import joint

parent = joint.pinned_module("quintet_conditioned_price", HERE.parent / "check.py",
                             "2272cc2397b907ae56bbbf9c53013b25302149f3")
DENOM_HASH = "a9ec77f6c5f399dbd7bd59b510b76b8c72f658948b083a5ff72ce3d0883d4a6a"
OLD_TRIAL = {"lower": Q("0.565703539413800"), "upper": Q("0.565746271349400")}
DONOR = {
    "lower": ("lower-wide.log", "b81568c75a4e57b94db733acc4d9deaac2a4eb568b7c498091b6654316d24ee3"),
    "upper": ("matched.log", "35b7573b85a005166032ff3596282d5ea3cbb07c4d77d1e38e2bbc8932c65e0d"),
}
H3_REPORT = "47c5e956fef39a9a1cd8b64a02e7aeb4aad2c006cdd6607f4019aedb7235f877"


def newly_available(side):
    name, digest = DONOR[side]
    raw = joint.read_pinned(QUINTET / "edges/old-donor" / name, digest).decode()
    section = raw.split(f"{side} floor ", 1)[1]
    donor = Q(re.search(r"^TOTAL (\d+/\d+) ", section, re.M)[1])
    report = joint.read_pinned(joint.R2 / "checks/quintet-h3-cohort-sol.md", H3_REPORT).decode()
    [h3] = re.findall(rf"^{side} (\d+/\d+)$", report, re.M)
    return donor, Q(h3)


def stationary_upper(side):
    receipt = joint.read_pinned(QUINTET / "stationary/joint24.log", DENOM_HASH).decode()
    [upper] = re.findall(rf"^{side} joint H savings.* stationary denominator in ticks \S+ (\S+)$",
                         receipt, re.M)
    return Q(upper)


def certify(lp, side, scale, num, den, original, fixed_lower, search):
    path = HERE / f"{side}.json"
    if search:
        optimum = parent.ratio_proposal(lp, num, den, scale)
        trial = Q((optimum*10**13).__floor__()-2, 10**13)
    else:
        saved = json.loads(path.read_bytes())
        assert saved["side"] == side
        trial = Q(saved["trial"])

    costs = [n-trial*d for n,d in zip(num, den, strict=True)]
    dual = parent.suggest_dual(lp, costs) if search else saved["integer_dual"]
    margin = parent.signed_margin(lp, costs, dual, fixed_lower)
    sign = 1 if side == "lower" else -1
    endpoint = Q(original["nominal"]) + sign*trial
    assert trial > 0 and margin > 0
    if search:
        path.write_text(json.dumps(dict(side=side, trial=str(trial), integer_dual=dual),
                                   separators=(",", ":"))+"\n")
    return endpoint, margin


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("work", type=Path)
    parser.add_argument("--side", choices=("lower","upper"), action="append")
    parser.add_argument("--solve", action="store_true")
    args = parser.parse_args()
    assert __debug__
    joint.read_pinned(Path(joint.__file__),
                      "a9f3677f0b00780ff284a4836042af3355134b2904a14483c7a2f16ce8f55bb8")
    original = joint.authenticated_receipt(args.work/"gold.json")
    lo, hi = [Q(original[s]["certified_exact"]) for s in ("lower", "upper")]
    kernel = joint.pinned_module("quintet_source_kernel", joint.R2/"checks/stationary-p4-sol.py",
                                 "b50cc2d2e51a6e7ab754e2d851abdd690cf20922")
    assert hi <= kernel.HI
    lp = joint.host.birth_model(lo, hi, 7).build(24, sided=True, kernel=kernel.class_kernel())
    assert (len(lp.rows),len(lp.cols)) == (613,5554)
    roots = parent.conditioned_roots()
    for side in args.side or ("lower","upper"):
        data = original[side]
        with closing(joint.Comparison(args.work/f"quintet_{side}.bin", side, data)) as weights:
            scale = weights.scale
            assert scale == 2**40
            sign = 1 if side == "lower" else -1
            donor, h3 = newly_available(side)
            assert min(donor, h3)>0
            upper = stationary_upper(side)
            edge_endpoint = OLD_TRIAL[side] + sign*scale*donor/upper
            both_endpoint = edge_endpoint + sign*scale*h3/upper
            cells = joint.host.census_cells(weights, 24, 7, data)
            prices = parent.edge_floors(scale)[side] + [scale*donor, scale*h3]
            num, den = parent.column_costs(lp, cells, data, weights.higher_cap, roots[side], prices)
            endpoint, margin = certify(lp, side, scale, num, den, data, lo, args.solve)
            assert sign*(endpoint-both_endpoint)>0
            print(side, "old donor only", joint.host.outward(edge_endpoint,side),
                  "both using stationary upper", joint.host.outward(both_endpoint,side))
            print(side, "new normalized floors", donor, h3, "endpoint",endpoint,
                  "outward",joint.host.outward(endpoint,side),"exact eta margin",margin,flush=True)
            weights.check_identity()
    print("PASS gold1 composition on every original cut-24 column")


if __name__ == "__main__":
    main()
