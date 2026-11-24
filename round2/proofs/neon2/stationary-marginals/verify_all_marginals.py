#!/usr/bin/env python3
"""
Batch zero-floating-point exact rational verifier for all 8 FT2_M5V1 stationary-marginal
observables (p_2, p_3, p_4, p_5, p_6, pi(J), pi(2,2), pi(L_o=0)) plus exact rational
corollaries (strict monotonicity p_2 > p_3 > p_4 > p_5 > p_6, degree-2 decomposition
consistency pi(2,1) + pi(2,2) = p_2, tail mass pi(D_o >= 7), bare degree>=3 mass, and
conditional neighbor-reciprocal expectations E_pi[1/D_parent | J] and
E_pi[1/D_v1 + 1/D_v2 | (2,2)]).
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import json
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
import cert_marginal as cm


OBS_LABELS = {
    1: ("p2_D2", "p_2 = pi(D_o = 2)"),
    2: ("p3_D3", "p_3 = pi(D_o = 3)"),
    3: ("p4_D4", "p_4 = pi(D_o = 4)"),
    4: ("p5_D5", "p_5 = pi(D_o = 5)"),
    5: ("p6_D6", "p_6 = pi(D_o = 6)"),
    6: ("pi_J_2_1", "pi(J) = pi(D_o = 2, c_o = 1)"),
    7: ("pi_bare_2_2", "pi(2, 2) = pi(D_o = 2, c_o = 2)"),
    8: ("pi_bare_L0_fin", "pi(L_o = 0) = pi(D_o = c_o)"),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-bin", required=True)
    ap.add_argument("--dir", required=True, help="Directory containing obs{1..8}_{lower,upper}.bin")
    ap.add_argument("--mu-lo", type=F, default=F(565713486333, 10**12))
    ap.add_argument("--mu-hi", type=F, default=F(565736986679, 10**12))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--skip-trees", action="store_true")
    ap.add_argument("--out", default=None)
    args = ap.parse_args()

    if not (F(0) < args.mu_lo <= args.mu_hi < F(1)):
        print(f"ERROR: invalid mu bracket [{args.mu_lo}, {args.mu_hi}] (require 0 < mu_lo <= mu_hi < 1)", file=sys.stderr)
        sys.exit(2)

    d = Path(args.dir)
    results = {}
    intervals = {}
    for oid in range(1, 9):
        lo_path = d / f"obs{oid}_lower.bin"
        hi_path = d / f"obs{oid}_upper.bin"
        lo = cm.combine_side(lo_path, args.verify_bin, threads=args.threads, expected_obs=oid,
                             expected_sign=-1, expected_mu_lo=args.mu_lo, expected_mu_hi=args.mu_hi,
                             return_raw=(not args.skip_trees))
        hi = cm.combine_side(hi_path, args.verify_bin, threads=args.threads, expected_obs=oid,
                             expected_sign=1, expected_mu_lo=args.mu_lo, expected_mu_hi=args.mu_hi,
                             return_raw=(not args.skip_trees))
        if not lo["ok"] or not hi["ok"]:
            lo.pop("_raw", None); hi.pop("_raw", None)
            print(json.dumps({"obs_id": oid, "lower": lo, "upper": hi}, indent=2), file=sys.stderr)
            sys.exit(2)
        lv, hv = F(lo["certified_exact"]), F(hi["certified_exact"])
        if lv > hv:
            print(f"ERROR: inverted interval for obs {oid}: {lv} > {hv}", file=sys.stderr)
            sys.exit(2)
        raw_lo = lo.pop("_raw", None)
        raw_hi = hi.pop("_raw", None)
        entry = {
            "obs_id": oid,
            "obs_name": OBS_LABELS[oid][0],
            "description": OBS_LABELS[oid][1],
            "lower_sha256": lo.pop("sha256"),
            "upper_sha256": hi.pop("sha256"),
            "certified_interval_exact": [str(lv), str(hv)],
            "certified_interval_dec": [cm.format_dec(lv, 8, round_up=False), cm.format_dec(hv, 8, round_up=True)],
            "certified_width_dec": cm.format_dec(hv - lv, 8, round_up=True),
            "lower": lo,
            "upper": hi,
        }
        if not args.skip_trees:
            entry["tree_oracle_lower"] = cm.run_tree_oracle(raw_lo, lo)
            entry["tree_oracle_upper"] = cm.run_tree_oracle(raw_hi, hi)
        results[str(oid)] = entry
        intervals[oid] = (lv, hv)

    # Exact rational corollaries
    mono_checks = {}
    all_mono = True
    for k in range(1, 5):
        gap = intervals[k][0] - intervals[k + 1][1]
        ok = gap > 0
        all_mono = all_mono and ok
        mono_checks[f"p{k+1}_gt_p{k+2}"] = {
            "holds": ok,
            "certified_gap_dec": cm.format_dec(gap, 8, round_up=False),
        }

    # Consistency of p_2 = pi(2,1) + pi(2,2)
    sum_2_lo = intervals[6][0] + intervals[7][0]
    sum_2_hi = intervals[6][1] + intervals[7][1]
    deg2_overlap = (max(sum_2_lo, intervals[1][0]) <= min(sum_2_hi, intervals[1][1]))
    if not deg2_overlap:
        print("ERROR: pi(2,1) + pi(2,2) interval does not overlap p_2 interval!", file=sys.stderr)
        sys.exit(2)

    # Tail mass pi(D_o >= 7) = 1 - (p_2 + ... + p_6)
    tail7_lo = F(1) - sum(intervals[k][1] for k in range(1, 6))
    tail7_hi = F(1) - sum(intervals[k][0] for k in range(1, 6))

    # Bare cores of degree >= 3: pi(L_o = 0, D_o >= 3) = pi(L_o = 0) - pi(2, 2)
    bare_ge3_lo = intervals[8][0] - intervals[7][1]
    bare_ge3_hi = intervals[8][1] - intervals[7][0]

    # Exact stationary identity at J = (2, 1): E_pi[1/D_parent | J] = mu * (1/pi(J) - 1) - 3/2
    j_lo, j_hi = intervals[6]
    if not (F(0) < j_lo <= j_hi < F(1)):
        print(f"ERROR: invalid pi(J) interval [{j_lo}, {j_hi}]", file=sys.stderr)
        sys.exit(2)
    r_parent_J_lo = args.mu_lo * (F(1) / j_hi - F(1)) - F(3, 2)
    r_parent_J_hi = args.mu_hi * (F(1) / j_lo - F(1)) - F(3, 2)

    # Exact stationary identity at (2, 2): E_pi[1/D_v1 + 1/D_v2 | (2, 2)] = pi(J)/(2*pi(2,2)) - mu
    b_lo, b_hi = intervals[7]
    if not (b_lo > 0):
        print(f"ERROR: nonpositive pi(2,2) lower bound {b_lo}", file=sys.stderr)
        sys.exit(2)
    r_nbrs_22_lo = j_lo / (2 * b_hi) - args.mu_hi
    r_nbrs_22_hi = j_hi / (2 * b_lo) - args.mu_lo

    for label, lo_c, hi_c in (
        ("pi_D_ge_7", tail7_lo, tail7_hi),
        ("pi_bare_D_ge_3", bare_ge3_lo, bare_ge3_hi),
        ("E_invD_parent_given_J", r_parent_J_lo, r_parent_J_hi),
        ("E_sum_invD_nbrs_given_2_2", r_nbrs_22_lo, r_nbrs_22_hi),
    ):
        if lo_c > hi_c:
            print(f"ERROR: inverted derived interval for {label}: {lo_c} > {hi_c}", file=sys.stderr)
            sys.exit(2)

    out = {
        "mu_bracket_exact": [str(args.mu_lo), str(args.mu_hi)],
        "mu_bracket_dec": [cm.format_dec(args.mu_lo, 12, round_up=False), cm.format_dec(args.mu_hi, 12, round_up=True)],
        "observables": results,
        "corollaries": {
            "strict_monotonicity_p2_to_p6": {"all_hold": all_mono, "pairs": mono_checks},
            "deg2_decomposition_overlap": {
                "holds": deg2_overlap,
                "pi_J_plus_pi_22_dec": [cm.format_dec(sum_2_lo, 8, round_up=False), cm.format_dec(sum_2_hi, 8, round_up=True)],
                "p2_dec": results["1"]["certified_interval_dec"],
            },
            "pi_D_ge_7": {
                "exact": [str(tail7_lo), str(tail7_hi)],
                "dec": [cm.format_dec(tail7_lo, 8, round_up=False), cm.format_dec(tail7_hi, 8, round_up=True)],
            },
            "pi_bare_D_ge_3": {
                "exact": [str(bare_ge3_lo), str(bare_ge3_hi)],
                "dec": [cm.format_dec(bare_ge3_lo, 8, round_up=False), cm.format_dec(bare_ge3_hi, 8, round_up=True)],
            },
            "E_invD_parent_given_J": {
                "formula": "mu / pi(J) - mu - 3/2",
                "exact": [str(r_parent_J_lo), str(r_parent_J_hi)],
                "dec": [cm.format_dec(r_parent_J_lo, 8, round_up=False), cm.format_dec(r_parent_J_hi, 8, round_up=True)],
            },
            "E_sum_invD_nbrs_given_2_2": {
                "formula": "pi(J) / (2 * pi(2, 2)) - mu",
                "exact": [str(r_nbrs_22_lo), str(r_nbrs_22_hi)],
                "dec": [cm.format_dec(r_nbrs_22_lo, 8, round_up=False), cm.format_dec(r_nbrs_22_hi, 8, round_up=True)],
            },
        },
    }
    text = json.dumps(out, indent=2)
    if args.out:
        Path(args.out).write_text(text + "\n")
    print(text)


if __name__ == "__main__":
    main()
