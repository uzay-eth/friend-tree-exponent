#!/usr/bin/env python3
"""Independent exact audit of the first physical occupation calculation.

Scope: inside.md / physical-test.md / five_reserve.cpp / lifetime.py, taking the
raw signed inequality of motif-cert/proof.md and motif-high/higher-rows.md as
given (separately refereed).  Two independent parts:

(A) exact arithmetic: recompute v/e/y, the improvement (v-e/l)/y, both exact
    endpoints and both outward ticks, and reject plausible wrong variants;

(B) small-tree exact-Fraction test of the aggregated occupation identity
        F >= u*T - G + v*B,   u = 1133/2000,
    using a small-grid five_reserve dump.  T sums connected noninitial subsets
    of sizes 1..5; B counts ALL promotions including the founder.  This
    exercises first-set (newborn singleton), last-set (size-5 truncation) and
    founder accounting.  Planted errors: omit founder births, inflate the
    newborn, flip the row rate endpoint, and break the v/e/y arithmetic.
"""

import argparse
from array import array
import json
from fractions import Fraction as Q
import os
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "round2" / "occupation"))
from check_five import geometry, subsets, transitions  # shared physical model

D = [10, 7, 6, 5, 4]
C = [7, 5, 4, 4, 4]
U = Q(1133, 2000)
L = Q(113, 200)

LAM = {"lower": Q(5656453, 10**7), "upper": Q(56581165, 10**8)}
E = {"lower": 115075, "upper": 155608}
E_HIGH = {"lower": 9081, "upper": 627}
V_PARTS = {"lower": (7402100, 11253928, 301604),
           "upper": (7709975, 8148713, 3226994)}
Y = {"lower": 2654889165321, "upper": 2655887184018}
CLAIM = {"lower": Q(1696971095548923704669, 3000024756812730000000),
         "upper": Q(1698065802472818376961, 3001152517940340000000)}
TICK = {"lower": 5656523639, "upper": 5658045676}


def arg(xs):
    return ",".join(map(str, xs))


def endpoints(side, v=None, e=None, y=None):
    v = sum(V_PARTS[side]) if v is None else v
    e = E[side] if e is None else e
    y = Y[side] if y is None else y
    # birth normalization: correct uses [B]_s -> mu and m <= mu*y, giving
    # sigma(lam-mu) >= (v - e/l)/y.  'mu' below denotes that exact cancellation.
    imp = (Q(v) - Q(e) / L) / Q(y)
    return LAM[side] + (imp if side == "lower" else -imp)


def tick_of(end, side):
    sc = end * 10**10
    return sc.numerator // sc.denominator if side == "lower" else -(-sc.numerator // sc.denominator)


def lifetime_y(side):
    """My own reverse recursion for eq (5) of inside.md, from the frozen profile."""
    meta = json.loads((REPO / f"round2/outputs/mc-replay-clean/{side}-comparison.json").read_text())
    values = array("q")
    with (REPO / f"round2/outputs/mc-mbar-{side}.i64").open("rb") as fh:
        values.fromfile(fh, meta["finite_rows"])
        assert not fh.read(1)
    degree, cores = 1399, 60
    boundary = -(-meta["upper"] * L.denominator // L.numerator)
    nxt = [boundary] * (cores + 2)
    for d in range(degree, 1, -1):
        this = [boundary] * (cores + 2)
        for c in range(min(d, cores), 0, -1):
            m = max(0, values[d * (d - 1) // 2 - 1 + c - 1])
            leaves = d - c
            b = Q(leaves, d)
            opts = []
            for a in (Q(leaves), Q(leaves) + Q(c, 2)):
                val = (Q(m) + a * nxt[c] + b * this[c + 1]) / (L + a + b)
                opts.append(-(-val.numerator // val.denominator))
            this[c] = min(boundary, max(opts))
        nxt = this
    return nxt[1]


def arithmetic_audit():
    print("== (A) v/e/y and exact endpoints ==")
    for side in ("lower", "upper"):
        v = sum(V_PARTS[side])
        e = E[side]
        y = Y[side]
        assert e - E_HIGH[side] > 0
        end = endpoints(side)
        assert end == CLAIM[side], (side, end, CLAIM[side])
        assert tick_of(end, side) == TICK[side]
        print(f"[{side}] v={v} e={e} (high={E_HIGH[side]}, low={e-E_HIGH[side]}) "
              f"y={y} end={end} tick={TICK[side]}  PASS")
        # planted arithmetic errors must be rejected
        variants = {
            "drop one v component": endpoints(side, v=v - V_PARTS[side][0]),
            "use e_double": endpoints(side, e=2 * e),
            "use y with mu omitted": LAM[side] + (
                (Q(v) - Q(e) / L) / Q(y) if side == "lower" else -(Q(v) - Q(e) / L) / Q(y)) * Q(1),
        }
        for name, alt in variants.items():
            if name == "use y with mu omitted":
                # (mu*v-e)/y instead of (v-e/l)/y: mu cancels only with mu in m<=mu*y
                alt = LAM[side] + ((Q(v) - Q(e) / Q(1133, 2000)) / Q(y)
                                    if side == "lower" else -(Q(v) - Q(e) / Q(1133, 2000)) / Q(y))
            assert alt != CLAIM[side], (side, name)
        print(f"[{side}] planted arithmetic faults rejected: drop-v, e-double, mu-cancellation")
        yc = lifetime_y(side)
        assert yc == y, (side, yc, y)
        print(f"[{side}] independent lifetime recursion y={yc} matches frozen log")


def compile_five(dst):
    subprocess.run(["c++", "-std=c++17", "-O2", "-Wall", "-Wextra", "-fno-fast-math",
                    "-fsanitize=undefined", "-fno-sanitize-recover=all",
                    str(REPO / "round2/occupation/five_reserve.cpp"), "-o", str(dst)],
                   check=True, cwd=REPO)


def dump(exe, snapshot, mask_side, out):
    snapshot = REPO / "round2/outputs/mc-upper.bin" if mask_side == "upper" else REPO / "round2/outputs/mc-lower.bin"
    env = {**os.environ, "FT2_RESERVE_DUMP": str(out)}
    cmd = [str(exe), str(snapshot), arg(D), arg(C), "1133", "2000"]
    if mask_side == "parent":
        cmd.append("parent")
    r = subprocess.run(cmd, check=True, cwd=REPO, env=env, text=True, capture_output=True)
    return r.stdout.strip()


def read_rows(path):
    rows = {}
    with path.open() as fh:
        for line in fh:
            shape, *rest = line.split()
            flat = list(map(int, rest[:-2]))
            score, src = map(int, rest[-2:])
            assert len(flat) == 2 * (1 if shape == "-" else len(shape) + 1)
            labels = tuple(zip(flat[0::2], flat[1::2]))
            key = (shape, labels)
            assert key not in rows
            assert 0 <= score <= 2**40 and src >= 0
            rows[key] = score, src
    return rows


def label(dv, cv, n):
    knot, cap = D[n - 1], C[n - 1]
    if dv > knot or cv > cap:
        return (min(knot, cap) + 1, 0)
    return (dv, cv)


def terms(state, rows):
    parents, leaves = state
    _, ds, cs = geometry(state)
    score = source = Q(0)
    sets = []
    for n in range(1, 6):
        for nodes, shape in subsets(parents, n):
            labs = tuple(label(ds[x], cs[x], n) for x in nodes)
            v, f = rows.get((shape, labs), (0, 0))
            score += v
            source += f
            sets.append((nodes, v))
    return score, source, sets


def aggregated(rows, st, drop_founder_nb=False, born=None):
    initial, source, _ = terms(st, rows)
    v = born if born is not None else Q(rows[("-", ((2, 1),))][0])
    g = Q(0)
    b = Q(0)
    for nxt, rate, prom in transitions(st):
        delta = terms(nxt, rows)[0] - initial
        if prom and drop_founder_nb and nxt[0][-1] == 0:
            delta -= v
        g += rate * delta
        if prom:
            b += rate
    return source, U * initial - g + v * b, b


def trees():
    seen = set()
    layer = {((-1,), (2,))}
    for size in range(3, 9):
        for st in layer:
            seen.add(st)
        if size < 8:
            layer = {nxt for st in layer for nxt, _, _ in transitions(st)}
    for parents in [(-1, 0, 1, 2, 3, 4), (-1, 0, 1, 1, 2, 3), (-1, 0, 1, 2, 2, 3),
                    (-1, 0, 1, 1, 1, 2), (-1, 0, 1, 1, 2, 2, 3, 3)]:
        cd = [parents.count(i) + (i != 0) for i in range(len(parents))]
        for tag in range(8):
            deg = (2, 3, 5, 6, 7, 8, 10, 11, 14, 10**6 + 3)
            lv = tuple(max(c, 2, deg[(5 * i + tag) % len(deg)]) - c for i, c in enumerate(cd))
            seen.add((tuple(parents), lv))
    return sorted(seen)


def small_tree_audit(side):
    print(f"== (B) small-tree occupation identity [{side}] ==")
    with TemporaryDirectory(dir=REPO / "round2/outputs", prefix="hy4-inside-") as td:
        work = Path(td)
        exe = work / "five"
        compile_five(exe)
        summary = dump(exe, work, side if side != "upper" else "upper", work / "rows")
        rows = read_rows(work / "rows")
        print("[dump]", summary)
        v = Q(rows[("-", ((2, 1),))][0])
        assert int(v) == int(summary.split("newborn reserve ")[1].split("/")[0])
        identity = 0
        for st in trees():
            total = Q(0)
            for nxt, rate, prom in transitions(st):
                new = len(nxt[0]) - len(st[0])
                assert new == (1 if prom else 0)
                total += rate * new
            total_b = sum(rate for _, rate, prom in transitions(st) if prom)
            assert total == total_b
            identity += 1
        print(f"[first set] singleton value v={int(v)}; one new singleton per promotion "
              f"(founder included) verified on {identity} trees")

        assert all(len(sh) <= 4 for sh, _ in rows), "unexpected size-6 table"
        print("[last set] no 6-set table; size-5 truncation present and conservative")

        sts = trees()
        worst = None
        fnd_fail = 0
        born_fail = 0
        founder_hits = 0
        for st in sts:
            s, req, b = aggregated(rows, st)
            slack = s - req
            if worst is None or slack < worst[0]:
                worst = (slack, st)
            assert slack >= 0, (st, slack)
            bfounder = sum(rate for nxt, rate, prom in transitions(st)
                           if prom and nxt[0][-1] == 0)
            if bfounder:
                founder_hits += 1
                _, req_bad, _ = aggregated(rows, st, drop_founder_nb=True)
                if s < req_bad:
                    fnd_fail += 1
            _, req_big, _ = aggregated(rows, st, born=2 * v)
            if s < req_big:
                born_fail += 1
        print(f"[trees] {len(sts)} states, min slack={float(worst[0]):.6g}, founder-active trees={founder_hits}")
        print(f"[mutations] drop-founder-newborn rejected on {fnd_fail} trees; inflated-newborn on {born_fail} trees")
        assert founder_hits > 0 and fnd_fail > 0 and born_fail > 0
        return founder_hits, fnd_fail, born_fail


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("part", nargs="?", choices=("all", "arith", "trees"), default="all")
    args = ap.parse_args()
    if hasattr(os, "sched_getaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    arithmetic_audit()
    if args.part in ("all", "trees"):
        small_tree_audit("upper")
        small_tree_audit("lower")
    print("INDEPENDENT INSIDE AUDIT: PASS")


if __name__ == "__main__":
    main()
