#!/usr/bin/env python3
"""Independent exact referee check of the reciprocal-degree occupation extension.

Derivation used here (written from the spec, not copied from solver.cpp/physical.py):

Physical model (shared process, imported from the separately reviewed first
calculation) gives, for a state, events (next_state, rate, is_promotion).

For a connected set of r noninitial cores with knot d_r and core cap m_r at each
position, a physical coordinate (D,c) is represented on the r-th table by

    c > m_r            -> score 0
    D < d_r            -> the single corner (D,c) with weight 1
    D >= d_r           -> mixture: mass t = d_r/D on the finite corner (d_r,c)
                          and mass 1-t on the infinity corner (d_r,c,inf).

The score of a set at a physical tree state is the tensor product of these
one-coordinate mixtures applied to the integer corner table.

Tail gain identity (derived here).  Put H = t*Hd + (1-t)*Hi, t = d/D, and after
a leaf gain t' = d/(D+1).  Then

    H' - H = (t'-t)(Hd-Hi) = d/(D(D+1)) (Hi-Hd) = t * [1/(D+1)] * (Hi-Hd).

The solver checks the corner row with a gain-rate interval containing
a/(D+1), giving the corner term min_{r in interval} r*(Hi-Hd) with the endpoint
chosen by the sign of Hi-Hd.  Since the physical coefficient is exactly
t * (a/(D+1)), the physical gain is >= t*(corner term), which is the affine
lower extension.  The aggregated identity tested on real trees is

    source >= u*T - G + v*B,   u = 1133/2000,  G = sum rate*(T(next)-T),
    B = sum rate over promotions, v = newborn score.

Injected faults: flat interpolation, shifted numerator, inflated newborn, and
a flipped rate-sign endpoint at the row level; each must be rejected.
"""

import argparse
import os
from fractions import Fraction as Q
from itertools import product
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(REPO / "round2" / "occupation"))
from check_five import geometry, subsets, transitions  # shared physical model

# A grid different from the author's physical.py grid, to keep the test independent.
D = [[32], [18, 14], [12, 9, 7], [10, 7, 5, 4]]
C = [[9], [6, 5], [4, 4, 3], [4, 3, 2, 2]]
U = Q(1133, 2000)


def arg(xs):
    return ",".join(str(v) for x in xs for v in x)


def compile_solver(dst):
    subprocess.run(
        ["c++", "-std=c++20", "-O2", "-Wall", "-Wextra", "-fno-fast-math",
         "-fsanitize=undefined", "-fno-sanitize-recover=all",
         str(REPO / "round2/occupation-tail/solver.cpp"), "-lgmpxx", "-lgmp", "-o", str(dst)],
        check=True, cwd=REPO,
    )


def dump(solver, snapshot, side, mask, out):
    env = {**os.environ, "FT2_TAIL_DUMP": str(out)}
    res = subprocess.run(
        [str(solver), str(snapshot), str(side), str(mask), "4", arg(D), arg(C)],
        check=True, cwd=REPO, env=env, text=True, capture_output=True,
    )
    return res.stdout.strip()


def read_rows(path):
    rows = {}
    with path.open() as fh:
        for line in fh:
            shape, *data = line.split()
            *flat, score, source = map(int, data)
            assert len(flat) % 3 == 0
            labels = tuple(zip(flat[0::3], flat[1::3], flat[2::3]))
            key = (shape, labels)
            assert key not in rows
            assert 0 <= score <= 2**40 and source >= 0
            rows[key] = score, source
    return rows


def interp(rows, shape, physical, corrupt=""):
    """My own tensor interpolation; corrupt in {'', 'flat', 'shift'}."""
    n = len(physical)
    knots, caps = D[n - 1], C[n - 1]
    coords = []
    for (dv, cv), (knot, cap) in zip(physical, zip(knots, caps)):
        if cv > cap:
            return Q(0), Q(0)
        if dv < knot:
            coords.append([((dv, cv, 0), Q(1))])
            continue
        if corrupt == "flat":
            w = Q(1)
        elif corrupt == "shift":
            w = Q(knot - 1, dv)
        else:
            w = Q(knot, dv)
        coords.append([((knot, cv, 1), w), ((knot, cv, 2), 1 - w)])

    total = [Q(0), Q(0)]
    for combo in product(*coords):
        w = Q(1)
        for _, m in combo:
            w *= m
        if not w:
            continue
        labels = tuple(s for s, _ in combo)
        vals = rows[(shape, labels)]
        for j in range(2):
            total[j] += w * vals[j]
    return total[0], total[1]


def terms(state, rows, corrupt=""):
    parents, leaves = state
    _, ds, cs = geometry(state)
    score = source = Q(0)
    for n in range(1, 5):
        for nodes, shape in subsets(parents, n):
            v, f = interp(rows, shape, tuple((ds[x], cs[x]) for x in nodes), corrupt)
            score += v
            source += f
    return score, source


def newborn(rows):
    return Q(rows[("-", ((2, 1, 0),))][0])


def trees():
    """Real small friend-trees plus decorated high-degree cores."""
    seen = set()
    layer = {((-1,), (2,))}
    for size in range(3, 9):
        for st in layer:
            seen.add(st)
        if size < 8:
            layer = {nxt for st in layer for nxt, _, _ in transitions(st)}
    for parents in [(-1, 0, 1, 2, 3, 4), (-1, 0, 1, 1, 2, 3),
                    (-1, 0, 1, 2, 2, 3), (-1, 0, 1, 1, 1, 2),
                    (-1, 0, 1, 1, 2, 2, 3, 3)]:
        cd = [parents.count(i) + (i != 0) for i in range(len(parents))]
        for tag in range(12):
            degrees = (2, 3, 7, 8, 9, 11, 18, 19, 20, 21, 36, 10**9 + 3)
            leaves = tuple(max(c, 2, degrees[(7 * i + tag + tag // (1 + i)) % len(degrees)]) - c
                           for i, c in enumerate(cd))
            seen.add((tuple(parents), leaves))
    return sorted(seen)


def aggregated(rows, st, corrupt="", born=None):
    initial, source = terms(st, rows, corrupt)
    v = newborn(rows) if born is None else born
    g = Q(0)
    b = Q(0)
    for nxt, rate, prom in transitions(st):
        g += rate * (terms(nxt, rows, corrupt)[0] - initial)
        if prom:
            b += rate
    required = U * initial - g + v * b
    return source, required, b


def tail_identity(rows):
    """Check the derived identity on tail rows and test two corrupted weights."""
    checked = {"": 0, "flat": 0, "shift": 0}
    broken = {"flat": 0, "shift": 0}
    for (shape, labels), (near, _) in rows.items():
        tail = [i for i, lab in enumerate(labels) if lab[2] == 1]
        if len(tail) != 1:
            continue
        pos = tail[0]
        knot, c, _ = labels[pos]
        far_labels = list(labels)
        far_labels[pos] = knot, c, 2
        far = rows[(shape, tuple(far_labels))][0]
        if near == far:
            continue
        dv_list = [knot, knot + 1, 2 * knot, 10**9 + 3]
        for dv in dv_list:
            rate = dv - c + Q(1, 3)
            base = [(di, ci) for di, ci, _ in labels]
            base[pos] = dv, c
            nxt = [(di, ci) for di, ci, _ in labels]
            nxt[pos] = dv + 1, c

            def measured(corrupt):
                a = interp(rows, shape, tuple(base), corrupt)[0]
                b = interp(rows, shape, tuple(nxt), corrupt)[0]
                return rate * (b - a)

            reference = Q(knot, dv) * rate / (dv + 1) * (far - near)
            if measured("") != reference:
                raise AssertionError("tail identity failed at the true weight")
            checked[""] += 1
            for corrupt in ("flat", "shift"):
                checked[corrupt] += 1
                if measured(corrupt) != reference:
                    broken[corrupt] += 1
    return checked, broken


def row_min_sign(rows, shape="-"):
    """Show the min-endpoint rule is necessary: flipping endpoints would accept
    an oversized score that the correct residual rejects."""
    bad = []
    samples = 0
    correct_ok = 0
    for (sh, labels), (v, src) in rows.items():
        if sh != shape:
            continue
        d, c, kind = labels[0]

        def val(dd, cc, kk):
            return rows.get((shape, ((dd, cc, kk),)), (0, 0))[0]

        ranges = []
        if kind == 0:
            L = d - c
            dg = val(d + 1, c, 0) if (d + 1, c, 0) in rows else val(d + 1, c, 1)
            ranges.append((Q(L), Q(L) + Q(c, 2), dg))
            if c < d:
                ranges.append((Q(d - c, d), Q(d - c, d), val(d, c + 1, 0)))
        elif kind == 1:
            ranges.append((Q(d - c, d + 1), Q(1), val(d, c, 2)))
            ranges.append((Q(d - c, d), Q(1), val(d, c + 1, 1)))
        else:
            ranges.append((Q(d - c, d), Q(1), val(d, c + 1, 2)))

        def residual(vv, flip):
            r = Q(src) - U * vv
            for lo, hi, dest in ranges:
                jump = dest - vv
                endpoint = (hi if jump >= 0 else lo) if flip else (lo if jump >= 0 else hi)
                r += endpoint * jump
            return r

        if not ranges:
            continue
        samples += 1
        if residual(v, False) < 0:
            continue
        correct_ok += 1
        for delta in range(1, 4000):
            vv = v + delta
            if residual(vv, False) < 0 <= residual(vv, True):
                bad.append((labels, v, vv))
                break
    return samples, correct_ok, bad


def run(side):
    snapshot = REPO / "round2/outputs/mc-upper.bin" if side == "upper" else REPO / "round2/outputs/mc-lower.bin"
    with TemporaryDirectory(dir=REPO / "round2/outputs", prefix="hy4-") as td:
        work = Path(td)
        solver = work / "solver"
        compile_solver(solver)
        summary = dump(solver, snapshot, 1 if side == "upper" else 2, 63, work / "rows")
        rows = read_rows(work / "rows")
        print(f"[{side}] solver: {summary}")

        # (1) aggregated physical inequality on real trees, exact Fractions.
        sts = trees()
        worst = None
        for st in sts:
            s, req, b = aggregated(rows, st)
            slack = s - req
            if worst is None or slack < worst[0]:
                worst = (slack, st)
            if slack < 0:
                raise AssertionError(f"aggregated inequality failed at {st}: {slack}")
        print(f"[{side}] aggregated inequality PASS on {len(sts)} real trees; "
              f"min slack={float(worst[0]):.6g} at {worst[1]}")

        # (2) derived tail identity + corrupted weights.
        checked, broken = tail_identity(rows)
        print(f"[{side}] tail identity: {checked['']} exact checks; "
              f"flat broken {broken['flat']}/{checked['flat']}, "
              f"shift broken {broken['shift']}/{checked['shift']}")
        assert checked[""] > 0 and broken["flat"] > 0 and broken["shift"] > 0

        # (3) injected faults against the TRUE source on real trees.
        faults = {"flat": 0, "shift": 0, "newborn": 0}
        v = newborn(rows)
        for st in sts:
            strue, _, _ = aggregated(rows, st)
            _, b = None, None
            b = sum(rate for _, rate, p in transitions(st) if p)
            for corrupt in ("flat", "shift"):
                _, req, _ = aggregated(rows, st, corrupt=corrupt)
                if strue < req:
                    faults[corrupt] += 1
            if b and strue < U * terms(st, rows)[0] - sum(
                    rate * (terms(nxt, rows)[0] - terms(st, rows)[0])
                    for nxt, rate, _ in transitions(st)) + 2 * v * b:
                faults["newborn"] += 1
        print(f"[{side}] injected faults rejected on trees: {faults}")
        assert faults["newborn"] > 0, faults
        assert broken["flat"] > 0 and broken["shift"] > 0, broken

        # (4) row-level min-sign endpoint rule.
        samples, correct_ok, bad = row_min_sign(rows)
        print(f"[{side}] row min-sign: {samples} one-vertex rows, correct residual >=0 on {correct_ok}; flipped-endpoint wrong accepts={len(bad)} e.g. {bad[:1]}")
        assert samples > 0 and len(bad) > 0
        return checked, broken, faults, len(bad)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("side", nargs="?", choices=("upper", "lower"), default="upper")
    args = ap.parse_args()
    if hasattr(os, "sched_getaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    run(args.side)
    print("INDEPENDENT REFEREE TEST: PASS")


if __name__ == "__main__":
    main()
