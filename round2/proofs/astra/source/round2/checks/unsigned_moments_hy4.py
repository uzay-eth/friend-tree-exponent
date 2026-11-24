#!/usr/bin/env python3
"""Standalone hy4 referee check of the V5 UNSIGNED per-newborn denominator.

Scope: the two rational dual sheets moments-v5-unsigned{,-fine}-{lower,upper},
their exact LP price certificates, the physical-tree pointwise potential bound,
the parent/child (Palm/oriented-edge) generator identity, the all-degree tail,
and negative controls that plant endpoint/sign/index faults.  No signed donor,
reserve, tip, star or lower-score file is opened here.

Run from the repository root (one core), e.g.
  python3 -B round2/checks/unsigned_moments_hy4.py unsigned-fine both
  python3 -B round2/checks/unsigned_moments_hy4.py unsigned lower
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
for sub in ("round2/moments", "round2/motif-cert", "round2/motif-high"):
    sys.path.insert(0, str(ROOT / sub))

from fractions import Fraction as F
import copy
import gzip
import json
from bisect import bisect_left
from random import Random
from types import SimpleNamespace

from occupation import Occupation
from marginal_lift import HighDegrees
from v5 import Application, score_inputs, OUT, CHILD, DEGREE, group
from test_occupation import observe, layout, fixture, dot, one_step
from low_oracle import Local
from oracle import Tree

LOW, HIGH = F(113, 200), F(1133, 2000)
UNIT = 1 << 40
CFG = {"unsigned": (5, 80), "unsigned-fine": (7, 144)}


def myverify(model, obj, w):
    """Independent Lagrangian dual lower bound for Ax=b, Gx>=h, 0<=x<=cap.

    eq prices are free; ge prices must be nonnegative; the capped-variable
    correction sum_j cap_j*min(0, c_j-(A^Tp+G^Ty)_j) handles integer rounding.
    """
    scale = F(w["dual_scale"])
    residual, bound = dict(obj), F()
    for tag, rows in (("eq", model.equal), ("ge", model.atleast)):
        seen = set()
        for i, n in w[tag]:
            assert type(i) is int and i not in seen and 0 <= i < len(rows), "bad dual index"
            assert tag == "eq" or n >= 0, "negative inequality price"
            seen.add(i)
            y = F(n, scale)
            _, row, rhs = rows[i]
            bound += rhs * y
            for j, a in row.items():
                residual[j] = residual.get(j, F()) - a * y
    for j, z in residual.items():
        if z < 0:
            bound += z * model.cap[j]
    return bound


def objective(num, den, sign, q):
    obj = {j: sign * v for j, v in num.items()}
    for j, v in den.items():
        obj[j] = obj.get(j, F()) - sign * q * v
    return obj


def build(form, side):
    K, M = CFG[form]
    m = Occupation(K, low=LOW, high=HIGH)
    HighDegrees(m, M)
    profile = OUT / f"moments-v5-{side}-k{K}.tsv.gz"
    app = Application(m, profile, score_inputs(side, profile), True,
                      score_path=OUT / f"moments-v5-{side}-score-k{K}.tsv.gz",
                      edge_score_path=OUT / f"moments-v5-{side}-edge-score-k{K}.tsv.gz",
                      score_only=True)
    sheet = json.load(gzip.open(OUT / f"moments-v5-{form}-{side}-dual.json.gz", "rt"))
    assert sheet["score_only"] and not sheet["use_prior"]
    assert sheet["K"] == K and sheet["M"] == M
    assert set(sheet["rows"]) == {"M-per-core", "M-per-birth"}
    return m, app, sheet


def full_vector(tree, m, idx):
    """Physical occupation x of a cyclic tree, including v5 tail cells and edge aggregates."""
    x, _ = observe(tree, m)
    kids, types = layout(tree)
    parent, leaves = tree
    for i, (d, c) in enumerate(types):
        k, L, b = c - 1, leaves[i], F(leaves[i], d)
        if d > m.lift.end:
            cell = (bisect_left(DEGREE, d), bisect_left(CHILD, k), 16 * L // d)
            ip = F(1, types[parent[i]][0])
            ic = sum((F(1, types[w][0]) for w in kids[i]), F())
            for name, val in [("mass", 1), ("slots", k), ("parent", ip), ("child", ic), ("birth", b)]:
                x[idx["v5-tail", cell, name]] += val
        a = group(types[parent[i]][0], types[parent[i]][1] - 1, m.cutoff)
        b2 = group(d, k, m.cutoff)
        x[idx["v5-edge", a, b2]] += 1
    return x, kids, types


def physical_census(m, app, side, K, M):
    """Pointwise phi <= den <= mbar on real trees, plus one-step generator identity."""
    cert = Local(OUT / f"mc-{side}.bin")
    idx = {k: j for j, k in enumerate(m.names)}
    rng = Random(7562501)
    checked = 0
    for deg in (2, 4, 16, 72, 1400, 1401, 2000003):
        p, l = fixture(rng, deg % 17, True)
        l[0] += max(0, deg - layout((p, l))[1][0][0])
        x, kids, types = full_vector((p, l), m, idx)
        states = [(li, len(ch)) for li, ch in zip(l, kids)]
        phi = sum((cert.f(s) for s in states), F())
        for fam, nodes in Tree.embeddings(SimpleNamespace(ch=kids, p=p)):
            phi += cert.weight(fam, tuple(states[i] for i in nodes))
        phi /= app.unit
        assert phi <= dot(app.den, x) <= dot(app.mbar, x), ("unsigned pointwise bound", deg)
        checked += 1
    base = Occupation(K, low=LOW, high=HIGH)  # plain model: one_step checks the generator, not v5 aggregates
    HighDegrees(base, M)
    for i in range(4):
        one_step(fixture(rng, i, True), base)  # asserts parent/child generator orientation
    return checked


def plant_tests(m, name, num, den, sign, w, q):
    base = myverify(m, objective(num, den, sign, q), w)
    assert base >= 0
    moved = myverify(m, objective(num, den, sign, F(2)), w)
    assert moved < 0, "endpoint plant not detected"
    flip = copy.deepcopy(w)
    flip["eq"][0][1] *= -1
    assert myverify(m, objective(num, den, sign, q), flip) < 0, "sign plant not detected"
    shifted = copy.deepcopy(w)
    shifted["ge"][0][0] += 1
    assert myverify(m, objective(num, den, sign, q), shifted) < 0, "index plant not detected"
    for bad in ("neg", "oob"):
        v = copy.deepcopy(w)
        if bad == "neg":
            v["ge"][0][1] = -abs(v["ge"][0][1])
        else:
            v["eq"][0][0] = 10 ** 9
        try:
            myverify(m, objective(num, den, sign, q), v)
        except AssertionError:
            continue
        raise AssertionError(f"{bad} plant not detected")


def endpoints(bounds, credits, losses, lam):
    out = {}
    for label, total in credits.items():
        ends = {}
        for side in ("lower", "upper"):
            q = (F(total[side], UNIT) - F(losses[side], UNIT) / LOW) / bounds[side]
            ends[side] = lam[side] + q if side == "lower" else lam[side] - q
        out[label] = ends
    return out


def run(form, side_arg):
    sides = ("lower", "upper") if side_arg == "both" else (side_arg,)
    bounds = {}
    print(f"== {form}: exact dual certificates (independent Fractions) ==")
    for side in sides:
        m, app, sheet = build(form, side)
        K, M = sheet["K"], sheet["M"]
        print(f" {side} K={K} M={M} rows-on-original-prior={LOW}..{HIGH}")
        for name in ("M-per-core", "M-per-birth"):
            num, den, sign = app.rows[name]
            w = sheet["duals"][name]
            q = F(w["bound"])
            margin = myverify(m, objective(num, den, sign, q), w)
            assert margin >= 0 and w["direction"] == sign, (side, name)
            print(f"   {name} = {float(q):.12f}  margin {float(margin):.3g}")
            if name == "M-per-birth":
                bounds[side] = q
                plant_tests(m, name, num, den, sign, w, q)
        n = physical_census(m, app, side, K, M)
        print(f"   physical census: {n} trees phi<=den<=mbar, 4 one-step generator identities")
    if sides == ("lower", "upper"):
        lam = {"lower": F(5656453, 10000000), "upper": F(11316233, 20000000)}
        losses = {"lower": 115075, "upper": 155608}
        credits = {"three lifetime reserves": {"lower": 18957632, "upper": 19085682},
                   "three reserves and affine source": {"lower": 34246660, "upper": 36544047}}
        print("== exact per-newborn substitution (ticks in 2^-40) ==")
        for label, ends in endpoints(bounds, credits, losses, lam).items():
            from math import ceil, floor
            print(f"   {label}: [0.{floor(ends['lower']*10**12):012d}, "
                  f"0.{ceil(ends['upper']*10**12):012d}]")
    print(f"PASS {form} {side_arg}")


if __name__ == "__main__":
    if len(sys.argv) != 3 or sys.argv[1] not in CFG or sys.argv[2] not in ("lower", "upper", "both"):
        raise SystemExit("usage: unsigned_moments_hy4.py {unsigned|unsigned-fine} {lower|upper|both}")
    run(sys.argv[1], sys.argv[2])
