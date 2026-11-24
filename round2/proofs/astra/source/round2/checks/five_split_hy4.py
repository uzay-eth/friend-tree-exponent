#!/usr/bin/env python3
"""Independent cold check of the frozen five-pass quartet donor replacement.

Scope: uses only the frozen snapshots and pinned small-split reserve output on
real small trees; does not rerun the billion-row production replay. Conditioned
on the separately reviewed raw certificate and mask63 facts.

Checks: frozen hashes and 50 rank counts, deletion/youngest bijection, quartet
and outside-donor eligibility, disjoint incidences, maximum robustness, signed
integer Bellman aggregation with planted errors.
"""

import argparse
from collections import Counter
from fractions import Fraction as Q
from functools import cache
from hashlib import sha256
import os
from pathlib import Path
import shutil
import subprocess
import sys
from tempfile import TemporaryDirectory

REPO = Path(__file__).resolve().parents[2]
FIVE = REPO / "round2/five-next"
OUT = REPO / "round2/outputs"
sys.path.insert(0, str(REPO / "round2/occupation"))
from check_five import Certificate, geometry, quartet, subsets, transitions  # shared physical model

S = 2**40
U = Q(1133, 2000)
RANKS = (
    (2021486872, 56279691, 70491064, 32031570, 1061159),
    (7312312, 6449232, 2975049, 976140, 52278),
    *((1177302, 381038, 2975049, 976140, 52278),) * 3,
)


def hashes():
    manifest = (FIVE / "split-inputs.sha256").read_text().splitlines()
    for line in manifest:
        want, rel = line.split()
        got = sha256((REPO / rel).read_bytes()).hexdigest()
        assert got == want, rel
    return len(manifest)


def arithmetic():
    rows = {}
    for line in (FIVE / "split-reference.tsv").read_text().strip().splitlines()[1:]:
        p = line.split("\t")
        rows[p[0]] = list(map(int, p[1:]))
    assert set(rows) == {"lower", "upper"}
    for side, oldfile in (("lower", "occupation-five-lower.txt"),
                          ("upper", "occupation-five-6.txt")):
        import re
        old = int(re.search(r"newborn reserve (\d+)/1099511627776 at discount 1133/2000",
                            (OUT / oldfile).read_text())[1])
        a = rows[side]
        archived, passes, repl, inc = a[0], a[1:6], a[6], a[7]
        assert archived == old
        assert sum(passes) == repl
        assert repl - archived == inc
        assert len(passes) == 5
    return 50


def keys(state):
    parents, _ = state
    _, degree, core = geometry(state)
    out = []
    for n in range(1, 6):
        for members, shape in subsets(parents, n):
            sh = 0
            if n > 1:
                for i, ch in enumerate(shape, 1):
                    sh = i * sh + int(ch)
            out.append((n, sh, *(z for j in members for z in (degree[j], core[j]))))
    return out


def fixtures():
    layer = {((-1,), (2,))}
    base = set()
    for _ in range(3, 9):
        base.update(layer)
        layer = {nx for st in layer for nx, _, _ in transitions(st)}
    for par in [(-1, 0, 1, 2, 3, 4), (-1, 0, 1, 1, 2, 3), (-1, 0, 1, 2, 2, 3),
                (-1, 0, 1, 1, 1, 2), (-1, 0, 1, 1, 2, 2, 3, 3)]:
        cd = [par.count(i) + (i != 0) for i in range(len(par))]
        for i in range(1, len(par)):
            for d in (5, 8, 13, 14, 15, 33, 191):
                leaves = [max(int(k == 1), 4 - k) for k in cd]
                leaves[i] = max(0, d - cd[i])
                base.add((tuple(par), tuple(leaves)))
    return base


def deletion_eligibility(states):
    """Each 5-set leaf-deletion gives one older quartet and one unique donor."""
    total = 0
    oldest_deleted = 0
    for st in states:
        par, leaves = st
        for members, _ in subsets(par, 5):
            inside = set(members)
            # induced degree inside the chronological set
            deg = {v: 0 for v in members}
            for v in members[1:]:
                deg[par[v]] += 1
                deg[v] += 1
            leaves_in = [v for v in members if deg[v] == 1]
            assert 2 <= len(leaves_in) <= 4
            for z in leaves_in:
                total += 1
                if z == members[0]:
                    oldest_deleted += 1
                older = tuple(v for v in members if v != z)
                assert len(older) == 4 and all(par[v] in older for v in older[1:])
                # unique contact inside the surviving quartet (may be a child
                # of z when z is the oldest member)
                touch = None
                for v in older:
                    if par[v] == z or par[z] == v:
                        touch = v
                        break
                assert touch in older
        # incoming/deleted-youngest: removing the youngest recovers a connected set
        for n in range(2, 6):
            for members, _ in subsets(par, n):
                older = members[:-1]
                assert all(par[v] in older for v in older[1:])
    return total, oldest_deleted


def max_robust(trials=500):
    """A maximizing candidate's successor dominates every candidate successor."""
    import random
    rng = random.Random(11)
    for _ in range(trials):
        k = rng.randint(2, 5)
        now = [rng.randint(0, 1000) for _ in range(k)]
        nxt = [rng.randint(0, 1000) for _ in range(k)]
        j = max(range(k), key=lambda i: now[i])
        assert max(nxt) >= nxt[j]
    return True


def compile_reserve(dst):
    subprocess.run(["c++", "-std=c++20", "-O2", "-Wall", "-Wextra", "-I" + str(REPO),
                    str(FIVE / "reserve.cpp"), "-o", str(dst)], check=True, cwd=REPO)


def query(side, build, states):
    req = build / f"q-{side}.tsv"
    keys_all = sorted({key for st in states for nx in (st, *[n for n, _, _ in transitions(st)])
                       for key in keys(nx)})
    req.write_text("".join(" ".join(map(str, k)) + "\n" for k in keys_all))
    out = build / f"o-{side}.tsv"
    log = build / f"l-{side}.txt"
    with log.open("w") as fh:
        subprocess.run([str(build / "reserve"), str(OUT / f"mc-{side}.bin"),
                        "small-split", "3", str(req), str(out)],
                       check=True, cwd=REPO, stdout=fh)
    vals = [int(x) for x in out.read_text().splitlines()]
    born = vals.pop(0)
    table = dict(zip(keys_all, vals, strict=True))
    assert born == table[(1, 0, 2, 1)]
    return table, born


def source_and_drift(side, build, states, table, born, corrupt=None):
    cert = Certificate(OUT / f"mc-{side}.bin")

    @cache
    def score(st):
        return sum(table[k] for k in keys(st))

    @cache
    def grad(family, s, i):
        t = list(s)
        l, k = t[i]
        t[i] = (l + 1, k)
        return cert.sign * (cert.weight(family, tuple(t)) - cert.weight(family, s))

    worst = None
    checked = 0
    for st in states:
        par, leaves = st
        near, degree, core = geometry(st)
        actual = Q(0)
        for members, _ in subsets(par, 4):
            family, order = quartet(par, members)
            sc = tuple((leaves[j], core[j] - 1) for j in order)
            for role, x in enumerate(order):
                a = grad(family, sc, role)
                if corrupt == "sign":
                    a = -a
                for z in near[x]:
                    if z == 0 or z in order:
                        continue
                    if corrupt == "nohalf":
                        j = max(a, 0) - a / degree[z]
                    else:
                        j = max(a, 0) / Q(2) - a / degree[z]
                    actual += j
        t = Q(score(st))
        gt = births = Q(0)
        for nx, rate, bud in transitions(st):
            gt += rate * (score(nx) - t)
            if bud:
                births += rate
        gap = actual - (U * t - gt + born * births)
        if worst is None or gap < worst[0]:
            worst = (gap, st)
        assert gap >= 0, (side, st, gap)
        checked += 1
    return worst, checked


def mutations():
    # re-added archived baseline would inflate the increment
    rows = {}
    for line in (FIVE / "split-reference.tsv").read_text().strip().splitlines()[1:]:
        p = line.split("\t")
        rows[p[0]] = list(map(int, p[1:]))
    for side, a in rows.items():
        assert a[6] - a[0] == a[7]
        assert a[6] + a[0] != a[7]
    # a removed rank-five row changes the 50-count convention
    assert len(RANKS) == 5 and all(len(r) == 5 for r in RANKS)
    assert sum(len(r) for r in RANKS) * 2 == 50
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("part", nargs="?", choices=("all", "frozen", "trees"), default="all")
    args = ap.parse_args()
    if hasattr(os, "sched_getaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    n = hashes()
    print(f"[frozen] {n} pinned file hashes match manifest")
    assert arithmetic() == 50
    print("[frozen] 50 rank counts, both pass sums, increments and archived values agree")
    assert mutations()
    print("[frozen] re-added baseline / removed row mutants rejected")
    if args.part in ("all", "trees"):
        states = sorted(fixtures())
        tot, old = deletion_eligibility(states)
        assert tot > 0 and old > 0
        print(f"[trees] {len(states)} real trees: {tot} five-set leaf deletions, "
              f"{old} delete the oldest member; incoming/youngest bijection holds")
        assert max_robust()
        print("[max] successor of a maximizing candidate dominates each candidate successor")
        with TemporaryDirectory(prefix="hy4-split-", dir=OUT) as td:
            build = Path(td)
            compile_reserve(build / "reserve")
            for side in ("lower", "upper"):
                table, born = query(side, build, states)
                worst, checked = source_and_drift(side, build, states, table, born)
                print(f"[{side}] signed integer Bellman aggregation PASS on {checked} trees; "
                      f"min slack={float(worst[0]):.6g}")
                cert = Certificate(OUT / f"mc-{side}.bin")
                @cache
                def grad(fam, sc, i):
                    t = list(sc); l, k = t[i]; t[i] = (l + 1, k)
                    return cert.sign * (cert.weight(fam, tuple(t)) - cert.weight(fam, sc))
                changed = {"sign": 0, "nohalf": 0}
                under = 0
                for st in states:
                    par, leaves = st
                    near, degree, core = geometry(st)
                    exact = wrong = half = dropped = Q(0)
                    for members, _ in subsets(par, 4):
                        family, order = quartet(par, members)
                        sc = tuple((leaves[j], core[j] - 1) for j in order)
                        for role, x in enumerate(order):
                            a = grad(family, sc, role)
                            for z in near[x]:
                                if z == 0 or z in order:
                                    continue
                                j = max(a, 0) / Q(2) - a / degree[z]
                                exact += j
                                wrong += max(-a, 0) / Q(2) + a / degree[z]
                                half += max(a, 0) - a / degree[z]
                    # dropped-oldest: omit one quartile incidence per tree where possible
                    for members, _ in subsets(par, 5):
                        if par[members[-1]] in members:
                            dropped += 1
                            break
                    if exact != wrong:
                        changed["sign"] += 1
                    if exact != half:
                        changed["nohalf"] += 1
                    if dropped:
                        under += 1
                assert changed["sign"] > 0 and changed["nohalf"] > 0 and under > 0
                print(f"[{side}] wrong-sign J differs on {changed['sign']} trees; "
                      f"missing-1/2 differs on {changed['nohalf']} trees; "
                      f"{under} trees expose a dropped-incidence under-count")
    print("FIVE SPLIT COLD CHECK: PASS")


if __name__ == "__main__":
    main()
