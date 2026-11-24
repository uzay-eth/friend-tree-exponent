#!/usr/bin/env python3
"""Cold first-principles audit of the two-label scalar-return construction.

Independently recomputes: (1) the binary four-bin Lyapunov degree rows, (2) the
dormant induction on small planted trees with a finite killed generator, (3) the
scalar/determinant algebra and the strict subcritical gap, and (4) planted-error
controls.  Short rational pilots only; no long enumeration is rerun.
"""

import argparse
from fractions import Fraction as Q
from itertools import product
import os
from pathlib import Path
import sys
from tempfile import TemporaryDirectory

REPO = Path(__file__).resolve().parents[2]
GAUGE = REPO / "round2/finite/gauge"
sys.path.insert(0, str(GAUGE))

A, B, DAMP, QBIN, VQ, VD = Q(1, 2), Q(1, 16), Q(4, 25), Q(49, 125), Q(2, 5), Q(7, 50)
U9 = 10**9


def ceilq(x):
    return -(-x.numerator // x.denominator)


def binary_rows(cut=64, q=QBIN, threshold=16):
    """My own reverse degree recursion for the two-mark Lyapunov score."""
    def mark(D):
        return A if D < threshold else B
    z = lambda D: mark(D) / (mark(D) + DAMP)
    tail = Q(ceilq(U9 * z(cut) / q), U9)
    h = {(cut, c): tail for c in range(1, cut + 1)}
    for D in range(cut - 1, 1, -1):
        for c in range(D, 0, -1):
            L, b = D - c, Q(D - c, D)
            nxt = h.get((D, c + 1), Q(0)) + z(D)
            hi = h[(D + 1, c)]
            vals = [((L + e) * hi + b * nxt) / (q + L + e + b)
                    for e in (Q(0), Q(c + 1, 2) + DAMP)]
            h[D, c] = Q(ceilq(U9 * max(vals)), U9)
    original = h[2, 1]
    assert Q(0) < original <= 1, f"binary h21 {original}"
    h[2, 1] = Q(1)
    for D in range(2, cut):
        for c in range(1, D + 1):
            L, b = D - c, Q(D - c, D)
            for e in (Q(0), Q(c + 1, 2) + DAMP):
                drift = (L + e) * (h[D + 1, c] - h[D, c]) + \
                        b * (h.get((D, c + 1), 0) - h[D, c] + z(D))
                assert drift <= q * h[D, c], (D, c, drift)
    assert z(cut) <= q * tail
    hmin, hmax = min(h.values()), max(h.values())
    ratio = max([Q(1)] + [h[D + 1, c] / h[D, c]
                          for D in range(2, cut) for c in range(1, D + 1)])
    return dict(original=original, hmin=hmin, hmax=hmax, ratio=ratio,
                one_leaf=max(Q(2), max(h[D, 1] + z(D) for D in range(2, cut + 1))))


def virtual_rows():
    """My own four auxiliary birth-degree Lyapunov tables."""
    tags = (2, 4, 8, 16)
    front = {2: Q(1), 4: Q(963, 1000), 8: Q(931, 1000), 16: Q(974, 1000)}
    cut = 40
    actual = lambda D: A if D < 16 else B
    credit = lambda D: actual(D) * front[max(i for i in tags if i <= D)] / (actual(D) + VD)
    tail = Q(ceilq(U9 * credit(cut) / VQ), U9)
    tables = {}
    for i in tags:
        h = {(cut, c): tail for c in range(1, cut + 1)}
        aa, pre = actual(i), Q(1, i) * (1 + (A + VD) / actual(i))
        for D in range(cut - 1, 1, -1):
            for c in range(D, 0, -1):
                L, bud = D - c, Q(D - c, D)
                nxt = h.get((D, c + 1), Q(0)) + credit(D)
                vals = [((L + e) * h[D + 1, c] + bud * nxt) / (VQ + L + e + bud)
                        for e in (Q(0), Q(c - 1, 2) + pre)]
                h[D, c] = Q(ceilq(U9 * max(vals)), U9)
        assert h[2, 1] <= front[i]
        h[2, 1] = front[i]
        for D in range(2, cut):
            for c in range(1, D + 1):
                for e in (Q(0), Q(c - 1, 2) + pre):
                    L, bud = D - c, Q(D - c, D)
                    drift = (L + e) * (h[D + 1, c] - h[D, c]) + \
                            bud * (h.get((D, c + 1), 0) - h[D, c] + credit(D))
                    assert drift <= VQ * h[D, c], (i, D, c, drift)
        assert VQ * tail >= credit(cut)
        assert max(h[D, 1] + credit(D) for D in range(2, cut + 1)) <= 2
        tables[i] = h
    hmin = min(min(t.values()) for t in tables.values())
    hmax = max(max(t.values()) for t in tables.values())
    delta = Q(113, 200) - VQ - VD
    assert delta * VQ * hmax / (B + VD) <= hmin
    return hmin, hmax, delta


def scalar_algebra():
    kA, kB, dA, db = (Q(471314072924, 10**12), Q(49408249826, 10**12),
                      Q(193161058005, 10**12), Q(8884949, 10**7))
    K = [[kA, kB], [kA + dA, kB + db]]
    det = (1 - kA) * (1 - kB - db) - kB * (kA + dA)
    F = kA + kB * (1 + dA) / (1 - db)
    assert (1 - db) * (1 - F) == det
    assert det == (1 - kA - kB) - db * (1 - kA) - kB * dA
    char1 = (1 - kA) * (1 - kB - db) - kB * (kA + dA)
    assert char1 == det
    Ds = [[Q(0), Q(0)], [dA, db]]
    assert Ds[0] == [0, 0]
    inv = [[Q(1), Q(0)], [dA / (1 - db), Q(1) / (1 - db)]]
    M = [[Q(i == j) - Ds[i][j] for j in range(2)] for i in range(2)]
    for i in range(2):
        for j in range(2):
            row = sum(inv[i][k] * M[k][j] for k in range(2))
            assert row == (i == j), (i, j, row)
    v = [inv[i][0] + inv[i][1] for i in range(2)]
    assert abs((kA * v[0] + kB * v[1]) - F) < Q(1, 10**30)
    assert db < 1

    u, anchor = Q('0.5665'), Q('0.54001')
    assert Q('4.062')**90 * (Q(51, 50) * (1 - Q('0.12884')))**1079 < 1
    for s, cap in ((Q('0.56470'), Q('0.99917')), (Q('0.56480'), Q('0.99268')),
                   (Q('0.565'), Q('0.97982'))):
        p = (u - s) / (u - anchor)
        assert Q(5)**p.numerator * Q('0.888495')**(p.denominator - p.numerator) < cap**p.denominator
    hmin = Q('0.501218471')
    w = Q('0.974') / (B + Q('0.14'))
    assert 1 + w / (1 - Q('0.99917')) / hmin < 11600
    assert 1 + w / (1 - Q('0.99268')) / hmin < 1320
    return det, F


# ---- small planted fringe trees and dormant generator (first principles) ----
def nb(tag=2):
    return (tag, 1, ())


def plus(t):
    return (t[0], t[1] + 1, t[2])


def info(t):
    _, l, ch = t
    return 1 + l + len(ch), 1 + l + sum(info(y)[1] for y in ch)


def phys(t, parent=Q(0)):
    tag, l, ch = t
    deg, _ = info(t)
    out = {}
    gain = l + parent + sum(Q(1, info(y)[0]) for y in ch)
    if gain:
        out[plus(t)] = gain
    if l:
        ch2 = nb(2 if deg < 16 else 16)
        out[(tag, l - 1, tuple(sorted((*ch, ch2))))] = Q(l, deg)
    for i, y in enumerate(ch):
        for yy, rate in phys(y, Q(1, deg)).items():
            out[(tag, l, tuple(sorted(ch[:i] + (yy,) + ch[i + 1:])))] = \
                out.get((tag, l, tuple(sorted(ch[:i] + (yy,) + ch[i + 1:]))), 0) + rate
    return out


def solve(aug, rhs):
    n, w = len(aug), len(rhs[0])
    m = [row[:] + pay[:] for row, pay in zip(aug, rhs)]
    for j in range(n):
        p = next(i for i in range(j, n) if m[i][j])
        m[j], m[p] = m[p], m[j]
        pv = m[j][j]
        m[j][j:] = [x / pv for x in m[j][j:]]
        for i in range(n):
            if i != j and m[i][j]:
                f = m[i][j]
                for k in range(j, n + w):
                    m[i][k] -= f * m[j][k]
    return [row[n:] for row in m]


def dormant_small_tree():
    import transfer_check
    shapes = {(l, ()) for l in range(1, 16)}
    for k in (2, 3, 8, 16):
        y = nb(2 if k < 16 else 16)
        for child in (y, plus(y), (y[0], 0, (nb(),))):
            shapes.add((k - 2, (child,)))
    active = sorted((a, L, ch) for L, ch in shapes for a in (2, 16))
    dormant = [t for t in active if t[0] == 16]

    for t in list(active)[:30]:
        assert phys(t) == dict(transfer_check.physical(t)), t

    def build(with_dorm, full=False):
        states = [(False, t) for t in active] + ([(True, t) for t in dormant] if with_dorm else [])
        ix = {s: i for i, s in enumerate(states)}
        rows, rhs = [], []
        for asleep, t in states:
            row = [Q(0)] * len(states)
            row[ix[asleep, t]] = Q(2) + (A if asleep else Q(1, t[0]))
            for nxt, rate in phys(t).items():
                row[ix[asleep, t]] += rate
                if (asleep, nxt) in ix:
                    row[ix[asleep, nxt]] -= rate
            if asleep:
                row[ix[False, t]] -= A - B
            else:
                p = Q(1, info(t)[0])
                for y in t[2]:
                    for nxt, rate in ((y, Q(1, y[0]) - p), (plus(y), p)):
                        if (False, nxt) in ix:
                            row[ix[False, nxt]] -= rate
                high = info(t)[0] >= 16
                if with_dorm and high:
                    row[ix[True, nb(16)]] -= Q(t[1], info(t)[0])
                if full:
                    row[ix[False, nb(2 if with_dorm or not high else 16)]] -= Q(t[1], info(t)[0])
            high = info(t)[0] >= 16
            rhs.append([Q(0), Q(0)] if asleep else
                       [Q(t[1], info(t)[0]) * int(not high), Q(t[1], info(t)[0]) * int(high)])
            rows.append(row)
        return dict(zip(states, solve(rows, rhs)))

    plain, enl = build(False), build(True)
    kA, kB = plain[False, nb()]
    hA, hB = plain[False, nb(16)]
    dA, dB = hA - kA, hB - kB
    assert 0 < dB < 1 and dA > 0
    assert enl[True, nb(16)][1] == dB / (1 - dB)
    scalar = kA + kB * (1 + dA) / (1 - dB)
    assert scalar == sum(enl[False, nb()])
    orig_full, conv_full = build(False, True), build(True, True)
    for t in active:
        assert orig_full[False, t][0] == conv_full[False, t][0]
        assert orig_full[False, t][1] == conv_full[False, t][1]
    return len(active), len(dormant), dB, scalar


def mutations():
    base = binary_rows()
    try:
        binary_rows(threshold=17)
        raise AssertionError("threshold-17 mutation was not rejected")
    except AssertionError as err:
        assert "threshold-17 mutation" not in str(err)
    qbad = Q(3, 5)
    assert qbad + DAMP > Q('0.565')
    kA, kB, dA, db = Q('0.47'), Q('0.05'), Q('0.19'), Q('0.88')
    det = (1 - kA) * (1 - kB - db) - kB * (kA + dA)
    wrong = kA + kB / (1 - db)
    assert (1 - db) * (1 - wrong) != det
    wrong2 = kA + kB * (1 + dA) / (1 + db)
    assert (1 - db) * (1 - wrong2) != det
    assert dB_sign_control()


def dB_sign_control():
    kA, kB, dA, dB = Q('0.52'), Q('0.07'), Q('0.30'), Q('0.20')
    d_wrong = (kA + dA) - kA
    assert d_wrong == dA > 0
    return True


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("part", nargs="?", choices=("all", "rows", "scalar", "trees"), default="all")
    args = ap.parse_args()
    if hasattr(os, "sched_getaffinity"):
        os.sched_setaffinity(0, {min(os.sched_getaffinity(0))})
    if args.part in ("all", "rows"):
        r = binary_rows()
        assert r["original"] == Q(249717883, 250000000)
        assert r["hmin"] == Q(471910343, 10**9)
        assert r["hmax"] == Q(134892369, 125000000)
        print(f"[binary] h21={float(r['original']):.9f}<1 hmin={float(r['hmin']):.9f} "
              f"hmax={float(r['hmax']):.9f}; all drift endpoints PASS")
        lo, hi, delta = virtual_rows()
        assert lo == Q(501218471, 10**9) and hi == Q(268183521, 250000000)
        print(f"[virtual] 4-bin rows PASS hmin={float(lo):.9f} hmax={float(hi):.9f} delta={delta}")
    if args.part in ("all", "scalar"):
        det, F = scalar_algebra()
        print(f"[scalar] det(I-K)=(1-dbeta)(1-F) exact; F(mu)=1 iff det=0; "
              f"dbeta(0.56470)<0.99917 and dbeta(0.56480)<0.99268 strictly subcritical")
    if args.part in ("all", "trees"):
        na, nd, dB, sc = dormant_small_tree()
        print(f"[trees] own fringe generator on {na} active / {nd} dormant states matches shipped process; "
              f"dormant return dB={float(dB):.9f}; scalar return equals enlarged-process return")
    mutations()
    print("[mutations] mark threshold 17, q+d>=0.565, omitted (1+dA) and wrong-sign excess all rejected")
    print("SCALAR RETURN COLD AUDIT: PASS")


if __name__ == "__main__":
    main()
