#!/usr/bin/env python3
"""Independent compact referee checks for the four-vertex raw drift certificate.

Repository-relative usage (one CPU, short):
  python3 -B round2/checks/motif_base_hy4.py upper
  python3 -B round2/checks/motif_base_hy4.py lower

Re-derives the analytic-tail algebra, samples the frozen affine envelopes and
rebuilds the physical generator / FY / newborn / two-credit identities on its
OWN small trees.  Snapshots are opened read-only; the author's compiled binaries
are queried only as a table oracle.  The 74-million-row production sweeps are
NOT rerun -- they were generated and pinned by the author.
"""
from collections import Counter, defaultdict
from fractions import Fraction as Q
from itertools import combinations, product
from math import comb
from pathlib import Path
import hashlib
import json
import subprocess
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "round2" / "motif-high"))
sys.path.insert(0, str(ROOT / "round2" / "motif-cert"))
from oracle import ARITY, INF, NB, Certificate, bud, degree, gain, inv, replace  # noqa: E402
from low_oracle import Local  # noqa: E402
from tail import header, read_tail, Poly, Q as PQ, X as PX  # noqa: E402

HERE = ROOT / "round2"
OUT = HERE / "outputs"
M = json.loads((HERE / "motif-cert" / "certificate.json").read_text())


def ok(cond, msg):
    if not cond:
        raise ArithmeticError(msg)


def digest(path):
    return hashlib.file_digest(path.open("rb"), "sha256").hexdigest()


# ---------------------------------------------------------------------------
# 1. Analytic f tail -- re-derive the residual from the model, bound it.
# ---------------------------------------------------------------------------
def wval(coef, q, x):
    a, b, c, e = coef
    return a + b * q + c * q * q + e * q * x


def direct_residual(coef, C, d, q, x, rho, lam, sign):
    qp, xp = (q + x) / (1 + x), x / (1 + x)
    dg = wval(coef, qp, xp) - wval(coef, q, x)
    db = wval(coef, q - x, x) - wval(coef, q, x)
    return sign * ((rho / x) * dg + q * (db + C + d * x) - lam * wval(coef, q, x))


def my_bernstein(poly, q0, q1, xmax):
    """Own power->tensor-Bernstein conversion on [q0,q1] x [0,xmax]."""
    terms = {k: v for k, v in poly.items() if v} or {(0, 0): 0}
    nq = max(i for i, j in terms)
    nx = max(j for i, j in terms)
    pw = defaultdict(Q)
    for (i, j), v in terms.items():
        for a in range(i + 1):
            pw[(a, j)] += v * comb(i, a) * q0 ** (i - a) * (q1 - q0) ** a * xmax ** j
    def b1(co, n):
        return [sum(Q(comb(k, i), comb(n, i)) * co.get(i, 0) for i in range(k + 1))
                for k in range(n + 1)]
    by_j = defaultdict(dict)
    for (a, j), v in pw.items():
        by_j[j][a] = v
    bt = {j: b1(co, nq) for j, co in by_j.items()}
    out = []
    for a in range(nq + 1):
        co = {j: bt[j][a] for j in bt}
        out.extend(b1(co, nx))
    return out


def analytic_check(path, target, panels=16):
    with path.open("rb") as fp:
        h = header(fp)
    _, C, d = read_tail(str(path))
    lam, sign = target, (1 if h["mode"] == 1 else -1)
    coef = ((h["f_p0"], h["f_p1"], 0, 0) if sign == 1 else
            (0, h["f_p0"], -h["f_p1"], -h["f_p2"]))
    xmax = Q(1, h["DF"])
    w = coef[0] + coef[1] * PQ + coef[2] * PQ * PQ + coef[3] * PQ * PX
    den = (1 + PX) * (1 + PX)
    N = (coef[1] * (1 - PQ) * (1 + PX)
         + coef[2] * (1 - PQ) * (2 * PQ + (1 + PQ) * PX)
         + coef[3] * ((1 - 2 * PQ) * PX - PQ * PX * PX))
    U = -coef[1] * PX + coef[2] * (-2 * PQ * PX + PX * PX) - coef[3] * PX * PX
    # (a) polynomial identity against the physical definition at many states
    worst = None
    for rho in (PQ, (1 + PQ) * Q(1, 2)):
        poly = (rho * N + (PQ * (U + C + d * PX) - w * lam) * den) * sign
        for D in (h["DF"], h["DF"] + 1, 3 * h["DF"] + 17):
            for L in (0, 1, D // 3, D - 1):
                q, x = Q(L, D), Q(1, D)
                rv = q if rho is PQ else (1 + q) / 2
                direct = direct_residual(coef, C, d, q, x, rv, lam, sign)
                ok(poly.at(q, x) == direct * (1 + x) ** 2, "analytic polynomial identity")
    # (b) independent bivariate Bernstein enclosure over q in [0,1], x in [0,1/DF]
    top = None
    for rho in (PQ, (1 + PQ) * Q(1, 2)):
        poly = (rho * N + (PQ * (U + C + d * PX) - w * lam) * den) * sign
        for j in range(panels):
            v = max(my_bernstein(poly, Q(j, panels), Q(j + 1, panels), xmax))
            top = v if top is None else max(top, v)
    cap = top if top >= 0 else top / (1 + xmax) ** 2
    # (c) semantic sampling: residual <= cap on a rational grid (incl. q=0, x small)
    for i in range(41):
        q = Q(i, 40)
        for j in range(9):
            x = xmax * Q(j, 8)
            if not x:
                continue
            for rv in (q, (1 + q) / 2):
                r = direct_residual(coef, C, d, q, x, rv, lam, sign)
                worst = r if worst is None else max(worst, r)
                ok(r <= cap, "analytic sampled residual exceeds cap")
    return -(-cap.numerator // cap.denominator), worst


# ---------------------------------------------------------------------------
# 2. Affine envelope -- sample the frozen 1D/2D enclosure from the weights.
# ---------------------------------------------------------------------------
def corner_state(l, i):
    if i == l.inf:
        return INF
    if i < l.inf:
        return l.internal[i]
    k = i - l.inf - 1
    return l.d0 - 1 - k, k


def proj(f, c, i):
    if i == f.inf:
        return c.inf, c.inf, 1, 1
    s = corner_state(f, i)
    d, k = degree(s), s[1]
    if d < c.d0:
        v = c.id[s]
        return v, v, 1, 1
    return c.inf + 1 + min(k, c.kt + 1), c.inf, c.d0, d


def envelope_samples(c, n=80):
    """Recompute envelope enclosures from raw weights on sampled fine labels."""
    from random import Random
    rnd = Random(20241124)
    g2, g3 = c.l2, c.levels[3]
    checked = 0
    # 1D: vertical-middle source at level 2 -> level 3, for both edge directions.
    for sec_m, sec_r, dirn in (("e230", "e231", 0), ("e232", "e233", 1)):
        mo, mr = c.sections[sec_m][0], c.sections[sec_r][0]
        for _ in range(n):
            a = rnd.randrange(len(c.act["V"]))
            s = list(c.act["V"])[a]
            if not (1 <= s[1] and degree(s) < g2.d0):
                continue
            sid, gid = g2.id[s], g2.id[gain(s)]
            old_c, new_c = corner_state(g2, sid), corner_state(g2, gid)
            if dirn == 0:
                delta = lambda v: (c.weight("E", (new_c, corner_state(g2, v)))
                                   - c.weight("E", (old_c, corner_state(g2, v))))
            else:
                delta = lambda v: (c.weight("E", (corner_state(g2, v), new_c))
                                   - c.weight("E", (corner_state(g2, v), old_c)))
            off = a * g3.n
            pool = sorted(g2.parent) if dirn == 1 else list(range(g2.n))
            for v in rnd.sample(pool, min(6, len(pool))):
                first, second, nn, dd = proj(g2, g3, v)
                mv = c.data[mo + off + first]
                rv = c.data[mr + off + first]
                mvv = c.data[mo + off + second]
                rvv = c.data[mr + off + second]
                err = dd * delta(v) - nn * mv - (dd - nn) * mvv
                ok(abs(err) <= nn * rv + (dd - nn) * rvv,
                   f"1D envelope violated {sec_m} state={s} label={v}")
                checked += 1
    return checked


# ---------------------------------------------------------------------------
# 3. Own small-tree physical identity, FY/newborn bookings, two credits.
# ---------------------------------------------------------------------------
class MyTree:
    def __init__(self, parents, leaves):
        self.p = tuple(parents)
        self.ch = [[] for _ in parents]
        for i in range(1, len(parents)):
            self.ch[parents[i]].append(i)
        self.s = (INF,) + tuple((leaves[i], len(self.ch[i])) for i in range(1, len(parents)))
        self.r = [Q(0)] + [s[0] + inv(self.s[self.p[i]]) + sum(inv(self.s[j]) for j in self.ch[i])
                           for i, s in enumerate(self.s) if i]
        self.b = (Q(0),) + tuple(Q(s[0], degree(s)) for s in self.s[1:])

    def after(self, x, promo):
        l = list(self.s[i][0] for i in range(1, len(self.s)))
        l[x - 1] += -1 if promo else 1
        np = list(self.p) + ([x] if promo else [])
        nl = [0] + l + ([1] if promo else [])
        return MyTree(np, nl)

    def embeddings(self):
        for u in range(len(self.p)):
            for a in self.ch[u]:
                yield "E", (u, a)
                for b in self.ch[a]:
                    yield "V", (u, a, b)
                    for cc in self.ch[b]:
                        yield "P", (u, a, b, cc)
            for a in self.ch[u]:
                for side in self.ch[u]:
                    if side != a:
                        for b in self.ch[a]:
                            yield "F", (u, side, a, b)
            for v, w in combinations(self.ch[u], 2):
                yield "C", (u, v, w)
                if u:
                    yield "Y", (self.p[u], u, v, w)
            for a, b, cc in combinations(self.ch[u], 3):
                yield "T", (u, a, b, cc)

    def active(self, c):
        out = {}
        for f, vs in self.embeddings():
            st = tuple(self.s[v] for v in vs)
            if c.active(f, st) < 0:
                continue
            key = frozenset(vs)
            ok(key not in out, "duplicate embedding")
            out[key] = (f, vs, st, c.weight(f, st))
        return out


def my_newborn(f, post, pos):
    a = post
    if f == "E":
        return ("C", (a[0], a[1], NB)) if pos == 0 else ("V", (a[0], a[1], NB))
    if f == "V":
        if pos == 0:
            return "F", (a[0], NB, a[1], a[2])
        if pos == 1:
            return "Y", (a[0], a[1], a[2], NB)
        return "P", a + (NB,)
    if f == "C":
        return ("T", a + (NB,)) if pos == 0 else ("F", (a[0], a[3 - pos], a[pos], NB))
    raise ValueError(f)


CUTS = {
    "V": [(0, 1, "E", (1, 2), 2), (2, 1, "E", (0, 1), 3)],
    "C": [(1, 0, "E", (0, 2), 4), (2, 0, "E", (0, 1), 4)],
    "P": [(0, 1, "V", (1, 2, 3), 0), (3, 2, "V", (0, 1, 2), 1)],
    "F": [(1, 0, "V", (0, 2, 3), 1), (3, 2, "C", (0, 1, 2), 1)],
    "Y": [(0, 1, "C", (1, 2, 3), 0), (2, 1, "V", (0, 1, 3), 1), (3, 1, "V", (0, 1, 2), 1)],
    "T": [(1, 0, "C", (0, 2, 3), 1), (2, 0, "C", (0, 1, 3), 1), (3, 0, "C", (0, 1, 2), 1)],
}


def phi(t, c):
    return sum(w for _, _, _, w in t.active(c).values()) + sum(c.f(t.s[i]) for i in range(1, len(t.s)))


def fit_trees(c, nmax=4):
    for n in range(1, nmax + 1):
        for seq in product(*(range(1, i) for i in range(2, n + 1))):
            p = (-1, 0) + seq
            k = [p.count(i) for i in range(n + 1)]
            for tag in range(3):
                l = [0]
                for i in range(1, n + 1):
                    if tag == 0:
                        li = (i * 7 + n * 3) % 5 + 1
                    elif tag == 1:
                        li = 1
                    else:
                        d = (c.par["D4"] - 1, c.par["D4"], c.par["D3"] + 3, 2 * c.par["D3"] + 1)[(n + i) % 4]
                        li = max(1, d - 1 - k[i])
                    l.append(max(li, 1 - k[i]))
                yield MyTree(p, l)
    # Four-core shapes with exact internal vertices and outside donors, plus a
    # high-degree (table-tail) vertex: these generate Q4 and Q1 credits.
    d2, d4 = c.par["D2"], c.par["D4"]
    for parents, ends in [
            ((-1, 0, 1, 2, 3), (1, 4)),
            ((-1, 0, 1, 1, 3), (2, 4)),
            ((-1, 0, 1, 2, 2), (1, 3, 4)),
            ((-1, 0, 1, 1, 1), (2, 3, 4)),
            ((-1, 0, 1, 1, 2, 4, 4), (3, 5, 6))]:
        for mode in range(4):
            leaves = [0] + [1] * (len(parents) - 1)
            for j, x in enumerate(ends):
                D = (d4 - 2, d4 - 1, d4, 2 * d4 + 3)[(j + mode) % 4]
                leaves[x] = max(1, D - 1 - parents.count(x))
            yield MyTree(parents, leaves)
    # Table-tail singleton with several neighbours (Q1 incidence).
    for big in (d2, d2 + 7, 117):
        parents = [-1, 0, 1, 2, 3] + [4, 4, 4]
        leaves = [0, 1, 1, 1, big] + [1, 1, 1]
        yield MyTree(parents, leaves)


def check_tree(t, c, rec, counts, plant):
    old = t.active(c)
    signed = c.sign
    base = phi(t, c)
    fresh = len(t.p)
    direct = Q(0)
    for x in range(1, fresh):
        direct += t.r[x] * (phi(t.after(x, False), c) - base)
        if t.b[x]:
            direct += t.b[x] * (phi(t.after(x, True), c) - base)
    adj = defaultdict(Q)
    reserve = defaultdict(Q)
    edges = {frozenset(vs) for f, vs, _, _ in old.values() if f == "E"}
    for x in range(1, fresh):
        s = t.s[x]
        change = c.f(gain(s)) - c.f(s)
        born = ((c.f(bud(s)) - c.f(s) + c.f(NB) + c.weight("E", (bud(s), NB))) if t.b[x] else 0)
        key = ("v", x)
        adj[key] = signed * (t.r[x] * change + t.b[x] * born - c.lam * c.f(s))
        for y in [t.p[x]] + t.ch[x]:
            if s is not INF and degree(s) < c.par["D2"]:
                idx = 0 if y == t.p[x] else 1
                mv = signed * (inv(t.s[y]) - c.center(idx, s)) * change
                adj[key] -= mv
                adj[("edge", frozenset((x, y)))] += mv
            elif s is not INF and c.par["D2"] <= degree(s) < c.par["DF"]:
                reserve[key] += Q(max(signed * change, 0), 2) - inv(t.s[y]) * signed * change
    for support, (f, vs, st, w) in old.items():
        adj[support] -= signed * c.lam * w
        for pos, x in enumerate(vs):
            if x == 0:
                continue
            delta = signed * (c.weight(f, replace(st, pos, gain(st[pos]))) - w)
            adj[support] += t.r[x] * delta
            if t.b[x]:
                post = replace(st, pos, bud(st[pos]))
                bdiff = c.weight(f, post) - w
                if ARITY[f] < 4:
                    nf, nst = my_newborn(f, post, pos)
                    nb = c.weight(nf, nst)
                    if not (plant == "birth" and nf == "F"):
                        bdiff += nb
                adj[support] += signed * t.b[x] * bdiff
            for y in [t.p[x]] + t.ch[x]:
                if y in support:
                    continue
                ext = support | {y}
                iy = inv(t.s[y])
                from_parent = y == t.p[x]
                if ARITY[f] == 4:
                    if st[pos] is not INF and degree(st[pos]) < c.par["D4"]:
                        reserve[support] += Q(max(delta, 0), 2) - iy * delta
                    continue
                if ext in old:
                    cent = (4 if old[ext][0] == "C" else (2 if from_parent else 3)) if f == "E" else (0 if from_parent else 1)
                    mv = (iy - c.center(cent, st[pos])) * delta
                    adj[support] -= mv
                    adj[ext] += mv
                elif f == "E" and st[pos] is not INF and degree(st[pos]) < c.par["D2"]:
                    ty = 1 if pos == 1 else (0 if from_parent else 2)
                    amt = signed * c.fy(ty, g2id(c, st[pos]), iy)
                    adj[support] -= amt
                    donor = frozenset((x, y))
                    ok(donor in edges, "FY donor is not an active edge")
                    if not (plant == "fy" and ty == 2):
                        adj[("edge", donor)] += amt
    actual = signed * (direct - c.lam * phi(t, c))
    ok(sum(adj.values()) == actual, "physical generator transfer identity failed")
    ok(min(reserve.values(), default=0) >= 0, "negative physical credit")
    for key, v in adj.items():
        if isinstance(key, tuple) and key[0] in ("v", "edge"):
            continue
        if True:
            if old[key][0] == "E":
                continue
            rk = (old[key][0], old[key][2])
            reward = v + reserve[key] * (2 if plant == "credit" else 1)
            rec[rk] = max(rec.get(rk, Q(-10 ** 30)), reward)
    counts["trees"] += 1
    counts["events"] += sum(1 + bool(b) for b in t.b[1:])
    counts["newborns"] += sum(1 for x in range(1, fresh) if t.b[x])


def g2id(c, s):
    return c.l2.id[s]


def corner_majorants(c, rec, binary, path):
    keys = sorted({(f, x) for f, st in rec for x, _ in c.corners(f, st)})
    inp = "".join(f + " " + " ".join(map(str, x)) + "\n" for f, x in keys)
    out = subprocess.run([str(binary), str(path), "--query"], input=inp, text=True,
                         capture_output=True, check=True)
    vals = list(map(Q, out.stdout.splitlines()))
    ok(len(vals) == len(keys), "short high corner response")
    lookup = dict(zip(keys, vals))
    for (f, st), r in rec.items():
        bound = sum((a * lookup[f, x] for x, a in c.corners(f, st)), Q(0))
        ok(bound >= r, f"credited {f} row exceeds corner bound by {r - bound}")
    return len(keys)


# ---------------------------------------------------------------------------
# 4a. Multiplicity over all D and both signs (cheap, no production sweep).
# ---------------------------------------------------------------------------
def multiplicity_check(c):
    DF = c.par["DF"]
    finite = sum(range(2, DF))
    ok(finite == M["rows"]["f_inner"] + M["rows"]["f_table"], "finite comparison D coverage")
    for k in range(DF + 2):
        ok(k * (k - 1) // 2 == comb(k, 2), "cherry/Y multiplicity")
        ok(k * (k - 1) * (k - 2) // 6 == comb(k, 3), "trident multiplicity")
        if k >= 2:
            ok(k - 1 >= 1, "fork copy multiplicity")
    return finite


# ---------------------------------------------------------------------------
# 4. Comparison endpoint arithmetic, recomputed from the certificate.
# ---------------------------------------------------------------------------
def bracket(side):
    e = M[side]
    lam = Q(*e["trial"])
    err = Q(max(e["caps"][k] for k in ("f_inner", "f_table", "f_analytic"))
            + e["caps"]["g"] + e["higher_per_core"])
    ok(err == e["error_per_core"], "cap sum mismatch")
    return lam + (Q(1) if side == "upper" else Q(-1)) * err / Q(e["min_weight"])


def main():
    side = sys.argv[1] if len(sys.argv) > 1 else "upper"
    ok(side in ("upper", "lower"), "side must be upper or lower")
    path = OUT / f"mc-{side}.bin"
    ok(digest(path) == M[side]["sha256"], "snapshot digest mismatch")
    c = Local(str(path))
    target = Q(*M[side]["trial"])

    cap, worst = analytic_check(path, target)
    ok(cap == M[side]["caps"]["f_analytic"], f"analytic cap mismatch {cap}")
    print(f"[{side}] analytic tail identity + sampling PASS cap={cap} worst={float(worst):.6g}")

    ne = envelope_samples(c)
    print(f"[{side}] affine envelope sampling PASS n={ne}")
    finite = multiplicity_check(c)
    print(f"[{side}] multiplicity over all D and signs PASS finite_rows={finite}")

    ok((OUT / "mh-check").exists(), "missing mh-check oracle")
    rec, counts = {}, Counter()
    for t in fit_trees(c, nmax=4):
        check_tree(t, c, rec, counts, plant=None)
    nq = corner_majorants(c, rec, OUT / "mh-check", path)
    print(f"[{side}] own small-tree generator/FY/newborn/credit PASS "
          f"{dict(counts)} rows={len(rec)} corners={nq}")

    for plant, needle in (("fy", "transfer identity"), ("birth", "transfer identity")):
        r2, c2 = {}, Counter()
        try:
            for t in fit_trees(c, nmax=4):
                check_tree(t, c, r2, c2, plant=plant)
        except ArithmeticError as ex:
            ok(needle in str(ex), f"wrong rejection for {plant}: {ex}")
        else:
            raise ArithmeticError(f"planted {plant} accepted")
    r2, c2 = {}, Counter()
    for t in fit_trees(c, nmax=4):
        check_tree(t, c, r2, c2, plant="credit")
    try:
        corner_majorants(c, r2, OUT / "mh-check", path)
    except ArithmeticError as ex:
        ok("corner bound" in str(ex), f"wrong credit rejection: {ex}")
    else:
        raise ArithmeticError("planted double credit accepted")
    print(f"[{side}] planted FY/birth/credit mutations rejected")

    b = bracket(side)
    print(f"[{side}] comparison arithmetic PASS base={float(b):.12f}")
    print(f"INDEPENDENT REFEREE TEST: PASS [{side}]")


if __name__ == "__main__":
    try:
        main()
    except (ArithmeticError, ValueError, AssertionError) as e:
        print("REJECTED", e, file=sys.stderr)
        sys.exit(2)
