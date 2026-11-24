#!/usr/bin/env python3
"""
Complete exact rational certificate verifier, independent Python Fraction concrete-tree oracle,
and planted negative-control mutation suite for FT2_U5V1 unoriented 5-vertex certificates.

ZERO floating-point arithmetic is used in any verification or bound deduction.
"""
from fractions import Fraction as F
from pathlib import Path
import argparse
import hashlib
import json
import os
import random
import struct
import subprocess
import sys
import tempfile

sys.path.insert(0, str(Path(__file__).resolve().parent))
import verify_tail_bernstein as vtb


def sha256_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def sha256_file(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as fp:
        for chunk in iter(lambda: fp.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def combine_side(bin_path, verify_bin, threads=14, override_pq=None, perturb=None, expected_mode=None, return_raw=False, expected_sha256=None):
    raw = Path(bin_path).read_bytes()
    h_raw = sha256_bytes(raw)
    if expected_sha256 is not None and h_raw != expected_sha256:
        return {
            "ok": False,
            "returncode": 2,
            "stderr": f"input binding mismatch: snapshot sha256={h_raw} != expected frozen base {expected_sha256}",
        }
    cmd = [str(verify_bin), str(bin_path), "--threads", str(threads)]
    if override_pq is not None:
        cmd += ["--override-pq", str(override_pq[0]), str(override_pq[1])]
    if perturb is not None:
        cmd += ["--perturb", perturb]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        return {"ok": False, "returncode": proc.returncode, "stderr": proc.stderr.strip()}
    v = json.loads(proc.stdout)
    # Cryptographically bind C++ verifier's in-memory buffer to Python's exact private bytes!
    if v.get("sha256") != h_raw:
        return {
            "ok": False,
            "returncode": 2,
            "stderr": f"input binding mismatch: verifier sha256={v.get('sha256')} != python sha256={h_raw}",
        }
    target = F(override_pq[0], override_pq[1]) if override_pq is not None else None
    t = vtb.verify_bytes(raw, target=target, panels=16)

    S = v["S"]
    alpha = v["alpha_ticks"]
    mode = v["mode"]
    if mode != t["mode"] or mode not in (1, 2):
        return {"ok": False, "returncode": 2, "stderr": f"mode mismatch: c++={mode}, py={t['mode']}"}
    if expected_mode is not None and mode != expected_mode:
        return {"ok": False, "returncode": 2, "stderr": f"expected mode {expected_mode}, got {mode}"}

    # Global residual cap per core vertex (valid for EVERY finite tree T with C(T) >= 1!):
    #   C = 1: residual <= cap_f <= cap_f + max(0, cap_g) + higher_per_core
    #   C >= 2: residual <= cap_f*C + cap_g*(C-1) + higher_per_core*C <= (cap_f + max(0, cap_g) + higher_per_core)*C
    cap_f = max(v["caps"]["f_in"], v["caps"]["f_tab"], t["cap_ticks"])
    cap_g_pos = max(0, v["caps"]["g"])
    higher = v["higher_per_core"]
    err_per_core = max(0, cap_f + cap_g_pos + higher)

    # Global positive comparability: 0 < m_min * C(T) <= S * Phi(T) <= m_max * C(T) for ALL C(T) >= 1!
    cf = v["comp_fin"]
    # Analytic tail (D >= DF): zero higher motifs, only f_ana + 1 edge g in [g_far_min, g_far_max]
    min_c0_ana = t["w_c0_min_ticks"]
    min_c1_ana = t["w_c1_min_ticks"] + cf["g_far_min"] - alpha
    min_c2_ana = t["w_c2_min_ticks"] + cf["g_far_min"]
    min_c3_ana = t["w_min_ticks"] + cf["g_far_min"] + alpha  # c >= 3 => alpha*(c-2) >= alpha
    m_min = min(cf["min_W_c0"], cf["min_W"], min_c0_ana, min_c1_ana, min_c2_ana, min_c3_ana)
    m_max = max(cf["max_W_c0"], cf["max_W"], t["w_max_ticks"] + cf["g_far_max"])
    if m_min <= 0 or m_max <= 0:
        return {"ok": False, "returncode": 2, "stderr": f"comparability failed: m_min={m_min}, m_max={m_max}"}

    lam_nom = F(v["P"], v["Q"])
    adj = F(err_per_core, m_min)
    certified = (lam_nom + adj) if mode == 1 else (lam_nom - adj)
    res = {
        "ok": True,
        "sha256": h_raw,
        "mode": mode,
        "nominal": str(lam_nom),
        "certified_exact": str(certified),
        "certified_dec": format_dec(certified, 12, round_up=(mode == 1)),
        "err_per_core": err_per_core,
        "cap_f": cap_f,
        "cap_f_ana": t["cap_ticks"],
        "cap_g": v["caps"]["g"],
        "higher_per_core": higher,
        "m_min": m_min,
        "m_max": m_max,
        "verify": v,
        "tail": t,
    }
    if return_raw:
        res["_raw"] = raw
    return res


def format_dec(q: F, digits=12, round_up=True):
    scale = 10 ** digits
    v = q * scale
    if round_up:
        n = -((-v.numerator) // v.denominator)
    else:
        n = v.numerator // v.denominator
    neg = n < 0
    n = abs(n)
    return f"{'-' if neg else ''}{n // scale}.{n % scale:0{digits}d}"


# === Independent Python Fraction Snapshot Reader & Concrete Unoriented Tree Oracle ===
class PySnapshot:
    """Independently parses FT2_U5V1 binary and evaluates S*Phi(T) in exact Fraction arithmetic."""

    def __init__(self, path=None, raw=None):
        if raw is None:
            raw = Path(path).read_bytes()
        self.h = vtb.read_header_bytes(raw)
        h = self.h
        self.S = h["S"]
        self.mode = h["mode"]
        self.lam = F(h["P"], h["Q"])

        # Build EndSpaces independently
        def build_end(dc, Kc):
            idx = {}
            ovf = {}
            nxt = 0
            for d in range(2, dc):
                for c in range(1, min(d, Kc) + 1):
                    idx[(d, c)] = nxt
                    nxt += 1
            for d in range(Kc + 1, dc):
                ovf[d] = nxt
                nxt += 1
            tail_base = nxt
            nxt += Kc + 1
            inf_id = nxt
            N = nxt + 1
            return {"dc": dc, "Kc": Kc, "idx": idx, "ovf": ovf, "tb": tail_base, "inf": inf_id, "N": N, "nin_ex": len(idx)}

        def build_inner(dc, mnC, mxC):
            idx = {}
            lst = []
            for d in range(max(2, mnC), dc):
                for c in range(mnC, min(d, mxC) + 1):
                    idx[(d, c)] = len(lst)
                    lst.append((d, c))
            return {"idx": idx, "lst": lst, "N": len(lst)}

        self.G2 = build_end(h["d2"], h["K2"])
        self.G3 = build_end(h["d3e"], h["K3e"])
        self.G4 = build_end(h["d4e"], h["K4e"])
        self.G5 = build_end(h["d5e"], h["K5e"])

        self.M3 = build_inner(h["d3m"], 2, h["K3m"])
        self.M4v = build_inner(h["d4i"], 2, h["K4i"])
        self.MS3 = build_inner(h["dS3"], 3, h["KS3"])
        self.I5_2 = build_inner(h["d5i"], 2, h["K5i"])
        self.I5_3 = build_inner(h["d5o"], 3, h["K5o"])
        self.I5_2c = build_inner(h["d5c"], 2, h["K5c"])
        self.MS4 = build_inner(h["dS4"], 4, h["KS4"])

        # P4 canonical blocks
        self.p4_blk = {}
        n_p4 = 0
        for b, (db, cb) in enumerate(self.M4v["lst"]):
            for c, (dc_, cc) in enumerate(self.M4v["lst"]):
                if b <= c and db + dc_ <= h["SP4"]:
                    self.p4_blk[(b, c)] = n_p4
                    n_p4 += (self.G4["N"] * (self.G4["N"] + 1) // 2) if b == c else (self.G4["N"] ** 2)

        # P5 canonical blocks
        self.p5_blk = {}
        n_p5 = 0
        for b, (db, cb) in enumerate(self.I5_2["lst"]):
            for c, (dc_, cc) in enumerate(self.I5_2["lst"]):
                if db + dc_ <= h["SP4"]:
                    for d, (dd, cd) in enumerate(self.I5_2["lst"]):
                        if b <= d and dc_ + dd <= h["SP4"] and db + dc_ + dd <= h["SP5"]:
                            self.p5_blk[(b, c, d)] = n_p5
                            n_p5 += (self.G5["N"] * (self.G5["N"] + 1) // 2) if b == d else (self.G5["N"] ** 2)

        # B5 pairs
        self.b5_map = {}
        nb5 = 0
        for o, (do, co) in enumerate(self.I5_3["lst"]):
            for c, (dc_, cc) in enumerate(self.I5_2c["lst"]):
                if do + dc_ <= h["SB5"]:
                    self.b5_map[(o, c)] = nb5
                    nb5 += 1

        # Triples in G4, quads in G5
        def tri_map(N):
            mp = {}
            k = 0
            for a in range(N):
                for b in range(a, N):
                    for c in range(b, N):
                        mp[(a, b, c)] = k
                        k += 1
            return mp, k

        def quad_map(N):
            mp = {}
            k = 0
            for a in range(N):
                for b in range(a, N):
                    for c in range(b, N):
                        for d in range(c, N):
                            mp[(a, b, c, d)] = k
                            k += 1
            return mp, k

        self.tri4, self.n_tri4 = tri_map(self.G4["N"])
        self.quad5, self.n_quad5 = quad_map(self.G5["N"])

        n_r1 = self.G2["nin_ex"]
        n_knots = n_r1 * h["nknots"]
        n_f = (h["DF"] - 2) * (h["DF"] + 3) // 2
        n_g = self.G2["N"] * (self.G2["N"] + 1) // 2
        n_h3 = self.M3["N"] * (self.G3["N"] * (self.G3["N"] + 1) // 2)
        n_s3 = self.MS3["N"] * self.n_tri4
        n_b5 = nb5 * (self.G5["N"] * (self.G5["N"] + 1) // 2) * self.G5["N"]
        n_s4 = self.MS4["N"] * self.n_quad5

        # Memory-slice payload int64s cleanly without unpacking huge arrays!
        self.buf = memoryview(raw)[512:].cast("q")
        cur = 0
        def take(n):
            nonlocal cur
            o = cur
            cur += n
            return o
        self.off_r1 = take(n_r1)
        self.off_knots = take(n_knots)
        self.off_f = take(n_f)
        self.off_g = take(n_g)
        self.off_h3 = take(n_h3)
        self.off_p4 = take(n_p4)
        self.off_s3 = take(n_s3)
        self.off_p5 = take(n_p5)
        self.off_b5 = take(n_b5)
        self.off_s4 = take(n_s4)

    @staticmethod
    def sym(a, b):
        if a > b:
            a, b = b, a
        return b * (b + 1) // 2 + a

    @staticmethod
    def wts(E, D, c):
        if D < E["dc"]:
            return [(E["idx"][(D, c)] if c <= E["Kc"] else E["ovf"][D], F(1))]
        t = F(E["dc"], D)
        return [(E["tb"] + min(c, E["Kc"] + 1) - 1, t), (E["inf"], 1 - t)]

    def eval_f(self, D, c):
        h = self.h
        if D < h["DF"]:
            return F(self.buf[self.off_f + (D - 2) * (D + 3) // 2 + c])
        q = F(D - c, D)
        if self.mode == 1:
            return h["f_p0"] + h["f_p1"] * q
        return h["f_p0"] * q - h["f_p1"] * q * q - h["f_p2"] * q * F(1, D)

    def eval_g(self, Du, cu, Dv, cv):
        res = F(0)
        for u, wu in self.wts(self.G2, Du, cu):
            for v, wv in self.wts(self.G2, Dv, cv):
                res += wu * wv * self.buf[self.off_g + self.sym(u, v)]
        return res

    def eval_h3(self, Da, ca, Dm, cm, Db, cb):
        m = self.M3["idx"].get((Dm, cm))
        if m is None:
            return F(0)
        p2 = self.G3["N"] * (self.G3["N"] + 1) // 2
        res = F(0)
        for a, wa in self.wts(self.G3, Da, ca):
            for b, wb in self.wts(self.G3, Db, cb):
                res += wa * wb * self.buf[self.off_h3 + m * p2 + self.sym(a, b)]
        return res

    def eval_p4(self, Da, ca, Db, cb, Dc, cc, Dd, cd):
        b = self.M4v["idx"].get((Db, cb))
        c = self.M4v["idx"].get((Dc, cc))
        if b is None or c is None or Db + Dc > self.h["SP4"]:
            return F(0)
        if b > c:
            b, c = c, b
            Da, ca, Dd, cd = Dd, cd, Da, ca
        base = self.p4_blk[(b, c)]
        N4 = self.G4["N"]
        res = F(0)
        for a, wa in self.wts(self.G4, Da, ca):
            for d, wd in self.wts(self.G4, Dd, cd):
                off = self.sym(a, d) if b == c else (a * N4 + d)
                res += wa * wd * self.buf[self.off_p4 + base + off]
        return res

    def eval_s3(self, Do, co, Da, ca, Db, cb, Dc, cc):
        o = self.MS3["idx"].get((Do, co))
        if o is None:
            return F(0)
        res = F(0)
        for a, wa in self.wts(self.G4, Da, ca):
            for b, wb in self.wts(self.G4, Db, cb):
                for c, wc in self.wts(self.G4, Dc, cc):
                    t = tuple(sorted((a, b, c)))
                    res += wa * wb * wc * self.buf[self.off_s3 + o * self.n_tri4 + self.tri4[t]]
        return res

    def eval_p5(self, Da, ca, Db, cb, Dc, cc, Dd, cd, De, ce):
        b = self.I5_2["idx"].get((Db, cb))
        c = self.I5_2["idx"].get((Dc, cc))
        d = self.I5_2["idx"].get((Dd, cd))
        if b is None or c is None or d is None:
            return F(0)
        if Db + Dc > self.h["SP4"] or Dc + Dd > self.h["SP4"] or Db + Dc + Dd > self.h["SP5"]:
            return F(0)
        if b > d:
            b, d = d, b
            Da, ca, De, ce = De, ce, Da, ca
        base = self.p5_blk[(b, c, d)]
        N5 = self.G5["N"]
        res = F(0)
        for a, wa in self.wts(self.G5, Da, ca):
            for e, we in self.wts(self.G5, De, ce):
                off = self.sym(a, e) if b == d else (a * N5 + e)
                res += wa * we * self.buf[self.off_p5 + base + off]
        return res

    def eval_b5(self, Do, co, Da, ca, Db, cb, Dc, cc, Dd, cd):
        o = self.I5_3["idx"].get((Do, co))
        c = self.I5_2c["idx"].get((Dc, cc))
        if o is None or c is None or Do + Dc > self.h["SB5"]:
            return F(0)
        bid = self.b5_map[(o, c)]
        N5 = self.G5["N"]
        p2 = N5 * (N5 + 1) // 2
        res = F(0)
        for a, wa in self.wts(self.G5, Da, ca):
            for b, wb in self.wts(self.G5, Db, cb):
                sab = self.sym(a, b)
                for d, wd in self.wts(self.G5, Dd, cd):
                    res += wa * wb * wd * self.buf[self.off_b5 + (bid * p2 + sab) * N5 + d]
        return res

    def eval_s4(self, Do, co, Da, ca, Db, cb, Dc, cc, Dd, cd):
        o = self.MS4["idx"].get((Do, co))
        if o is None:
            return F(0)
        res = F(0)
        for a, wa in self.wts(self.G5, Da, ca):
            for b, wb in self.wts(self.G5, Db, cb):
                for c, wc in self.wts(self.G5, Dc, cc):
                    for d, wd in self.wts(self.G5, Dd, cd):
                        q = tuple(sorted((a, b, c, d)))
                        res += wa * wb * wc * wd * self.buf[self.off_s4 + o * self.n_quad5 + self.quad5[q]]
        return res

    def eval_phi(self, cadj, D):
        n = len(cadj)
        c = [len(cadj[u]) for u in range(n)]
        tot = F(0)
        for u in range(n):
            tot += self.eval_f(D[u], c[u])
        for u in range(n):
            for v in cadj[u]:
                if u < v:
                    tot += self.eval_g(D[u], c[u], D[v], c[v])
        for m in range(n):
            nb = cadj[m]
            for i in range(len(nb)):
                for j in range(i + 1, len(nb)):
                    tot += self.eval_h3(D[nb[i]], c[nb[i]], D[m], c[m], D[nb[j]], c[nb[j]])
        for b in range(n):
            for cc in cadj[b]:
                if b < cc:
                    for a in cadj[b]:
                        if a != cc:
                            for d in cadj[cc]:
                                if d != b:
                                    tot += self.eval_p4(D[a], c[a], D[b], c[b], D[cc], c[cc], D[d], c[d])
        for o in range(n):
            nb = cadj[o]
            k = len(nb)
            for i in range(k):
                for j in range(i + 1, k):
                    for l in range(j + 1, k):
                        tot += self.eval_s3(D[o], c[o], D[nb[i]], c[nb[i]], D[nb[j]], c[nb[j]], D[nb[l]], c[nb[l]])
        for m in range(n):
            nb = cadj[m]
            for i in range(len(nb)):
                for j in range(i + 1, len(nb)):
                    b, d = nb[i], nb[j]
                    for a in cadj[b]:
                        if a != m:
                            for e in cadj[d]:
                                if e != m:
                                    tot += self.eval_p5(D[a], c[a], D[b], c[b], D[m], c[m], D[d], c[d], D[e], c[e])
        for o in range(n):
            nb = cadj[o]
            k = len(nb)
            for ic in range(k):
                cc = nb[ic]
                for i in range(k):
                    if i == ic:
                        continue
                    for j in range(i + 1, k):
                        if j == ic:
                            continue
                        for d in cadj[cc]:
                            if d != o:
                                tot += self.eval_b5(D[o], c[o], D[nb[i]], c[nb[i]], D[nb[j]], c[nb[j]], D[cc], c[cc], D[d], c[d])
        for o in range(n):
            nb = cadj[o]
            k = len(nb)
            for i in range(k):
                for j in range(i + 1, k):
                    for l in range(j + 1, k):
                        for r in range(l + 1, k):
                            tot += self.eval_s4(D[o], c[o], D[nb[i]], c[nb[i]], D[nb[j]], c[nb[j]], D[nb[l]], c[nb[l]], D[nb[r]], c[nb[r]])
        return tot

    def exact_tree_check(self, cadj, D):
        n = len(cadj)
        phi0 = self.eval_phi(cadj, D)
        Lphi = F(0)
        for u in range(n):
            Lv = D[u] - len(cadj[u])
            gr = F(Lv) + sum(F(1, D[v]) for v in cadj[u])
            if gr > 0:
                D2 = list(D)
                D2[u] += 1
                Lphi += gr * (self.eval_phi(cadj, D2) - phi0)
            if Lv > 0:
                br = F(Lv, D[u])
                cadj2 = [list(x) for x in cadj] + [[u]]
                cadj2[u].append(n)
                D2 = list(D) + [2]
                Lphi += br * (self.eval_phi(cadj2, D2) - phi0)
        sgn = 1 if self.mode == 1 else -1
        return phi0, sgn * (Lphi - self.lam * phi0)


def run_tree_oracle(bin_path, side_res, raw=None):
    snap = PySnapshot(bin_path, raw=raw)
    DF, d2 = snap.h["DF"], snap.h["d2"]
    m_min, m_max, err = side_res["m_min"], side_res["m_max"], side_res["err_per_core"]
    trees = []
    # 1. Single-core stars C=1 (including D=DF-1, D=DF, D=DF+17!)
    for D0 in (2, 3, 5, d2, DF - 1, DF, DF + 17):
        trees.append(("star_C1_D" + str(D0), [[]], [D0]))
    # 2. Bare chains & tail-end chains
    for ln in (2, 3, 4, 5, 6, 7):
        adj = [[] for _ in range(ln)]
        for i in range(ln - 1):
            adj[i].append(i + 1)
            adj[i + 1].append(i)
        trees.append((f"chain_{ln}", adj, [2] * ln))
        D2 = [2] * ln
        D2[0] = d2 + 5
        D2[-1] = DF + 11
        trees.append((f"chain_tail_{ln}", adj, D2))
    # 3. High-arity stars (including bare kind=3 overflow stars + brooms!)
    for arms in (3, 4, 5, 8):
        adj = [[i for i in range(1, arms + 1)]] + [[0] for _ in range(arms)]
        trees.append((f"bare_star_{arms}", adj, [arms] + [2] * arms))
        trees.append((f"leaf_star_{arms}", adj, [arms + 1] + [2] * arms))
        adj_b = [list(x) for x in adj] + [[1]]
        adj_b[1].append(arms + 1)
        trees.append((f"broom_{arms}", adj_b, [arms + 2, 3] + [2] * (arms - 1) + [d2 + 3]))
    # 4. Double-star (two adjacent high-arity hubs)
    a1, a2 = 5, 5
    adj = [[] for _ in range(2 + a1 + a2)]
    adj[0].append(1); adj[1].append(0)
    for i in range(a1):
        adj[0].append(2 + i); adj[2 + i].append(0)
    for i in range(a2):
        adj[1].append(2 + a1 + i); adj[2 + a1 + i].append(1)
    trees.append(("double_star_ovf", adj, [a1 + 1, a2 + 2] + [2] * (a1 + a2)))
    # 5. 15 random friend trees of sizes 6..28
    rng = random.Random(20261124)
    for tr in range(15):
        n = rng.randint(6, 28)
        full = [[1, 2], [0], [0]] + [[] for _ in range(n - 3)]
        for k in range(3, n):
            v = rng.randrange(k)
            w = full[v][rng.randrange(len(full[v]))]
            full[w].append(k); full[k].append(w)
        cid = {}
        for i in range(n):
            if len(full[i]) >= 2:
                cid[i] = len(cid)
        cadj = [[] for _ in range(len(cid))]
        D = [0] * len(cid)
        for i, u in cid.items():
            D[u] = len(full[i])
            for j in full[i]:
                if j in cid:
                    cadj[u].append(cid[j])
        trees.append((f"rand_{tr}", cadj, D))

    worst_ratio = None
    for name, cadj, D in trees:
        C = len(cadj)
        phi0, res = snap.exact_tree_check(cadj, D)
        if not (m_min * C <= phi0 <= m_max * C):
            raise AssertionError(f"Comparability failed on {name}: phi0={phi0}, bounds=[{m_min*C}, {m_max*C}]")
        if res > err * C:
            raise AssertionError(f"Global drift inequality FAILED on {name}: res={res} > err*C={err*C}")
        r = res / C
        if worst_ratio is None or r > worst_ratio:
            worst_ratio = r
    return {"trees_checked": len(trees), "worst_res_per_core": str(worst_ratio), "err_bound": err}


def run_mutations(bin_path, verify_bin, base_res, expected_sha256):
    results = {}
    # 1. Hard fail-closed structural & envelope mutations (MUST fail with returncode == 2 and specific diagnostic stderr!)
    expected_diag = {
        "ginf": ("g[INF,INF]",),
        "env23": ("env23", "envelope radius negative"),
        "env34p": ("env34p", "envelope radius negative"),
        "env34s": ("env34s", "envelope radius negative"),
        "env45p": ("env45p", "envelope radius negative"),
        "env45bp": ("env45bp", "envelope radius negative"),
        "env45bs": ("env45bs", "envelope radius negative"),
        "env45s4": ("env45s4", "envelope radius negative"),
    }
    for tgt, kws in expected_diag.items():
        r = combine_side(bin_path, verify_bin, threads=1, perturb=tgt, expected_sha256=expected_sha256)
        if r["ok"] or r.get("returncode") != 2 or not r.get("stderr"):
            raise AssertionError(f"Planted mutation {tgt} was NOT rejected cleanly (rc=2): {r}")
        if not any(kw in r["stderr"] for kw in kws):
            raise AssertionError(f"Planted mutation {tgt} unexpected stderr: {r['stderr']}")
        results[tgt] = f"REJECTED (rc={r['returncode']}: {r['stderr']})"

    # 2. Weight/center/knot perturbations (MUST reject or drastically blow up err_per_core on the SAME frozen bytes!)
    base_err = base_res["err_per_core"]
    for tgt in ("f", "g", "h3", "p4", "s3", "p5", "b5", "s4", "r1", "knots"):
        r = combine_side(bin_path, verify_bin, threads=1, perturb=tgt, expected_sha256=expected_sha256)
        if r["ok"]:
            if r["verify"]["sha256"] != expected_sha256:
                raise AssertionError(f"Planted perturbation {tgt} ran on wrong snapshot digest!")
            if r["err_per_core"] <= base_err + 100000:
                raise AssertionError(f"Planted perturbation {tgt} did not increase err_per_core: {r['err_per_core']} vs {base_err}")
        results[tgt] = f"CAUGHT (err_per_core={r.get('err_per_core', 'REJECTED')} vs base {base_err})"

    # 3. Tightening target exponent P/Q by 1e-4 in the forbidden direction MUST be caught!
    P, Q = base_res["verify"]["P"], base_res["verify"]["Q"]
    bad = F(P, Q) + (F(-1, 10000) if base_res["mode"] == 1 else F(1, 10000))
    r_bad = combine_side(bin_path, verify_bin, threads=1, override_pq=(bad.numerator, bad.denominator), expected_sha256=expected_sha256)
    if r_bad["ok"]:
        if r_bad["verify"]["sha256"] != expected_sha256:
            raise AssertionError("Tighter P/Q probe ran on wrong snapshot digest!")
        if r_bad["err_per_core"] <= base_err + 10**7:
            raise AssertionError("Tighter P/Q mutation was NOT caught!")
    results["tighter_pq"] = f"CAUGHT (err={r_bad.get('err_per_core', 'REJECTED')}, cert={r_bad.get('certified_dec', 'N/A')})"
    return results


def run_wrapper_and_header_controls(frozen_lo, frozen_hi, raw_lo, lo_sha, verify_bin):
    """Fail-closed regression tests for wrapper mode binding, digest binding & malformed binary headers."""
    results = {}
    self_script = str(Path(__file__).resolve())
    for label, a_lo, a_hi in (
        ("dup_upper", frozen_hi, frozen_hi),
        ("dup_lower", frozen_lo, frozen_lo),
        ("swapped_sides", frozen_hi, frozen_lo),
    ):
        proc = subprocess.run(
            [sys.executable, "-B", self_script, "--verify-bin", str(verify_bin),
             "--lower", str(a_lo), "--upper", str(a_hi), "--threads", "1", "--skip-trees"],
            capture_output=True, text=True,
        )
        if proc.returncode != 2:
            raise AssertionError(f"Wrapper control {label} was NOT rejected with rc=2 (got rc={proc.returncode})")
        results[label] = f"REJECTED (rc={proc.returncode}: {proc.stderr.strip()})"

    # One-field header attacks derived strictly from frozen raw_lo bytes:
    # (a) f_p0 = INT64_MIN (offset 56): must fail header validation (code 30), never silently narrow!
    # (b) knot_num[1] = 630 + 2^60 (offset 312): must fail range validation before product (code 28)!
    # (c) P = -7 (offset 16) even WITH --override-pq 1 2: original header MUST fail validation (code 3)!
    for label, off, val, extra_args in (
        ("header_fp0_int64min", 56, -(1 << 63), []),
        ("header_knot1_overflow", 312, 630 + (1 << 60), []),
        ("header_badP_override", 16, -7, ["--override-pq", "1", "2"]),
    ):
        mut = bytearray(raw_lo)
        struct.pack_into("<q", mut, off, val)
        fd, tmp_name = tempfile.mkstemp(prefix=f"q5_{label}_", suffix=".bin")
        try:
            os.write(fd, mut)
            os.close(fd)
            proc = subprocess.run([str(verify_bin), tmp_name, "--threads", "1", *extra_args], capture_output=True, text=True)
            if proc.returncode != 2:
                raise AssertionError(f"Header control {label} was NOT rejected with rc=2 (got rc={proc.returncode})")
            results[label] = f"REJECTED (rc={proc.returncode}: {proc.stderr.strip()})"
        finally:
            if os.path.exists(tmp_name):
                os.unlink(tmp_name)

    # Strict CLI integer parsing: leading newline must be rejected!
    proc = subprocess.run([str(verify_bin), str(frozen_lo), "--threads", "\n1"], capture_output=True, text=True)
    if proc.returncode == 0:
        raise AssertionError("CLI control bad_threads_newline was NOT rejected")
    results["cli_threads_newline"] = f"REJECTED (rc={proc.returncode}: {proc.stderr.strip()})"

    # Frozen-digest binding control: presenting a different snapshot under lo_sha MUST be rejected!
    r_mis = combine_side(frozen_hi, verify_bin, threads=1, expected_sha256=lo_sha)
    if r_mis["ok"] or r_mis.get("returncode") != 2 or "input binding mismatch" not in r_mis.get("stderr", ""):
        raise AssertionError(f"Digest-binding control was NOT rejected cleanly: {r_mis}")
    results["digest_binding_mismatch"] = f"REJECTED (rc={r_mis['returncode']}: {r_mis['stderr']})"
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-bin", required=True)
    ap.add_argument("--lower", required=True)
    ap.add_argument("--upper", required=True)
    ap.add_argument("--threads", type=int, default=14)
    ap.add_argument("--run-mutations", action="store_true")
    ap.add_argument("--skip-trees", action="store_true")
    args = ap.parse_args()

    need_raw = (not args.skip_trees) or args.run_mutations
    lo = combine_side(args.lower, args.verify_bin, threads=args.threads, expected_mode=2, return_raw=need_raw)
    hi = combine_side(args.upper, args.verify_bin, threads=args.threads, expected_mode=1, return_raw=need_raw)
    if not lo["ok"] or not hi["ok"]:
        lo.pop("_raw", None)
        hi.pop("_raw", None)
        print(json.dumps({"lower": lo, "upper": hi}, indent=2), file=sys.stderr)
        sys.exit(2)

    width = F(hi["certified_exact"]) - F(lo["certified_exact"])
    if width < 0:
        lo.pop("_raw", None)
        hi.pop("_raw", None)
        print(json.dumps({"error": f"inverted certified interval: {lo['certified_exact']} > {hi['certified_exact']}"}), file=sys.stderr)
        sys.exit(2)

    raw_lo = lo.pop("_raw", None)
    raw_hi = hi.pop("_raw", None)
    lo_sha = lo.pop("sha256")
    hi_sha = hi.pop("sha256")
    out = {
        "lower_sha256": lo_sha,
        "upper_sha256": hi_sha,
        "certified_interval_dec": [lo["certified_dec"], hi["certified_dec"]],
        "certified_width_ppm": format_dec(width * 10**6, 4, round_up=True),
        "certifies_mu_gt_0_5657": (F(lo["certified_exact"]) > F(5657, 10000)),
        "lower": lo,
        "upper": hi,
    }
    if not args.skip_trees:
        out["tree_oracle_lower"] = run_tree_oracle(args.lower, lo, raw=raw_lo)
        out["tree_oracle_upper"] = run_tree_oracle(args.upper, hi, raw=raw_hi)
    if args.run_mutations:
        # Execute ALL optional controls strictly on private frozen copies of the verified base bytes
        # (never reopening the external --lower/--upper pathnames after base verification)!
        with tempfile.TemporaryDirectory(prefix="q5_frozen_") as td:
            frozen_lo = Path(td) / "lower.bin"
            frozen_hi = Path(td) / "upper.bin"
            frozen_lo.write_bytes(raw_lo)
            frozen_hi.write_bytes(raw_hi)
            out["mutations_lower"] = run_mutations(frozen_lo, args.verify_bin, lo, lo_sha)
            out["mutations_upper"] = run_mutations(frozen_hi, args.verify_bin, hi, hi_sha)
            out["wrapper_header_controls"] = run_wrapper_and_header_controls(frozen_lo, frozen_hi, raw_lo, lo_sha, args.verify_bin)
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
