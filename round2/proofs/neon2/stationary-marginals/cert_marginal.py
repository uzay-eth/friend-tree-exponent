#!/usr/bin/env python3
"""
Exact zero-floating-point rational verifier, independent concrete-tree Fraction oracle,
and planted negative-control suite for FT2_M5V1 stationary-marginal certificates.

Proves certified two-sided rational enclosures for stationary marginals pi(psi) of the
friend-tree eigenmeasure pi, conditional ONLY on a certified mu enclosure [mu_lo, mu_hi]
(taken from Neon-2 Claim 6 gold2: [565713486333/10^12, 565736986679/10^12]).
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
import verify_marginal_tail as vmt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "quintet-cert"))
import verify_tail_bernstein as vtb


def sha256_bytes(raw):
    return hashlib.sha256(raw).hexdigest()


def format_dec(q: F, digits=8, round_up=True):
    scale = 10 ** digits
    v = q * scale
    n = (-((-v.numerator) // v.denominator)) if round_up else (v.numerator // v.denominator)
    neg = n < 0
    n = abs(n)
    return f"{'-' if neg else ''}{n // scale}.{n % scale:0{digits}d}"


def eval_obs(obs_id, D, c, DF):
    if D >= DF:
        return 0
    if obs_id == 1:
        return 1 if D == 2 else 0
    if obs_id == 2:
        return 1 if D == 3 else 0
    if obs_id == 3:
        return 1 if D == 4 else 0
    if obs_id == 4:
        return 1 if D == 5 else 0
    if obs_id == 5:
        return 1 if D == 6 else 0
    if obs_id == 6:
        return 1 if (D == 2 and c == 1) else 0
    if obs_id == 7:
        return 1 if (D == 2 and c == 2) else 0
    if obs_id == 8:
        return 1 if D == c else 0
    raise ValueError(f"unknown obs_id={obs_id}")


def bare_core_tail_upper(DF: int, mu_lo: F) -> F:
    """
    Exact analytic upper bound on pi(L_o = 0, D_o >= DF) via the 4-line Lyapunov function
    h(D, c) = 0 (D <= 5), 2^(-(DF-D)) * 4^(-(D-c)) (6 <= D <= DF), 4^(-(D-c)) (D >= DF):
      (L H - mu H)(T) + mu * #{D_u >= DF, L_u = 0} <= (5/8) * 2^(-(DF - 6)) * C(T) = 40 * 2^(-DF) * C(T)
    => pi(L_o = 0, D_o >= DF) <= (40 / mu_lo) * 2^(-DF).
    """
    if DF < 6 or mu_lo <= 0:
        raise ValueError("DF must be >= 6 and mu_lo > 0")
    return F(40, 1) / (mu_lo * (1 << DF))


def combine_side(bin_path, verify_bin, threads=8, override_theta=None,
                 expected_obs=None, expected_sign=None, expected_mu_lo=None, expected_mu_hi=None,
                 return_raw=False, expected_sha256=None):
    if expected_mu_lo is not None and not (F(0) < expected_mu_lo < F(1)):
        return {"ok": False, "returncode": 2, "stderr": f"invalid expected_mu_lo={expected_mu_lo}"}
    if expected_mu_hi is not None and not (F(0) < expected_mu_hi < F(1)):
        return {"ok": False, "returncode": 2, "stderr": f"invalid expected_mu_hi={expected_mu_hi}"}
    if expected_mu_lo is not None and expected_mu_hi is not None and expected_mu_lo > expected_mu_hi:
        return {"ok": False, "returncode": 2, "stderr": f"reversed expected mu bracket: {expected_mu_lo} > {expected_mu_hi}"}
    raw = Path(bin_path).read_bytes()
    h_raw = sha256_bytes(raw)
    if expected_sha256 is not None and h_raw != expected_sha256:
        return {"ok": False, "returncode": 2, "stderr": f"input binding mismatch: {h_raw} != {expected_sha256}"}
    cmd = [str(verify_bin), str(bin_path), "--threads", str(threads)]
    if override_theta is not None:
        cmd += ["--override-theta", str(override_theta[0]), str(override_theta[1])]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        return {"ok": False, "returncode": proc.returncode, "stderr": proc.stderr.strip()}
    v = json.loads(proc.stdout)
    if v.get("perturbed"):
        return {"ok": False, "returncode": 2, "stderr": "refusing to certify perturbed verifier output"}
    if v.get("sha256") != h_raw:
        return {"ok": False, "returncode": 2, "stderr": "SHA-256 mismatch between C++ verifier and Python"}

    obs_id = v["obs_id"]
    obs_sign = v["obs_sign"]
    if expected_obs is not None and obs_id != expected_obs:
        return {"ok": False, "returncode": 2, "stderr": f"expected obs_id {expected_obs}, got {obs_id}"}
    if expected_sign is not None and obs_sign != expected_sign:
        return {"ok": False, "returncode": 2, "stderr": f"expected obs_sign {expected_sign}, got {obs_sign}"}

    mu_lo = F(v["mu_lo_P"], v["mu_lo_Q"])
    mu_hi = F(v["mu_hi_P"], v["mu_hi_Q"])
    if not (F(0) < mu_lo <= mu_hi < F(1)):
        return {"ok": False, "returncode": 2, "stderr": f"invalid snapshot mu bracket [{mu_lo}, {mu_hi}]"}
    if expected_mu_lo is not None and mu_lo > expected_mu_lo:
        return {"ok": False, "returncode": 2, "stderr": f"snapshot mu_lo={mu_lo} > certified mu_lo={expected_mu_lo}"}
    if expected_mu_hi is not None and mu_hi < expected_mu_hi:
        return {"ok": False, "returncode": 2, "stderr": f"snapshot mu_hi={mu_hi} < certified mu_hi={expected_mu_hi}"}

    th_ov = F(override_theta[0], override_theta[1]) if override_theta is not None else None
    t = vmt.verify_bytes(raw, mu_lo=mu_lo, mu_hi=mu_hi, theta_override=th_ov, panels=16)

    S = v["S"]
    # Exact per-core error at BOTH endpoints mu_lo and mu_hi (by affinity in lambda, max covers ALL lam in [mu_lo, mu_hi]!)
    cap_f_lo = max(v["caps_lo"]["f_in"], v["caps_lo"]["f_tab"], t["cap_ticks_lo"])
    cap_f_hi = max(v["caps_hi"]["f_in"], v["caps_hi"]["f_tab"], t["cap_ticks_hi"])
    err_lo = max(0, cap_f_lo + max(0, v["caps_lo"]["g"]) + v["caps_lo"]["higher_per_core"])
    err_hi = max(0, cap_f_hi + max(0, v["caps_hi"]["g"]) + v["caps_hi"]["higher_per_core"])
    err_max = max(err_lo, err_hi)

    # Exact linear boundedness |S * Phi(T)| <= M_ticks * C(T) for ALL finite trees C(T) >= 1:
    cf = v["comp_fin"]
    min_ana = t["w_min_ticks"] + cf["g_far_min"]
    max_ana = t["w_max_ticks"] + cf["g_far_max"]
    m_lo = min(cf["min_W_c0"], cf["min_W"], t["w_min_ticks"], min_ana)
    m_hi = max(cf["max_W_c0"], cf["max_W"], t["w_max_ticks"], max_ana)
    M_ticks = max(abs(m_lo), abs(m_hi))

    theta_nom = F(v["theta_P"], v["theta_Q"])
    eps = F(err_max, S)
    # Inequality: obs_sign * (pi(psi_fin) - theta_nom) <= eps
    certified_fin = theta_nom + obs_sign * eps
    tail_add = F(0)
    if obs_id == 8 and obs_sign == 1:
        tail_add = bare_core_tail_upper(v["DF"], mu_lo)
    certified_full = certified_fin + tail_add

    res = {
        "ok": True,
        "sha256": h_raw,
        "obs_id": obs_id,
        "obs_name": v["obs_name"],
        "obs_sign": obs_sign,
        "mu_bracket": [str(mu_lo), str(mu_hi)],
        "theta_nominal": str(theta_nom),
        "err_lo_ticks": err_lo,
        "err_hi_ticks": err_hi,
        "err_max_ticks": err_max,
        "eps_exact": str(eps),
        "M_bound_ticks": M_ticks,
        "certified_fin_exact": str(certified_fin),
        "bare_tail_add_exact": str(tail_add),
        "certified_exact": str(certified_full),
        "certified_dec": format_dec(certified_full, 8, round_up=(obs_sign == 1)),
        "verify": v,
        "tail": t,
    }
    if override_theta is not None:
        res["theta_overridden"] = True
    if return_raw:
        res["_raw"] = raw
    return res


def probe_mutation(bin_path, verify_bin, perturb, expected_sha256):
    """Test-only mutation probe; never returns a certification receipt."""
    raw = Path(bin_path).read_bytes()
    h_raw = sha256_bytes(raw)
    if h_raw != expected_sha256:
        raise RuntimeError(f"frozen probe input mismatch: {h_raw} != {expected_sha256}")
    cmd = [str(verify_bin), str(bin_path), "--threads", "1", "--perturb", perturb]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0:
        return {"test_probe": True, "rejected": True, "returncode": proc.returncode, "stderr": proc.stderr.strip()}
    v = json.loads(proc.stdout)
    if not v.get("perturbed") or v.get("perturb_target") != perturb or v.get("unperturbed_sha256") != h_raw:
        raise RuntimeError(f"unexpected unperturbed metadata for probe {perturb}")
    if v.get("sha256") == h_raw:
        raise RuntimeError(f"perturbed buffer digest matched unperturbed digest for probe {perturb}")
    mu_lo = F(v["mu_lo_P"], v["mu_lo_Q"])
    mu_hi = F(v["mu_hi_P"], v["mu_hi_Q"])
    t = vmt.verify_bytes(raw, mu_lo=mu_lo, mu_hi=mu_hi, panels=16)
    cap_f_lo = max(v["caps_lo"]["f_in"], v["caps_lo"]["f_tab"], t["cap_ticks_lo"])
    cap_f_hi = max(v["caps_hi"]["f_in"], v["caps_hi"]["f_tab"], t["cap_ticks_hi"])
    err_lo = max(0, cap_f_lo + max(0, v["caps_lo"]["g"]) + v["caps_lo"]["higher_per_core"])
    err_hi = max(0, cap_f_hi + max(0, v["caps_hi"]["g"]) + v["caps_hi"]["higher_per_core"])
    return {
        "test_probe": True,
        "rejected": False,
        "returncode": 0,
        "perturbed_sha256": v["sha256"],
        "err_max_ticks": max(err_lo, err_hi),
    }


class PyMarginalSnapshot:
    """Independently parses FT2_M5V1 binary and checks exact generator inequality on concrete trees in Fraction."""

    def __init__(self, raw):
        self.mh = vmt.read_marginal_header_bytes(raw)
        self.h = self.mh["base"]
        h = self.h
        self.S = h["S"]
        self.obs_id = self.mh["obs_id"]
        self.obs_sign = self.mh["obs_sign"]
        self.theta = F(self.mh["theta_P"], self.mh["theta_Q"])
        self.mu_lo = F(self.mh["mu_lo_P"], self.mh["mu_lo_Q"])
        self.mu_hi = F(self.mh["mu_hi_P"], self.mh["mu_hi_Q"])

        def build_end(dc, Kc):
            idx, ovf, nxt = {}, {}, 0
            for d in range(2, dc):
                for c in range(1, min(d, Kc) + 1):
                    idx[(d, c)] = nxt; nxt += 1
            for d in range(Kc + 1, dc):
                ovf[d] = nxt; nxt += 1
            tb = nxt; nxt += Kc + 1
            inf_id = nxt
            return {"dc": dc, "Kc": Kc, "idx": idx, "ovf": ovf, "tb": tb, "inf": inf_id, "N": nxt + 1, "nin_ex": len(idx)}

        def build_inner(dc, mnC, mxC):
            idx, lst = {}, []
            for d in range(max(2, mnC), dc):
                for c in range(mnC, min(d, mxC) + 1):
                    idx[(d, c)] = len(lst); lst.append((d, c))
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

        self.p4_blk, n_p4 = {}, 0
        for b, (db, _) in enumerate(self.M4v["lst"]):
            for c, (dc_, _) in enumerate(self.M4v["lst"]):
                if b <= c and db + dc_ <= h["SP4"]:
                    self.p4_blk[(b, c)] = n_p4
                    n_p4 += (self.G4["N"] * (self.G4["N"] + 1) // 2) if b == c else (self.G4["N"] ** 2)

        self.p5_blk, n_p5 = {}, 0
        for b, (db, _) in enumerate(self.I5_2["lst"]):
            for c, (dc_, _) in enumerate(self.I5_2["lst"]):
                if db + dc_ <= h["SP4"]:
                    for d, (dd, _) in enumerate(self.I5_2["lst"]):
                        if b <= d and dc_ + dd <= h["SP4"] and db + dc_ + dd <= h["SP5"]:
                            self.p5_blk[(b, c, d)] = n_p5
                            n_p5 += (self.G5["N"] * (self.G5["N"] + 1) // 2) if b == d else (self.G5["N"] ** 2)

        self.b5_map, nb5 = {}, 0
        for o, (do, _) in enumerate(self.I5_3["lst"]):
            for c, (dc_, _) in enumerate(self.I5_2c["lst"]):
                if do + dc_ <= h["SB5"]:
                    self.b5_map[(o, c)] = nb5; nb5 += 1

        def tri_map(N):
            mp, k = {}, 0
            for a in range(N):
                for b in range(a, N):
                    for c in range(b, N):
                        mp[(a, b, c)] = k; k += 1
            return mp, k

        def quad_map(N):
            mp, k = {}, 0
            for a in range(N):
                for b in range(a, N):
                    for c in range(b, N):
                        for d in range(c, N):
                            mp[(a, b, c, d)] = k; k += 1
            return mp, k

        self.tri4, self.n_tri4 = tri_map(self.G4["N"])
        self.quad5, self.n_quad5 = quad_map(self.G5["N"])

        self.buf = memoryview(raw)[640:].cast("q")
        cur = 0
        def take(n):
            nonlocal cur
            o = cur; cur += n; return o
        self.off_r1 = take(self.G2["nin_ex"])
        self.off_knots = take(self.G2["nin_ex"] * h["nknots"])
        self.off_f = take((h["DF"] - 2) * (h["DF"] + 3) // 2)
        self.off_g = take(self.G2["N"] * (self.G2["N"] + 1) // 2)
        self.off_h3 = take(self.M3["N"] * (self.G3["N"] * (self.G3["N"] + 1) // 2))
        self.off_p4 = take(n_p4)
        self.off_s3 = take(self.MS3["N"] * self.n_tri4)
        self.off_p5 = take(n_p5)
        self.off_b5 = take(nb5 * (self.G5["N"] * (self.G5["N"] + 1) // 2) * self.G5["N"])
        self.off_s4 = take(self.MS4["N"] * self.n_quad5)

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
        q, x = F(D - c, D), F(1, D)
        return h["f_p0"] + h["f_p1"] * q + h["f_p2"] * q * q + self.mh["f_p3"] * q * x

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
        b = self.M4v["idx"].get((Db, cb)); c = self.M4v["idx"].get((Dc, cc))
        if b is None or c is None or Db + Dc > self.h["SP4"]:
            return F(0)
        if b > c:
            b, c = c, b; Da, ca, Dd, cd = Dd, cd, Da, ca
        base = self.p4_blk[(b, c)]; N4 = self.G4["N"]; res = F(0)
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
                    res += wa * wb * wc * self.buf[self.off_s3 + o * self.n_tri4 + self.tri4[tuple(sorted((a, b, c)))]]
        return res

    def eval_p5(self, Da, ca, Db, cb, Dc, cc, Dd, cd, De, ce):
        b = self.I5_2["idx"].get((Db, cb)); c = self.I5_2["idx"].get((Dc, cc)); d = self.I5_2["idx"].get((Dd, cd))
        if b is None or c is None or d is None:
            return F(0)
        if Db + Dc > self.h["SP4"] or Dc + Dd > self.h["SP4"] or Db + Dc + Dd > self.h["SP5"]:
            return F(0)
        if b > d:
            b, d = d, b; Da, ca, De, ce = De, ce, Da, ca
        base = self.p5_blk[(b, c, d)]; N5 = self.G5["N"]; res = F(0)
        for a, wa in self.wts(self.G5, Da, ca):
            for e, we in self.wts(self.G5, De, ce):
                off = self.sym(a, e) if b == d else (a * N5 + e)
                res += wa * we * self.buf[self.off_p5 + base + off]
        return res

    def eval_b5(self, Do, co, Da, ca, Db, cb, Dc, cc, Dd, cd):
        o = self.I5_3["idx"].get((Do, co)); c = self.I5_2c["idx"].get((Dc, cc))
        if o is None or c is None or Do + Dc > self.h["SB5"]:
            return F(0)
        bid = self.b5_map[(o, c)]; N5 = self.G5["N"]; p2 = N5 * (N5 + 1) // 2; res = F(0)
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
                        res += wa * wb * wc * wd * self.buf[self.off_s4 + o * self.n_quad5 + self.quad5[tuple(sorted((a, b, c, d)))]]
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
            nb = cadj[o]; k = len(nb)
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
            nb = cadj[o]; k = len(nb)
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
            nb = cadj[o]; k = len(nb)
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
        obs_sum = 0
        for u in range(n):
            cu = len(cadj[u])
            obs_sum += eval_obs(self.obs_id, D[u], cu, self.h["DF"])
            Lv = D[u] - cu
            gr = F(Lv) + sum(F(1, D[v]) for v in cadj[u])
            if gr > 0:
                D2 = list(D); D2[u] += 1
                Lphi += gr * (self.eval_phi(cadj, D2) - phi0)
            if Lv > 0:
                br = F(Lv, D[u])
                cadj2 = [list(x) for x in cadj] + [[u]]
                cadj2[u].append(n)
                D2 = list(D) + [2]
                Lphi += br * (self.eval_phi(cadj2, D2) - phi0)
        forcing = self.obs_sign * (F(obs_sum) - self.theta * n) * self.S
        res_lo = Lphi - self.mu_lo * phi0 + forcing
        res_hi = Lphi - self.mu_hi * phi0 + forcing
        return phi0, res_lo, res_hi


def run_tree_oracle(raw, side_res):
    if side_res.get("theta_overridden"):
        raise ValueError("run_tree_oracle uses nominal header theta and cannot run on theta_overridden results")
    snap = PyMarginalSnapshot(raw)
    DF, d2 = snap.h["DF"], snap.h["d2"]
    M_bnd = side_res["M_bound_ticks"]
    err_lo, err_hi = side_res["err_lo_ticks"], side_res["err_hi_ticks"]
    trees = []
    for D0 in (2, 3, 5, d2, DF - 1, DF, DF + 17):
        trees.append((f"star_C1_D{D0}", [[]], [D0]))
    for ln in (2, 3, 4, 5, 6):
        adj = [[] for _ in range(ln)]
        for i in range(ln - 1):
            adj[i].append(i + 1); adj[i + 1].append(i)
        trees.append((f"chain_{ln}", adj, [2] * ln))
        D2 = [2] * ln; D2[0] = d2 + 5; D2[-1] = DF + 11
        trees.append((f"chain_tail_{ln}", adj, D2))
    for arms in (3, 4, 6):
        adj = [[i for i in range(1, arms + 1)]] + [[0] for _ in range(arms)]
        trees.append((f"bare_star_{arms}", adj, [arms] + [2] * arms))
        trees.append((f"leaf_star_{arms}", adj, [arms + 1] + [2] * arms))
    rng = random.Random(20261124)
    for tr in range(12):
        n = rng.randint(6, 24)
        full = [[1, 2], [0], [0]] + [[] for _ in range(n - 3)]
        for k in range(3, n):
            v = rng.randrange(k)
            w = full[v][rng.randrange(len(full[v]))]
            full[w].append(k); full[k].append(w)
        cid = {i: idx for idx, i in enumerate(i for i in range(n) if len(full[i]) >= 2)}
        cadj = [[cid[j] for j in full[i] if j in cid] for i in cid]
        D = [len(full[i]) for i in cid]
        trees.append((f"rand_{tr}", cadj, D))

    worst_lo = None
    worst_hi = None
    for name, cadj, D in trees:
        C = len(cadj)
        phi0, r_lo, r_hi = snap.exact_tree_check(cadj, D)
        if abs(phi0) > M_bnd * C:
            raise AssertionError(f"Linear boundedness failed on {name}: |phi0|={abs(phi0)} > {M_bnd*C}")
        if r_lo > err_lo * C:
            raise AssertionError(f"Generator bound FAILED at mu_lo on {name}: {r_lo} > {err_lo*C}")
        if r_hi > err_hi * C:
            raise AssertionError(f"Generator bound FAILED at mu_hi on {name}: {r_hi} > {err_hi*C}")
        worst_lo = (r_lo / C) if (worst_lo is None or r_lo / C > worst_lo) else worst_lo
        worst_hi = (r_hi / C) if (worst_hi is None or r_hi / C > worst_hi) else worst_hi
    return {"trees_checked": len(trees), "worst_lo_per_core": str(worst_lo), "worst_hi_per_core": str(worst_hi)}


def run_mutations(bin_path, verify_bin, base_res, expected_sha256):
    results = {}
    env_neg = "ERROR: envelope radius negative or magnitude out of range"
    struct_expect = {
        "ginf":    {"ERROR: g[INF,INF] != 0"},
        "env23":   {"ERROR: env23 verification failed", env_neg},
        "env34p":  {"ERROR: env34p verification failed", env_neg},
        "env34s":  {"ERROR: env34s verification failed", env_neg},
        "env45p":  {"ERROR: env45p verification failed", env_neg},
        "env45bp": {"ERROR: env45bp/env45bs verification failed", env_neg},
        "env45bs": {"ERROR: env45bp/env45bs verification failed", env_neg},
        "env45s4": {"ERROR: env45s4 verification failed", env_neg},
    }
    for tgt, allowed in struct_expect.items():
        r = probe_mutation(bin_path, verify_bin, tgt, expected_sha256)
        if (not r["rejected"]) or r["returncode"] != 2 or r["stderr"] not in allowed:
            raise AssertionError(f"Planted structural mutation {tgt} was NOT rejected cleanly with expected diagnostic: {r}")
        results[tgt] = f"REJECTED (rc={r['returncode']}: {r['stderr']})"
    base_err = base_res["err_max_ticks"]
    allowed_range_reject = {"ERROR: weight out of range", "ERROR: r1 out of [0, S_r/2]"}
    for tgt in ("f", "g", "h3", "p4", "s3", "p5", "b5", "s4", "r1", "knots"):
        r = probe_mutation(bin_path, verify_bin, tgt, expected_sha256)
        if r["rejected"]:
            if r["returncode"] != 2 or r["stderr"] not in allowed_range_reject:
                raise AssertionError(f"Planted perturbation {tgt} failed unexpectedly (not range-clip or widened): {r}")
            results[tgt] = f"REJECTED (rc=2: {r['stderr']})"
        else:
            if r["err_max_ticks"] <= base_err + 100000:
                raise AssertionError(f"Planted perturbation {tgt} did not widen err_max: {r['err_max_ticks']} vs {base_err}")
            results[tgt] = f"WIDENED (err_max={r['err_max_ticks']} vs base {base_err})"
    th = F(base_res["verify"]["theta_P"], base_res["verify"]["theta_Q"])
    bad_th = th + (-F(1, 1000) if base_res["obs_sign"] == 1 else F(1, 1000))
    bad_th = max(F(0), min(F(1), bad_th))
    r_bad = combine_side(bin_path, verify_bin, threads=1,
                         override_theta=(bad_th.numerator, bad_th.denominator),
                         expected_sha256=expected_sha256)
    if (not r_bad["ok"]) or r_bad["err_max_ticks"] <= base_err + 10**8:
        raise AssertionError(f"Tighter theta probe did not cleanly widen err_max: {r_bad}")
    results["tighter_theta"] = f"WIDENED (err_max={r_bad['err_max_ticks']} vs base {base_err}, cert={r_bad['certified_dec']})"
    return results


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--verify-bin", required=True)
    ap.add_argument("--lower", required=True)
    ap.add_argument("--upper", required=True)
    ap.add_argument("--mu-lo", type=F, default=F(565713486333, 10**12))
    ap.add_argument("--mu-hi", type=F, default=F(565736986679, 10**12))
    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--run-mutations", action="store_true")
    ap.add_argument("--skip-trees", action="store_true")
    args = ap.parse_args()

    if not (F(0) < args.mu_lo <= args.mu_hi < F(1)):
        print(json.dumps({"error": f"invalid mu bracket [{args.mu_lo}, {args.mu_hi}]"}), file=sys.stderr)
        sys.exit(2)

    need_raw = (not args.skip_trees) or args.run_mutations
    lo = combine_side(args.lower, args.verify_bin, threads=args.threads, expected_sign=-1,
                      expected_mu_lo=args.mu_lo, expected_mu_hi=args.mu_hi, return_raw=need_raw)
    hi = combine_side(args.upper, args.verify_bin, threads=args.threads, expected_sign=1,
                      expected_mu_lo=args.mu_lo, expected_mu_hi=args.mu_hi, return_raw=need_raw)
    if not lo["ok"] or not hi["ok"]:
        lo.pop("_raw", None); hi.pop("_raw", None)
        print(json.dumps({"lower": lo, "upper": hi}, indent=2), file=sys.stderr)
        sys.exit(2)
    if lo["obs_id"] != hi["obs_id"]:
        print(json.dumps({"error": f"obs_id mismatch: {lo['obs_id']} vs {hi['obs_id']}"}), file=sys.stderr)
        sys.exit(2)
    lo_val, hi_val = F(lo["certified_exact"]), F(hi["certified_exact"])
    if lo_val > hi_val:
        print(json.dumps({"error": f"inverted interval: {lo_val} > {hi_val}"}), file=sys.stderr)
        sys.exit(2)

    raw_lo = lo.pop("_raw", None)
    raw_hi = hi.pop("_raw", None)
    out = {
        "obs_id": lo["obs_id"],
        "obs_name": lo["obs_name"],
        "mu_bracket_checked": [str(args.mu_lo), str(args.mu_hi)],
        "lower_sha256": lo.pop("sha256"),
        "upper_sha256": hi.pop("sha256"),
        "certified_interval_dec": [lo["certified_dec"], hi["certified_dec"]],
        "certified_width_dec": format_dec(hi_val - lo_val, 8, round_up=True),
        "lower": lo,
        "upper": hi,
    }
    if not args.skip_trees:
        out["tree_oracle_lower"] = run_tree_oracle(raw_lo, lo)
        out["tree_oracle_upper"] = run_tree_oracle(raw_hi, hi)
    if args.run_mutations:
        with tempfile.TemporaryDirectory(prefix="m5_frozen_") as td:
            fl = Path(td) / "lo.bin"; fh = Path(td) / "hi.bin"
            fl.write_bytes(raw_lo); fh.write_bytes(raw_hi)
            out["mutations_lower"] = run_mutations(fl, args.verify_bin, lo, out["lower_sha256"])
            out["mutations_upper"] = run_mutations(fh, args.verify_bin, hi, out["upper_sha256"])
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
