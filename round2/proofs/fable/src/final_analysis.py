#!/usr/bin/env python3
"""Final numbers from the archived per-chunk moments (moments_reduce.py).  Nothing here is a proof.

Every estimate is a ratio  R = MU_REF + sum e / sum d  over independent runs (see moments_reduce.py).  Two standard
errors are printed: "run" = delete-one-run jackknife in its linearised form, sqrt(n/(n-1) sum psi_r^2) with
psi_r = (e_r - delta d_r)/sum d, delta = sum e / sum d; "chunk" = delete-one-chunk jackknife.  The run-level error
is the one quoted.  Paired differences use the same runs for both estimates.
Usage: final_analysis.py MOMENTS.npz [--value-file data/final_value.txt]
"""
import sys, json, argparse
import numpy as np

ap = argparse.ArgumentParser(); ap.add_argument("npz"); ap.add_argument("--value-file", default=None)
args = ap.parse_args()
Z = np.load(args.npz, allow_pickle=False); META = json.loads(str(Z["meta"]))
NAMES5 = ["plain", "sane", "full", "AN", "ANPQ", "AR"]; NAMES2 = ["plain", "sane", "full"]
MU_LOW = 0.565                      # proved lower bound 113/200, used only to round the contraction factor upward
RHO_UP = float(np.exp(-4 * MU_LOW))     # 0.10435048...; printed rounded upward

def has(kind): return kind in META
def names(kind): return NAMES5 if META[kind]["nth"] == 5 else (NAMES2 if META[kind]["nth"] == 2 else ["plain"] + [f"th{i+1}" for i in range(META[kind]["nth"])])
def vecs(kind, obs, ta, tb):
    """selector vectors (ce, cd) for observable obs over the window ta -> tb of sample kind"""
    m = META[kind]; cps = m["checkpoints"]; nth = m["nth"]; k = names(kind).index(obs); ia, ib = cps.index(float(ta)), cps.index(float(tb))
    dim = 2 * (nth + 1) * (len(cps) - 1); ce = np.zeros(dim); cd = np.zeros(dim)
    for i in range(ia, ib): ce[(i * (nth + 1) + k) * 2] = 1.0; cd[(i * (nth + 1) + k) * 2 + 1] = 1.0
    return ce, cd
class Est:
    """ratio estimate with its per-chunk pieces and the linear functional v with psi_r = v.u_r"""
    def __init__(self, kind, obs, ta, tb):
        self.kind = kind; S = Z[kind + "_sum"]; ce, cd = vecs(kind, obs, ta, tb)
        self.Ec = S @ ce; self.Dc = S @ cd; self.E = self.Ec.sum(); self.D = self.Dc.sum(); self.delta = self.E / self.D
        self.R = META[kind]["mu_ref"] + self.delta; self.v = (ce - self.delta * cd) / self.D
        self.n = int(Z[kind + "_n"].sum()); self.G = len(self.Ec)
    def run_se(self): return quad_se(self.kind, self.v)
    def chunk_se(self):
        loo = (self.E - self.Ec) / (self.D - self.Dc); return float(np.sqrt((self.G - 1) / self.G * ((loo - loo.mean())**2).sum()))
def quad_se(kind, v):
    Q = Z[kind + "_sq"].sum(axis=0); n = int(Z[kind + "_n"].sum()); S = Z[kind + "_sum"].sum(axis=0)
    ss = float(v @ Q @ v) - float(v @ S)**2 / n          # sum psi^2 about the mean (the mean is zero up to rounding)
    return float(np.sqrt(max(ss, 0.0) * n / (n - 1)))
def paired(a, b):
    """difference R_a - R_b of two estimates on the same sample, run-level and chunk-level errors"""
    assert a.kind == b.kind
    d = a.R - b.R; se = quad_se(a.kind, a.v - b.v)
    loo = (a.E - a.Ec) / (a.D - a.Dc) - (b.E - b.Ec) / (b.D - b.Dc); G = a.G
    return d, se, float(np.sqrt((G - 1) / G * ((loo - loo.mean())**2).sum()))
def pooled(specs):
    """one ratio from several samples (sum of numerators over sum of denominators); specs = [(kind, obs, ta, tb), ...]"""
    parts = []
    for (kind, obs, ta, tb) in specs:
        ce, cd = vecs(kind, obs, ta, tb); S = Z[kind + "_sum"].sum(axis=0); parts.append((kind, ce, cd, float(S @ ce), float(S @ cd)))
    E = sum(x[3] for x in parts); D = sum(x[4] for x in parts); delta = E / D; n = sum(int(Z[x[0] + "_n"].sum()) for x in parts); ss = 0.0
    for (kind, ce, cd, _, _) in parts:
        u = (ce - delta * cd) / D; ss += float(u @ Z[kind + "_sq"].sum(axis=0) @ u)
    return META[parts[0][0]]["mu_ref"] + delta, float(np.sqrt(ss * n / (n - 1))), n
def fmt(x, s): return f"{x:.10f} +- {s:.2e}"
def up(x, nd=6):
    """round upward at nd decimals, exactly (coefficients of inequalities are never rounded down)"""
    from decimal import Decimal, ROUND_CEILING
    return float(Decimal.from_float(float(x)).quantize(Decimal(1).scaleb(-nd), rounding=ROUND_CEILING))
def upsig(x, sig=3):
    """round a positive bound upward at sig significant digits, exactly"""
    from decimal import Decimal, ROUND_CEILING
    if x <= 0: return x
    d = Decimal.from_float(float(x)); e = d.adjusted() - (sig - 1)
    return float(d.quantize(Decimal(1).scaleb(e), rounding=ROUND_CEILING))
def line(label, e): print(f"   {label:30s} {e.R:.10f}   run SE {e.run_se():.2e}   chunk SE {e.chunk_se():.2e}   (runs {e.n}, chunks {e.G})")

print(f"moments file {args.npz}; samples: " + ", ".join(f"{k} ({int(Z[k+'_n'].sum())} runs)" for k in META))
print(f"contraction factor rounded upward: rho_up = exp(-4*{MU_LOW}) = {RHO_UP:.7f}")
final = None

if has("main"):
    print("\n== 1. PRIMARY (fixed at 04:13 UTC): observable full, window 16 -> 24, runs made after all weights were frozen (sample main)")
    P = Est("main", "full", 16, 24); line("R_full(16,24)", P); final = (P.R, P.run_se())
    for o in ("sane", "AN", "ANPQ", "AR", "plain"): line(f"R_{o}(16,24)", Est("main", o, 16, 24))
    print("   windows of the same sample, observable full:")
    for (a, b) in ((12, 16), (16, 20), (20, 24), (12, 24)): line(f"R_full({a},{b})", Est("main", "full", a, b))
    print("   paired differences on the same runs (run SE; chunk SE):")
    for (a, b) in ((16, 24), (20, 24), (16, 20), (12, 16)):
        for o in ("sane", "AN", "ANPQ", "AR"):
            d, s, sc = paired(Est("main", "full", a, b), Est("main", o, a, b)); print(f"      full - {o:5s} at {a}->{b}: {d:+.3e} +- {s:.2e} ({d/s:+.1f} sigma; chunk SE {sc:.2e})")
    d, s, sc = paired(Est("main", "full", 20, 24), Est("main", "full", 16, 20)); print(f"      full(20,24) - full(16,20): {d:+.3e} +- {s:.2e} ({d/s:+.1f} sigma)")
    if has("main1"):
        print("   sensitivity: all runs, i.e. including the 3000 phase-1 runs (the first 600 of them prompted the choice of the window; all 3000 predate the three control observables) (samples main1 + main):")
        for o in ("full", "sane"):
            r, s, n = pooled([("main1", o, 16, 24), ("main", o, 16, 24)]); print(f"      R_{o}(16,24) all runs: {fmt(r, s)}  (runs {n})")

    print("\n== 2. INDEPENDENT RUNGS (phase-2 runs)")
    R0 = Est("r1216", "full", 12, 16) if has("r1216") else None; R1 = Est("r1620", "full", 16, 20) if has("r1620") else None
    for kind, (a, b) in (("r1216", (12, 16)), ("r1620", (16, 20))):
        if not has(kind): continue
        for o in ("full", "sane", "AN", "ANPQ", "AR"): line(f"{kind}: R_{o}({a},{b})", Est(kind, o, a, b))
        for o in ("sane", "AN", "ANPQ", "AR"):
            d, s, sc = paired(Est(kind, "full", a, b), Est(kind, o, a, b)); print(f"      full - {o:5s} at {a}->{b}: {d:+.3e} +- {s:.2e} ({d/s:+.1f} sigma; chunk SE {sc:.2e})")

    if has("b1216"):
        print("   fresh 12->16 cohort of the block test (stream 9210; not part of the registered contrast):")
        for o in ("full", "sane", "ANPQ"): line(f"b1216: R_{o}(12,16)", Est("b1216", o, 12, 16))
        r, s_, n_ = pooled([("r1216", "full", 12, 16), ("b1216", "full", 12, 16)]); print(f"      pooled R_full(12,16), r1216 + b1216: {fmt(r, s_)}  (runs {n_})")
        d, s_, sc = paired(Est("b1216", "full", 12, 16), Est("b1216", "ANPQ", 12, 16)); print(f"      full - ANPQ  at 12->16 (b1216): {d:+.3e} +- {s_:.2e} ({d/s_:+.1f} sigma)")
    print("\n== 3. CONDITIONAL SYSTEMATIC BOUNDS for the primary (each holds only under the premise named on its line; bounds rounded upward)")
    S = Z["main_sum"]; ce16, cd16 = vecs("main", "full", 16, 20); ce20, cd20 = vecs("main", "full", 20, 24)
    D16 = float(S.sum(axis=0) @ cd16); D20 = float(S.sum(axis=0) @ cd20); w = D16 / (D16 + D20)
    se_w = quad_se("main", (cd16 - w * (cd16 + cd20)) / (D16 + D20)); w_up = w + 1.96 * se_w
    print(f"   w = D(16,20)/[D(16,20)+D(20,24)] = {w:.9f} +- {se_w:.2e};  upper limit used {w_up:.9f};  growth of the window denominator D(20,24)/D(16,20) = {D20/D16:.5f}")
    for rho, tag in ((RHO_UP, "premise (C): |b(t+4)| <= rho |b(t)| at t = 12, 16 with rho = exp(-4 mu)"), (float(np.exp(-4 * (MU_LOW - 0.30))), "premise (C'): the same finite-window inequalities with the larger factor rho' = exp(-4(mu - 0.30))")):
        f = w_up + (1 - w_up) * rho; k2 = up(f * rho / (1 - f * rho)); k1 = up(f / (1 - f))
        print(f"   {tag}: factor <= {up(rho, 7):.7f}, f <= {up(f):.6f}, two-step coefficient f rho/(1 - f rho) <= {k2:.6f}, one-step coefficient f/(1-f) <= {k1:.6f} (rounded upward)")
        if R0 is not None:
            c = R0.R - P.R; sc = float(np.hypot(R0.run_se(), P.run_se())); print(f"      two-step: R_full(12,16)[r1216] - primary = {c:+.3e} +- {sc:.2e};  |primary - mu| <= {k2:.6f} * (|contrast| + 2.2414 SE) = {upsig(k2 * (abs(c) + 2.2414 * sc)):.2e}")
        if R1 is not None:
            c = R1.R - P.R; sc = float(np.hypot(R1.run_se(), P.run_se())); print(f"      one-step: R_full(16,20)[r1620] - primary = {c:+.3e} +- {sc:.2e};  |primary - mu| <= {k1:.6f} * (|contrast| + 2.2414 SE) = {upsig(k1 * (abs(c) + 2.2414 * sc)):.2e}")

if has("late"):
    print("\n== 4. LATE COHORT (pre-registered 06:00 UTC; checkpoints 20, 24, 28; independent of the samples above)")
    L1 = Est("late", "full", 24, 28); line("L1 = R_full(24,28)", L1)
    for o in ("sane", "AN", "ANPQ", "AR", "plain"): line(f"R_{o}(24,28)", Est("late", o, 24, 28))
    line("R_full(20,24) in the cohort", Est("late", "full", 20, 24)); line("R_full(20,28) in the cohort", Est("late", "full", 20, 28))
    d, s, sc = paired(L1, Est("late", "full", 20, 24)); print(f"   L2 = R_full(24,28) - R_full(20,24), paired: {d:+.3e} +- {s:.2e} ({d/s:+.1f} sigma)")
    if has("main"):
        c = L1.R - P.R; sc = float(np.hypot(L1.run_se(), P.run_se())); print(f"   L3 = L1 - primary: {c:+.3e} +- {sc:.2e} ({c/sc:+.1f} sigma)   [reading rule: agreement if within 2.5 sigma]")
        M20 = Est("main", "full", 20, 24); c2 = M20.R - L1.R; sc2 = float(np.hypot(M20.run_se(), L1.run_se())); kk = up(RHO_UP / (1 - RHO_UP))
        print(f"   bias bound for L1 under one contraction with factor rho at t = 20: |L1 - mu| <= rho/(1-rho) * (|R_full(20,24)[main] - L1| + 2.2414 SE) = {kk:.5f} * ({abs(c2):.2e} + 2.2414*{sc2:.2e}) = {upsig(kk * (abs(c2) + 2.2414 * sc2)):.2e}")
        if has("r1216"):
            R0l = Est("r1216", "full", 12, 16); c3 = R0l.R - L1.R; sc3 = float(np.hypot(R0l.run_se(), L1.run_se()))
            for rho, tag in ((RHO_UP, "premise (C+), factor rho"), (float(np.exp(-4 * (MU_LOW - 0.30))), "premise (C'+), factor rho'")):
                k3 = up(rho**3 / (1 - rho**3)); print(f"   bias bound for L1 under three contractions, {tag}: |L1 - mu| <= r^3/(1-r^3) * (|R_full(12,16)[r1216] - L1| + 2.2414 SE) = {k3:.6f} * ({abs(c3):.2e} + 2.2414*{sc3:.2e}) = {upsig(k3 * (abs(c3) + 2.2414 * sc3)):.2e}")
        wgt = 1 / P.run_se()**2, 1 / L1.run_se()**2; comb = (wgt[0] * P.R + wgt[1] * L1.R) / sum(wgt); print(f"   secondary, inverse-variance mean of primary and L1 (independent samples): {comb:.10f} +- {1/np.sqrt(sum(wgt)):.2e}")
    for o in ("sane", "AN", "ANPQ", "AR"):
        d, s, sc = paired(L1, Est("late", o, 24, 28)); print(f"   L4({o}) = full - {o} at 24->28: {d:+.3e} +- {s:.2e} ({d/s:+.1f} sigma)")
    for o in ("sane", "AN", "ANPQ", "AR"):
        d, s, sc = paired(Est("late", "full", 20, 24), Est("late", o, 20, 24)); print(f"   cohort, full - {o} at 20->24: {d:+.3e} +- {s:.2e} ({d/s:+.1f} sigma)")

if has("mainA"):
    print("\n== 4b. A-TEST 1 (registered 06:25 UTC): paired full - ANPQ at 20->24 on main chunks numbered >= 550 only")
    for o in ("ANPQ", "sane", "AN", "AR"):
        d, s, sc = paired(Est("mainA", "full", 20, 24), Est("mainA", o, 20, 24)); print(f"   full - {o:5s} at 20->24, chunks >= 550: {d:+.3e} +- {s:.2e} ({d/s:+.1f} sigma; runs {Est('mainA', 'full', 20, 24).n})")
    line("R_full(16,24), chunks >= 550", Est("mainA", "full", 16, 24))
ref = final[0] if final else 0.5657236
print(f"\n== 5. EARLY WINDOWS AND KNOWN-TRUTH LADDERS (offsets are estimate minus reference)")
for kind, label, refv in (("e0408", "friend tree 4->8 (stream 9410)", ref), ("e0812", "friend tree 8->12 (stream 9420)", ref), ("e0408five", "friend tree 4->8 (stream 9440)", ref), ("e0812five", "friend tree 8->12 (stream 9430)", ref),
                          ("ir0early", "IR-0 (stream 7303)", 0.4822669183233687), ("ir0ladder", "IR-0 (stream 7302)", 0.4822669183233687)):
    if not has(kind): continue
    cps = META[kind]["checkpoints"]
    for i in range(len(cps) - 1):
        out = []
        for o in names(kind)[1:]:
            e = Est(kind, o, cps[i], cps[i + 1]); out.append(f"{o}: {e.R - refv:+.3e} +- {e.run_se():.1e}")
        print(f"   {label}, window {cps[i]:g}->{cps[i+1]:g}, reference {refv:.10f}: " + " | ".join(out))

if final and args.value_file:
    open(args.value_file, "w").write(f"{final[0]:.10f} {final[1]:.2e}\n"); print(f"\nwrote {args.value_file}: {final[0]:.10f} {final[1]:.2e}")
