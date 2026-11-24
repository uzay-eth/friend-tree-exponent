#!/usr/bin/env python3
"""Inverse-symbolic search around a numerical value mu0 +- sigma.

For each family of candidate closed forms we count the members that fall inside the window
[mu0 - K sigma, mu0 + K sigma] and compare with the count in NCTRL control windows of the same width centred at
mu0 + j * step (j = +-1, +-2, ...), which measures how many hits are expected by chance for a family of that size.
A family "has no evidence" when the hit count at mu0 is within the range of the control counts.

Families
  rat      p/q, 1 <= q <= Q
  alg d    roots of integer polynomials of degree d with coefficients |a_i| <= H_d (d = 2..6)
  surd     (a + b sqrt2 + c sqrt13 + d sqrt26)/e, |a|,|b|,|c|,|d| <= H, 1 <= e <= E
  hkrz     P(mu) + Q(mu) F(mu) = 0, F(mu) = int_0^1 x^mu e^(x-1) dx, deg P, deg Q <= 2, |coeff| <= H
           (the equation family of the solvable DAN / PAN / VAN models: R(mu) = F(mu) with R rational of degree <= 2)
  hkrzg    same with F replaced by G_z(mu) = int_0^1 x^mu e^(z(x-1)) dx for z in {2, 1/2, -1} and by Gamma(1+mu)
  named    a list of named constants and simple expressions

Usage: isc_search.py MU0 SIGMA [--K 3] [--nctrl 200] [--quick]
"""
import argparse, itertools, math, sys
import numpy as np
from mpmath import mp, mpf, quad, exp, gamma, sqrt, lambertw, log, pi, euler, zeta, findroot

ap = argparse.ArgumentParser(); ap.add_argument("mu0", type=float); ap.add_argument("sigma", type=float)
ap.add_argument("--K", type=float, default=3.0); ap.add_argument("--nctrl", type=int, default=200); ap.add_argument("--quick", action="store_true")
args = ap.parse_args()
mp.dps = 30
mu0, sig, K = args.mu0, args.sigma, args.K
w = K * sig
step = max(50 * w, 2e-6)                      # control windows well separated from the target window
ctrl = [mu0 + j * step for j in range(-args.nctrl // 2, args.nctrl // 2 + 1) if j != 0]
print(f"target {mu0:.10f} +- {sig:.1e}; window half-width {w:.2e} ({K} sigma); {len(ctrl)} control windows at spacing {step:.1e}")

def report(name, size, hits, ctrl_counts, examples):
    cc = np.array(ctrl_counts, float)
    print(f"[{name}] candidates searched ~{size:.3g}; hits in target window: {hits}; control windows: mean {cc.mean():.2f}, max {int(cc.max())}")
    for e in examples[:12]: print("      ", e)

# ---------- rationals ----------
def rat_hits(x, Q):
    out = []
    for q in range(1, Q + 1):
        p = round(x * q)
        if abs(p / q - x) <= w and math.gcd(p, q) == 1: out.append((p, q))
    return out
Q = 3000 if args.quick else 20000
h = rat_hits(mu0, Q)
report(f"rat q<={Q}", 0.3 * Q * Q, len(h), [len(rat_hits(c, Q)) for c in ctrl], [f"{p}/{q} = {p/q:.10f}" for p, q in h])

# ---------- algebraic numbers ----------
def alg_hits(x, d, H, collect=False):
    """polynomials sum a_i x^i, a_d != 0 (a_d > 0), |a_i| <= H: choose a_1..a_d, solve a_0 = round(-sum); accept if a root lies in window"""
    rng = np.arange(-H, H + 1)
    pw = np.array([x ** i for i in range(1, d + 1)])
    grids = np.meshgrid(*([rng] * (d - 1) + [np.arange(1, H + 1)]), indexing="ij")
    coef = np.stack([gg.ravel() for gg in grids], axis=1).astype(np.int64)            # a_1..a_d
    val = coef @ pw
    a0 = -np.rint(val)
    ok = np.abs(a0) <= H
    resid = val + a0                                                                   # P(x)
    deriv = coef @ np.array([i * x ** (i - 1) for i in range(1, d + 1)])
    with np.errstate(divide="ignore", invalid="ignore"):
        dx = np.abs(resid / deriv)                                                     # distance to the root (Newton)
    sel = ok & (dx <= w) & (np.abs(deriv) > 1e-6)
    idx = np.where(sel)[0]
    res = set()
    for i in idx:
        a = (int(a0[i]),) + tuple(int(t) for t in coef[i])
        g = 0
        for t in a: g = math.gcd(g, abs(t))
        if g > 1: continue                                                             # non-primitive: counted at smaller height
        res.add(a)
    return res if collect else len(res)
plan = [(2, 60 if args.quick else 400), (3, 20 if args.quick else 60), (4, 10 if args.quick else 22), (5, 6 if args.quick else 11), (6, 4 if args.quick else 7)]
for d, H in plan:
    hs = alg_hits(mu0, d, H, collect=True)
    # drop polynomials that are reducible with a lower-degree factor already in the window is not attempted; listed as found
    cc = [alg_hits(c, d, H) for c in ctrl[:: max(1, len(ctrl) // 40)]]
    report(f"alg deg {d}, |coeff|<={H}", (2 * H + 1) ** d * H, len(hs), cc, [f"coeffs a0..a{d} = {a}" for a in sorted(hs, key=lambda t: max(abs(v) for v in t))])

# ---------- Q(sqrt2, sqrt13) ----------
def surd_hits(x, H, E, collect=False):
    s2, s13, s26 = math.sqrt(2), math.sqrt(13), math.sqrt(26)
    rng = np.arange(-H, H + 1)
    B, Cc, Dd = np.meshgrid(rng, rng, rng, indexing="ij")
    base = (B * s2 + Cc * s13 + Dd * s26).ravel(); bb = B.ravel(); cc_ = Cc.ravel(); dd = Dd.ravel()
    res = set()
    for e in range(1, E + 1):
        a = np.rint(e * x - base)
        ok = (np.abs(a) <= H) & (np.abs((a + base) / e - x) <= w)
        for i in np.where(ok)[0]:
            t = (int(a[i]), int(bb[i]), int(cc_[i]), int(dd[i]), e)
            g = 0
            for v in t: g = math.gcd(g, abs(v))
            if g == 1 and (t[1] or t[2] or t[3]): res.add(t)
    return res if collect else len(res)
H, E = (8, 12) if args.quick else (20, 40)
hs = surd_hits(mu0, H, E, collect=True)
report(f"surd (a+b√2+c√13+d√26)/e, |a..d|<={H}, e<={E}", (2 * H + 1) ** 3 * E, len(hs), [surd_hits(c, H, E) for c in ctrl[:: max(1, len(ctrl) // 50)]],
       [f"(a,b,c,d,e) = {t}" for t in sorted(hs, key=lambda t: max(abs(v) for v in t))])

# ---------- HKRZ family ----------
def Fz(mu, z=1):
    return quad(lambda x: x ** mu * exp(z * (x - 1)), [0, 1])
def dnum(f, x, hh=1e-6):
    return float((f(mpf(x) + hh) - f(mpf(x) - hh)) / (2 * hh))
def trans_hits(x, f, df, H, collect=False):
    """P(mu) + Q(mu) T(mu) = 0 with P = p0 + p1 mu + p2 mu^2, Q = q0 + q1 mu + q2 mu^2, integers |.| <= H, Q != 0."""
    rng = np.arange(-H, H + 1)
    P1, P2, Q0, Q1, Q2 = np.meshgrid(rng, rng, rng, rng, rng, indexing="ij")
    P1, P2, Q0, Q1, Q2 = [a.ravel().astype(np.float64) for a in (P1, P2, Q0, Q1, Q2)]
    val = P1 * x + P2 * x * x + (Q0 + Q1 * x + Q2 * x * x) * f
    p0 = -np.rint(val); resid = val + p0
    deriv = P1 + 2 * P2 * x + (Q1 + 2 * Q2 * x) * f + (Q0 + Q1 * x + Q2 * x * x) * df
    with np.errstate(divide="ignore", invalid="ignore"): dx = np.abs(resid / deriv)
    ok = (np.abs(p0) <= H) & (dx <= w) & (np.abs(deriv) > 1e-6) & ((Q0 != 0) | (Q1 != 0) | (Q2 != 0))
    res = set()
    for i in np.where(ok)[0]:
        t = tuple(int(v) for v in (p0[i], P1[i], P2[i], Q0[i], Q1[i], Q2[i]))
        g = 0
        for v in t: g = math.gcd(g, abs(v))
        if g != 1: continue
        lead = next(v for v in t if v != 0)
        if lead < 0: t = tuple(-v for v in t)
        res.add(t)
    return res if collect else len(res)
H = 5 if args.quick else 9
fams = [("hkrz F(mu)=int_0^1 x^mu e^(x-1)dx", lambda m: Fz(m, 1)), ("G_2(mu)=int x^mu e^(2(x-1))", lambda m: Fz(m, 2)),
        ("G_1/2(mu)", lambda m: Fz(m, mpf(1) / 2)), ("G_-1(mu)=int x^mu e^(1-x)", lambda m: Fz(m, -1)), ("Gamma(1+mu)", lambda m: gamma(1 + m))]
for name, f in fams:
    f0 = float(f(mpf(mu0))); d0 = dnum(f, mu0)
    hs = trans_hits(mu0, f0, d0, H, collect=True)
    cc = []
    for c in ctrl[:: max(1, len(ctrl) // 50)]:
        cc.append(trans_hits(c, float(f(mpf(c))), dnum(f, c), H))
    report(f"{name}: P+Q*T=0, deg<=2, |coeff|<={H}", (2 * H + 1) ** 5, len(hs), cc,
           [f"(p0,p1,p2 | q0,q1,q2) = {t}" for t in sorted(hs, key=lambda t: max(abs(v) for v in t))])
# the two exact HKRZ equations themselves, evaluated at mu0
F0 = Fz(mpf(mu0))
print(f"   F(mu0) = {float(F0):.10f};  1-mu0 = {1-mu0:.10f} (PAN eq.), 1-mu0^2 = {1-mu0*mu0:.10f} (DAN eq.)")

# ---------- VAN(omega) family: (1-mu)(mu+omega)/(omega+(1-omega)mu) = F(mu): which omega would give mu0? ----------
F0f = float(F0)
# solve (1-mu)(mu+w) = F (w + (1-w) mu)  ->  w [(1-mu) - F(1-mu)] = F mu - mu(1-mu)
w_van = (F0f * mu0 - mu0 * (1 - mu0)) / ((1 - mu0) * (1 - F0f))
F1f = float(Fz(mpf(mu0 + sig)))        # sensitivity of omega to mu: the full formula re-evaluated at mu0 + sigma
dw = abs((F1f * (mu0 + sig) - (mu0 + sig) * (1 - mu0 - sig)) / ((1 - mu0 - sig) * (1 - F1f)) - w_van)
print(f"[van] the VAN(omega) model of HKRZ has exponent mu0 for omega = {w_van:.8f} +- {dw:.1e} (PAN: 0, DAN: 1); nearest rationals with q <= 50: "
      + ", ".join(f"{p}/{q}" for q in range(1, 51) for p in [round(w_van * q)] if abs(p / q - w_van) < 3 * K * dw and math.gcd(p, q) == 1) + " (none listed = none)")

# ---------- elementary transcendental families ----------
def small_family(name, funcs, H):
    """a*f1(mu) + b*f2(mu) + c*f3(mu) + ... = 0 with integer coefficients |.| <= H (first non-zero positive)"""
    def count(x):
        vals = np.array([float(f(mpf(x))) for f in funcs]); der = np.array([dnum(f, x) for f in funcs]); n = len(funcs)
        rng = np.arange(-H, H + 1); grids = np.meshgrid(*([rng] * (n - 1)), indexing="ij")
        co = np.stack([gg.ravel() for gg in grids], axis=1).astype(float)            # coefficients of funcs[1:]
        rest = co @ vals[1:]; a0 = -np.rint(rest / vals[0])                             # coefficient of funcs[0]
        resid = a0 * vals[0] + rest; deriv = a0 * der[0] + co @ der[1:]
        with np.errstate(divide="ignore", invalid="ignore"): dx = np.abs(resid / deriv)
        ok = (np.abs(a0) <= H) & (dx <= w) & (np.abs(deriv) > 1e-9) & (np.abs(co).sum(axis=1) > 0) & (a0 != 0)
        res = set()
        for i in np.where(ok)[0]:
            t = (int(a0[i]),) + tuple(int(v) for v in co[i]); g = 0
            for v in t: g = math.gcd(g, abs(v))
            if g != 1: continue
            lead = next(v for v in t if v != 0)
            res.add(t if lead > 0 else tuple(-v for v in t))
        return res
    hs = count(mu0); cc = [len(count(c)) for c in ctrl[:: max(1, len(ctrl) // 50)]]
    report(f"{name}, |coeff|<={H}", (2 * H + 1) ** (len(funcs) - 1), len(hs), cc, [str(t) for t in sorted(hs, key=lambda t: max(abs(v) for v in t))])
one = lambda m: mpf(1)
H3 = 12 if args.quick else 40
small_family("a + b mu + c log(mu) + d mu log(mu) = 0", [one, lambda m: m, lambda m: log(m), lambda m: m * log(m)], H3)
small_family("a + b mu + c e^mu + d mu e^mu = 0", [one, lambda m: m, lambda m: exp(m), lambda m: m * exp(m)], H3)
# (the family with e^(-mu) and mu e^(-mu) is the previous one multiplied by e^(-mu); it is not searched twice)
small_family("a + b mu + c e^mu + d e^(-mu) = 0", [one, lambda m: m, lambda m: exp(m), lambda m: exp(-m)], H3)
small_family("a + b mu + c 2^mu + d 3^mu = 0", [one, lambda m: m, lambda m: mpf(2) ** m, lambda m: mpf(3) ** m], H3)
small_family("a + b cos(pi mu) + c sin(pi mu) + d mu = 0", [one, lambda m: mp.cos(pi * m), lambda m: mp.sin(pi * m), lambda m: m], H3)
small_family("a + b mu + c Gamma(mu) + d Gamma(1+mu)^2 = 0", [one, lambda m: m, lambda m: gamma(m), lambda m: gamma(1 + m) ** 2], H3)
small_family("a + b mu + c psi(mu) + d psi(1+mu) mu = 0", [one, lambda m: m, lambda m: mp.digamma(m), lambda m: m * mp.digamma(1 + m)], H3)
small_family("a + b mu + c M(1,2+mu,1) + d mu M(1,2+mu,1) = 0  (Kummer, z=+1)", [one, lambda m: m, lambda m: mp.hyp1f1(1, 2 + m, 1), lambda m: m * mp.hyp1f1(1, 2 + m, 1)], H3)
small_family("a + b mu + c U(1,2+mu,1) + d mu U(1,2+mu,1) = 0  (Tricomi)", [one, lambda m: m, lambda m: mp.hyperu(1, 2 + m, 1), lambda m: m * mp.hyperu(1, 2 + m, 1)], H3)

# ---------- integer relations with classical constants (PSLQ-type, small coefficients) ----------
basis = {"1": mpf(1), "pi": pi, "e": exp(1), "log2": log(2), "gamma": euler, "sqrt2": sqrt(2), "sqrt3": sqrt(3), "sqrt5": sqrt(5), "sqrt13": sqrt(13), "zeta3": zeta(3), "1/pi": 1 / pi, "1/e": exp(-1)}
names = list(basis); Hc = 4 if args.quick else 8
def lin_hits(x):
    cnt = 0; ex = []
    for i, j in itertools.combinations(range(1, len(names)), 2):
        ci, cj = float(basis[names[i]]), float(basis[names[j]])
        rng = np.arange(-Hc, Hc + 1); Bm, Cm, Dm = np.meshgrid(rng, rng, np.arange(1, Hc + 1), indexing="ij")
        b, c, d = Bm.ravel(), Cm.ravel(), Dm.ravel()                       # d*mu = a + b*ci + c*cj
        a = np.rint(d * x - b * ci - c * cj); val = (a + b * ci + c * cj) / d
        ok = (np.abs(a) <= Hc) & (np.abs(val - x) <= w) & ((b != 0) | (c != 0))
        for t in np.where(ok)[0]:
            g = 0
            for v in (int(a[t]), int(b[t]), int(c[t]), int(d[t])): g = math.gcd(g, abs(v))
            if g == 1:
                cnt += 1
                if len(ex) < 12: ex.append(f"({int(a[t])} + {int(b[t])}*{names[i]} + {int(c[t])}*{names[j]})/{int(d[t])} = {val[t]:.10f}")
    return cnt, ex
hc, hex_ = lin_hits(mu0); cc = [lin_hits(c)[0] for c in ctrl[:: max(1, len(ctrl) // 50)]]
report(f"(a + b X + c Y)/d with X, Y two of {names[1:]}, |a|,|b|,|c|,d <= {Hc}", 66 * (2 * Hc + 1) ** 2 * Hc, hc, cc, hex_)

# ---------- named constants ----------
named = {
    "17/30 (Strogatz 2022)": mpf(17) / 30, "W(1), omega constant": lambertw(1), "1/sqrt(pi)": 1 / sqrt(pi), "exp(-gamma)": exp(-euler),
    "2 sqrt(2)/5": 2 * sqrt(2) / 5, "mu_PAN (HKRZ)": findroot(lambda m: 1 - m - Fz(m), 0.55), "mu_DAN (HKRZ)": findroot(lambda m: 1 - m * m - Fz(m), 0.77),
    "IR-0 (no core-to-core sends)": mpf("0.4822669183233687"), "(sqrt(13)-3)/6 (AB et al. lower)": (sqrt(13) - 3) / 6, "sqrt(2)-1/2 (AB et al. upper)": sqrt(2) - mpf(1) / 2,
    "113/200 (certified lower bound)": mpf(113) / 200, "1133/2000 (certified upper bound)": mpf(1133) / 2000, "13/23": mpf(13) / 23, "43/76": mpf(43) / 76, "56/99": mpf(56) / 99, "log(2)^(3/2)": log(2) ** mpf(1.5), "pi/sqrt(31)": pi / sqrt(31),
}
print("[named] distance from mu0 in units of sigma:")
for k, v in named.items():
    print(f"       {k:34s} = {float(v):.10f}   ({(float(v) - mu0) / sig:+.1f} sigma)")
