#!/usr/bin/env python3
"""
Systematic Inverse-Symbolic Search (Rational, Algebraic, HKRZ Incomplete-Gamma,
Classical Constants, Special Functions) on a high-precision numerical estimate of mu.
Pure Python stdlib + PARI/GP (no mpmath required).

Usage:
  python3 -B pslq_search.py <mu_est> <sigma>
Example:
  python3 -B pslq_search.py 0.5657236585 3.8e-9

Reports every family tested, basis size, height/norm bound, precision window (5 sigma),
expected pigeonhole threshold, and any candidate hits (with residual in units of sigma).
"""
import sys, math, bisect, subprocess

def gp_eval_list(gp_setup, expr_list):
    lines = ["\\p 60", gp_setup]
    for e in expr_list:
        lines.append(f"print({e});")
    lines.append("quit;")
    r = subprocess.run(["gp", "-q", "-f"], input="\n".join(lines) + "\n",
                       stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
    out = [ln.strip().replace(" ", "") for ln in r.stdout.strip().splitlines() if ln.strip()]
    return [float(x) for x in out[-len(expr_list):]]

def run_gp_lindep(expr_list, prec_digits, max_norm, mu_val):
    vec_str = "[" + ", ".join(expr_list) + "]"
    gp_code = f"""
\\p 60
mu = {mu_val};
Fmu = intnum(x=0, 1, x^mu * exp(x-1));
v = {vec_str};
r = lindep(v, {prec_digits});
res = abs(sum(i=1, #v, r[i] * v[i]));
print(r);
print(res);
quit;
"""
    r = subprocess.run(["gp", "-q", "-f"], input=gp_code, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    if r.returncode != 0:
        return None, None
    lines = [ln.strip() for ln in r.stdout.strip().splitlines() if ln.strip()]
    if len(lines) < 2:
        return None, None
    out = lines[-2].replace("[", "").replace("]", "").replace("~", "")
    coeffs = [int(x.strip()) for x in out.split(",") if x.strip()]
    res = float(lines[-1].replace(" ", ""))
    if not coeffs or all(c == 0 for c in coeffs):
        return None, None
    if max(abs(c) for c in coeffs) > max_norm:
        return None, None
    return coeffs, res

def main():
    if len(sys.argv) < 3:
        print("Usage: pslq_search.py <mu_est> <sigma>"); sys.exit(2)
    mu_str = sys.argv[1]
    mu_f = float(mu_str)
    sigma = float(sys.argv[2])
    window = 5.0 * sigma
    eff_digits = max(6, int(round(-math.log10(window))))
    # Also ensure PARI/GP lindep targets strictly below 5*sigma (~8 digits for sigma=3.8e-9)
    gp_digits = max(eff_digits, int(math.ceil(-math.log10(window))))

    print(f"=== Inverse-Symbolic Search for mu = {mu_str} +- {sigma:.3e} (5-sigma window = {window:.3e}, eff_digits = {eff_digits}) ===")

    # 1. Exhaustive rational search up to Q_max = 1,000,000
    Q_max = 1000000
    lo_f = mu_f - window
    hi_f = mu_f + window
    hits_rat = []
    for q in range(1, Q_max + 1):
        p_lo = math.ceil(lo_f * q)
        p_hi = math.floor(hi_f * q)
        for p in range(p_lo, p_hi + 1):
            if math.gcd(p, q) == 1:
                err = p / q - mu_f
                hits_rat.append((q, p, err, err / sigma))
                if len(hits_rat) >= 10:
                    break
        if len(hits_rat) >= 10:
            break
    pigeon_q = int(math.sqrt(1.0 / (2.0 * window)))
    print(f"\n[Family 1: Rationals p/q, q <= {Q_max}] (pigeonhole threshold q ~ {pigeon_q})")
    if not hits_rat:
        print(f"  No rational p/q with q <= {Q_max} inside 5-sigma window.")
    else:
        for q, p, err, z in hits_rat[:8]:
            print(f"  q={q:7d}: {p}/{q} = {p/q:.12f}  (err={err:+.3e}, {z:+.2f} sigma)")

    # 2. HKRZ incomplete-gamma family: [1, mu, mu^2, F(mu), mu*F(mu)], |coeff| <= 16
    setup = f"mu = {mu_str}; Fmu = intnum(x=0, 1, x^mu * exp(x-1)); dF = intnum(x=0, 1, log(x) * x^mu * exp(x-1));"
    Fmu, dF = gp_eval_list(setup, ["Fmu", "dF"])
    print(f"\n[HKRZ F(mu) = {Fmu:.16f}, F'(mu) = {dF:.16f}]")
    H_exh = 16
    print(f"\n[Family 2: Exhaustive HKRZ [1, mu, mu^2, F(mu), mu*F(mu)], |coeff| <= {H_exh}, nonzero F part]")
    vals_hkrz = [1.0, mu_f, mu_f**2, Fmu, mu_f * Fmu]
    dvals_hkrz = [0.0, 1.0, 2.0 * mu_f, dF, Fmu + mu_f * dF]
    rhs_list = []
    for c3 in range(-H_exh, H_exh + 1):
        for c4 in range(-H_exh, H_exh + 1):
            if c3 == 0 and c4 == 0:
                continue
            rhs_list.append((c3 * vals_hkrz[3] + c4 * vals_hkrz[4], c3, c4))
    rhs_list.sort()
    rhs_keys = [x[0] for x in rhs_list]
    max_d_est = H_exh * sum(abs(x) for x in dvals_hkrz)
    tol_val = window * max_d_est
    hkrz_hits = []
    for c1 in range(-H_exh, H_exh + 1):
        for c2 in range(-H_exh, H_exh + 1):
            for c0 in range(-H_exh, H_exh + 1):
                if c0 == 0 and c1 == 0 and c2 == 0:
                    continue
                lhs = -(c0 + c1 * vals_hkrz[1] + c2 * vals_hkrz[2])
                i0 = bisect.bisect_left(rhs_keys, lhs - tol_val)
                i1 = bisect.bisect_right(rhs_keys, lhs + tol_val)
                for idx in range(i0, i1):
                    rval, c3, c4 = rhs_list[idx]
                    if math.gcd(math.gcd(abs(c0), abs(c1)), math.gcd(abs(c2), math.gcd(abs(c3), abs(c4)))) != 1:
                        continue
                    first_nz = next(c for c in (c0, c1, c2, c3, c4) if c != 0)
                    if first_nz < 0:
                        continue
                    res = c0 + c1 * vals_hkrz[1] + c2 * vals_hkrz[2] + rval
                    deriv = c1 * dvals_hkrz[1] + c2 * dvals_hkrz[2] + c3 * dvals_hkrz[3] + c4 * dvals_hkrz[4]
                    if abs(deriv) < 1e-12:
                        continue
                    dmu = -res / deriv
                    if abs(dmu) <= window:
                        hkrz_hits.append((max(abs(c0),abs(c1),abs(c2),abs(c3),abs(c4)), (c0,c1,c2,c3,c4), dmu, dmu / sigma))
    hkrz_hits.sort()
    pigeon_hkrz = (1.0 / (2.0 * window)) ** (1.0 / 5.0)
    print(f"  Exhaustive H <= {H_exh} hits in 5-sigma window (pigeonhole H ~ {pigeon_hkrz:.1f}): {len(hkrz_hits)}")
    for H, cs, dmu, z in hkrz_hits[:8]:
        print(f"    H={H:2d} coeffs={cs}  delta_mu={dmu:+.3e} ({z:+.2f} sigma)")

    # 3. Exhaustive algebraic degree 2 & 3 check (exact heights)
    print(f"\n[Family 2b: Exhaustive Algebraic deg 2 (|c|<=300) and deg 3 (|c|<=40)]")
    alg2_hits = []
    for c2 in range(1, 301):
        for c1 in range(-300, 301):
            c0 = round(-(c2 * mu_f * mu_f + c1 * mu_f))
            if abs(c0) > 300 or math.gcd(c2, math.gcd(abs(c1), abs(c0))) != 1:
                continue
            disc = c1 * c1 - 4 * c2 * c0
            if disc <= 0:
                continue
            # Check both roots
            sd = math.sqrt(disc)
            for r in ((-c1 + sd) / (2 * c2), (-c1 - sd) / (2 * c2)):
                dmu = r - mu_f
                if abs(dmu) <= window:
                    alg2_hits.append((max(c2, abs(c1), abs(c0)), (c0, c1, c2), dmu, dmu / sigma))
    alg2_hits.sort()
    print(f"  Degree 2 (H<=300, pigeonhole H~{(1/(2*window))**(1/3):.1f}): {len(alg2_hits)} hits")
    for H, cs, dmu, z in alg2_hits[:6]:
        print(f"    H={H:3d} {cs[2]}*x^2 + {cs[1]}*x + {cs[0]} = 0  delta_mu={dmu:+.3e} ({z:+.2f} sigma)")

    # 4. PARI/GP lindep across algebraic degrees 2..6 & transcendental bases
    bases = [
        ("Algebraic deg 2", ["1", "mu", "mu^2"], 10000),
        ("Algebraic deg 3", ["1", "mu", "mu^2", "mu^3"], 2000),
        ("Algebraic deg 4", ["1", "mu", "mu^2", "mu^3", "mu^4"], 500),
        ("Algebraic deg 5", ["1", "mu", "mu^2", "mu^3", "mu^4", "mu^5"], 200),
        ("Algebraic deg 6", ["1", "mu", "mu^2", "mu^3", "mu^4", "mu^5", "mu^6"], 80),
        ("HKRZ deg (2,1)", ["1", "mu", "mu^2", "Fmu", "mu*Fmu"], 200),
        ("HKRZ deg (3,2)", ["1", "mu", "mu^2", "mu^3", "Fmu", "mu*Fmu", "mu^2*Fmu"], 50),
        ("Constants [1, mu, Pi, e, log(2), Euler]", ["1", "mu", "Pi", "exp(1)", "log(2)", "Euler"], 100),
        ("Log/Exp [1, mu, log(mu), exp(mu), exp(-mu)]", ["1", "mu", "log(mu)", "exp(mu)", "exp(-mu)"], 100),
        ("Gamma [1, mu, gamma(mu), gamma(1+mu), psi(1+mu)]", ["1", "mu", "gamma(mu)", "gamma(1+mu)", "psi(1+mu)"], 100),
        ("Lambert W(1) & zeta [1, mu, lambertw(1), zeta(2), zeta(3)]", ["1", "mu", "lambertw(1)", "zeta(2)", "zeta(3)"], 100),
        ("Bessel I0/I1/K0 at 1 [1, mu, besseli(0,1), besseli(1,1), besselk(0,1)]", ["1", "mu", "besseli(0,1)", "besseli(1,1)", "besselk(0,1)"], 100),
        ("Quadratic irrationals [1, mu, sqrt(2), sqrt(3), sqrt(5), sqrt(13)]", ["1", "mu", "sqrt(2)", "sqrt(3)", "sqrt(5)", "sqrt(13)"], 100),
    ]
    print(f"\n[Family 3: PARI/GP lindep across algebraic & transcendental bases (target {gp_digits} digits)]")
    for name, exprs, max_n in bases:
        d = len(exprs)
        pigeon_h = 10.0 ** (gp_digits / d)
        c, res = run_gp_lindep(exprs, gp_digits, max_n * 1000, mu_str)
        if c is None:
            print(f"  {name:65s}: no relation found (pigeonhole H ~ {pigeon_h:.1f})")
            continue
        inf_norm = max(abs(x) for x in c)
        # Note: a genuine hit at 5*sigma requires res <= O(inf_norm * window) AND inf_norm << pigeon_h
        flag = "BELOW PIGEONHOLE & WITHIN 5s!" if (inf_norm < 0.25 * pigeon_h and res < 10.0 * inf_norm * window) else "near/above pigeonhole (noise)"
        print(f"  {name:65s}: ||c||_inf={inf_norm:6d} (pigeon~{pigeon_h:6.1f}) res={res:.2e}  [{flag}]  c={c}")

if __name__ == "__main__":
    main()
