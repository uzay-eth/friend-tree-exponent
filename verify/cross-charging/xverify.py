#!/usr/bin/env python3
"""Exact verification of a cross-charged drift certificate (proofs/xcharge.md).

Input: a pickle written by xcharge.py (float LP solution: values of f on all states, G, kappa, nu).
The values are converted to exact rationals (a float is an exact dyadic rational), clipped to the sign
constraints of the statement, and every inequality is then checked in exact rational arithmetic:

 upper mode (value lambda):  f >= 1 everywhere, G >= 1, kappa(TOP) >= 0, Kp := max(0, max_b kappa(b)), and
     for every box state and both endpoints of every 'L > Lmax' family
         main + pay * c * Kp + sum_j max(0, hi_j * e_j) [+ coarse: max(0, c/2 * e1, c*tau * e2)]  <=  lambda * f(x);
     beyond the box:  f(newborn of a TOP creator) <= lambda * G  and  Kp <= lambda * G.
   Conclusion:  E[C_n] <= f(T_3 centre) * prod_{m=3}^{n-1} (1 + lambda/m).
 lower mode (value delta):  f >= 0, kappa(b) >= 0, kappa(TOP) = 0, f = 0 beyond the box, and for every row
         main + sum_j min(0, hi_j * e_j) [+ coarse: min(0, c/2 * e1)]  >=  delta * f(x).
   Conclusion:  E[C_n] >= (f(T_3 centre)/sup f) * prod_{m=3}^{n-1} (1 + delta/m).

usage: xverify.py CERT.pkl VALUE [--expect-fail]       VALUE as a fraction, e.g. 577/1000
"""
import sys, pickle, time
from fractions import Fraction as Fr
from xcharge import Scheme


def load(path):
    d = pickle.load(open(path, 'rb'))
    S = Scheme(*d['params'])
    mode = d['mode']
    val = {k: Fr(v) for k, v in d['sol'].items()}
    return S, mode, val, d


def verify(S, mode, val, x, verbose=True):
    upper = (mode == 'upper')
    R = lambda p, q: Fr(p, q)
    TOP = S.TOP
    # sign constraints (clip)
    for k in list(val):
        if k == 'Kp':
            continue
        if k == 'G' or k[0] == 'f':
            val[k] = max(val[k], Fr(1) if upper else Fr(0))
        elif k[0] == 'k':
            if upper:
                if k[1] == TOP:
                    val[k] = max(val[k], Fr(0))
            else:
                val[k] = Fr(0) if k[1] == TOP else max(val[k], Fr(0))
    if upper:
        Kp = max([Fr(0)] + [val[('k', b)] for b in range(S.ncls)])
    else:
        val['G'] = Fr(0)
    def ev(d):
        return sum(co * val[k] for k, co in d.items())
    nchecks = 0; worst = None
    tau = Fr(1, S.top_min_deg)
    for name, r in S.all_rows(R):
        t = ev(r['main'])
        c = name[2] if name[0] == 'box' else name[1]
        if upper:
            if r['pay']:
                t += r['pay'] * c * Kp
            if r['top'] is not None:
                t += max(Fr(0), r['top'][0] * ev(r['top'][1]))
            for (a, hi, e) in r['nbr']:
                t += max(Fr(0), hi * ev(e))
            if r['coarse'] is not None:
                cc, e1, e2 = r['coarse']
                t += max(Fr(0), Fr(cc, 2) * ev(e1), cc * tau * ev(e2))
            slack = x * val[r['fx']] - t
        else:
            if r['top'] is not None:
                t += min(Fr(0), r['top'][0] * ev(r['top'][1]))
            for (a, hi, e) in r['nbr']:
                t += min(Fr(0), hi * ev(e))
            if r['coarse'] is not None:
                cc, e1, e2 = r['coarse']
                t += min(Fr(0), Fr(cc, 2) * ev(e1))
            slack = t - x * val[r['fx']]
        nchecks += 1
        if worst is None or slack < worst[0]:
            worst = (slack, name)
        if slack < 0:
            return False, nchecks, (float(slack), name)
    if upper:
        nb = val[('f', 1, 1, (TOP,))]
        for lhs in (nb, Kp):
            nchecks += 1
            if lhs > x * val['G']:
                return False, nchecks, ('beyond', float(lhs), float(x * val['G']))
    return True, nchecks, (float(worst[0]), worst[1])


def main():
    path = sys.argv[1]; x = Fr(sys.argv[2]); expect_fail = '--expect-fail' in sys.argv
    t0 = time.time()
    S, mode, val, d = load(path)
    ok, n, info = verify(S, mode, val, x)
    Dcap, c1, Lmax, cmax = d['params']
    fvals = [v for k, v in val.items() if k == 'G' or (k != 'Kp' and k[0] == 'f')]
    x0 = val[('f', 2, 0, ())]
    print(f"{mode} certificate {path}: Dcap={Dcap} c1={c1} Lmax={Lmax} cmax={cmax}, {len(fvals)} values of f")
    if not ok:
        print(f"  FAILS at {x} = {float(x):.6f}: {info}  (after {n} exact checks)")
        sys.exit(0 if expect_fail else 1)
    if expect_fail:
        print(f"  UNEXPECTED PASS at {x}"); sys.exit(1)
    from mpmath import mp, gamma, mpf
    mp.dps = 30
    xm = mpf(x.numerator) / x.denominator
    if mode == 'upper':
        const = mpf(x0.numerator) / x0.denominator * 2 / gamma(3 + xm)
        print(f"  PASS: all {n} exact inequalities hold at lambda = {x} = {float(x):.6f}; smallest slack {info[0]:.3e} at {info[1]}")
        print(f"  min f = {float(min(fvals)):.6f} (>= 1), f(T_3 centre) = {float(x0):.6f}")
        print(f"  => E[C_n] <= {float(x0):.6f} * prod_(m=3..n-1)(1 + lambda/m) = {float(const):.6f} * Gamma(n+lambda)/Gamma(n)   (n >= 3)")
    else:
        sup = max(fvals)
        pos = [v for v in fvals if v > 0]
        nbs = sorted({('f', 1, 1, (b,)) for b in range(S.ncls)})
        eta = min(val[k] for k in nbs)
        kap = x0 / sup
        const = mpf(kap.numerator) / kap.denominator * 2 / gamma(3 + xm)
        print(f"  PASS: all {n} exact inequalities hold at delta = {x} = {float(x):.6f}; smallest slack {info[0]:.3e} at {info[1]}")
        print(f"  sup f = {float(sup):.6f}, f(T_3 centre) = {float(x0):.6f}, smallest positive f = {float(min(pos)):.3e}, "
              f"min newborn value = {float(eta):.6f}")
        print(f"  => E[C_n] >= {float(kap):.6f} * prod_(m=3..n-1)(1 + delta/m) = {float(const):.6f} * Gamma(n+delta)/Gamma(n)   (n >= 3)")
    print(f"  ({time.time()-t0:.0f}s)")


if __name__ == '__main__':
    main()
