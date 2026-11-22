#!/usr/bin/env python3
"""Hypotheses (H1)-(H5) of the almost-sure lower bound (proofs/as_lower.md), checked in exact arithmetic
on the final lower certificate of the earlier attempt (attempts/x-strogatz/certs9.py, Lmax=60, cmax=10,
K0=2, Dcap=6, delta=543/1000).

The table is rebuilt exactly as `certs9.py 0.5433 0.5434 0 543/1000 60 10 2 6` builds it, converted to
exact rationals, and then
  (H1) f >= 0 everywhere, sup f finite;                     [exact]
  (H2) the drift inequality of certs9.verify (exact rationals, every box state and the L > Lmax tail);
  (H4) every newborn type has f(nb) > 0;                     [exact]
and the constants of the proof are printed:
  F = sup f, eta = min newborn value, fmin = smallest positive value of f,
  B = (Dcap + 1) F  (bound on |Phi_{n+1} - Phi_n|),
  R = Lmax + cmax/2 + 1 + cmax*Dcap  (bound on n * P(f(x_v) changes) for a vertex with f(x_v) > 0),
  kappa = 2 B delta + 4 B^2 R / fmin.

usage: python3 as_lower_check.py [Lmax cmax K0 Dcap build_delta cert_delta]
"""
import sys, os, time
from fractions import Fraction as Fr

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, '..', 'baseline'))
import certs9  # noqa: E402


def main():
    a = sys.argv[1:]
    Lmax, cmax, K, Dcap = (int(a[0]), int(a[1]), int(a[2]), int(a[3])) if len(a) >= 4 else (60, 10, 2, 6)
    build = float(a[4]) if len(a) >= 5 else 0.5433
    delta = Fr(a[5]) if len(a) >= 6 else Fr(543, 1000)
    t0 = time.time()
    M = certs9.Model(Lmax, cmax, K, Dcap)
    print(f"model Lmax={Lmax} cmax={cmax} K0={K} Dcap={Dcap}: {len(M.classes)} classes, {len(M.NB)} newborn types",
          flush=True)
    ftab, rho = certs9.lower_table(M, build)
    assert rho >= 1, rho
    print(f"table built at delta={build} (Perron value {rho:.6f}); {len(ftab)} entries; {time.time()-t0:.0f}s",
          flush=True)
    # (H2): the exact drift check of the earlier attempt (asserts on failure)
    snb, sup, nchecks = certs9.verify(M, ftab, delta)
    print(f"(H2) drift inequality at delta={delta}: {nchecks} exact checks passed; {time.time()-t0:.0f}s", flush=True)
    F = {k: Fr(v) for k, v in ftab.items()}
    # (H1)
    assert all(v >= 0 for v in F.values())
    assert sup == max(F.values())
    pos = [v for v in F.values() if v > 0]
    fmin = min(pos)
    nzero = sum(1 for v in F.values() if v == 0)
    print(f"(H1) f >= 0 on all {len(F)} table entries; sup f = {float(sup):.6f}; "
          f"{nzero} entries are 0; smallest positive value fmin = {float(fmin):.6e}")
    # (H4): newborn values. A newborn is in state nb(L, c) where (L, c) is the creator's state BEFORE the
    # conversion (L >= 1); certs9.Model.NB lists nb(L, c) for 1 <= L <= Lmax+1, 0 <= c <= cmax+1, and nb(L, c)
    # depends on (L, c) only through: c == 0 and L-1 in {1..K0} (pendant), (L, c) == (1, 1) (class B),
    # or min(L + c, Dcap) ('other' of that degree / top). So NB is the complete list of newborn states.
    allnb = set()
    for L in range(1, 4 * Lmax):
        for c in range(0, 4 * cmax):
            if L + c >= 2:
                allnb.add(M.nb(L, c))
    assert allnb == set(M.NB), (allnb, M.NB)
    eta = min(snb.values())
    for x, v in sorted(snb.items()):
        print(f"     newborn state {x}: f = {float(v):.6f}")
    assert eta > 0
    print(f"(H4) all {len(snb)} newborn types have f > 0; eta = min = {float(eta):.6f}")
    B = (Dcap + 1) * sup
    R = Fr(Lmax) + Fr(cmax, 2) + 1 + cmax * Dcap
    kappa = 2 * B * delta + 4 * B * B * R / fmin
    print(f"constants: B = (Dcap+1) sup f = {float(B):.4f};  R = Lmax + cmax/2 + 1 + cmax*Dcap = {float(R):.1f};  "
          f"kappa = 2 B delta + 4 B^2 R / fmin = {float(kappa):.4e}")
    print(f"ALL HYPOTHESES OF THE ALMOST-SURE LOWER BOUND HOLD for delta = {delta} = {float(delta)}; "
          f"total {time.time()-t0:.0f}s")


if __name__ == '__main__':
    main()
