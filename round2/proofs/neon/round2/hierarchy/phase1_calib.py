#!/usr/bin/env python3
"""
Phase 1 Calibration & Exact Hierarchy Script:
1. Re-derive DAN, PAN, VAN(omega) to 30 digits via both GF/CF and hitting-time Laplace transform.
2. Solve Model IR-0 (undirected promotion rate L/(L+c), core-to-core credit rho = 0) via exact D-recursion and check closed form.
3. Solve 1-vertex mean-field closures:
   - MF-sym: rho_u = alpha * c_u (self-consistent alpha = E_edge[1/D])
   - MF-dir: rho_u = alpha_par * h_u + alpha_ch * k_u (parent/child self-consistent!)
   - MF-IR-eff: also compare with exact small-n and simulated values.
"""
import math
import mpmath as mp

mp.mp.dps = 40

def F_hkrz(mu):
    return mp.quad(lambda x: x**mu * mp.exp(x - 1), [0, 1])

print("=== 1. Calibration: HKRZ DAN / PAN / VAN(omega) ===")
mu_dan = mp.findroot(lambda m: 1 - m**2 - F_hkrz(m), mp.mpf('0.77'))
mu_pan = mp.findroot(lambda m: 1 - m - F_hkrz(m), mp.mpf('0.55'))
print("DAN mu =", mp.nstr(mu_dan, 30))
print("PAN mu =", mp.nstr(mu_pan, 30))

# Verify hitting-time Laplace transform phi_1(mu) via backward continued fraction
def phi1_cf(mu, N=200):
    # (l + 1 + mu) phi_l = l phi_{l+1} + phi_{l-1}, ratio r_l = phi_l / phi_{l-1}
    # r_l = 1 / (l + 1 + mu - l * r_{l+1})
    r = mp.mpf(0)
    for l in range(N, 0, -1):
        r = 1 / (l + 1 + mu - l * r)
    return r

print("Check |F(mu_pan) - phi1_cf(mu_pan)| =", mp.nstr(abs(F_hkrz(mu_pan) - phi1_cf(mu_pan)), 10))
print("Check |F(mu_dan) - phi1_cf(mu_dan)| =", mp.nstr(abs(F_hkrz(mu_dan) - phi1_cf(mu_dan)), 10))

print("\n=== 2. Model IR-0 (exact L/(L+c), rho = 0) ===")
# Exact forward recursion in D = 2, 3, ... for stationary unnormalized q_{L, c} with source = mu at (1,1)
# Wait: for a single newborn starting at (D,c)=(2,1), let M(mu) be expected discounted offspring.
# Since D only increases (2 -> 3 -> ...), backward DP from D_max down to 2 gives V(D, c) = E_{(D,c)}[sum e^{-mu tau_birth}]!
# Wait: let's compute BOTH forward sum_{L,c} (L/(L+c)) q_{L,c}(mu) / mu AND backward V(2,1)!
def ir0_forward_M(mu, Dmax=2000):
    # Normalize source to 1 at (L=1, c=1), so M(mu) = sum_{L,c} (L/(L+c)) q_{L,c}
    # Equation for D >= 2, L = D-1 down to 0:
    # For L = D-1 (c = 1):
    #   (mu + L + L/D) q_{L, 1} = (L-1) q_{L-1, 1} + delta_{D,2} * 1
    # For 0 <= L <= D-2 (c = D - L >= 2):
    #   (mu + L + L/D) q_{L, c} = (L-1) q_{L-1, c} + ((L+1)/D) q_{L+1, c-1}
    # Notice at each D, we compute L = D-1 down to 0 using previous degree D-1 (which has L-1, c)!
    prev = {} # map L -> q_{L, (D-1)-L}
    total_birth = mp.mpf(0)
    total_mass = mp.mpf(0)
    tail_last = mp.mpf(0)
    for D in range(2, Dmax + 1):
        curr = {}
        for L in range(D - 1, -1, -1):
            c = D - L
            rhs = mp.mpf(0)
            if D == 2 and L == 1:
                rhs += 1
            if L >= 1 and (L - 1) in prev:
                rhs += (L - 1) * prev[L - 1]
            if c >= 2:
                rhs += mp.mpf(L + 1) / D * curr[L + 1]
            denom = mu + L + mp.mpf(L) / D
            curr[L] = rhs / denom
        bD = sum(mp.mpf(L) / D * curr[L] for L in range(1, D))
        mD = sum(curr.values())
        total_birth += bD
        total_mass += mD
        tail_last = bD
        prev = curr
    # Wait: tail bD ~ C * D^{-1-mu}, so sum_{D>Dmax} bD ~ b_{Dmax} * Dmax / mu
    return total_birth, total_mass, tail_last

# Wait: for fast exact high-precision IR-0, let's check the asymptotic tail of bD!
# Even better: solve the ODE for h_D(u) or compute M(mu) with Richardson extrapolation in Dmax!
for Dmax in [100, 200, 400]:
    b, m, last = ir0_forward_M(mp.mpf('0.5'), Dmax)
    print(f"IR-0 at mu=0.5, Dmax={Dmax}: birth={mp.nstr(b,12)}, mass={mp.nstr(m,12)}, b_last={mp.nstr(last,8)}")

# Find mu_0 using Dmax = 200, 400, 800
def find_mu_ir0(Dmax):
    # Include asymptotic tail correction b_last * Dmax / mu (first-order)
    def f(mu):
        b, m, last = ir0_forward_M(mu, Dmax)
        # More accurate Euler-Maclaurin tail: b_D ~ A * D^{-(1+mu)}
        return b + last * Dmax / mu - 1
    return mp.findroot(f, mp.mpf('0.45'))

mu0_100 = find_mu_ir0(100)
mu0_200 = find_mu_ir0(200)
mu0_400 = find_mu_ir0(400)
print("IR-0 mu0 (Dmax=100, tail-corr):", mp.nstr(mu0_100, 15))
print("IR-0 mu0 (Dmax=200, tail-corr):", mp.nstr(mu0_200, 15))
print("IR-0 mu0 (Dmax=400, tail-corr):", mp.nstr(mu0_400, 15))
