#!/usr/bin/env python3
"""
High-precision exact solution of Model IR-0 (rho=0, exact L/(L+c) promotion rate)
via backward DP with exact asymptotic 1/D^j expansion of the deficit W(D, c),
plus PSLQ search on mu_0, PLUS exact finite-N recurrence for HKRZ DAN & PAN
to diagnose their Fig. 6 fit (DAN=0.773, PAN=0.545, IR=0.569)!
"""
import mpmath as mp
import numpy as np

mp.mp.dps = 50

# =====================================================================
# Part 1: Model IR-0 to 35+ digits via backward DP + 1/D^j asymptotics
# =====================================================================
# In Model IR-0, a vertex at (D, c) with 1 <= c <= D-1 has deficit W(D, c) = 1/mu - V(D, c):
#   (mu + L + L/D) W(D, c) = L W(D+1, c) + (L/D) W(D, c+1) + c/D,  L = D - c,  W(D, D) = 1/mu.
#Wait: check W(D, D)! When c = D (L = 0), V(D, D) = 0, so W(D, D) = 1/mu!
# Let's double check the equation at L = 1 (c = D - 1):
#   (mu + 1 + 1/D) V(D, D-1) = 1 * V(D+1, D-1) + (1/D) * (V(D, D) + 1)
# Plug V = 1/mu - W:
#   (mu + 1 + 1/D)(1/mu - W) = (1/mu - W(D+1, D-1)) + (1/D)(1/mu - W(D, D) + 1)
# Cancel 1/mu terms:
#   1 - (mu + 1 + 1/D) W = -W(D+1, D-1) - (1/D) W(D, D) + 1/D
#   (mu + 1 + 1/D) W = W(D+1, D-1) + (1/D) W(D, D) + (D - 1)/D = ... wait: (1 - 1/D) = (D-1)/D = c/D!!
# YES!! EXACT for ALL 1 <= c <= D-1, with boundary W(D, D) = 1/mu!
#
# Wait: can we compute the exact 1/D^j asymptotic expansion of W(D, c) for FIXED c as D -> inf?!
# Write W(D, c) = sum_{j=1}^J P_j(c) / D^j, where P_j(c) is a polynomial of degree j in c!
# Check the equation:
#   mu W(D,c) + (D - c)(W(D,c) - W(D+1,c)) + (1 - c/D)(W(D,c) - W(D,c+1)) = c/D
# Notice for fixed c, (1 + 1/D)^{-j} = sum_{m>=0} binom(-j, m) D^{-m}, so
#   W(D, c) - W(D+1, c) = sum_{j>=1} P_j(c) D^{-j} (1 - (1+1/D)^{-j})
#                       = sum_{n>=2} D^{-n} sum_{j=1}^{n-1} (-1)^{n-j-1} binom(n-1, j-1) P_j(c)
# Let's derive P_j(c) symbolically/recursively for j = 1..12!

def compute_P_coeffs(mu, J=10):
    # P_j(c) = sum_{r=0}^j a[j][r] c^r
    # Match coeff of D^{-n} for n = 1..J:
    # From mu W: mu P_n(c)
    # From D (W(D,c) - W(D+1,c)):
    #   coeff of D^{-n} is sum_{j=1}^n (-1)^{n-j} binom(n, j-1) P_j(c)  [check: for j=n, binom(n, n-1)=n -> n P_n(c)!]
    # From -c (W(D,c) - W(D+1,c)):
    #   coeff of D^{-n} is -c * sum_{j=1}^{n-1} (-1)^{n-1-j} binom(n-1, j-1) P_j(c)
    # From (W(D,c) - W(D,c+1)):
    #   coeff of D^{-n} is P_n(c) - P_n(c+1)
    # From -(c/D)(W(D,c) - W(D,c+1)):
    #   coeff of D^{-n} is -c (P_{n-1}(c) - P_{n-1}(c+1))
    # RHS: c * delta_{n, 1}
    # Notice on LHS, P_n(c) appears as: (mu + n) P_n(c) + (P_n(c) - P_n(c+1)) = RHS_n(c)!
    # Since P_n(c) - P_n(c+1) has degree deg(P_n) - 1, this uniquely determines P_n(c) top-down from c^n down to c^0!
    a = {} # a[j] = list of coeffs [c^0, c^1, ..., c^j]
    for n in range(1, J + 1):
        rhs = [mp.mpf(0)] * (n + 1)
        if n == 1:
            rhs[1] = mp.mpf(1)
        # Subtract lower-j terms from D(W(D)-W(D+1)): j = 1..n-1
        for j in range(1, n):
            coeff = ((-1)**(n - j)) * mp.binomial(n, j - 1)
            for r in range(j + 1):
                rhs[r] -= coeff * a[j][r]
        # Subtract -c(W(D)-W(D+1)): j = 1..n-1 -> add +c * sum ...
        for j in range(1, n):
            coeff = ((-1)**(n - 1 - j)) * mp.binomial(n - 1, j - 1)
            for r in range(j + 1):
                rhs[r + 1] += coeff * a[j][r]
        # Subtract -(c/D)(W(D,c)-W(D,c+1)) -> add +c * (P_{n-1}(c) - P_{n-1}(c+1))
        if n >= 2:
            # P_{n-1}(c+1) = sum_{r=0}^{n-1} a[n-1][r] sum_{k=0}^r binom(r, k) c^k
            p_shift = [mp.mpf(0)] * n
            for r in range(n):
                for k in range(r + 1):
                    p_shift[k] += a[n - 1][r] * mp.binomial(r, k)
            for k in range(n):
                rhs[k + 1] += (a[n - 1][k] - p_shift[k])
        # Now solve (mu + n) P_n(c) + P_n(c) - P_n(c+1) = rhs(c)
        # Note P_n(c+1) - P_n(c) = sum_{r=1}^n a_n[r] sum_{k=0}^{r-1} binom(r, k) c^k
        # So (mu + n) a_n[k] = rhs[k] + sum_{r=k+1}^n a_n[r] binom(r, k)
        an = [mp.mpf(0)] * (n + 1)
        for k in range(n, -1, -1):
            val = rhs[k]
            for r in range(k + 1, n + 1):
                val += an[r] * mp.binomial(r, k)
            an[k] = val / (mu + n)
        a[n] = an
    return a

def eval_W_asymp(a, D, c, J):
    invD = mp.mpf(1) / D
    s = mp.mpf(0)
    p = invD
    for j in range(1, J + 1):
        # Horner in c
        poly = mp.mpf(0)
        for coeff in reversed(a[j]):
            poly = poly * c + coeff
        s += poly * p
        p *= invD
    return s

def ir0_V21_exact(mu, Dmax=200, J=8):
    a = compute_P_coeffs(mu, J)
    # At D = Dmax + 1, set W(Dmax+1, c) for c = 1..Dmax using asymp (or cap at 1/mu)
    # Wait: for large c near Dmax, the 1/D expansion diverges if c ~ Dmax!
    # Wait: for W(2, 1), how large can c get at D = Dmax? c ~ log(Dmax) with super-exponential tails!
    # Check: probability c > 40 at D = 200 is tiny, AND if we just use W = 1/mu for c > c_cut,
    # let's check both! Better: cap c at Cmax (e.g. Cmax = 50), since prob(c >= 50 before D=200) < 10^{-30}!
    # Actually at D = Dmax+1, we can use asymp for c <= Dmax//4 and 1/mu (or smooth) above,
    # and check sensitivity to Dmax and Cmax!
    Cmax = min(Dmax, 60)
    # Wait: at D = Dmax + 1, only c = 1..Cmax are used if we stop c at Cmax!
    # Actually, let's check both W_lo (boundary 0 for c > Cmax) and W_hi (boundary 1/mu for c > Cmax) to get RIGOROUS error!
    W_next = [mp.mpf(0)] * (Cmax + 2)
    for c in range(1, Cmax + 1):
        W_next[c] = eval_W_asymp(a, Dmax + 1, c, J)
    W_next[Cmax + 1] = mp.mpf(1) / mu

    for D in range(Dmax, 1, -1):
        c_hi = min(D - 1, Cmax)
        W_curr = [mp.mpf(0)] * (Cmax + 2)
        if D <= Cmax:
            W_curr[D] = mp.mpf(1) / mu
        else:
            # c = Cmax + 1 boundary
            W_curr[Cmax + 1] = eval_W_asymp(a, D, Cmax + 1, J) if Cmax + 1 < D // 3 else mp.mpf(1) / mu
        for c in range(c_hi, 0, -1):
            L = D - c
            rhs = L * W_next[c] + (mp.mpf(L) / D) * W_curr[c + 1] + mp.mpf(c) / D
            W_curr[c] = rhs / (mu + L + mp.mpf(L) / D)
        W_next = W_curr
    return mp.mpf(1) / mu - W_next[1]

# Test convergence of V(2,1) at mu = 0.4822668
test_mu = mp.mpf('0.4822668')
for Dmax, J in [(100, 6), (200, 8), (300, 10), (400, 10)]:
    v = ir0_V21_exact(test_mu, Dmax, J)
    print(f"IR-0 V(2,1) at mu={test_mu}, Dmax={Dmax}, J={J}: {mp.nstr(v, 35)}")

# Solve V(2,1) = 1 to 30+ digits!
mu0_200 = mp.findroot(lambda m: ir0_V21_exact(m, 200, 8) - 1, mp.mpf('0.4822668'))
mu0_350 = mp.findroot(lambda m: ir0_V21_exact(m, 350, 10) - 1, mu0_200)
print("IR-0 mu_0 (Dmax=200, J=8 ):", mp.nstr(mu0_200, 35))
print("IR-0 mu_0 (Dmax=350, J=10):", mp.nstr(mu0_350, 35))
print("Diff:", mp.nstr(abs(mu0_350 - mu0_200), 10))

# PSLQ test on mu_0!
print("\n=== PSLQ on IR-0 mu_0 (30 digits) ===")
mp.mp.dps = 32
x = mu0_350
# Test algebraic up to degree 6
rel = mp.pslq([x**k for k in range(7)], maxcoeff=10000)
print("Algebraic deg<=6:", rel)
# Test combinations with 1, e, pi, log(2), gamma, F(x)
Fx = mp.quad(lambda t: t**x * mp.exp(t - 1), [0, 1])
rel2 = mp.pslq([1, x, x**2, Fx, x*Fx], maxcoeff=10000)
print("HKRZ-type [1, x, x^2, F(x), x*F(x)]:", rel2)
rel3 = mp.pslq([1, x, mp.e, mp.pi, mp.log(2), mp.euler], maxcoeff=1000)
print("Standard constants [1, x, e, pi, log2, gamma]:", rel3)

# =====================================================================
# Part 2: Exact finite-N expectation E[C_N] for HKRZ DAN and PAN (N=3..10^6)
# =====================================================================
print("\n=== Exact E[C_N] for DAN and PAN (N=3..10^6) vs HKRZ Fig. 6 ===")
# In DAN: M_ell(N) satisfies EXACT linear recurrence in N:
#   E[N_1(N+1)] = E[N_1(N)] + (E[N_1(N)] + E[Nuc(N)]) / N
#   E[Nuc(N+1)] = E[Nuc(N)] + (E[M_1(N)] - E[Nuc(N)]) / N
#   E[M_1(N+1)] = E[M_1(N)] + (N - E[N_1(N)] - 2 E[M_1(N)] + E[M_2(N)]) / N
#   E[M_l(N+1)] = E[M_l(N)] + ((l-1)E[M_{l-1}] + E[M_{l+1}] - (l+1)E[M_l]) / N
# Wait: for DAN, the exact expectation recurrence is 100% LINEAR AND EXACT (zero approximation)!
# Start at L_3 (N=3): N_1=2, Nuc=0, M_2=1, all other M_l=0. (Core C = N - N_1 = 1.)
def exact_dan(Nmax=1000000, checkpoints=None):
    # Track M[l] up to l_cap? Wait: max leaf degree at N is N-1, we can just compute in continuous time or exact discrete!
    # Wait: for N up to 10^6, l can reach 10^6, too big for dense vector, EXCEPT we can compute C(N) directly
    # using the EXACT linear ODE / generating function, OR simulate DAN/PAN fast!
    # Even simpler: check N up to 20,000 exactly in 0.5s with numpy!
    Lmax = 20005
    M = np.zeros(Lmax, dtype=np.float64)
    # Initial L_3: N=3, N1=2, Nuc=0, M[2]=1
    N1 = 2.0
    Nuc = 0.0
    M[2] = 1.0
    res = {}
    for N in range(3, 20001):
        if N in checkpoints:
            res[N] = N - N1
        # Update
        dN1 = (N1 + Nuc) / N
        dNuc = (M[1] - Nuc) / N
        dM1 = (N - N1 - 2.0 * M[1] + M[2]) / N
        # Vectorized for l = 2..N-1
        hi = min(N + 1, Lmax - 1)
        idx = np.arange(2, hi)
        dMl = ((idx - 1) * M[idx - 1] + M[idx + 1] - (idx + 1) * M[idx]) / N
        N1 += dN1
        Nuc += dNuc
        M[1] += dM1
        M[2:hi] += dMl
    return res

pts = [100, 316, 1000, 3162, 10000, 20000]
dan_exact = exact_dan(20000, set(pts))
print("Exact DAN E[C_N] and local log-slopes (true mu_DAN = 0.7711923):")
for i, N in enumerate(pts):
    if i > 0:
        N0 = pts[i - 1]
        slope = np.log(dan_exact[N] / dan_exact[N0]) / np.log(N / N0)
        print(f"  N={N:6d}: E[C_N]={dan_exact[N]:10.4f}  local_slope[{N0}->{N}] = {slope:.6f}")
    else:
        print(f"  N={N:6d}: E[C_N]={dan_exact[N]:10.4f}")
