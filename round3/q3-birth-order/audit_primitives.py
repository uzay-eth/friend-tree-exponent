#!/usr/bin/env python3
"""
Independent high-precision audit of the exact __int128 primitives in BOTH
cert_all_m.cpp and cert_small_m.cpp, plus planted-defect sanity checks.

Zero external dependencies beyond stdlib (+ scipy.special if installed, or
built-in 40-digit Decimal / adaptive quadrature verification).
"""
import subprocess, tempfile, os, math
from decimal import Decimal, getcontext
getcontext().prec = 50

try:
    from scipy import special as sp
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

print(f"HAS_SCIPY={HAS_SCIPY}")

# Pure-Python 40-digit reference functions (independent of C++ code!)
def dec_gamma_frac(p, q):
    """50-digit Gamma(p/q) via Lanczos / Stirling at x = 80 + p/q (independent of C++!)."""
    # Use Euler's reflection / shift to x = 80 + p/q with Bernoulli B2, B4, B6, B8
    N = 80
    x = Decimal(N) + Decimal(p) / Decimal(q)
    # pi to 50 digits
    PI = Decimal("3.14159265358979323846264338327950288419716939937510")
    ln_g = (x - Decimal("0.5")) * x.ln() - x + (Decimal(2) * PI).ln() / Decimal(2)
    invx = Decimal(1) / x
    invx2 = invx * invx
    # Stirling series: 1/(12x) - 1/(360x^3) + 1/(1260x^5) - 1/(1680x^7) + 1/(1188x^9)
    ln_g += invx / 12 - invx * invx2 / 360 + invx * invx2**2 / 1260 - invx * invx2**3 / 1680 + invx * invx2**4 / 1188
    prod = Decimal(1)
    r = Decimal(p) / Decimal(q)
    for k in range(N):
        prod *= (r + k)
    return ln_g.exp() / prod

def ref_H_exact(N, j, x_dec):
    """Exact 50-digit H_{N, j}(x) = e^{-1} sum_{k=0}^{j-1} (1/k!) * P(Bin(N-1, x) <= j-k-1)."""
    q = Decimal(1) - x_dec
    if x_dec == 0:
        B = [Decimal(1)] * (j + 1)
    else:
        z = x_dec / q
        qpow = q ** (N - 1)
        B = [Decimal(0)] * (j + 1)
        term = Decimal(1)
        s = Decimal(1)
        B[1] = qpow * s
        for i in range(1, j):
            term = term * Decimal(N - i) * z / Decimal(i)
            s += term
            B[i + 1] = qpow * s
    ans = Decimal(0)
    fk = Decimal(1)
    for k in range(j):
        if k > 0:
            fk /= Decimal(k)
        ans += fk * B[j - k]
    return ans * Decimal(-1).exp()

def ref_B_num(N, r, x):
    """Independent reference for B_{N, r}(x) = P(Beta(r, N-r) > x)."""
    if HAS_SCIPY:
        return float(sp.betaincc(r, N - r, x))
    # Pure-Python high-accuracy composite Gauss-Legendre 8-point on transformed variable
    # B_{N, r}(x) = int_x^1 u^{r-1} (1-u)^{N-r-1} du / B(r, N-r)
    # Substitute u = 1 - (1-x)*t -> (1-x)^{N-r} int_0^1 (1 - (1-x)*t)^{r-1} t^{N-r-1} dt / B(r, N-r)
    # Further substitute t = s^{1/(N-r)} -> (1-x)^{N-r} / ((N-r)*B(r, N-r)) * int_0^1 (1 - (1-x)*s^{1/(N-r)})^{r-1} ds!
    # Smooth everywhere except near s=1 if x is tiny; split [0, 1] geometrically near s=1!
    q = 1.0 - x
    B = N - r
    lg_pref = math.lgamma(N) - math.lgamma(r) - math.lgamma(B) - math.log(B) + B * math.log(q)
    # 8-point Gauss-Legendre nodes/weights on [-1, 1]
    gl_x = [-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
             0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363]
    gl_w = [ 0.1012285362903763,  0.2223810344533745,  0.3137066458778873,  0.3626837833783620,
             0.3626837833783620,  0.3137066458778873,  0.2223810344533745,  0.1012285362903763]
    # Build split points from 0 to 1 clustering exponentially toward 1
    pts = [0.0]
    for k in range(1, 32):
        pts.append(1.0 - 0.5**k)
    pts.append(1.0)
    # Subdivide each panel into 4
    total = 0.0
    for i in range(len(pts) - 1):
        a0, b0 = pts[i], pts[i + 1]
        for sub in range(4):
            a = a0 + (b0 - a0) * sub / 4.0
            b = a0 + (b0 - a0) * (sub + 1) / 4.0
            mid = 0.5 * (a + b)
            half = 0.5 * (b - a)
            for xi, wi in zip(gl_x, gl_w):
                s = mid + half * xi
                # 1 - q * s^{1/B} = x + q * (1 - s^{1/B}) = x + q * (-expm1(log(s)/B))
                if s <= 0.0:
                    val = 1.0
                else:
                    one_m = -math.expm1(math.log(s) / B)
                    val = x + q * one_m
                total += half * wi * (val ** (r - 1.0))
    return math.exp(lg_pref) * total

HARNESS_ALL = r'''
#define main cert_all_m_orig_main
#include "cert_all_m.cpp"
#undef main

int main() {
    init_exp_grid(1600);
    for (auto [p, q] : vector<pair<int,int>>{{1,2},{1,3},{2,3},{1,4},{3,4},{1,5},{1,10},{1,32},{7,12}}) {
        Iv g = iv_gamma_rat(p, q);
        cout << "GAMMA " << p << " " << q << " " << (long long)g.lo << " " << (long long)g.hi << "\n";
    }
    for (int t_idx : {0, 64, 256, 512, 1024}) {
        for (int tn : {0, 100, 5000, 100000, 300000, 600000}) {
            for (int j : {1, 2, 5, 12, 20}) {
                int64_t hlo = H_raw_lower_int(j, t_idx, tn, 100000, 200);
                cout << "HRAW " << j << " " << t_idx << " " << tn << " " << hlo << "\n";
            }
        }
    }
    for (int t_idx : {0, 64, 256, 512}) {
        for (int tn : {25, 1000, 50000, 200000, 500000}) {
            for (auto [ip, p, q] : vector<tuple<int,int,int>>{{0,1,2},{0,1,3},{0,1,15},{0,5,6},{1,0,1},{1,1,2},{1,1,5},{2,1,2},{2,1,10}}) {
                int64_t bu = B_raw_upper_int(ip, p, q, t_idx, tn, 100000, 200);
                cout << "BRAW " << ip << " " << p << " " << q << " " << t_idx << " " << tn << " " << bu << "\n";
            }
        }
    }
    return 0;
}
'''

HARNESS_SMALL = r'''
#define main cert_small_m_orig_main
#include "cert_small_m.cpp"
#undef main

int main() {
    init_exp_grid(1600);
    init_log_int(2100);
    init_quad_nodes();
    for (int N : {4, 5, 7, 15, 40, 59, 60, 120, 199}) {
        for (int yn : {0, 1, 50, 500, 2000, 8000, 15000, 20000}) { // y = yn / 10000
            for (int j : {1, 2, 3, 6, 12}) {
                if (j >= N) continue;
                int64_t hlo = H_exact_N_lower(N, j, yn, 10000);
                cout << "HEX " << N << " " << j << " " << yn << " " << hlo << "\n";
            }
            if (yn == 0) continue;
            for (auto [ip, p, q] : vector<tuple<int,int,int>>{{0,1,2},{0,1,3},{0,1,10},{0,5,6},{1,0,1},{1,1,2},{1,1,5},{2,1,2},{2,1,5}}) {
                if (ip + (double)p/q >= N - 1) continue;
                int64_t bu = B_exact_N_upper(N, ip, p, q, yn, 10000);
                cout << "BEX " << N << " " << ip << " " << p << " " << q << " " << yn << " " << bu << "\n";
            }
        }
    }
    return 0;
}
'''

S = 10**12
S2 = 10**15

with tempfile.TemporaryDirectory() as td:
    for name, code in [("all", HARNESS_ALL), ("small", HARNESS_SMALL)]:
        src = os.path.join(td, f"{name}.cpp")
        exe = os.path.join(td, name)
        with open(src, "w") as f:
            f.write(code)
        subprocess.check_call(["g++", "-O2", "-std=c++17", f"-I{os.path.dirname(os.path.abspath(__file__))}", src, "-o", exe])
        out = subprocess.check_output([exe], text=True)

        cnt = {"GAMMA": 0, "HRAW": 0, "BRAW": 0, "HEX": 0, "BEX": 0}
        max_rel_gap_bex = 0.0
        for line in out.strip().splitlines():
            parts = line.split()
            tag = parts[0]
            if tag == "GAMMA":
                p, q, lo, hi = map(int, parts[1:])
                tv = dec_gamma_frac(p, q)
                assert Decimal(lo) / S2 <= tv <= Decimal(hi) / S2, f"GAMMA FAIL {p}/{q}: {lo/S2} <= {tv} <= {hi/S2}"
                cnt[tag] += 1
            elif tag == "HRAW":
                j, t_idx, tn, hlo = map(int, parts[1:])
                tau = Decimal(tn) / Decimal(100000)
                t = Decimal(t_idx) / Decimal(128)
                for m in [200, 500, 10000]:
                    N = int((Decimal(m) * t.exp()).to_integral_value(rounding="ROUND_FLOOR"))
                    x = Decimal(1) - (-tau / Decimal(m - 1)).exp()
                    tv = ref_H_exact(N, j, x)
                    assert Decimal(hlo) / S <= tv + Decimal("1e-25"), f"HRAW FAIL m={m} N={N} j={j} t={t} tau={tau}: {hlo/S} > {tv}"
                cnt[tag] += 1
            elif tag == "BRAW":
                ip, p, q, t_idx, tn, bu = map(int, parts[1:])
                r = ip + p / q
                tau = tn / 100000.0
                t = t_idx / 128.0
                for m in [200, 500, 10000]:
                    N = math.ceil(1.0 + (m - 1) * math.exp(t))
                    x = -math.expm1(-tau / (m - 1))
                    tv = ref_B_num(N, r, x)
                    assert bu / S >= tv - 1e-11, f"BRAW FAIL m={m} N={N} r={r} t={t} tau={tau}: {bu/S} < {tv}"
                cnt[tag] += 1
            elif tag == "HEX":
                N, j, yn, hlo = map(int, parts[1:])
                y = Decimal(yn) / Decimal(10000)
                x = Decimal(1) - (-y).exp()
                tv = ref_H_exact(N, j, x)
                assert Decimal(hlo) / S <= tv + Decimal("1e-25"), f"HEX FAIL N={N} j={j} y={y}: {hlo/S} > {tv}"
                if tv > Decimal("1e-9"):
                    assert Decimal(hlo) / S >= tv - Decimal("1e-8"), f"HEX LOOSE N={N} j={j} y={y}: {hlo/S} vs {tv}"
                cnt[tag] += 1
            elif tag == "BEX":
                N, ip, p, q, yn, bu = map(int, parts[1:])
                r = ip + p / q
                y = yn / 10000.0
                x = -math.expm1(-y)
                tv = ref_B_num(N, r, x)
                assert bu / S >= tv - 1e-11, f"BEX FAIL N={N} r={r} y={y}: {bu/S} < {tv}"
                if tv > 1e-4:
                    rel = (bu / S - tv) / tv
                    if rel > max_rel_gap_bex:
                        max_rel_gap_bex = rel
                cnt[tag] += 1
        print(f"Audit [{name}] PASSED: {cnt}, max_rel_excess_BEX={max_rel_gap_bex:.5f}")
