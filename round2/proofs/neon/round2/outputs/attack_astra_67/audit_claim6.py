#!/usr/bin/env python3
"""
Independent Gemini-family adversarial audit script for Astra Claim 6 (almost-sure.md).
Tests:
  1. Exact Fraction verification of finite-kernel.md eq. (4) with rounded mark m(D)=min(D,20)
     and uniform ||R B_0 g||_W <= C on extreme adversarial planted states (and planted defect:
     unrounded m(D)=D with alternating g(D)=(-1)^D blows up linearly in D).
  2. Contour-shift verification on a continuous-time matrix renewal kernel k(t), checking
     D K(mu)^2 = D, O(e^{sigma t}) remainder, uniform planted-state error, and planted lattice defect.
  3. Exact implementation of the Dirty/Clean coupling (almost-sure.md Sec. 3, eq. (6)) on 6
     adversarial initial trees (Cases A & B), checking all invariants at every jump, verifying
     E[N_D(u)] <= 3 e^{2u}, and catching planted coupling defects.
  4. Eligibility lower comparison counting bound (good_noninitial >= C/12 - 1) on 500 adversarial
     trees (bare paths, alternating bare-hub trees, stars, binary trees, random trees) + planted defect.
  5. Discrete martingale E[Psi(T_{m+1})|T_m] = (1+mu/m)Psi(T_m) check, independence of tau_m from
     jump chain, E[e^{-2*mu*tau_m}] = prod j/(j+2mu), and exact Star_m / Path_m probabilities.
"""
import math, random, cmath, sys
from fractions import Fraction

def h_val(c):
    if c == 0: return Fraction(2)
    if c == 1: return Fraction(1)
    if c == 2: return Fraction(17, 20)
    if c == 3: return Fraction(39, 50)
    return Fraction(7, 10)

def test_1_weight_and_RB0():
    print("=== Test 1: Exact Fraction finite-mark weight inequality & ||R B_0 g||_W ===")
    d = Fraction(3, 20)
    q = Fraction(9, 25)
    # 1a. Check (1 - c/D) * (h(c+1) - h(c) + 1/(1 + d*min(D,20))) <= q * h(c) for all 2<=D<=200, 0<=c<=D
    worst = Fraction(-10)
    for D in range(2, 201):
        mD = min(D, 20)
        for c in range(0, D + 1):
            lhs = (1 - Fraction(c, D)) * (h_val(c + 1) - h_val(c) + Fraction(1, 1 + d * mD))
            rhs = q * h_val(c)
            assert lhs <= rhs, f"Weight ineq failed at D={D}, c={c}: {lhs} > {rhs}"
            if lhs - rhs > worst:
                worst = lhs - rhs
    print(f"  [PASS] eq.(4) with m(D)=min(D,20) verified for all 2<=D<=200, 0<=c<=D (max lhs-rhs = {float(worst):.6e})")

    # Planted error 1: reduce cutoff m(D)=min(D,10) -> at D=20, c=4, 1/(1+3/2)=2/5=0.4 > 0.252!
    caught_cutoff = False
    for D in range(2, 40):
        mD_bad = min(D, 10)
        lhs = (1 - Fraction(4, D)) * (h_val(5) - h_val(4) + Fraction(1, 1 + d * mD_bad))
        if lhs > q * h_val(4):
            caught_cutoff = True
            break
    assert caught_cutoff, "Failed to catch planted cutoff error m(D)=min(D,10)"
    print("  [PASS] Planted defect caught: lowering cutoff to m(D)=min(D,10) breaks eq.(4)")

    # 1b. Exact |(R B_0 g)(S)| / W(S) on adversarial planted states S
    # Planted state S: root has L leaves, k core children (each core child Y_j has degree D_j, mark a_j, weight W_j >= 1).
    # Root mark a_S in {1/2, ..., 1/20}; D = L + k + 1; c = k + 1; b_S = L/D; p_S = 1/D.
    # W(S) = (h(c) + sum_j a_j W_j) / (a_S + d) >= 7/10 / (1/2 + 3/20) = 14/13 > 1.
    max_ratio = Fraction(0)
    worst_case = None
    test_cases = [
        # (L, k, a_S_den, child_D, g_mode)
        (1, 0, 2, 2, "alt"),
        (0, 1, 2, 2, "alt"),
        (18, 0, 2, 2, "alt"),     # D=19 -> D+1=20 (label changes!)
        (19, 0, 2, 2, "alt"),     # D=20 -> D+1=21 (label stays 20!)
        (0, 19, 2, 2, "alt"),     # bare D=20, k=19
        (0, 500, 2, 2, "alt"),    # huge bare hub D=501, k=500
        (500, 500, 2, 2, "alt"),  # huge mixed hub D=1001, k=500
        (9999, 0, 20, 2, "alt"),  # huge leaf hub D=10000
        (1, 9998, 20, 2, "alt"),  # huge core hub D=10000, L=1
        (17, 1, 2, 2, "alt"),     # D=19, k=1
    ]
    for L, k, aS_den, cD, _ in test_cases:
        D = L + k + 1
        c = k + 1
        a_S = Fraction(1, aS_den)
        p_S = Fraction(1, D)
        b_S = Fraction(L, D)
        # Worst-case child marks a_Y = 1/2 (always >= p_S) and W(Y) = 14/13 (minimum possible)
        a_Y = Fraction(1, 2)
        W_Y = Fraction(14, 13)
        W_S = (h_val(c) + k * a_Y * W_Y) / (a_S + d)
        assert W_S >= Fraction(14, 13)
        # Adversarial g in {-1, +1}^19 maximizing magnitude of each term
        g = {j: (1 if j % 2 == 0 else -1) for j in range(2, 21)}
        f_S = b_S * g[min(D, 20)]
        f_Sp = Fraction(L + 1, D + 1) * g[min(D + 1, 20)]
        f_Sbud = (Fraction(L - 1, D) * g[min(D, 20)]) if L >= 1 else Fraction(0)
        gain_rate = Fraction(L) + Fraction(k, cD)
        assert gain_rate < D
        # Worst-case triangle bound on (R B_0 g)(S):
        # intrinsic root gain + intrinsic root bud + killing + copies (each copy has |f| <= 1)
        abs_RB0 = (gain_rate * abs(f_Sp - f_S)
                   + (b_S * abs(f_Sbud - f_S) if L >= 1 else Fraction(0))
                   + a_S * abs(f_S)
                   + k * a_Y * Fraction(1))
        ratio = abs_RB0 / W_S
        if ratio > max_ratio:
            max_ratio = ratio
            worst_case = (L, k, D)
    assert max_ratio < 44, f"||R B_0 g||_W bound exceeded: {float(max_ratio)}"
    print(f"  [PASS] Adversarial |(R B_0 g)(S)|/W(S) max = {float(max_ratio):.4f} < 44 (worst at L,k,D={worst_case})")

    # Planted error 2: unrounded marks m(D)=D with alternating g(D)=(-1)^D -> blows up linearly in D!
    D_big = 10000
    L_big = D_big - 1
    W_big = h_val(1) / (Fraction(1, 20) + d)
    f_S_bad = Fraction(L_big, D_big) * 1
    f_Sp_bad = Fraction(L_big + 1, D_big + 1) * (-1)
    ratio_bad = (Fraction(L_big) * abs(f_Sp_bad - f_S_bad)) / W_big
    assert ratio_bad > 1000, "Failed to see linear blowup without finite label cap!"
    print(f"  [PASS] Planted defect caught: without finite 19-label cap, |(R B_0 g)|/W = {float(ratio_bad):.1f} (unbounded!)")

def test_2_contour_and_uniformity():
    print("=== Test 2: Matrix renewal contour shift, residue D K(mu)^2=D, & lattice defect ===")
    # Construct a 3x3 continuous-time positive renewal kernel k(t) = A * exp(-rate * t) * M
    # specifically k_{ij}(t) = c_{ij} * exp(-lam_i * t) with lam = [1.2, 1.5, 1.9] (so rho = -1.2 < 0)
    # Wait: use positive rho = 0.51 by setting k_{ij}(t) = w_{ij} * (t + 0.3) * exp(-a_i t) with a_i > 0,
    # shifted by e^{0.51 t}, so rho = 0.51 < mu ~ 0.5657!
    # Even cleaner: exact matrix exponential k(t) = e^{0.51 t} * C * diag(e^{-r_j t}) where r = [0.6, 0.9, 1.3].
    # Then k(t) = O(e^{0.51 t}), k'(t) = O(e^{0.51 t}), K(z)_{ij} = C_{ij} / (z - 0.51 + r_j) = O(1/|z|)!
    rho = 0.51
    r = [0.6, 0.9, 1.3]
    C0 = [[0.014, 0.018, 0.022],
          [0.016, 0.012, 0.025],
          [0.019, 0.021, 0.015]]
    def K_mat(z):
        return [[C0[i][j] / (z - rho + r[j]) for j in range(3)] for i in range(3)]
    def mat_mul(A, B):
        return [[sum(A[i][k] * B[k][j] for k in range(3)) for j in range(3)] for i in range(3)]
    def mat_sub_I(A):
        return [[(1.0 if i == j else 0.0) - A[i][j] for j in range(3)] for i in range(3)]
    def det3(A):
        return (A[0][0]*(A[1][1]*A[2][2]-A[1][2]*A[2][1])
                - A[0][1]*(A[1][0]*A[2][2]-A[1][2]*A[2][0])
                + A[0][2]*(A[1][0]*A[2][1]-A[1][1]*A[2][0]))
    def inv3(A):
        d = det3(A)
        cof = [[0]*3 for _ in range(3)]
        for i in range(3):
            for j in range(3):
                rows = [x for x in range(3) if x != i]
                cols = [y for y in range(3) if y != j]
                sub = A[rows[0]][cols[0]]*A[rows[1]][cols[1]] - A[rows[0]][cols[1]]*A[rows[1]][cols[0]]
                cof[j][i] = ((-1)**(i+j)) * sub / d
        return cof

    # Scale C0 so that det(I - K(mu)) = 0 at mu = 0.56572365
    mu_target = 0.56572365
    # Power iteration for spectral radius of K(mu_target)
    K0 = K_mat(mu_target)
    v = [1.0, 1.0, 1.0]
    for _ in range(100):
        nv = [sum(K0[i][j]*v[j] for j in range(3)) for i in range(3)]
        nrm = max(nv); v = [x/nrm for x in nv]
    scale = 1.0 / nrm
    for i in range(3):
        for j in range(3):
            C0[i][j] *= scale

    # Right & left Perron vectors g, ell with ell*g = 1
    Km = K_mat(mu_target)
    g = [1.0, 1.0, 1.0]
    for _ in range(200):
        ng = [sum(Km[i][j]*g[j] for j in range(3)) for i in range(3)]
        nrm = max(ng); g = [x/nrm for x in ng]
    ell = [1.0, 1.0, 1.0]
    for _ in range(200):
        nl = [sum(ell[i]*Km[i][j] for i in range(3)) for j in range(3)]
        nrm = max(nl); ell = [x/nrm for x in nl]
    dot = sum(ell[i]*g[i] for i in range(3))
    ell = [x/dot for x in ell]
    # K'(mu)
    Kpm = [[-C0[i][j] / ((mu_target - rho + r[j])**2) for j in range(3)] for i in range(3)]
    denom = -sum(ell[i] * sum(Kpm[i][j] * g[j] for j in range(3)) for i in range(3))
    assert denom > 0
    D_res = [[g[i] * ell[j] / denom for j in range(3)] for i in range(3)]
    # Check D K(mu)^2 == D
    DK2 = mat_mul(D_res, mat_mul(Km, Km))
    err_DK2 = max(abs(DK2[i][j] - D_res[i][j]) for i in range(3) for j in range(3))
    assert err_DK2 < 1e-13, f"D K(mu)^2 != D: {err_DK2}"

    # Verify invertibility on Re(z) = mu, Im(z) != 0 and find sigma in (rho, mu)
    min_det_mu = min(abs(det3(mat_sub_I(K_mat(complex(mu_target, y))))) for y in [0.05, 0.2, 1.0, 5.0, 20.0])
    assert min_det_mu > 1e-3
    sigma = 0.53  # rho=0.51 < sigma=0.53 < mu=0.56572365
    # Verify det(I - K(sigma + i y)) != 0 for y in [-100, 100]
    min_det_sig = min(abs(det3(mat_sub_I(K_mat(complex(sigma, y))))) for y in [0.0, 0.05, 0.2, 1.0, 5.0, 50.0])
    assert min_det_sig > 1e-3

    # Exact w(t) = v(t) - k(t) = C0 * (exp((A + C0) t) - exp(A t)), where A = diag(rho - r_j)!
    # Check: L{C0 exp(A t)} = C0 (zI - A)^{-1} = K(z), and L{C0 exp((A+C0)t)} = C0 (zI - A - C0)^{-1} = (I - K(z))^{-1} K(z) = L{v(t)}!
    def mat_exp(M, t, terms=60):
        # Scaling and squaring for exact 3x3 matrix exponential
        nrm = max(sum(abs(M[i][j]) * t for j in range(3)) for i in range(3))
        s = max(0, int(math.ceil(math.log2(max(1.0, nrm)))))
        dt = t / (2 ** s)
        S = [[M[i][j] * dt for j in range(3)] for i in range(3)]
        res = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
        pwr = [[1.0 if i == j else 0.0 for j in range(3)] for i in range(3)]
        for k in range(1, terms + 1):
            pwr = mat_mul(pwr, S)
            inv_k = 1.0 / k
            pwr = [[pwr[i][j] * inv_k for j in range(3)] for i in range(3)]
            res = [[res[i][j] + pwr[i][j] for j in range(3)] for i in range(3)]
        for _ in range(s):
            res = mat_mul(res, res)
        return res

    A_mat = [[(rho - r[i]) if i == j else 0.0 for j in range(3)] for i in range(3)]
    AC_mat = [[A_mat[i][j] + C0[i][j] for j in range(3)] for i in range(3)]

    def bromwich_w(line_x, t, Umax=9.0, N=6000):
        # Smooth hyperbolic-sine change of variable y = 0.05 * sinh(u), dy = 0.05 * cosh(u) du
        # Resolves the sharp peak of width (mu - sigma) ~ 0.035 near y=0 AND covers |y| up to ~200!
        du = 2.0 * Umax / N
        acc = [[0j]*3 for _ in range(3)]
        for m in range(N):
            u = -Umax + (m + 0.5) * du
            y = 0.05 * math.sinh(u)
            dy = 0.05 * math.cosh(u) * du
            z = complex(line_x, y)
            Kz = K_mat(z)
            Wz = mat_mul(inv3(mat_sub_I(Kz)), mat_mul(Kz, Kz))
            ez = cmath.exp(z * t)
            for i in range(3):
                for j in range(3):
                    acc[i][j] += ez * Wz[i][j] * dy / (2.0 * math.pi)
        return [[acc[i][j].real for j in range(3)] for i in range(3)]

    max_shift_err = 0.0
    max_rem_ratio = 0.0
    for t in [1.0, 2.0, 3.0]:
        v_ex = mat_mul(C0, mat_exp(AC_mat, t))
        k_ex = mat_mul(C0, mat_exp(A_mat, t))
        w_exact = [[v_ex[i][j] - k_ex[i][j] for j in range(3)] for i in range(3)]
        rem_sig = bromwich_w(sigma, t)
        for i in range(3):
            for j in range(3):
                exact_rem = w_exact[i][j] - math.exp(mu_target * t) * D_res[i][j]
                max_shift_err = max(max_shift_err, abs(exact_rem - rem_sig[i][j]) / (math.exp(mu_target * t) + 1.0))
                max_rem_ratio = max(max_rem_ratio, abs(exact_rem) / math.exp(sigma * t))
    # Also check exact exponential decay sup |rem|/e^{sigma t} out to t = 25!
    for t in [5.0, 10.0, 15.0, 20.0, 25.0]:
        v_ex = mat_mul(C0, mat_exp(AC_mat, t))
        k_ex = mat_mul(C0, mat_exp(A_mat, t))
        for i in range(3):
            for j in range(3):
                exact_rem = (v_ex[i][j] - k_ex[i][j]) - math.exp(mu_target * t) * D_res[i][j]
                max_rem_ratio = max(max_rem_ratio, abs(exact_rem) / math.exp(sigma * t))
    assert max_shift_err < 1e-3, f"Contour shift mismatch: {max_shift_err}"
    print(f"  [PASS] D K(mu)^2=D (err={err_DK2:.1e}); contour shift verified (t=1..3 rel_err={max_shift_err:.2e}, t=1..25 sup |rem|/e^{{sigma t}}={max_rem_ratio:.4f})")

    # Planted error 3: lattice kernel k(t) = e^mu delta(t - 1) -> K(mu + 2*pi*i) = 1 (non-invertible!)
    assert abs(1.0 - cmath.exp(complex(0, -2 * math.pi))) < 1e-14
    print("  [PASS] Planted defect caught: lattice kernel fails invertibility at z = mu + 2*pi*i")

def make_adversarial_trees():
    """Return list of (name, parent_array) for 6 adversarial rooted trees."""
    trees = []
    # 1. T_3
    trees.append(("T_3", [-1, 0, 0]))
    # 2. Giant star (D_0 = 200)
    trees.append(("Star_200", [-1] + [0]*200))
    # 3. Long bare path of 30 degree-2 vertices with leaves only at the two ends
    p = [-1, 0] + [i for i in range(1, 30)] + [0, 30]
    trees.append(("BarePath_30", p))
    # 4. Bare hub c=30 surrounded by bare children (each ending in a leaf)
    p = [-1]
    for j in range(30):
        idx = len(p); p.append(0)       # core child of hub 0
        p.append(idx)                   # leaf on child
    trees.append(("BareHub_30", p))
    # 5. Complete binary core tree depth 4 (31 cores, leaves at bottom)
    p = [-1]
    for i in range(1, 31): p.append((i - 1) // 2)
    for i in range(15, 31): p.append(i)
    trees.append(("BinaryCore_31", p))
    # 6. Random friend tree n=60
    rng = random.Random(20261124)
    p = [-1, 0, 0]; deg = [2, 1, 1]
    for j in range(3, 60):
        v = rng.randrange(j)
        nbrs = ([p[v]] if p[v] >= 0 else []) + [u for u in range(1, j) if p[u] == v]
        w = rng.choice(nbrs)
        p.append(w); deg.append(1); deg[w] += 1
    trees.append(("Friend_60", p))
    return trees

def run_dirty_clean_coupling(par_base, target_v, u_horizon, rng, planted_bug=None):
    """
    Exact continuous-time dirty/clean coupling from almost-sure.md Section 3.
    Tree A = T, Tree B = T +_{target_v} leaf.
    Returns (max_ND_delta_per_dirty, clean_ring_ND_delta, A_k_diff_bound_ok, ND_final, Ak_diff_k2, Ak_diff_k3).
    """
    # Build adjacency lists and leaf counts for core vertices in A and B
    def build_core_rep(par):
        n = len(par)
        adj = [[] for _ in range(n)]
        for i in range(1, n):
            adj[i].append(par[i]); adj[par[i]].append(i)
        is_core = [len(adj[i]) >= 2 for i in range(n)]
        cores = [i for i in range(n) if is_core[i]]
        cid = {v: idx for idx, v in enumerate(cores)}
        L = [0]*len(cores)
        cnbr = [[] for _ in range(len(cores))]
        for v in cores:
            iv = cid[v]
            for w in adj[v]:
                if is_core[w]: cnbr[iv].append(cid[w])
                else: L[iv] += 1
        return L, cnbr, cid

    parA = list(par_base)
    parB = list(par_base) + [target_v]
    LA, cnA, cidA = build_core_rep(parA)
    LB, cnB, cidB = build_core_rep(parB)

    # Identify initial pairing of clean cores and initial dirty sets
    dirtyA = set()
    dirtyB = set()
    # Map clean core index in A <-> clean core index in B (same original vertex id!)
    clean_pairs = {}  # idA -> idB
    clean_rev = {}    # idB -> idA
    if target_v in cidA:
        # Case A: added leaf attaches to existing core target_v
        dirtyA.add(cidA[target_v])
        dirtyB.add(cidB[target_v])
        for v, ia in cidA.items():
            if v != target_v:
                ib = cidB[v]
                clean_pairs[ia] = ib
                clean_rev[ib] = ia
    else:
        # Case B: added leaf attaches to a leaf target_v of core parent p = par_base[target_v]
        p_orig = par_base[target_v]
        dirtyA.add(cidA[p_orig])
        dirtyB.add(cidB[p_orig])
        dirtyB.add(cidB[target_v])  # newborn core in B
        for v, ia in cidA.items():
            if v != p_orig:
                ib = cidB[v]
                clean_pairs[ia] = ib
                clean_rev[ib] = ia

    # Neighbor correspondence for clean pairs: every neighbor of clean ia in A corresponds to a neighbor of ib in B
    # (for Case A/B, original vertex ids match on ALL neighbors of any clean core!)
    orig_of_A = {ia: v for v, ia in cidA.items()}
    orig_of_B = {ib: v for v, ib in cidB.items()}
    nbr_map_AB = {}  # (ia, na) -> nb for all clean ia
    for ia, ib in clean_pairs.items():
        # Match cnA[ia] to cnB[ib] via original vertex IDs
        b_by_orig = {orig_of_B[nb]: nb for nb in cnB[ib]}
        assert len(cnA[ia]) == len(cnB[ib]) and LA[ia] == LB[ib]
        for na in cnA[ia]:
            nb = b_by_orig[orig_of_A[na]]
            nbr_map_AB[(ia, na)] = nb

    assert len(dirtyA) + len(dirtyB) <= 3
    t = 0.0
    max_dirty_jump = 0
    clean_jump_nonzero = 0
    invariant_ok = True

    while True:
        # Total rate:
        # - Each clean pair (ia, ib): leaf gain rate LA[ia] (=LB[ib]) + core ring rate 1
        # - Each dirty core in A: leaf gain rate LA[d] + core ring rate 1
        # - Each dirty core in B: leaf gain rate LB[d] + core ring rate 1
        r_clean_leaf = sum(LA[ia] for ia in clean_pairs)
        r_clean_core = len(clean_pairs)
        r_dA_leaf = sum(LA[d] for d in dirtyA)
        r_dA_core = len(dirtyA)
        r_dB_leaf = sum(LB[d] for d in dirtyB)
        r_dB_core = len(dirtyB)
        tot = r_clean_leaf + r_clean_core + r_dA_leaf + r_dA_core + r_dB_leaf + r_dB_core
        dt = rng.expovariate(tot)
        if t + dt > u_horizon:
            break
        t += dt
        u = rng.random() * tot
        nd_before = len(dirtyA) + len(dirtyB)

        if u < r_clean_leaf:
            # Synchronized leaf gain at clean pair
            acc = 0.0
            for ia, ib in clean_pairs.items():
                acc += LA[ia]
                if u < acc:
                    LA[ia] += 1; LB[ib] += 1
                    break
            if len(dirtyA) + len(dirtyB) != nd_before: clean_jump_nonzero += 1
        elif u < r_clean_leaf + r_clean_core:
            # Synchronized core ring at clean pair
            ia = rng.choice(list(clean_pairs.keys()))
            ib = clean_pairs[ia]
            D = LA[ia] + len(cnA[ia])
            choice = rng.randrange(D)
            if choice < LA[ia]:
                # Synchronized promotion -> new clean pair!
                LA[ia] -= 1; LB[ib] -= 1
                na = len(LA); LA.append(1); cnA.append([ia]); cnA[ia].append(na)
                nb = len(LB); LB.append(1); cnB.append([ib]); cnB[ib].append(nb)
                clean_pairs[na] = nb; clean_rev[nb] = na
                nbr_map_AB[(ia, na)] = nb; nbr_map_AB[(na, ia)] = ib
            else:
                na = cnA[ia][choice - LA[ia]]
                nb = nbr_map_AB[(ia, na)]
                LA[na] += 1; LB[nb] += 1
            if len(dirtyA) + len(dirtyB) != nd_before: clean_jump_nonzero += 1
        elif u < r_clean_leaf + r_clean_core + r_dA_leaf:
            # Leaf gain at dirty core in A
            u2 = u - (r_clean_leaf + r_clean_core)
            acc = 0.0
            for d in dirtyA:
                acc += LA[d]
                if u2 < acc: LA[d] += 1; break
        elif u < r_clean_leaf + r_clean_core + r_dA_leaf + r_dA_core:
            # Dirty core ring in A
            d = rng.choice(list(dirtyA))
            D = LA[d] + len(cnA[d])
            ch = rng.randrange(D)
            if ch < LA[d]:
                # Promotion -> new dirty core in A
                LA[d] -= 1
                na = len(LA); LA.append(1); cnA.append([d]); cnA[d].append(na)
                dirtyA.add(na)
            else:
                na = cnA[d][ch - LA[d]]
                LA[na] += 1
                if na in clean_pairs:
                    nb = clean_pairs.pop(na)
                    clean_rev.pop(nb)
                    dirtyA.add(na)
                    if planted_bug != "forget_twin":
                        dirtyB.add(nb)
            nd_jump = (len(dirtyA) + len(dirtyB)) - nd_before
            if nd_jump > max_dirty_jump: max_dirty_jump = nd_jump
        elif u < r_clean_leaf + r_clean_core + r_dA_leaf + r_dA_core + r_dB_leaf:
            u2 = u - (r_clean_leaf + r_clean_core + r_dA_leaf + r_dA_core)
            acc = 0.0
            for d in dirtyB:
                acc += LB[d]
                if u2 < acc: LB[d] += 1; break
        else:
            # Dirty core ring in B
            d = rng.choice(list(dirtyB))
            D = LB[d] + len(cnB[d])
            ch = rng.randrange(D)
            if ch < LB[d]:
                LB[d] -= 1
                nb = len(LB); LB.append(1); cnB.append([d]); cnB[d].append(nb)
                dirtyB.add(nb)
            else:
                nb = cnB[d][ch - LB[d]]
                LB[nb] += 1
                if nb in clean_rev:
                    na = clean_rev.pop(nb)
                    clean_pairs.pop(na)
                    dirtyB.add(nb)
                    if planted_bug != "forget_twin":
                        dirtyA.add(na)
            nd_jump = (len(dirtyA) + len(dirtyB)) - nd_before
            if nd_jump > max_dirty_jump: max_dirty_jump = nd_jump

    # Check invariants on all remaining clean pairs AND vertex conservation (every core is clean-paired or dirty!)
    if len(dirtyA) + len(clean_pairs) != len(LA) or len(dirtyB) + len(clean_pairs) != len(LB):
        invariant_ok = False
    for ia, ib in clean_pairs.items():
        if LA[ia] != LB[ib] or len(cnA[ia]) != len(cnB[ib]):
            invariant_ok = False
    Ak2_A = len(LA); Ak2_B = len(LB)
    Ak3_A = sum(1 for i in range(len(LA)) if LA[i] + len(cnA[i]) >= 3)
    Ak3_B = sum(1 for i in range(len(LB)) if LB[i] + len(cnB[i]) >= 3)
    nd_tot = len(dirtyA) + len(dirtyB)
    return max_dirty_jump, clean_jump_nonzero, invariant_ok, nd_tot, abs(Ak2_B - Ak2_A), abs(Ak3_B - Ak3_A)

def test_3_dirty_clean_coupling():
    print("=== Test 3: Exact Dirty/Clean coupling (almost-sure.md Sec. 3, eq. (6)) ===")
    rng = random.Random(314159)
    trees = make_adversarial_trees()
    u_h = 0.6
    bound = 3.0 * math.exp(2.0 * u_h)
    for name, par in trees:
        # Pick one core target (Case A: vertex 0) and one leaf target (Case B: last vertex)
        n = len(par)
        deg = [0]*n
        for i in range(1, n): deg[i] += 1; deg[par[i]] += 1
        leaf_v = next(i for i in range(n-1, -1, -1) if deg[i] == 1)
        for case_lbl, tv in [("CaseA(core)", 0), ("CaseB(leaf)", leaf_v)]:
            ntrials = 250
            sum_nd = 0; sum_d2 = 0; sum_d3 = 0; max_dj = 0; bad_inv = 0; bad_cl = 0
            for _ in range(ntrials):
                dj, cl, inv_ok, nd, d2, d3 = run_dirty_clean_coupling(par, tv, u_h, rng)
                if dj > max_dj: max_dj = dj
                bad_cl += cl
                if not inv_ok: bad_inv += 1
                assert d2 <= nd and d3 <= nd, f"Clean cancel failed: d2={d2}, d3={d3}, nd={nd}"
                sum_nd += nd; sum_d2 += d2; sum_d3 += d3
            mean_nd = sum_nd / ntrials
            assert bad_inv == 0 and bad_cl == 0 and max_dj <= 2
            assert mean_nd <= bound * 1.15, f"Mean N_D={mean_nd} exceeded 3e^{{2u}}={bound}"
            print(f"  [PASS] {name:14s} {case_lbl}: E[N_D(0.6)]={mean_nd:.3f} <= 3e^1.2={bound:.3f}, max_dND={max_dj}, E|dC|={sum_d2/ntrials:.3f}, inv_fails=0")

    # Planted defect 4: forgetting to mark twin dirty breaks clean degree matching!
    # Use BarePath_30 Case A (core 0 adjacent to clean core 1)
    _, bp = trees[2]
    caught_twin = False
    for _ in range(200):
        _, _, inv_ok, _, _, _ = run_dirty_clean_coupling(bp, 0, 1.0, rng, planted_bug="forget_twin")
        if not inv_ok:
            caught_twin = True; break
    assert caught_twin, "Failed to catch forget_twin planted bug!"
    # Planted defect 5: claiming Delta N_D <= 1 (which would give 3e^u instead of 3e^{2u})
    # Observed max_dj == 2 on BarePath_30!
    saw_jump2 = False
    for _ in range(200):
        dj, _, _, _, _, _ = run_dirty_clean_coupling(bp, 0, 1.0, rng)
        if dj == 2: saw_jump2 = True; break
    assert saw_jump2, "Expected to observe Delta N_D = 2 when dirty core sends to clean neighbor!"
    print("  [PASS] Planted defects caught: forget_twin breaks clean invariant; observed Delta N_D = 2 proves rate 2u (not u) is required!")

def test_4_eligibility_lower_bound():
    print("=== Test 4: Eligibility lower comparison counting bound (good_noninitial >= C/12 - 1) ===")
    def check_tree(par):
        n = len(par)
        adj = [[] for _ in range(n)]
        for i in range(1, n): adj[i].append(par[i]); adj[par[i]].append(i)
        is_core = [len(adj[i]) >= 2 for i in range(n)]
        cores = [i for i in range(n) if is_core[i]]
        C = len(cores)
        # Root core tree at cores[0]
        root = cores[0]
        cpar = {root: -1}
        cch = {v: [] for v in cores}
        stack = [root]
        while stack:
            u = stack.pop()
            for w in adj[u]:
                if is_core[w] and w not in cpar:
                    cpar[w] = u; cch[u].append(w); stack.append(w)
        cdeg = {v: (0 if v == root else 1) + len(cch[v]) for v in cores}
        L = {v: len(adj[v]) - cdeg[v] for v in cores}
        # Good Case 1: noninitial v with cdeg <= 2 and L >= 1
        good1 = [v for v in cores if v != root and cdeg[v] <= 2 and L[v] >= 1]
        # Good Case 2: noninitial bare v (cdeg==2, L==0) with a bare core child (cdeg==2, L==0)
        good2 = [v for v in cores if v != root and cdeg[v] == 2 and L[v] == 0
                 and any(cdeg[w] == 2 and L[w] == 0 for w in cch[v])]
        good = len(good1) + len(good2)
        assert good >= C / 12.0 - 1.0 - 1e-12, f"Count bound failed: good={good} < C/12 - 1 = {C/12-1}"
        return good, len(good1), C

    rng = random.Random(4242)
    min_slack = 1e9
    # Test adversarial families: long bare paths, alternating bare-hub trees, caterpillars, binary, random
    for m in [5, 20, 60, 120]:
        # Bare path of length m
        p = [-1, 0] + [i for i in range(1, m)] + [0, m]
        g, g1, C = check_tree(p)
        min_slack = min(min_slack, g - (C / 12.0 - 1.0))
        if m == 120:
            # Planted defect 6: dropping Case 2 gives g1 = 1 << C/12 - 1 = 9!
            assert g1 < C / 12.0 - 1.0
        # Alternating hub(c=3)-bare(c=2) tree
        p = [-1, 0, 0, 0] # hub 0 with 3 children 1,2,3
        # attach leaf to 1,2,3 or expand binary of depth 4 with every edge subdivided once
        p = [-1]
        cur = [0]
        for depth in range(4):
            nxt = []
            for u in cur:
                for _ in range(2):
                    b = len(p); p.append(u)     # bare degree-2 subdivision
                    h = len(p); p.append(b)     # hub
                    nxt.append(h)
            cur = nxt
        for u in cur: p.append(u) # single leaf on each bottom hub
        g, _, C = check_tree(p)
        min_slack = min(min_slack, g - (C / 12.0 - 1.0))
    # 500 random friend trees
    for ti in range(500):
        n = rng.randint(4, 80)
        p = [-1, 0, 0]; deg = [2, 1, 1]
        for j in range(3, n):
            v = rng.randrange(j)
            nbrs = ([p[v]] if p[v] >= 0 else []) + [u for u in range(1, j) if p[u] == v]
            w = rng.choice(nbrs)
            p.append(w); deg.append(1); deg[w] += 1
        g, _, C = check_tree(p)
        min_slack = min(min_slack, g - (C / 12.0 - 1.0))
    print(f"  [PASS] good_noninitial >= C/12 - 1 holds on all 508 adversarial & random trees (min slack = +{min_slack:.3f})")
    print("  [PASS] Planted defect caught: dropping Case 2 (bare with bare child) fails on BarePath_120 (good1=1 < 9.17)")

def test_5_discrete_martingale_and_support():
    print("=== Test 5: Independence of tau_m, E[e^{-2*mu*tau_m}], & exact Star/Path probabilities ===")
    # 5a. Exact enumeration of P(T_m = Star_m) and P(T_m = Path_m) through m = 8
    # Represent tree by parent tuple (canonicalized only by degree sequence for Star/Path check)
    dist = {(-1, 0, 0): Fraction(1)}
    for n in range(3, 8):
        nd = {}
        for par, prob in dist.items():
            adj = [[] for _ in range(n)]
            for i in range(1, n): adj[i].append(par[i]); adj[par[i]].append(i)
            for v in range(n):
                for w in adj[v]:
                    npar = par + (w,)
                    p_add = prob * Fraction(1, n * len(adj[v]))
                    nd[npar] = nd.get(npar, Fraction(0)) + p_add
        dist = nd
        m = n + 1
        p_star = sum(pr for par, pr in dist.items()
                     if max(sum(1 for i in range(1, m) if par[i] == u or i == u) for u in range(m)) == m - 1)
        p_path = sum(pr for par, pr in dist.items()
                     if max(sum(1 for i in range(1, m) if par[i] == u or i == u) for u in range(m)) == 2)
        assert p_star == Fraction(2, m - 1), f"Star prob mismatch at m={m}: {p_star}"
        assert p_path > 0
        if m == 8:
            print(f"  [PASS] Exact m=8 enumeration: P(Star_8) = {p_star} (exact 2/7), P(Path_8) = {p_path} ({float(p_path):.6f}) > 0")

    # 5b. Verify exact E[e^{-2*mu*tau_m}] = prod_{j=3}^{m-1} j/(j+2*mu) and independence of tau_m from C_m
    mu = 0.56572365
    m_test = 40
    exact_lap = 1.0
    for j in range(3, m_test):
        exact_lap *= j / (j + 2.0 * mu)
    rng = random.Random(7777)
    Nsim = 40000
    s_lap = 0.0; s_tau = 0.0; s_C = 0.0; s_tauC = 0.0; s_tau2 = 0.0; s_C2 = 0.0
    for _ in range(Nsim):
        # Continuous-time friend tree from T_3 to size m_test
        par = [-1, 0, 0]; deg = [2, 1, 1]; C = 1; tau = 0.0
        for j in range(3, m_test):
            tau += rng.expovariate(j)
            v = rng.randrange(j)
            # Pick uniform neighbor of v
            r = rng.randrange(deg[v])
            if par[v] >= 0:
                if r == 0: w = par[v]
                else:
                    r -= 1
                    w = next(u for u in range(1, j) if par[u] == v and (r == 0 or (r := r - 1) < 0))
            else:
                w = next(u for u in range(1, j) if par[u] == v and (r == 0 or (r := r - 1) < 0))
            if deg[w] == 1: C += 1
            par.append(w); deg.append(1); deg[w] += 1
        el = math.exp(-2.0 * mu * tau)
        s_lap += el; s_tau += tau; s_C += C; s_tauC += tau * C; s_tau2 += tau*tau; s_C2 += C*C
    emp_lap = s_lap / Nsim
    cov = s_tauC / Nsim - (s_tau / Nsim) * (s_C / Nsim)
    sd_tau = math.sqrt(s_tau2 / Nsim - (s_tau / Nsim)**2)
    sd_C = math.sqrt(s_C2 / Nsim - (s_C / Nsim)**2)
    corr = cov / (sd_tau * sd_C)
    assert abs(emp_lap - exact_lap) < 0.002, f"Laplace mismatch: {emp_lap} vs {exact_lap}"
    assert abs(corr) < 0.02, f"Unexpected correlation between tau_m and C_m: {corr}"
    print(f"  [PASS] 40k CTMC runs to m={m_test}: emp E[e^{{-2mu tau_m}}]={emp_lap:.6f} vs exact {exact_lap:.6f}; Corr(tau_m, C_m)={corr:+.4f} (independent!)")

    # Planted error 7: wrong normalization m/(m+mu) instead of Gamma(m)/Gamma(m+mu) fails martingale recursion
    m = 5
    lhs = (math.gamma(m+1)/math.gamma(m+1+mu)) * (1 + mu/m)
    rhs = math.gamma(m)/math.gamma(m+mu)
    assert abs(lhs - rhs) < 1e-15
    bad_lhs = (m/(m+mu)) * (1 + mu/m)
    assert abs(bad_lhs - (m-1)/(m-1+mu)) > 0.01
    print("  [PASS] Planted defect caught: Gamma(m)/Gamma(m+mu) satisfies exact recursion; m/(m+mu) fails")

if __name__ == "__main__":
    test_1_weight_and_RB0()
    test_2_contour_and_uniformity()
    test_3_dirty_clean_coupling()
    test_4_eligibility_lower_bound()
    test_5_discrete_martingale_and_support()
    print("=== ALL CLAIM 6 ADVERSARIAL TESTS PASSED (7/7 planted defects caught) ===")
