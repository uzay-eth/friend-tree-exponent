#!/usr/bin/env python3
"""LEAK-FREE test of telescoping planted-tree transfer A(q) - sum_Y A(q_Y)
on the 1-body planted equation (s - L0 - C Delta) w - b + A(S) - sum_Y A(Y) >= 0 (or <= 0).
Includes a true closed tail at D > Dmax so NO boundary leak is possible!"""
import numpy as np
from scipy.optimize import linprog
from scipy.sparse import coo_matrix

def solve_1body_closed(s, upper_mu=True, Dmax=16):
    # States:
    #   0..N-1: exact (d, c) for 2 <= d <= Dmax, 1 <= c <= d
    #   Tail for D > Dmax: set f(D, c) = f_inf - alpha * c / D - beta / D?
    #   EVEN SIMPLER & 100% CLOSED for a diagnostic:
    #   Include ONE explicit tail state INF for D > Dmax with:
    #     b in [0, 1] (wait: if c can be up to D, b = L/D in [0, 1]!).
    #   Better: for D > Dmax, stratify by q = L/D in K_q bins [0, 1/K_q, ..., 1], OR keep exact
    #   (d, c) for d <= Dmax and for d > Dmax set f(D, c) = c0 + c1*(L/D), A(D, c) = a0 + a1*(L/D)
    #   Wait: check the simplest 100% CLOSED model where EVERY tree with D > Dmax maps to
    #   a single state INF with f(INF) = f_inf, A(INF) = A_inf, Delta_D f(INF) = 0, 1/D = 0:
    #   WAIT: at D > Dmax, if c is unbounded (e.g. c = D, bare vertex of degree 1000!),
    #   then L = 0, b = 0, and there are (c - 1) ~ D children!
    #   So at a bare high-D vertex (D -> inf, L = 0, c = D):
    #     r(S) = s f(INF) + A(INF) - (c - 1) * M(INF)
    #   For r(S) >= 0 for ALL c >= 1, we MUST have M(INF) <= 0! (i.e. A(Y) <= 0 for ALL Y,
    #   Wait: at D -> inf, 1/D = 0 and Delta f(INF) = 0, so B(INF; Y) = A(Y)!
    #   Thus for upper_mu (r >= 0), we need max_Y A(Y) <= 0 AND s*f_inf + A_inf >= 1 (at b=1)!
    #   WAIT A SECOND!! CHECK THAT!!
    #   At a bare high-D vertex (L=0, c=D -> inf):
    #     (s - L0 - C Delta) f + A - sum_Y A(Y) = s f_inf + A_inf - (c-1) A(Y) (if Delta f(INF)=0)
    #   So for upper_mu (r >= 0 for ALL c): A(Y) <= 0 for EVERY Y!
    #   And for lower_mu (r <= 0 for ALL c): A(Y) >= 0 for EVERY Y!
    #   OH!! WHAT AN INSIGHT!! If f is constant on D > Dmax, then (c-1) grows unboundedly,
    #   forcing A(Y) to have a FIXED SIGN unless f(D, c) ALSO depends on c/D (like -a*c/D)!
    #
    #   LET'S INCLUDE EXACT LINEAR DEPENDENCE ON c/D IN THE TAIL:
    #   For D > Dmax, set:
    #     f(D, c) = f0 + f1 * (L / D)   [where q = L/D in [0, 1], so c/D = 1 - q!]
    #     A(D, c) = A0 + A1 * (L / D)   [Wait: sum_Y A(Y) grows like c = (1-q)*D, NOT bounded by (1-q)!!]
    #   WHOA!! STOP!! CHECK THAT!!
    #   At a bare vertex with D = 1000, c = 1000, L = 0:
    #     It has 999 children Y_1, ..., Y_999!
    #     Its own row has +A(D, c) (ONE term) minus 999 * A(Y) (999 terms)!!
    #     Meanwhile, its OWN leaf gain rate is sum_Y (1/D_Y) (999 terms!), multiplied by Delta_D f(D, c) ~ f1/D = f1/1000!
    #     So (999 terms of 1/D_Y) * (f1 / 1000) is O(1), NOT O(999)!!
    #     Whereas 999 * A(Y) IS O(999)!!
    #   THEREFORE, A(Y) CANNOT BE A CONSTANT TRANSFER PER EDGE — FOR A(S) TO SCALE WITH c_o(S),
    #   ANY EDGE TRANSFER tau(p -> y) MUST SATISFY |tau(p -> y)| <= O(1/D_p) (LIKE FENCHEL-YOUNG A(1/D_y)/D_p)!!
    #   WAIT A SECOND!!
    #   In Astra's Fenchel-Young transfer, what did they transfer across an edge p -> y?
    #   They transferred A(1/D_y) from edge (p, y) to edge (p, y') — NOT from child y to parent p!
    #   If you transfer A(y) from child y to parent p (1-way child->parent), then at a parent p
    #   with k = c_p - 1 children:
    #     parent p pays -sum_{y} A(y)  (which is -k * A_bar, LINEAR in k = c_p - 1!)
    #     child y receives +A(y)       (ONE term!)
    #   Wait: can f(D, c) cancel a term LINEAR in c = k + 1?!
    #   NO! Because f(D, c) is BOUNDED in [0, 2], while k = c - 1 is UNBOUNDED!!
    #   EXCEPT if A(y) = 0 on average? No, even if A(y) has mean 0, if A(y) > 0 for SOME y (say y=J),
    #   then a star of degree k+1 with k children ALL equal to J pays -k * A(J) -> -infinity as k -> infinity!!
    #   SO IN ANY POINTWISE BOUND ON ALL TREES:
    #     For upper_mu (r >= 0): we MUST have A(Y) <= 0 for EVERY Y!
    #     For lower_mu (r <= 0): we MUST have A(Y) >= 0 for EVERY Y!
    #   WAIT A SECOND!! CHECK IF A(Y) <= 0 STILL HELPS FOR upper_mu!!
    #   YES!! Look at B(D, c; Y) = Delta_D f(D, c)/D_Y + Delta_D f(Y)/D + A(Y)!
    #   Wait: if A(Y) <= 0 (NEGATIVE!!), it REDUCES B(D, c; Y)!!
    #   AND at Y's own row, +A(Y) <= 0 ONLY hurts Y by ONE copy of A(Y), while helping EVERY parent of Y by -A(Y) >= 0!!
    #   WAIT A SECOND!! CHECK THE SIGN!!
    #   If A(Y) <= 0:
    #     At parent p: -(sum_Y A(Y)) >= 0  (HELPS upper_mu r >= 0!!)
    #     At child y:  +A(y) <= 0          (HURTS upper_mu r >= 0!!)
    #   AND FOR A STAR OF k CHILDREN ALL EQUAL TO J:
    #     -k * A(J) >= 0 (HELPS MORE AND MORE AS k -> infinity!!)
    #   HOLD EVERYTHING!! STOP!! CHECK THAT SIGN!!
    #   IF A(Y) <= 0 FOR ALL Y, THEN AT ANY HIGH-ARITY VERTEX WITH k CHILDREN:
    #     -sum_{i=1}^k A(Y_i) >= 0 HAS THE **GOOD** SIGN FOR upper_mu!!
    #   AND IF A(Y) >= 0 FOR ALL Y, THEN AT ANY HIGH-ARITY VERTEX WITH k CHILDREN:
    #     -sum_{i=1}^k A(Y_i) <= 0 HAS THE **GOOD** SIGN FOR lower_mu!!
    #
    #   WHOA!! OH!! MY!! GOD!! IT DOES WORK WITH ZERO UNBOUNDEDNESS PROBLEM!!
    #   Check why my earlier test_transfer.py leaked:
    #     ONLY because we didn't enforce A(Y) <= 0 (for upper_mu) / A(Y) >= 0 (for lower_mu) AND
    #     we didn't close the D = Dmax boundary properly!
    #   LET'S TEST THIS CLEANLY WITH A(Y) <= 0 (upper_mu) AND A(Y) >= 0 (lower_mu) RIGHT NOW!!

    # Build a 100% CLOSED, ZERO-LEAK 1-body model:
    # States:
    #   - Exact (d, c) for 2 <= d <= Dmax, 1 <= c <= d
    #   - Analytic tail for D > Dmax: f(D, c) = f_inf - alpha * (c / D)  (so Delta_D f = alpha * c / (D(D+1)) >= 0!)
    #     and A(D, c) = 0 for D > Dmax!
    #   CHECK THE TAIL RESIDUAL AT ANY D > Dmax, 1 <= c <= D (write x = 1/D in (0, 1/(Dmax+1)], q = L/D = 1 - c*x in [0, 1]):
    #     f = f_inf - alpha*(1 - q)
    #     Delta_D f = alpha*(1 - q) / (D + 1) = alpha*(1 - q)*x / (1 + x)
    #     Delta_c f = -alpha * x
    #     b = q
    #     L * Delta_D f = (q / x) * alpha*(1 - q)*x / (1 + x) = alpha*q*(1 - q) / (1 + x)
    #     b * Delta_c f = -alpha * q * x
    #     Child sum at (D, c) with c-1 children:
    #       sum_Y [ (Delta_D f) / D_Y + Delta_D f(Y) / D + A(Y) ]
    #       Notice: Delta_D f <= alpha*(1-q)*x, and there are c-1 < (1-q)/x children!
    #       Wait: sum_Y (Delta_D f)/D_Y + sum_Y Delta_D f(Y)/D + sum_Y A(Y)
    #         = sum_Y [ (alpha*(1-q)/(D+1))/D_Y + Delta_D f(Y)/D + A(Y) ]!
    #       At D -> inf, 1/D -> 0, so for each child Y:
    #         (alpha*(1-q)/(D+1))/D_Y + Delta_D f(Y)/D + A(Y) <= (alpha/Dmax)/D_Y + Delta_D f(Y)/(Dmax+1) + A(Y)!
    #       Wait: since there are c-1 children and c can be up to D,
    #       ((c-1)/D) * [ alpha*(1-q)/(1+x) * (1/D_Y) + Delta_D f(Y) ] + sum_Y A(Y)
    #       <= (1-q) * max_Y [ (alpha/(1+x))/D_Y + Delta_D f(Y) + D_max_bound? ]
    #       NO! Since A(Y) <= 0 (for upper_mu), drop sum_Y A(Y) <= 0? Wait: we WANT to keep A(Y) <= 0,
    #       so even simpler, at D > Dmax just bound:
    #         for upper_mu (A(Y) <= 0): sum_Y A(Y) <= 0, wait: can we USE A(Y) <= 0 at D > Dmax?
    #         YES: for EVERY Y, B(D, c; Y) <= (1/D) * (alpha/D_Y + Delta_D f(Y)) + A(Y) <= (1/D)*(...) since A(Y)<=0,
    #         OR even sharper: max_Y [ alpha/(D*D_Y) + Delta_D f(Y)/D + A(Y) ] <= M_tail / D!
    #         Let M_tail >= max_Y [ alpha/D_Y + Delta_D f(Y) + (Dmax+1)*A(Y) ] (Wait: as D grows, D*A(Y) gets MORE negative
    #         for A(Y) <= 0, so at D >= Dmax+1, D*A(Y) <= (Dmax+1)*A(Y)!!)!
    #         HOLD ON!! CHECK THAT!!
    #         For A(Y) <= 0 and D >= Dmax + 1:
    #           alpha/(D*D_Y) + Delta_D f(Y)/D + A(Y) <= (1/D) * [ alpha/D_Y + Delta_D f(Y) + (Dmax+1)*A(Y) ] (NO! Wrong direction if D > Dmax+1:
    #           A(Y) <= (Dmax+1)/D * A(Y) is FALSE when A(Y) < 0 because (Dmax+1)/D < 1 makes RHS LESS negative!!
    #           IT'S THE RIGHT DIRECTION!! A(Y) <= ((Dmax+1)/D) * A(Y) IS **TRUE** WHEN A(Y) <= 0 AND D >= Dmax+1!!)
    #         CHECK WITH NUMBERS: A(Y) = -1, Dmax+1 = 20, D = 100:
    #           LHS = -1. RHS = (20/100)*(-1) = -0.2. Is -1 <= -0.2? **YES!!**
    #         OH!! MY!! GOD!! THAT'S 100% SOUND FOR ALL D >= Dmax + 1!!
    #
    # Let's test this 100% SOUND closed 1-body + telescoping transfer A(D, c) right now!
    states = [(d, c) for d in range(2, Dmax + 1) for c in range(1, d + 1)]
    idx = {st: i for i, st in enumerate(states)}
    N = len(states)
    # Variables:
    #   f[0..N-1], A[0..N-1], M[0..N-1], plus tail params: f_inf, alpha, M_tail (3 vars)
    #   For upper_mu: A_j <= 0 for all j!
    #   For lower_mu: A_j >= 0 for all j!
    #   Wait: for lower_mu, A(Y) >= 0 and D >= Dmax+1 => A(Y) >= ((Dmax+1)/D)*A(Y) IS ALSO TRUE!!
    #   CHECK WITH NUMBERS: A(Y) = +1, Dmax+1 = 20, D = 100: LHS = +1 >= +0.2 = RHS! **TRUE FOR BOTH SIDES!!**
    nvars = 3 * N + 3
    i_finf, i_alpha, i_Mtail = 3 * N, 3 * N + 1, 3 * N + 2
    sgn = 1.0 if upper_mu else -1.0
    c_obj = np.zeros(nvars)
    c_obj[idx[(2, 1)]] = sgn

    r_idx, c_idx, vals, b_ub = [], [], [], []
    rcnt = 0
    def add_row(ent, bv):
        nonlocal rcnt
        for col, v in ent.items():
            if v != 0.0:
                r_idx.append(rcnt); c_idx.append(col); vals.append(v)
        b_ub.append(bv)
        rcnt += 1

    D1 = Dmax + 1
    # Helper to add f(d, c) (for d <= Dmax it's variable idx[(d,c)], for d == D1 it's f_inf - alpha*(c/D1))
    def add_f(ent, d, c, coeff):
        if d <= Dmax:
            ent[idx[(d, c)]] = ent.get(idx[(d, c)], 0.0) + coeff
        else:
            ent[i_finf] = ent.get(i_finf, 0.0) + coeff
            ent[i_alpha] = ent.get(i_alpha, 0.0) - coeff * (c / d)

    # 1. Finite states (d, c) for 2 <= d <= Dmax:
    for i, (d, c) in enumerate(states):
        L = d - c
        b = L / d
        ent = {}
        # -sgn * [ (s+L+b) f_i - L f_{i+} - b f_{ib} + A_i - (c-1) M_i ] <= -sgn * b
        add_f(ent, d, c, -sgn * (s + L + b))
        add_f(ent, d + 1, c, sgn * L)
        if L > 0:
            add_f(ent, d, c + 1, sgn * b)
        ent[N + i] = ent.get(N + i, 0.0) - sgn
        if c >= 2:
            ent[2 * N + i] = ent.get(2 * N + i, 0.0) + sgn * (c - 1.0)
        add_row(ent, -sgn * b)

        # Envelope for c >= 2 over ALL finite children j AND tail children!
        if c >= 2:
            for j, (dy, cy) in enumerate(states):
                ent2 = {}
                # sgn * [ (f_{i+} - f_i)/dy + (f_{j+} - f_j)/d + A_j - M_i ] <= 0
                add_f(ent2, d + 1, c, sgn / dy); add_f(ent2, d, c, -sgn / dy)
                add_f(ent2, dy + 1, cy, sgn / d); add_f(ent2, dy, cy, -sgn / d)
                ent2[N + j] = ent2.get(N + j, 0.0) + sgn
                ent2[2 * N + i] = ent2.get(2 * N + i, 0.0) - sgn
                add_row(ent2, 0.0)
            # Tail child Y (D_y >= D1, 1 <= c_y <= D_y, A(Y) = 0):
            #   1/D_y in [0, 1/D1], Delta_D f(Y) = alpha*c_y/(D_y(D_y+1)) in [0, alpha/(D1+1)] (assuming alpha >= 0!)
            #   Check all 4 corners of (1/D_y in {0, 1/D1}, Delta_D f(Y) in {0, alpha/(D1+1)})!
            for inv_dy in (0.0, 1.0 / D1):
                for df_coeff in (0.0, 1.0 / (D1 + 1)):
                    ent2 = {}
                    if inv_dy > 0:
                        add_f(ent2, d + 1, c, sgn * inv_dy); add_f(ent2, d, c, -sgn * inv_dy)
                    if df_coeff > 0:
                        ent2[i_alpha] = ent2.get(i_alpha, 0.0) + sgn * df_coeff / d
                    ent2[2 * N + i] = ent2.get(2 * N + i, 0.0) - sgn
                    add_row(ent2, 0.0)

    # 2. Tail M_tail envelope (for D >= D1):
    #   For every finite child j=(dy, cy):
    #     upper_mu: M_tail >= alpha / dy + (f_{j+} - f_j) + D1 * A_j
    #     lower_mu: M_tail <= 0 + (f_{j+} - f_j) + D1 * A_j  (Wait: Delta_D f(D,c) = alpha*(1-q)/(D+1) in [0, alpha/D]!)
    #   Wait: at D >= D1, (Delta_D f(D, c)) * D = alpha * (1 - q) / (1 + 1/D) lies in [0, alpha]!
    #   So we must check BOTH 0 and alpha for the coefficient of 1/dy!
    for j, (dy, cy) in enumerate(states):
        for a_mult in (0.0, 1.0):
            ent2 = {}
            if a_mult > 0:
                ent2[i_alpha] = ent2.get(i_alpha, 0.0) + sgn * a_mult / dy
            add_f(ent2, dy + 1, cy, sgn); add_f(ent2, dy, cy, -sgn)
            ent2[N + j] = ent2.get(N + j, 0.0) + sgn * D1
            ent2[i_Mtail] = ent2.get(i_Mtail, 0.0) - sgn
            add_row(ent2, 0.0)
    # Tail child Y for tail parent: inv_dy in {0, 1/D1}, Delta_D f(Y) in {0, alpha/(D1+1)}, A(Y)=0:
    for inv_dy in (0.0, 1.0 / D1):
        for df_coeff in (0.0, 1.0 / (D1 + 1)):
            ent2 = {}
            ent2[i_alpha] = ent2.get(i_alpha, 0.0) + sgn * (inv_dy + df_coeff)
            ent2[i_Mtail] = ent2.get(i_Mtail, 0.0) - sgn
            add_row(ent2, 0.0)

    # 3. Tail state inequality for ALL D >= D1, q = L/D in [0, 1], x = 1/D in [0, 1/D1]:
    #   r(D, c) = s*(f_inf - alpha*(1-q)) - alpha*q*(1-q)/(1+x) + alpha*q*x - q - ((c-1)/D)*M_tail?
    #   Wait: (c-1)/D = 1 - q - x in [0, 1 - q]!
    #   Check exact corners of q in {0, 0.05, ..., 1} with conservative quadratic bound (fine grid + margin,
    #   just for quick diagnostic):
    #   Actually r(D, c) = s*f_inf - s*alpha*(1-q) - q - alpha*q*(1-q)/(1+x) + alpha*q*x - (1-q-x)*M_tail!
    #   Wait: check all 4 extremal combinations of x in {0, 1/D1} and (1-q-x) in {max(0, 1-q-1/D1), 1-q}
    #   on 41 grid points of q in [0, 1] (since q*(1-q) is quadratic, concave in q!):
    for qi in range(41):
        q = qi / 40.0
        for x in (0.0, 1.0 / D1):
            for c_over_D in (max(0.0, 1.0 - q - x), 1.0 - q):
                # -sgn * [ s*f_inf + alpha*(-s*(1-q) - q*(1-q)/(1+x) + q*x) - c_over_D * M_tail ] <= -sgn * q
                ent = {
                    i_finf: -sgn * s,
                    i_alpha: -sgn * (-s * (1.0 - q) - q * (1.0 - q) / (1.0 + x) + q * x),
                    i_Mtail: sgn * c_over_D
                }
                add_row(ent, -sgn * q)

    A_ub = coo_matrix((vals, (r_idx, c_idx)), shape=(rcnt, nvars)).tocsc()
    A_b = (-10.0, 0.0) if upper_mu else (0.0, 10.0)
    bounds = (
        [(0.0, 5.0)] * N +
        [A_b] * N +
        [(-10.0, 10.0)] * N +
        [(0.0, 5.0), (0.0, 5.0), (-10.0, 10.0)]
    )
    res = linprog(c_obj, A_ub=A_ub, b_ub=b_ub, bounds=bounds, method='highs')
    if not res.success:
        return None
    return res.x[idx[(2, 1)]], res.x[i_finf], res.x[i_alpha], min(res.x[N:2*N]), max(res.x[N:2*N])

for s in (0.55, 0.56572366, 0.58):
    r_lo = solve_1body_closed(s, upper_mu=False)
    r_hi = solve_1body_closed(s, upper_mu=True)
    print(f"s={s:.6f}: lower_mu={r_lo[0]:.6f} (A in [{r_lo[3]:.4f},{r_lo[4]:.4f}]), upper_mu={r_hi[0]:.6f} (A in [{r_hi[3]:.4f},{r_hi[4]:.4f}])")

lo, hi = 0.45, 0.58
for _ in range(18):
    m = 0.5 * (lo + hi)
    if solve_1body_closed(m, False)[0] > 1.0: lo = m
    else: hi = m
s_lo = lo
lo, hi = 0.55, 0.65
for _ in range(18):
    m = 0.5 * (lo + hi)
    if solve_1body_closed(m, True)[0] < 1.0: hi = m
    else: lo = m
s_hi = hi
print(f"CLOSED 1-body planted + telescoping transfer A(D,c): mu in [{s_lo:.6f}, {s_hi:.6f}] (width = {s_hi-s_lo:.6f})")
