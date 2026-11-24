// feat2.hpp -- dense-indexed local features of core vertices, with the exact action of the generator on them.
//
// SPECIFICATION (brute force, `feat_bf` + `gen_bf`): F(u) is a sparse vector depending on the neighbourhood of
// radius <= 2 of core vertex u in the core tree with leaf counts.  For every u:
//     b_u = F(u)
//     a_u = sum_{events e} rate(e) [F_after_e(u) - F(u)]  +  (L_u/D_u) F_after(newborn child of u)
// with the elementary events of the continuous-time friend tree
//     G_v : v gains a leaf            rate g_v = L_v + sum_{core nbrs j} 1/D_j
//     K_v : v converts a leaf         rate k_v = L_v / D_v     (new core vertex w, (L_w,c_w) = (1,1), child of v;
//                                                               v -> (L_v - 1, c_v + 1))
// so that  d/dt E[sum_u theta.F(u)] = E[sum_u theta.a_u]  for every fixed theta.
// `gen_fast` computes the same a_u, b_u analytically (no event replay); test_feat2.cpp checks equality.
//
// Feature groups (all indices dense, layout in Layout):
//  A  own state exact:      L < Lfine:  1[(L, min(c,ccap))]
//  H  own state, hubs:      L >= Lfine: for cc = min(c,ccap), k = 0,1,2:  (Lfine/L)^k
//  N  neighbour, low deg:   own class s1 = (min(L,L1cap), min(c,c1cap)); nbr j with D_j < Dc:
//                           count of type (D_j, min(L_j,Lc), rel)      rel = 1 if j is the parent of u
//  B  neighbour, high deg:  nbr j with D_j >= Dc:  sum of (Dc/D_j)^k, k = 0,1,2, per (s1, rel)
//  N2,B2                    the same for L >= L1cap with the extra factors x and x^2 (two blocks), x = L1cap/L,
//                           per (min(c,c1cap), ...).  Options: topx=1 removes the unweighted N,B terms for
//                           L >= L1cap; Bconst=0 removes the k = 0 (constant) term of B and B2.  With both, the
//                           value of the neighbour features of a vertex tends to 0 as its own L -> infinity and
//                           the contribution of a neighbour tends to 0 as that neighbour's degree -> infinity.
//  P  pairs:                L <= L2cap, 2 <= c <= c2cap: number of unordered pairs of nbrs by coarse types
//                           ct(j) = (min(D_j,Dp), 1[L_j > 0])
//  Q  paths:                L <= L2cap, 1 <= c <= c2cap: number of paths u-j-k (k != u) with D_j < Dp,
//                           by (ct(j), min(D_k, Dk))
#pragma once
#include "lump.hpp"
#include <string>

struct F2Cfg {
    int Lfine = 12, ccap = 7;
    int L1cap = 4, c1cap = 4, Dc = 8, Lc = 2;
    int pairs = 1, paths = 1, L2cap = 2, c2cap = 3, Dp = 5, Dk = 5;
    int Rpow = 0, L3cap = 6, c3cap = 8, Rsum = 1;
    int Bconst = 1;                          // 0: no constant term per high-degree neighbour (features vanish as D_j -> inf)
    int topx = 0;                            // 1: for L >= L1cap neighbour features appear only with weights x, x^2 (x = L1cap/L)
    bool set(const std::string &k, int v) {
        if (k == "Lfine") Lfine = v; else if (k == "ccap") ccap = v; else if (k == "L1cap") L1cap = v; else if (k == "c1cap") c1cap = v;
        else if (k == "Dc") Dc = v; else if (k == "Lc") Lc = v; else if (k == "pairs") pairs = v; else if (k == "paths") paths = v;
        else if (k == "L2cap") L2cap = v; else if (k == "Rpow") Rpow = v; else if (k == "L3cap") L3cap = v; else if (k == "c3cap") c3cap = v; else if (k == "Rsum") Rsum = v; else if (k == "Bconst") Bconst = v; else if (k == "topx") topx = v; else if (k == "c2cap") c2cap = v; else if (k == "Dp") Dp = v; else if (k == "Dk") Dk = v;
        else return false;
        return true;
    }
};

struct Layout {
    F2Cfg G; int A0, H0, N0, B0, N20, B20, P0, Q0, R0, total; int NS1, NT, NC, ND, NR; int ex[128][3];
    explicit Layout(const F2Cfg &g) : G(g) {
        NS1 = (G.L1cap + 1) * (G.c1cap + 1); NT = (G.Dc - 2) * (G.Lc + 1) * 2; NC = 2 * (G.Dp - 1); ND = G.Dk - 1;
        A0 = 0; H0 = A0 + G.Lfine * (G.ccap + 1); N0 = H0 + 3 * (G.ccap + 1); B0 = N0 + NS1 * NT; N20 = B0 + NS1 * 6;
        B20 = N20 + 2 * (G.c1cap + 1) * NT; P0 = B20 + 2 * (G.c1cap + 1) * 6;
        Q0 = P0 + (G.pairs ? (G.L2cap + 1) * (G.c2cap - 1) * NC * NC : 0);
        R0 = Q0 + (G.paths ? (G.L2cap + 1) * G.c2cap * NC * ND : 0);
        // group R: all monomials of total degree 2..Rpow in the Rsum summaries x_m = sum_{nbrs j of class m} 1/D_j
        //   Rsum = 1: one class (x = cred);  Rsum = 2: classes L_j = 0 / L_j > 0;  Rsum = 3: L_j = 0 / 1 / >= 2
        NR = 0;
        if (G.Rpow >= 2) for (int e0 = 0; e0 <= G.Rpow; e0++) for (int e1 = 0; e1 <= (G.Rsum >= 2 ? G.Rpow : 0); e1++) for (int e2 = 0; e2 <= (G.Rsum >= 3 ? G.Rpow : 0); e2++) {
            int d = e0 + e1 + e2; if (d < 2 || d > G.Rpow) continue; ex[NR][0] = e0; ex[NR][1] = e1; ex[NR][2] = e2; NR++;
        }
        total = R0 + (G.L3cap + 1) * (G.c3cap + 1) * NR;
    }
    inline int ct(uint64_t D, uint64_t L) const { return (int)((std::min<uint64_t>(D, (uint64_t)G.Dp) - 2) * 2 + (L > 0 ? 1 : 0)); }
    inline int dk(uint64_t D) const { return (int)(std::min<uint64_t>(D, (uint64_t)G.Dk) - 2); }
    // ---- emitters: sink.add(index, value) ----
    template <class S> inline void own(uint64_t L, uint64_t c, double w, S &s) const {
        uint64_t cc = std::min<uint64_t>(c, (uint64_t)G.ccap);
        if (L < (uint64_t)G.Lfine) s.add(A0 + (int)L * (G.ccap + 1) + (int)cc, w);
        else { double x = (double)G.Lfine / (double)L; int b = H0 + 3 * (int)cc; s.add(b, w); s.add(b + 1, w * x); s.add(b + 2, w * x * x); }
    }
    // contribution of one neighbour (Dj, Lj, rel) to a vertex with own state (L, c)
    // top-class (L >= L1cap) terms: weights x and x^2 with x = L1cap/L, stored at N20/B20 + (xp*(c1cap+1) + c1)*...
    template <class S> inline void nbr(uint64_t L, uint64_t c, uint64_t Dj, uint64_t Lj, int rel, double w, S &s) const {
        int c1 = (int)std::min<uint64_t>(c, (uint64_t)G.c1cap);
        int s1 = (int)std::min<uint64_t>(L, (uint64_t)G.L1cap) * (G.c1cap + 1) + c1;
        bool top = L >= (uint64_t)G.L1cap; double x = top ? (double)G.L1cap / (double)L : 0.0;
        bool plainterm = !(top && G.topx);
        if (Dj < (uint64_t)G.Dc) {
            int tau = (int)((Dj - 2) * (G.Lc + 1) + std::min<uint64_t>(Lj, (uint64_t)G.Lc)) * 2 + rel;
            if (plainterm) s.add(N0 + s1 * NT + tau, w);
            if (top) { s.add(N20 + c1 * NT + tau, w * x); s.add(N20 + ((G.c1cap + 1) + c1) * NT + tau, w * x * x); }
        } else {
            double y = (double)G.Dc / (double)Dj; int b = B0 + (s1 * 2 + rel) * 3;
            if (plainterm) { if (G.Bconst) s.add(b, w); s.add(b + 1, w * y); s.add(b + 2, w * y * y); }
            if (top) for (int xp = 0; xp < 2; xp++) {
                double wx = xp == 0 ? w * x : w * x * x; int b2 = B20 + ((xp * (G.c1cap + 1) + c1) * 2 + rel) * 3;
                if (G.Bconst) s.add(b2, wx); s.add(b2 + 1, wx * y); s.add(b2 + 2, wx * y * y);
            }
        }
    }
    inline bool pair_ok(uint64_t L, uint64_t c) const { return G.pairs && L <= (uint64_t)G.L2cap && c >= 2 && c <= (uint64_t)G.c2cap; }
    inline bool path_ok(uint64_t L, uint64_t c) const { return G.paths && L <= (uint64_t)G.L2cap && c >= 1 && c <= (uint64_t)G.c2cap; }
    template <class S> inline void pair(uint64_t L, uint64_t c, int x, int y, double w, S &s) const {   // caller checks pair_ok
        if (x > y) std::swap(x, y);
        s.add(P0 + (((int)L * (G.c2cap - 1) + (int)c - 2) * NC + x) * NC + y, w);
    }
    template <class S> inline void path(uint64_t L, uint64_t c, int ctj, int dkk, double w, S &s) const { // caller checks path_ok
        s.add(Q0 + (((int)L * G.c2cap + (int)c - 1) * NC + ctj) * ND + dkk, w);
    }
    // group R emitters.  x[3] = the summaries (unused entries 0).
    inline bool r_ok(uint64_t L, uint64_t c) const { return NR > 0 && L <= (uint64_t)G.L3cap && c <= (uint64_t)G.c3cap; }
    inline int rcls(uint64_t Lj) const { return G.Rsum == 1 ? 0 : (G.Rsum == 2 ? (Lj > 0 ? 1 : 0) : (int)std::min<uint64_t>(Lj, 2)); }
    inline void mono(const long double *x, long double *out) const {
        long double pw[3][9];
        for (int m = 0; m < 3; m++) { pw[m][0] = 1.0L; for (int k = 1; k <= G.Rpow && k < 9; k++) pw[m][k] = pw[m][k - 1] * x[m]; }
        for (int i = 0; i < NR; i++) out[i] = pw[0][ex[i][0]] * pw[1][ex[i][1]] * pw[2][ex[i][2]];
    }
    template <class S> inline void rmono(uint64_t L, uint64_t c, const long double *x, double w, S &s) const {
        if (!r_ok(L, c)) return;
        long double v[128]; mono(x, v); int b = R0 + ((int)L * (G.c3cap + 1) + (int)c) * NR;
        for (int i = 0; i < NR; i++) s.add(b + i, (double)(w * v[i]));
    }
    // w * [mono(x1) - mono(x0)], difference formed in extended precision
    template <class S> inline void rmono_diff(uint64_t L, uint64_t c, const long double *x0, const long double *x1, double w, S &s) const {
        if (!r_ok(L, c)) return;
        long double v0[128], v1[128]; mono(x0, v0); mono(x1, v1); int b = R0 + ((int)L * (G.c3cap + 1) + (int)c) * NR;
        for (int i = 0; i < NR; i++) s.add(b + i, (double)((long double)w * (v1[i] - v0[i])));
    }
    // ---- exact differences for large arguments (avoid adding and subtracting huge equal terms) ----
    // own(L+1,c) - own(L,c), weight w, valid for L >= Lfine
    template <class S> inline void own_gain_diff_hub(uint64_t L, uint64_t c, double w, S &s) const {
        int b = H0 + 3 * (int)std::min<uint64_t>(c, (uint64_t)G.ccap); double Ld = (double)L, f = (double)G.Lfine;
        s.add(b + 1, -w * f / (Ld * (Ld + 1.0))); s.add(b + 2, -w * f * f * (2.0 * Ld + 1.0) / (Ld * Ld * (Ld + 1.0) * (Ld + 1.0)));
    }
    // nbr(L+1,c,j) - nbr(L,c,j), weight w, valid for L >= L1cap: only the x- and x^2-weighted terms change
    template <class S> inline void nbr_gain_diff_top(uint64_t L, uint64_t c, uint64_t Dj, uint64_t Lj, int rel, double w, S &s) const {
        int c1 = (int)std::min<uint64_t>(c, (uint64_t)G.c1cap); double Ld = (double)L, f = (double)G.L1cap;
        double wx[2] = {-w * f / (Ld * (Ld + 1.0)), -w * f * f * (2.0 * Ld + 1.0) / (Ld * Ld * (Ld + 1.0) * (Ld + 1.0))};
        for (int xp = 0; xp < 2; xp++) {
            if (Dj < (uint64_t)G.Dc) { int tau = (int)((Dj - 2) * (G.Lc + 1) + std::min<uint64_t>(Lj, (uint64_t)G.Lc)) * 2 + rel; s.add(N20 + (xp * (G.c1cap + 1) + c1) * NT + tau, wx[xp]); }
            else { double y = (double)G.Dc / (double)Dj; int b2 = B20 + ((xp * (G.c1cap + 1) + c1) * 2 + rel) * 3; if (G.Bconst) s.add(b2, wx[xp]); s.add(b2 + 1, wx[xp] * y); s.add(b2 + 2, wx[xp] * y * y); }
        }
    }
    // nbr(L,c; Dj+1) - nbr(L,c; Dj), weight w, valid for Dj >= Dc (leaf count of a high-degree nbr is not resolved)
    template <class S> inline void nbr_big_diff(uint64_t L, uint64_t c, uint64_t Dj, int rel, double w, S &s) const {
        int c1 = (int)std::min<uint64_t>(c, (uint64_t)G.c1cap);
        int s1 = (int)std::min<uint64_t>(L, (uint64_t)G.L1cap) * (G.c1cap + 1) + c1;
        double Dd = (double)Dj, dc = (double)G.Dc;
        double d1 = -w * dc / (Dd * (Dd + 1.0)), d2 = -w * dc * dc * (2.0 * Dd + 1.0) / (Dd * Dd * (Dd + 1.0) * (Dd + 1.0));
        bool top = L >= (uint64_t)G.L1cap;
        if (!(top && G.topx)) { int b = B0 + (s1 * 2 + rel) * 3; s.add(b + 1, d1); s.add(b + 2, d2); }
        if (top) { double x = (double)G.L1cap / (double)L; for (int xp = 0; xp < 2; xp++) { double xx = xp == 0 ? x : x * x; int b2 = B20 + ((xp * (G.c1cap + 1) + c1) * 2 + rel) * 3; s.add(b2 + 1, d1 * xx); s.add(b2 + 2, d2 * xx); } }
    }
};

struct VecSink { double *v; inline void add(int i, double x) { v[i] += x; } };

template <class Fn> static inline void for_nbrs(const Snap &T, uint32_t u, Fn fn) {
    if (T.parent[u] != NONE) fn(T.parent[u]);
    for (uint32_t k = T.first_child[u]; k != NONE; k = T.next_sib[k]) fn(k);
}

// event rates of vertex v: gain = rate of G_v, conv = rate of K_v, for the three models (see ftsim.hpp)
static inline void rates(const Snap &T, uint32_t v, int model, double &gain, double &conv) {
    double L = (double)T.L[v];
    if (model == 0) { gain = L + T.cred[v]; conv = L / (L + (double)T.c[v]); }
    else if (model == 3) { gain = L; conv = L / (L + (double)T.c[v]); }
    else { gain = L + ((model == 2 && T.L[v] == 0) ? 1.0 : 0.0); conv = T.L[v] > 0 ? 1.0 : 0.0; }
}
// ------------------------------------------------------------------ brute force (the specification)
template <class S> static void feat_bf(const Layout &Y, const Snap &T, uint32_t u, double w, S &s) {
    uint64_t L = T.L[u], c = T.c[u];
    Y.own(L, c, w, s);
    std::vector<uint32_t> nb; for_nbrs(T, u, [&](uint32_t j) { nb.push_back(j); });
    for (uint32_t j : nb) Y.nbr(L, c, T.L[j] + T.c[j], T.L[j], T.parent[u] == j ? 1 : 0, w, s);
    if (Y.r_ok(L, c)) { long double x[3] = {0, 0, 0}; for (uint32_t j : nb) x[Y.rcls(T.L[j])] += 1.0L / (long double)(T.L[j] + T.c[j]); Y.rmono(L, c, x, w, s); }
    if (Y.pair_ok(L, c))
        for (size_t i = 0; i < nb.size(); i++) for (size_t k = i + 1; k < nb.size(); k++)
            Y.pair(L, c, Y.ct(T.L[nb[i]] + T.c[nb[i]], T.L[nb[i]]), Y.ct(T.L[nb[k]] + T.c[nb[k]], T.L[nb[k]]), w, s);
    if (Y.path_ok(L, c))
        for (uint32_t j : nb) {
            uint64_t Dj = T.L[j] + T.c[j]; if (Dj >= (uint64_t)Y.G.Dp) continue;
            for_nbrs(T, j, [&](uint32_t k) { if (k != u) Y.path(L, c, Y.ct(Dj, T.L[j]), Y.dk(T.L[k] + T.c[k]), w, s); });
        }
}
// a_u and b_u by replaying every event within distance 2 of u.  T.cred must be current.
template <class S> static void gen_bf(const Layout &Y, Snap &T, uint32_t u, S &sa, S &sb, int model = 0) {
    feat_bf(Y, T, u, 1.0, sb);
    Scheme dummy; Lumper lp(dummy, T); std::vector<uint32_t> bl; lp.ball(u, 2, bl);
    for (uint32_t v : bl) {
        double gain, conv; rates(T, v, model, gain, conv);
        if (gain > 0) { T.L[v] += 1; feat_bf(Y, T, u, gain, sa); T.L[v] -= 1; feat_bf(Y, T, u, -gain, sa); }
        if (conv > 0) {
            uint32_t w = (uint32_t)T.L.size();
            T.L[v] -= 1; T.c[v] += 1; T.L.push_back(1); T.c.push_back(1); T.parent.push_back(v);
            T.first_child.push_back(NONE); T.next_sib.push_back(T.first_child[v]); T.first_child[v] = w;
            feat_bf(Y, T, u, conv, sa);
            if (v == u) feat_bf(Y, T, w, conv, sa);
            T.first_child[v] = T.next_sib[w];
            T.L.pop_back(); T.c.pop_back(); T.parent.pop_back(); T.first_child.pop_back(); T.next_sib.pop_back();
            T.L[v] += 1; T.c[v] -= 1;
            feat_bf(Y, T, u, -conv, sa);
        }
    }
}

// ------------------------------------------------------------------ theta-aware evaluation of group R
// Production only.  For a fixed weight vector theta the R-part of phi at own state s = (L,c) is the cubic polynomial
//   P_s(x) = sum_i theta[R0 + s*NR + i] * mono_i(x)   in the summaries x (Rpow <= 3 required).
// RFast accumulates theta.a (into A) and theta.b (into B) for group R directly, in double precision; differences for
// a small decrement of one summary use the exact third-order Taylor expansion (no cancellation).
// ftprod checks on request (rcheck=1) that this path reproduces the generic emitter path.
struct RFast {
    const Layout &Y; int nth, NS, NR, c3; std::vector<double> coef; double A[8], B[8];
    RFast(const Layout &y, int n, const double *const *th) : Y(y), nth(n) {
        NR = Y.NR; c3 = Y.G.c3cap + 1; NS = (Y.G.L3cap + 1) * c3; coef.resize((size_t)nth * NS * NR);
        for (int k = 0; k < nth; k++) for (int q = 0; q < NS; q++) for (int i = 0; i < NR; i++) coef[((size_t)k * NS + q) * NR + i] = th[k][Y.R0 + q * NR + i];
        reset();
    }
    void reset() { for (int k = 0; k < 8; k++) A[k] = B[k] = 0.0; }
    inline const double *cf(int k, uint64_t L, uint64_t c) const { return &coef[((size_t)k * NS + (size_t)L * c3 + (size_t)c) * NR]; }
    inline void mono(const double *x, double *v) const {
        double pw[3][4];
        for (int m = 0; m < 3; m++) { pw[m][0] = 1.0; pw[m][1] = x[m]; pw[m][2] = x[m] * x[m]; pw[m][3] = pw[m][2] * x[m]; }
        for (int i = 0; i < NR; i++) v[i] = pw[0][Y.ex[i][0]] * pw[1][Y.ex[i][1]] * pw[2][Y.ex[i][2]];
    }
    inline void addP(double *acc, uint64_t L, uint64_t c, const double *x, double w) const {
        if (!Y.r_ok(L, c)) return;
        double v[64]; mono(x, v);
        for (int k = 0; k < nth; k++) { const double *q = cf(k, L, c); double t = 0; for (int i = 0; i < NR; i++) t += q[i] * v[i]; acc[k] += w * t; }
    }
    // same with precomputed monomials v = mono(x)
    inline void addPv(double *acc, uint64_t L, uint64_t c, const double *v, double w) const {
        if (!Y.r_ok(L, c)) return;
        for (int k = 0; k < nth; k++) { const double *q = cf(k, L, c); double t = 0; for (int i = 0; i < NR; i++) t += q[i] * v[i]; acc[k] += w * t; }
    }
    // w * [P_s(x1) - P_s(x0)] with v0 = mono(x0) precomputed
    inline void addPdiffv(double *acc, uint64_t L, uint64_t c, const double *v0, const double *x1, double w) const {
        if (!Y.r_ok(L, c)) return;
        double v1[64]; mono(x1, v1);
        for (int k = 0; k < nth; k++) { const double *q = cf(k, L, c); double t = 0; for (int i = 0; i < NR; i++) t += q[i] * (v1[i] - v0[i]); acc[k] += w * t; }
    }
    // Taylor tables for class m at x: t1[k], t2[k], t3[k] = first/second/third pure derivatives of P_s in x_m
    inline void taylor_tab(uint64_t L, uint64_t c, const double *x, int m, double *t1, double *t2, double *t3) const {
        double pw[3][4];
        for (int a = 0; a < 3; a++) { pw[a][0] = 1.0; pw[a][1] = x[a]; pw[a][2] = x[a] * x[a]; pw[a][3] = pw[a][2] * x[a]; }
        int m1 = (m + 1) % 3, m2 = (m + 2) % 3;
        for (int k = 0; k < nth; k++) t1[k] = t2[k] = t3[k] = 0.0;
        for (int i = 0; i < NR; i++) {
            int e = Y.ex[i][m]; if (e == 0) continue;
            double other = pw[m1][Y.ex[i][m1]] * pw[m2][Y.ex[i][m2]];
            double a1 = e * pw[m][e - 1] * other, a2 = e >= 2 ? e * (e - 1) * pw[m][e - 2] * other : 0.0, a3 = e >= 3 ? 6.0 * other : 0.0;
            for (int k = 0; k < nth; k++) { double q = cf(k, L, c)[i]; t1[k] += q * a1; t2[k] += q * a2; t3[k] += q * a3; }
        }
    }
    inline void addTaylor(double *acc, const double *t1, const double *t2, const double *t3, double d, double w) const {
        for (int k = 0; k < nth; k++) acc[k] += w * d * (-t1[k] + d * (0.5 * t2[k] - d * t3[k] / 6.0));
    }
    // w * [P_s(x1) - P_s(x0)]
    inline void addPdiff(double *acc, uint64_t L, uint64_t c, const double *x0, const double *x1, double w) const {
        if (!Y.r_ok(L, c)) return;
        double v0[64], v1[64]; mono(x0, v0); mono(x1, v1);
        for (int k = 0; k < nth; k++) { const double *q = cf(k, L, c); double t = 0; for (int i = 0; i < NR; i++) t += q[i] * (v1[i] - v0[i]); acc[k] += w * t; }
    }
    // w * [P_s(x - d e_m) - P_s(x)] = w * ( -d P_m + d^2/2 P_mm - d^3/6 P_mmm )   (exact for cubic P)
    inline void addPtaylor(double *acc, uint64_t L, uint64_t c, const double *x, int m, double d, double w) const {
        if (!Y.r_ok(L, c)) return;
        double pw[3][4];
        for (int a = 0; a < 3; a++) { pw[a][0] = 1.0; pw[a][1] = x[a]; pw[a][2] = x[a] * x[a]; pw[a][3] = pw[a][2] * x[a]; }
        int m1 = (m + 1) % 3, m2 = (m + 2) % 3;
        double d1[64], d2[64], d3[64];
        for (int i = 0; i < NR; i++) {
            int e = Y.ex[i][m]; double other = pw[m1][Y.ex[i][m1]] * pw[m2][Y.ex[i][m2]];
            d1[i] = e >= 1 ? e * pw[m][e - 1] * other : 0.0;
            d2[i] = e >= 2 ? e * (e - 1) * pw[m][e - 2] * other : 0.0;
            d3[i] = e >= 3 ? 6.0 * other : 0.0;
        }
        for (int k = 0; k < nth; k++) {
            const double *q = cf(k, L, c); double t1 = 0, t2 = 0, t3 = 0;
            for (int i = 0; i < NR; i++) { t1 += q[i] * d1[i]; t2 += q[i] * d2[i]; t3 += q[i] * d3[i]; }
            acc[k] += w * d * (-t1 + d * (0.5 * t2 - d * t3 / 6.0));
        }
    }
};

// ------------------------------------------------------------------ analytic (fast)
struct NbInfo { uint32_t id; uint64_t D, L; int rel; double g, k; };
struct FastGen {
    const Layout &Y; const Snap &T; std::vector<double> gain, conv; std::vector<NbInfo> nb; std::vector<int> cts;
    FastGen(const Layout &y, const Snap &t, int model = 0) : Y(y), T(t) {
        size_t C = T.size(); gain.resize(C); conv.resize(C);
        for (size_t v = 0; v < C; v++) rates(T, (uint32_t)v, model, gain[v], conv[v]);
    }
    // pairs among the current neighbour types cts (plus optionally one extra type), for own state (L,c)
    template <class S> inline void pairs_all(uint64_t L, uint64_t c, int extra, double w, S &s) {
        if (!Y.pair_ok(L, c)) return;
        size_t n = cts.size();
        for (size_t i = 0; i < n; i++) for (size_t k = i + 1; k < n; k++) Y.pair(L, c, cts[i], cts[k], w, s);
        if (extra >= 0) for (size_t i = 0; i < n; i++) Y.pair(L, c, cts[i], extra, w, s);
    }
    template <class S> inline void paths_all(uint64_t L, uint64_t c, double w, S &s) {
        if (!Y.path_ok(L, c)) return;
        for (auto &j : nb) {
            if (j.D >= (uint64_t)Y.G.Dp) continue;
            int ctj = Y.ct(j.D, j.L);
            for_nbrs(T, j.id, [&](uint32_t k) { if (k != cur) Y.path(L, c, ctj, Y.dk(T.L[k] + T.c[k]), w, s); });
        }
    }
    uint32_t cur = NONE;
    RFast *rf = nullptr;      // if set, group R is evaluated by rf (theta-aware) instead of being emitted
    double rv0[64], rt1[3][8], rt2[3][8], rt3[3][8]; bool rtab[3] = {false, false, false}, rok = false;
    template <class S> void compute(uint32_t u, S &sa, S &sb) {
        cur = u;
        uint64_t L = T.L[u], c = T.c[u], Du = L + c; double g = gain[u], k = conv[u];
        nb.clear(); cts.clear();
        for_nbrs(T, u, [&](uint32_t j) { nb.push_back({j, T.L[j] + T.c[j], T.L[j], T.parent[u] == j ? 1 : 0, gain[j], conv[j]}); });
        bool anyP = Y.G.pairs && L <= (uint64_t)Y.G.L2cap + 1 && c <= (uint64_t)Y.G.c2cap;
        if (anyP || Y.G.pairs) for (auto &j : nb) cts.push_back(Y.ct(j.D, j.L));
        // ---- b_u ----
        Y.own(L, c, 1.0, sb);
        for (auto &j : nb) Y.nbr(L, c, j.D, j.L, j.rel, 1.0, sb);
        pairs_all(L, c, -1, 1.0, sb); paths_all(L, c, 1.0, sb);
        long double xs[3] = {0, 0, 0}; double xd[3] = {0, 0, 0};
        if (Y.NR > 0) {
            if (rf) {
                for (auto &j : nb) xd[Y.rcls(j.L)] += 1.0 / (double)j.D;
                rok = Y.r_ok(L, c) || Y.r_ok(L + 1, c);
                if (rok) rf->mono(xd, rv0);
                rf->addPv(rf->B, L, c, rv0, 1.0);
                rtab[0] = rtab[1] = rtab[2] = false;
            }
            else { for (auto &j : nb) xs[Y.rcls(j.L)] += 1.0L / (long double)j.D; Y.rmono(L, c, xs, 1.0, sb); }
        }
        // ---- G_u ----
        if (g > 0) {
            if (L >= (uint64_t)Y.G.Lfine) Y.own_gain_diff_hub(L, c, g, sa); else { Y.own(L + 1, c, g, sa); Y.own(L, c, -g, sa); }
            if (L >= (uint64_t)Y.G.L1cap) for (auto &j : nb) Y.nbr_gain_diff_top(L, c, j.D, j.L, j.rel, g, sa);
            else for (auto &j : nb) { Y.nbr(L + 1, c, j.D, j.L, j.rel, g, sa); Y.nbr(L, c, j.D, j.L, j.rel, -g, sa); }
            pairs_all(L + 1, c, -1, g, sa); pairs_all(L, c, -1, -g, sa);
            paths_all(L + 1, c, g, sa); paths_all(L, c, -g, sa);
            if (rf) { if (rok) { rf->addPv(rf->A, L + 1, c, rv0, g); rf->addPv(rf->A, L, c, rv0, -g); } }
            else { Y.rmono(L + 1, c, xs, g, sa); Y.rmono(L, c, xs, -g, sa); }
        }
        // ---- K_u (and the birth) ----
        if (k > 0) {
            Y.own(L - 1, c + 1, k, sa); Y.own(L, c, -k, sa);
            for (auto &j : nb) { Y.nbr(L - 1, c + 1, j.D, j.L, j.rel, k, sa); Y.nbr(L, c, j.D, j.L, j.rel, -k, sa); }
            Y.nbr(L - 1, c + 1, 2, 1, 0, k, sa);                                   // the new child seen from u
            pairs_all(L - 1, c + 1, Y.ct(2, 1), k, sa); pairs_all(L, c, -1, -k, sa);
            paths_all(L - 1, c + 1, k, sa); paths_all(L, c, -k, sa);              // the new child has no further nbrs
            // newborn w: own (1,1), single neighbour u (parent) in state (D_u, L-1); paths w-u-k for all old nbrs k of u
            if (Y.NR > 0 && rf) {
                double x2[3] = {xd[0], xd[1], xd[2]}; x2[Y.rcls(1)] += 0.5;
                rf->addP(rf->A, L - 1, c + 1, x2, k); if (Y.r_ok(L, c)) rf->addPv(rf->A, L, c, rv0, -k);
                double xw[3] = {0, 0, 0}; xw[Y.rcls(L - 1)] = 1.0 / (double)Du;
                rf->addP(rf->A, 1, 1, xw, k);
            } else if (Y.NR > 0) {
                long double x2[3] = {xs[0], xs[1], xs[2]}; x2[Y.rcls(1)] += 0.5L;               // new child (D=2, L=1)
                Y.rmono(L - 1, c + 1, x2, k, sa); Y.rmono(L, c, xs, -k, sa);
                long double xw[3] = {0, 0, 0}; xw[Y.rcls(L - 1)] = 1.0L / (long double)Du;     // newborn sees u = (D_u, L-1)
                Y.rmono(1, 1, xw, k, sa);
            }
            Y.own(1, 1, k, sa); Y.nbr(1, 1, Du, L - 1, 1, k, sa);
            if (Y.path_ok(1, 1) && Du < (uint64_t)Y.G.Dp) { int ctu = Y.ct(Du, L - 1); for (auto &j : nb) Y.path(1, 1, ctu, Y.dk(j.D), k, sa); }
        }
        // ---- neighbour events ----
        bool pk = Y.pair_ok(L, c), qk = Y.path_ok(L, c);
        for (size_t i = 0; i < nb.size(); i++) {
            auto &j = nb[i];
            if (j.g > 0) {                                                         // G_j: (D_j, L_j) -> (D_j+1, L_j+1)
                if (Y.r_ok(L, c) && rf) {
                    int c0 = Y.rcls(j.L), c1 = Y.rcls(j.L + 1);
                    if (c0 == c1) {
                        if (!rtab[c0]) { rf->taylor_tab(L, c, xd, c0, rt1[c0], rt2[c0], rt3[c0]); rtab[c0] = true; }
                        rf->addTaylor(rf->A, rt1[c0], rt2[c0], rt3[c0], 1.0 / ((double)j.D * ((double)j.D + 1.0)), j.g);
                    } else { double x2[3] = {xd[0], xd[1], xd[2]}; x2[c0] -= 1.0 / (double)j.D; x2[c1] += 1.0 / ((double)j.D + 1.0); rf->addPdiffv(rf->A, L, c, rv0, x2, j.g); }
                } else if (Y.r_ok(L, c)) {
                    long double x2[3] = {xs[0], xs[1], xs[2]}; int c0 = Y.rcls(j.L), c1 = Y.rcls(j.L + 1);
                    if (c0 == c1) x2[c0] -= 1.0L / ((long double)j.D * ((long double)j.D + 1.0L));
                    else { x2[c0] -= 1.0L / (long double)j.D; x2[c1] += 1.0L / ((long double)j.D + 1.0L); }
                    Y.rmono_diff(L, c, xs, x2, j.g, sa);
                }
                if (j.D >= (uint64_t)Y.G.Dc) Y.nbr_big_diff(L, c, j.D, j.rel, j.g, sa);
                else { Y.nbr(L, c, j.D + 1, j.L + 1, j.rel, j.g, sa); Y.nbr(L, c, j.D, j.L, j.rel, -j.g, sa); }
                if (pk) { int a = cts[i], b = Y.ct(j.D + 1, j.L + 1); if (a != b) for (size_t m = 0; m < nb.size(); m++) if (m != i) { Y.pair(L, c, b, cts[m], j.g, sa); Y.pair(L, c, a, cts[m], -j.g, sa); } }
                if (qk && j.D < (uint64_t)Y.G.Dp) {
                    int a = Y.ct(j.D, j.L); bool stay = j.D + 1 < (uint64_t)Y.G.Dp; int b = stay ? Y.ct(j.D + 1, j.L + 1) : -1;
                    for_nbrs(T, j.id, [&](uint32_t kk) { if (kk == u) return; int d = Y.dk(T.L[kk] + T.c[kk]); Y.path(L, c, a, d, -j.g, sa); if (stay) Y.path(L, c, b, d, j.g, sa); });
                }
            }
            if (j.k > 0) {                                                         // K_j: (D_j, L_j) -> (D_j, L_j-1), new child of j
                if (Y.r_ok(L, c)) {
                    int c0 = Y.rcls(j.L), c1 = Y.rcls(j.L - 1);
                    if (c0 != c1 && rf) { double x2[3] = {xd[0], xd[1], xd[2]}; x2[c0] -= 1.0 / (double)j.D; x2[c1] += 1.0 / (double)j.D; rf->addPdiffv(rf->A, L, c, rv0, x2, j.k); }
                    else if (c0 != c1) { long double x2[3] = {xs[0], xs[1], xs[2]}; x2[c0] -= 1.0L / (long double)j.D; x2[c1] += 1.0L / (long double)j.D; Y.rmono_diff(L, c, xs, x2, j.k, sa); }
                }
                if (j.D < (uint64_t)Y.G.Dc) { Y.nbr(L, c, j.D, j.L - 1, j.rel, j.k, sa); Y.nbr(L, c, j.D, j.L, j.rel, -j.k, sa); }
                if (pk) { int a = cts[i], b = Y.ct(j.D, j.L - 1); if (a != b) for (size_t m = 0; m < nb.size(); m++) if (m != i) { Y.pair(L, c, b, cts[m], j.k, sa); Y.pair(L, c, a, cts[m], -j.k, sa); } }
                if (qk && j.D < (uint64_t)Y.G.Dp) {
                    int a = Y.ct(j.D, j.L), b = Y.ct(j.D, j.L - 1);
                    if (a != b) for_nbrs(T, j.id, [&](uint32_t kk) { if (kk == u) return; int d = Y.dk(T.L[kk] + T.c[kk]); Y.path(L, c, a, d, -j.k, sa); Y.path(L, c, b, d, j.k, sa); });
                    Y.path(L, c, b, Y.dk(2), j.k, sa);                             // new far vertex of degree 2
                }
            }
            if (qk && j.D < (uint64_t)Y.G.Dp) {                                    // G_k for far vertices k
                int a = Y.ct(j.D, j.L);
                for_nbrs(T, j.id, [&](uint32_t kk) {
                    if (kk == u) return; uint64_t Dk = T.L[kk] + T.c[kk]; int d0 = Y.dk(Dk), d1 = Y.dk(Dk + 1);
                    if (d0 != d1 && gain[kk] > 0) { Y.path(L, c, a, d1, gain[kk], sa); Y.path(L, c, a, d0, -gain[kk], sa); }
                });
            }
        }
    }
};
