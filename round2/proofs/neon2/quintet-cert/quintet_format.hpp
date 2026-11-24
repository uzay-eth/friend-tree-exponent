#pragma once
#include <algorithm>
#include <array>
#include <cassert>
#include <climits>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <utility>
#include <vector>

#if defined(__BYTE_ORDER__) && defined(__ORDER_LITTLE_ENDIAN__)
static_assert(__BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__, "FT2_U5V1 requires a little-endian host");
#endif

__extension__ typedef __int128 q5_i128;
__extension__ typedef unsigned __int128 q5_u128;

/*
 * Frozen Binary Format & Exact Deterministic Indexing for Neon-2 Unoriented 5-Vertex
 * (Diameter-4) Certificate ("FT2_U5V1").
 *
 * ZERO floating-point arithmetic is used anywhere in this header. Both the untrusted
 * builder (build_quintet.cpp) and the standalone zero-float verifier (verify_quintet.cpp)
 * share ONLY this frozen layout specification.
 */

static constexpr uint64_t Q5_MAGIC = 0x4654325F55355631ULL; // "FT2_U5V1"
static constexpr int64_t  Q5_WEIGHT_SCALE = INT64_C(1) << 40;
static constexpr int64_t  Q5_CENTER_SCALE = INT64_C(1) << 30;
static constexpr int64_t  Q5_MAX_W = 16 * Q5_WEIGHT_SCALE;
static constexpr int64_t  Q5_MAX_E = 32 * Q5_WEIGHT_SCALE;
static constexpr int64_t  Q5_MAX_Q = INT64_C(1) << 40;
static constexpr uint64_t Q5_MAX_TOTAL_I64 = (uint64_t)1 << 31; // 16 GiB payload cap
static constexpr uint64_t Q5_MAX_DENSE_MAP = (uint64_t)1 << 26; // 64M dense table entries cap
static constexpr int      Q5_MAX_KNOTS = 24;

struct QuintetHeader {
    uint64_t magic;
    // 0..8 (9)
    int64_t mode;         // 1 = UPPER (L Phi <= lam Phi), 2 = LOWER (L Phi >= lam Phi)
    int64_t P;            // target lambda = P / Q
    int64_t Q;
    int64_t S;            // weight tick scale (2^40)
    int64_t S_r;          // centering scale (2^30)
    int64_t alpha_lo_num; // 0 for UPPER, >0 (e.g. S_r/4) for LOWER comparison shift alpha*(c-2)
    int64_t f_p0;         // analytic f tail params in weight ticks S
    int64_t f_p1;
    int64_t f_p2;
    // 9..34 (26)
    int64_t DF;
    int64_t d2, K2;
    int64_t d3m, K3m, d3e, K3e;
    int64_t d4i, K4i, SP4, dS3, KS3, d4e, K4e;
    int64_t d5i, K5i, SP5, d5o, K5o, d5c, K5c, SB5, dS4, KS4, d5e, K5e;
    // 35..36 (2)
    int64_t Q_K;          // exact common denominator of FY knots
    int64_t nknots;       // number of active FY knots (3..24)
    // 37..60 (24)
    int64_t knot_num[Q5_MAX_KNOTS];
    // 61..62 (2)
    int64_t reserved0;
    int64_t reserved1;
};
static_assert(sizeof(QuintetHeader) == 512, "QuintetHeader must be exactly 512 bytes");

// Exact signed 128-bit helpers (den MUST be > 0; checked narrowing to int64_t!)
inline q5_i128 q5_ceil_div(q5_i128 num, q5_i128 den) {
    if (den <= 0) { std::fprintf(stderr, "ERROR: non-positive denominator in q5_ceil_div\n"); std::exit(2); }
    if (num >= 0) return (num + den - 1) / den;
    return -((-num) / den);
}
inline q5_i128 q5_floor_div(q5_i128 num, q5_i128 den) {
    if (den <= 0) { std::fprintf(stderr, "ERROR: non-positive denominator in q5_floor_div\n"); std::exit(2); }
    if (num >= 0) return num / den;
    return -(((-num) + den - 1) / den);
}
inline q5_i128 q5_abs128(q5_i128 x) { return x < 0 ? -x : x; }
inline int64_t q5_to_i64(q5_i128 x, const char* ctx) {
    if (x < (q5_i128)INT64_MIN || x > (q5_i128)INT64_MAX) {
        std::fprintf(stderr, "ERROR: int64_t narrowing overflow in %s\n", ctx);
        std::exit(2);
    }
    return (int64_t)x;
}

inline q5_u128 q5_end_count(int64_t dc, int64_t Kc, q5_u128* out_nin_ex = nullptr) {
    q5_u128 nin = 0;
    for (int64_t d = 2; d < dc; ++d) nin += (q5_u128)std::min(d, Kc);
    if (out_nin_ex) *out_nin_ex = nin;
    q5_u128 novf = (dc > Kc + 1) ? (q5_u128)(dc - 1 - Kc) : 0;
    return nin + novf + (q5_u128)(Kc + 1) + 1;
}

inline q5_u128 q5_inner_count(int64_t dc, int64_t mnC, int64_t mxC) {
    q5_u128 cnt = 0;
    for (int64_t d = std::max<int64_t>(2, mnC); d < dc; ++d) {
        int64_t hi = std::min(d, mxC);
        if (hi >= mnC) cnt += (q5_u128)(hi - mnC + 1);
    }
    return cnt;
}

inline int q5_validate_header(const QuintetHeader& h) {
    if (h.magic != Q5_MAGIC) return 1;
    if (h.mode != 1 && h.mode != 2) return 2;
    if (h.P <= 0 || h.Q <= 0 || h.P >= h.Q || h.Q > Q5_MAX_Q) return 3;
    if (h.S != Q5_WEIGHT_SCALE || h.S_r != Q5_CENTER_SCALE) return 4;
    if (h.reserved0 != 0 || h.reserved1 != 0) return 5;
    if (h.mode == 1 && (h.f_p2 != 0 || h.alpha_lo_num != 0)) return 6;
    if (h.mode == 2 && (h.alpha_lo_num < 0 || h.alpha_lo_num > h.S_r)) return 7;
    if (h.f_p0 < -Q5_MAX_W || h.f_p0 > Q5_MAX_W ||
        h.f_p1 < -Q5_MAX_W || h.f_p1 > Q5_MAX_W ||
        h.f_p2 < -Q5_MAX_W || h.f_p2 > Q5_MAX_W) return 30;
    if (h.DF < 20 || h.DF > 4096) return 8;
    if (h.d2 < 6 || h.d2 >= h.DF || h.d2 > 256) return 9;
    if (h.K2 < 2 || h.K2 >= h.d2 || h.K2 > 64) return 10;
    if (h.d3m < 4 || h.d3m > h.d2 || h.d3m > 128 || h.K3m < 2 || h.K3m >= h.d3m || h.K3m > h.K2 || h.K3m > 32) return 11;
    if (h.d3e < 4 || h.d3e > h.d2 || h.d3e > 128 || h.K3e < 2 || h.K3e >= h.d3e || h.K3e > h.K2 || h.K3e > 32) return 12;
    if (h.d4i < 3 || h.d4i > std::min(h.d3m, h.d3e) || h.d4i > 64 ||
        h.K4i < 2 || h.K4i >= h.d4i || h.K4i > std::min(h.K3m, h.K3e) || h.K4i > 24) return 13;
    if (h.SP4 < 4 || h.SP4 > 2 * (h.d4i - 1)) return 14;
    if (h.dS3 < 4 || h.dS3 > h.d3m || h.dS3 > 64 || h.KS3 < 3 || h.KS3 >= h.dS3 || h.KS3 > h.K3m || h.KS3 > 24) return 15;
    if (h.d4e < 3 || h.d4e > h.d3e || h.d4e > 64 || h.K4e < 2 || h.K4e >= h.d4e || h.K4e > h.K3e || h.K4e > 24) return 16;
    if (h.d5i < 3 || h.d5i > std::min(h.d4i, h.d4e) || h.d5i > 32 ||
        h.K5i < 2 || h.K5i >= h.d5i || h.K5i > std::min(h.K4i, h.K4e) || h.K5i > 16) return 17;
    if (h.SP5 < 6 || h.SP5 > 3 * (h.d5i - 1)) return 18;
    if (h.d5o < 4 || h.d5o > std::min(h.d4i, h.dS3) || h.d5o > 32 ||
        h.K5o < 3 || h.K5o >= h.d5o || h.K5o > std::min(h.K4i, h.KS3) || h.K5o > 16) return 19;
    if (h.d5c < 3 || h.d5c > std::min(h.d4i, h.d4e) || h.d5c > 32 ||
        h.K5c < 2 || h.K5c >= h.d5c || h.K5c > std::min(h.K4i, h.K4e) || h.K5c > 16) return 20;
    if (h.SB5 < 5 || h.SB5 > h.SP4 || h.SB5 > (h.d5o - 1) + (h.d5c - 1)) return 21;
    if (h.dS4 < 5 || h.dS4 > h.dS3 || h.dS4 > 32 || h.KS4 < 4 || h.KS4 >= h.dS4 || h.KS4 > h.KS3 || h.KS4 > 16) return 22;
    if (h.d5e < 3 || h.d5e > h.d4e || h.d5e > 32 || h.K5e < 2 || h.K5e >= h.d5e || h.K5e > h.K4e || h.K5e > 16) return 23;
    if (h.Q_K <= 0 || h.Q_K > 65536 || h.nknots < 3 || h.nknots > Q5_MAX_KNOTS) return 24;
    if (h.knot_num[0] != 0) return 25;
    // Range-check ALL knot numerators against [0, Q_K] BEFORE any multiplication (zero signed overflow UB!)
    for (int i = 1; i < h.nknots; ++i)
        if (h.knot_num[i] <= h.knot_num[i - 1] || h.knot_num[i] > h.Q_K) return 28;
    for (int i = (int)h.nknots; i < Q5_MAX_KNOTS; ++i)
        if (h.knot_num[i] != 0) return 29;
    if (h.knot_num[1] * h.d2 != h.Q_K) return 26;           // xi_1 == 1/d2 EXACTLY!
    if (2 * h.knot_num[h.nknots - 1] != h.Q_K) return 27;   // xi_{m-1} == 1/2 EXACTLY!
    return 0;
}

// Checked pre-allocation exact byte-size & resource-bound check (zero allocation, pure q5_u128!)
inline int q5_compute_expected_bytes(const QuintetHeader& h, uint64_t& out_bytes) {
    int vh = q5_validate_header(h);
    if (vh != 0) return vh;
    q5_u128 g2_in = 0;
    q5_u128 nG2 = q5_end_count(h.d2, h.K2, &g2_in);
    q5_u128 nG3 = q5_end_count(h.d3e, h.K3e);
    q5_u128 nG4 = q5_end_count(h.d4e, h.K4e);
    q5_u128 nG5 = q5_end_count(h.d5e, h.K5e);
    if (nG2 > 32768 || nG3 > 8192 || nG4 > 400 || nG5 > 90) return 31;
    if (nG4 * nG4 * nG4 > Q5_MAX_DENSE_MAP) return 32;
    if (nG5 * nG5 * nG5 * nG5 > Q5_MAX_DENSE_MAP) return 33;

    q5_u128 nM3  = q5_inner_count(h.d3m, 2, h.K3m);
    q5_u128 nM4v = q5_inner_count(h.d4i, 2, h.K4i);
    q5_u128 nMS3 = q5_inner_count(h.dS3, 3, h.KS3);
    q5_u128 nI52 = q5_inner_count(h.d5i, 2, h.K5i);
    q5_u128 nI53 = q5_inner_count(h.d5o, 3, h.K5o);
    q5_u128 nI5c = q5_inner_count(h.d5c, 2, h.K5c);
    q5_u128 nMS4 = q5_inner_count(h.dS4, 4, h.KS4);
    if (nM4v * nM4v > Q5_MAX_DENSE_MAP || nI52 * nI52 * nI52 > Q5_MAX_DENSE_MAP || nI53 * nI5c > Q5_MAX_DENSE_MAP) return 34;

    // Count ordered P4, P5, B5 blocks exactly without allocating!
    q5_u128 g4_p2 = nG4 * (nG4 + 1) / 2, g5_p2 = nG5 * (nG5 + 1) / 2;
    q5_u128 nP4_ord = 0, n_p4 = 0;
    for (int64_t db = 2; db < h.d4i; ++db) {
        int64_t cb_cnt = std::min(db, h.K4i) - 1;
        if (cb_cnt <= 0) continue;
        for (int64_t dc = 2; dc < h.d4i; ++dc) {
            if (db + dc > h.SP4) continue;
            int64_t cc_cnt = std::min(dc, h.K4i) - 1;
            if (cc_cnt <= 0) continue;
            nP4_ord += (q5_u128)cb_cnt * (q5_u128)cc_cnt;
        }
    }
    // Canonical P4 blocks: b == c gives g4_p2, b < c gives nG4^2; total pairs = nP4_ord, diagonal pairs = count(2*db <= SP4)
    q5_u128 nP4_diag = 0;
    for (int64_t db = 2; db < h.d4i; ++db)
        if (2 * db <= h.SP4) {
            int64_t cb_cnt = std::min(db, h.K4i) - 1;
            if (cb_cnt > 0) nP4_diag += (q5_u128)cb_cnt;
        }
    n_p4 = nP4_diag * g4_p2 + ((nP4_ord - nP4_diag) / 2) * (nG4 * nG4);

    q5_u128 nP5_ord = 0, nP5_diag = 0;
    for (int64_t db = 2; db < h.d5i; ++db) {
        int64_t cb_cnt = std::min(db, h.K5i) - 1;
        if (cb_cnt <= 0) continue;
        for (int64_t dc = 2; dc < h.d5i; ++dc) {
            if (db + dc > h.SP4) continue;
            int64_t cc_cnt = std::min(dc, h.K5i) - 1;
            if (cc_cnt <= 0) continue;
            for (int64_t dd = 2; dd < h.d5i; ++dd) {
                if (dc + dd > h.SP4 || db + dc + dd > h.SP5) continue;
                int64_t cd_cnt = std::min(dd, h.K5i) - 1;
                if (cd_cnt <= 0) continue;
                nP5_ord += (q5_u128)cb_cnt * (q5_u128)cc_cnt * (q5_u128)cd_cnt;
                if (db == dd) nP5_diag += (q5_u128)cb_cnt * (q5_u128)cc_cnt;
            }
        }
    }
    q5_u128 n_p5 = nP5_diag * g5_p2 + ((nP5_ord - nP5_diag) / 2) * (nG5 * nG5);

    q5_u128 nB5 = 0;
    for (int64_t do_ = 3; do_ < h.d5o; ++do_) {
        int64_t co_cnt = std::min(do_, h.K5o) - 2;
        if (co_cnt <= 0) continue;
        for (int64_t dc = 2; dc < h.d5c; ++dc) {
            if (do_ + dc > h.SB5) continue;
            int64_t cc_cnt = std::min(dc, h.K5c) - 1;
            if (cc_cnt <= 0) continue;
            nB5 += (q5_u128)co_cnt * (q5_u128)cc_cnt;
        }
    }

    q5_u128 nG4_tri  = nG4 * (nG4 + 1) * (nG4 + 2) / 6;
    q5_u128 nG5_tri  = nG5 * (nG5 + 1) * (nG5 + 2) / 6;
    q5_u128 nG5_quad = nG5 * (nG5 + 1) * (nG5 + 2) * (nG5 + 3) / 24;
    if (nP4_ord > INT32_MAX / 2 || nP5_ord > INT32_MAX / 2 || nB5 > INT32_MAX / 2 ||
        nG4_tri > INT32_MAX / 2 || nG5_tri > INT32_MAX / 2 || nG5_quad > INT32_MAX / 2) return 35;

    q5_u128 tot = 0;
    tot += g2_in;                                          // r1
    tot += g2_in * (q5_u128)h.nknots;                      // knots
    tot += (q5_u128)(h.DF - 2) * (q5_u128)(h.DF + 3) / 2;  // f
    tot += nG2 * (nG2 + 1) / 2;                            // g
    tot += nM3 * (nG3 * (nG3 + 1) / 2);                    // h3
    tot += n_p4;                                           // p4
    tot += nMS3 * nG4_tri;                                 // s3
    tot += n_p5;                                           // p5
    tot += nB5 * g5_p2 * nG5;                              // b5
    tot += nMS4 * nG5_quad;                                // s4
    tot += 2 * (nM3 * nG3);                                // e23
    tot += 2 * (nP4_ord * nG4);                            // e34p
    tot += 2 * (nMS3 * g4_p2);                             // e34s
    tot += 2 * (nP5_ord * nG5);                            // e45p
    tot += 2 * (nB5 * nG5 * nG5);                          // e45bp
    tot += 2 * (nB5 * g5_p2);                              // e45bs
    tot += 2 * (nMS4 * nG5_tri);                           // e45s4

    if (tot > Q5_MAX_TOTAL_I64) return 36;
    out_bytes = (uint64_t)sizeof(QuintetHeader) + (uint64_t)tot * sizeof(int64_t);
    return 0;
}

struct Q5EndSpace {
    int dcut = 0, Kcut = 0, N = 0, INF = 0, maxR = 0, Nin_ex = 0, tail_base = 0;
    std::vector<int> kind, D, C, L;
    std::vector<std::vector<int>> idx;
    std::vector<int> ovf_idx;
    std::vector<std::vector<int>> by_r;

    void build(int dc, int Kc) {
        dcut = dc; Kcut = Kc;
        kind.clear(); D.clear(); C.clear(); L.clear();
        idx.assign(dc + 2, std::vector<int>(Kc + 2, -1));
        ovf_idx.assign(dc + 2, -1);
        for (int d = 2; d < dc; ++d) {
            for (int c = 1; c <= std::min(d, Kc); ++c) {
                int id = (int)kind.size();
                kind.push_back(0); D.push_back(d); C.push_back(c); L.push_back(d - c);
                idx[d][c] = id;
            }
        }
        Nin_ex = (int)kind.size();
        for (int d = Kc + 1; d < dc; ++d) {
            int id = (int)kind.size();
            kind.push_back(3); D.push_back(d); C.push_back(Kc + 1); L.push_back(0);
            ovf_idx[d] = id;
        }
        tail_base = (int)kind.size();
        for (int c = 1; c <= Kc + 1; ++c) {
            kind.push_back(1); D.push_back(dc); C.push_back(c); L.push_back(std::max(0, dc - c));
        }
        INF = (int)kind.size();
        kind.push_back(2); D.push_back(0); C.push_back(0); L.push_back(0);
        N = (int)kind.size();
        maxR = dc + 1;
        by_r.assign(maxR + 1, {});
        for (int d = 2; d < dc; ++d) {
            for (int c = 1; c <= std::min(d, Kc); ++c) by_r[d].push_back(idx[d][c]);
            if (ovf_idx[d] >= 0) by_r[d].push_back(ovf_idx[d]);
        }
        for (int c = 1; c <= Kc + 1; ++c) by_r[dc].push_back(t_id(c));
        by_r[maxR].push_back(INF);
    }
    inline bool is_inner(int u) const { return kind[u] == 0 || kind[u] == 3; }
    inline int t_id(int c) const { return tail_base + std::min(c, Kcut + 1) - 1; }
    inline int up(int u) const {
        if (kind[u] == 0) return (D[u] + 1 < dcut) ? idx[D[u] + 1][C[u]] : t_id(C[u]);
        if (kind[u] == 3) return (D[u] + 1 < dcut) ? ovf_idx[D[u] + 1] : t_id(Kcut + 1);
        return INF;
    }
    inline int bud(int u) const {
        if (kind[u] == 0) {
            if (L[u] < 1) return u;
            return (C[u] + 1 <= Kcut) ? idx[D[u]][C[u] + 1] : ovf_idx[D[u]];
        }
        if (kind[u] == 1 && C[u] <= Kcut) return t_id(C[u] + 1);
        return u;
    }
    // Exact integer projection of fine state fw in F -> (ce, num, den):
    // weight on ce is num/den, weight on INF is (den-num)/den.
    inline void proj_int(const Q5EndSpace& F, int fw, int& ce, int& num, int& den) const {
        if (F.kind[fw] == 2) { ce = INF; num = 1; den = 1; }
        else if (F.is_inner(fw) && F.D[fw] < dcut) {
            ce = (F.kind[fw] == 0 && F.C[fw] <= Kcut) ? idx[F.D[fw]][F.C[fw]] : ovf_idx[F.D[fw]];
            num = 1; den = 1;
        } else {
            ce = t_id(F.C[fw]); num = dcut; den = F.D[fw];
        }
    }
};

struct Q5InnerSpace {
    int dcut = 0, minC = 0, maxC = 0, N = 0;
    std::vector<int> D, C, L;
    std::vector<std::vector<int>> idx;
    void build(int dc, int mnC, int mxC) {
        dcut = dc; minC = mnC; maxC = mxC;
        D.clear(); C.clear(); L.clear();
        idx.assign(dc + 2, std::vector<int>(mxC + 2, -1));
        for (int d = std::max(2, mnC); d < dc; ++d)
            for (int c = mnC; c <= std::min(d, mxC); ++c) {
                int id = (int)D.size();
                D.push_back(d); C.push_back(c); L.push_back(d - c);
                idx[d][c] = id;
            }
        N = (int)D.size();
    }
    inline int up(int u) const { return (D[u] + 1 < dcut) ? idx[D[u] + 1][C[u]] : -1; }
    inline int bud(int u) const { return (L[u] >= 1 && C[u] + 1 <= maxC) ? idx[D[u]][C[u] + 1] : -1; }
};

struct Q5ProjMap {
    int NC = 0, NF = 0, C_INF = 0, F_INF = 0;
    std::vector<int> exact_fine;
    struct FineEntry { int fw, num, den; }; // weight on ce = num/den, on INF = (den-num)/den
    std::vector<std::vector<FineEntry>> buckets;
    std::vector<int> non_exact;

    void build(const Q5EndSpace& F, const Q5EndSpace& Ce) {
        NC = Ce.N; NF = F.N; C_INF = Ce.INF; F_INF = F.INF;
        exact_fine.assign(NC, -1);
        buckets.assign(NC, {});
        non_exact.clear();
        for (int e = 0; e < NC; ++e) {
            if (Ce.kind[e] == 0) exact_fine[e] = F.idx[Ce.D[e]][Ce.C[e]];
            else if (Ce.kind[e] == 3 || Ce.kind[e] == 1) non_exact.push_back(e);
        }
        for (int fw = 0; fw < NF; ++fw) {
            if (fw == F_INF) continue;
            int ce, num, den;
            Ce.proj_int(F, fw, ce, num, den);
            if (Ce.kind[ce] == 3 || Ce.kind[ce] == 1) buckets[ce].push_back({fw, num, den});
        }
        for (int e : non_exact) assert(!buckets[e].empty());
    }
};

struct Q5Layout {
    QuintetHeader h;
    Q5EndSpace G2, G3, G4, G5;
    Q5ProjMap PM23, PM34, PM45;
    Q5InnerSpace M3, M4v, MS3, I5_2, I5_3, I5_2c, MS4;

    // Ordered P4 internal pairs (for 1D envelopes env34_p) and canonical P4 storage (b<=c)!
    int NP4_ord = 0;
    std::vector<std::pair<int, int>> p4_ord;
    std::vector<int> p4_ord_map;
    std::vector<size_t> p4_block_off; // indexed by (b * M4v.N + c) for b <= c!
    size_t n_p4 = 0;

    // Ordered P5 internal triples (for 1D envelopes env45_p) and canonical P5 storage (b<=d)!
    int NP5_ord = 0;
    std::vector<std::array<int, 3>> p5_ord;
    std::vector<int> p5_ord_map;
    std::vector<size_t> p5_block_off; // indexed by ((b * N52 + c) * N52 + d) for b <= d!
    size_t n_p5 = 0;

    // B5 internal pairs (o, c)
    int NB5_int = 0;
    std::vector<std::pair<int, int>> b5_int;
    std::vector<int> b5_int_map;

    int N_G4_tri = 0, N_G5_tri = 0, N_G5_quad = 0;
    std::vector<std::array<int, 3>> tri_G4, tri_G5;
    std::vector<int> tri_G4_map, tri_G5_map;
    std::vector<std::array<int, 4>> quad_G5;
    std::vector<int> quad_G5_map;

    // Payload sizes (in int64_t counts) & byte offsets
    size_t n_r1=0, n_knots=0, n_f=0, n_g=0, n_h3=0, n_s3=0, n_b5=0, n_s4=0;
    size_t n_e23=0, n_e34p=0, n_e34s=0, n_e45p=0, n_e45bp=0, n_e45bs=0, n_e45s4=0;
    size_t off_r1=0, off_knots=0, off_f=0, off_g=0, off_h3=0, off_p4=0, off_s3=0, off_p5=0, off_b5=0, off_s4=0;
    size_t off_e23_m=0, off_e23_r=0, off_e34p_m=0, off_e34p_r=0, off_e34s_m=0, off_e34s_r=0;
    size_t off_e45p_m=0, off_e45p_r=0, off_e45bp_m=0, off_e45bp_r=0, off_e45bs_m=0, off_e45bs_r=0, off_e45s4_m=0, off_e45s4_r=0;
    size_t total_i64 = 0, total_bytes = 0;

    static inline size_t sym_off(int a, int b) {
        if (a > b) std::swap(a, b);
        return (size_t)b * (b + 1) / 2 + a;
    }
    inline int f_id(int D, int c) const { return (D - 2) * (D + 3) / 2 + c; }
    inline size_t g_off(int u, int v) const { return sym_off(u, v); }
    inline size_t h3_off(int m, int a, int b) const {
        return (size_t)m * ((size_t)G3.N * (G3.N + 1) / 2) + sym_off(a, b);
    }
    inline size_t p4_off(int b, int c, int a, int d) const {
        if (b > c) { std::swap(b, c); std::swap(a, d); }
        size_t base = p4_block_off[(size_t)b * M4v.N + c];
        return base + ((b == c) ? sym_off(a, d) : ((size_t)a * G4.N + d));
    }
    inline size_t s3_off(int o, int a, int b, int c) const {
        return (size_t)o * N_G4_tri + tri_G4_map[((size_t)a * G4.N + b) * G4.N + c];
    }
    inline size_t p5_off(int b, int c, int d, int a, int e) const {
        if (b > d) { std::swap(b, d); std::swap(a, e); }
        size_t base = p5_block_off[((size_t)b * I5_2.N + c) * I5_2.N + d];
        return base + ((b == d) ? sym_off(a, e) : ((size_t)a * G5.N + e));
    }
    inline size_t b5_off(int bid, int a, int b, int d) const {
        size_t n2 = (size_t)G5.N * (G5.N + 1) / 2;
        return ((size_t)bid * n2 + sym_off(a, b)) * G5.N + d;
    }
    inline size_t s4_off(int o, int a, int b, int c, int d) const {
        size_t N5 = G5.N;
        return (size_t)o * N_G5_quad + quad_G5_map[(((size_t)a * N5 + b) * N5 + c) * N5 + d];
    }

    void init(const QuintetHeader& hdr) {
        h = hdr;
        assert(q5_validate_header(h) == 0);
        G2.build((int)h.d2, (int)h.K2);
        G3.build((int)h.d3e, (int)h.K3e);
        G4.build((int)h.d4e, (int)h.K4e);
        G5.build((int)h.d5e, (int)h.K5e);
        PM23.build(G2, G3);
        PM34.build(G3, G4);
        PM45.build(G4, G5);

        M3.build((int)h.d3m, 2, (int)h.K3m);
        M4v.build((int)h.d4i, 2, (int)h.K4i);
        MS3.build((int)h.dS3, 3, (int)h.KS3);
        I5_2.build((int)h.d5i, 2, (int)h.K5i);
        I5_3.build((int)h.d5o, 3, (int)h.K5o);
        I5_2c.build((int)h.d5c, 2, (int)h.K5c);
        MS4.build((int)h.dS4, 4, (int)h.KS4);

        // P4 ordered & canonical palindrome blocks
        p4_ord.clear();
        p4_ord_map.assign((size_t)M4v.N * M4v.N, -1);
        p4_block_off.assign((size_t)M4v.N * M4v.N, (size_t)-1);
        n_p4 = 0;
        for (int b = 0; b < M4v.N; ++b)
            for (int c = 0; c < M4v.N; ++c)
                if (M4v.D[b] + M4v.D[c] <= h.SP4) {
                    p4_ord_map[(size_t)b * M4v.N + c] = (int)p4_ord.size();
                    p4_ord.push_back({b, c});
                    if (b <= c) {
                        p4_block_off[(size_t)b * M4v.N + c] = n_p4;
                        n_p4 += (b == c) ? ((size_t)G4.N * (G4.N + 1) / 2) : ((size_t)G4.N * G4.N);
                    }
                }
        NP4_ord = (int)p4_ord.size();

        // P5 ordered & canonical palindrome blocks
        int N52 = I5_2.N;
        p5_ord.clear();
        p5_ord_map.assign((size_t)N52 * N52 * N52, -1);
        p5_block_off.assign((size_t)N52 * N52 * N52, (size_t)-1);
        n_p5 = 0;
        for (int b = 0; b < N52; ++b)
            for (int c = 0; c < N52; ++c)
                if (I5_2.D[b] + I5_2.D[c] <= h.SP4)
                    for (int d = 0; d < N52; ++d)
                        if (I5_2.D[c] + I5_2.D[d] <= h.SP4 && I5_2.D[b] + I5_2.D[c] + I5_2.D[d] <= h.SP5) {
                            size_t flat = ((size_t)b * N52 + c) * N52 + d;
                            p5_ord_map[flat] = (int)p5_ord.size();
                            p5_ord.push_back({b, c, d});
                            if (b <= d) {
                                p5_block_off[flat] = n_p5;
                                n_p5 += (b == d) ? ((size_t)G5.N * (G5.N + 1) / 2) : ((size_t)G5.N * G5.N);
                            }
                        }
        NP5_ord = (int)p5_ord.size();

        // B5 internal pairs
        b5_int.clear();
        b5_int_map.assign((size_t)I5_3.N * I5_2c.N, -1);
        for (int o = 0; o < I5_3.N; ++o)
            for (int c = 0; c < I5_2c.N; ++c)
                if (I5_3.D[o] + I5_2c.D[c] <= h.SB5) {
                    b5_int_map[(size_t)o * I5_2c.N + c] = (int)b5_int.size();
                    b5_int.push_back({o, c});
                }
        NB5_int = (int)b5_int.size();

        auto build_tri = [](int N, std::vector<std::array<int, 3>>& lst, std::vector<int>& mp) {
            lst.clear(); mp.assign((size_t)N * N * N, -1);
            for (int a = 0; a < N; ++a) for (int b = a; b < N; ++b) for (int c = b; c < N; ++c) {
                int id = (int)lst.size();
                lst.push_back({a, b, c});
                int p[3] = {a, b, c};
                do { mp[((size_t)p[0] * N + p[1]) * N + p[2]] = id; } while (std::next_permutation(p, p + 3));
            }
        };
        build_tri(G4.N, tri_G4, tri_G4_map); N_G4_tri = (int)tri_G4.size();
        build_tri(G5.N, tri_G5, tri_G5_map); N_G5_tri = (int)tri_G5.size();

        int N5 = G5.N;
        quad_G5.clear(); quad_G5_map.assign((size_t)N5 * N5 * N5 * N5, -1);
        for (int a = 0; a < N5; ++a) for (int b = a; b < N5; ++b) for (int c = b; c < N5; ++c) for (int d = c; d < N5; ++d) {
            int id = (int)quad_G5.size();
            quad_G5.push_back({a, b, c, d});
            int p[4] = {a, b, c, d};
            do { quad_G5_map[(((size_t)p[0] * N5 + p[1]) * N5 + p[2]) * N5 + p[3]] = id; } while (std::next_permutation(p, p + 4));
        }
        N_G5_quad = (int)quad_G5.size();

        size_t g4_p2 = (size_t)G4.N * (G4.N + 1) / 2;
        size_t g5_p2 = (size_t)G5.N * (G5.N + 1) / 2;
        n_r1    = (size_t)G2.Nin_ex;
        n_knots = (size_t)G2.Nin_ex * (size_t)h.nknots;
        n_f     = (size_t)f_id((int)h.DF, 0);
        n_g     = (size_t)G2.N * (G2.N + 1) / 2;
        n_h3    = (size_t)M3.N * ((size_t)G3.N * (G3.N + 1) / 2);
        n_s3    = (size_t)MS3.N * N_G4_tri;
        n_b5    = (size_t)NB5_int * g5_p2 * G5.N;
        n_s4    = (size_t)MS4.N * N_G5_quad;

        n_e23   = (size_t)M3.N * G3.N;
        n_e34p  = (size_t)NP4_ord * G4.N;
        n_e34s  = (size_t)MS3.N * g4_p2;
        n_e45p  = (size_t)NP5_ord * G5.N;
        n_e45bp = (size_t)NB5_int * G5.N * G5.N;
        n_e45bs = (size_t)NB5_int * g5_p2;
        n_e45s4 = (size_t)MS4.N * N_G5_tri;

        size_t cur = 0;
        auto alloc = [&](size_t sz) -> size_t { size_t o = cur; cur += sz; return o; };
        off_r1      = alloc(n_r1);
        off_knots   = alloc(n_knots);
        off_f       = alloc(n_f);
        off_g       = alloc(n_g);
        off_h3      = alloc(n_h3);
        off_p4      = alloc(n_p4);
        off_s3      = alloc(n_s3);
        off_p5      = alloc(n_p5);
        off_b5      = alloc(n_b5);
        off_s4      = alloc(n_s4);
        off_e23_m   = alloc(n_e23);   off_e23_r   = alloc(n_e23);
        off_e34p_m  = alloc(n_e34p);  off_e34p_r  = alloc(n_e34p);
        off_e34s_m  = alloc(n_e34s);  off_e34s_r  = alloc(n_e34s);
        off_e45p_m  = alloc(n_e45p);  off_e45p_r  = alloc(n_e45p);
        off_e45bp_m = alloc(n_e45bp); off_e45bp_r = alloc(n_e45bp);
        off_e45bs_m = alloc(n_e45bs); off_e45bs_r = alloc(n_e45bs);
        off_e45s4_m = alloc(n_e45s4); off_e45s4_r = alloc(n_e45s4);
        total_i64 = cur;
        total_bytes = sizeof(QuintetHeader) + total_i64 * sizeof(int64_t);
    }
};
