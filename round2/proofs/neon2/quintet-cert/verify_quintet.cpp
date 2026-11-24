#include "quintet_format.hpp"
#include <charconv>
#include <cinttypes>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <iostream>
#include <omp.h>
#include <string>
#include <vector>

using namespace std;

/*
 * Standalone ZERO-FLOATING-POINT Exact 128-bit Verifier for FT2_U5V1 Certificates.
 *
 * Reads ONLY the binary snapshot (read-only), performs zero floating-point arithmetic,
 * computes an exact FIPS 180-4 SHA-256 digest over the exact in-memory header+payload bytes
 * it verifies, and verifies:
 *   1. Header invariants (validated on original header before any --override-pq), checked
 *      pre-allocation resource/size bounds, exact file size, g[INF, INF] == 0, magnitude & radius bounds.
 *   2. ALL 7 cross-level envelope families (1D, 2D gen/sym, 3D sym) over ALL fine states/tuples.
 *   3. ALL finite rows of f (including D=DF-1 analytic seam), g (including exact FY primal & dual),
 *      H3, P4, S3, P5, B5, S4 in directed sub-tick 128-bit arithmetic (K_SUB = 1024) with checked
 *      int64_t narrowing.
 *   4. Per-core multiplicity aggregation & global comparability bounds (D < DF + far-tail g).
 *
 * Also supports --override-pq P Q and --perturb <target> for fail-closed negative-control testing.
 */

static constexpr q5_i128 K_SUB = 1024;
static constexpr int64_t MAX_W = Q5_MAX_W;
static constexpr int64_t MAX_E = Q5_MAX_E;

static bool parse_i64_strict(const char* s, int64_t& out) {
    if (!s || s[0] == '\0') return false;
    size_t len = strlen(s);
    int64_t v = 0;
    auto res = from_chars(s, s + len, v, 10);
    if (res.ec != errc() || res.ptr != s + len) return false;
    out = v;
    return true;
}

// Self-contained FIPS 180-4 SHA-256 (zero external dependencies, zero floats)
struct Q5Sha256 {
    uint32_t st[8] = {
        0x6a09e667u, 0xbb67ae85u, 0x3c6ef372u, 0xa54ff53au,
        0x510e527fu, 0x9b05688cu, 0x1f83d9abu, 0x5be0cd19u
    };
    uint64_t bits = 0;
    uint8_t blk[64];
    size_t blen = 0;

    static inline uint32_t rotr(uint32_t x, uint32_t n) { return (x >> n) | (x << (32u - n)); }

    void compress(const uint8_t b[64]) {
        static constexpr uint32_t K[64] = {
            0x428a2f98u,0x71374491u,0xb5c0fbcfu,0xe9b5dba5u,0x3956c25bu,0x59f111f1u,0x923f82a4u,0xab1c5ed5u,
            0xd807aa98u,0x12835b01u,0x243185beu,0x550c7dc3u,0x72be5d74u,0x80deb1feu,0x9bdc06a7u,0xc19bf174u,
            0xe49b69c1u,0xefbe4786u,0x0fc19dc6u,0x240ca1ccu,0x2de92c6fu,0x4a7484aau,0x5cb0a9dcu,0x76f988dau,
            0x983e5152u,0xa831c66du,0xb00327c8u,0xbf597fc7u,0xc6e00bf3u,0xd5a79147u,0x06ca6351u,0x14292967u,
            0x27b70a85u,0x2e1b2138u,0x4d2c6dfcu,0x53380d13u,0x650a7354u,0x766a0abbu,0x81c2c92eu,0x92722c85u,
            0xa2bfe8a1u,0xa81a664bu,0xc24b8b70u,0xc76c51a3u,0xd192e819u,0xd6990624u,0xf40e3585u,0x106aa070u,
            0x19a4c116u,0x1e376c08u,0x2748774cu,0x34b0bcb5u,0x391c0cb3u,0x4ed8aa4au,0x5b9cca4fu,0x682e6ff3u,
            0x748f82eeu,0x78a5636fu,0x84c87814u,0x8cc70208u,0x90befffau,0xa4506cebu,0xbef9a3f7u,0xc67178f2u
        };
        uint32_t w[64];
        for (int i = 0; i < 16; ++i)
            w[i] = ((uint32_t)b[4*i] << 24) | ((uint32_t)b[4*i+1] << 16) | ((uint32_t)b[4*i+2] << 8) | (uint32_t)b[4*i+3];
        for (int i = 16; i < 64; ++i) {
            uint32_t s0 = rotr(w[i-15], 7) ^ rotr(w[i-15], 18) ^ (w[i-15] >> 3);
            uint32_t s1 = rotr(w[i-2], 17) ^ rotr(w[i-2], 19) ^ (w[i-2] >> 10);
            w[i] = w[i-16] + s0 + w[i-7] + s1;
        }
        uint32_t a=st[0], bb=st[1], c=st[2], d=st[3], e=st[4], f=st[5], g=st[6], h=st[7];
        for (int i = 0; i < 64; ++i) {
            uint32_t S1 = rotr(e, 6) ^ rotr(e, 11) ^ rotr(e, 25);
            uint32_t ch = (e & f) ^ ((~e) & g);
            uint32_t t1 = h + S1 + ch + K[i] + w[i];
            uint32_t S0 = rotr(a, 2) ^ rotr(a, 13) ^ rotr(a, 22);
            uint32_t maj = (a & bb) ^ (a & c) ^ (bb & c);
            uint32_t t2 = S0 + maj;
            h = g; g = f; f = e; e = d + t1; d = c; c = bb; bb = a; a = t1 + t2;
        }
        st[0]+=a; st[1]+=bb; st[2]+=c; st[3]+=d; st[4]+=e; st[5]+=f; st[6]+=g; st[7]+=h;
    }

    void update(const void* data, size_t n) {
        const uint8_t* p = static_cast<const uint8_t*>(data);
        bits += (uint64_t)n * 8u;
        if (blen > 0) {
            size_t take = min(64 - blen, n);
            memcpy(blk + blen, p, take);
            blen += take; p += take; n -= take;
            if (blen == 64) { compress(blk); blen = 0; }
        }
        while (n >= 64) { compress(p); p += 64; n -= 64; }
        if (n > 0) { memcpy(blk, p, n); blen = n; }
    }

    string hex_digest() {
        blk[blen++] = 0x80u;
        if (blen > 56) { memset(blk + blen, 0, 64 - blen); compress(blk); blen = 0; }
        memset(blk + blen, 0, 56 - blen);
        for (int i = 0; i < 8; ++i) blk[56 + i] = (uint8_t)(bits >> (56 - 8 * i));
        compress(blk);
        char hex[65];
        static constexpr char H[] = "0123456789abcdef";
        for (int i = 0; i < 8; ++i)
            for (int j = 0; j < 8; ++j)
                hex[8 * i + j] = H[(st[i] >> (28 - 4 * j)) & 0xfu];
        hex[64] = '\0';
        return string(hex);
    }
};

static inline void add_rat(q5_i128& acc, q5_i128 num, q5_i128 den) {
    acc += q5_ceil_div(num * K_SUB, den);
}

static inline void add_ext_rate(q5_i128& acc, int sgn, int64_t dv, q5_i128 lo_n, q5_i128 lo_d, q5_i128 hi_n, q5_i128 hi_d) {
    q5_i128 sdv = (q5_i128)sgn * dv;
    if (sdv > 0) { if (hi_n > 0) add_rat(acc, sdv * hi_n, hi_d); }
    else if (sdv < 0) { if (lo_n > 0) add_rat(acc, sdv * lo_n, lo_d); }
}

static inline void add_src_env(q5_i128& acc, int sgn, const Q5EndSpace& E, int don, int64_t r1_int, int64_t S_r, int64_t mid, int64_t rad) {
    q5_i128 nf, df;
    if (E.kind[don] == 2) { nf = -(q5_i128)r1_int; df = S_r; }
    else { nf = (q5_i128)S_r - (q5_i128)E.D[don] * r1_int; df = (q5_i128)E.D[don] * S_r; }
    add_rat(acc, (q5_i128)sgn * nf * mid + q5_abs128(nf) * rad, df);
}

// Add exact tail-endpoint terms (gamma, promotion, and adjacent-inner quadratic credit) for kind==1!
static inline void add_tail_end(q5_i128& acc, int sgn, const Q5EndSpace& E, int u, int D_nbr,
                                int64_t val, int64_t v_up, int64_t v_bud, int64_t d0_nbr, int64_t nbr_up) {
    int dc = E.dcut, c = E.C[u], Kc = E.Kcut;
    // 1. Leaf gain: gamma in [g_lo, 1]
    q5_i128 glo_n = 0, glo_d = 1;
    if (c <= Kc) {
        if (D_nbr > 0) { glo_n = (q5_i128)(dc - c) * D_nbr + 1; glo_d = (q5_i128)(dc + 1) * D_nbr; }
        else { glo_n = dc - c; glo_d = dc + 1; }
        if (glo_n < 0) glo_n = 0;
    }
    add_ext_rate(acc, sgn, v_up - val, glo_n, glo_d, 1, 1);
    // 2. Promotion: b in [max(0, dc-c)/dc, 1] for c <= Kc (0 for overflow c == Kc+1)
    if (c <= Kc) {
        q5_i128 blo_n = (dc > c) ? (dc - c) : 0;
        add_ext_rate(acc, sgn, v_bud - val, blo_n, dc, 1, 1);
    }
    // 3. Adjacent inner neighbor quadratic credit: (1/dc) * max(sgn*d0_nbr, sgn*(nbr_up - val))
    if (D_nbr > 0) {
        q5_i128 s0 = (q5_i128)sgn * d0_nbr;
        q5_i128 s1 = (q5_i128)sgn * (nbr_up - val);
        add_rat(acc, (s1 > s0) ? s1 : s0, dc);
    }
}

static inline void add_ovf_end(q5_i128& acc, int sgn, const Q5EndSpace& E, int u, int64_t v_up, int64_t val) {
    // kind == 3: extra rate in [0, (2*D - Kcut - 2) / 2]
    q5_i128 hi_n = 2 * E.D[u] - E.Kcut - 2;
    add_ext_rate(acc, sgn, v_up - val, 0, 1, hi_n, 2);
}

// Exact FY primal evaluation A(inv_v) in exact rational num/den
static inline void eval_fy_primal(const QuintetHeader& h, const Q5EndSpace& G2, int v, const int64_t* ku, q5_i128& num, q5_i128& den) {
    int nk = (int)h.nknots;
    if (G2.kind[v] == 2) { num = ku[0]; den = 1; return; }
    if (G2.kind[v] == 1) { num = ku[1]; den = 1; return; } // inv = 1/d2 == xi_1 EXACTLY!
    int Dv = G2.D[v];
    int64_t QK = h.Q_K;
    // Find unique k with knot_num[k]*Dv <= QK <= knot_num[k+1]*Dv
    int k = 1;
    while (k + 1 < nk - 1 && h.knot_num[k + 1] * Dv < QK) ++k;
    int64_t n0 = h.knot_num[k], n1 = h.knot_num[k + 1];
    assert(n0 * Dv <= QK && QK <= n1 * Dv);
    q5_i128 dseg = (q5_i128)(n1 - n0) * Dv;
    num = (q5_i128)ku[k] * dseg + (q5_i128)(ku[k + 1] - ku[k]) * (QK - (q5_i128)n0 * Dv);
    den = dseg;
}

// Exact 1D projection verification
static bool check_env_1d(const Q5ProjMap& PM, const Q5EndSpace& F, const Q5EndSpace& Ce,
                         const int64_t* df, const int64_t* mid, const int64_t* rad, size_t& checks) {
    int CINF = Ce.INF;
    for (int fw = 0; fw < F.N; ++fw) {
        int ce, num, den;
        Ce.proj_int(F, fw, ce, num, den);
        q5_i128 w1 = num, w0 = den - num;
        q5_i128 pred = w1 * mid[ce] + w0 * mid[CINF];
        q5_i128 bnd  = w1 * rad[ce] + w0 * rad[CINF];
        if (q5_abs128((q5_i128)den * df[fw] - pred) > bnd) return false;
        ++checks;
    }
    (void)PM;
    return true;
}

// Exact 2D projection verification (sym=true means mid/rad stored as sym_off!)
static bool check_env_2d(const Q5EndSpace& F, const Q5EndSpace& Ce, const vector<int64_t>& Mf,
                         const int64_t* mid, const int64_t* rad, bool sym, size_t& checks) {
    int NF = F.N, NC = Ce.N, CINF = Ce.INF;
    auto getM = [&](const int64_t* arr, int a, int b) -> int64_t {
        return sym ? arr[Q5Layout::sym_off(a, b)] : arr[(size_t)a * NC + b];
    };
    for (int fw1 = 0; fw1 < NF; ++fw1) {
        int c1, n1, d1; Ce.proj_int(F, fw1, c1, n1, d1);
        int e1[2] = {CINF, c1}; q5_i128 w1[2] = {d1 - n1, n1};
        int fw2_start = sym ? fw1 : 0;
        for (int fw2 = fw2_start; fw2 < NF; ++fw2) {
            int c2, n2, d2; Ce.proj_int(F, fw2, c2, n2, d2);
            int e2[2] = {CINF, c2}; q5_i128 w2[2] = {d2 - n2, n2};
            q5_i128 pred = 0, bnd = 0;
            for (int a = 0; a < 2; ++a) if (w1[a]) for (int b = 0; b < 2; ++b) if (w2[b]) {
                q5_i128 w = w1[a] * w2[b];
                pred += w * getM(mid, e1[a], e2[b]);
                bnd  += w * getM(rad, e1[a], e2[b]);
            }
            if (q5_abs128((q5_i128)d1 * d2 * Mf[(size_t)fw1 * NF + fw2] - pred) > bnd) return false;
            ++checks;
        }
    }
    return true;
}

// Exact 3D symmetric projection verification
static bool check_env_3d(const Q5EndSpace& F, const Q5EndSpace& Ce,
                         const vector<array<int, 3>>& tri_F, const vector<int>& tri_C_map,
                         const vector<int64_t>& Tf, const int64_t* mid, const int64_t* rad, size_t& checks) {
    int NC = Ce.N, CINF = Ce.INF;
    for (size_t tid = 0; tid < tri_F.size(); ++tid) {
        auto [fw1, fw2, fw3] = tri_F[tid];
        int c1, n1, d1, c2, n2, d2, c3, n3, d3;
        Ce.proj_int(F, fw1, c1, n1, d1); Ce.proj_int(F, fw2, c2, n2, d2); Ce.proj_int(F, fw3, c3, n3, d3);
        int e1[2]={CINF,c1}, e2[2]={CINF,c2}, e3[2]={CINF,c3};
        q5_i128 w1[2]={d1-n1,n1}, w2[2]={d2-n2,n2}, w3[2]={d3-n3,n3};
        q5_i128 pred = 0, bnd = 0;
        for (int a=0;a<2;++a) if (w1[a]) for (int b=0;b<2;++b) if (w2[b]) for (int c=0;c<2;++c) if (w3[c]) {
            q5_i128 w = w1[a] * w2[b] * w3[c];
            int idc = tri_C_map[(e1[a] * NC + e2[b]) * NC + e3[c]];
            pred += w * mid[idc];
            bnd  += w * rad[idc];
        }
        if (q5_abs128((q5_i128)d1 * d2 * d3 * Tf[tid] - pred) > bnd) return false;
        ++checks;
    }
    return true;
}

static inline int64_t binom2(int c) { return (c >= 2) ? (int64_t)c * (c - 1) / 2 : 0; }
static inline int64_t binom3(int c) { return (c >= 3) ? (int64_t)c * (c - 1) * (c - 2) / 6 : 0; }
static inline int64_t binom4(int c) { return (c >= 4) ? (int64_t)c * (c - 1) * (c - 2) * (c - 3) / 24 : 0; }

int main(int argc, char** argv) {
    if (argc < 2) {
        cerr << "Usage: verify_quintet <snapshot.bin> [--override-pq P Q] [--perturb <target>] [--threads N]\n";
        return 1;
    }
    string path = argv[1];
    bool has_override = false;
    int64_t ov_P = -1, ov_Q = -1;
    string perturb = "";
    int nt = 14;
    for (int i = 2; i < argc; ++i) {
        string a = argv[i];
        if (a == "--override-pq") {
            if (i + 2 >= argc || !parse_i64_strict(argv[i + 1], ov_P) || !parse_i64_strict(argv[i + 2], ov_Q) ||
                ov_P <= 0 || ov_Q <= 0 || ov_P >= ov_Q || ov_Q > Q5_MAX_Q) {
                cerr << "ERROR: invalid --override-pq arguments\n"; return 2;
            }
            has_override = true;
            i += 2;
        } else if (a == "--perturb") {
            if (i + 1 >= argc) { cerr << "ERROR: missing --perturb target\n"; return 1; }
            perturb = argv[++i];
        } else if (a == "--threads") {
            int64_t ntv = 0;
            if (i + 1 >= argc || !parse_i64_strict(argv[i + 1], ntv) || ntv < 1 || ntv > 256) {
                cerr << "ERROR: invalid --threads value (must be 1..256)\n"; return 1;
            }
            nt = (int)ntv;
            ++i;
        } else { cerr << "Unknown arg: " << a << "\n"; return 1; }
    }

    FILE* fp = fopen(path.c_str(), "rb");
    if (!fp) { perror("fopen"); return 1; }
    uint8_t hdr_raw[sizeof(QuintetHeader)];
    if (fread(hdr_raw, 1, sizeof(hdr_raw), fp) != sizeof(hdr_raw)) { fclose(fp); cerr << "ERROR: short header\n"; return 2; }
    QuintetHeader h;
    memcpy(&h, hdr_raw, sizeof(h));
    // Validate the ORIGINAL on-disk header BEFORE any --override-pq replacement!
    uint64_t expected_bytes = 0;
    int vh = q5_compute_expected_bytes(h, expected_bytes);
    if (vh != 0) { fclose(fp); cerr << "ERROR: header/resource validation failed code=" << vh << "\n"; return 2; }
    if (has_override) {
        h.P = ov_P; h.Q = ov_Q;
        int vh2 = q5_validate_header(h);
        if (vh2 != 0) { fclose(fp); cerr << "ERROR: override header validation failed code=" << vh2 << "\n"; return 2; }
    }

    if (fseek(fp, 0, SEEK_END) != 0) { fclose(fp); cerr << "ERROR: fseek(END) failed\n"; return 2; }
    long long fsz_ll = ftell(fp);
    if (fsz_ll < 0 || (uint64_t)fsz_ll != expected_bytes) {
        fclose(fp);
        cerr << "ERROR: file size " << fsz_ll << " != expected " << expected_bytes << "\n";
        return 2;
    }
    if (fseek(fp, (long)sizeof(h), SEEK_SET) != 0) { fclose(fp); cerr << "ERROR: fseek(SET) failed\n"; return 2; }

    Q5Layout L; L.init(h);
    if ((uint64_t)L.total_bytes != expected_bytes) {
        fclose(fp); cerr << "ERROR: internal layout size mismatch\n"; return 2;
    }
    vector<int64_t> buf(L.total_i64);
    if (fread(buf.data(), sizeof(int64_t), L.total_i64, fp) != L.total_i64) {
        fclose(fp); cerr << "ERROR: short payload\n"; return 2;
    }
    fclose(fp);

    // Cryptographically authenticate the EXACT in-memory bytes (original header + payload) verified below!
    Q5Sha256 sha;
    sha.update(hdr_raw, sizeof(hdr_raw));
    sha.update(buf.data(), L.total_i64 * sizeof(int64_t));
    string input_sha256 = sha.hex_digest();

    auto check_struct_and_bounds = [&]() -> int {
        if (buf[L.off_g + L.g_off(L.G2.INF, L.G2.INF)] != 0) { cerr << "ERROR: g[INF,INF] != 0\n"; return 2; }
        for (size_t i = 0; i < L.n_r1; ++i)
            if (buf[L.off_r1 + i] < 0 || buf[L.off_r1 + i] > h.S_r / 2) { cerr << "ERROR: r1 out of [0, S_r/2]\n"; return 2; }
        for (size_t i = L.off_knots; i < L.off_e23_m; ++i)
            if (buf[i] < -MAX_W || buf[i] > MAX_W) { cerr << "ERROR: weight out of range\n"; return 2; }
        auto check_env_bounds = [&](size_t om, size_t orad, size_t n) -> bool {
            for (size_t i = 0; i < n; ++i) {
                if (buf[om + i] < -MAX_E || buf[om + i] > MAX_E) return false;
                if (buf[orad + i] < 0 || buf[orad + i] > MAX_E) return false;
            }
            return true;
        };
        if (!check_env_bounds(L.off_e23_m, L.off_e23_r, L.n_e23) ||
            !check_env_bounds(L.off_e34p_m, L.off_e34p_r, L.n_e34p) ||
            !check_env_bounds(L.off_e34s_m, L.off_e34s_r, L.n_e34s) ||
            !check_env_bounds(L.off_e45p_m, L.off_e45p_r, L.n_e45p) ||
            !check_env_bounds(L.off_e45bp_m, L.off_e45bp_r, L.n_e45bp) ||
            !check_env_bounds(L.off_e45bs_m, L.off_e45bs_r, L.n_e45bs) ||
            !check_env_bounds(L.off_e45s4_m, L.off_e45s4_r, L.n_e45s4)) {
            cerr << "ERROR: envelope radius negative or magnitude out of range\n"; return 2;
        }
        return 0;
    };

    // === 1a. Validate base snapshot BEFORE applying any optional negative-control mutation ===
    if (check_struct_and_bounds() != 0) return 2;

    // Optional single-entry mutation control (for fail-closed negative tests ONLY!)
    if (!perturb.empty()) {
        int64_t delta = (h.mode == 1) ? -INT64_C(1000000000) : INT64_C(1000000000);
        auto add_checked = [&](size_t idx, int64_t d) {
            if (idx >= buf.size()) { cerr << "ERROR: perturb index out of bounds\n"; std::exit(2); }
            buf[idx] = q5_to_i64((q5_i128)buf[idx] + (q5_i128)d, "perturb");
        };
        if      (perturb == "f")       add_checked(L.off_f + (size_t)L.f_id(2, 1), delta);
        else if (perturb == "g")       add_checked(L.off_g + L.g_off(0, 0), delta);
        else if (perturb == "h3")      add_checked(L.off_h3 + 0, delta);
        else if (perturb == "p4")      add_checked(L.off_p4 + 0, delta);
        else if (perturb == "s3")      add_checked(L.off_s3 + 0, delta);
        else if (perturb == "p5")      add_checked(L.off_p5 + 0, delta);
        else if (perturb == "b5")      add_checked(L.off_b5 + 0, delta);
        else if (perturb == "s4")      add_checked(L.off_s4 + 0, delta);
        else if (perturb == "r1")      add_checked(L.off_r1 + 0, INT64_C(50000000));
        else if (perturb == "knots") {
            int fy_u = -1;
            for (int u = 0; u < L.G2.Nin_ex; ++u) {
                int Du = L.G2.D[u], cu = L.G2.C[u];
                bool is_m3 = (Du < h.d3m && cu >= 2 && cu <= h.K3m);
                if (cu >= 2 && !is_m3) { fy_u = u; break; }
            }
            if (fy_u < 0) { cerr << "ERROR: no FY state available for knots perturb\n"; return 2; }
            add_checked(L.off_knots + (size_t)fy_u * (size_t)h.nknots + 2, delta);
        }
        else if (perturb == "env23")   add_checked(L.off_e23_r + (size_t)0 * L.G3.N + L.G3.t_id(1), -1);
        else if (perturb == "env34p")  add_checked(L.off_e34p_r + (size_t)0 * L.G4.N + L.G4.t_id(1), -1);
        else if (perturb == "env34s")  add_checked(L.off_e34s_r + Q5Layout::sym_off(L.G4.t_id(1), L.G4.t_id(1)), -1);
        else if (perturb == "env45p")  add_checked(L.off_e45p_r + (size_t)0 * L.G5.N + L.G5.t_id(1), -1);
        else if (perturb == "env45bp") add_checked(L.off_e45bp_r + (size_t)L.G5.t_id(1) * L.G5.N + L.G5.t_id(1), -1);
        else if (perturb == "env45bs") add_checked(L.off_e45bs_r + Q5Layout::sym_off(L.G5.t_id(1), L.G5.t_id(1)), -1);
        else if (perturb == "env45s4") add_checked(L.off_e45s4_r + L.tri_G5_map[((size_t)L.G5.t_id(1)*L.G5.N + L.G5.t_id(1))*L.G5.N + L.G5.t_id(1)], -1);
        else if (perturb == "ginf")    buf[L.off_g + L.g_off(L.G2.INF, L.G2.INF)] = 1;
        else { cerr << "Unknown perturb target: " << perturb << "\n"; return 1; }
        // === 1b. Re-validate structural & magnitude bounds AFTER mutation ===
        if (check_struct_and_bounds() != 0) return 2;
    }

    const int64_t* r1    = &buf[L.off_r1];
    const int64_t* knots = &buf[L.off_knots];
    const int64_t* f     = &buf[L.off_f];
    const int64_t* g     = &buf[L.off_g];
    const int64_t* h3    = &buf[L.off_h3];
    const int64_t* p4    = &buf[L.off_p4];
    const int64_t* s3    = &buf[L.off_s3];
    const int64_t* p5    = &buf[L.off_p5];
    const int64_t* b5    = &buf[L.off_b5];
    const int64_t* s4    = &buf[L.off_s4];
    const int64_t* e23_m = &buf[L.off_e23_m],   *e23_r = &buf[L.off_e23_r];
    const int64_t* e34p_m= &buf[L.off_e34p_m],  *e34p_r= &buf[L.off_e34p_r];
    const int64_t* e34s_m= &buf[L.off_e34s_m],  *e34s_r= &buf[L.off_e34s_r];
    const int64_t* e45p_m= &buf[L.off_e45p_m],  *e45p_r= &buf[L.off_e45p_r];
    const int64_t* e45bp_m=&buf[L.off_e45bp_m], *e45bp_r=&buf[L.off_e45bp_r];
    const int64_t* e45bs_m=&buf[L.off_e45bs_m], *e45bs_r=&buf[L.off_e45bs_r];
    const int64_t* e45s4_m=&buf[L.off_e45s4_m], *e45s4_r=&buf[L.off_e45s4_r];

    int sgn = (h.mode == 1) ? 1 : -1;
    int NG = L.G2.N, N3 = L.G3.N, N4 = L.G4.N, N5 = L.G5.N;
    size_t g4_p2 = (size_t)N4 * (N4 + 1) / 2, g5_p2 = (size_t)N5 * (N5 + 1) / 2;
    int J_g = L.G2.idx[2][1], J_e3 = L.G3.idx[2][1], J_e4 = L.G4.idx[2][1], J_e5 = L.G5.idx[2][1];
    int eff_nt = (L.n_p4 < 500000) ? 1 : nt;

    // === 2. Verify ALL 7 cross-level envelopes in exact q5_i128 ===
    size_t c_env1 = 0, c_env2 = 0, c_env3 = 0;
    int env_fail = 0;
    #pragma omp parallel for schedule(static) num_threads(eff_nt) reduction(+:c_env1,env_fail) if(eff_nt > 1)
    for (int m = 0; m < L.M3.N; ++m) {
        int gu = L.G2.idx[L.M3.D[m]][L.M3.C[m]], gup = L.G2.up(gu);
        vector<int64_t> dg(NG);
        for (int w = 0; w < NG; ++w) dg[w] = g[L.g_off(gup, w)] - g[L.g_off(gu, w)];
        if (!check_env_1d(L.PM23, L.G2, L.G3, dg.data(), &e23_m[(size_t)m * N3], &e23_r[(size_t)m * N3], c_env1)) env_fail += 1;
    }
    if (env_fail) { cerr << "ERROR: env23 verification failed\n"; return 2; }

    #pragma omp parallel for schedule(static) num_threads(eff_nt) reduction(+:c_env1,env_fail) if(eff_nt > 1)
    for (int pid = 0; pid < L.NP4_ord; ++pid) {
        auto [b, cc] = L.p4_ord[pid];
        int b_e3 = L.G3.idx[L.M4v.D[b]][L.M4v.C[b]], b_e3up = L.G3.up(b_e3);
        int c_m3 = L.M3.idx[L.M4v.D[cc]][L.M4v.C[cc]];
        vector<int64_t> dh(N3);
        for (int d3 = 0; d3 < N3; ++d3) dh[d3] = h3[L.h3_off(c_m3, b_e3up, d3)] - h3[L.h3_off(c_m3, b_e3, d3)];
        if (!check_env_1d(L.PM34, L.G3, L.G4, dh.data(), &e34p_m[(size_t)pid * N4], &e34p_r[(size_t)pid * N4], c_env1)) env_fail += 1;
    }
    if (env_fail) { cerr << "ERROR: env34p verification failed\n"; return 2; }

    #pragma omp parallel for schedule(static) num_threads(eff_nt) reduction(+:c_env2,env_fail) if(eff_nt > 1)
    for (int o = 0; o < L.MS3.N; ++o) {
        int om3 = L.M3.idx[L.MS3.D[o]][L.MS3.C[o]], om3_up = L.M3.up(om3);
        vector<int64_t> Mf((size_t)N3 * N3);
        for (int u = 0; u < N3; ++u) for (int v = 0; v < N3; ++v)
            Mf[(size_t)u * N3 + v] = ((om3_up >= 0) ? h3[L.h3_off(om3_up, u, v)] : 0) - h3[L.h3_off(om3, u, v)];
        if (!check_env_2d(L.G3, L.G4, Mf, &e34s_m[(size_t)o * g4_p2], &e34s_r[(size_t)o * g4_p2], true, c_env2)) env_fail += 1;
    }
    if (env_fail) { cerr << "ERROR: env34s verification failed\n"; return 2; }

    #pragma omp parallel for schedule(static) num_threads(eff_nt) reduction(+:c_env1,env_fail) if(eff_nt > 1)
    for (int pid5 = 0; pid5 < L.NP5_ord; ++pid5) {
        auto [b, c, d] = L.p5_ord[pid5];
        int b_e4 = L.G4.idx[L.I5_2.D[b]][L.I5_2.C[b]], b_e4up = L.G4.up(b_e4);
        int c_m4 = L.M4v.idx[L.I5_2.D[c]][L.I5_2.C[c]], d_m4 = L.M4v.idx[L.I5_2.D[d]][L.I5_2.C[d]];
        vector<int64_t> dp(N4);
        for (int e4 = 0; e4 < N4; ++e4) dp[e4] = p4[L.p4_off(c_m4, d_m4, b_e4up, e4)] - p4[L.p4_off(c_m4, d_m4, b_e4, e4)];
        if (!check_env_1d(L.PM45, L.G4, L.G5, dp.data(), &e45p_m[(size_t)pid5 * N5], &e45p_r[(size_t)pid5 * N5], c_env1)) env_fail += 1;
    }
    if (env_fail) { cerr << "ERROR: env45p verification failed\n"; return 2; }

    #pragma omp parallel for schedule(static) num_threads(eff_nt) reduction(+:c_env2,env_fail) if(eff_nt > 1)
    for (int bid = 0; bid < L.NB5_int; ++bid) {
        auto [o, c] = L.b5_int[bid];
        int Do = L.I5_3.D[o], co = L.I5_3.C[o], Dc = L.I5_2c.D[c], cc = L.I5_2c.C[c];
        int o_m4 = L.M4v.idx[Do][co], o_m4up = L.M4v.up(o_m4), c_m4 = L.M4v.idx[Dc][cc];
        bool has_up = (o_m4up >= 0 && Do + 1 + Dc <= h.SP4);
        vector<int64_t> Mf((size_t)N4 * N4);
        for (int u = 0; u < N4; ++u) for (int v = 0; v < N4; ++v)
            Mf[(size_t)u * N4 + v] = (has_up ? p4[L.p4_off(o_m4up, c_m4, u, v)] : 0) - p4[L.p4_off(o_m4, c_m4, u, v)];
        if (!check_env_2d(L.G4, L.G5, Mf, &e45bp_m[(size_t)bid * N5 * N5], &e45bp_r[(size_t)bid * N5 * N5], false, c_env2)) env_fail += 1;

        int o_ms3 = L.MS3.idx[Do][co], c_e4 = L.G4.idx[Dc][cc], c_e4up = L.G4.up(c_e4);
        for (int u = 0; u < N4; ++u) for (int v = 0; v < N4; ++v)
            Mf[(size_t)u * N4 + v] = s3[L.s3_off(o_ms3, u, v, c_e4up)] - s3[L.s3_off(o_ms3, u, v, c_e4)];
        if (!check_env_2d(L.G4, L.G5, Mf, &e45bs_m[(size_t)bid * g5_p2], &e45bs_r[(size_t)bid * g5_p2], true, c_env2)) env_fail += 1;
    }
    if (env_fail) { cerr << "ERROR: env45bp/env45bs verification failed\n"; return 2; }

    #pragma omp parallel for schedule(static) num_threads(eff_nt) reduction(+:c_env3,env_fail) if(eff_nt > 1)
    for (int o = 0; o < L.MS4.N; ++o) {
        int oms3 = L.MS3.idx[L.MS4.D[o]][L.MS4.C[o]], oms3_up = L.MS3.up(oms3);
        vector<int64_t> Tf(L.N_G4_tri);
        for (int tid = 0; tid < L.N_G4_tri; ++tid)
            Tf[tid] = ((oms3_up >= 0) ? s3[(size_t)oms3_up * L.N_G4_tri + tid] : 0) - s3[(size_t)oms3 * L.N_G4_tri + tid];
        if (!check_env_3d(L.G4, L.G5, L.tri_G4, L.tri_G5_map, Tf, &e45s4_m[(size_t)o * L.N_G5_tri], &e45s4_r[(size_t)o * L.N_G5_tri], c_env3)) env_fail += 1;
    }
    if (env_fail) { cerr << "ERROR: env45s4 verification failed\n"; return 2; }

    // === 3. Verify ALL finite rows of f (D < DF, including exact analytic seam at D = DF-1!) ===
    int64_t cap_f_in = INT64_MIN, cap_f_tab = INT64_MIN;
    int64_t f_J = f[L.f_id(2, 1)];
    size_t rows_f_in = 0, rows_f_tab = 0;
    for (int D = 2; D < h.DF; ++D) {
        for (int c = 0; c <= D; ++c) {
            int Lv = D - c;
            int64_t f0 = f[L.f_id(D, c)];
            // Gain jump df = f_up - f0 as exact rational (df_n / df_d):
            q5_i128 df_n, df_d;
            if (D + 1 < h.DF) { df_n = f[L.f_id(D + 1, c)] - f0; df_d = 1; }
            else {
                // Exact analytic f_ana(DF, c) at seam D+1 == DF!
                q5_i128 DF = h.DF, L_up = DF - c;
                if (h.mode == 1) {
                    df_n = (q5_i128)h.f_p0 * DF + (q5_i128)h.f_p1 * L_up - (q5_i128)f0 * DF;
                    df_d = DF;
                } else {
                    df_n = (q5_i128)h.f_p0 * L_up * DF - (q5_i128)h.f_p1 * L_up * L_up - (q5_i128)h.f_p2 * L_up - (q5_i128)f0 * DF * DF;
                    df_d = DF * DF;
                }
            }
            // Promotion & newborn edge g((D, c+1), J)
            q5_i128 bud_n = 0, bud_d = 1;
            if (Lv > 0) {
                int64_t f_bd = f[L.f_id(D, c + 1)];
                if (D < h.d2) {
                    int id = (c + 1 <= h.K2) ? L.G2.idx[D][c + 1] : L.G2.ovf_idx[D];
                    bud_n = (q5_i128)Lv * (f_bd - f0 + f_J + g[L.g_off(id, J_g)]);
                    bud_d = D;
                } else {
                    q5_i128 g_inf = g[L.g_off(L.G2.INF, J_g)];
                    q5_i128 g_t   = g[L.g_off(L.G2.t_id(c + 1), J_g)];
                    q5_i128 nb_n  = (q5_i128)D * (f_bd - f0 + f_J) + (q5_i128)(D - h.d2) * g_inf + (q5_i128)h.d2 * g_t;
                    bud_n = (q5_i128)Lv * nb_n;
                    bud_d = (q5_i128)D * D;
                }
            }
            q5_i128 acc = 0;
            add_rat(acc, -(q5_i128)sgn * h.P * f0, h.Q);
            if (Lv > 0) add_rat(acc, (q5_i128)sgn * bud_n, bud_d);
            bool centered = (D < h.d2 && c >= 1 && c <= h.K2);
            if (centered) {
                q5_i128 ur_n = (q5_i128)Lv * h.S_r + (q5_i128)c * r1[L.G2.idx[D][c]];
                add_rat(acc, (q5_i128)sgn * ur_n * df_n, h.S_r * df_d);
                int64_t rcap = q5_to_i64(q5_ceil_div(acc, K_SUB), "rcap");
                if (rcap > cap_f_in) cap_f_in = rcap;
                ++rows_f_in;
            } else {
                // Extremize r in [Lv, Lv + c/2] = [2*Lv/2, (2*Lv + c)/2]
                q5_i128 sdf = (q5_i128)sgn * df_n;
                q5_i128 r_num = (sdf > 0) ? (2 * Lv + c) : (2 * Lv);
                if (r_num > 0) add_rat(acc, sdf * r_num, 2 * df_d);
                int64_t rcap = q5_to_i64(q5_ceil_div(acc, K_SUB), "rcap");
                if (rcap > cap_f_tab) cap_f_tab = rcap;
                ++rows_f_tab;
            }
        }
    }

    // === 4. Verify ALL rows of g (u <= v, excluding (INF, INF)) ===
    int64_t cap_g = INT64_MIN;
    size_t rows_g = 0;
    vector<int64_t> dLf(L.G2.Nin_ex);
    for (int u = 0; u < L.G2.Nin_ex; ++u) dLf[u] = f[L.f_id(L.G2.D[u] + 1, L.G2.C[u])] - f[L.f_id(L.G2.D[u], L.G2.C[u])];
    #pragma omp parallel for schedule(dynamic) num_threads(eff_nt) reduction(max:cap_g) reduction(+:rows_g) if(eff_nt > 1)
    for (int u = 0; u < NG; ++u) {
        int uk = L.G2.kind[u], Du = L.G2.D[u], cu = L.G2.C[u], Lu = L.G2.L[u];
        bool u_in = L.G2.is_inner(u);
        int u_up = L.G2.up(u), u_bud = L.G2.bud(u);
        bool u_m3 = (uk == 0 && Du < h.d3m && cu >= 2 && cu <= h.K3m);
        bool u_fy = (uk == 0 && cu >= 2 && !u_m3);

        for (int v = u; v < NG; ++v) {
            if (u == L.G2.INF && v == L.G2.INF) continue;
            int vk = L.G2.kind[v], Dv = L.G2.D[v], cv = L.G2.C[v], Lv = L.G2.L[v];
            bool v_in = L.G2.is_inner(v);
            int v_up = L.G2.up(v), v_bud = L.G2.bud(v);
            bool v_m3 = (vk == 0 && Dv < h.d3m && cv >= 2 && cv <= h.K3m);
            bool v_fy = (vk == 0 && cv >= 2 && !v_m3);

            int64_t val = g[L.g_off(u, v)];
            int64_t gu_up = g[L.g_off(u_up, v)], gv_up = g[L.g_off(u, v_up)];
            int64_t gu_bd = g[L.g_off(u_bud, v)], gv_bd = g[L.g_off(u, v_bud)];
            q5_i128 acc = 0;
            add_rat(acc, -(q5_i128)sgn * h.P * val, h.Q);

            // Incoming transfers from f + FY primal incoming + internal transitions on u and v
            auto do_end = [&](int x, int xk, int Dx, int cx, int Lx, int x_bud, bool x_m3, bool x_fy,
                              int64_t gx_up, int64_t gx_bd, int y, bool y_in, int Dy) {
                if (xk == 0) {
                    // f->g source: (inv_y - r1[x]) * dLf[x]
                    add_src_env(acc, sgn, L.G2, y, r1[x], h.S_r, dLf[x], 0);
                    // Own leaf + neighbor + centered M3 gain
                    q5_i128 dg = gx_up - val;
                    if (Lx > 0) add_rat(acc, (q5_i128)sgn * Lx * dg, 1);
                    if (y_in) add_rat(acc, (q5_i128)sgn * dg, Dy);
                    if (x_m3) add_rat(acc, (q5_i128)sgn * (cx - 1) * r1[x] * dg, h.S_r);
                    // FY incoming A(inv_y) and dual (cx-1)*A*(dg)
                    if (x_fy) {
                        const int64_t* kx = &knots[(size_t)x * h.nknots];
                        q5_i128 fn, fd; eval_fy_primal(h, L.G2, y, kx, fn, fd);
                        add_rat(acc, (q5_i128)sgn * (cx - 1) * fn, fd);
                        q5_i128 best = -((q5_i128)1 << 120);
                        for (int kk = 0; kk < h.nknots; ++kk) {
                            q5_i128 cval = (q5_i128)sgn * (h.knot_num[kk] * dg - h.Q_K * kx[kk]);
                            if (cval > best) best = cval;
                        }
                        add_rat(acc, (q5_i128)(cx - 1) * best, h.Q_K);
                    }
                    // Promotion + H3 newborn
                    if (Lx > 0) {
                        add_rat(acc, (q5_i128)sgn * Lx * (gx_bd - val), Dx);
                        if (L.G2.kind[x_bud] == 0 && L.G2.D[x_bud] < h.d3m && L.G2.C[x_bud] >= 2 && L.G2.C[x_bud] <= h.K3m) {
                            int m = L.M3.idx[L.G2.D[x_bud]][L.G2.C[x_bud]], cy, ny, dy;
                            L.G3.proj_int(L.G2, y, cy, ny, dy);
                            q5_i128 nb_n = (q5_i128)ny * h3[L.h3_off(m, cy, J_e3)] + (q5_i128)(dy - ny) * h3[L.h3_off(m, L.G3.INF, J_e3)];
                            add_rat(acc, (q5_i128)sgn * Lx * nb_n, (q5_i128)Dx * dy);
                        }
                    }
                } else if (xk == 3) {
                    if (y_in) add_rat(acc, (q5_i128)sgn * (gx_up - val), Dy);
                    add_ovf_end(acc, sgn, L.G2, x, gx_up, val);
                } else if (xk == 1) {
                    int64_t d0 = y_in ? (g[L.g_off(L.G2.INF, L.G2.up(y))] - g[L.g_off(L.G2.INF, y)]) : 0;
                    int64_t y_up_v = y_in ? g[L.g_off(x, L.G2.up(y))] : 0;
                    add_tail_end(acc, sgn, L.G2, x, y_in ? Dy : 0, val, gx_up, gx_bd, d0, y_up_v);
                }
            };
            do_end(u, uk, Du, cu, Lu, u_bud, u_m3, u_fy, gu_up, gu_bd, v, v_in, Dv);
            do_end(v, vk, Dv, cv, Lv, v_bud, v_m3, v_fy, gv_up, gv_bd, u, u_in, Du);

            int64_t rcap = q5_to_i64(q5_ceil_div(acc, K_SUB), "rcap");
            if (rcap > cap_g) cap_g = rcap;
            ++rows_g;
        }
    }

    // === 5. Verify ALL rows of H3 ===
    vector<int64_t> cap_h3_m(L.M3.N, INT64_MIN);
    size_t rows_h3 = 0;
    #pragma omp parallel for schedule(dynamic) num_threads(eff_nt) reduction(+:rows_h3) if(eff_nt > 1)
    for (int m = 0; m < L.M3.N; ++m) {
        int rm = L.M3.D[m], cm = L.M3.C[m], Lm = rm - cm;
        int m_up = L.M3.up(m), m_bud = L.M3.bud(m);
        int64_t r_m = r1[L.G2.idx[rm][cm]];
        int m_m4 = (rm < h.d4i && cm <= h.K4i) ? L.M4v.idx[rm][cm] : -1;
        bool m_in_s3 = (rm < h.dS3 && cm >= 3 && cm <= h.KS3);
        int64_t loc_cap = INT64_MIN;

        for (int a = 0; a < N3; ++a) {
            int ak = L.G3.kind[a], Da = L.G3.D[a], ca = L.G3.C[a], La = L.G3.L[a];
            bool a_in = L.G3.is_inner(a);
            int a_up = L.G3.up(a), a_bud = L.G3.bud(a);
            bool a_p4 = (m_m4 >= 0 && ak == 0 && Da < h.d4i && ca >= 2 && ca <= h.K4i && Da + rm <= h.SP4);
            int64_t ma = (ak == 0) ? r1[L.G2.idx[Da][ca]] : 0;

            for (int b = a; b < N3; ++b) {
                int bk = L.G3.kind[b], Db = L.G3.D[b], cb = L.G3.C[b], Lb = L.G3.L[b];
                bool b_in = L.G3.is_inner(b);
                int b_up = L.G3.up(b), b_bud = L.G3.bud(b);
                bool b_p4 = (m_m4 >= 0 && bk == 0 && Db < h.d4i && cb >= 2 && cb <= h.K4i && Db + rm <= h.SP4);
                int64_t mb = (bk == 0) ? r1[L.G2.idx[Db][cb]] : 0;

                int64_t val = h3[L.h3_off(m, a, b)];
                int64_t ha_up = h3[L.h3_off(m, a_up, b)], hb_up = h3[L.h3_off(m, a, b_up)];
                int64_t hm_up = (m_up >= 0) ? h3[L.h3_off(m_up, a, b)] : 0;
                int64_t ha_bd = h3[L.h3_off(m, a_bud, b)], hb_bd = h3[L.h3_off(m, a, b_bud)];
                int64_t hm_bd = (m_bud >= 0) ? h3[L.h3_off(m_bud, a, b)] : 0;

                q5_i128 acc = 0;
                add_rat(acc, -(q5_i128)sgn * h.P * val, h.Q);
                add_src_env(acc, sgn, L.G3, a, r_m, h.S_r, e23_m[(size_t)m * N3 + b], e23_r[(size_t)m * N3 + b]);
                add_src_env(acc, sgn, L.G3, b, r_m, h.S_r, e23_m[(size_t)m * N3 + a], e23_r[(size_t)m * N3 + a]);

                // Middle m transitions
                int64_t dm = hm_up - val;
                if (Lm > 0) add_rat(acc, (q5_i128)sgn * Lm * dm, 1);
                if (a_in) add_rat(acc, (q5_i128)sgn * dm, Da);
                if (b_in) add_rat(acc, (q5_i128)sgn * dm, Db);
                if (m_in_s3) add_rat(acc, (q5_i128)sgn * (cm - 2) * r_m * dm, h.S_r);
                else if (cm > 2) add_ext_rate(acc, sgn, dm, 0, 1, cm - 2, 2);
                if (Lm > 0) {
                    add_rat(acc, (q5_i128)sgn * Lm * (hm_bd - val), rm);
                    if (rm < h.dS3 && cm + 1 >= 3 && cm + 1 <= h.KS3) {
                        int mo = L.MS3.idx[rm][cm + 1], c1, n1, d1, c2, n2, d2;
                        L.G4.proj_int(L.G3, a, c1, n1, d1); L.G4.proj_int(L.G3, b, c2, n2, d2);
                        int e1[2]={L.G4.INF,c1}, e2[2]={L.G4.INF,c2}; q5_i128 w1[2]={d1-n1,n1}, w2[2]={d2-n2,n2};
                        q5_i128 nb = 0;
                        for (int p=0;p<2;++p) if (w1[p]) for (int q=0;q<2;++q) if (w2[q])
                            nb += w1[p] * w2[q] * s3[L.s3_off(mo, e1[p], e2[q], J_e4)];
                        add_rat(acc, (q5_i128)sgn * Lm * nb, (q5_i128)rm * d1 * d2);
                    }
                }
                // Endpoints a and b
                auto do_h3_end = [&](int x, int xk, int Dx, int cx, int Lx, int x_bud, bool x_p4, int64_t rx,
                                     int64_t hx_up, int64_t hx_bd, int y) {
                    if (xk == 0) {
                        int64_t dx = hx_up - val;
                        if (Lx > 0) add_rat(acc, (q5_i128)sgn * Lx * dx, 1);
                        add_rat(acc, (q5_i128)sgn * dx, rm);
                        if (x_p4) add_rat(acc, (q5_i128)sgn * (cx - 1) * rx * dx, h.S_r);
                        else if (cx > 1) add_ext_rate(acc, sgn, dx, 0, 1, cx - 1, 2);
                        if (Lx > 0) {
                            add_rat(acc, (q5_i128)sgn * Lx * (hx_bd - val), Dx);
                            if (L.G3.kind[x_bud] == 0 && Dx < h.d4i && cx + 1 >= 2 && cx + 1 <= h.K4i && m_m4 >= 0 && Dx + rm <= h.SP4) {
                                int mx = L.M4v.idx[Dx][cx + 1], cy, ny, dy;
                                L.G4.proj_int(L.G3, y, cy, ny, dy);
                                q5_i128 nb = (q5_i128)ny * p4[L.p4_off(mx, m_m4, J_e4, cy)] + (q5_i128)(dy - ny) * p4[L.p4_off(mx, m_m4, J_e4, L.G4.INF)];
                                add_rat(acc, (q5_i128)sgn * Lx * nb, (q5_i128)Dx * dy);
                            }
                        }
                    } else if (xk == 3) {
                        add_rat(acc, (q5_i128)sgn * (hx_up - val), rm);
                        add_ovf_end(acc, sgn, L.G3, x, hx_up, val);
                    } else if (xk == 1) {
                        int64_t d0 = ((m_up >= 0) ? h3[L.h3_off(m_up, L.G3.INF, y)] : 0) - h3[L.h3_off(m, L.G3.INF, y)];
                        add_tail_end(acc, sgn, L.G3, x, rm, val, hx_up, hx_bd, d0, hm_up);
                    }
                };
                do_h3_end(a, ak, Da, ca, La, a_bud, a_p4, ma, ha_up, ha_bd, b);
                do_h3_end(b, bk, Db, cb, Lb, b_bud, b_p4, mb, hb_up, hb_bd, a);

                int64_t rcap = q5_to_i64(q5_ceil_div(acc, K_SUB), "rcap");
                if (rcap > loc_cap) loc_cap = rcap;
                ++rows_h3;
            }
        }
        cap_h3_m[m] = loc_cap;
    }

    // === 6. Verify ALL rows of P4 (canonical b <= c; if b == c, a <= d) ===
    vector<int64_t> cap_p4_ord(L.NP4_ord, INT64_MIN);
    size_t rows_p4 = 0;
    #pragma omp parallel for schedule(dynamic) num_threads(eff_nt) reduction(+:rows_p4) if(eff_nt > 1)
    for (int pid = 0; pid < L.NP4_ord; ++pid) {
        auto [b, c] = L.p4_ord[pid];
        if (b > c) continue; // checked canonically with b <= c!
        int rb = L.M4v.D[b], cb = L.M4v.C[b], Lb = rb - cb;
        int rc = L.M4v.D[c], cc = L.M4v.C[c], Lc = rc - cc;
        int b_up = L.M4v.up(b), c_up = L.M4v.up(c), b_bud = L.M4v.bud(b), c_bud = L.M4v.bud(c);
        int pid_rev = L.p4_ord_map[(size_t)c * L.M4v.N + b];
        bool has_bup = (b_up >= 0 && rb + 1 + rc <= h.SP4);
        bool has_cup = (c_up >= 0 && rb + rc + 1 <= h.SP4);
        int64_t r_b = r1[L.G2.idx[rb][cb]], r_c = r1[L.G2.idx[rc][cc]];
        bool cov_b = (rb < h.d5o && cb >= 3 && cb <= h.K5o && rc < h.d5c && cc <= h.K5c && rb + rc <= h.SB5);
        bool cov_c = (rc < h.d5o && cc >= 3 && cc <= h.K5o && rb < h.d5c && cb <= h.K5c && rb + rc <= h.SB5);
        bool bc_p5 = (rb < h.d5i && cb <= h.K5i && rc < h.d5i && cc <= h.K5i);
        int64_t loc_cap = INT64_MIN;

        for (int a = 0; a < N4; ++a) {
            int ak = L.G4.kind[a], Da = L.G4.D[a], ca = L.G4.C[a], La = L.G4.L[a];
            bool a_in = L.G4.is_inner(a);
            int a_up = L.G4.up(a), a_bud = L.G4.bud(a);
            bool cov_a = (bc_p5 && ak == 0 && Da < h.d5i && ca >= 2 && ca <= h.K5i && Da + rb <= h.SP4 && Da + rb + rc <= h.SP5);
            int64_t r_a = (ak == 0) ? r1[L.G2.idx[Da][ca]] : 0;

            int d_start = (b == c) ? a : 0;
            for (int d = d_start; d < N4; ++d) {
                int dk = L.G4.kind[d], Dd = L.G4.D[d], cd = L.G4.C[d], Ld = L.G4.L[d];
                bool d_in = L.G4.is_inner(d);
                int d_up = L.G4.up(d), d_bud = L.G4.bud(d);
                bool cov_d = (bc_p5 && dk == 0 && Dd < h.d5i && cd >= 2 && cd <= h.K5i && Dd + rc <= h.SP4 && rb + rc + Dd <= h.SP5);
                int64_t r_d = (dk == 0) ? r1[L.G2.idx[Dd][cd]] : 0;

                int64_t val = p4[L.p4_off(b, c, a, d)];
                int64_t pa_up = p4[L.p4_off(b, c, a_up, d)], pd_up = p4[L.p4_off(b, c, a, d_up)];
                int64_t pb_up = has_bup ? p4[L.p4_off(b_up, c, a, d)] : 0;
                int64_t pc_up = has_cup ? p4[L.p4_off(b, c_up, a, d)] : 0;
                int64_t pa_bd = p4[L.p4_off(b, c, a_bud, d)], pd_bd = p4[L.p4_off(b, c, a, d_bud)];
                int64_t pb_bd = (b_bud >= 0) ? p4[L.p4_off(b_bud, c, a, d)] : 0;
                int64_t pc_bd = (c_bud >= 0) ? p4[L.p4_off(b, c_bud, a, d)] : 0;

                q5_i128 acc = 0;
                add_rat(acc, -(q5_i128)sgn * h.P * val, h.Q);
                add_src_env(acc, sgn, L.G4, a, r_b, h.S_r, e34p_m[(size_t)pid * N4 + d], e34p_r[(size_t)pid * N4 + d]);
                add_src_env(acc, sgn, L.G4, d, r_c, h.S_r, e34p_m[(size_t)pid_rev * N4 + a], e34p_r[(size_t)pid_rev * N4 + a]);

                // Internal b and c
                auto do_mid = [&](int rx, int cx, int Lx, int ry, int cy, bool cov_x, int64_t r_x,
                                  int64_t px_up, int64_t px_bd, bool near_in, int D_near, int e_near, int e_far) {
                    int64_t dx = px_up - val;
                    if (Lx > 0) add_rat(acc, (q5_i128)sgn * Lx * dx, 1);
                    if (near_in) add_rat(acc, (q5_i128)sgn * dx, D_near);
                    add_rat(acc, (q5_i128)sgn * dx, ry);
                    if (cov_x) add_rat(acc, (q5_i128)sgn * (cx - 2) * r_x * dx, h.S_r);
                    else if (cx > 2) add_ext_rate(acc, sgn, dx, 0, 1, cx - 2, 2);
                    if (Lx > 0) {
                        add_rat(acc, (q5_i128)sgn * Lx * (px_bd - val), rx);
                        if (rx < h.d5o && cx + 1 >= 3 && cx + 1 <= h.K5o && ry < h.d5c && cy <= h.K5c && rx + ry <= h.SB5) {
                            int bid = L.b5_int_map[(size_t)L.I5_3.idx[rx][cx + 1] * L.I5_2c.N + L.I5_2c.idx[ry][cy]];
                            int c1, n1, d1, c2, n2, d2;
                            L.G5.proj_int(L.G4, e_near, c1, n1, d1); L.G5.proj_int(L.G4, e_far, c2, n2, d2);
                            int u1[2]={L.G5.INF,c1}, u2[2]={L.G5.INF,c2}; q5_i128 w1[2]={d1-n1,n1}, w2[2]={d2-n2,n2};
                            q5_i128 nb = 0;
                            for (int p=0;p<2;++p) if (w1[p]) for (int q=0;q<2;++q) if (w2[q])
                                nb += w1[p] * w2[q] * b5[L.b5_off(bid, u1[p], J_e5, u2[q])];
                            add_rat(acc, (q5_i128)sgn * Lx * nb, (q5_i128)rx * d1 * d2);
                        }
                    }
                };
                do_mid(rb, cb, Lb, rc, cc, cov_b, r_b, pb_up, pb_bd, a_in, Da, a, d);
                do_mid(rc, cc, Lc, rb, cb, cov_c, r_c, pc_up, pc_bd, d_in, Dd, d, a);

                // Endpoints a and d
                auto do_end4 = [&](int x, int xk, int Dx, int cx, int Lx, int x_bud, bool cov_x, int64_t r_x,
                                   int64_t px_up, int64_t px_bd, int r_n, int c_n, int r_f, int c_f,
                                   int e_far, int64_t pn_up, int64_t d0_n) {
                    if (xk == 0) {
                        int64_t dx = px_up - val;
                        if (Lx > 0) add_rat(acc, (q5_i128)sgn * Lx * dx, 1);
                        add_rat(acc, (q5_i128)sgn * dx, r_n);
                        if (cov_x) add_rat(acc, (q5_i128)sgn * (cx - 1) * r_x * dx, h.S_r);
                        else if (cx > 1) add_ext_rate(acc, sgn, dx, 0, 1, cx - 1, 2);
                        if (Lx > 0) {
                            add_rat(acc, (q5_i128)sgn * Lx * (px_bd - val), Dx);
                            if (L.G4.kind[x_bud] == 0 && Dx < h.d5i && cx + 1 >= 2 && cx + 1 <= h.K5i &&
                                r_n < h.d5i && c_n <= h.K5i && r_f < h.d5i && c_f <= h.K5i &&
                                Dx + r_n <= h.SP4 && r_n + r_f <= h.SP4 && Dx + r_n + r_f <= h.SP5) {
                                int ix = L.I5_2.idx[Dx][cx + 1], in = L.I5_2.idx[r_n][c_n], ifr = L.I5_2.idx[r_f][c_f];
                                int cf, nf, df; L.G5.proj_int(L.G4, e_far, cf, nf, df);
                                q5_i128 nb = (q5_i128)nf * p5[L.p5_off(ix, in, ifr, J_e5, cf)] + (q5_i128)(df - nf) * p5[L.p5_off(ix, in, ifr, J_e5, L.G5.INF)];
                                add_rat(acc, (q5_i128)sgn * Lx * nb, (q5_i128)Dx * df);
                            }
                        }
                    } else if (xk == 3) {
                        add_rat(acc, (q5_i128)sgn * (px_up - val), r_n);
                        add_ovf_end(acc, sgn, L.G4, x, px_up, val);
                    } else if (xk == 1) {
                        add_tail_end(acc, sgn, L.G4, x, r_n, val, px_up, px_bd, d0_n, pn_up);
                    }
                };
                int64_t d0_ab = (ak == 1) ? ((has_bup ? p4[L.p4_off(b_up, c, L.G4.INF, d)] : 0) - p4[L.p4_off(b, c, L.G4.INF, d)]) : 0;
                int64_t d0_dc = (dk == 1) ? ((has_cup ? p4[L.p4_off(b, c_up, a, L.G4.INF)] : 0) - p4[L.p4_off(b, c, a, L.G4.INF)]) : 0;
                do_end4(a, ak, Da, ca, La, a_bud, cov_a, r_a, pa_up, pa_bd, rb, cb, rc, cc, d, pb_up, d0_ab);
                do_end4(d, dk, Dd, cd, Ld, d_bud, cov_d, r_d, pd_up, pd_bd, rc, cc, rb, cb, a, pc_up, d0_dc);

                int64_t rcap = q5_to_i64(q5_ceil_div(acc, K_SUB), "rcap");
                if (rcap > loc_cap) loc_cap = rcap;
                ++rows_p4;
            }
        }
        cap_p4_ord[pid] = cap_p4_ord[pid_rev] = loc_cap;
    }

    // === 7. Verify ALL rows of S3 ===
    vector<int64_t> cap_s3_o(L.MS3.N, INT64_MIN);
    size_t rows_s3 = 0;
    #pragma omp parallel for schedule(dynamic) num_threads(eff_nt) reduction(+:rows_s3) if(eff_nt > 1)
    for (int o = 0; o < L.MS3.N; ++o) {
        int ro = L.MS3.D[o], co = L.MS3.C[o], Lo = ro - co;
        int o_up = L.MS3.up(o), o_bud = L.MS3.bud(o);
        int64_t r_o = r1[L.G2.idx[ro][co]];
        bool cov_o = (ro < h.dS4 && co >= 4 && co <= h.KS4);
        int64_t loc_cap = INT64_MIN;

        for (int tid = 0; tid < L.N_G4_tri; ++tid) {
            auto [a, b, c] = L.tri_G4[tid];
            int e[3] = {a, b, c};
            int64_t val = s3[(size_t)o * L.N_G4_tri + tid];
            int64_t ho_up = (o_up >= 0) ? s3[(size_t)o_up * L.N_G4_tri + tid] : 0;
            int64_t ho_bd = (o_bud >= 0) ? s3[(size_t)o_bud * L.N_G4_tri + tid] : 0;

            q5_i128 acc = 0;
            add_rat(acc, -(q5_i128)sgn * h.P * val, h.Q);
            for (int k = 0; k < 3; ++k) {
                int d1 = e[(k+1)%3], d2 = e[(k+2)%3];
                size_t so = (size_t)o * g4_p2 + Q5Layout::sym_off(d1, d2);
                add_src_env(acc, sgn, L.G4, e[k], r_o, h.S_r, e34s_m[so], e34s_r[so]);
            }
            int64_t do_ = ho_up - val;
            if (Lo > 0) add_rat(acc, (q5_i128)sgn * Lo * do_, 1);
            if (cov_o) add_rat(acc, (q5_i128)sgn * (co - 3) * r_o * do_, h.S_r);
            else if (co > 3) add_ext_rate(acc, sgn, do_, 0, 1, co - 3, 2);
            if (Lo > 0) {
                add_rat(acc, (q5_i128)sgn * Lo * (ho_bd - val), ro);
                if (ro < h.dS4 && co + 1 >= 4 && co + 1 <= h.KS4) {
                    int o4 = L.MS4.idx[ro][co + 1];
                    int ck[3], nk_[3], dk_[3];
                    for (int k=0;k<3;++k) L.G5.proj_int(L.G4, e[k], ck[k], nk_[k], dk_[k]);
                    int u0[2]={L.G5.INF,ck[0]}, u1[2]={L.G5.INF,ck[1]}, u2[2]={L.G5.INF,ck[2]};
                    q5_i128 w0[2]={dk_[0]-nk_[0],nk_[0]}, w1[2]={dk_[1]-nk_[1],nk_[1]}, w2[2]={dk_[2]-nk_[2],nk_[2]};
                    q5_i128 nb = 0;
                    for (int p=0;p<2;++p) if (w0[p]) for (int q=0;q<2;++q) if (w1[q]) for (int r=0;r<2;++r) if (w2[r])
                        nb += w0[p] * w1[q] * w2[r] * s4[L.s4_off(o4, u0[p], u1[q], u2[r], J_e5)];
                    add_rat(acc, (q5_i128)sgn * Lo * nb, (q5_i128)ro * dk_[0] * dk_[1] * dk_[2]);
                }
            }
            for (int k = 0; k < 3; ++k) {
                int u = e[k], u1 = e[(k+1)%3], u2 = e[(k+2)%3];
                int uk = L.G4.kind[u], Du = L.G4.D[u], cu = L.G4.C[u], Lu = L.G4.L[u];
                int u_up = L.G4.up(u), u_bud = L.G4.bud(u);
                int64_t hu_up = s3[L.s3_off(o, u_up, u1, u2)];
                int64_t hu_bd = s3[L.s3_off(o, u_bud, u1, u2)];
                if (uk == 0) {
                    add_rat(acc, (q5_i128)sgn * do_, Du);
                    int64_t du = hu_up - val;
                    if (Lu > 0) add_rat(acc, (q5_i128)sgn * Lu * du, 1);
                    add_rat(acc, (q5_i128)sgn * du, ro);
                    bool cov_u = (ro < h.d5o && co <= h.K5o && Du < h.d5c && cu >= 2 && cu <= h.K5c && ro + Du <= h.SB5);
                    if (cov_u) add_rat(acc, (q5_i128)sgn * (cu - 1) * r1[L.G2.idx[Du][cu]] * du, h.S_r);
                    else if (cu > 1) add_ext_rate(acc, sgn, du, 0, 1, cu - 1, 2);
                    if (Lu > 0) {
                        add_rat(acc, (q5_i128)sgn * Lu * (hu_bd - val), Du);
                        if (L.G4.kind[u_bud] == 0 && ro < h.d5o && co <= h.K5o && Du < h.d5c && cu + 1 >= 2 && cu + 1 <= h.K5c && ro + Du <= h.SB5) {
                            int bid = L.b5_int_map[(size_t)L.I5_3.idx[ro][co] * L.I5_2c.N + L.I5_2c.idx[Du][cu + 1]];
                            int c1, n1, d1, c2, n2, d2;
                            L.G5.proj_int(L.G4, u1, c1, n1, d1); L.G5.proj_int(L.G4, u2, c2, n2, d2);
                            int z1[2]={L.G5.INF,c1}, z2[2]={L.G5.INF,c2}; q5_i128 w1[2]={d1-n1,n1}, w2[2]={d2-n2,n2};
                            q5_i128 nb = 0;
                            for (int p=0;p<2;++p) if (w1[p]) for (int q=0;q<2;++q) if (w2[q])
                                nb += w1[p] * w2[q] * b5[L.b5_off(bid, z1[p], z2[q], J_e5)];
                            add_rat(acc, (q5_i128)sgn * Lu * nb, (q5_i128)Du * d1 * d2);
                        }
                    }
                } else if (uk == 3) {
                    add_rat(acc, (q5_i128)sgn * do_, Du);
                    add_rat(acc, (q5_i128)sgn * (hu_up - val), ro);
                    add_ovf_end(acc, sgn, L.G4, u, hu_up, val);
                } else if (uk == 1) {
                    int64_t d0 = ((o_up >= 0) ? s3[L.s3_off(o_up, L.G4.INF, u1, u2)] : 0) - s3[L.s3_off(o, L.G4.INF, u1, u2)];
                    add_tail_end(acc, sgn, L.G4, u, ro, val, hu_up, hu_bd, d0, ho_up);
                }
            }
            int64_t rcap = q5_to_i64(q5_ceil_div(acc, K_SUB), "rcap");
            if (rcap > loc_cap) loc_cap = rcap;
            ++rows_s3;
        }
        cap_s3_o[o] = loc_cap;
    }

    // === 8. Verify ALL rows of P5 (canonical b <= d; if b == d, a <= e) ===
    vector<int64_t> cap_p5_ord(L.NP5_ord, INT64_MIN);
    size_t rows_p5 = 0;
    int N52 = L.I5_2.N;
    #pragma omp parallel for schedule(dynamic) num_threads(eff_nt) reduction(+:rows_p5) if(eff_nt > 1)
    for (int pid = 0; pid < L.NP5_ord; ++pid) {
        auto [b, c, d] = L.p5_ord[pid];
        if (b > d) continue;
        int rb = L.I5_2.D[b], cb = L.I5_2.C[b], Lb = rb - cb;
        int rc = L.I5_2.D[c], cc = L.I5_2.C[c], Lc = rc - cc;
        int rd = L.I5_2.D[d], cd = L.I5_2.C[d], Ld = rd - cd;
        int b_up = L.I5_2.up(b), c_up = L.I5_2.up(c), d_up = L.I5_2.up(d);
        int b_bd = L.I5_2.bud(b), c_bd = L.I5_2.bud(c), d_bd = L.I5_2.bud(d);
        int pid_rev = L.p5_ord_map[((size_t)d * N52 + c) * N52 + b];
        bool has_bup = (b_up >= 0 && rb+1+rc <= h.SP4 && rb+1+rc+rd <= h.SP5);
        bool has_cup = (c_up >= 0 && rb+rc+1 <= h.SP4 && rc+1+rd <= h.SP4 && rb+rc+1+rd <= h.SP5);
        bool has_dup = (d_up >= 0 && rc+rd+1 <= h.SP4 && rb+rc+rd+1 <= h.SP5);
        int64_t r_b = r1[L.G2.idx[rb][cb]], r_d = r1[L.G2.idx[rd][cd]];
        int64_t loc_cap = INT64_MIN;

        for (int a = 0; a < N5; ++a) {
            int ak = L.G5.kind[a], Da = L.G5.D[a], ca = L.G5.C[a], La = L.G5.L[a];
            bool a_in = L.G5.is_inner(a);
            int a_up = L.G5.up(a), a_bud = L.G5.bud(a);
            int e_start = (b == d) ? a : 0;
            for (int e = e_start; e < N5; ++e) {
                int ek = L.G5.kind[e], De = L.G5.D[e], ce = L.G5.C[e], Le = L.G5.L[e];
                bool e_in = L.G5.is_inner(e);
                int e_up = L.G5.up(e), e_bud = L.G5.bud(e);

                int64_t val = p5[L.p5_off(b, c, d, a, e)];
                int64_t pa_up = p5[L.p5_off(b, c, d, a_up, e)], pe_up = p5[L.p5_off(b, c, d, a, e_up)];
                int64_t pb_up = has_bup ? p5[L.p5_off(b_up, c, d, a, e)] : 0;
                int64_t pc_up = has_cup ? p5[L.p5_off(b, c_up, d, a, e)] : 0;
                int64_t pd_up = has_dup ? p5[L.p5_off(b, c, d_up, a, e)] : 0;
                int64_t pa_bd = p5[L.p5_off(b, c, d, a_bud, e)], pe_bd = p5[L.p5_off(b, c, d, a, e_bud)];
                int64_t pb_bd = (b_bd >= 0) ? p5[L.p5_off(b_bd, c, d, a, e)] : 0;
                int64_t pc_bd = (c_bd >= 0) ? p5[L.p5_off(b, c_bd, d, a, e)] : 0;
                int64_t pd_bd = (d_bd >= 0) ? p5[L.p5_off(b, c, d_bd, a, e)] : 0;

                q5_i128 acc = 0;
                add_rat(acc, -(q5_i128)sgn * h.P * val, h.Q);
                add_src_env(acc, sgn, L.G5, a, r_b, h.S_r, e45p_m[(size_t)pid * N5 + e], e45p_r[(size_t)pid * N5 + e]);
                add_src_env(acc, sgn, L.G5, e, r_d, h.S_r, e45p_m[(size_t)pid_rev * N5 + a], e45p_r[(size_t)pid_rev * N5 + a]);

                // Internal b, c, d
                int64_t db = pb_up - val, dc_ = pc_up - val, dd = pd_up - val;
                if (Lb > 0) { add_rat(acc, (q5_i128)sgn * Lb * db, 1); add_rat(acc, (q5_i128)sgn * Lb * (pb_bd - val), rb); }
                if (a_in) add_rat(acc, (q5_i128)sgn * db, Da);
                add_rat(acc, (q5_i128)sgn * db, rc);
                if (cb > 2) add_ext_rate(acc, sgn, db, 0, 1, cb - 2, 2);

                if (Lc > 0) { add_rat(acc, (q5_i128)sgn * Lc * dc_, 1); add_rat(acc, (q5_i128)sgn * Lc * (pc_bd - val), rc); }
                add_rat(acc, (q5_i128)sgn * dc_, rb); add_rat(acc, (q5_i128)sgn * dc_, rd);
                if (cc > 2) add_ext_rate(acc, sgn, dc_, 0, 1, cc - 2, 2);

                if (Ld > 0) { add_rat(acc, (q5_i128)sgn * Ld * dd, 1); add_rat(acc, (q5_i128)sgn * Ld * (pd_bd - val), rd); }
                add_rat(acc, (q5_i128)sgn * dd, rc);
                if (e_in) add_rat(acc, (q5_i128)sgn * dd, De);
                if (cd > 2) add_ext_rate(acc, sgn, dd, 0, 1, cd - 2, 2);

                // Endpoints a and e
                auto do_end5 = [&](int x, int xk, int Dx, int cx, int Lx, int64_t px_up, int64_t px_bd, int r_n, int64_t pn_up, int64_t d0_n) {
                    if (xk == 0) {
                        int64_t dx = px_up - val;
                        if (Lx > 0) { add_rat(acc, (q5_i128)sgn * Lx * dx, 1); add_rat(acc, (q5_i128)sgn * Lx * (px_bd - val), Dx); }
                        add_rat(acc, (q5_i128)sgn * dx, r_n);
                        if (cx > 1) add_ext_rate(acc, sgn, dx, 0, 1, cx - 1, 2);
                    } else if (xk == 3) {
                        add_rat(acc, (q5_i128)sgn * (px_up - val), r_n);
                        add_ovf_end(acc, sgn, L.G5, x, px_up, val);
                    } else if (xk == 1) {
                        add_tail_end(acc, sgn, L.G5, x, r_n, val, px_up, px_bd, d0_n, pn_up);
                    }
                };
                int64_t d0_ab = (ak == 1) ? ((has_bup ? p5[L.p5_off(b_up, c, d, L.G5.INF, e)] : 0) - p5[L.p5_off(b, c, d, L.G5.INF, e)]) : 0;
                int64_t d0_ed = (ek == 1) ? ((has_dup ? p5[L.p5_off(b, c, d_up, a, L.G5.INF)] : 0) - p5[L.p5_off(b, c, d, a, L.G5.INF)]) : 0;
                do_end5(a, ak, Da, ca, La, pa_up, pa_bd, rb, pb_up, d0_ab);
                do_end5(e, ek, De, ce, Le, pe_up, pe_bd, rd, pd_up, d0_ed);

                int64_t rcap = q5_to_i64(q5_ceil_div(acc, K_SUB), "rcap");
                if (rcap > loc_cap) loc_cap = rcap;
                ++rows_p5;
            }
        }
        cap_p5_ord[pid] = cap_p5_ord[pid_rev] = loc_cap;
    }

    // === 9. Verify ALL rows of B5 ===
    vector<int64_t> cap_b5_id(L.NB5_int, INT64_MIN);
    size_t rows_b5 = 0;
    #pragma omp parallel for schedule(dynamic) num_threads(eff_nt) reduction(+:rows_b5) if(eff_nt > 1)
    for (int bid = 0; bid < L.NB5_int; ++bid) {
        auto [o, c] = L.b5_int[bid];
        int ro = L.I5_3.D[o], co = L.I5_3.C[o], Lo = ro - co;
        int rc = L.I5_2c.D[c], cc = L.I5_2c.C[c], Lc = rc - cc;
        int o_up = L.I5_3.up(o), c_up = L.I5_2c.up(c), o_bd = L.I5_3.bud(o), c_bd = L.I5_2c.bud(c);
        int b_oup = (o_up >= 0) ? L.b5_int_map[(size_t)o_up * L.I5_2c.N + c] : -1;
        int b_cup = (c_up >= 0) ? L.b5_int_map[(size_t)o * L.I5_2c.N + c_up] : -1;
        int b_obd = (o_bd >= 0) ? L.b5_int_map[(size_t)o_bd * L.I5_2c.N + c] : -1;
        int b_cbd = (c_bd >= 0) ? L.b5_int_map[(size_t)o * L.I5_2c.N + c_bd] : -1;
        int64_t r_o = r1[L.G2.idx[ro][co]], r_c = r1[L.G2.idx[rc][cc]];
        int64_t loc_cap = INT64_MIN;

        for (int a = 0; a < N5; ++a) {
            int ak = L.G5.kind[a], Da = L.G5.D[a], ca = L.G5.C[a], La = L.G5.L[a];
            bool a_in = L.G5.is_inner(a);
            int a_up = L.G5.up(a), a_bud = L.G5.bud(a);
            for (int b = a; b < N5; ++b) {
                int bk = L.G5.kind[b], Db = L.G5.D[b], cb = L.G5.C[b], Lb = L.G5.L[b];
                bool b_in = L.G5.is_inner(b);
                int b_up = L.G5.up(b), b_bud = L.G5.bud(b);
                size_t obs = (size_t)bid * g5_p2 + Q5Layout::sym_off(a, b);
                for (int d = 0; d < N5; ++d) {
                    int dk = L.G5.kind[d], Dd = L.G5.D[d], cd = L.G5.C[d], Ld = L.G5.L[d];
                    bool d_in = L.G5.is_inner(d);
                    int d_up = L.G5.up(d), d_bud = L.G5.bud(d);

                    int64_t val = b5[L.b5_off(bid, a, b, d)];
                    int64_t bo_up = (b_oup >= 0) ? b5[L.b5_off(b_oup, a, b, d)] : 0;
                    int64_t bc_up = (b_cup >= 0) ? b5[L.b5_off(b_cup, a, b, d)] : 0;
                    int64_t ba_up = b5[L.b5_off(bid, a_up, b, d)];
                    int64_t bb_up = b5[L.b5_off(bid, a, b_up, d)];
                    int64_t bd_up = b5[L.b5_off(bid, a, b, d_up)];
                    int64_t bo_bd = (b_obd >= 0) ? b5[L.b5_off(b_obd, a, b, d)] : 0;
                    int64_t bc_bd = (b_cbd >= 0) ? b5[L.b5_off(b_cbd, a, b, d)] : 0;
                    int64_t ba_bd = b5[L.b5_off(bid, a_bud, b, d)];
                    int64_t bb_bd = b5[L.b5_off(bid, a, b_bud, d)];
                    int64_t bd_bd = b5[L.b5_off(bid, a, b, d_bud)];

                    q5_i128 acc = 0;
                    add_rat(acc, -(q5_i128)sgn * h.P * val, h.Q);
                    add_src_env(acc, sgn, L.G5, b, r_o, h.S_r, e45bp_m[((size_t)bid*N5 + a)*N5 + d], e45bp_r[((size_t)bid*N5 + a)*N5 + d]);
                    add_src_env(acc, sgn, L.G5, a, r_o, h.S_r, e45bp_m[((size_t)bid*N5 + b)*N5 + d], e45bp_r[((size_t)bid*N5 + b)*N5 + d]);
                    add_src_env(acc, sgn, L.G5, d, r_c, h.S_r, e45bs_m[obs], e45bs_r[obs]);

                    int64_t do_ = bo_up - val, dc_ = bc_up - val;
                    if (Lo > 0) { add_rat(acc, (q5_i128)sgn * Lo * do_, 1); add_rat(acc, (q5_i128)sgn * Lo * (bo_bd - val), ro); }
                    if (a_in) add_rat(acc, (q5_i128)sgn * do_, Da);
                    if (b_in) add_rat(acc, (q5_i128)sgn * do_, Db);
                    add_rat(acc, (q5_i128)sgn * do_, rc);
                    if (co > 3) add_ext_rate(acc, sgn, do_, 0, 1, co - 3, 2);

                    if (Lc > 0) { add_rat(acc, (q5_i128)sgn * Lc * dc_, 1); add_rat(acc, (q5_i128)sgn * Lc * (bc_bd - val), rc); }
                    add_rat(acc, (q5_i128)sgn * dc_, ro);
                    if (d_in) add_rat(acc, (q5_i128)sgn * dc_, Dd);
                    if (cc > 2) add_ext_rate(acc, sgn, dc_, 0, 1, cc - 2, 2);

                    auto do_e = [&](int x, int xk, int Dx, int cx, int Lx, int64_t vx_up, int64_t vx_bd, int r_n, int64_t vn_up, int64_t d0_n) {
                        if (xk == 0) {
                            int64_t dx = vx_up - val;
                            if (Lx > 0) { add_rat(acc, (q5_i128)sgn * Lx * dx, 1); add_rat(acc, (q5_i128)sgn * Lx * (vx_bd - val), Dx); }
                            add_rat(acc, (q5_i128)sgn * dx, r_n);
                            if (cx > 1) add_ext_rate(acc, sgn, dx, 0, 1, cx - 1, 2);
                        } else if (xk == 3) {
                            add_rat(acc, (q5_i128)sgn * (vx_up - val), r_n);
                            add_ovf_end(acc, sgn, L.G5, x, vx_up, val);
                        } else if (xk == 1) {
                            add_tail_end(acc, sgn, L.G5, x, r_n, val, vx_up, vx_bd, d0_n, vn_up);
                        }
                    };
                    int64_t d0_ao = (ak == 1) ? (((b_oup >= 0) ? b5[L.b5_off(b_oup, L.G5.INF, b, d)] : 0) - b5[L.b5_off(bid, L.G5.INF, b, d)]) : 0;
                    int64_t d0_bo = (bk == 1) ? (((b_oup >= 0) ? b5[L.b5_off(b_oup, a, L.G5.INF, d)] : 0) - b5[L.b5_off(bid, a, L.G5.INF, d)]) : 0;
                    int64_t d0_dc = (dk == 1) ? (((b_cup >= 0) ? b5[L.b5_off(b_cup, a, b, L.G5.INF)] : 0) - b5[L.b5_off(bid, a, b, L.G5.INF)]) : 0;
                    do_e(a, ak, Da, ca, La, ba_up, ba_bd, ro, bo_up, d0_ao);
                    do_e(b, bk, Db, cb, Lb, bb_up, bb_bd, ro, bo_up, d0_bo);
                    do_e(d, dk, Dd, cd, Ld, bd_up, bd_bd, rc, bc_up, d0_dc);

                    int64_t rcap = q5_to_i64(q5_ceil_div(acc, K_SUB), "rcap");
                    if (rcap > loc_cap) loc_cap = rcap;
                    ++rows_b5;
                }
            }
        }
        cap_b5_id[bid] = loc_cap;
    }

    // === 10. Verify ALL rows of S4 ===
    vector<int64_t> cap_s4_o(L.MS4.N, INT64_MIN);
    size_t rows_s4 = 0;
    #pragma omp parallel for schedule(dynamic) num_threads(eff_nt) reduction(+:rows_s4) if(eff_nt > 1)
    for (int o = 0; o < L.MS4.N; ++o) {
        int ro = L.MS4.D[o], co = L.MS4.C[o], Lo = ro - co;
        int o_up = L.MS4.up(o), o_bd = L.MS4.bud(o);
        int64_t r_o = r1[L.G2.idx[ro][co]];
        int64_t loc_cap = INT64_MIN;

        for (int qid = 0; qid < L.N_G5_quad; ++qid) {
            auto [a, b, c, d] = L.quad_G5[qid];
            int e[4] = {a, b, c, d};
            int64_t val = s4[(size_t)o * L.N_G5_quad + qid];
            int64_t ho_up = (o_up >= 0) ? s4[(size_t)o_up * L.N_G5_quad + qid] : 0;
            int64_t ho_bd = (o_bd >= 0) ? s4[(size_t)o_bd * L.N_G5_quad + qid] : 0;

            q5_i128 acc = 0;
            add_rat(acc, -(q5_i128)sgn * h.P * val, h.Q);
            for (int k = 0; k < 4; ++k) {
                int u1 = e[(k+1)%4], u2 = e[(k+2)%4], u3 = e[(k+3)%4];
                int tid5 = L.tri_G5_map[((size_t)u1 * N5 + u2) * N5 + u3];
                add_src_env(acc, sgn, L.G5, e[k], r_o, h.S_r, e45s4_m[(size_t)o * L.N_G5_tri + tid5], e45s4_r[(size_t)o * L.N_G5_tri + tid5]);
            }
            int64_t do_ = ho_up - val;
            if (Lo > 0) { add_rat(acc, (q5_i128)sgn * Lo * do_, 1); add_rat(acc, (q5_i128)sgn * Lo * (ho_bd - val), ro); }
            if (co > 4) add_ext_rate(acc, sgn, do_, 0, 1, co - 4, 2);
            for (int k = 0; k < 4; ++k) {
                int u = e[k], u1 = e[(k+1)%4], u2 = e[(k+2)%4], u3 = e[(k+3)%4];
                int uk = L.G5.kind[u], Du = L.G5.D[u], cu = L.G5.C[u], Lu = L.G5.L[u];
                int u_up = L.G5.up(u), u_bud = L.G5.bud(u);
                int64_t hu_up = s4[L.s4_off(o, u_up, u1, u2, u3)];
                int64_t hu_bd = s4[L.s4_off(o, u_bud, u1, u2, u3)];
                if (uk == 0) {
                    add_rat(acc, (q5_i128)sgn * do_, Du);
                    int64_t du = hu_up - val;
                    if (Lu > 0) { add_rat(acc, (q5_i128)sgn * Lu * du, 1); add_rat(acc, (q5_i128)sgn * Lu * (hu_bd - val), Du); }
                    add_rat(acc, (q5_i128)sgn * du, ro);
                    if (cu > 1) add_ext_rate(acc, sgn, du, 0, 1, cu - 1, 2);
                } else if (uk == 3) {
                    add_rat(acc, (q5_i128)sgn * do_, Du);
                    add_rat(acc, (q5_i128)sgn * (hu_up - val), ro);
                    add_ovf_end(acc, sgn, L.G5, u, hu_up, val);
                } else if (uk == 1) {
                    int64_t d0 = ((o_up >= 0) ? s4[L.s4_off(o_up, L.G5.INF, u1, u2, u3)] : 0) - s4[L.s4_off(o, L.G5.INF, u1, u2, u3)];
                    add_tail_end(acc, sgn, L.G5, u, ro, val, hu_up, hu_bd, d0, ho_up);
                }
            }
            int64_t rcap = q5_to_i64(q5_ceil_div(acc, K_SUB), "rcap");
            if (rcap > loc_cap) loc_cap = rcap;
            ++rows_s4;
        }
        cap_s4_o[o] = loc_cap;
    }

    // === 11. Aggregate positive higher-family caps per core vertex (strict upper bound!) ===
    auto max_vec = [](const vector<int64_t>& v) { int64_t r = INT64_MIN; for (int64_t x : v) if (x > r) r = x; return r; };
    int64_t max_h3 = max_vec(cap_h3_m), max_p4 = max_vec(cap_p4_ord), max_s3 = max_vec(cap_s3_o);
    int64_t max_p5 = max_vec(cap_p5_ord), max_b5 = max_vec(cap_b5_id), max_s4 = max_vec(cap_s4_o);

    // Precompute per-vertex-state positive higher cap for all 2 <= D < d3m, 1 <= c <= min(D, K3m)!
    int64_t higher_per_core = 0;
    // Also precompute per-vertex-state value extrema for global comparability!
    vector<int64_t> g_min_u(NG, INT64_MAX), g_max_u(NG, INT64_MIN);
    for (int u = 0; u < NG; ++u)
        for (int p = 0; p < NG; ++p) {
            int64_t v = g[L.g_off(p, u)];
            if (v < g_min_u[u]) g_min_u[u] = v;
            if (v > g_max_u[u]) g_max_u[u] = v;
        }
    vector<int64_t> h3_min_m(L.M3.N, INT64_MAX), h3_max_m(L.M3.N, INT64_MIN);
    for (int m = 0; m < L.M3.N; ++m)
        for (int a = 0; a < N3; ++a) for (int b = a; b < N3; ++b) {
            int64_t v = h3[L.h3_off(m, a, b)];
            if (v < h3_min_m[m]) h3_min_m[m] = v;
            if (v > h3_max_m[m]) h3_max_m[m] = v;
        }
    vector<int64_t> p4_min_o(L.NP4_ord, INT64_MAX), p4_max_o(L.NP4_ord, INT64_MIN);
    for (int pid = 0; pid < L.NP4_ord; ++pid) {
        auto [b, c] = L.p4_ord[pid];
        int64_t mn = INT64_MAX, mx = INT64_MIN;
        for (int a = 0; a < N4; ++a) for (int d = 0; d < N4; ++d) {
            int64_t v = p4[L.p4_off(b, c, a, d)];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        p4_min_o[pid] = mn; p4_max_o[pid] = mx;
    }
    vector<int64_t> s3_min_o(L.MS3.N, INT64_MAX), s3_max_o(L.MS3.N, INT64_MIN);
    for (int o = 0; o < L.MS3.N; ++o)
        for (int tid = 0; tid < L.N_G4_tri; ++tid) {
            int64_t v = s3[(size_t)o * L.N_G4_tri + tid];
            if (v < s3_min_o[o]) s3_min_o[o] = v;
            if (v > s3_max_o[o]) s3_max_o[o] = v;
        }
    vector<int64_t> p5_min_o(L.NP5_ord, INT64_MAX), p5_max_o(L.NP5_ord, INT64_MIN);
    for (int pid = 0; pid < L.NP5_ord; ++pid) {
        auto [b, c, d] = L.p5_ord[pid];
        int64_t mn = INT64_MAX, mx = INT64_MIN;
        for (int a = 0; a < N5; ++a) for (int e = 0; e < N5; ++e) {
            int64_t v = p5[L.p5_off(b, c, d, a, e)];
            if (v < mn) mn = v;
            if (v > mx) mx = v;
        }
        p5_min_o[pid] = mn; p5_max_o[pid] = mx;
    }
    vector<int64_t> b5_min_id(L.NB5_int, INT64_MAX), b5_max_id(L.NB5_int, INT64_MIN);
    for (int bid = 0; bid < L.NB5_int; ++bid)
        for (int a = 0; a < N5; ++a) for (int b = a; b < N5; ++b) for (int d = 0; d < N5; ++d) {
            int64_t v = b5[L.b5_off(bid, a, b, d)];
            if (v < b5_min_id[bid]) b5_min_id[bid] = v;
            if (v > b5_max_id[bid]) b5_max_id[bid] = v;
        }
    vector<int64_t> s4_min_o(L.MS4.N, INT64_MAX), s4_max_o(L.MS4.N, INT64_MIN);
    for (int o = 0; o < L.MS4.N; ++o)
        for (int qid = 0; qid < L.N_G5_quad; ++qid) {
            int64_t v = s4[(size_t)o * L.N_G5_quad + qid];
            if (v < s4_min_o[o]) s4_min_o[o] = v;
            if (v > s4_max_o[o]) s4_max_o[o] = v;
        }

    // Scan finite vertices D < DF for both higher_per_core AND comparability [min_W_fin, max_W_fin]!
    // Note: alpha shift (alpha*(c - 2)) is in ticks: alpha_ticks = (alpha_lo_num * S) / S_r (exact since S_r = 2^30 divides S = 2^40!)
    assert(h.S % h.S_r == 0);
    int64_t alpha_ticks = h.alpha_lo_num * (h.S / h.S_r);
    int64_t min_W_c0_fin = INT64_MAX, max_W_c0_fin = INT64_MIN;
    int64_t min_W_fin = INT64_MAX, max_W_fin = INT64_MIN;

    for (int D = 2; D < h.DF; ++D) {
        // c == 0 (single core tree C = 1)
        int64_t fc0 = f[L.f_id(D, 0)];
        if (fc0 < min_W_c0_fin) min_W_c0_fin = fc0;
        if (fc0 > max_W_c0_fin) max_W_c0_fin = fc0;

        for (int c = 1; c <= D; ++c) {
            q5_i128 hcap = 0;
            q5_i128 lo_w = f[L.f_id(D, c)], hi_w = f[L.f_id(D, c)];

            // Edge g(p_u, u) allocated to child u != r (0 if u == r, so include 0!)
            int64_t g_lo_v, g_hi_v;
            if (D < h.d2) {
                int id = (c <= h.K2) ? L.G2.idx[D][c] : L.G2.ovf_idx[D];
                g_lo_v = min<int64_t>(0, g_min_u[id]);
                g_hi_v = max<int64_t>(0, g_max_u[id]);
            } else {
                int tid = L.G2.t_id(c);
                q5_i128 mn_i = ((q5_i128)(D - h.d2) * g_min_u[L.G2.INF] + (q5_i128)h.d2 * g_min_u[tid]);
                q5_i128 mx_i = ((q5_i128)(D - h.d2) * g_max_u[L.G2.INF] + (q5_i128)h.d2 * g_max_u[tid]);
                g_lo_v = min<int64_t>(0, q5_to_i64(q5_floor_div(mn_i, D), "g_lo_v"));
                g_hi_v = max<int64_t>(0, q5_to_i64(q5_ceil_div(mx_i, D), "g_hi_v"));
            }
            lo_w += g_lo_v; hi_w += g_hi_v;

            // H3 at middle m = (D, c)
            if (D < h.d3m && c >= 2 && c <= h.K3m) {
                int m = L.M3.idx[D][c];
                int64_t b2 = binom2(c);
                if (cap_h3_m[m] > 0) hcap += (q5_i128)b2 * cap_h3_m[m];
                lo_w += (q5_i128)b2 * min<int64_t>(0, h3_min_m[m]);
                hi_w += (q5_i128)b2 * max<int64_t>(0, h3_max_m[m]);
            }
            // P4: 1/2 allocated to internal vertex b = (D, c)
            // For each neighbor y of b (c choices of y, c-1 choices of a, c_y-1 choices of d):
            // 1/2 * sum_{y~b} (c-1)*(c_y-1)*val <= 1/2 * c*(c-1) * max_y (c_y-1)*val (including 0)!
            if (D < h.d4i && c >= 2 && c <= h.K4i) {
                int b = L.M4v.idx[D][c];
                int64_t mult = (int64_t)c * (c - 1); // divide by 2 via floor/ceil!
                q5_i128 c_max = 0, v_mn = 0, v_mx = 0;
                for (int y = 0; y < L.M4v.N; ++y) {
                    int pid = L.p4_ord_map[(size_t)b * L.M4v.N + y];
                    if (pid < 0) continue;
                    int cy = L.M4v.C[y];
                    if (cap_p4_ord[pid] > 0) {
                        q5_i128 t = (q5_i128)(cy - 1) * cap_p4_ord[pid];
                        if (t > c_max) c_max = t;
                    }
                    q5_i128 tmn = (q5_i128)(cy - 1) * p4_min_o[pid];
                    q5_i128 tmx = (q5_i128)(cy - 1) * p4_max_o[pid];
                    if (tmn < v_mn) v_mn = tmn;
                    if (tmx > v_mx) v_mx = tmx;
                }
                hcap += q5_ceil_div(mult * c_max, 2);
                lo_w += q5_floor_div(mult * v_mn, 2);
                hi_w += q5_ceil_div(mult * v_mx, 2);
            }
            // S3 at hub o = (D, c)
            if (D < h.dS3 && c >= 3 && c <= h.KS3) {
                int o = L.MS3.idx[D][c];
                int64_t b3 = binom3(c);
                if (cap_s3_o[o] > 0) hcap += (q5_i128)b3 * cap_s3_o[o];
                lo_w += (q5_i128)b3 * min<int64_t>(0, s3_min_o[o]);
                hi_w += (q5_i128)b3 * max<int64_t>(0, s3_max_o[o]);
            }
            // P5 at center c_mid = (D, c): binom2(c) unordered pairs {b, d} of neighbors of c_mid
            if (D < h.d5i && c >= 2 && c <= h.K5i) {
                int cmid = L.I5_2.idx[D][c];
                int64_t b2 = binom2(c);
                q5_i128 c_max = 0, v_mn = 0, v_mx = 0;
                for (int b = 0; b < N52; ++b) for (int d = b; d < N52; ++d) {
                    int pid = L.p5_ord_map[((size_t)b * N52 + cmid) * N52 + d];
                    if (pid < 0) continue;
                    int64_t ed = (int64_t)(L.I5_2.C[b] - 1) * (L.I5_2.C[d] - 1);
                    if (cap_p5_ord[pid] > 0) { q5_i128 t = (q5_i128)ed * cap_p5_ord[pid]; if (t > c_max) c_max = t; }
                    q5_i128 tmn = (q5_i128)ed * p5_min_o[pid], tmx = (q5_i128)ed * p5_max_o[pid];
                    if (tmn < v_mn) v_mn = tmn;
                    if (tmx > v_mx) v_mx = tmx;
                }
                hcap += (q5_i128)b2 * c_max;
                lo_w += (q5_i128)b2 * v_mn;
                hi_w += (q5_i128)b2 * v_mx;
            }
            // B5 at hub o = (D, c): c choices of arm y, binom2(c-1) pairs {a,b}, (c_y-1) choices of d
            if (D < h.d5o && c >= 3 && c <= h.K5o) {
                int o = L.I5_3.idx[D][c];
                int64_t mult = (int64_t)c * binom2(c - 1);
                q5_i128 c_max = 0, v_mn = 0, v_mx = 0;
                for (int y = 0; y < L.I5_2c.N; ++y) {
                    int bid = L.b5_int_map[(size_t)o * L.I5_2c.N + y];
                    if (bid < 0) continue;
                    int cy = L.I5_2c.C[y];
                    if (cap_b5_id[bid] > 0) { q5_i128 t = (q5_i128)(cy - 1) * cap_b5_id[bid]; if (t > c_max) c_max = t; }
                    q5_i128 tmn = (q5_i128)(cy - 1) * b5_min_id[bid], tmx = (q5_i128)(cy - 1) * b5_max_id[bid];
                    if (tmn < v_mn) v_mn = tmn;
                    if (tmx > v_mx) v_mx = tmx;
                }
                hcap += (q5_i128)mult * c_max;
                lo_w += (q5_i128)mult * v_mn;
                hi_w += (q5_i128)mult * v_mx;
            }
            // S4 at hub o = (D, c)
            if (D < h.dS4 && c >= 4 && c <= h.KS4) {
                int o = L.MS4.idx[D][c];
                int64_t b4 = binom4(c);
                if (cap_s4_o[o] > 0) hcap += (q5_i128)b4 * cap_s4_o[o];
                lo_w += (q5_i128)b4 * min<int64_t>(0, s4_min_o[o]);
                hi_w += (q5_i128)b4 * max<int64_t>(0, s4_max_o[o]);
            }

            int64_t hcap_i64 = q5_to_i64(hcap, "higher_per_core");
            if (hcap_i64 > higher_per_core) higher_per_core = hcap_i64;
            int64_t lo_shifted_i64 = q5_to_i64(lo_w + (q5_i128)alpha_ticks * (c - 2), "min_W_fin");
            if (lo_shifted_i64 < min_W_fin) min_W_fin = lo_shifted_i64;
            int64_t hi_w_i64 = q5_to_i64(hi_w, "max_W_fin");
            if (hi_w_i64 > max_W_fin) max_W_fin = hi_w_i64;
        }
    }

    // Far-tail edge bounds for D >= DF (affine in t = d2/D in [0, d2/DF])!
    int64_t g_tail_mn = INT64_MAX, g_tail_mx = INT64_MIN;
    for (int kap = 1; kap <= h.K2 + 1; ++kap) {
        int tid = L.G2.t_id(kap);
        if (g_min_u[tid] < g_tail_mn) g_tail_mn = g_min_u[tid];
        if (g_max_u[tid] > g_tail_mx) g_tail_mx = g_max_u[tid];
    }
    int64_t g_inf_mn = g_min_u[L.G2.INF], g_inf_mx = g_max_u[L.G2.INF];
    int64_t g_far_min = min<int64_t>(0, min<int64_t>(g_inf_mn, q5_to_i64(q5_floor_div((q5_i128)(h.DF - h.d2) * g_inf_mn + (q5_i128)h.d2 * g_tail_mn, h.DF), "g_far_min")));
    int64_t g_far_max = max<int64_t>(0, max<int64_t>(g_inf_mx, q5_to_i64(q5_ceil_div((q5_i128)(h.DF - h.d2) * g_inf_mx + (q5_i128)h.d2 * g_tail_mx, h.DF), "g_far_max")));

    // Output strict JSON record (including SHA-256 of the exact in-memory input bytes verified)!
    printf("{\n");
    printf("  \"format\": \"FT2_U5V1\",\n");
    printf("  \"sha256\": \"%s\",\n", input_sha256.c_str());
    printf("  \"mode\": %" PRId64 ",\n", h.mode);
    printf("  \"P\": %" PRId64 ", \"Q\": %" PRId64 ",\n", h.P, h.Q);
    printf("  \"S\": %" PRId64 ", \"S_r\": %" PRId64 ", \"alpha_ticks\": %" PRId64 ",\n", h.S, h.S_r, alpha_ticks);
    printf("  \"envelope_checks\": {\"1d\": %zu, \"2d\": %zu, \"3d\": %zu},\n", c_env1, c_env2, c_env3);
    printf("  \"rows\": {\"f_in\": %zu, \"f_tab\": %zu, \"g\": %zu, \"h3\": %zu, \"p4\": %zu, \"s3\": %zu, \"p5\": %zu, \"b5\": %zu, \"s4\": %zu},\n",
           rows_f_in, rows_f_tab, rows_g, rows_h3, rows_p4, rows_s3, rows_p5, rows_b5, rows_s4);
    printf("  \"caps\": {\"f_in\": %" PRId64 ", \"f_tab\": %" PRId64 ", \"g\": %" PRId64 ", \"h3\": %" PRId64 ", \"p4\": %" PRId64 ", \"s3\": %" PRId64 ", \"p5\": %" PRId64 ", \"b5\": %" PRId64 ", \"s4\": %" PRId64 "},\n",
           cap_f_in, cap_f_tab, cap_g, max_h3, max_p4, max_s3, max_p5, max_b5, max_s4);
    printf("  \"higher_per_core\": %" PRId64 ",\n", higher_per_core);
    printf("  \"comp_fin\": {\"min_W_c0\": %" PRId64 ", \"max_W_c0\": %" PRId64 ", \"min_W\": %" PRId64 ", \"max_W\": %" PRId64 ", \"g_far_min\": %" PRId64 ", \"g_far_max\": %" PRId64 "}\n",
           min_W_c0_fin, max_W_c0_fin, min_W_fin, max_W_fin, g_far_min, g_far_max);
    printf("}\n");
    return 0;
}
