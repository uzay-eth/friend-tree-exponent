// Read-only view of the shared frozen FT2_4V01 certificate.
#pragma once
#include "../motif/cert4v_format.h"
#include <algorithm>
#include <array>
#include <cerrno>
#include <climits>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

#if __BYTE_ORDER__ != __ORDER_LITTLE_ENDIAN__
#error "FT2_4V01 uses little-endian integers"
#endif
namespace mh {
static_assert(sizeof(size_t) == 8 && sizeof(long) == 8, "64-bit LP64 build required");
using I = int64_t;
using Wide = __int128_t;
using Labels = std::array<int, 4>;
enum Family { Edge, V, C, P, F, Y, T, NFamily };
constexpr const char *names[] = {"E", "V", "C", "P", "F", "Y", "T"};

// Order is also the coordinate order accepted by the exact row-query interface.
// Unordered branches describe DISTINCT physical children, even if their labels agree.
struct Shape { int n; Labels parent; };
constexpr Shape shapes[NFamily] = {
    {2, {-1, 0}}, {3, {-1, 0, 1}}, {3, {-1, 0, 0}},
    {4, {-1, 0, 1, 2}}, {4, {-1, 0, 0, 2}}, {4, {-1, 0, 1, 1}},
    {4, {-1, 0, 0, 0}}
};

inline void require(bool ok, const std::string &why) {
    if (!ok) throw std::runtime_error(why);
}

struct Mass {
    int a, b, na, nb;  // weights na/(na+nb), nb/(na+nb); nb=0 means a point mass
};
struct State {
    int d, k, l, type, up, bud, center;
    Mass coarse;
};

struct Snapshot {
    static_assert(sizeof(Cert4vHeader) == 248, "unexpected certificate ABI");
    Cert4vHeader h{};
    Cert4vLayout c{};
    Cert4vSizes z{};
    size_t bytes{};
    void *mapping{MAP_FAILED};
    const I *cent[5], *weight[NFamily];
    const I *e23[4], *ep[4], *ef[4], *eyv[2], *eyc[2];
    int sign{};
    std::vector<State> s3, s4;
    std::vector<int> hu3, cu3, hu4, cu4, yi4, ti4;
    std::vector<int> hu_state, cu_state;

    static Mass project(const Cert4vLevel &fine, const Cert4vLevel &coarse, int i) {
        if (i == fine.INF_ID) return {coarse.INF_ID, 0, 1, 0};
        int d = fine.gD[i], k = fine.gK[i];
        if (d < coarse.D0) return {c4_lvl_inner(&coarse, fine.gL[i], k), 0, 1, 0};
        int tail = c4_lvl_tail(&coarse, k);
        if (d == coarse.D0) return {tail, 0, 1, 0};
        return {tail, coarse.INF_ID, coarse.D0, d - coarse.D0};
    }

    explicit Snapshot(const char *path) {
        const int fd = open(path, O_RDONLY);
        require(fd >= 0, std::string("open ") + path + " failed");
        struct stat st{};
        if (fstat(fd, &st) != 0 || st.st_size < (off_t)sizeof h) {
            close(fd); throw std::runtime_error("truncated header");
        }
        bytes = (size_t)st.st_size;
        if (bytes > (size_t(5) << 30)) {
            close(fd); throw std::runtime_error("certificate exceeds workspace budget");
        }
        mapping = mmap(nullptr, bytes, PROT_READ, MAP_PRIVATE, fd, 0);
        close(fd);
        require(mapping != MAP_FAILED, "mmap failed");
        std::memcpy(&h, mapping, sizeof h);
        // Bounds before invoking shared C routines that use signed int indexing.
        // These bounds keep all intermediate layout products and fixed-point operands defined.
        require(h.DF >= 4 && h.DF <= 65535 && h.KHU >= 0 && h.KHU <= 128 &&
                h.KT_2 >= 0 && h.KT_2 <= 128 && h.KT_3 >= 0 && h.KT_3 <= 128 &&
                h.KT_4 >= 0 && h.KT_4 <= 128, "unsupported layout limits");
        require(c4_validate_header(&h) == 0, "frozen header failed validation");
        require(h.S_r <= (I(1) << 40), "centering scale exceeds arithmetic range");
        // Conservative bounds for the two pair maps, their inverse lists, the
        // fourth level's triangular index tables and the smaller level arrays.
        size_t hu = 0, cu = 0;
        for (int k = 1; k <= h.KHU; ++k) hu += std::max(0,h.DHU-1-k);
        for (int k = 2; k <= h.KCU; ++k) cu += std::max(0,h.DCU-1-k);
        size_t n4 = size_t(h.D0_4)*(h.D0_4-1)/2 + h.KT_4+2;
        size_t scratch = 16*(hu*hu+hu*cu) + 8*n4*n4 + (size_t(64) << 20);
        require(bytes+scratch <= (size_t(5) << 30), "layout and payload exceed workspace budget");
        c4_init_layout(&c, &h);
        c4_compute_sizes(&c, &z);
        require(bytes == z.total_bytes, "payload size mismatch");
        const I *q = reinterpret_cast<const I *>(static_cast<const char *>(mapping) + sizeof h);
        auto take = [&](size_t n) { const I *p = q; q += n; return p; };
        for (auto &p : cent) p = take(c.L2.INNER_N);
        take(z.n_fy); take(z.n_f);
        weight[Edge] = take(z.n_g);
        weight[V] = take(z.n_hV); weight[C] = take(z.n_hC);
        weight[P] = take(z.n_hP4); weight[F] = take(z.n_hF4);
        weight[Y] = take(z.n_hY4); weight[T] = take(z.n_hT4);
        for (auto &p : e23) p = take(z.n_env_23 / 4);
        for (auto &p : ep) p = take(z.n_env_p4 / 4);
        for (auto &p : ef) p = take(z.n_env_f4 / 4);
        for (auto &p : eyv) p = take(z.n_env_y4_hv / 2);
        for (auto &p : eyc) p = take(z.n_env_y4_hc / 2);
        require(reinterpret_cast<const char *>(q) == static_cast<const char *>(mapping) + bytes,
                "payload cursor mismatch");
        const I lim = I(1) << 52, clim = I(1) << 40;
        for (int i = 0; i < 5*c.L2.INNER_N; ++i)
            require(cent[0][i] >= -clim && cent[0][i] <= clim, "centering value outside arithmetic range");
        // This single scan includes every generated weight used below (g supplies separate source checks).
        for (const I *p = weight[V]; p < q; ++p)
            if (*p < -lim || *p > lim) throw std::runtime_error("weight or envelope outside arithmetic range");
        auto radii = [&](const I **table, int n, size_t len) {
            for (int j = 1; j < n; j += 2)
                for (size_t i = 0; i < len; ++i)
                    if (table[j][i] < 0) throw std::runtime_error("negative source envelope radius");
        };
        radii(e23,4,z.n_env_23/4); radii(ep,4,z.n_env_p4/4); radii(ef,4,z.n_env_f4/4);
        radii(eyv,2,z.n_env_y4_hv/2); radii(eyc,2,z.n_env_y4_hc/2);
        sign = h.mode == 1 ? 1 : -1;
        s3 = states(c.L3); s4 = states(c.L4);
        init_active(c.L3, hu3, cu3);
        init_active(c.L4, hu4, cu4);
        yi4.assign(c.L4.NG, -1); ti4.assign(c.L4.NG, -1);
        for (int i = 0; i < c.NY4_ACT; ++i) yi4[c.y4_u[i]] = i;
        for (int i = 0; i < c.NT4_ACT; ++i) ti4[c.t4_u[i]] = i;
        hu_state.assign(c.NHU_ACT, -1); cu_state.assign(c.NCU_ACT, -1);
        for (int i = 0; i < c.L3.INNER_N; ++i) {
            if (hu3[i] >= 0) hu_state[hu3[i]] = i;
            if (cu3[i] >= 0) cu_state[cu3[i]] = i;
        }
    }
    Snapshot(const Snapshot &) = delete;
    ~Snapshot() { c4_free_layout(&c); munmap(mapping, bytes); }

    std::vector<State> states(const Cert4vLevel &g) const {
        std::vector<State> a;
        for (int i = 0; i < g.NG; ++i) {
            int t = g.gKind[i], l = g.gL[i], k = g.gK[i], d = g.gD[i];
            int up = t != 0 ? g.INF_ID :
                     (d+1 < g.D0 ? c4_lvl_inner(&g,l+1,k) : c4_lvl_tail(&g,k));
            int bud = t == 0 && l > 0 ? c4_lvl_inner(&g,l-1,k+1) :
                      (t == 2 && k <= g.KT ? c4_lvl_tail(&g,k+1) : i);
            // All centered contacts have degree below their support cutoff, hence below L2.D0.
            int center = t == 0 ? c4_lvl_inner(&c.L2,l,k) : -1;
            a.push_back({d,k,l,t,up,bud,center,project(g,c.L4,i)});
        }
        return a;
    }
    void init_active(const Cert4vLevel &g, std::vector<int> &hu, std::vector<int> &cu) const {
        hu.assign(g.NG, -1); cu.assign(g.NG, -1);
        for (int i = 0; i < g.INNER_N; ++i) {
            int a = c4_hu_idx(&c,g.gL[i],g.gK[i]), b = c4_cu_idx(&c,g.gL[i],g.gK[i]);
            if (a >= 0) hu[i] = c.hu_to_act[a];
            if (b >= 0) cu[i] = c.cu_to_act[b];
        }
    }
    const Cert4vLevel &level(Family f) const { return f <= C ? c.L3 : c.L4; }
    const std::vector<State> &states(Family f) const { return f <= C ? s3 : s4; }

    int active(Family f, const Labels &x) const {
        switch (f) {
        case V: return hu3[x[1]];
        case C: return cu3[x[0]];
        case P: {
            int a = hu4[x[1]], b = hu4[x[2]];
            return a < 0 || b < 0 ? -1 : c.p4_map[a*c.NHU_ACT+b];
        }
        case F: {
            int a = cu4[x[0]], b = hu4[x[2]];
            return a < 0 || b < 0 ? -1 : c.f4_map[a*c.NHU_ACT+b];
        }
        case Y: return yi4[x[1]];
        case T: return ti4[x[0]];
        default: throw std::logic_error("unexpected shape");
        }
    }
    I value(Family f, int act, const Labels &x) const {
        if (act < 0) return 0;  // all jumps across the exact support boundary land at zero
        size_t a = act, n = level(f).NG, off{};
        switch (f) {
        case V: off = (size_t(c.L3.p_to_pi[x[0]])*c.NHU_ACT+a)*n+x[2]; break;
        case C: off = a*c.L3.pair_n+c4_lvl_pair(&c.L3,x[1],x[2]); break;
        case P: off = (a*c.L4.npars+c.L4.p_to_pi[x[0]])*n+x[3]; break;
        case F: off = (a*n+x[1])*n+x[3]; break;
        case Y: off = (a*c.L4.npars+c.L4.p_to_pi[x[0]])*c.L4.pair_n+c4_lvl_pair(&c.L4,x[2],x[3]); break;
        case T: off = a*c.L4.tri_n+c4_lvl_tri(&c.L4,x[1],x[2],x[3]); break;
        default: throw std::logic_error("unexpected table");
        }
        return sign * weight[f][off];
    }
    I value(Family f, const Labels &x) const { return value(f, active(f,x), x); }
};
} // namespace mh
