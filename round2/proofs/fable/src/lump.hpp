// lump.hpp -- exact per-vertex "flux" accumulation for a lumped (Markov-state-model) description of the
// core of the continuous-time friend tree.
//
// A cell scheme is a function  cell(u) = 64-bit code computed from the neighbourhood of core vertex u.
// For a snapshot tree we accumulate, exactly (no sampling of events):
//     n[s]        = number of core vertices in cell s
//     F[s][s']    = sum over vertices u in cell s of the rate at which u's own cell changes to s'
//                   (off-diagonal, >= 0; the matching loss is put on the diagonal F[s][s]),
//                   plus the rate at which u gives birth to a new core vertex whose cell is s'.
// Elementary events of the model (each core vertex v):
//     G_v  "v gains a leaf":      rate L_v + sum_{core nbrs j} 1/D_j     (leaves of v ring; core nbrs ring and pick v)
//     K_v  "v converts a leaf":   rate L_v / D_v                          (v rings and picks a leaf; that leaf
//                                                                          becomes a core vertex with (L,c)=(1,1))
// For the stationary (tilted) law pi of the cell of a typical core vertex, d/dt E[n] = E[n M] with
// M[s][s'] = E_pi[flux | s], so  pi M = mu pi : mu is the Perron root of the lumped matrix, for ANY cell scheme.
#pragma once
#include "ftsim.hpp"
#include <algorithm>
#include <unordered_map>

static inline uint64_t mix64(uint64_t x) {
    x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33;
    return x;
}
static inline uint64_t hcomb(uint64_t h, uint64_t v) { return mix64(h ^ (v + 0x9e3779b97f4a7c15ULL + (h << 6) + (h >> 2))); }

struct Scheme {
    // --- centre ---
    int Lfine = 16;        // L < Lfine: exact; above: geometric bins with ratio 2^(1/Lsub)
    int Lsub = 1;
    int ccap = 8;          // c >= ccap pooled
    // --- neighbours (depth 1) ---
    int R = 0;             // radius: 0, 1 or 2
    int Dcap1 = 6;         // neighbour degree resolved exactly below Dcap1, pooled into log2 bins up to Dbig1, then "BIG"
    int Dbig1 = 6;         // neighbours with D >= Dbig1 are one class (if Dbig1 > Dcap1: log2 bins between)
    int Lnb1 = 0;          // neighbour leaf count resolved as min(L_j, Lnb1) when D_j < Dcap1 (0 = not resolved)
    int Lcut = 1 << 30;    // centre vertices with L >= Lcut get no neighbour detail (cell = own state only)
    int ccut = 1 << 30;    // centre vertices with c >= ccut get no neighbour detail
    // --- second shell (depth 2), used when R == 2: for a depth-1 neighbour j with D_j < Dcap2 its other
    //     core neighbours k are resolved as min(D_k, Dk2) ---
    int Dcap2 = 0;
    int Dk2 = 4;
    int relpar = 0;        // 1: neighbour code also records whether the neighbour is the parent of the centre
};

struct Snap {
    std::vector<uint64_t> L;
    std::vector<uint32_t> c, parent, first_child, next_sib;
    std::vector<double> cred;   // sum over core neighbours of 1/D
    size_t size() const { return L.size(); }
    void build_children() {
        size_t C = L.size();
        first_child.assign(C, NONE); next_sib.assign(C, NONE);
        for (size_t u = C; u-- > 1;) { uint32_t p = parent[u]; next_sib[u] = first_child[p]; first_child[p] = (uint32_t)u; }
    }
    void build_cred() {
        size_t C = L.size(); cred.assign(C, 0.0);
        for (size_t u = 1; u < C; u++) {
            uint32_t p = parent[u];
            cred[u] += 1.0 / (double)(L[p] + c[p]);
            cred[p] += 1.0 / (double)(L[u] + c[u]);
        }
    }
};

struct Lumper {
    const Scheme &S; Snap &T;
    Lumper(const Scheme &s, Snap &t) : S(s), T(t) {}
    inline uint64_t binL(uint64_t L) const {
        if (L < (uint64_t)S.Lfine) return L;
        double x = log2((double)L / (double)S.Lfine) * S.Lsub;
        return (uint64_t)S.Lfine + (uint64_t)floor(x);
    }
    inline uint64_t own(uint32_t u) const {
        uint64_t cc = T.c[u] < (uint32_t)S.ccap ? T.c[u] : (uint32_t)S.ccap;
        return binL(T.L[u]) * 64 + cc;
    }
    // code of neighbour j as seen from u (depth-1 message)
    uint64_t nbcode(uint32_t j, uint32_t u) const {
        uint64_t D = T.L[j] + T.c[j];
        uint64_t code;
        if (D >= (uint64_t)S.Dbig1) code = 1000003;                      // BIG
        else if (D >= (uint64_t)S.Dcap1) { uint64_t b = 0; uint64_t x = D; while (x > 1) { x >>= 1; b++; } code = 2000000 + b; }
        else {
            uint64_t l = S.Lnb1 > 0 ? std::min<uint64_t>(T.L[j], (uint64_t)S.Lnb1) : 0;
            code = D * 256 + l;
            if (S.R >= 2 && D < (uint64_t)S.Dcap2) {
                // second shell: multiset of capped degrees of j's other core neighbours
                uint64_t cnt[16] = {0};
                auto add = [&](uint32_t k) { uint64_t Dk = T.L[k] + T.c[k]; if (Dk > (uint64_t)S.Dk2) Dk = S.Dk2; cnt[Dk < 16 ? Dk : 15]++; };
                if (T.parent[j] != NONE && T.parent[j] != u) add(T.parent[j]);
                for (uint32_t k = T.first_child[j]; k != NONE; k = T.next_sib[k]) if (k != u) add(k);
                uint64_t h = code;
                for (int i = 2; i <= S.Dk2 && i < 16; i++) h = h * 31 + cnt[i];
                code = hcomb(0x1234, h);
            }
        }
        if (S.relpar) code = code * 2 + (T.parent[u] == j ? 1 : 0);
        return code;
    }
    uint64_t cell(uint32_t u) const {
        uint64_t o = own(u);
        if (S.R == 0 || T.L[u] >= (uint64_t)S.Lcut || T.c[u] >= (uint32_t)S.ccut) return mix64(o * 2 + 1);
        uint64_t buf[64]; std::vector<uint64_t> big; uint64_t *a = buf; size_t n = 0;
        size_t cu = T.c[u];
        if (cu > 64) { big.resize(cu); a = big.data(); }
        if (T.parent[u] != NONE) a[n++] = nbcode(T.parent[u], u);
        for (uint32_t k = T.first_child[u]; k != NONE; k = T.next_sib[k]) a[n++] = nbcode(k, u);
        std::sort(a, a + n);
        uint64_t h = mix64(o * 2);
        for (size_t i = 0; i < n; i++) h = hcomb(h, a[i]);
        return h;
    }
    // collect vertices within distance r of v (including v)
    void ball(uint32_t v, int r, std::vector<uint32_t> &out) const {
        out.clear(); out.push_back(v);
        size_t lo = 0;
        for (int d = 0; d < r; d++) {
            size_t hi = out.size();
            for (size_t i = lo; i < hi; i++) {
                uint32_t x = out[i];
                // neighbours of x not already present: tree => only the vertex we came from is a repeat
                if (T.parent[x] != NONE) out.push_back(T.parent[x]);
                for (uint32_t k = T.first_child[x]; k != NONE; k = T.next_sib[k]) out.push_back(k);
            }
            lo = hi;
        }
        std::sort(out.begin(), out.end()); out.erase(std::unique(out.begin(), out.end()), out.end());
    }
};

struct PairHash { size_t operator()(const std::pair<uint64_t, uint64_t> &p) const { return (size_t)hcomb(p.first, p.second); } };
typedef std::unordered_map<std::pair<uint64_t, uint64_t>, double, PairHash> FluxMap;
typedef std::unordered_map<uint64_t, double> CountMap;

// Accumulate n and F for one snapshot.  Also returns plain sums: C, S = sum L/D.
static void accumulate(const Scheme &S, Snap &T, CountMap &n, FluxMap &F, double &Csum, double &Ssum) {
    Lumper lp(S, T);
    size_t C = T.size();
    T.build_cred();
    std::vector<uint64_t> cellv(C);
    for (size_t u = 0; u < C; u++) { cellv[u] = lp.cell((uint32_t)u); n[cellv[u]] += 1.0; }
    Csum += (double)C;
    std::vector<uint32_t> bl;
    int R = S.R;
    for (size_t v = 0; v < C; v++) {
        double D = (double)(T.L[v] + T.c[v]);
        double gain = (double)T.L[v] + T.cred[v];
        double conv = (double)T.L[v] / D;
        Ssum += conv;
        lp.ball((uint32_t)v, R, bl);
        // --- event G_v ---
        if (gain > 0) {
            T.L[v] += 1;
            for (uint32_t u : bl) {
                uint64_t s2 = lp.cell(u);
                if (s2 != cellv[u]) { F[{cellv[u], s2}] += gain; F[{cellv[u], cellv[u]}] -= gain; }
            }
            T.L[v] -= 1;
        }
        // --- event K_v ---
        if (conv > 0) {
            uint32_t w = (uint32_t)T.L.size();
            T.L[v] -= 1; T.c[v] += 1;
            T.L.push_back(1); T.c.push_back(1); T.parent.push_back((uint32_t)v);
            T.first_child.push_back(NONE); T.next_sib.push_back(T.first_child[v]); T.first_child[v] = w;
            for (uint32_t u : bl) {
                uint64_t s2 = lp.cell(u);
                if (s2 != cellv[u]) { F[{cellv[u], s2}] += conv; F[{cellv[u], cellv[u]}] -= conv; }
            }
            F[{cellv[v], lp.cell(w)}] += conv;   // birth, attributed to the parent's cell
            // revert
            T.first_child[v] = T.next_sib[w];
            T.L.pop_back(); T.c.pop_back(); T.parent.pop_back(); T.first_child.pop_back(); T.next_sib.pop_back();
            T.L[v] += 1; T.c[v] -= 1;
        }
    }
}
