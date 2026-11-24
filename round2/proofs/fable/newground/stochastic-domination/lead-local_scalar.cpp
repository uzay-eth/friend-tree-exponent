// Lead's LOCAL scalar certificate of stochastic domination between the adjacent labels u, u+1 of the random friend tree.
//
// Criterion (proved in COUPLING.md, section 4.2).  For a tree T let r(w) = sum over neighbours x of w of 1/D_x and
// rho(w) = r(w)/(D_w (D_w + 1)); for a leaf of a core p this is 1/(2 D_p).  For K >= 2 put
//    dem_K(T) = sum over cores x with D_x >= K of [ sum over core neighbours z of (rho(z) - rho(x))^+  +  L_x (1/(2 D_x) - rho(x))^+ ],
//    sup_K(T) = sum over cores p with L_p - 1 >= K of L_p (rho(p) - 1/(2 D_p))^+
//             + sum over core-core edges {x, z} with rho(z) > rho(x) and L_z >= K of (rho(z) - rho(x)).
// If E dem_K(T_(u-1)) <= E sup_K(T_(u-1)) for every K >= 2, then Z_u stochastically dominates Z_(u+1).
// The program computes both expectations EXACTLY (law of the unrooted tree T_(u-1) on all free trees of that size,
// integer masses over a common denominator; the local terms as integers over lcm(1..n) lcm(1..n+1)).
//
// Usage: local_scalar ULO UHI   (EXT=1 in the environment: extend the law of T_(u-2) by one step instead of storing T_(u-1))
// Build: g++ -O2 -std=c++17 [-DNLIMBS=7] local_scalar.cpp -lgmpxx -lgmp -o local_scalar
//        (NLIMBS=7 is needed when the tree size u-1 is 24 or more)
#include <bits/stdc++.h>
#include <gmpxx.h>
using namespace std;
typedef uint64_t u64; typedef unsigned __int128 u128;
#ifndef NLIMBS
#define NLIMBS 6
#endif
static const int MAXN = 28, LIMBS = NLIMBS, K = 40;

struct Big { u64 l[LIMBS]; };
static inline Big bzero() { Big a; for (int i = 0; i < LIMBS; i++) a.l[i] = 0; return a; }
static inline void badd(Big& a, const Big& b) { u128 c = 0; for (int i = 0; i < LIMBS; i++) { c += (u128)a.l[i] + b.l[i]; a.l[i] = (u64)c; c >>= 64; } }
static inline Big bmul(const Big& a, u64 m) { Big r; u128 c = 0; for (int i = 0; i < LIMBS; i++) { c += (u128)a.l[i] * m; r.l[i] = (u64)c; c >>= 64; } return r; }
static inline bool bneg(const Big& a) { return a.l[LIMBS - 1] >> 63; }
static inline bool biszero(const Big& a) { for (int i = 0; i < LIMBS; i++) if (a.l[i]) return false; return true; }
static inline Big bnegate(const Big& a) { Big r; u128 c = 1; for (int i = 0; i < LIMBS; i++) { c += (u128)(~a.l[i]); r.l[i] = (u64)c; c >>= 64; } return r; }
static mpz_class tompz(const Big& a) { bool n = bneg(a); Big b = n ? bnegate(a) : a; mpz_class z; mpz_import(z.get_mpz_t(), LIMBS, -1, 8, 0, 0, b.l); return n ? mpz_class(-z) : z; }

struct Table {
    vector<u64> keys; vector<Big> vals; size_t cap = 0, cnt = 0; int shards = 1, shard = 0;      // optional filter: keep keys of one hash class only
    void init(size_t want) { cap = 1024; while (cap * 3 < want * 4) cap <<= 1; keys.assign(cap, 0); vals.assign(cap, bzero()); cnt = 0; }
    static inline u64 h(u64 x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
    void add(u64 k, const Big& v) {
        if (shards > 1 && (int)((h(k) >> 40) % (u64)shards) != shard) return;
        size_t i = h(k) & (cap - 1);
        while (true) {
            if (keys[i] == k) { badd(vals[i], v); return; }
            if (keys[i] == 0) { keys[i] = k; vals[i] = v; if (++cnt * 10 > cap * 8) grow(); return; }
            i = (i + 1) & (cap - 1);
        }
    }
    void grow() {
        vector<u64> ok; vector<Big> ov; ok.swap(keys); ov.swap(vals); size_t oc = cap; cap <<= 1; keys.assign(cap, 0); vals.assign(cap, bzero()); cnt = 0;
        for (size_t i = 0; i < oc; i++) if (ok[i]) { size_t j = h(ok[i]) & (cap - 1); while (keys[j]) j = (j + 1) & (cap - 1); keys[j] = ok[i]; vals[j] = ov[i]; cnt++; }
    }
    void clear() { vector<u64>().swap(keys); vector<Big>().swap(vals); cap = cnt = 0; }
    mpz_class total() const { mpz_class s = 0; for (size_t i = 0; i < cap; i++) if (keys[i]) s += tompz(vals[i]); return s; }
    size_t nonzero() const { size_t c = 0; for (size_t i = 0; i < cap; i++) if (keys[i] && !biszero(vals[i])) c++; return c; }
};

struct Tree { int n; int deg[MAXN]; int adj[MAXN][MAXN]; };
static void decode(u64 code, int n, Tree& t) {
    t.n = n; for (int i = 0; i <= n; i++) t.deg[i] = 0;
    int st[MAXN], sp = 0, nxt = 0;
    for (int pos = 2 * n - 1; pos >= 0; --pos) {
        if ((code >> pos) & 1) { int v = nxt++; if (sp > 0) { int p = st[sp - 1]; t.adj[p][t.deg[p]++] = v; t.adj[v][t.deg[v]++] = p; } st[sp++] = v; }
        else sp--;
    }
    if (nxt != n || sp != 0) { fprintf(stderr, "decode error\n"); exit(1); }
}
struct Code { u64 bits; int len; };
static Code canon(const Tree& t, int v, int p) {
    u64 kids[MAXN]; int k = 0;
    for (int i = 0; i < t.deg[v]; i++) { int w = t.adj[v][i]; if (w == p) continue; Code c = canon(t, w, v); kids[k++] = ((u64)c.len << 56) | c.bits; }
    sort(kids, kids + k);
    u64 bits = 1; int len = 1;
    for (int i = 0; i < k; i++) { int cl = (int)(kids[i] >> 56); bits = (bits << cl) | (kids[i] & ((1ULL << 56) - 1)); len += cl; }
    return { bits << 1, len + 1 };
}
// ---- unrooted trees (carrier of the law of T_n): code of the tree rooted at its centre; for two centres a virtual
// root above both (two more bits, so the two kinds cannot collide) ----
static u64 ucanon(const Tree& t) {
    int n = t.n, d[MAXN], q[MAXN], lay[MAXN], h = 0, e = 0, left = n;
    for (int v = 0; v < n; v++) { d[v] = t.deg[v]; lay[v] = -1; if (d[v] <= 1) { q[e++] = v; lay[v] = 0; } }
    int last = 0;
    while (h < e) { int v = q[h++]; last = lay[v]; left--; for (int j = 0; j < t.deg[v]; j++) { int w = t.adj[v][j]; if (lay[w] < 0 && --d[w] == 1) { lay[w] = lay[v] + 1; q[e++] = w; } } }
    (void)left;
    int c1 = -1, c2 = -1;                                          // vertices of the last layer: one or two centres
    for (int v = 0; v < n; v++) if (lay[v] == last) { if (c1 < 0) c1 = v; else c2 = v; }
    if (c2 < 0) return canon(t, c1, -1).bits;
    Code a = canon(t, c1, c2), b = canon(t, c2, c1);
    u64 x = ((u64)a.len << 56) | a.bits, y = ((u64)b.len << 56) | b.bits; if (x > y) { Code tmp = a; a = b; b = tmp; }
    u64 bits = 1; bits = (bits << a.len) | a.bits; bits = (bits << b.len) | b.bits; return bits << 1;
}
static void udecode(u64 code, int n, Tree& t) {
    if ((code >> (2 * n)) == 0) { decode(code, n, t); return; }   // one centre: an ordinary rooted word of 2n bits
    Tree r; decode(code, n + 1, r);                                // two centres: vertex 0 is the virtual root
    if (r.deg[0] != 2) { fprintf(stderr, "udecode error\n"); exit(1); }
    t.n = n; for (int v = 0; v < n; v++) t.deg[v] = 0;
    for (int v = 1; v <= n; v++) for (int j = 0; j < r.deg[v]; j++) { int w = r.adj[v][j]; if (w != 0) t.adj[v - 1][t.deg[v - 1]++] = w - 1; }
    int a = r.adj[0][0] - 1, b = r.adj[0][1] - 1; t.adj[a][t.deg[a]++] = b; t.adj[b][t.deg[b]++] = a;
}
static u64 lcm_upto(int m) { u64 L = 1; for (int j = 2; j <= m; j++) L = L / std::__gcd(L, (u64)j) * j; return L; }

// one growth step from size n to n+1; rootnew = false keeps the root, true roots the result at the new leaf
// one growth step from size n to n+1.  srcun: the source keys are unrooted codes.  mode 0: keep the root (rooted source
// only); mode 1: root the result at the new leaf; mode 2: unrooted result
static void propagate(const Table& src, int n, Table& dst, int mode, bool negate = false, bool srcun = false) {
    u64 J = lcm_upto(n - 1); Tree t; u64 inv[MAXN];
    for (size_t i = 0; i < src.cap; i++) {
        if (!src.keys[i] || biszero(src.vals[i])) continue;
        if (srcun) udecode(src.keys[i], n, t); else decode(src.keys[i], n, t);
        for (int v = 0; v < n; v++) inv[v] = J / t.deg[v];
        for (int w = 0; w < n; w++) {
            u64 wt = 0; for (int j = 0; j < t.deg[w]; j++) wt += inv[t.adj[w][j]];
            t.adj[w][t.deg[w]++] = n; t.adj[n][0] = w; t.deg[n] = 1; t.n = n + 1;
            u64 key = mode == 2 ? ucanon(t) : canon(t, mode == 1 ? n : 0, -1).bits;
            Big m = bmul(src.vals[i], wt); dst.add(key, negate ? bnegate(m) : m);
            t.deg[w]--; t.n = n;
        }
    }
}


typedef __int128 i128;
// local terms of one tree with n vertices, scaled by L1 L2: demand by level k = D_x, supply by level j
static void local_terms(const Tree& t, int n, u64 L1, u64 L2, i128* demv, i128* supv) {
    int Lf[MAXN]; i128 rho[MAXN];
    for (int b = 0; b <= n + 1; b++) { demv[b] = 0; supv[b] = 0; }
    for (int v = 0; v < n; v++) { int c = 0; for (int j = 0; j < t.deg[v]; j++) c += (t.deg[t.adj[v][j]] == 1); Lf[v] = c; }
    for (int v = 0; v < n; v++) { u64 r = 0; for (int j = 0; j < t.deg[v]; j++) r += L1 / t.deg[t.adj[v][j]]; rho[v] = (i128)r * (L2 / ((u64)t.deg[v] * (t.deg[v] + 1))); }
    for (int x = 0; x < n; x++) {
        int D = t.deg[x]; if (D < 2) continue;                                         // cores only
        i128 leafrho = (i128)(L1 / D) * (L2 / 2);                                       // rho of a leaf of x, times L1 L2
        if (Lf[x] > 0) { if (leafrho > rho[x]) demv[D] += (i128)Lf[x] * (leafrho - rho[x]); else supv[Lf[x] - 1] += (i128)Lf[x] * (rho[x] - leafrho); }
        for (int j = 0; j < D; j++) { int z = t.adj[x][j]; if (t.deg[z] < 2 || z < x) continue;      // each core-core edge once
            if (rho[z] > rho[x]) { demv[D] += rho[z] - rho[x]; supv[Lf[z]] += rho[z] - rho[x]; }
            else if (rho[x] > rho[z]) { demv[t.deg[z]] += rho[x] - rho[z]; supv[Lf[x]] += rho[x] - rho[z]; } }
    }
}
int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s ULO UHI   (environment EXT=1: do not store the law of T_(u-1); extend the law of T_(u-2) by one step instead)\n", argv[0]); return 2; }
    int ulo = max(atoi(argv[1]), 4), uhi = atoi(argv[2]); bool ext = getenv("EXT") != nullptr;
    vector<mpz_class> Den(MAXN + 1); Den[2] = 1; for (int n = 2; n < MAXN; n++) Den[n + 1] = Den[n] * n * mpz_class((unsigned long)lcm_upto(n - 1));
    static const size_t FR[] = {0, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159, 7741, 19320, 48629, 123867, 317955, 823065, 2144505, 5623756, 14828074, 39299897, 104636890};
    Table cur; cur.init(4); { Big one = bzero(); one.l[0] = 1; cur.add(0b110100ULL, one); } int have = 2;
    for (int u = ulo; u <= uhi; u++) {
        int n = u - 1, stored = ext ? n - 1 : n;
        if (Den[n] >= (mpz_class(1) << (64 * LIMBS - 2))) { fprintf(stderr, "integers too short for size %d: rebuild with -DNLIMBS=7\n", n); return 1; }
        while (have < stored) {
            Table nx; nx.init(FR[have + 1]); propagate(cur, have, nx, 2, false, true);
            if (nx.total() != Den[have + 1] || nx.cnt != FR[have + 1]) { fprintf(stderr, "mass or census check failed at size %d\n", have + 1); return 1; }
            cur.clear(); cur = std::move(nx); have++;
        }
        if (have != stored) { fprintf(stderr, "sizes must increase\n"); return 1; }
        u64 L1 = lcm_upto(n), L2 = lcm_upto(n + 1), J = lcm_upto(stored - 1);
        vector<mpz_class> demLo(n + 2, 0), demHi(n + 2, 0), supLo(n + 2, 0), supHi(n + 2, 0); Tree t; mpz_class mass, total = 0;
        auto addto = [&](vector<mpz_class>& lo, vector<mpz_class>& hi, int b, i128 v) {
            if (v == 0) return; unsigned __int128 w = (unsigned __int128)v;
            mpz_addmul_ui(lo[b].get_mpz_t(), mass.get_mpz_t(), (unsigned long)(u64)w); mpz_addmul_ui(hi[b].get_mpz_t(), mass.get_mpz_t(), (unsigned long)(u64)(w >> 64)); };
        i128 demv[MAXN + 2], supv[MAXN + 2], dacc[MAXN + 2], sacc[MAXN + 2];
        for (size_t i = 0; i < cur.cap; i++) {
            if (!cur.keys[i]) continue;
            udecode(cur.keys[i], stored, t); mass = tompz(cur.vals[i]);
            if (!ext) { local_terms(t, n, L1, L2, dacc, sacc); total += mass; }
            else {                                                                      // one more step: all attachments
                for (int b = 0; b <= n + 1; b++) { dacc[b] = 0; sacc[b] = 0; }
                u64 inv[MAXN]; for (int v = 0; v < stored; v++) inv[v] = J / t.deg[v]; u64 wsum = 0;
                for (int w = 0; w < stored; w++) {
                    u64 wt = 0; for (int j = 0; j < t.deg[w]; j++) wt += inv[t.adj[w][j]]; wsum += wt;
                    t.adj[w][t.deg[w]++] = stored; t.adj[stored][0] = w; t.deg[stored] = 1; t.n = n;
                    local_terms(t, n, L1, L2, demv, supv);
                    for (int b = 0; b <= n; b++) { dacc[b] += demv[b] * (i128)wt; sacc[b] += supv[b] * (i128)wt; }
                    t.deg[w]--; t.n = stored;
                }
                total += mass * mpz_class((unsigned long)wsum);
            }
            for (int b = 0; b <= n; b++) { addto(demLo, demHi, b, dacc[b]); addto(supLo, supHi, b, sacc[b]); }
        }
        if (total != Den[n]) { fprintf(stderr, "total mass check failed for size %d\n", n); return 1; }
        mpz_class cd = 0, cs = 0; bool ok = true; mpq_class worst(0); int worstK = 0;
        for (int Kk = n; Kk >= 2; Kk--) {
            cd += demLo[Kk] + (demHi[Kk] << 64); cs += supLo[Kk] + (supHi[Kk] << 64);
            if (cd > cs) ok = false;
            if (cd > 0) { if (cs == 0) { worst = -1; worstK = Kk; } else { mpq_class r(cd, cs); r.canonicalize(); if (worst >= 0 && r > worst) { worst = r; worstK = Kk; } } }
        }
        mpz_class scale = Den[n] * mpz_class((unsigned long)L1) * mpz_class((unsigned long)L2);
        mpq_class dem2(cd, scale), sup2(cs, scale); dem2.canonicalize(); sup2.canonicalize();
        printf("pair (%d,%d): law of T_%d %s; E dem_2 = %.9f, E sup_2 = %.9f (in units of 1/(n(n+1))); largest ratio E dem_K / E sup_K over K >= 2: %s at K = %d: %s\n",
               u, u + 1, n, ext ? "by one step from all free trees of the size below" : "on all free trees of that size", dem2.get_d(), sup2.get_d(), worst >= 0 ? to_string(worst.get_d()).c_str() : "infinite", worstK,
               ok ? "PASS (stochastic domination certified by the local criterion)" : "local criterion not satisfied");
        fflush(stdout);
    }
    return 0;
}
