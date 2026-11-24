// Third program for board claim 3(a), written by the lead independently of the child's unrooted-orbit program.
//
// For adjacent labels u, u+1 of the random friend tree it propagates the SIGNED difference of the two rooted laws
//     q_u(s) - q_{u+1}(s),   q_v(s) = P( (T_n, v) is the rooted shape s ),
// exactly (integers over a common denominator) from size u+1 upwards, and at each size evaluates
//     p_u - p_{u+1} >= sum_s [ (q_u - q_{u+1})^+ lb(s) - (q_u - q_{u+1})^- ub(s) ]          (hub probabilities)
//     m_u - m_{u+1} >= sum_s [ (q_u - q_{u+1})^+ max(L-1,0)/n - (q_u - q_{u+1})^- D/n ]      (means of Z)
// with the lead's baseline classifier bounds lb, ub of hub_order_exact.py (same formulas, re-implemented here).
//
// Representation: a rooted shape is its balanced-parenthesis word with children sorted (2 bits per vertex, one machine
// word); the root is the marked label.  No automorphism orbits are needed: a step attaches a leaf at each of the n
// vertices of a representative and re-canonicalises; equal shapes merge in a hash table.
// The law of the newborn label v is obtained from the law of the unrooted tree T_{v-1} (keys: the tree rooted at its
// centre, with a virtual root above two centres) by attaching the new leaf and rooting the tree at it.  Labels 1 and 2
// have the same law (swap of the two initial vertices).  (Version 1 of this program carried the law of (T_n, label 1)
// instead; both versions give identical output.)
//
// Masses are 384-bit two's-complement integers over Den_n = prod_{j=2}^{n-1} j * lcm(1..j-1); classifier values are
// exact GMP rationals.  Total masses are checked against Den_n at every step.
//
// Usage: signed_orbit_dp ULO UHI EXTRA [NMAXSIZE]
//   for each u in [ULO, UHI]: sizes u+1 .. min(u+1+EXTRA, NMAXSIZE)
//   environment: IMPROVED=1 strengthened bounds; NOEDGE=1 without the edge theorem; SHARDS=k (with EXTRA = 0) builds
//   the table of size u+1 in k passes over disjoint hash classes of the shapes, to save memory
// Build: g++ -O2 -std=c++17 signed_orbit_dp.cpp -lgmpxx -lgmp -o signed_orbit_dp
#include <bits/stdc++.h>
#include <gmpxx.h>
using namespace std;
typedef uint64_t u64; typedef unsigned __int128 u128;
static const int MAXN = 26, LIMBS = 6, K = 40;

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

// ---------------- classifier (exact rationals, interned) ----------------
static vector<mpq_class> VAL; static map<string, int> VIDX;
static int intern(mpq_class q) { q.canonicalize(); string s = q.get_str(); auto it = VIDX.find(s); if (it != VIDX.end()) return it->second; VAL.push_back(q); return VIDX[s] = (int)VAL.size() - 1; }
static int ZERO, ONE;
static map<pair<int, int>, int> Bc, Fc, MAXc; static map<int, int> OMc; static map<vector<int>, int> BAREc; static map<int, vector<mpq_class>> FT;
static int Bid(int l, int c) {
    if (l == 0) return ZERO;
    auto key = make_pair(l, c); auto it = Bc.find(key); if (it != Bc.end()) return it->second;
    vector<mpq_class> a; a.push_back(1); for (int k = 1; k < 40; k++) a.push_back(a.back() / (c + k));
    mpq_class tail = a.back() / (c + 40) * mpq_class(c + 41, c + 40), part = 0, all = 0;
    for (int k = 0; k < 40; k++) { all += a[k]; if (k < l) part += a[k]; }
    return Bc[key] = intern(part / (all + tail));
}
static const vector<mpq_class>& Ftab(int c) {
    auto it = FT.find(c); if (it != FT.end()) return it->second;
    vector<mpq_class> al(K + 2), be(K + 2), F(K + 2);
    for (int k = 1; k <= K; k++) {
        mpq_class d = k * (k + c) + k, up = k * (k + c), lo = k - 1;
        if (k == 1) { al[k] = 1 / d; be[k] = up / d; }
        else { mpq_class den = d - lo * be[k - 1]; al[k] = (1 + lo * al[k - 1]) / den; be[k] = up / den; }
    }
    F[K + 1] = mpq_class(1, K);
    for (int k = K; k >= 1; k--) F[k] = al[k] + be[k] * F[k + 1];
    for (int k = 1; k <= K; k++) if (!(F[k] >= F[k + 1] && F[k] <= 1 && F[k + 1] >= 0)) { fprintf(stderr, "F table not monotone, c=%d\n", c); exit(1); }
    return FT[c] = F;
}
static int Fid(int k, int c) {
    auto key = make_pair(k, c); auto it = Fc.find(key); if (it != Fc.end()) return it->second;
    return Fc[key] = intern(k <= K ? Ftab(c)[k] : mpq_class(1, k - 1));
}
// ---- strengthened classifier (IMPROVED=1): four one-step refinements, each re-derived by the lead (see README) ----
static bool IMPROVED = false; static const int CUT = 24;
static map<pair<int, int>, mpq_class> FPt, LPt; static map<int, vector<mpq_class>> HT; static bool tables_ready = false;
static mpq_class Bq(int l, int c) { return VAL[Bid(l, c)]; }
static void build_tables() {
    if (tables_ready) return;
    tables_ready = true;
    // leaf root, upper bound for "ever promoted": parent state (k leaves incl. the marked one, c core neighbours);
    // a loss (k-1, c+1) keeps the parent's degree.  Stop at degree CUT with the baseline F.
    for (int c = CUT - 1; c >= 0; c--) for (int k = CUT - 1 - c; k >= 1; k--) {
        mpq_class up = k * (k + c), nextk = (k == CUT - 1 - c) ? Ftab(c)[k + 1] : FPt[{k + 1, c}], loss = (k == 1) ? mpq_class(0) : FPt[{k - 1, c + 1}];
        mpq_class v = (nextk * up + 1 + (k - 1) * loss) / (up + k); if (Ftab(c)[k] < v) v = Ftab(c)[k];
        if (!(v >= nextk) || !((up + k) * v >= up * nextk + 1 + (k - 1) * loss) || v > 1) { fprintf(stderr, "FP check failed k=%d c=%d\n", k, c); exit(1); }
        FPt[{k, c}] = v;
    }
    // core root, lower bound for "hub": next change of (k, c) is a birth with probability in [d/(d+1), 1], d = k + c
    for (int d = CUT - 1; d >= 2; d--) for (int k = 1; k <= d; k++) {
        int c = d - k; mpq_class birth = (k + 1 == CUT) ? mpq_class(0) : (d == CUT - 1 ? Bq(k + 1, c) : LPt[{k + 1, c}]), loss = (k == 1) ? mpq_class(0) : LPt[{k - 1, c + 1}];
        mpq_class nx = (d * birth + loss) / (d + 1), v = nx < birth ? nx : birth; if (Bq(k, c) > v) v = Bq(k, c);
        LPt[{k, c}] = v;
    }
}
static const vector<mpq_class>& Htab(int c) {            // bare core: one neighbour with k leaves and c >= 1 core neighbours selects it
    auto it = HT.find(c); if (it != HT.end()) return it->second;
    vector<mpq_class> al(K + 2), be(K + 2), H(K + 2); al[0] = 1; be[0] = 0;
    for (int k = 1; k <= K; k++) { mpq_class up = k * (k + c), den = up + k + 1 - k * be[k - 1]; al[k] = (1 + k * al[k - 1]) / den; be[k] = up / den; }
    H[K + 1] = mpq_class(1, K); for (int k = K; k >= 1; k--) H[k] = al[k] + be[k] * H[k + 1]; H[0] = 1;
    for (int k = 0; k <= K; k++) if (!(H[k] >= H[k + 1] && H[k] <= 1)) { fprintf(stderr, "H table not monotone, c=%d\n", c); exit(1); }
    for (int k = 1; k <= K; k++) if ((k * (k + c) + k + 1) * H[k] != k * (k + c) * H[k + 1] + k * H[k - 1] + 1) { fprintf(stderr, "H equation failed\n"); exit(1); }
    return HT[c] = H;
}
static int oneminus(int id) { auto it = OMc.find(id); if (it != OMc.end()) return it->second; return OMc[id] = intern(1 - VAL[id]); }
static int maxid(int a, int b) { if (a == b) return a; auto key = make_pair(min(a, b), max(a, b)); auto it = MAXc.find(key); if (it != MAXc.end()) return it->second; return MAXc[key] = (VAL[a] >= VAL[b] ? a : b); }
static int leafcnt[MAXN]; static map<pair<int, int>, int> FPc, LPc; static map<int, int> LEAFc;
static int ubid(const Tree& t, int v) {
    if (t.deg[v] == 1) {                                                                            // leaf: promotion bound through its parent's stock
        int w = t.adj[v][0], k = leafcnt[w], c = t.deg[w] - leafcnt[w];
        if (IMPROVED && k + c < CUT) { auto key = make_pair(k, c); auto it = FPc.find(key); if (it != FPc.end()) return it->second; return FPc[key] = intern(FPt.at(key)); }
        return Fid(k, c);
    }
    if (leafcnt[v] >= 1) return ONE;
    vector<int> key;                                                                                // bare core: must receive a leaf first
    for (int j = 0; j < t.deg[v]; j++) { int w = t.adj[v][j]; key.push_back(IMPROVED ? leafcnt[w] * 64 + (t.deg[w] - leafcnt[w]) : max(1, leafcnt[w] - 1)); }
    sort(key.begin(), key.end()); auto it = BAREc.find(key); if (it != BAREc.end()) return it->second;
    mpq_class s = 0;
    for (int x : key) {
        if (!IMPROVED) { s += mpq_class(1, x); continue; }
        int l = x / 64, c = x % 64; mpq_class y(1, max(1, l - 1)); if (l <= K && Htab(c)[l] < y) y = Htab(c)[l]; s += y;
    }
    if (s > 1) s = 1;
    return BAREc[key] = intern(s);
}
static int lbid(const Tree& t, int v, bool edge) {
    int best;
    if (t.deg[v] >= 2) {
        int k = leafcnt[v], c = t.deg[v] - leafcnt[v];
        if (IMPROVED && k >= 1 && k + c < CUT) { auto key = make_pair(k, c); auto it = LPc.find(key); best = it != LPc.end() ? it->second : (LPc[key] = intern(LPt.at(key))); }
        else best = Bid(k, c);
    } else if (IMPROVED) {                                                                          // leaf: selected at its parent's first ring, then B(1,1)
        int d = t.deg[t.adj[v][0]]; auto it = LEAFc.find(d); best = it != LEAFc.end() ? it->second : (LEAFc[d] = intern(Bq(1, 1) / (2 * d)));
    } else best = ZERO;
    if (edge) for (int j = 0; j < t.deg[v]; j++) best = maxid(best, oneminus(ubid(t, t.adj[v][j])));
    return best;
}
static string dec(const mpq_class& q) { char buf[64]; snprintf(buf, sizeof buf, "%+.12f", q.get_d()); return buf; }

struct Acc {
    vector<Big> hubLo, hubHi, mPosL, mNegD, mPosD, mNegL, hallP, hallN; Big universal = bzero(); bool hallok = true; size_t pos = 0, neg = 0;
    Acc() : mPosL(MAXN, bzero()), mNegD(MAXN, bzero()), mPosD(MAXN, bzero()), mNegL(MAXN, bzero()), hallP(MAXN, bzero()), hallN(MAXN, bzero()) {}
};
static void accumulate(const Table& diff, int n, bool edge, Acc& a) {
    Tree t;
    for (size_t i = 0; i < diff.cap; i++) {
        if (!diff.keys[i] || biszero(diff.vals[i])) continue;
        const Big& m = diff.vals[i]; decode(diff.keys[i], n, t);
        for (int v = 0; v < n; v++) { int c = 0; for (int j = 0; j < t.deg[v]; j++) c += (t.deg[t.adj[v][j]] == 1); leafcnt[v] = c; }
        int lb = lbid(t, 0, edge), ub = ubid(t, 0);
        if (a.hubLo.size() < VAL.size()) { a.hubLo.resize(VAL.size(), bzero()); a.hubHi.resize(VAL.size(), bzero()); }
        int L1 = max(leafcnt[0] - 1, 0), D = t.deg[0];
        if (!bneg(m)) { a.pos++; badd(a.hubLo[lb], m); badd(a.hubHi[ub], m); badd(a.mPosL[L1], m); badd(a.mPosD[D], m); }
        else { a.neg++; badd(a.hubLo[ub], m); badd(a.hubHi[lb], m); badd(a.mNegD[D], m); badd(a.mNegL[L1], m); }
        // scalar transport test for stochastic domination (see newground/stochastic-domination): a leaf whose parent has
        // k other neighbours is >= a leaf whose parent has at least k other leaf neighbours; a root with a leaf child is
        // >= every leaf-rooted shape
        if (D >= 2) { if (bneg(m) || leafcnt[0] == 0) a.hallok = false; else badd(a.universal, m); }
        else { int par = t.adj[0][0]; if (!bneg(m)) badd(a.hallP[t.deg[par] - 1], m); else badd(a.hallN[leafcnt[par] - 1], bnegate(m)); }
    }
}
static void finalize(const Acc& a, int n, int u, const mpz_class& Den) {
    mpq_class lo = 0, hi = 0; mpz_class mlo = 0, mhi = 0, pmass = 0;
    for (size_t id = 0; id < a.hubLo.size(); id++) { lo += VAL[id] * tompz(a.hubLo[id]); hi += VAL[id] * tompz(a.hubHi[id]); }
    for (int j = 0; j < MAXN; j++) { mlo += j * (tompz(a.mPosL[j]) + tompz(a.mNegD[j])); mhi += j * (tompz(a.mPosD[j]) + tompz(a.mNegL[j])); pmass += tompz(a.mPosD[j]); }
    lo /= Den; hi /= Den; mpq_class qmlo(mlo, Den * n), qmhi(mhi, Den * n), tvq(pmass, Den); qmlo.canonicalize(); qmhi.canonicalize(); tvq.canonicalize();
    printf("H %d %d %s %s %s %s\n", u, n, dec(lo).c_str(), dec(hi).c_str(), lo.get_str().c_str(), hi.get_str().c_str());
    printf("M %d %d %s %s %s %s\n", u, n, dec(qmlo).c_str(), dec(qmhi).c_str(), qmlo.get_str().c_str(), qmhi.get_str().c_str());
    printf("# u=%d size=%d shapes+ %zu shapes- %zu  total variation %s  hub gap %s  mean gap %s\n", u, n, a.pos, a.neg, dec(tvq).c_str(),
           lo > 0 ? "CERTIFIED > 0" : "not decided", qmlo > 0 ? "CERTIFIED > 0" : "not decided");
    // nested Hall condition for the scalar rule: for every K, positive leaf-rooted mass with k >= K <= negative mass with j >= K
    { mpz_class cp = 0, cn = 0, uni = tompz(a.universal), tot = 0; bool ok = a.hallok; mpq_class worst(-1); int worstK = -1;
      for (int K = MAXN - 1; K >= 0; K--) { cp += tompz(a.hallP[K]); cn += tompz(a.hallN[K]); if (cp > cn) ok = false;
          if (cp > 0) { mpq_class r(cp, cn > 0 ? cn : mpz_class(1)); r.canonicalize(); if (cn == 0 || r > worst) { worst = r; worstK = K; } } }
      tot = cp; mpq_class frac(tot, uni > 0 ? uni : mpz_class(1)); frac.canonicalize();
      printf("T %d %d scalar transport rule: positive leaf-rooted mass %s, negative mass %s, universal %s (ratio to universal %s); largest ratio of cumulative masses %s at K = %d: %s\n",
             u, n, tot.get_str().c_str(), cn.get_str().c_str(), uni.get_str().c_str(), dec(frac).c_str(), worst >= 0 ? dec(worst).c_str() : "none", worstK,
             (ok && cn == tot + uni) ? "PASS (stochastic domination certified by the scalar rule)" : "not certified by the scalar rule"); }
    fflush(stdout);
}

int main(int argc, char** argv) {
    if (argc < 4) { fprintf(stderr, "usage: %s ULO UHI EXTRA [NMAXSIZE]\n", argv[0]); return 2; }
    int ulo = atoi(argv[1]), uhi = atoi(argv[2]), extra = atoi(argv[3]), nmaxsize = argc > 4 ? atoi(argv[4]) : 99;
    bool edge = getenv("NOEDGE") == nullptr;
    if (ulo < 2 || uhi < ulo || min(uhi + 1 + extra, nmaxsize) > MAXN - 3) { fprintf(stderr, "range\n"); return 2; }
    ZERO = intern(0); ONE = intern(1); IMPROVED = getenv("IMPROVED") != nullptr; if (IMPROVED) build_tables();
    printf("# classifier: %s%s\n", IMPROVED ? "strengthened" : "baseline", edge ? "" : ", edge theorem not used"); 
    vector<mpz_class> Den(MAXN + 1); Den[2] = 1; for (int n = 2; n < MAXN; n++) Den[n + 1] = Den[n] * n * mpz_class((unsigned long)lcm_upto(n - 1));
    static const size_t R[] = {0, 1, 1, 2, 4, 9, 20, 48, 115, 286, 719, 1842, 4766, 12486, 32973, 87811, 235381, 634847, 1721159, 4688676, 12826228, 35221832, 97055181, 268282855, 743724984};   // rooted trees (OEIS A000081), size hints only
    static const size_t FR[] = {0, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159, 7741, 19320, 48629, 123867, 317955, 823065, 2144505, 5623756, 14828074, 39299897};   // free trees (OEIS A000055), size hints only
    int S = getenv("SHARDS") ? atoi(getenv("SHARDS")) : 1; if (S < 1 || (S > 1 && extra > 0)) { fprintf(stderr, "SHARDS needs EXTRA = 0\n"); return 2; }
    // law of the unrooted tree T_n, computed as needed
    int top = max(uhi, 3); vector<Table> U(top + 2); int have = 2;
    U[2].init(4); { Big one = bzero(); one.l[0] = 1; U[2].add(0b110100ULL, one); }
    auto need = [&](int n) {
        while (have < n) {
            U[have + 1].init(FR[have + 1]); propagate(U[have], have, U[have + 1], 2, false, true);
            if (U[have + 1].total() != Den[have + 1]) { fprintf(stderr, "mass check failed at size %d\n", have + 1); exit(1); }
            have++; fprintf(stderr, "law of T_%d: %zu unrooted shapes (expected %zu)\n", have, U[have].cnt, FR[have]);
            if (U[have].cnt != FR[have]) { fprintf(stderr, "census of free trees failed\n"); exit(1); }
        }
    };
    for (int u = ulo; u <= uhi; u++) {
        Table Qu;
        if (u >= 3) {
            need(u - 1); Qu.init(R[u - 1]); propagate(U[u - 1], u - 1, Qu, 1, false, true);          // newborn u at size u
            if (Qu.total() != Den[u]) { fprintf(stderr, "newborn mass check failed (u=%d)\n", u); return 1; }
            if (Qu.cnt != R[u - 1]) { fprintf(stderr, "census of rooted trees failed (u=%d)\n", u); return 1; }
        }
        need(u);
        Acc acc; mpz_class newbornTotal = 0, allTotal = 0; Table A; int n = u + 1;
        for (int sh = 0; sh < S; sh++) {
            A.clear(); A.init((R[u] + R[u - 1]) / S + (S > 1 ? (R[u] + R[u - 1]) / (8 * S) : 0) + 1024); A.shards = S; A.shard = sh;
            propagate(U[u], u, A, 1, true, true);                                                   // minus newborn u+1 at size u+1
            newbornTotal += A.total();
            if (u == 2) { Big one = bzero(); one.l[0] = 1; A.add(0b110100ULL, one); A.add(0b111000ULL, one); }   // label 2 at size 3 has the law of label 1: centre or end of the path, 1/2 each
            else propagate(Qu, u, A, 0);                                                            // plus label u at size u+1
            allTotal += A.total();
            if (S > 1) { fprintf(stderr, "pair %d: shard %d of %d, %zu shapes\n", u, sh + 1, S, A.cnt); accumulate(A, n, edge, acc); }
        }
        Qu.clear();
        if (newbornTotal != -Den[u + 1] || allTotal != 0) { fprintf(stderr, "mass check failed (pair %d)\n", u); return 1; }
        if (S > 1) { finalize(acc, n, u, Den[n]); continue; }
        A.shards = 1; { Acc a1; accumulate(A, n, edge, a1); finalize(a1, n, u, Den[n]); }
        while (n < u + 1 + extra && n < nmaxsize) {
            Table Nx; Nx.init(min(R[n + 1], A.nonzero() * (size_t)n)); propagate(A, n, Nx, 0); A.clear(); A = std::move(Nx); n++;
            if (A.total() != 0) { fprintf(stderr, "difference mass not zero\n"); return 1; }
            Acc a1; accumulate(A, n, edge, a1); finalize(a1, n, u, Den[n]);
        }
    }
    return 0;
}
