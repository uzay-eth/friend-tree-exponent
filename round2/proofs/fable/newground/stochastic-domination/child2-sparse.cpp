// Lead's independent exact certificate of stochastic domination between the adjacent labels u, u+1 of the random friend
// tree, at the birth size N = u+1 of the later label (see COUPLING.md for the order ">=" and the coupling lemma).
//
// 1. The signed measure d = q_u - q_(u+1) on rooted shapes of size N is computed exactly (same code as
//    lead-signed_orbit_dp.cpp: unrooted law of T_n, newborn rooted at the new leaf, one more step for the older label).
// 2. Shapes with d > 0 whose root has degree >= 2 have a leaf child (the later label was attached to the older one);
//    they are >= every leaf-rooted shape, and every shape with d < 0 is leaf-rooted.  So the coupling exists as soon as
//    ALL the positive mass on leaf-rooted shapes can be transported, along the order, into the negative mass.
// 3. Flow network on leaf-rooted shapes of sizes <= N: an arc s -> t when t is s with one whole child branch deleted
//    at a parent of even depth (s >= t), and an arc t' -> t when t' is t with one whole child branch deleted at a
//    parent of odd depth (t' >= t).  The order is transitive, so every path joins comparable shapes.  Sources: shapes
//    of size N with d > 0 (capacity d); sinks: shapes of size N with d < 0 (capacity -d); inner arcs unbounded.
//    Exact integer max-flow (Dinic).  PASS iff the flow equals the total leaf-rooted positive mass.
//
// Usage: transport_sparse ULO UHI        Build: g++ -O2 -std=c++17 transport_sparse.cpp -lgmpxx -lgmp -o transport_sparse
#include <bits/stdc++.h>
#include <gmpxx.h>
using namespace std;
typedef uint64_t u64; typedef unsigned __int128 u128;
#ifndef SDG_LIMBS
#define SDG_LIMBS 6
#endif
#ifndef SDG_INIT
#define SDG_INIT 75
#endif
static const int MAXN = 26, LIMBS = SDG_LIMBS, K = 40;

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
    void init(size_t want) { cap = 1024; while (cap * SDG_INIT < want * 100) cap <<= 1; keys.assign(cap, 0); vals.assign(cap, bzero()); cnt = 0; }
    static inline u64 h(u64 x) { x ^= x >> 33; x *= 0xff51afd7ed558ccdULL; x ^= x >> 33; x *= 0xc4ceb9fe1a85ec53ULL; x ^= x >> 33; return x; }
    void add(u64 k, const Big& v) {
        if (shards > 1 && (int)((h(k) >> 40) % (u64)shards) != shard) return;
        size_t i = h(k) & (cap - 1);
        while (true) {
            if (keys[i] == k) { badd(vals[i], v); return; }
            if (keys[i] == 0) { keys[i] = k; vals[i] = v; if (++cnt * 10 > cap * 9) grow(); return; }
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

static void rawtest(const Table &a,int n,const mpz_class &den) {
    u64 J=lcm_upto(n-1), W=lcm_upto(n);
    // Consecutive denominators d,d+1 are coprime and both divide W.
    mpz_class bigW=(unsigned long)W;bigW*=(unsigned long)J;
    // Aggregates require more bits than a prior law entry. Batch one-limb-rate
    // products in a private 512-bit accumulator even when law entries use 7 limbs.
    struct Wide {
        array<u64,8> m{};
        void add(const Big &a,u128 x) {
            auto word=[&](u64 y,int k) {u128 carry=0;
                for(int i=0;i<LIMBS;i++) {carry+=(u128)a.l[i]*y+m[i+k];m[i+k]=carry;carry>>=64;}
                for(int i=LIMBS+k;i<8 && carry;i++) {carry+=m[i];m[i]=carry;carry>>=64;}assert(!carry);
            };
            word((u64)x,0);if(x>>64){assert(LIMBS<8);word(x>>64,1);}
        }
        mpz_class mp() const {mpz_class z;mpz_import(z.get_mpz_t(),8,-1,8,0,0,m.data());return z;}
    };
    if(den*bigW*n*n >= (mpz_class(1) << 511)) {fprintf(stderr,"raw bits exceeded\n");exit(1);}
    vector<Wide> P(n),Q(n);Wide U;
    for(size_t k=0;k<a.cap;k++) if(a.keys[k]) {
        Tree t;udecode(a.keys[k],n,t);
        u128 rho[MAXN];int l[MAXN]={};
        for(int v=0;v<n;v++) {
            u64 r=0;int d=t.deg[v];
            for(int j=0;j<d;j++) {int w=t.adj[v][j];r+=J/t.deg[w];l[v]+=t.deg[w]==1;}
            assert(W%(d*(d+1))==0);
            rho[v]=(u128)(W/(d*(d+1)))*r;
            U.add(a.vals[k],(u128)(W/(d+1))*r);
        }
        for(int v=0;v<n;v++) for(int j=0;j<t.deg[v];j++) {
            int w=t.adj[v][j];if(rho[w]<=rho[v])continue;
            u128 x=rho[w]-rho[v];
            if(t.deg[v]>1)P[t.deg[v]].add(a.vals[k],x);
            if(t.deg[w]>1)Q[l[w]-(t.deg[v]==1)].add(a.vals[k],x);
        }
    }
    mpz_class slack=0,p=0,q=0;int pivot=-1;
    for(int j=n-1;j>=2;j--) {
        mpz_class pj=P[j].mp();p+=pj;q+=Q[j].mp();mpz_class g=q-p;
        if(pj!=0 && (pivot<0 || g<slack)) {slack=g;pivot=j;}
    }
    // Retain the legacy J^2 W output scaling for identical replay logs.
    p*=J;slack*=J;mpz_class uu=U.mp()*J;
    printf("RAW N=%d supply=%s U=%s pivot=%d slack=%s slack/U=%.12g\n",n+2,p.get_str().c_str(),uu.get_str().c_str(),pivot,slack.get_str().c_str(),slack.get_d()/uu.get_d());
    fflush(stdout);
}

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


// ---------------- transport ----------------
static inline bool bless(const Big& a, const Big& b) { for (int i = LIMBS - 1; i >= 0; i--) if (a.l[i] != b.l[i]) return a.l[i] < b.l[i]; return false; }   // nonnegative values
static inline void bsub(Big& a, const Big& b) { Big nb = bnegate(b); badd(a, nb); }
static inline int size_of_key(u64 k) { return (64 - __builtin_clzll(k)) / 2; }
static int skipv;
static Code canon_skip(const Tree& t, int v, int p) {
    u64 kids[MAXN]; int k = 0;
    for (int i = 0; i < t.deg[v]; i++) { int w = t.adj[v][i]; if (w == p || w == skipv) continue; Code c = canon_skip(t, w, v); kids[k++] = ((u64)c.len << 56) | c.bits; }
    sort(kids, kids + k);
    u64 bits = 1; int len = 1;
    for (int i = 0; i < k; i++) { int cl = (int)(kids[i] >> 56); bits = (bits << cl) | (kids[i] & ((1ULL << 56) - 1)); len += cl; }
    return { bits << 1, len + 1 };
}
struct Net {
    unordered_map<u64, int> id; vector<u64> key; vector<char> seenP, seenM;
    vector<int> tail, head;                                   // inner arcs tail -> head (tail >= head in the order)
    int node(u64 k) { auto it = id.find(k); if (it != id.end()) return it->second; int i = (int)key.size(); id.emplace(k, i); key.push_back(k); seenP.push_back(0); seenM.push_back(0); return i; }
    // all shapes obtained from the shape of node i by deleting one whole child branch at a parent of the given parity
    void moves(int i, int parity, vector<u64>& out) {
        out.clear(); u64 k = key[i]; int n = size_of_key(k); Tree t; decode(k, n, t);
        int par[MAXN], dep[MAXN]; par[0] = -1; dep[0] = 0;       // decode numbers the vertices in preorder: parents come first
        for (int v = 0; v < n; v++) for (int j = 0; j < t.deg[v]; j++) { int w = t.adj[v][j]; if (w != par[v]) { par[w] = v; dep[w] = dep[v] + 1; } }
        for (int v = 1; v < n; v++) {
            if ((dep[par[v]] & 1) != parity) continue;
            if (par[v] == 0) continue;                         // deleting the only branch of the leaf root leaves a single vertex: useless
            skipv = v; out.push_back(canon_skip(t, 0, -1).bits);
        }
        sort(out.begin(), out.end()); out.erase(unique(out.begin(), out.end()), out.end());
    }
};
int main(int argc, char** argv) {
    if (argc < 3) { fprintf(stderr, "usage: %s ULO UHI\n", argv[0]); return 2; }
    int ulo = atoi(argv[1]), uhi = atoi(argv[2]);
    vector<mpz_class> Den(MAXN + 1); Den[2] = 1; for (int n = 2; n < MAXN; n++) Den[n + 1] = Den[n] * n * mpz_class((unsigned long)lcm_upto(n - 1));
    static const size_t R[] = {0, 1, 1, 2, 4, 9, 20, 48, 115, 286, 719, 1842, 4766, 12486, 32973, 87811, 235381, 634847, 1721159, 4688676, 12826228, 35221832, 97055181, 268282855, 743724984};
    static const size_t FR[] = {0, 1, 1, 1, 2, 3, 6, 11, 23, 47, 106, 235, 551, 1301, 3159, 7741, 19320, 48629, 123867, 317955, 823065, 2144505, 5623756, 14828074, 39299897};
    int top = max(uhi, 3); vector<Table> U(top + 2); int have = 2;
    U[2].init(4); { Big one = bzero(); one.l[0] = 1; U[2].add(0b110100ULL, one); }
    bool rawmode=argc>=4 && string(argv[3])=="raw";
    auto need = [&](int n) { while (have < n) { U[have + 1].init(FR[have + 1]); propagate(U[have], have, U[have + 1], 2, false, true);
            if (U[have + 1].total() != Den[have + 1] || U[have + 1].cnt != FR[have + 1]) { fprintf(stderr, "mass or census check failed at size %d\n", have + 1); exit(1); }
            if(rawmode) U[have].clear();
            have++; } };
    for (int u = max(ulo, 2); u <= uhi; u++) {
        if(rawmode) {
            if(u>=3) {need(u-1);rawtest(U[u-1],u-1,Den[u-1]);}
            continue;
        }
        int N = u + 1; Table A; A.init(R[u] + R[u - 1]);
        need(u); propagate(U[u], u, A, 1, true, true);                                              // minus newborn u+1
        if (u == 2) { Big one = bzero(); one.l[0] = 1; A.add(0b110100ULL, one); A.add(0b111000ULL, one); }
        else { Table Qu; Qu.init(R[u - 1]); propagate(U[u - 1], u - 1, Qu, 1, false, true); if (Qu.total() != Den[u]) { fprintf(stderr, "newborn mass check failed\n"); return 1; } propagate(Qu, u, A, 0); }
        if (A.total() != 0) { fprintf(stderr, "difference mass not zero\n"); return 1; }
        if (Den[N] >= (mpz_class(1) << (64 * LIMBS - 2))) { fprintf(stderr, "integers too short for size %d\n", N); return 1; }
        // sources, sinks, universal mass
        bool classonly = argc >= 4;
        vector<mpz_class> supply(N), demand(N);
        Net net; vector<Big> capS, capT; vector<int> srcNode, snkNode; mpz_class universal = 0, want = 0, sinkmass = 0; Tree t;
        size_t cpos=0,cneg=0;
        for (size_t i = 0; i < A.cap; i++) {
            if (!A.keys[i] || biszero(A.vals[i])) continue;
            decode(A.keys[i], N, t); const Big& m = A.vals[i];
            if (t.deg[0] >= 2) {                                   // root with a leaf child: universal
                bool leafchild = false; for (int j = 0; j < t.deg[0]; j++) leafchild |= (t.deg[t.adj[0][j]] == 1);
                if (bneg(m) || !leafchild) { fprintf(stderr, "unexpected non-leaf-rooted shape\n"); return 1; }
                universal += tompz(m); continue;
            }
            int p = t.adj[0][0], leaves = 0;
            for (int j=0;j<t.deg[p];j++) {int w=t.adj[p][j]; leaves += w!=0 && t.deg[w]==1;}
            if (!bneg(m)) { supply[t.deg[p]-1] += tompz(m); cpos++; }
            else { demand[leaves] += tompz(bnegate(m)); cneg++; }
            if (classonly) continue;
            int v = net.node(A.keys[i]);
            if (!bneg(m)) { srcNode.push_back(v); capS.push_back(m); want += tompz(m); }
            else { snkNode.push_back(v); capT.push_back(bnegate(m)); sinkmass += tompz(bnegate(m)); }
        }
        A.clear();
        mpz_class ps=0, ns=0, slack=0, deficit=0; int pivot=-1;
        for(int j=N-1;j>=1;j--) {
            ps+=supply[j]; ns+=demand[j];
            mpz_class gap=ns-ps;
            if(supply[j]>0 && (pivot<0 || gap<slack)) {slack=gap; pivot=j;}
            if(supply[j]>0 && gap<0) deficit=max(deficit,mpz_class(-gap));
        }
        assert(ps+universal==ns+demand[0]);
        mpz_class classflow=ps-deficit;
        printf("LEAVES N=%d +%zu -%zu supply=%s transported=%s U=%s pivot=%d slack=%s slack/U=%.12g\n",N,cpos,cneg,ps.get_str().c_str(),classflow.get_str().c_str(),universal.get_str().c_str(),pivot,slack.get_str().c_str(),slack.get_d()/universal.get_d());
        fflush(stdout);
        if(classonly) continue;
        if (want + universal != sinkmass) { fprintf(stderr, "mass balance failed\n"); return 1; }
        size_t nsrc = srcNode.size(), nsnk = snkNode.size();
        // closures: downwards from the sources by even-parent deletions, downwards from the sinks by odd-parent deletions
        vector<u64> mv; vector<int> queue;
        for (int v : srcNode) { net.seenP[v] = 1; queue.push_back(v); }
        for (size_t h = 0; h < queue.size(); h++) { int x = queue[h]; net.moves(x, 0, mv);
            for (u64 k : mv) { int y = net.node(k); net.tail.push_back(x); net.head.push_back(y); if (!net.seenP[y]) { net.seenP[y] = 1; queue.push_back(y); } } }
        queue.clear(); for (int v : snkNode) { net.seenM[v] = 1; queue.push_back(v); }
        for (size_t h = 0; h < queue.size(); h++) { int x = queue[h]; net.moves(x, 1, mv);
            for (u64 k : mv) { int y = net.node(k); net.tail.push_back(y); net.head.push_back(x); if (!net.seenM[y]) { net.seenM[y] = 1; queue.push_back(y); } } }
        int V = (int)net.key.size(); size_t E = net.tail.size();
        // adjacency (CSR): arc code 2e for tail -> head (unbounded), 2e+1 for head -> tail (capacity = current flow)
        vector<int> deg(V + 1, 0); for (size_t e = 0; e < E; e++) { deg[net.tail[e]]++; deg[net.head[e]]++; }
        vector<size_t> st(V + 1, 0); for (int v = 0; v < V; v++) st[v + 1] = st[v] + deg[v];
        vector<unsigned> adj(2 * E); { vector<size_t> pos(st.begin(), st.end() - 1); for (size_t e = 0; e < E; e++) { adj[pos[net.tail[e]]++] = (unsigned)(2 * e); adj[pos[net.head[e]]++] = (unsigned)(2 * e + 1); } }
        vector<Big> flow(E, bzero()); vector<int> sinkIdx(V, -1), level(V); for (size_t j = 0; j < nsnk; j++) sinkIdx[snkNode[j]] = (int)j;
        vector<size_t> it(V); Big INF = bzero(); INF.l[LIMBS - 1] = 1ULL << 62;
        function<Big(int, Big)> dfs = [&](int x, Big lim) -> Big {
            if (sinkIdx[x] >= 0 && !biszero(capT[sinkIdx[x]])) { Big g = bless(lim, capT[sinkIdx[x]]) ? lim : capT[sinkIdx[x]]; bsub(capT[sinkIdx[x]], g); return g; }
            for (; it[x] < st[x + 1]; it[x]++) {
                unsigned a = adj[it[x]]; size_t e = a >> 1; int y = (a & 1) ? net.tail[e] : net.head[e];
                if (level[y] != level[x] + 1) continue;
                Big cap = (a & 1) ? flow[e] : INF; if (biszero(cap)) continue;
                Big g = dfs(y, bless(lim, cap) ? lim : cap);
                if (!biszero(g)) { if (a & 1) bsub(flow[e], g); else badd(flow[e], g); return g; }
            }
            return bzero();
        };
        int phases = 0;
        while (true) {
            fill(level.begin(), level.end(), -1); vector<int> bq; bool reach = false;
            for (size_t i = 0; i < nsrc; i++) if (!biszero(capS[i])) { level[srcNode[i]] = 0; bq.push_back(srcNode[i]); }
            for (size_t h = 0; h < bq.size(); h++) { int x = bq[h]; if (sinkIdx[x] >= 0 && !biszero(capT[sinkIdx[x]])) reach = true;
                for (size_t p = st[x]; p < st[x + 1]; p++) { unsigned a = adj[p]; size_t e = a >> 1; int y = (a & 1) ? net.tail[e] : net.head[e];
                    if (level[y] >= 0) continue; if ((a & 1) && biszero(flow[e])) continue; level[y] = level[x] + 1; bq.push_back(y); } }
            if (!reach) break;
            phases++; for (int v = 0; v < V; v++) it[v] = st[v];
            for (size_t i = 0; i < nsrc; i++) while (!biszero(capS[i])) { Big g = dfs(srcNode[i], capS[i]); if (biszero(g)) break; bsub(capS[i], g); }
        }
        mpz_class left = 0; for (size_t i = 0; i < nsrc; i++) left += tompz(capS[i]);
        mpz_class room = 0; for (size_t j = 0; j < nsnk; j++) room += tompz(capT[j]);
        mpz_class flowed = want - left; bool pass = (left == 0) && (room == universal);
        printf("pair (%d,%d) size %d: leaf-rooted shapes + %zu, - %zu; network %d shapes, %zu arcs, %d phases; leaf-rooted positive mass %s, transported %s, universal %s, denominator bits %zu: %s\n",
               u, u + 1, N, nsrc, nsnk, V, E, phases, want.get_str().c_str(), flowed.get_str().c_str(), universal.get_str().c_str(), mpz_sizeinbase(Den[N].get_mpz_t(), 2),
               pass ? "PASS (full stochastic domination certified)" : "NOT CERTIFIED");
        fflush(stdout);
    }
    return 0;
}
