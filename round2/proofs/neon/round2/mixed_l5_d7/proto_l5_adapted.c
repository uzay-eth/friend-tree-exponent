/* Included after the B3-FullH helpers in diag_slack.c. Pure-inner L5 DP.
 * Tables are indexed by rooted shape, apex h, and multisets of equal subtrees.
 * Coarse modes also support one or two large endpoint buckets. */
#include <assert.h>

typedef struct {
    int n, a, child[4], ns[2], state[2][64], ix[2][64*64];
    size_t len[2];
    double *q[2];
} L5Shape;

typedef struct { int n, p[5], g[5]; } L5Tree;
typedef struct { uint64_t key; double val; } CoarseCache;
typedef struct {
    CoarseCache *cache;
    int nf[3][2], *fiber[3]; /* apex (both h), apex h=1, descendant h=1; inner then tail */
    int d, count, path4, top4, bot4, star4;
    size_t entries, visited;
    const B3FLayout *c;
    L5Shape s[17];  /* rooted shapes of size 1..5: 1+1+2+4+9 */
    const double *p4, *yt, *yb, *s3;
    double lam;
    int mode;
} L5;

static size_t l5_choose(size_t n, int k) {
    if (n < (size_t)k) return 0;
    size_t r = 1;
    for (int j = 1; j <= k; ++j) r = r * (n - j + 1) / j;
    return r;
}

static int l5_find(const L5 *x, int a, const int *child) {
    for (int i = 0; i < x->count; ++i)
        if (x->s[i].a == a && !memcmp(x->s[i].child, child, (size_t)a*sizeof(int))) return i;
    return -1;
}

static void l5_shapes(L5 *x, int n, int a, int first, int left, int kids[4]) {
    if (!left) {
        assert(x->count < 17);
        L5Shape *s = &x->s[x->count++];
        s->n = n; s->a = a; memcpy(s->child, kids, (size_t)a*sizeof(int));
        return;
    }
    for (int i = first; i < x->count && x->s[i].n <= left; ++i) {
        kids[a] = i;
        l5_shapes(x, n, a+1, i, left-x->s[i].n, kids);
    }
}

/* -1/-2 = BIG_INNER (or unsplit BIG), -3/-4 = TAIL; parity encodes index h. */
static inline int coarse_id(int h, int b) { return -1-h-2*b; }
static inline int coarse_bucket(int g) { return (-g-1)/2; }
static inline int l5_h(const L5 *x, int g) { return g < 0 ? (-g-1)%2 : x->c->H[g]; }
static inline int l5_degree(const L5Tree *t, int u) {
    int d = !!u;
    for (int i = u+1; i < t->n; ++i) d += t->p[i] == u;
    return d;
}
static int coarse_selected(const L5Tree *t, int u) {
    return coarse_mask & (l5_degree(t, u) > 1 ? 1 : 2);
}
static int coarse_project4(const L5 *x, L5Tree *t) {
    for (int i=0; i<4; ++i) {
        int g = t->g[i];
        if (x->c->kind[g] == 0 && x->c->D[g] < coarse) continue;
        if (l5_degree(t, i) != 1) return 0;
        /* INF is shared by both apex h's, and q5 could depend on h. Its
           own bud is excluded; at other slots evaluate using a common BIG
           apex bucket (see l5_h override in the shape index below). */
        t->g[i] = coarse_id(i != 0, coarse_buckets == 2 && x->c->kind[g] != 0);
    }
    return 1;
}

static L5 *l5_init(const B3FLayout *c, int d) {
    L5 *x = b3f_xcalloc(1, sizeof(*x));
    x->c = c; x->d = d;
    int kids[4];
    for (int n = 1; n <= 5; ++n) l5_shapes(x, n, 0, 0, n-1, kids);
    assert(x->count == 17 && c->INNER_N <= 64*64 && d <= c->D0);
    for (int i = 0; i < x->count; ++i) {
        L5Shape *s = &x->s[i];
        for (int h = 0; h < 2; ++h) {
            for (int g = 0; g < c->INNER_N; ++g) s->ix[h][g] = -1;
            for (int g = 0; g < c->INNER_N; ++g) if (c->D[g] < d && c->H[g] == h && c->K[g] >= s->a) {
                assert(s->ns[h] < 64);
                s->ix[h][g] = s->ns[h]; s->state[h][s->ns[h]++] = g;
            }
            /* A BIG apex uses h=0 in the index; h never affects a BIG
               transition or bound. A nonapex is necessarily h=1. */
            if (coarse && ((h==1 && s->a==0) || (h==0 && s->a==1)))
                for (int b=0; b<coarse_buckets; ++b) {
                    assert(s->ns[h] < 64);
                    s->state[h][s->ns[h]++] = coarse_id(h,b);
                }
            size_t z = s->ns[h];
            for (int j = 0, end; j < s->a; j = end) {
                for (end = j+1; end < s->a && s->child[end] == s->child[j]; ++end) {}
                z *= l5_choose(x->s[s->child[j]].len[1] + end-j-1, end-j);
            }
            s->len[h] = z;
            if (s->n == 5) { s->q[h] = b3f_xmalloc(z*sizeof(double)); x->entries += z; }
        }
    }
    if (coarse) {
        x->cache = b3f_xcalloc(1u<<19, sizeof(*x->cache));
        for (int r=0; r<3; ++r) {
            x->fiber[r] = b3f_xmalloc(c->NG * sizeof(int));
            for (int b=0, at=0; b<coarse_buckets; ++b) for (int g=0; g<c->NG; ++g) {
                int k=c->kind[g];
                int valid = r==2 ? c->w_to_wi[g]>=0 : c->p_to_pi[g]>=0;
                if (coarse_buckets==2 && b!=(k!=0)) continue;
                if (valid && (k==1 || ((r==0 || c->H[g]==1) && (k || c->D[g]>=coarse)))) {
                    x->fiber[r][at++] = g;
                    ++x->nf[r][b];
                }
            }
        }
    }
    if (coarse) {
        size_t bytes = x->entries*sizeof(double) + sizeof(*x)
                     + (1u<<19)*sizeof(*x->cache) + 3*c->NG*sizeof(int);
        printf("COARSE allocated_payload_bytes=%zu coarse=%d mask=%d buckets=%d apex_fibers=%d/%d child_fibers=%d/%d\n",
               bytes, coarse, coarse_mask, coarse_buckets,
               x->nf[0][0], x->nf[0][1], x->nf[2][0], x->nf[2][1]);
    }
    printf("L5 D5=%d entries=%zu bytes=%zu\n", d, x->entries, x->entries*sizeof(double)); fflush(stdout);
    return x;
}

static void l5_free(L5 *x) {
    if (!x) return;
    for (int i = 0; i < x->count; ++i) for (int h = 0; h < 2; ++h) free(x->s[i].q[h]);
    free(x->cache);
    for (int r=0; r<3; ++r) free(x->fiber[r]);
    free(x);
}

static int l5_encode(const L5 *x, const L5Tree *t, int u, size_t *row) {
    int a = 0, kid[4]; size_t key[4];
    for (int v = u+1; v < t->n; ++v) if (t->p[v] == u) {
        size_t r; int id = l5_encode(x, t, v, &r);
        int j = a++;
        while (j && (kid[j-1] > id || (kid[j-1] == id && key[j-1] > r))) {
            kid[j] = kid[j-1]; key[j] = key[j-1]; --j;
        }
        kid[j] = id; key[j] = r;
    }
    int id = l5_find(x, a, kid), g = t->g[u], h = l5_h(x,g);
    assert(id >= 0);
    int idx = g < 0 ? x->s[id].ns[h]-coarse_buckets+coarse_bucket(g) : x->s[id].ix[h][g];
    assert(idx >= 0 && x->s[id].state[h][idx]==g);
    size_t r = idx;
    for (int j = 0, end; j < a; j = end) {
        size_t z = 0;
        for (end = j; end < a && kid[end] == kid[j]; ++end)
            z += l5_choose(key[end] + end-j, end-j+1);
        r = r*l5_choose(x->s[kid[j]].len[1]+end-j-1, end-j)+z;
    }
    assert(r < x->s[id].len[h]);
    *row = r;
    return id;
}

static void l5_decode(const L5 *x, L5Tree *t, int id, int h, int p, size_t row) {
    const L5Shape *s = &x->s[id];
    size_t keys[4];
    for (int j = s->a, first; j > 0; j = first) {
        int k = s->child[j-1];
        for (first = j-1; first && s->child[first-1] == k; --first) {}
        size_t max = x->s[k].len[1] + j-first-1;
        size_t radix = l5_choose(max, j-first), m = row % radix;
        row /= radix;
        for (int i = j-first; i > 0; --i) {
            size_t lo = i-1, hi = max;
            while (lo+1 < hi) {
                size_t mid = lo+(hi-lo)/2;
                if (l5_choose(mid, i) <= m) lo = mid; else hi = mid;
            }
            m -= l5_choose(lo, i); max = lo;
            keys[first+i-1] = lo-(i-1);
        }
        assert(!m);
    }
    int u = t->n++;
    assert(row < (size_t)s->ns[h]);
    t->p[u] = p; t->g[u] = s->state[h][row];
    for (int j = 0; j < s->a; ++j) l5_decode(x, t, s->child[j], 1, u, keys[j]);
}

static double l5_read(L5 *x, L5Tree *t);

static double l5_q4(const L5 *x, const L5Tree *t) {
    const B3FLayout *c = x->c;
    int n[4] = {0}, child[4][3], g[4];
    for (int i = 0; i < 4; ++i) g[i] = t->g[i];
    for (int i = 1; i < 4; ++i) { int p = t->p[i]; child[p][n[p]++] = i; }
    if (n[0] == 3) {
        int u = c->s3_lut[c->H[g[0]]][c->K[g[0]]][c->L[g[0]]];
        return u < 0 ? 0 : x->s3[b3f_s3_off(c, u, g[1], g[2], g[3])];
    }
    int b = child[0][0];
    if (n[0] == 2) {
        int a = child[0][1];
        if (n[a]) { int z = b; b = a; a = z; }
        assert(n[b] == 1 && !n[a]);
        int u = c->ybu_lut[c->H[g[0]]][c->K[g[0]]][c->L[g[0]]];
        int w = c->ybw_lut[c->K[g[b]]][c->L[g[b]]];
        return u < 0 || w < 0 ? 0 : x->yb[b3f_yb_off(c, u, w, g[a], g[child[b][0]])];
    }
    if (n[b] == 2) {
        int u = c->yt_lut[c->K[g[b]]][c->L[g[b]]];
        return u < 0 ? 0 : x->yt[b3f_yt_off(c, g[0], u, g[child[b][0]], g[child[b][1]])];
    }
    int mid = child[b][0];
    int b4 = c->p4_lut[c->K[g[b]]][c->L[g[b]]];
    int c4 = c->p4_lut[c->K[g[mid]]][c->L[g[mid]]];
    assert(n[b] == 1 && n[mid] == 1);
    if (b4 < 0 || c4 < 0) return 0;
    return x->p4[b3f_p4_off(c, g[0], b4, c4, g[child[mid][0]])];
}

#include "coarse_source.c"

static double l5_read(L5 *x, L5Tree *t) {
    const B3FLayout *c = x->c;
    size_t row;
    int id = l5_encode(x, t, 0, &row);
    double *out = &x->s[id].q[l5_h(x,t->g[0])][row];
    if (!isnan(*out)) return *out;
    int neighbors[5] = {0};
    for (int i = 1; i < 5; ++i) ++neighbors[i], ++neighbors[t->p[i]];
    double num = 0, den = x->lam, up[5], slack[5];
    for (int z = 0; z < 5; ++z) if (neighbors[z] == 1) {
        int u = z ? t->p[z] : 1;
        if (t->g[z]<0) continue;
        double w = 1.0/c->D[t->g[z]] - 1.0/x->d;
        L5Tree m = {0};
        for (int v = 0; v < 5; ++v) if (v != z) {
            int p = t->p[v];
            m.g[m.n] = t->g[v]; m.p[m.n++] = p == z ? -1 : p - (p > z);
        }
        int v = u-(u>z), g = m.g[v];
        if (coarse) {
            if (coarse_selected(&m, v)) num += w * coarse_source(x, &m, v);
            continue;
        }
        double before = l5_q4(x, &m);
        m.g[v] = b3f_next_up(c, g);
        /* At D5-1, q4 can remain active after the gain. */
        num += w * (l5_q4(x, &m) - before);
    }
    for (int i = 0; i < 5; ++i) {
        int g = t->g[i];
        if (g<0) {
            up[i]=slack[i]=0;
            if (coarse_buckets==2 && coarse_bucket(g)==0) {
                /* Only D0-1 can cross; its rate L+rho <= D0-1.
                   No crossing otherwise. Promotion does not change bucket. */
                t->g[i] = coarse_id(l5_h(x,g),1);
                up[i] = l5_read(x,t);
                t->g[i] = g;
                slack[i] = c->D0-1;
            }
            continue;
        }
        double r = c->L[g], extra = 0;
        for (int j=0; j<5; ++j) if ((i && t->p[i]==j) || t->p[j]==i) {
            if (t->g[j]>=0) r+=1.0/c->D[t->g[j]];
            else extra+=1.0/(coarse_bucket(t->g[j]) ? c->D0 : coarse);
        }
        if (c->D[g]+1 < x->d) {
            t->g[i] = b3f_next_up(c, g); up[i] = l5_read(x, t); t->g[i] = g;
        } else if (coarse && neighbors[i]==1) {
            t->g[i] = coarse_id(i != 0,0); up[i] = l5_read(x, t); t->g[i]=g;
        } else up[i] = 0;
        num += r*up[i]; den += r;
        slack[i] = extra + 0.5*(c->K[g] + c->H[g] - neighbors[i]);
        assert(slack[i] >= 0);
        if (c->L[g]) {
            double b = (double)c->L[g]/c->D[g];
            t->g[i] = b3f_next_bud(c, g);
            num += b*l5_read(x, t); den += b;
            t->g[i] = g;
        }
    }
    double v = num/den;
    for (int tries = 0, previous = -1; tries < 32; ++tries) {
        int bits = 0; double n = num, d = den;
        for (int i = 0; i < 5; ++i) if (slack[i] && (x->mode == 1 ? up[i] > v : up[i] < v)) {
            bits |= 1<<i; n += slack[i]*up[i]; d += slack[i];
        }
        v = n/d;
        if (bits == previous) { ++x->visited; return *out = v; }
        previous = bits;
    }
    abort();
}

static void l5_sweep(L5 *x, int mode, double lam, const double *p4, const double *yt, const double *yb, const double *s3) {
    x->mode = mode; x->lam = lam; x->p4 = p4; x->yt = yt; x->yb = yb; x->s3 = s3;
    x->visited = 0;
    if (coarse) memset(x->cache, 0, (1u<<19)*sizeof(*x->cache));
    for (int i = 0; i < x->count; ++i) if (x->s[i].n == 5) for (int h = 0; h < 2; ++h)
        for (size_t j = 0; j < x->s[i].len[h]; ++j) x->s[i].q[h][j] = NAN;
    for (int i = 0; i < x->count; ++i) if (x->s[i].n == 5) for (int h = 0; h < 2; ++h)
        for (size_t j = 0; j < x->s[i].len[h]; ++j) if (isnan(x->s[i].q[h][j])) {
            L5Tree t = {0}; size_t row;
            l5_decode(x, &t, i, h, -1, j);
            assert(l5_encode(x, &t, 0, &row) == i && row == j);
            l5_read(x, &t);
            assert(!isnan(x->s[i].q[h][j]));
        }
    assert(x->visited == x->entries);
}

static double l5_birth(const L5 *x, L5Tree t, int at, int newborn) {
    const B3FLayout *c = x->c;
    if (!x || !x->visited || !c->L[t.g[at]]) return 0;
    t.g[at] = b3f_next_bud(c, t.g[at]);
    t.n = 5; t.p[4] = at; t.g[4] = newborn;
    size_t row;
    int id = l5_encode(x, &t, 0, &row);
    return x->s[id].q[l5_h(x,t.g[0])][row];
}

static void l5_scale(L5 *x, double rate) {
    for (int i = 0; i < x->count; ++i) if (x->s[i].n == 5) for (int h = 0; h < 2; ++h)
        for (size_t j = 0; j < x->s[i].len[h]; ++j) x->s[i].q[h][j] *= rate;
}

static void l5_selfcheck(L5 *x) {
    const B3FLayout *c = x->c;
    /* Constant zero L4 yields exactly zero L5, independent of all event rates. */
    double *p4 = b3f_xcalloc(c->p4_len, sizeof(double)), *yt = b3f_xcalloc(c->yt_len, sizeof(double));
    double *yb = b3f_xcalloc(c->yb_len, sizeof(double)), *s3 = b3f_xcalloc(c->s3_len, sizeof(double));
    l5_sweep(x, 1, 0.5657237, p4, yt, yb, s3);
    for (int i = 0; i < x->count; ++i) if (x->s[i].n == 5) for (int h = 0; h < 2; ++h)
        for (size_t j = 0; j < x->s[i].len[h]; ++j) assert(x->s[i].q[h][j] == 0);
    /* A varying table must create a nonzero response; check both directions. */
    for (size_t j = 0; j < c->p4_len; ++j) p4[j] = 0.001*(j%17);
    for (int mode = 1; mode <= 2; ++mode) {
        l5_sweep(x, mode, 0.5657237, p4, yt, yb, s3);
        size_t nz = 0;
        for (int i = 0; i < x->count; ++i) if (x->s[i].n == 5) for (int h = 0; h < 2; ++h)
            for (size_t j = 0; j < x->s[i].len[h]; ++j) {
                assert(isfinite(x->s[i].q[h][j]));
                nz += x->s[i].q[h][j] != 0;
            }
        assert(nz);
        printf("L5 check mode=%d rows=%zu nonzero=%zu\n", mode, x->visited, nz);
    }
    free(p4); free(yt); free(yb); free(s3);
}
