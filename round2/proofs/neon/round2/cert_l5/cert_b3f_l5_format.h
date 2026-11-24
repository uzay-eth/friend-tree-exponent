#ifndef CERT_B3F_L5_FORMAT_H
#define CERT_B3F_L5_FORMAT_H

#include <assert.h>
#include <limits.h>
#include <stdio.h>
#include "../cert/cert_b3f_format.h"

_Static_assert(CHAR_BIT == 8 && sizeof(long) == 8 && sizeof(void *) == 8 &&
               sizeof(size_t) == 8 && sizeof(int) == 4, "v2 requires LP64");

typedef struct { uint32_t dc, buckets, mask, reserved; uint64_t rows; } L5Header;
enum { L5_EXT_BYTES = 24, L5_SHAPES = 17, L5_CACHE = 1 << 19 };

static inline int l5_header_io(FILE *fp, L5Header *h, int writing) {
    if (!b3f_is_little_endian()) return 0;
#define FIELD(f) (writing ? fwrite(&h->f, sizeof h->f, 1, fp) : fread(&h->f, sizeof h->f, 1, fp))
    return FIELD(dc) == 1 && FIELD(buckets) == 1 && FIELD(mask) == 1 &&
           FIELD(reserved) == 1 && FIELD(rows) == 1;
#undef FIELD
}

static inline int b3f_l5_validate_header(const CertB3FHeader *h, const L5Header *e) {
    CertB3FHeader v1 = *h;
    v1.version = 1;
    int err = b3f_validate_header(&v1);
    if (err) return err;
    if (h->version != 2 || e->reserved || e->mask < 1 || e->mask > 3 ||
        e->buckets < 1 || e->buckets > 2 || e->dc < 4 || e->dc > CERT_B3F_MAX_D0) return -20;
    int d = (int)e->dc;
    if (d > h->D0 || (e->buckets == 2 && d == h->D0) ||
        d > h->DP4 || d > h->DYT || d > h->DYB || d > h->DS3 ||
        h->KP4 < d-1 || h->KYT < d-1 || h->KYB < d-1 || h->KS3 < d-1) return -21;
    return 0;
}

static inline size_t l5_add(size_t a, size_t b) {
    if (a > SIZE_MAX-b) return SIZE_MAX;
    return a+b;
}
static inline size_t l5_mul(size_t a, size_t b) {
    if (b && a > SIZE_MAX/b) return SIZE_MAX;
    return a*b;
}
static inline void *b3f_l5_xmalloc(size_t bytes) {
    void *p = bytes == SIZE_MAX ? NULL : malloc(bytes ? bytes : 1);
    if (!p) { fprintf(stderr,"Allocation failed (%zu bytes)\n",bytes); exit(1); }
    return p;
}
static inline void *b3f_l5_xcalloc(size_t n, size_t width) {
    size_t bytes = l5_mul(n,width);
    void *p = bytes == SIZE_MAX ? NULL : calloc(bytes ? bytes : 1,1);
    if (!p) { fprintf(stderr,"Allocation failed (%zu entries of %zu bytes)\n",n,width); exit(1); }
    return p;
}
/* The inherited v1 layout only allocates small, header-bounded index lists.
 * Use overflow-checked allocation for all v2 arrays and row workspaces. */
#define b3f_xmalloc b3f_l5_xmalloc
#define b3f_xcalloc b3f_l5_xcalloc
static inline void *l5_array(size_t n, size_t size) {
    return b3f_xmalloc(l5_mul(n,size));
}

typedef struct { int n, p[5], g[5]; } L5Tree;
typedef struct {
    int n, a, child[4], ns[2];
    int state[2][CERT_B3F_MAX_D*CERT_B3F_MAX_D];
    int ix[2][CERT_B3F_MAX_D*CERT_B3F_MAX_D];
    size_t len[2], off[2];
} L5Shape;
typedef struct {
    const B3FLayout *c;
    int dc, buckets, mask, count, unary[L5_SHAPES];
    size_t rows;
    L5Shape s[L5_SHAPES];
    int nf[2][2], *fiber[2][2]; /* apex (both physical h), child (h=1) */
} L5;

static inline size_t l5_choose(size_t n, int k) {
    if (n < (size_t)k) return 0;
    size_t r = 1;
    for (int i = 1; i <= k; ++i) {
        r = l5_mul(r, n-i+1);
        if (r == SIZE_MAX) return r;
        r /= i;
    }
    return r;
}
static inline int l5_big(int h, int b) { return -1-h-2*b; }
static inline int l5_bucket(int g) { return (-g-1)/2; }
static inline int l5_h(const L5 *x, int g) { return g < 0 ? (-g-1)%2 : x->c->H[g]; }
static inline int l5_degree(const L5Tree *t, int u) {
    int d = !!u;
    for (int i = u+1; i < t->n; ++i) d += t->p[i] == u;
    return d;
}
static inline int l5_charge(const L5Tree *t) {
    return l5_degree(t,0) > 1 ? 0 : 1;
}
static inline unsigned long l5_charge_cap(int branches, int dc) {
    return (unsigned long)branches*l5_choose(branches-1+4*(dc-2),3)/4;
}
/* Upper bound for embeddings of this rooted shape into a tree with at
 * most k children at its root and Dc-2 children at each supported
 * non-root interior. Identical subtree shapes use unordered choices. */
static inline unsigned long l5_shape_charge_cap(const L5 *x, int id, int k) {
    const L5Shape *s = &x->s[id];
    if (s->a > k) return 0;
    unsigned long result = 1;
    for (int j = 0, end; j < s->a; j = end) {
        for (end = j+1; end < s->a && s->child[end] == s->child[j]; ++end) {}
        result *= l5_choose(k,end-j); k -= end-j;
        for (int i = j; i < end; ++i) result *= l5_shape_charge_cap(x,s->child[i],x->dc-2);
    }
    return result;
}
static inline int l5_selected(const L5 *x, const L5Tree *t, int u) {
    return !!(x->mask & (l5_degree(t, u) > 1 ? 1 : 2));
}
static inline void l5_shapes(L5 *x, int n, int a, int first, int left, int kids[4]) {
    if (!left) {
        assert(x->count < L5_SHAPES);
        if (a == 1) x->unary[kids[0]] = x->count;
        L5Shape *s = &x->s[x->count++];
        s->n = n; s->a = a; memcpy(s->child, kids, (size_t)a*sizeof(int));
        return;
    }
    if (a >= 4) abort();
    for (int i = first; i < x->count && x->s[i].n <= left; ++i) {
        kids[a] = i;
        l5_shapes(x, n, a+1, i, left-x->s[i].n, kids);
    }
}
static inline L5 *l5_init(const B3FLayout *c, const L5Header *e) {
    L5 *x = b3f_xcalloc(1, sizeof *x);
    x->c = c; x->dc = (int)e->dc; x->buckets = (int)e->buckets; x->mask = (int)e->mask;
    int kids[4];
    for (int n = 1; n <= 5; ++n) l5_shapes(x, n, 0, 0, n-1, kids);
    for (int i = 0; i < x->count; ++i) for (int h = 0; h < 2; ++h) {
        L5Shape *s = &x->s[i];
        memset(s->ix[h], -1, sizeof s->ix[h]);
        for (int g = 0; g < c->INNER_N; ++g)
            if (c->D[g] < x->dc && c->H[g] == h && c->K[g] >= s->a) {
                s->ix[h][g] = s->ns[h]; s->state[h][s->ns[h]++] = g;
            }
        /* Coarse apex collapses both physical h values to h=0. */
        if (s->a == 1-h) for (int b = 0; b < x->buckets; ++b)
            s->state[h][s->ns[h]++] = l5_big(h,b);
        size_t z = s->ns[h];
        for (int j = 0, end; j < s->a; j = end) {
            for (end = j+1; end < s->a && s->child[end] == s->child[j]; ++end) {}
            z = l5_mul(z, l5_choose(l5_add(x->s[s->child[j]].len[1], end-j-1), end-j));
        }
        s->len[h] = z;
        if (s->n == 5) { s->off[h] = x->rows; x->rows = l5_add(x->rows, z); }
    }
    return x;
}

static inline size_t l5_expected_bytes(const L5 *x) {
    return l5_add(l5_add(b3f_expected_bytes(x->c), L5_EXT_BYTES), l5_mul(x->rows, sizeof(int64_t)));
}

static inline void l5_init_fibers(L5 *x) {
    const B3FLayout *c = x->c;
    for (int r = 0; r < 2; ++r) for (int b = 0; b < x->buckets; ++b) {
        int *f = x->fiber[r][b] = l5_array(c->NG, sizeof *f);
        for (int g = 0; g < c->NG; ++g) {
            int k = c->kind[g];
            if (x->buckets == 2 && b != (k != 0)) continue;
            if ((r ? c->w_to_wi[g] : c->p_to_pi[g]) >= 0 &&
                (k == 1 || ((!r || c->H[g] == 1) && (k || c->D[g] >= x->dc))))
                f[x->nf[r][b]++] = g;
        }
        assert(x->nf[r][b]);
    }
}
static inline void l5_free(L5 *x) {
    for (int r = 0; r < 2; ++r) for (int b = 0; b < 2; ++b) free(x->fiber[r][b]);
    free(x);
}
static inline int l5_shape_id(const L5 *x, int a, const int kid[4]) {
    if (a == 0) return 0;
    if (a == 1) { assert(x->unary[kid[0]] > 0); return x->unary[kid[0]]; }
    int id = a == 4 ? 8 : a == 3 ? 4 : 2;
    while (id < x->count && (x->s[id].a != a || memcmp(x->s[id].child,kid,(size_t)a*sizeof(int)))) ++id;
    assert(id < x->count);
    return id;
}
static inline int l5_encode(const L5 *x, const L5Tree *t, int u, size_t *row) {
    int a = 0, kid[4] = {0}; size_t key[4];
    for (int v = u+1; v < t->n; ++v) if (t->p[v] == u) {
        size_t r; int id = l5_encode(x, t, v, &r);
        int j = a++;
        while (j && (kid[j-1] > id || (kid[j-1] == id && key[j-1] > r))) {
            kid[j] = kid[j-1]; key[j] = key[j-1]; --j;
        }
        kid[j] = id; key[j] = r;
    }
    int id = l5_shape_id(x,a,kid);
    const L5Shape *s = &x->s[id];
    int g = t->g[u], h = l5_h(x,g);
    int idx = g < 0 ? s->ns[h]-x->buckets+l5_bucket(g) : s->ix[h][g];
    assert(idx >= 0 && s->state[h][idx] == g);
    size_t r = (size_t)idx;
    for (int j = 0, end; j < a; j = end) {
        size_t z = 0;
        for (end = j; end < a && kid[end] == kid[j]; ++end)
            z += l5_choose(key[end]+end-j, end-j+1);
        r = r*l5_choose(x->s[kid[j]].len[1]+end-j-1, end-j)+z;
    }
    assert(r < s->len[h]); *row = r;
    return id;
}
static inline size_t l5_index(const L5 *x, const L5Tree *t) {
    size_t r; int id = l5_encode(x, t, 0, &r);
    return x->s[id].off[l5_h(x,t->g[0])] + r;
}
static inline void l5_decode(const L5 *x, L5Tree *t, int id, int h, int p, size_t row) {
    const L5Shape *s = &x->s[id];
    size_t keys[4];
    for (int j = s->a, first; j > 0; j = first) {
        int k = s->child[j-1];
        for (first = j-1; first && s->child[first-1] == k; --first) {}
        size_t max = x->s[k].len[1]+j-first-1;
        size_t radix = l5_choose(max, j-first), m = row%radix;
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

/* Map a row of physical controls (including INF) to coarse endpoints. */
static inline int l5_project(const L5 *x, L5Tree *t) {
    for (int i = 0; i < t->n; ++i) {
        int g = t->g[i], kind = x->c->kind[g];
        if (!kind && x->c->D[g] < x->dc) continue;
        if (l5_degree(t,i) != 1) return 0;
        t->g[i] = l5_big(i != 0, x->buckets == 2 && kind != 0);
    }
    return 1;
}
static inline size_t l5_newborn(const L5 *x, L5Tree t, int at, int newborn) {
    int g = t.g[at];
    if (g < 0 || !x->c->L[g]) return SIZE_MAX;
    t.g[at] = b3f_next_bud(x->c,g);
    t.n = 5; t.p[4] = at; t.g[4] = newborn;
    return l5_index(x,&t);
}
static inline void l5_couple(const L5 *x, L5Tree t, int newborn, size_t birth[4], int ext_den[4]) {
    int supported = l5_project(x,&t);
    for (int i = 0; i < 4; ++i) {
        birth[i] = supported ? l5_newborn(x,t,i,newborn) : SIZE_MAX;
        ext_den[i] = supported && t.g[i] >= 0 && l5_selected(x,&t,i) ? x->dc : 2;
    }
}
static inline int l5_delete(const L5Tree *t, int z, L5Tree *m) {
    int u = z ? t->p[z] : 1;
    m->n = 0;
    for (int v = 0; v < 5; ++v) if (v != z) {
        int p = t->p[v];
        m->g[m->n] = t->g[v]; m->p[m->n++] = p == z ? -1 : p-(p > z);
    }
    return u-(u > z);
}

/* L4 addresses for the two controls of a deletion. Returning SIZE_MAX for
 * a cutoff miss makes the after-value exactly zero in either implementation. */
static inline size_t l5_l4_index(const B3FLayout *c, const L5Tree *t, int *table) {
    const int *g = t->g;
    int n[4] = {0}, child[4][3];
    for (int i = 1; i < 4; ++i) { int p = t->p[i]; child[p][n[p]++] = i; }
    for (int i = 0; i < 4; ++i)
        if (n[i]+!!i > 1 && c->kind[g[i]]) return SIZE_MAX;
    if (n[0] == 3) {
        *table = 3;
        int u = c->s3_lut[c->H[g[0]]][c->K[g[0]]][c->L[g[0]]];
        return u < 0 ? SIZE_MAX : b3f_s3_off(c,u,g[1],g[2],g[3]);
    }
    int b = child[0][0];
    if (n[0] == 2) {
        *table = 2;
        int a = child[0][1];
        if (n[a]) { int z = b; b = a; a = z; }
        int u = c->ybu_lut[c->H[g[0]]][c->K[g[0]]][c->L[g[0]]];
        int w = c->ybw_lut[c->K[g[b]]][c->L[g[b]]];
        return u < 0 || w < 0 ? SIZE_MAX : b3f_yb_off(c,u,w,g[a],g[child[b][0]]);
    }
    if (n[b] == 2) {
        *table = 1;
        int u = c->yt_lut[c->K[g[b]]][c->L[g[b]]];
        return u < 0 ? SIZE_MAX : b3f_yt_off(c,g[0],u,g[child[b][0]],g[child[b][1]]);
    }
    *table = 0;
    int mid = child[b][0];
    int b4 = c->p4_lut[c->K[g[b]]][c->L[g[b]]];
    int c4 = c->p4_lut[c->K[g[mid]]][c->L[g[mid]]];
    return b4 < 0 || c4 < 0 ? SIZE_MAX : b3f_p4_off(c,g[0],b4,c4,g[child[mid][0]]);
}
static inline uint64_t l5_source_key(const L5Tree *m, int u) {
    uint64_t key = (uint64_t)u+1;
    for (int i = 1; i < 4; ++i) key = (key<<2)|(unsigned)m->p[i];
    for (int i = 0; i < 4; ++i) {
        assert(m->g[i] >= -4 && m->g[i]+4 < 8192);
        key = (key<<13)|(unsigned)(m->g[i]+4);
    }
    return key;
}
static inline size_t l5_cache_index(uint64_t h) {
    h = (h^(h>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    h = (h^(h>>27))*UINT64_C(0x94d049bb133111eb);
    return (h^(h>>31)) & (L5_CACHE-1);
}

#endif
