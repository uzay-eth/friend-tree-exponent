/* In-place triangular solve. Within a shape a gain increases degree/bucket
 * and a promotion increases k; all dependencies preserve the rooted shape. */
typedef struct { uint64_t key; double value; } L5FloatCache;
#define L5_CORNER(rate,up,slot) do { \
    if ((rate) > 0 && (mode == 1 ? (up) > v : (up) < v)) { \
        double cap = 2.0*(rate)/ed[slot]; num += cap*(up); den += cap; m |= 1<<(slot); \
    } \
} while (0)
typedef struct {
    const L5 *x;
    const double *l4[4];
    double *q, lam, eps;
    int mode;
    L5FloatCache *cache;
} L5Solve;

static void l5_float_fiber(const L5Solve *s, L5Tree *m, int u, int i, double *ext) {
    while (i < 4 && m->g[i] >= 0) ++i;
    if (i == 4) {
        int g = m->g[u], a, b;
        size_t before = l5_l4_index(s->x->c,m,&a);
        m->g[u] = b3f_next_up(s->x->c,g);
        size_t after = l5_l4_index(s->x->c,m,&b);
        m->g[u] = g;
        double d = (after == SIZE_MAX ? 0 : s->l4[b][after]) -
                   (before == SIZE_MAX ? 0 : s->l4[a][before]);
        if (s->mode == 1 ? d > *ext : d < *ext) *ext = d;
        return;
    }
    int big = m->g[i], b = l5_bucket(big), r = !!i;
    for (int j = 0; j < s->x->nf[r][b]; ++j) {
        m->g[i] = s->x->fiber[r][b][j];
        l5_float_fiber(s,m,u,i+1,ext);
    }
    m->g[i] = big;
}

static double l5_solve_row(L5Solve *s, L5Tree *t, size_t row) {
    if (!isnan(s->q[row])) return s->q[row];
    const L5 *x = s->x;
    const B3FLayout *c = x->c;
    double num = s->eps, den = s->lam, up[5] = {0}, slack[5] = {0};
    int degree[5];
    for (int i = 0; i < 5; ++i) degree[i] = l5_degree(t,i);
    for (int z = 0; z < 5; ++z) if (degree[z] == 1 && t->g[z] >= 0) {
        L5Tree m; int u = l5_delete(t,z,&m);
        if (!l5_selected(x,&m,u)) continue;
        uint64_t key = l5_source_key(&m,u);
        L5FloatCache *e = s->cache+l5_cache_index(key);
        if (e->key != key) {
            double v = s->mode == 1 ? -INFINITY : INFINITY;
            l5_float_fiber(s,&m,u,0,&v);
            *e = (L5FloatCache){key,v};
        }
        num += (1.0/c->D[t->g[z]] - 1.0/x->dc)*e->value;
    }
    for (int i = 0; i < 5; ++i) {
        int g = t->g[i];
        if (g < 0) {
            if (x->buckets == 2 && !l5_bucket(g)) {
                t->g[i] = l5_big(i != 0,1);
                up[i] = l5_solve_row(s,t,l5_index(x,t));
                t->g[i] = g;
                slack[i] = c->D0-1;
            }
            continue;
        }
        double rate = c->L[g];
        slack[i] = 0.5*(c->K[g]+c->H[g]-degree[i]);
        for (int j = 0; j < 5; ++j) if (t->p[i] == j || t->p[j] == i) {
            int w = t->g[j];
            if (w >= 0) rate += 1.0/c->D[w];
            else slack[i] += 1.0/(l5_bucket(w) ? c->D0 : x->dc);
        }
        if (c->D[g]+1 < x->dc || degree[i] == 1) {
            t->g[i] = c->D[g]+1 < x->dc ? b3f_next_up(c,g) : l5_big(i != 0,0);
            up[i] = l5_solve_row(s,t,l5_index(x,t));
            t->g[i] = g;
        }
        num += rate*up[i]; den += rate;
        if (c->L[g]) {
            double b = (double)c->L[g]/c->D[g];
            t->g[i] = b3f_next_bud(c,g);
            num += b*l5_solve_row(s,t,l5_index(x,t)); den += b;
            t->g[i] = g;
        }
    }
    /* Greedily include optional rates in order of their target value. The
     * updated weighted average lies between its previous value and the
     * target, so no previously included rate needs reconsideration. */
    for (int j = 0; j < 5; ++j) {
        int best = -1;
        for (int i = 0; i < 5; ++i) if (slack[i] &&
            (best < 0 || (s->mode == 1 ? up[i] > up[best] : up[i] < up[best]))) best = i;
        if (best < 0 || (s->mode == 1 ? up[best] <= num/den : up[best] >= num/den)) break;
        num += slack[best]*up[best]; den += slack[best]; slack[best] = 0;
    }
    return s->q[row] = num/den;
}

static void l5_sweep(L5Solve base, int threads) {
    #pragma omp parallel for num_threads(threads)
    for (size_t j = 0; j < base.x->rows; ++j) base.q[j] = NAN;
    #pragma omp parallel num_threads(threads)
    {
        L5Solve s = base;
        s.cache = b3f_xcalloc(L5_CACHE,sizeof *s.cache);
        #pragma omp for schedule(dynamic,1)
        for (int i = 0; i < s.x->count; ++i) if (s.x->s[i].n == 5) {
            for (int h = 0; h < 2; ++h) for (size_t j = 0; j < s.x->s[i].len[h]; ++j) {
                size_t row = s.x->s[i].off[h]+j;
                if (isnan(s.q[row])) {
                    L5Tree t = {0};
                    l5_decode(s.x,&t,i,h,-1,j);
                    l5_solve_row(&s,&t,row);
                }
            }
        }
        free(s.cache);
    }
}
