/* Only a fully evaluated exact extremum is cached. The difference of two
 * validated entries is within 2*MAX_VAL, so storing it as int64 is lossless. */
typedef struct { uint64_t key; int64_t value; } L5ExactCache;
typedef struct {
    const L5 *x;
    const int64_t *l4[4], *q;
    int mode;
    L5ExactCache *cache;
    mpz_t sum, unit, lambda, diff, ext, tmp;
    mpz_t reciprocals[CERT_B3F_MAX_D0+1];
    int have;
} L5Exact;

static void l5_exact_fiber(L5Exact *s, L5Tree *m, int u, int i) {
    while (i < 4 && m->g[i] >= 0) ++i;
    if (i == 4) {
        int g = m->g[u], a, b;
        size_t before = l5_l4_index(s->x->c,m,&a);
        m->g[u] = b3f_next_up(s->x->c,g);
        size_t after = l5_l4_index(s->x->c,m,&b);
        m->g[u] = g;
        mpz_set_si(s->diff,after == SIZE_MAX ? 0 : s->l4[b][after]);
        mpz_set_si(s->tmp,before == SIZE_MAX ? 0 : s->l4[a][before]);
        mpz_sub(s->diff,s->diff,s->tmp);
        int cmp = mpz_cmp(s->diff,s->ext);
        if (!s->have || (s->mode == 1 ? cmp > 0 : cmp < 0)) mpz_set(s->ext,s->diff);
        s->have = 1;
        return;
    }
    int big = m->g[i], b = l5_bucket(big), r = !!i;
#ifdef L5_PLANT_BUCKET
    if (b == 1) b = 0;
#endif
    for (int j = 0; j < s->x->nf[r][b]; ++j) {
        m->g[i] = s->x->fiber[r][b][j];
        l5_exact_fiber(s,m,u,i+1);
    }
    m->g[i] = big;
}

static int64_t l5_exact_source(L5Exact *s, L5Tree *m, int u) {
    uint64_t key = l5_source_key(m,u);
    L5ExactCache *e = s->cache+l5_cache_index(key);
    if (e->key != key) {
        s->have = 0;
        l5_exact_fiber(s,m,u,0);
        assert(s->have && mpz_fits_slong_p(s->ext));
        *e = (L5ExactCache){key,mpz_get_si(s->ext)};
    }
    return e->value;
}

static inline int64_t l5_qint(const L5 *x, const int64_t *q, const L5Tree *t) {
    return q[l5_index(x,t)];
}
static inline void l5_addmul(L5Exact *s, const mpz_t unit, int64_t k) {
    if (k >= 0) mpz_addmul_ui(s->sum,unit,(unsigned long)k);
    else mpz_submul_ui(s->sum,unit,(unsigned long)-k);
}
static inline int l5_corner(const L5Exact *s, int64_t jump) {
    return s->mode == 1 ? jump > 0 : jump < 0;
}

/* Accumulate the numerator over Q*lcm(1,...,D0). Splitting the transfer
 * coefficient into 1/d - 1/Dc avoids a needless product denominator. */
static int l5_exact_row(L5Exact *s, L5Tree *t, size_t row) {
    const L5 *x = s->x;
    const B3FLayout *c = x->c;
    int64_t v = s->q[row];
    mpz_mul_si(s->sum,s->lambda,v);
    int degree[5];
    for (int i = 0; i < 5; ++i) degree[i] = l5_degree(t,i);
    for (int z = 0; z < 5; ++z) if (degree[z] == 1 && t->g[z] >= 0) {
        L5Tree m; int u = l5_delete(t,z,&m);
        if (!l5_selected(x,&m,u)) continue;
        int64_t ext = l5_exact_source(s,&m,u);
#ifdef L5_PLANT_SOURCE
        ext = -ext;
#endif
        l5_addmul(s,s->reciprocals[c->D[t->g[z]]],ext);
        l5_addmul(s,s->reciprocals[x->dc],-ext);
    }
    for (int i = 0; i < 5; ++i) {
        int g = t->g[i];
        if (g < 0) {
            if (x->buckets == 2 && !l5_bucket(g)) {
                t->g[i] = l5_big(i != 0,1);
                int64_t delta = l5_qint(x,s->q,t)-v;
                t->g[i] = g;
                if (l5_corner(s,delta)) l5_addmul(s,s->unit,(c->D0-1)*delta);
            }
            continue;
        }
        int d = c->D[g], l = c->L[g];
        int64_t delta = -v;
        if (d+1 < x->dc || degree[i] == 1) {
            t->g[i] = d+1 < x->dc ? b3f_next_up(c,g) : l5_big(i != 0,0);
            delta += l5_qint(x,s->q,t);
            t->g[i] = g;
        }
        l5_addmul(s,s->unit,l*delta);
        int corner = l5_corner(s,delta);
        if (corner) l5_addmul(s,s->reciprocals[2],(c->K[g]+c->H[g]-degree[i])*delta);
        for (int j = 0; j < 5; ++j) if (t->p[i] == j || t->p[j] == i) {
            int w = t->g[j];
            if (w >= 0) l5_addmul(s,s->reciprocals[c->D[w]],delta);
            else if (corner) l5_addmul(s,s->reciprocals[l5_bucket(w) ? c->D0 : x->dc],delta);
        }
        if (l) {
            t->g[i] = b3f_next_bud(c,g);
            l5_addmul(s,s->reciprocals[d],l*(l5_qint(x,s->q,t)-v));
            t->g[i] = g;
        }
    }
    return s->mode == 1 ? mpz_sgn(s->sum) > 0 : mpz_sgn(s->sum) < 0;
}

static long long verify_l5(const CertData *d, int threads) {
    const L5 *x = d->X;
    long long fails = 0;
    #pragma omp parallel num_threads(threads) reduction(+:fails)
    {
        L5Exact s = {.x=x,.l4={d->qP4,d->qYt,d->qYb,d->qS3},.q=d->q5,.mode=d->H->mode};
        s.cache = b3f_xcalloc(L5_CACHE,sizeof *s.cache);
        mpz_inits(s.sum,s.unit,s.lambda,s.diff,s.ext,s.tmp,NULL);
        mpz_set_ui(s.unit,1);
        for (int i = 2; i <= x->c->D0; ++i) mpz_lcm_ui(s.unit,s.unit,(unsigned long)i);
        mpz_mul_si(s.lambda,s.unit,-d->H->P);
        mpz_mul_si(s.unit,s.unit,d->H->Q);
        for (int i = 1; i <= x->c->D0; ++i) {
            mpz_init(s.reciprocals[i]); mpz_divexact_ui(s.reciprocals[i],s.unit,(unsigned long)i);
        }
        for (int i = 0; i < x->count; ++i) if (x->s[i].n == 5) for (int h = 0; h < 2; ++h) {
            #pragma omp for schedule(static)
            for (size_t j = 0; j < x->s[i].len[h]; ++j) {
                L5Tree t = {0};
                l5_decode(x,&t,i,h,-1,j);
                fails += l5_exact_row(&s,&t,x->s[i].off[h]+j);
            }
        }
        for (int i = 1; i <= x->c->D0; ++i) mpz_clear(s.reciprocals[i]);
        mpz_clears(s.sum,s.unit,s.lambda,s.diff,s.ext,s.tmp,NULL);
        free(s.cache);
    }
    return fails;
}

static inline int64_t newborn_int(const CertData *d, size_t i) {
#ifdef L5_PLANT_NEWBORN
    (void)d; (void)i;
    return 0;
#else
    return i == SIZE_MAX ? 0 : d->q5[i];
#endif
}
