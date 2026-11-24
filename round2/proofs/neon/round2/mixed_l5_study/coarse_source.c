/* Convex-hull extrema over L4 control values. The marked receiver is exact
 * in every L5 state: attach the deleted leaf and it is an interior. */
static void coarse_fiber(const L5 *x, L5Tree *m, int u, int first, double *ext) {
    int i = first;
    while (i<4 && m->g[i]>=0) ++i;
    if (i==4) {
        int g=m->g[u];
        double before=l5_q4(x,m);
        m->g[u]=b3f_next_up(x->c,g);
        double diff=l5_q4(x,m)-before;
        m->g[u]=g;
        if (x->mode==1 ? diff>*ext : diff<*ext) *ext=diff;
        return;
    }
    int big=m->g[i], r=i ? 2 : 0;
    for (int j=0; j<x->nf[r]; ++j) {
        m->g[i]=x->fiber[r][j];
        coarse_fiber(x,m,u,i+1,ext);
    }
    m->g[i]=big;
}

static double coarse_source(L5 *x, L5Tree *m, int u) {
    assert(m->g[u]>=0);
    /* Full ordered tree (including marked slot); collision means recompute. */
    uint64_t key=(uint64_t)u+1;
    for (int i=1; i<4; ++i) key=(key<<2)|(unsigned)m->p[i];
    for (int i=0; i<4; ++i) {
        assert(m->g[i]+2<8192);
        key=(key<<13)|(unsigned)(m->g[i]+2);
    }
    uint64_t h=key;
    h=(h^(h>>30))*UINT64_C(0xbf58476d1ce4e5b9);
    h=(h^(h>>27))*UINT64_C(0x94d049bb133111eb);
    h^=h>>31;
    CoarseCache *e=x->cache+(h & ((1u<<19)-1));
    if (e->key==key) return e->val;
    double ext=x->mode==1 ? -INFINITY : INFINITY;
    coarse_fiber(x,m,u,0,&ext);
    assert(isfinite(ext));
    e->key=key; e->val=ext;
    return ext;
}
