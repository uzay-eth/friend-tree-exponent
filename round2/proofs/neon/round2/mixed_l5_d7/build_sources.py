#!/usr/bin/env python3
"""Apply the experiment delta to the previous study's generated solver."""
from pathlib import Path

HERE = Path(__file__).resolve().parent
OLD = HERE.parent / 'mixed_l5_study'


def edit(name, changes):
    text = (OLD / name).read_text()
    for before, after in changes:
        assert text.count(before) == 1, (name, before)
        text = text.replace(before, after)
    (HERE / name).write_text(text)


edit('diag_slack.c', [
    ('coarse_mask = 3;', 'coarse_mask = 3, coarse_buckets = 1;'),
    ('        else if (!strcmp(argv[i], "--coarse-mask")',
     '        else if (!strcmp(argv[i], "--coarse-buckets") && i + 1 < argc) coarse_buckets = atoi(argv[++i]);\n'
     '        else if (!strcmp(argv[i], "--coarse-mask")'),
    ('    if (coarse && (proto_d5 != 6 || coarse != 6 || split || slot_mask || coarse_mask < 0 || coarse_mask > 3 || H.D0 <= coarse || H.D0 > 20)) return 2;',
     '    if (coarse_buckets < 1 || coarse_buckets > 2 || (!coarse && coarse_buckets != 1)) return 2;\n'
     '    if (coarse && (proto_d5 != coarse || coarse < 6 || coarse > 8 || split || slot_mask || coarse_mask < 0 || coarse_mask > 3 || H.D0 <= coarse || H.D0 > 24)) return 2;'),
])

edit('proto_l5_adapted.c', [
    (' * q5 is zero outside D<D5, and the residual row uses the original 1/2 cap. */',
     ' * Coarse modes also support one or two large endpoint buckets. */'),
    ('    int nf[3], *fiber[3]; /* apex h=0, apex h=1, descendant h=1 */',
     '    int nf[3][2], *fiber[3]; /* apex (both h), apex h=1, descendant h=1; inner then tail */'),
    ('''/* BIG uses negative ids with actual apex h encoded in the sentinel. */
static inline int l5_h(const L5 *x, int g) { return g < 0 ? -g-1 : x->c->H[g]; }''',
     '''/* -1/-2 = BIG_INNER (or unsplit BIG), -3/-4 = TAIL; parity encodes index h. */
static inline int coarse_id(int h, int b) { return -1-h-2*b; }
static inline int coarse_bucket(int g) { return (-g-1)/2; }
static inline int l5_h(const L5 *x, int g) { return g < 0 ? (-g-1)%2 : x->c->H[g]; }'''),
    ('''        t->g[i] = -1;
        if (i != 0) t->g[i] = -2;''',
     '''        t->g[i] = coarse_id(i != 0, coarse_buckets == 2 && x->c->kind[g] != 0);'''),
    ('''            if (coarse && ((h==1 && s->a==0) || (h==0 && s->a==1)))
                s->state[h][s->ns[h]++] = -1-h;''',
     '''            if (coarse && ((h==1 && s->a==0) || (h==0 && s->a==1)))
                for (int b=0; b<coarse_buckets; ++b) {
                    assert(s->ns[h] < 64);
                    s->state[h][s->ns[h]++] = coarse_id(h,b);
                }'''),
    ('''            for (int g=0; g<c->NG; ++g) {
                int k=c->kind[g];
                int valid = r==2 ? c->w_to_wi[g]>=0 : c->p_to_pi[g]>=0;
                if (valid && (k==1 || (r==0 || (c->H[g]==1)) && (k || c->D[g]>=coarse)))
                    x->fiber[r][x->nf[r]++] = g;
            }''',
     '''            for (int b=0, at=0; b<coarse_buckets; ++b) for (int g=0; g<c->NG; ++g) {
                int k=c->kind[g];
                int valid = r==2 ? c->w_to_wi[g]>=0 : c->p_to_pi[g]>=0;
                if (coarse_buckets==2 && b!=(k!=0)) continue;
                if (valid && (k==1 || ((r==0 || c->H[g]==1) && (k || c->D[g]>=coarse)))) {
                    x->fiber[r][at++] = g;
                    ++x->nf[r][b];
                }
            }'''),
    ('''        printf("COARSE allocated_payload_bytes=%zu coarse=%d mask=%d fibers=%d/%d/%d\\n",
               bytes, coarse, coarse_mask, x->nf[0], x->nf[1], x->nf[2]);''',
     '''        printf("COARSE allocated_payload_bytes=%zu coarse=%d mask=%d buckets=%d apex_fibers=%d/%d child_fibers=%d/%d\\n",
               bytes, coarse, coarse_mask, coarse_buckets,
               x->nf[0][0], x->nf[0][1], x->nf[2][0], x->nf[2][1]);'''),
    ('g < 0 ? x->s[id].ns[h]-1 :', 'g < 0 ? x->s[id].ns[h]-coarse_buckets+coarse_bucket(g) :'),
    ('''        if (g<0) { up[i]=slack[i]=0; continue; }''',
     '''        if (g<0) {
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
        }'''),
    ('else extra+=1.0/coarse;', 'else extra+=1.0/(coarse_bucket(t->g[j]) ? c->D0 : coarse);'),
    ('t->g[i] = i ? -2 : -1;', 't->g[i] = coarse_id(i != 0,0);'),
])

edit('coarse_source.c', [
    ('''    for (int j=0; j<x->nf[r]; ++j) {''',
     '''    int b=coarse_bucket(big), start=b ? x->nf[r][0] : 0;
    for (int j=start; j<start+x->nf[r][b]; ++j) {'''),
    ('''        assert(m->g[i]+2<8192);
        key=(key<<13)|(unsigned)(m->g[i]+2);''',
     '''        assert(m->g[i]>=-4 && m->g[i]+4<8192);
        key=(key<<13)|(unsigned)(m->g[i]+4);'''),
])
