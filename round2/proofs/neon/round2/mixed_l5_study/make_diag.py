#!/usr/bin/env python3
"""Make the slot attribution diagnostic from the previous coupled solver."""
from pathlib import Path
import re

root = Path(__file__).resolve().parent
s = (root / '../level5_study/diag_slack.c').read_text()
s = s.replace('#include "proto_l5.c"', '#include "proto_l5_adapted.c"')
s = s.replace('static int scenario, in_h, proto_d5, split;', '''static int scenario, in_h, proto_d5, split, slot_mask, slot_denom, coarse, coarse_mask = 3;
static inline double slot_scale(int bit) {
    return (slot_mask & bit) ? (slot_denom ? 2.0 / slot_denom : 0.0) : 1.0;
}''')
s = s.replace('return split == 2 ? 0 : l4_cap(n);', 'return split == 2 ? 0 : l4_cap(n) * slot_scale(4);')
s = s.replace('else if (!strcmp(argv[i], "--count"))', '''else if (!strcmp(argv[i], "--slot-mask") && i + 1 < argc) slot_mask = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--slot-denom") && i + 1 < argc) slot_denom = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--coarse") && i + 1 < argc) coarse = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--coarse-mask") && i + 1 < argc) coarse_mask = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--count"))''')
s = s.replace('if (split < 0', 'if (coarse && (proto_d5 != 6 || coarse != 6 || split || slot_mask || coarse_mask < 0 || coarse_mask > 3 || H.D0 <= coarse || H.D0 > 20)) return 2;\n    if (split < 0')
parts = re.split(r'(?=        /\* 5[abcd]\. Solve )', s)
assert len(parts) == 5
for j, bits in enumerate(((2,1,1,2), (2,1,2,2), (1,2,1,2), (1,2,2,2)), 1):
    i = iter(bits)
    parts[j], n = re.subn(r'(double r\w+_ext = r\w+_ext_base);',
                          lambda m: f'{m[1]} * slot_scale({next(i)});', parts[j])
    assert n == 4
    # A/B/C are partitions of the remaining mixed slack only.
    p = parts[j]
    caps = re.findall(r'double (r\w+_ext) = r\w+_ext_base', p)
    assert len(caps) == 4
    parts[j], n = re.subn(r'(double scale = eligible \? \(split == 1 \? 0 : 2\.0/proto_d5\) : \(split == 2 \? 0 : 1\);)',
                          r'\1\n                                if (eligible && slot_mask) {\n' +
                          ''.join(f'                                    {m} = {m}_base;\n' for m in caps) +
                          '                                }', parts[j])
    assert n == 1
    p = parts[j]
    m = re.search(r'L5Tree t = (\{4,[^;]+);', p)
    assert m
    buds = re.findall(r'(q\w+_b) \+= l5_birth\(l5, t, (\d), iw_nb\)', p)
    assert len(buds) == 4
    body = f'''if (l5 && coarse) {{
                                L5Tree t = {m[1]};
                                if (coarse_project4(l5, &t)) {{
'''
    for (bud, idx), cap in zip(buds, caps):
        body += f'''                                    if (t.g[{idx}] >= 0) {{
                                        {bud} += l5_birth(l5, t, {idx}, iw_nb);
                                        if (coarse_selected(&t, {idx})) {cap} *= 2.0/proto_d5;
                                    }}
'''
    body += '''                                }
                            } else if (l5) {'''
    assert p.count('if (l5) {') == 1
    parts[j] = p.replace('if (l5) {', body)
s = ''.join(parts)
(root / 'diag_slack.c').write_text(s)

p = (root / '../level5_study/proto_l5.c').read_text()

def sub(old, new):
    global p
    assert p.count(old) == 1, old
    p = p.replace(old, new)

sub('typedef struct {\n    int d, count', '''typedef struct { uint64_t key; double val; } CoarseCache;
typedef struct {
    CoarseCache *cache;
    int nf[3], *fiber[3]; /* apex h=0, apex h=1, descendant h=1 */
    int d, count''')
sub('static L5 *l5_init', '''/* BIG uses negative ids with actual apex h encoded in the sentinel. */
static inline int l5_h(const L5 *x, int g) { return g < 0 ? -g-1 : x->c->H[g]; }
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
        t->g[i] = -1;
        if (i != 0) t->g[i] = -2;
    }
    return 1;
}

static L5 *l5_init''')
sub('            size_t z = s->ns[h];', '''            /* A BIG apex uses h=0 in the index; h never affects a BIG
               transition or bound. A nonapex is necessarily h=1. */
            if (coarse && ((h==1 && s->a==0) || (h==0 && s->a==1)))
                s->state[h][s->ns[h]++] = -1-h;
            size_t z = s->ns[h];''')
sub('    printf("L5 D5=', '''    if (coarse) {
        x->cache = b3f_xcalloc(1u<<19, sizeof(*x->cache));
        for (int r=0; r<3; ++r) {
            x->fiber[r] = b3f_xmalloc(c->NG * sizeof(int));
            for (int g=0; g<c->NG; ++g) {
                int k=c->kind[g];
                int valid = r==2 ? c->w_to_wi[g]>=0 : c->p_to_pi[g]>=0;
                if (valid && (k==1 || (r==0 || (c->H[g]==1)) && (k || c->D[g]>=coarse)))
                    x->fiber[r][x->nf[r]++] = g;
            }
        }
    }
    if (coarse) {
        size_t bytes = x->entries*sizeof(double) + sizeof(*x)
                     + (1u<<19)*sizeof(*x->cache) + 3*c->NG*sizeof(int);
        printf("COARSE allocated_payload_bytes=%zu coarse=%d mask=%d fibers=%d/%d/%d\\n",
               bytes, coarse, coarse_mask, x->nf[0], x->nf[1], x->nf[2]);
    }
    printf("L5 D5=''')
sub('    free(x);', '    free(x->cache);\n    for (int r=0; r<3; ++r) free(x->fiber[r]);\n    free(x);')
sub('int id = l5_find(x, a, kid), g = t->g[u], h = x->c->H[g];\n    assert(id >= 0 && x->s[id].ix[h][g] >= 0);\n    size_t r = x->s[id].ix[h][g];', '''int id = l5_find(x, a, kid), g = t->g[u], h = l5_h(x,g);
    assert(id >= 0);
    int idx = g < 0 ? x->s[id].ns[h]-1 : x->s[id].ix[h][g];
    assert(idx >= 0 && x->s[id].state[h][idx]==g);
    size_t r = idx;''')
sub('static double l5_read(L5 *x, L5Tree *t) {', '''#include "coarse_source.c"

static double l5_read(L5 *x, L5Tree *t) {''')
sub('    double *out = &x->s[id].q[c->H[t->g[0]]][row];', '    double *out = &x->s[id].q[l5_h(x,t->g[0])][row];')
sub('        double w = 1.0/c->D[t->g[z]] - 1.0/x->d;', '''        if (t->g[z]<0) continue;
        double w = 1.0/c->D[t->g[z]] - 1.0/x->d;''')
sub('        double before = l5_q4(x, &m);', '''        if (coarse) {
            if (coarse_selected(&m, v)) num += w * coarse_source(x, &m, v);
            continue;
        }
        double before = l5_q4(x, &m);''')
sub('        double r = c->L[g];\n        if (i) r += 1.0/c->D[t->g[t->p[i]]];\n        for (int j = i+1; j < 5; ++j) if (t->p[j] == i) r += 1.0/c->D[t->g[j]];', '''        if (g<0) { up[i]=slack[i]=0; continue; }
        double r = c->L[g], extra = 0;
        for (int j=0; j<5; ++j) if ((i && t->p[i]==j) || t->p[j]==i) {
            if (t->g[j]>=0) r+=1.0/c->D[t->g[j]];
            else extra+=1.0/coarse;
        }''')
sub('        } else up[i] = 0;', '''        } else if (coarse && neighbors[i]==1) {
            t->g[i] = i ? -2 : -1; up[i] = l5_read(x, t); t->g[i]=g;
        } else up[i] = 0;''')
sub('slack[i] = 0.5*', 'slack[i] = extra + 0.5*')
sub('    x->visited = 0;', '    x->visited = 0;\n    if (coarse) memset(x->cache, 0, (1u<<19)*sizeof(*x->cache));')
sub('    return x->s[id].q[c->H[t.g[0]]][row];', '    return x->s[id].q[l5_h(x,t.g[0])][row];')
(root / 'proto_l5_adapted.c').write_text(p)
