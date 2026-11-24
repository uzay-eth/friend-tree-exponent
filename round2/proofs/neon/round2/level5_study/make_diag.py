#!/usr/bin/env python3
"""Generate diag_slack.c from the pinned local builder, with checked substitutions."""
from pathlib import Path
import re

root = Path(__file__).resolve().parent
s = (root / '../cert/cert_b3f_build.c').read_text()


def sub(old, new, n=1):
    global s
    assert s.count(old) == n, (old, s.count(old), n)
    s = s.replace(old, new)


sub('#include "cert_b3f_format.h"', '''#include "../cert/cert_b3f_format.h"

/* One ablation per invocation. These hypothetical rows are NOT certificates. */
enum { BASE, CENTER, D5_6, D5_7, D5_8, H_CENTER };
static int scenario, in_h, proto_d5, split;
static inline double l4_cap(int n) {
    return scenario == CENTER ? 0.0 : n / (double)(scenario >= D5_6 && scenario <= D5_8 ? scenario + 4 : 2);
}
static inline double l4_tail_cap(int n) {
    return split == 2 ? 0 : l4_cap(n);
}
static inline double l4_fixed(int np, int nc, double rp, double rc) {
    return scenario == CENTER ? np * rp + nc * rc : 0.0;
}
static inline void add_fixed(double r, double v, double *num, double *den) {
    *num += r * v; *den += r;
}''')
sub('double Rbar_ext, double Y_ref,', 'double Rbar_ext, double Y_ref, double ext_lo,')
sub('double rexts[2] = { 0.0, rho_ext_max };\n            int nrext = (rho_ext_max > 1e-15) ? 2 : 1;', '''double rexts[2] = { ext_lo, rho_ext_max };
            int nrext = (rho_ext_max - ext_lo > 1e-15) ? 2 : 1;
            if (in_h) {
                rins[0] = (r_in_lo + r_in_hi) * 0.5; nrin = 1;
                rexts[0] = rho_ext_max * 0.5; nrext = 1;
            }''')
sub('double rn = BN[j] / BB[j] - shift, rd = BD[j] / BB[j], val = rn - rd * v;', '''/* At t=tmax/2 use the actual rational row, not a hull extremum. */
                    if (in_h) {
                        const int w[5] = {1, 4, 6, 4, 1};
                        double nn = 0, dd = 0, bb = 0;
                        for (int l = 0; l < 5; ++l) { nn += w[l]*BN[l]; dd += w[l]*BD[l]; bb += w[l]*BB[l]; }
                        best_n = nn / bb - shift; best_d = dd / bb;
                        break;
                    }
                    double rn = BN[j] / BB[j] - shift, rd = BD[j] / BB[j], val = rn - rd * v;''')
sub('double a_lo = -((double)(c + 1) - rho_min) / D0 * d1;', '''if (in_h) rho_min = rho_max = 0.5 * (rho_min + rho_max);
            double a_lo = -((double)(c + 1) - rho_min) / D0 * d1;''')
sub('if (kx > KT && ((mode == 1 && F0 > v)', 'if (in_h && kx > KT) { *num += 0.5 * F0; *den += 0.5; }\n        else if (kx > KT && ((mode == 1 && F0 > v)')
sub('if (kx > KT && ((mode == 1 && F1 > v)', 'if (in_h && kx > KT) { *num += F1; *den += 1; }\n        else if (kx > KT && ((mode == 1 && F1 > v)')
sub('double best_val = 0, best_n = 0, best_d = 0; int bid = 0;', '''if (in_h) {
                double bb = BB2[0] + 2*BB2[1] + BB2[2];
                *num += 0.5*(BN2[0] + 2*BN2[1] + BN2[2])/bb;
                *den += 0.5*(BD2[0] + 2*BD2[1] + BD2[2])/bb;
                if (n_nbr_in > 0) {
                    *num += 0.5*(d0_nbr - 3*d1_nbr + 3*d2_nbr + F3_up_sum)/dD0;
                    *den += 0.5*n_nbr_in/dD0;
                }
                return;
            }
            double best_val = 0, best_n = 0, best_d = 0; int bid = 0;''')
sub('/* 3. Solve FULL hV */','in_h = scenario == H_CENTER;\n        /* 3. Solve FULL hV */')
sub('/* 5a. Solve P4 */','in_h = 0;\n        /* 5a. Solve P4 */')
sub('static void write_i64_array', '#include "proto_l5.c"\n\nstatic void write_i64_array')

# (row state, external parent count, external child count, up value, tail-centering index)
slots = [
    [('a','ha','ka - 1','qa_up','b3f_tidx(&C, ka, ha)'), ('b_g','0','kb - 1','qb_up',None), ('c_g','0','kc - 1','qc_up',None), ('d','0','kd','qd_up','b3f_tidx(&C, kd, 1)')],
    [('p','hp','kp - 1','qp_up','b3f_tidx(&C, kp, hp)'), ('u_g','0','ku - 2','qu_up',None), ('w1','0','k1','q1_up','b3f_tidx(&C, k1, 1)'), ('w2','0','k2','q2_up','b3f_tidx(&C, k2, 1)')],
    [('u_g','hu','ku - 2','qu_up',None), ('w1','0','k1','q1_up','b3f_tidx(&C, k1, 1)'), ('w2_g','0','k2 - 1','q2_up',None), ('c','0','kc','qc_up','b3f_tidx(&C, kc, 1)')],
    [('u_g','hu','ku - 3','qu_up',None), ('w1','0','k1','q1_up','b3f_tidx(&C, k1, 1)'), ('w2','0','k2','q2_up','b3f_tidx(&C, k2, 1)'), ('w3','0','k3','q3_up','b3f_tidx(&C, k3, 1)')],
]

parts = re.split(r'(?=        /\* 5[abcd]\. Solve )', s)
assert len(parts) == 5
parts[0], n = re.subn(r'(b3_tail_var\(mode[^\n]+?,\s*(?:p_Rbar, p_Yref|w_Rbar, w_Yref|0, 0)), ', r'\1, 0, ', parts[0]); assert n == 8, n
for part, spec in enumerate(slots, 1):
    b = parts[part]
    for state, np, nc, up, tail in spec:
        ext = re.search(r'if \((r\w+_ext) > 0 && \(\(mode == 1 && ' + up + r' > v\)', b)[1]
        b, n = re.subn(r'(\b'+ext+r' = (?:\(\w+ == 0\) \? )?)0\.5 \* (\([^)]*\)|\w+)', rf'\1l4_cap(\2)', b); assert n == 1
        fixed = f'l4_fixed({np}, {nc}, rpar[{state}], rch[{state}])'
        anchor = 'double v = num0 / den0; int mp = -1;'
        assert b.count(anchor) == 1
        b = b.replace(anchor, f'if (scenario == CENTER && C.kind[{state}] == 0) add_fixed({fixed}, {up}, &num0, &den0);\n                            ' + anchor)
        if tail:
            # Known-k tail: narrow the rate corners in T2 and T3 as well.
            cap = f'l4_tail_cap(({np}) + ({nc}))'
            tf = f'l4_fixed({np}, {nc}, trpar[{tail}], trch[{tail}])'
            lines = b.splitlines(keepends=True)
            # The tail call stores the *neighbor's* up value, so match by external slot kind.
            kind = {'a':'ak','d':'dk','p':'pk','w1':'k1k','w2':'k2k','c':'ck','w3':'k3k'}[state]
            for op in ('const', 'var'):
                indices = [i for i,l in enumerate(lines) if re.search(r'b3_tail_'+op+r'\(mode, D0, KT, '+kind+r',', l)]
                assert len(indices) == 1
                i = indices[0]; l = lines[i]
                line_rate = r'0\.5 \* (\([^)]*\)|\w+)'
                val = f'({kind} == 3' + ('' if op == 'const' else f' || {kind} == 4') + f' ? {tf} : 0.0)'
                if op == 'const':
                    l = re.sub(r'(,\s*(?:ib|ic|iu|i2)), (ib|ic|iu|i2) \+ '+line_rate, rf'\1 + {val}, \2 + {val} + {cap}', l)
                else:
                    l = re.sub(line_rate+r',\s*0, 0,', f'{val} + {cap}, 0, 0, {val},', l)
                assert l != lines[i], (part, state, op, l)
                lines[i] = l
            b = ''.join(lines)
    # Cut off ONLY fully eligible rows. Mixed rows need the original bound, including tails.
    # Each promotion of an old leaf is charged to the unique preimage M=N\{newborn}.
    order, parents, bud_names = {
        1: ('a,b_g,c_g,d', '-1,0,1,2', ('qa_b','qb_b','qc_b','qd_b')),
        2: ('p,u_g,w1,w2', '-1,0,1,1', ('qp_b','qu_b','q1_b','q2_b')),
        3: ('u_g,w1,w2_g,c', '-1,0,0,2', ('qu_b','q1_b','q2_b','qc_b')),
        4: ('u_g,w1,w2,w3', '-1,0,0,0', ('qu_b','q1_b','q2_b','q3_b')),
    }[part]
    gs = order.split(',')
    eligibility = ' && '.join(f'C.kind[{g}]==0 && C.D[{g}]<proto_d5' for g in gs)
    # Only the mixed oracle also shrinks the residuals in known-k L4 tail rows.
    body = f'''if (l5) {{
                                int eligible = {eligibility};
                                if (eligible) {{
                                    L5Tree t = {{4, {{{parents}}}, {{{order}}}}};
'''
    for i, bud in enumerate(bud_names):
        body += f'                                    {bud} += l5_birth(l5, t, {i}, iw_nb);\n'
    body += '''                                }
                                double scale = eligible ? (split == 1 ? 0 : 2.0/proto_d5) : (split == 2 ? 0 : 1);
'''
    caps = re.findall(r'if \((r\w+_ext) > 0', b)
    assert len(caps) == 4
    body += ''.join(f'                                {cap} *= scale;\n' for cap in caps)
    body += '                            }\n                            '
    # Existing r_ext variables in outer loops must not be changed across rows.
    for cap in caps:
        b = re.sub(r'\b' + cap + r'\b', cap+'_base', b)
    for cap in caps:
        b = b.replace(f'if ({cap}_base > 0', f'if ({cap} > 0')
        b = b.replace(f'num += {cap}_base *', f'num += {cap} *').replace(f'den += {cap}_base;', f'den += {cap};')
    declarations = ''.join(f'double {cap} = {cap}_base;\n                            ' for cap in caps)
    b = b.replace('double num0 =', declarations+body+'double num0 =', 1)
    parts[part] = b
s = ''.join(parts)

# Each layer has exactly one last reader before its rebuild, including the final phase.
sub('int nit = 22, nthreads = 4;', 'int nit = 22, nthreads = 1, count = 0, fixed = 0, l5_check = 0;')
sub('else if (!strcmp(argv[i], "--out")', '''else if (!strcmp(argv[i], "--scenario") && i + 1 < argc) {
            const char *names[] = {"base", "center", "d5-6", "d5-7", "d5-8", "h-center"};
            const char *name = argv[++i];
            for (scenario = 0; scenario < 6 && strcmp(names[scenario], name); ++scenario) {}
            if (scenario == 6) return 2;
        }
        else if (!strcmp(argv[i], "--count")) count = 1;
        else if (!strcmp(argv[i], "--l5") && i + 1 < argc) proto_d5 = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--l5-split") && i + 1 < argc) split = atoi(argv[++i]);
        else if (!strcmp(argv[i], "--l5-check")) l5_check = 1;
        else if (!strcmp(argv[i], "--out")''')
sub('if (!out_path) { fprintf(stderr, "Missing --out <path>\\n"); return 2; }', '''if (nthreads != 1 || nit < 4 || (out_path && scenario != BASE)) return 2;
    if (proto_d5 && (proto_d5 < 6 || proto_d5 > 8 || proto_d5 >= H.D0 || scenario != BASE || out_path ||
        H.DP4 < proto_d5 || H.DYT < proto_d5 || H.DYB < proto_d5 || H.DS3 < proto_d5 ||
        H.KP4 < proto_d5-2 || H.KYT < proto_d5-2 || H.KYB < proto_d5-1 || H.KS3 < proto_d5-1)) {
        fprintf(stderr, "L5 requires every pure-inner four-vertex preimage in L4 (check D5 and L4 cutoffs)\\n");
        return 2;
    }
    if (split < 0 || split > 2 || ((split || l5_check) && !proto_d5)) return 2;
    fixed = out_path != NULL;''')
sub('H.LF      = 12000; H.KF = 160;', 'H.LF      = 6000; H.KF = 100;')
sub('H.D0      = 16;    H.KT = 15;', 'H.D0      = 12;    H.KT = 11;')
for field, old, new in [('P4','9,6','10,8'), ('YT','8,5','9,7'), ('YB','9,6','10,8'), ('S3','7,5','8,6')]:
    a,b = old.split(','); c,d = new.split(',')
    sub(f'H.D{field}     = {a};     H.K{field} = {b};', f'H.D{field}     = {c};     H.K{field} = {d};')
sub('    fflush(stdout);\n\n    int64_t *rpar_q', '''    size_t pot = C.NV + C.g_len + C.hv_len + C.hc_len + C.p4_len + C.yt_len + C.yb_len + C.s3_len;
    printf("COUNTS f=%d g=%zu hv=%zu hc=%zu p4=%zu yt=%zu yb=%zu s3=%zu potential=%zu file_bytes=%zu\\n",
           C.NV, C.g_len, C.hv_len, C.hc_len, C.p4_len, C.yt_len, C.yb_len, C.s3_len, pot, b3f_expected_bytes(&C));
    fflush(stdout);
    L5 *l5 = proto_d5 ? l5_init(&C, proto_d5) : NULL;
    if (count) { l5_free(l5); b3f_layout_free(&C); return 0; }
    if (l5_check) {
        l5_selfcheck(l5);
        l5_free(l5); b3f_layout_free(&C); return 0;
    }

    int64_t *rpar_q''')
for a, n, length in [('g','gn','g_len'), ('hV','hVn','hv_len'), ('hC','hCn','hc_len'), ('qP4','qP4n','p4_len'), ('qYt','qYtn','yt_len'), ('qYb','qYbn','yb_len'), ('qS3','qS3n','s3_len')]:
    s, k = re.subn(r'\*'+n+r'\s*= \(double\*\)calloc\(C\.'+length+r', sizeof\(double\)\)', f'*{n} = fixed ? (double*)calloc(C.{length}, sizeof(double)) : {a}', s)
    assert k == 1
    sub(f'*{n} = fixed ? (double*)calloc(C.{length}, sizeof(double)) : {a}', f'*{n} = {a}')
    sub(f'memcpy({a}, {n}, C.{length} * sizeof(double));\n', '')
    sub(f'free({a}); free({n});', f'free({a});')
sub('int total_it = nit + 5;', 'int total_it = nit + (fixed ? 5 : 0);')
sub('        double rho = fn[v_nb] / f[v_nb];', '''        if (l5) l5_sweep(l5, mode, lam, qP4n, qYtn, qYbn, qS3n);
        double rho = fn[v_nb] / f[v_nb];''')
sub('double inv_s = 1.0 / fn[v_nb];', '''double inv_s = 1.0 / fn[v_nb];
            if (l5) l5_scale(l5, inv_s);
            printf("SWEEP scenario=%d mode=%d it=%d lam=%.17g rho=%.17g elapsed=%.3f\\n", scenario, mode, it+1, lam, rho, omp_get_wtime()-start);
            fflush(stdout);''')
sub('double lam = 0.5657237;', 'double start = omp_get_wtime(), lam = 0.5657237;')
sub('    FILE *fp = fopen(out_path, "wb");', '    printf("FLOAT scenario=%d mode=%d D0=%d KT=%d lam=%.17g\\n", scenario, mode, D0, KT, lam);\n    if (fixed) {\n    FILE *fp = fopen(out_path, "wb");')
sub('    free(rpar_q); free(rch_q);', '    }\n    free(rpar_q); free(rch_q);')
sub('    b3f_layout_free(&C);\n    return 0;', '    l5_free(l5);\n    b3f_layout_free(&C);\n    return 0;')
# P4/Yb rank buckets can exceed 4096 at large cutoffs. Size for all valid inputs.
for expr in ('C.NP4 * (size_t)C.NP4', 'C.NYBU * (size_t)C.NYBW'):
    s = s.replace('int pairs[4096][2], npairs = 0;', f'int pairs[1 + {expr}][2], npairs = 0;', 1)
(root / 'diag_slack.c').write_text(s)
