/* Untrusted, deterministic proposal generator. All proof checks read its integer output. */
#define main vintage_main
#include "../motif/exp_v5.c"
#undef main
#include "../motif/cert4v_format.h"
#include <unistd.h>
#include <errno.h>
#include <inttypes.h>

static const int64_t weight_unit = INT64_C(1) << 40;
static const int64_t center_unit = INT64_C(1) << 30;

static void require(int ok, const char *what) {
    if (!ok) { fprintf(stderr, "export error: %s\n", what); exit(1); }
}
static int64_t quantize(double x, int64_t scale) {
    require(isfinite(x) && fabs(x) < 16., "candidate magnitude outside checker range");
    return llround(x * scale);
}
static void write_array(FILE *out, const double *a, size_t n, int64_t unit) {
    int64_t buf[4096];
    for (size_t start = 0; start < n;) {
        size_t len = n - start;
        if (len > 4096) len = 4096;
        for (size_t i = 0; i < len; ++i) buf[i] = quantize(a[start+i], unit);
        require(fwrite(buf, sizeof(*buf), len, out) == len, "writing payload");
        start += len;
    }
}
static void write_snapshot(Model *m, Cert4vHeader *h, const Cert4vSizes *z, const char *file) {
    int upper = m->mode == 1;
    h->f_p0 = quantize(upper ? m->A_up : m->B1_lo, h->S);
    h->f_p1 = quantize(upper ? m->B_up : m->B2_lo, h->S);
    h->f_p2 = upper ? 0 : quantize(m->CD_lo, h->S);
    require(!c4_validate_header(h), "frozen header");
    FILE *out = fopen(file, "wbx");
    require(out != NULL, "cannot exclusively create snapshot staging file");
    require(fwrite(h, sizeof(*h), 1, out) == 1, "writing header");
    const double *cent[] = {m->rpar,m->rch,m->xpar,m->xch,m->xcc};
    for (int i = 0; i < 5; ++i) write_array(out, cent[i], z->n_cent/5, h->S_r);
    double (*fy[])[NKNOTS] = {m->k_pc,m->k_cp,m->k_cc};
    for (int i = 0; i < 3; ++i)
        for (size_t u = 0; u < z->n_cent/5; ++u)
            write_array(out, fy[i][u], h->nknots, h->S);
    const double *w[] = {m->f,m->g,m->hV,m->hC,m->hP4,m->hF4,m->hY4,m->hT4};
    const size_t n[] = {z->n_f,z->n_g,z->n_hV,z->n_hC,z->n_hP4,z->n_hF4,z->n_hY4,z->n_hT4};
    for (int i = 0; i < 8; ++i) write_array(out, w[i], n[i], h->S);
    require(fflush(out) == 0 && ftruncate(fileno(out), z->total_bytes) == 0 && fclose(out) == 0, "sizing envelope section");
    printf("staging=%s bytes=%zu header=%zu analytic=(%" PRId64 ",%" PRId64 ",%" PRId64 ")\n",
           file,z->total_bytes,sizeof(*h),h->f_p0,h->f_p1,h->f_p2);
}

int main(int argc, char **argv) {
    require(argc == 4, "usage: export {upper|lower} {full|small|fixture} FILE.part");
    int upper = !strcmp(argv[1],"upper"), small = !strcmp(argv[2],"small"), fixture = !strcmp(argv[2],"fixture");
    require(upper || !strcmp(argv[1],"lower"), "mode");
    require(small || fixture || !strcmp(argv[2],"full"), "geometry");
    Model m = {0};
    m.mode = upper ? 1 : 2;
    Cert4vHeader h = {.magic=CERT4V_MAGIC,.mode=m.mode,.P=upper?5658147:5656423,.Q=10000000,
        .S=weight_unit,.S_r=center_unit,.alpha_lo_num=upper?0:center_unit/4,
        .DF=fixture?12:1400,.D0_2=fixture?8:72,.KT_2=fixture?3:24,
        .D0_3=fixture?6:(small?14:28),.KT_3=fixture?3:12,
        .DHU=fixture?6:(small?12:24),.KHU=fixture?4:12,.DCU=fixture?6:(small?10:22),.KCU=fixture?4:10,
        .D0_4=fixture?5:(small?10:14),.KT_4=fixture?3:8,
        .DP4_1=fixture?5:(small?10:11),.DP4_2=fixture?5:(small?10:11),.DP4_S=fixture?9:16,
        .KP4_1=fixture?3:7,.KP4_2=fixture?3:7,
        .DF4_1=fixture?5:10,.DF4_2=fixture?5:10,.DF4_S=fixture?9:15,.KF4_1=fixture?3:6,.KF4_2=fixture?3:6,
        .DY4=fixture?5:10,.KY4=fixture?3:6,.DT4=fixture?5:8,.KT4=fixture?3:5,
        .Q_K=fixture?8:10080,.nknots=fixture?3:NKNOTS};
    int grid[NKNOTS] = {0,140,210,315,420,560,720,1008,1440,2016,2520,3360,4200,4620,5040};
    if (fixture) { h.knot_num[1]=1; h.knot_num[2]=4; }
    else memcpy(h.knot_num,grid,sizeof(grid));
    h.f_p0=h.f_p1=1;
    require(!c4_validate_header(&h),"initial geometry");
    Cert4vLayout setup;
    c4_init_layout(&setup,&h);
    Cert4vSizes z;
    c4_compute_sizes(&setup,&z);
    size_t nf = setup.L2.INNER_N;
    c4_free_layout(&setup);
    m.rpar=calloc(nf,sizeof(double)); m.rch=calloc(nf,sizeof(double));
    m.xpar=calloc(nf,sizeof(double)); m.xch=calloc(nf,sizeof(double)); m.xcc=calloc(nf,sizeof(double));
    /* Fixture uses an abbreviated knot grid and zero arrays. */
    m.k_pc=calloc(nf,sizeof(*m.k_pc)); m.k_cp=calloc(nf,sizeof(*m.k_cp)); m.k_cc=calloc(nf,sizeof(*m.k_cc));
    m.f=calloc(z.n_f,sizeof(double)); m.g=calloc(z.n_g,sizeof(double));
    m.hV=calloc(z.n_hV+1,sizeof(double)); m.hC=calloc(z.n_hC+1,sizeof(double));
    m.hP4=calloc(z.n_hP4+1,sizeof(double)); m.hF4=calloc(z.n_hF4+1,sizeof(double));
    m.hY4=calloc(z.n_hY4+1,sizeof(double)); m.hT4=calloc(z.n_hT4+1,sizeof(double));
    for (size_t i=0;i<z.n_f;i++) m.f[i]=1.;
    if (fixture) {
        m.A_up=1.; m.B_up=0.0625; m.B1_lo=1.; m.B2_lo=0.0625;
    } else {
        init_layout(&m.C,h.DF,h.D0_2,h.D0_3,h.DHU,h.DCU,h.D0_4);
        for (int i=0;i<NKNOTS;i++) m.KNOTS[i]=(double)h.knot_num[i]/h.Q_K;
        m.emp_edge=calloc(z.n_g,sizeof(double));
        run_empirical(&m,500000,30);
        /* Work close to the normalized root. Signs and stopping criteria have no proof role. */
        double l0 = upper?0.5657900:0.5656500;
        double rho=eval_lambda(&m,l0,15,7);
        double l1=l0+0.5*(rho-1.);
        double rho1=eval_lambda(&m,l1,8,2);
        double root=l1;
        if (fabs(rho1-rho)>1e-12) root=l1-(rho1-1.)*(l1-l0)/(rho1-rho);
        printf("proposal root=%.12f, verification P/Q=%d/%d\n",root,h.P,h.Q);
        eval_lambda(&m,root,14,0);
    }
    write_snapshot(&m,&h,&z,argv[3]);
    return 0;
}
