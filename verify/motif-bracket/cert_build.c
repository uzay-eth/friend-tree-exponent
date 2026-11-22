/* Exact-Rational Certificate Builder for the Random Friend Tree Exponent
 *
 * Builds deterministic quantized exact-rational Analytic 1/D-Tail 3-Vertex + Cherry
 * drift certificates (f + g + h_V + h_C) for:
 *   - mode 1 (UPPER): proves E[C_n] <= A_up * prod_{j=3}^{n-1} (1 + (P/Q)/j)
 *     Default: P/Q = 1133/2000 = 0.5665 < 17/30 = 0.566666... (REFUTES Strogatz 17/30!)
 *   - mode 2 (LOWER): proves E[C_n] >= A_lo * prod_{j=3}^{n-1} (1 + (P/Q)/j)
 *     Default: P/Q = 113/200 = 0.5650
 *
 * Output: compact binary certificate file (~135 MB) containing ONLY exact int64_t tables
 *         verified independently in arbitrary-precision GMP arithmetic by cert_verify_gmp.c.
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>
#include <omp.h>
#include "cert_format.h"

static uint64_t xs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t next64(void) { uint64_t r=rotl(xs[1]*5,7)*9, t=xs[1]<<17; xs[2]^=xs[0]; xs[3]^=xs[1]; xs[1]^=xs[2]; xs[0]^=xs[3]; xs[2]^=t; xs[3]=rotl(xs[3],45); return r; }
static uint64_t sm(uint64_t *x) { uint64_t z=(*x+=0x9e3779b97f4a7c15ULL); z=(z^(z>>30))*0xbf58476d1ce4e5b9ULL; z=(z^(z>>27))*0x94d049bb133111ebULL; return z^(z>>31); }
static inline uint32_t unif(uint32_t m) { uint64_t x=next64()>>32, prod=x*(uint64_t)m; uint32_t low=(uint32_t)prod; if(low<m){uint32_t thr=(uint32_t)(-m)%m;while(low<thr){x=next64()>>32;prod=x*(uint64_t)m;low=(uint32_t)prod;}} return (uint32_t)(prod>>32); }

static void run_empirical(const CertLayout *C, int NMAX, int RUNS,
                          double *emp_rpar, double *emp_rch, double *emp_edge) {
    int N = C->NG;
    double *cp = calloc(N, sizeof(double)), *cc = calloc(N, sizeof(double));
    uint64_t seed = 20261122ULL; for (int i=0; i<4; i++) xs[i]=sm(&seed);
    uint32_t *par=malloc(NMAX*sizeof(uint32_t)), *deg=malloc(NMAX*sizeof(uint32_t));
    uint32_t *kcnt=calloc(NMAX,sizeof(uint32_t)), *kcap=calloc(NMAX,sizeof(uint32_t));
    uint32_t **kids=calloc(NMAX,sizeof(uint32_t*)); int *st=malloc(NMAX*sizeof(int));
    for (int run=0; run<RUNS; run++) {
        for (int i=0; i<NMAX; i++) { kcnt[i]=0; if (kcap[i]>64) { free(kids[i]); kids[i]=0; kcap[i]=0; } }
        par[0]=UINT32_MAX; deg[0]=2; par[1]=0; deg[1]=1; par[2]=0; deg[2]=1;
        if (kcap[0]<2) { kids[0]=realloc(kids[0],2*sizeof(uint32_t)); kcap[0]=2; }
        kids[0][0]=1; kids[0][1]=2; kcnt[0]=2;
        for (uint32_t n=3; n<(uint32_t)NMAX; n++) {
            uint32_t V=unif(n), d=deg[V], r=unif(d);
            uint32_t W=(par[V]!=UINT32_MAX)?((r==0)?par[V]:kids[V][r-1]):kids[V][r];
            par[n]=W; deg[n]=1; deg[W]++;
            if (kcnt[W]==kcap[W]) { uint32_t nc=kcap[W]?2*kcap[W]:2; kids[W]=realloc(kids[W],nc*sizeof(uint32_t)); kcap[W]=nc; }
            kids[W][kcnt[W]++]=n;
        }
        for (int u=0; u<NMAX; u++) {
            if (deg[u]<2) { st[u]=-1; continue; }
            int L=0; for (uint32_t i=0; i<kcnt[u]; i++) if (deg[kids[u][i]]==1) L++;
            int h=(par[u]!=UINT32_MAX)?1:0, k=(int)deg[u]-L-h;
            st[u] = ((int)deg[u] < C->D0) ? cl_inner_idx(C,L,k,h) : cl_t1_idx(C,k,h);
        }
        for (int w=1; w<NMAX; w++) {
            if (deg[w]<2) continue;
            int p=(int)par[w], sp=st[p], sw=st[w];
            emp_edge[sp*N+sw] += 1.0;
            emp_rpar[sw] += 1.0/deg[p]; cp[sw] += 1.0;
            emp_rch[sp]  += 1.0/deg[w]; cc[sp] += 1.0;
        }
    }
    for (int u=0; u<C->INNER_N; u++) {
        int k=C->gK[u], L=C->gL[u];
        emp_rpar[u] = (cp[u]>=5.0) ? (emp_rpar[u]/cp[u]) : (C->gH[u] ? (0.25/(1+0.15*L+0.35*k)) : 0.0);
        emp_rch[u]  = (cc[u]>=5.0) ? (emp_rch[u]/cc[u])  : (k>0   ? (0.36/(1+0.04*L+0.08*k)) : 0.0);
    }
    for (int i=0; i<NMAX; i++) free(kids[i]);
    free(par); free(deg); free(kcnt); free(kcap); free(kids); free(st); free(cp); free(cc);
}

static inline double eval_g_p(const CertLayout *C, const double *g, int Lp, int kp, int hp, int iw) {
    int Dp = Lp + kp + hp;
    if (Dp < C->D0) return g[cl_inner_idx(C, Lp, kp, hp) * C->NG + iw];
    double tp = (double)C->D0 / Dp;
    return (1.0 - tp) * g[C->INF_ID * C->NG + iw] + tp * g[cl_t1_idx(C, kp, hp) * C->NG + iw];
}

static inline double eval_pl(const double *KNOTS, const double *vals, double x) {
    if (x <= KNOTS[0]) return vals[0];
    if (x >= KNOTS[CERT_NKNOTS - 1]) return vals[CERT_NKNOTS - 1];
    for (int i = 0; i < CERT_NKNOTS - 1; i++) {
        if (x <= KNOTS[i + 1]) {
            double t = (x - KNOTS[i]) / (KNOTS[i + 1] - KNOTS[i]);
            return vals[i] * (1.0 - t) + vals[i + 1] * t;
        }
    }
    return vals[CERT_NKNOTS - 1];
}

static void project_convex(const double *KNOTS, double *a, double sgn) {
    int m=CERT_NKNOTS-1, bs[CERT_NKNOTS], be[CERT_NKNOTS], nb=0; double sl[CERT_NKNOTS], wt[CERT_NKNOTS];
    for (int i=0; i<m; i++) {
        double dx=KNOTS[i+1]-KNOTS[i]; sl[nb]=sgn*(a[i+1]-a[i])/dx; wt[nb]=dx; bs[nb]=i; be[nb]=i; nb++;
        while (nb>=2 && sl[nb-2]>sl[nb-1]) {
            double nw=wt[nb-2]+wt[nb-1];
            sl[nb-2]=(sl[nb-2]*wt[nb-2]+sl[nb-1]*wt[nb-1])/nw;
            wt[nb-2]=nw; be[nb-2]=be[nb-1]; nb--;
        }
    }
    double cur=a[0];
    for (int b=0; b<nb; b++)
        for (int i=bs[b]; i<=be[b]; i++) { cur+=sgn*sl[b]*(KNOTS[i+1]-KNOTS[i]); a[i+1]=cur; }
}

static void opt_knots(const double *KNOTS, double *a, double sgn,
                      const double *xs_v, const double *wxs, int nx,
                      const double *ys_v, const double *wys, int ny) {
    double ym=0, wyt=0, wxt=0;
    for (int j=0; j<ny; j++) { ym+=wys[j]*ys_v[j]; wyt+=wys[j]; }
    for (int i=0; i<nx; i++) wxt+=wxs[i];
    if (wyt<=0||wxt<=0) { for(int k=0;k<CERT_NKNOTS;k++) a[k]=0; return; }
    ym/=wyt; for (int k=0; k<CERT_NKNOTS; k++) a[k]=ym*KNOTS[k];
    double mu[CERT_NKNOTS]={0};
    for (int i=0; i<nx; i++) {
        double x=xs_v[i], w=wxs[i]/wxt;
        if (x<=KNOTS[0]) { mu[0]+=w; continue; }
        if (x>=KNOTS[CERT_NKNOTS-1]) { mu[CERT_NKNOTS-1]+=w; continue; }
        for (int k=0; k<CERT_NKNOTS-1; k++)
            if (x<=KNOTS[k+1]) { double t=(x-KNOTS[k])/(KNOTS[k+1]-KNOTS[k]); mu[k]+=w*(1-t); mu[k+1]+=w*t; break; }
    }
    double yvar=0; for (int j=0; j<ny; j++) { double d=ys_v[j]-ym; yvar+=(wys[j]/wyt)*d*d; }
    double ysd=sqrt(yvar+1e-12), best_a[CERT_NKNOTS], best_obj=1e30;
    memcpy(best_a,a,sizeof(best_a));
    for (int step=0; step<40; step++) {
        double grad[CERT_NKNOTS], obj=0;
        for (int k=0; k<CERT_NKNOTS; k++) { grad[k]=sgn*mu[k]; obj+=sgn*mu[k]*a[k]; }
        for (int j=0; j<ny; j++) {
            double y=ys_v[j], w=wys[j]/wyt, bv=-1e30; int bk=0;
            for (int k=0; k<CERT_NKNOTS; k++) { double c=sgn*(KNOTS[k]*y-a[k]); if (c>bv) { bv=c; bk=k; } }
            obj+=w*bv; grad[bk]-=sgn*w;
        }
        if (obj<best_obj) { best_obj=obj; memcpy(best_a,a,sizeof(best_a)); }
        double lr=0.25*ysd/sqrt(1.0+0.5*step);
        for (int k=0; k<CERT_NKNOTS; k++) a[k]-=lr*grad[k];
        double c0=a[0]; for (int k=0; k<CERT_NKNOTS; k++) a[k]-=c0;
        project_convex(KNOTS, a, sgn);
    }
    memcpy(a, best_a, sizeof(best_a));
}

static inline double gamma_min_t1(int D0, int KT, int kap, int h) {
    if (kap > KT) return 0.0;
    int m = D0 - kap - h;
    return (m > 0) ? ((double)m / (D0 + 1.0)) : 0.0;
}
static inline double b_min_t1(int D0, int KT, int kap, int h) {
    if (kap > KT) return 0.0;
    int m = D0 - kap - h;
    return (m > 0) ? ((double)m / (double)D0) : 0.0;
}

static void solve_g_hv_hc(
    const CertLayout *C, const double *KNOTS, int mode, double lam_e,
    double rew_g, double rew_h,
    const double *gInvD, const double *rpar, const double *rch,
    const double *xpar, const double *xch, const double *xcc,
    const double *dLf_g, const double *A_in, const double *A_out,
    double (*k_pc)[CERT_NKNOTS], double (*k_cp)[CERT_NKNOTS], double (*k_cc)[CERT_NKNOTS],
    int **p_by_r, const int *p_rcnt, int **c_by_r, const int *c_rcnt,
    const double *hV_old, const double *hC_old,
    double *gn, double *hVn, double *hCn)
{
    int N = C->NG, D0 = C->D0, KT = C->KT, max_rank = D0 + 1;
    int iw_nb = cl_inner_idx(C, 1, 0, 1);
    double sgn = (mode == 1) ? 1.0 : -1.0;

    /* 1. Solve g in exact DAG order */
    memset(gn, 0, (size_t)N * N * sizeof(double));
    for (int rsum = 2 * max_rank; rsum >= 4; rsum--) {
        int rp_lo = rsum - max_rank; if (rp_lo < 2) rp_lo = 2;
        int rp_hi = rsum - 2; if (rp_hi > max_rank) rp_hi = max_rank;
        #pragma omp parallel for schedule(dynamic) num_threads(3)
        for (int rp = rp_lo; rp <= rp_hi; rp++) {
            int rw = rsum - rp;
            if (rp >= D0 && rw >= D0) continue;
            for (int ip = 0; ip < p_rcnt[rp]; ip++) {
                int p = p_by_r[rp][ip], pk = C->gKind[p], Lp = C->gL[p], kp = C->gK[p], hp = C->gH[p];
                int p_up = (pk == 0) ? ((C->gD[p]+1 < D0) ? cl_inner_idx(C,Lp+1,kp,hp) : cl_t1_idx(C,kp,hp)) : C->INF_ID;
                int p_bud = (pk == 0 && Lp>=1) ? cl_inner_idx(C,Lp-1,kp+1,hp) : ((pk == 2 && kp <= KT) ? cl_t1_idx(C,kp+1,hp) : p);
                double bp_lo = (pk == 0) ? ((Lp>=1)?(double)Lp/C->gD[p]:0.0) : ((pk==2 && kp<=KT) ? b_min_t1(D0,KT,kp,hp) : 0.0);
                double bp_hi = (pk == 0) ? bp_lo : ((pk==2 && kp<=KT) ? 1.0 : 0.0);
                double gp_lo = (pk == 2) ? gamma_min_t1(D0,KT,kp,hp) : 0.0, gp_hi = (pk == 2) ? 1.0 : 0.0;
                int p_in_hu = (pk == 0 && cl_hu_idx(C,Lp,kp,hp)>=0), p_in_cu = (pk == 0 && cl_cu_idx(C,Lp,kp,hp)>=0);
                double inv_p = gInvD[p];

                for (int iw = 0; iw < c_rcnt[rw]; iw++) {
                    int w = c_by_r[rw][iw], wk = C->gKind[w], Lw = C->gL[w], kw = C->gK[w];
                    int w_up = (wk == 0) ? ((C->gD[w]+1 < D0) ? cl_inner_idx(C,Lw+1,kw,1) : cl_t1_idx(C,kw,1)) : C->INF_ID;
                    int w_bud = (wk == 0 && Lw>=1) ? cl_inner_idx(C,Lw-1,kw+1,1) : ((wk == 2 && kw <= KT) ? cl_t1_idx(C,kw+1,1) : w);
                    double bw_lo = (wk == 0) ? ((Lw>=1)?(double)Lw/C->gD[w]:0.0) : ((wk==2 && kw<=KT) ? b_min_t1(D0,KT,kw,1) : 0.0);
                    double bw_hi = (wk == 0) ? bw_lo : ((wk==2 && kw<=KT) ? 1.0 : 0.0);
                    double gw_lo = (wk == 2) ? gamma_min_t1(D0,KT,kw,1) : 0.0, gw_hi = (wk == 2) ? 1.0 : 0.0;
                    int w_in_hu = (wk == 0 && cl_hu_idx(C,Lw,kw,1)>=0);
                    double inv_w = gInvD[w];

                    double a_const = A_out[p*N + w] + A_in[p*N + w];
                    double gp_plus = gn[p_up*N + w], gw_plus = gn[p*N + w_up];
                    double gp_bud_v = gn[p_bud*N + w], gw_bud_v = gn[p*N + w_bud];

                    int cu_pb = (pk == 0 && Lp >= 1) ? cl_cu_idx(C, Lp-1, kp+1, hp) : -1;
                    double p_hcnb = (cu_pb >= 0) ? hC_old[cl_hc_off(C, cu_pb, w, iw_nb)] : 0.0;
                    int hu_wb = (wk == 0 && Lw >= 1) ? cl_hu_idx(C, Lw-1, kw+1, 1) : -1;
                    double w_hvnb = (hu_wb >= 0) ? hV_old[cl_hv_off(C, p, hu_wb, iw_nb)] : 0.0;

                    int has_pc = (pk == 0 && hp == 1 && !p_in_hu);
                    int has_cc = (pk == 0 && kp >= 2 && !p_in_cu);
                    int has_cp = (wk == 0 && kw >= 1 && !w_in_hu);
                    double d0_w = (pk == 2) ? (gn[C->INF_ID*N + w_up] - gn[C->INF_ID*N + w]) : 0.0;
                    double d0_p = (wk == 2) ? (gn[p_up*N + C->INF_ID] - gn[p*N + C->INF_ID]) : 0.0;

                    double num0 = a_const + sgn * rew_g, den0 = lam_e;
                    if (pk == 0) {
                        num0 += (inv_w - rch[p]) * dLf_g[p];
                        double ur = Lp + (wk == 0 ? inv_w : 0.0) + (p_in_hu ? xpar[p] : 0.0) + (p_in_cu ? (kp-1)*xcc[p] : 0.0);
                        num0 += ur * gp_plus; den0 += ur;
                        if (bp_lo > 0) { num0 += bp_lo * (gp_bud_v + p_hcnb); den0 += bp_lo; }
                    }
                    if (wk == 0) {
                        num0 += (inv_p - rpar[w]) * dLf_g[w];
                        double ur = Lw + (pk == 0 ? inv_p : 0.0) + (w_in_hu ? kw*xch[w] : 0.0);
                        num0 += ur * gw_plus; den0 += ur;
                        if (bw_lo > 0) { num0 += bw_lo * (gw_bud_v + w_hvnb); den0 += bw_lo; }
                    }
                    double v = num0 / den0;
                    if (has_pc || has_cc || has_cp || pk == 2 || wk == 2) {
                        int i1p=-1, i2p=-1, i3p=-1, mask_p=-1;
                        for (int ns = 0; ns < 16; ns++) {
                            double num = num0, den = den0; int i1=0, i2=0, i3=0, mask=0;
                            if (pk == 2) {
                                double gr = ((mode==1 && gp_plus > v) || (mode==2 && gp_plus < v)) ? (mask|=1, gp_hi) : gp_lo;
                                num += gr * gp_plus; den += gr;
                                if (bp_hi > 0) { double br = ((mode==1 && gp_bud_v > v) || (mode==2 && gp_bud_v < v)) ? (mask|=2, bp_hi) : bp_lo; num += br * gp_bud_v; den += br; }
                                double d1 = gw_plus - v;
                                if ((mode==1 && d1 > d0_w) || (mode==2 && d1 < d0_w)) { mask |= 4; num += inv_p * gw_plus; den += inv_p; }
                                else num += inv_p * d0_w;
                            }
                            if (wk == 2) {
                                double gr = ((mode==1 && gw_plus > v) || (mode==2 && gw_plus < v)) ? (mask|=8, gw_hi) : gw_lo;
                                num += gr * gw_plus; den += gr;
                                if (bw_hi > 0) { double br = ((mode==1 && gw_bud_v > v) || (mode==2 && gw_bud_v < v)) ? (mask|=16, bw_hi) : bw_lo; num += br * gw_bud_v; den += br; }
                                double d1 = gp_plus - v;
                                if ((mode==1 && d1 > d0_p) || (mode==2 && d1 < d0_p)) { mask |= 32; num += inv_w * gp_plus; den += inv_w; }
                                else num += inv_w * d0_p;
                            }
                            if (has_pc) { double yp=gp_plus-v, bv=(mode==1)?-1e30:1e30; for (int kk=0;kk<CERT_NKNOTS;kk++){double c=KNOTS[kk]*yp-k_pc[p][kk]; if((mode==1&&c>bv)||(mode==2&&c<bv)){bv=c;i1=kk;}} num+=KNOTS[i1]*gp_plus-k_pc[p][i1]; den+=KNOTS[i1]; }
                            if (has_cc) { double yp=gp_plus-v, bv=(mode==1)?-1e30:1e30; for (int kk=0;kk<CERT_NKNOTS;kk++){double c=KNOTS[kk]*yp-k_cc[p][kk]; if((mode==1&&c>bv)||(mode==2&&c<bv)){bv=c;i2=kk;}} num+=(kp-1)*(KNOTS[i2]*gp_plus-k_cc[p][i2]); den+=(kp-1)*KNOTS[i2]; }
                            if (has_cp) { double yw=gw_plus-v, bv=(mode==1)?-1e30:1e30; for (int kk=0;kk<CERT_NKNOTS;kk++){double c=KNOTS[kk]*yw-k_cp[w][kk]; if((mode==1&&c>bv)||(mode==2&&c<bv)){bv=c;i3=kk;}} num+=kw*(KNOTS[i3]*gw_plus-k_cp[w][i3]); den+=kw*KNOTS[i3]; }
                            double vn = num / den;
                            if (i1==i1p && i2==i2p && i3==i3p && mask==mask_p) { v=vn; break; }
                            i1p=i1; i2p=i2; i3p=i3; mask_p=mask; v=vn;
                        }
                    }
                    gn[p*N + w] = v;
                }
            }
        }
    }

    /* 2. Solve h_V in exact DAG order */
    size_t hv_len = (size_t)C->npars * C->NHU_ACT * C->nchs;
    memset(hVn, 0, hv_len * sizeof(double));
    int ru_max = C->DHU - 1;
    for (int rsum = 2 * max_rank + ru_max; rsum >= 6; rsum--) {
        int tasks[16384][3], nt = 0;
        for (int ru = 2; ru <= ru_max; ru++)
            for (int rp = 2; rp <= max_rank; rp++) {
                int rw = rsum - ru - rp;
                if (rw >= 2 && rw <= max_rank && p_rcnt[rp] > 0 && c_rcnt[rw] > 0) {
                    tasks[nt][0]=rp; tasks[nt][1]=ru; tasks[nt][2]=rw; nt++;
                }
            }
        #pragma omp parallel for schedule(dynamic) num_threads(3)
        for (int ti = 0; ti < nt; ti++) {
            int rp = tasks[ti][0], ru = tasks[ti][1], rw = tasks[ti][2];
            for (int ku = C->KHU; ku >= 1; ku--) {
                int Lu = ru - 1 - ku; if (Lu < 0 || Lu >= C->DHU) continue;
                int ui = cl_hu_idx(C, Lu, ku, 1); if (ui < 0) continue;
                int gu = cl_inner_idx(C, Lu, ku, 1), gup = (ru + 1 < D0) ? cl_inner_idx(C, Lu+1, ku, 1) : cl_t1_idx(C, ku, 1);
                double iu = 1.0/ru, bu = (double)Lu/ru, ru_hi = 0.5*(ku - 1), xp_u = xpar[gu], xw_u = xch[gu];
                int ui_up = cl_hu_idx(C, Lu+1, ku, 1), ui_bud = (Lu>=1) ? cl_hu_idx(C, Lu-1, ku+1, 1) : -1;

                for (int ip = 0; ip < p_rcnt[rp]; ip++) {
                    int p = p_by_r[rp][ip], pk = C->gKind[p], Lp = C->gL[p], kp = C->gK[p], hp = C->gH[p];
                    int p_up = (pk == 0) ? ((C->gD[p]+1 < D0) ? cl_inner_idx(C,Lp+1,kp,hp) : cl_t1_idx(C,kp,hp)) : C->INF_ID;
                    int p_bud = (pk == 0 && Lp>=1) ? cl_inner_idx(C,Lp-1,kp+1,hp) : ((pk == 2 && kp <= KT) ? cl_t1_idx(C,kp+1,hp) : p);
                    double bp_lo = (pk == 0) ? ((Lp>=1)?(double)Lp/C->gD[p]:0.0) : ((pk==2 && kp<=KT) ? b_min_t1(D0,KT,kp,hp) : 0.0);
                    double bp_hi = (pk == 0) ? bp_lo : ((pk==2 && kp<=KT) ? 1.0 : 0.0);
                    double gp_lo = (pk == 2) ? gamma_min_t1(D0,KT,kp,hp) : 0.0, gp_hi = (pk == 2) ? 1.0 : 0.0;
                    double rp_hi = (pk == 0) ? 0.5*(hp + kp - 1) : 0.0, inv_p = gInvD[p], dgpu = gn[p*N + gup] - gn[p*N + gu];

                    for (int iw = 0; iw < c_rcnt[rw]; iw++) {
                        int w = c_by_r[rw][iw], wk = C->gKind[w], Lw = C->gL[w], kw = C->gK[w];
                        int w_up = (wk == 0) ? ((C->gD[w]+1 < D0) ? cl_inner_idx(C,Lw+1,kw,1) : cl_t1_idx(C,kw,1)) : C->INF_ID;
                        int w_bud = (wk == 0 && Lw>=1) ? cl_inner_idx(C,Lw-1,kw+1,1) : ((wk == 2 && kw <= KT) ? cl_t1_idx(C,kw+1,1) : w);
                        double bw_lo = (wk == 0) ? ((Lw>=1)?(double)Lw/C->gD[w]:0.0) : ((wk==2 && kw<=KT) ? b_min_t1(D0,KT,kw,1) : 0.0);
                        double bw_hi = (wk == 0) ? bw_lo : ((wk==2 && kw<=KT) ? 1.0 : 0.0);
                        double gw_lo = (wk == 2) ? gamma_min_t1(D0,KT,kw,1) : 0.0, gw_hi = (wk == 2) ? 1.0 : 0.0;
                        double rw_hi = (wk == 0) ? 0.5*kw : 0.0, inv_w = gInvD[w], dguw = gn[gup*N + w] - gn[gu*N + w];

                        double hp_up = hVn[cl_hv_off(C, p_up, ui, w)];
                        double hu_up = (ui_up>=0) ? hVn[cl_hv_off(C, p, ui_up, w)] : 0.0;
                        double hw_up = hVn[cl_hv_off(C, p, ui, w_up)];
                        double hp_bud = (bp_hi>0) ? hVn[cl_hv_off(C, p_bud, ui, w)] : 0.0;
                        double hu_bud = (ui_bud>=0) ? hVn[cl_hv_off(C, p, ui_bud, w)] : 0.0;
                        double hw_bud = (bw_hi>0) ? hVn[cl_hv_off(C, p, ui, w_bud)] : 0.0;

                        double num0 = (inv_p - xp_u)*dguw + (inv_w - xw_u)*dgpu + sgn * rew_h, den0 = lam_e;
                        if (pk == 0) { double r = Lp + iu; num0 += r*hp_up; den0 += r; if (bp_lo>0) { num0 += bp_lo*hp_bud; den0 += bp_lo; } }
                        if (wk == 0) { double r = Lw + iu; num0 += r*hw_up; den0 += r; if (bw_lo>0) { num0 += bw_lo*hw_bud; den0 += bw_lo; } }
                        { double ru_b = Lu + (pk == 0 ? inv_p : 0.0) + (wk == 0 ? inv_w : 0.0); num0 += ru_b*hu_up; den0 += ru_b; if (bu > 0) { num0 += bu*hu_bud; den0 += bu; } }

                        double d0_pu = (pk == 2) ? (((ui_up>=0)?hVn[cl_hv_off(C,C->INF_ID,ui_up,w)]:0.0) - hVn[cl_hv_off(C,C->INF_ID,ui,w)]) : 0.0;
                        double d0_wu = (wk == 2) ? (((ui_up>=0)?hVn[cl_hv_off(C,p,ui_up,C->INF_ID)]:0.0) - hVn[cl_hv_off(C,p,ui,C->INF_ID)]) : 0.0;

                        double v = num0 / den0; int mask_p = -1;
                        for (int ns = 0; ns < 12; ns++) {
                            double num = num0, den = den0; int mask = 0;
                            if (rp_hi > 0 && ((mode==1 && hp_up>v)||(mode==2 && hp_up<v))) { num += rp_hi*hp_up; den += rp_hi; mask|=1; }
                            if (ru_hi > 0 && ((mode==1 && hu_up>v)||(mode==2 && hu_up<v))) { num += ru_hi*hu_up; den += ru_hi; mask|=2; }
                            if (rw_hi > 0 && ((mode==1 && hw_up>v)||(mode==2 && hw_up<v))) { num += rw_hi*hw_up; den += rw_hi; mask|=4; }
                            if (pk == 2) {
                                double gr = ((mode==1 && hp_up>v)||(mode==2 && hp_up<v)) ? (mask|=8, gp_hi) : gp_lo; num += gr*hp_up; den += gr;
                                if (bp_hi > 0) { double br = ((mode==1 && hp_bud>v)||(mode==2 && hp_bud<v)) ? (mask|=16, bp_hi) : bp_lo; num += br*hp_bud; den += br; }
                                double d1 = hu_up - v;
                                if ((mode==1 && d1 > d0_pu) || (mode==2 && d1 < d0_pu)) { mask|=32; num += inv_p*hu_up; den += inv_p; } else num += inv_p*d0_pu;
                            }
                            if (wk == 2) {
                                double gr = ((mode==1 && hw_up>v)||(mode==2 && hw_up<v)) ? (mask|=64, gw_hi) : gw_lo; num += gr*hw_up; den += gr;
                                if (bw_hi > 0) { double br = ((mode==1 && hw_bud>v)||(mode==2 && hw_bud<v)) ? (mask|=128, bw_hi) : bw_lo; num += br*hw_bud; den += br; }
                                double d1 = hu_up - v;
                                if ((mode==1 && d1 > d0_wu) || (mode==2 && d1 < d0_wu)) { mask|=256; num += inv_w*hu_up; den += inv_w; } else num += inv_w*d0_wu;
                            }
                            double vn = num / den; if (mask == mask_p) { v = vn; break; } mask_p = mask; v = vn;
                        }
                        hVn[cl_hv_off(C, p, ui, w)] = v;
                    }
                }
            }
        }
    }

    /* 3. Solve h_C in exact DAG order */
    size_t hc_len = (size_t)C->NCU_ACT * C->nchs * (C->nchs + 1) / 2;
    memset(hCn, 0, hc_len * sizeof(double));
    int rcu_max = C->DCU - 1;
    for (int rsum = rcu_max + 2 * max_rank; rsum >= 6; rsum--) {
        int tasks[16384][3], nt = 0;
        for (int ru = 2; ru <= rcu_max; ru++)
            for (int r1 = 2; r1 <= max_rank; r1++) {
                int r2 = rsum - ru - r1;
                if (r2 >= r1 && r2 <= max_rank && c_rcnt[r1] > 0 && c_rcnt[r2] > 0) {
                    tasks[nt][0]=ru; tasks[nt][1]=r1; tasks[nt][2]=r2; nt++;
                }
            }
        #pragma omp parallel for schedule(dynamic) num_threads(3)
        for (int ti = 0; ti < nt; ti++) {
            int ru=tasks[ti][0], r1=tasks[ti][1], r2=tasks[ti][2];
            for (int hu=1; hu>=0; hu--) for (int ku=C->KCU; ku>=2; ku--) {
                int Lu = ru - hu - ku; if (Lu < 0 || Lu >= C->DCU) continue;
                int ci = cl_cu_idx(C, Lu, ku, hu); if (ci < 0) continue;
                int gu = cl_inner_idx(C, Lu, ku, hu), gup = (ru + 1 < D0) ? cl_inner_idx(C, Lu+1, ku, hu) : cl_t1_idx(C, ku, hu);
                double iu = 1.0/ru, bu = (double)Lu/ru, ru_hi = 0.5*(hu + ku - 2), xc_u = xcc[gu];
                int ci_up = cl_cu_idx(C, Lu+1, ku, hu), ci_bud = (Lu>=1) ? cl_cu_idx(C, Lu-1, ku+1, hu) : -1;

                for (int i1 = 0; i1 < c_rcnt[r1]; i1++) {
                    int w1 = c_by_r[r1][i1], k1_kind = C->gKind[w1], L1 = C->gL[w1], k1 = C->gK[w1];
                    int w1_up = (k1_kind==0) ? ((C->gD[w1]+1<D0)?cl_inner_idx(C,L1+1,k1,1):cl_t1_idx(C,k1,1)) : C->INF_ID;
                    int w1_bud = (k1_kind==0 && L1>=1) ? cl_inner_idx(C,L1-1,k1+1,1) : ((k1_kind==2 && k1<=KT)?cl_t1_idx(C,k1+1,1):w1);
                    double b1_lo = (k1_kind==0) ? ((L1>=1)?(double)L1/C->gD[w1]:0.0) : ((k1_kind==2&&k1<=KT)?b_min_t1(D0,KT,k1,1):0.0);
                    double b1_hi = (k1_kind==0) ? b1_lo : ((k1_kind==2&&k1<=KT)?1.0:0.0);
                    double g1_lo = (k1_kind==2) ? gamma_min_t1(D0,KT,k1,1) : 0.0, g1_hi = (k1_kind==2) ? 1.0 : 0.0;
                    double r1_hi = (k1_kind==0) ? 0.5*k1 : 0.0, inv1 = gInvD[w1], dg1 = gn[gup*N + w1] - gn[gu*N + w1];

                    int j_start = (r1 == r2) ? i1 : 0;
                    for (int i2 = j_start; i2 < c_rcnt[r2]; i2++) {
                        int w2 = c_by_r[r2][i2], k2_kind = C->gKind[w2], L2 = C->gL[w2], k2 = C->gK[w2];
                        int w2_up = (k2_kind==0) ? ((C->gD[w2]+1<D0)?cl_inner_idx(C,L2+1,k2,1):cl_t1_idx(C,k2,1)) : C->INF_ID;
                        int w2_bud = (k2_kind==0 && L2>=1) ? cl_inner_idx(C,L2-1,k2+1,1) : ((k2_kind==2 && k2<=KT)?cl_t1_idx(C,k2+1,1):w2);
                        double b2_lo = (k2_kind==0) ? ((L2>=1)?(double)L2/C->gD[w2]:0.0) : ((k2_kind==2&&k2<=KT)?b_min_t1(D0,KT,k2,1):0.0);
                        double b2_hi = (k2_kind==0) ? b2_lo : ((k2_kind==2&&k2<=KT)?1.0:0.0);
                        double g2_lo = (k2_kind==2) ? gamma_min_t1(D0,KT,k2,1) : 0.0, g2_hi = (k2_kind==2) ? 1.0 : 0.0;
                        double r2_hi = (k2_kind==0) ? 0.5*k2 : 0.0, inv2 = gInvD[w2], dg2 = gn[gup*N + w2] - gn[gu*N + w2];

                        double cu_up  = (ci_up>=0)  ? hCn[cl_hc_off(C, ci_up, w1, w2)] : 0.0;
                        double cu_bud = (ci_bud>=0) ? hCn[cl_hc_off(C, ci_bud, w1, w2)] : 0.0;
                        double c1_up  = hCn[cl_hc_off(C, ci, w1_up, w2)];
                        double c2_up  = hCn[cl_hc_off(C, ci, w1, w2_up)];
                        double c1_bud = (b1_hi>0) ? hCn[cl_hc_off(C, ci, w1_bud, w2)] : 0.0;
                        double c2_bud = (b2_hi>0) ? hCn[cl_hc_off(C, ci, w1, w2_bud)] : 0.0;

                        double num0 = (inv1 - xc_u)*dg2 + (inv2 - xc_u)*dg1 + sgn * rew_h, den0 = lam_e;
                        if (k1_kind==0) { double r=L1+iu; num0+=r*c1_up; den0+=r; if (b1_lo>0){num0+=b1_lo*c1_bud;den0+=b1_lo;} }
                        if (k2_kind==0) { double r=L2+iu; num0+=r*c2_up; den0+=r; if (b2_lo>0){num0+=b2_lo*c2_bud;den0+=b2_lo;} }
                        { double ru_b = Lu + (k1_kind==0 ? inv1 : 0.0) + (k2_kind==0 ? inv2 : 0.0); num0 += ru_b*cu_up; den0 += ru_b; if (bu>0) { num0 += bu*cu_bud; den0 += bu; } }

                        double d0_1 = (k1_kind==2) ? (((ci_up>=0)?hCn[cl_hc_off(C,ci_up,C->INF_ID,w2)]:0.0) - hCn[cl_hc_off(C,ci,C->INF_ID,w2)]) : 0.0;
                        double d0_2 = (k2_kind==2) ? (((ci_up>=0)?hCn[cl_hc_off(C,ci_up,w1,C->INF_ID)]:0.0) - hCn[cl_hc_off(C,ci,w1,C->INF_ID)]) : 0.0;

                        double v = num0 / den0; int mask_p = -1;
                        for (int ns=0; ns<12; ns++) {
                            double num=num0, den=den0; int mask=0;
                            if (ru_hi>0 && ((mode==1&&cu_up>v)||(mode==2&&cu_up<v))) { num+=ru_hi*cu_up; den+=ru_hi; mask|=1; }
                            if (r1_hi>0 && ((mode==1&&c1_up>v)||(mode==2&&c1_up<v))) { num+=r1_hi*c1_up; den+=r1_hi; mask|=2; }
                            if (r2_hi>0 && ((mode==1&&c2_up>v)||(mode==2&&c2_up<v))) { num+=r2_hi*c2_up; den+=r2_hi; mask|=4; }
                            if (k1_kind==2) {
                                double gr = ((mode==1&&c1_up>v)||(mode==2&&c1_up<v)) ? (mask|=8, g1_hi) : g1_lo; num+=gr*c1_up; den+=gr;
                                if (b1_hi>0) { double br=((mode==1&&c1_bud>v)||(mode==2&&c1_bud<v))?(mask|=16,b1_hi):b1_lo; num+=br*c1_bud; den+=br; }
                                double d1 = cu_up - v;
                                if ((mode==1 && d1>d0_1)||(mode==2 && d1<d0_1)) { mask|=32; num+=inv1*cu_up; den+=inv1; } else num += inv1*d0_1;
                            }
                            if (k2_kind==2) {
                                double gr = ((mode==1&&c2_up>v)||(mode==2&&c2_up<v)) ? (mask|=64, g2_hi) : g2_lo; num+=gr*c2_up; den+=gr;
                                if (b2_hi>0) { double br=((mode==1&&c2_bud>v)||(mode==2&&c2_bud<v))?(mask|=128,b2_hi):b2_lo; num+=br*c2_bud; den+=br; }
                                double d1 = cu_up - v;
                                if ((mode==1 && d1>d0_2)||(mode==2 && d1<d0_2)) { mask|=256; num+=inv2*cu_up; den+=inv2; } else num += inv2*d0_2;
                            }
                            double vn = num / den; if (mask==mask_p) { v=vn; break; } mask_p=mask; v=vn;
                        }
                        hCn[cl_hc_off(C, ci, w1, w2)] = v;
                    }
                }
            }
        }
    }
}

static void build_cert(int mode, int P, int Q, const char *out_path) {
    CertHeader H;
    memset(&H, 0, sizeof(H));
    H.magic = CERT_MAGIC;
    H.mode = mode; H.P = P; H.Q = Q;
    H.S = 10000000000LL; /* 10^10 */
    H.S_r = 1000000LL;   /* 10^6 */
    H.LF = 800; H.KF = 150; H.D0 = 24; H.KT = 12;
    H.DHU = 16; H.KHU = 8; H.DCU = 12; H.KCU = 6;
    H.Q_K = 240; H.nknots = CERT_NKNOTS;
    int kn[CERT_NKNOTS] = {0, 10, 15, 25, 40, 60, 80, 100, 120};
    for (int i = 0; i < CERT_NKNOTS; i++) H.knot_num[i] = kn[i];

    CertLayout C;
    cl_init(&C, &H);
    int N = C.NG, v_nb = cl_vidx(&C, 1, 0, 1), iw_nb = cl_inner_idx(&C, 1, 0, 1);
    double KNOTS[CERT_NKNOTS];
    for (int i = 0; i < CERT_NKNOTS; i++) KNOTS[i] = (double)H.knot_num[i] / H.Q_K;

    double *gInvD = calloc(N, sizeof(double));
    for (int u = 0; u < N; u++) gInvD[u] = (C.gKind[u] == 1) ? 0.0 : (1.0 / C.gD[u]);

    double *emp_rpar = calloc(N, sizeof(double)), *emp_rch = calloc(N, sizeof(double));
    double *emp_edge = calloc((size_t)N * N, sizeof(double));
    run_empirical(&C, 500000, 35, emp_rpar, emp_rch, emp_edge);

    /* Immediately quantize rpar, rch, xpar, xch, xcc to exact integers / S_r! */
    int64_t *rpar_q = calloc(C.INNER_N, sizeof(int64_t)), *rch_q = calloc(C.INNER_N, sizeof(int64_t));
    int64_t *xpar_q = calloc(C.INNER_N, sizeof(int64_t)), *xch_q = calloc(C.INNER_N, sizeof(int64_t)), *xcc_q = calloc(C.INNER_N, sizeof(int64_t));
    double *rpar = calloc(N, sizeof(double)), *rch = calloc(N, sizeof(double));
    double *xpar = calloc(N, sizeof(double)), *xch = calloc(N, sizeof(double)), *xcc = calloc(N, sizeof(double));
    for (int u = 0; u < C.INNER_N; u++) {
        rpar_q[u] = llround(emp_rpar[u] * H.S_r); rpar[u] = (double)rpar_q[u] / H.S_r;
        rch_q[u]  = llround(emp_rch[u]  * H.S_r); rch[u]  = (double)rch_q[u]  / H.S_r;
        if (cl_hu_idx(&C, C.gL[u], C.gK[u], C.gH[u]) >= 0) {
            xpar_q[u] = rpar_q[u]; xpar[u] = rpar[u];
            xch_q[u]  = rch_q[u];  xch[u]  = rch[u];
        }
        if (cl_cu_idx(&C, C.gL[u], C.gK[u], C.gH[u]) >= 0) {
            xcc_q[u] = rch_q[u]; xcc[u] = rch[u];
        }
    }

    int max_rank = C.D0 + 1;
    int *p_by_r[128], p_rcnt[128] = {0}, *c_by_r[128], c_rcnt[128] = {0};
    for (int r = 0; r <= max_rank; r++) { p_by_r[r] = malloc(N * sizeof(int)); c_by_r[r] = malloc(N * sizeof(int)); }
    for (int ip = 0; ip < C.npars; ip++) { int p = C.pars[ip], r = (C.gKind[p]==1)?(C.D0+1):((C.gKind[p]==2)?C.D0:C.gD[p]); p_by_r[r][p_rcnt[r]++] = p; }
    for (int iw = 0; iw < C.nchs;  iw++) { int w = C.chs[iw],  r = (C.gKind[w]==1)?(C.D0+1):((C.gKind[w]==2)?C.D0:C.gD[w]); c_by_r[r][c_rcnt[r]++] = w; }
    for (int r = 2; r <= C.D0; r++) {
        for (int i=0;i<p_rcnt[r];i++) for (int j=i+1;j<p_rcnt[r];j++) if (C.gK[p_by_r[r][j]]>C.gK[p_by_r[r][i]]) { int t=p_by_r[r][i];p_by_r[r][i]=p_by_r[r][j];p_by_r[r][j]=t; }
        for (int i=0;i<c_rcnt[r];i++) for (int j=i+1;j<c_rcnt[r];j++) if (C.gK[c_by_r[r][j]]>C.gK[c_by_r[r][i]]) { int t=c_by_r[r][i];c_by_r[r][i]=c_by_r[r][j];c_by_r[r][j]=t; }
    }

    size_t hv_len = (size_t)C.npars * C.NHU_ACT * C.nchs;
    size_t hc_len = (size_t)C.NCU_ACT * C.nchs * (C.nchs + 1) / 2;
    double *f = calloc(C.NV, sizeof(double)), *fn = calloc(C.NV, sizeof(double));
    double *g = calloc((size_t)N*N, sizeof(double)), *gn = calloc((size_t)N*N, sizeof(double));
    double *hV = calloc(hv_len, sizeof(double)), *hVn = calloc(hv_len, sizeof(double));
    double *hC = calloc(hc_len, sizeof(double)), *hCn = calloc(hc_len, sizeof(double));
    double *dLf_g = calloc(N, sizeof(double)), *A_in = calloc((size_t)N*N, sizeof(double)), *A_out = calloc((size_t)N*N, sizeof(double));
    double (*k_pc)[CERT_NKNOTS] = calloc(C.INNER_N, sizeof(*k_pc));
    double (*k_cp)[CERT_NKNOTS] = calloc(C.INNER_N, sizeof(*k_cp));
    double (*k_cc)[CERT_NKNOTS] = calloc(C.INNER_N, sizeof(*k_cc));

    double lam_target = (double)P / Q;
    /* Solve f at lam_f slightly stricter than lam_target to guarantee strict vertex slack on EVERY f(u) > 0! */
    double lam_f = (mode == 1) ? (lam_target - 2.5e-5) : (lam_target + 1.8e-5);
    double rew_g = 2e-7, rew_h = 5e-8;
    double sgn = (mode == 1) ? 1.0 : -1.0;

    for (int u = 0; u < C.NV; u++) f[u] = 1.0;
    if (mode == 2) f[C.VSTAR0] = f[C.VSTAR1] = 0.0; else f[C.VSTAR0] = f[C.VSTAR1] = 1.0 / lam_f;

    printf("Building [%s] certificate for %d/%d = %.6f (lam_f=%.6f, S=%lld, S_r=%lld)...\n",
           mode==1?"UPPER":"LOWER", P, Q, lam_target, lam_f, (long long)H.S, (long long)H.S_r);
    fflush(stdout);

    double rho = 1.0;
    for (int it = 0; it < 16; it++) {
        if (mode == 2) { fn[C.VSTAR0] = fn[C.VSTAR1] = 0.0; }
        else {
            double g0_inf = g[C.INF_ID*N + iw_nb];
            double g0_hi  = (1.0 - (double)C.D0/(C.KF+1)) * g0_inf + ((double)C.D0/(C.KF+1)) * g[cl_t1_idx(&C,C.KT+1,0)*N + iw_nb];
            double g1_hi  = (1.0 - (double)C.D0/(C.KF+2)) * g0_inf + ((double)C.D0/(C.KF+2)) * g[cl_t1_idx(&C,C.KT+1,1)*N + iw_nb];
            fn[C.VSTAR0] = fmax(0.0, f[v_nb] + fmax(g0_inf, g0_hi)) / lam_f;
            fn[C.VSTAR1] = fmax(0.0, f[v_nb] + fmax(g0_inf, g1_hi)) / lam_f;
        }
        for (int hb=1; hb>=0; hb--) for (int k=C.KF; k>=0; k--) for (int L=C.LF; L>=0; L--) {
            int u = cl_vidx(&C,L,k,hb), D = L+k+hb; if (D < 2) { fn[u] = 0; continue; }
            double bL = (double)L / D;
            if (L == C.LF) {
                int ub = cl_vidx(&C,L-1,k+1,hb), ubt = cl_vidx(&C,C.LF,k+1,hb);
                double snb = f[v_nb] + eval_g_p(&C, g, L-1, k+1, hb, iw_nb);
                double t_hi = (double)C.D0 / (C.LF + 1 + k + hb);
                double g_t0 = g[C.INF_ID*N + iw_nb];
                double g_t1 = (1.0 - t_hi)*g_t0 + t_hi*g[cl_t1_idx(&C, k+1, hb)*N + iw_nb];
                double snbt = f[v_nb] + ((mode==1) ? fmax(g_t0, g_t1) : fmin(g_t0, g_t1));
                double blo = (double)(C.LF+1) / (C.LF + 1 + k + hb);
                double v1=(bL*(fn[ub]+snb))/(lam_f+bL), v2=(blo*(fn[ubt]+snbt))/(lam_f+blo), v3=(fn[ubt]+snbt)/(lam_f+1.0);
                fn[u] = (mode==1) ? fmax(v1, fmax(v2, v3)) : fmin(v1, fmin(v2, v3));
                continue;
            }
            double f_up = fn[cl_vidx(&C,L+1,k,hb)];
            int ub = (L>=1) ? cl_vidx(&C,L-1,k+1,hb) : u;
            double bv = (L>=1) ? (fn[ub] + f[v_nb] + eval_g_p(&C, g, L-1, k+1, hb, iw_nb)) : 0.0;
            if (D < C.D0) {
                int gu = cl_inner_idx(&C,L,k,hb);
                double up = L + hb*rpar[gu] + k*rch[gu];
                fn[u] = (up*f_up + (L>=1?bL*bv:0.0)) / (lam_f + up + (L>=1?bL:0.0));
            } else {
                double r0=L, r1=L+0.5*(hb+k);
                double v0=(r0*f_up+(L>=1?bL*bv:0))/(lam_f+r0+(L>=1?bL:0));
                double v1=(r1*f_up+(L>=1?bL*bv:0))/(lam_f+r1+(L>=1?bL:0));
                fn[u] = (mode==1) ? fmax(v0, v1) : fmin(v0, v1);
            }
        }
        memset(dLf_g, 0, N*sizeof(double));
        for (int u=0; u<C.INNER_N; u++) dLf_g[u] = fn[cl_vidx(&C,C.gL[u]+1,C.gK[u],C.gH[u])] - fn[cl_vidx(&C,C.gL[u],C.gK[u],C.gH[u])];

        memset(A_in, 0, (size_t)N*N*sizeof(double)); memset(A_out, 0, (size_t)N*N*sizeof(double));
        int do_opt = (it >= 1 && it <= 10);
        #pragma omp parallel for schedule(dynamic) num_threads(3)
        for (int u=0; u<C.INNER_N; u++) {
            int Lu=C.gL[u], ku=C.gK[u], hu=C.gH[u];
            int up = (C.gD[u]+1 < C.D0) ? cl_inner_idx(&C,Lu+1,ku,hu) : cl_t1_idx(&C,ku,hu);
            int in_hu = (cl_hu_idx(&C,Lu,ku,hu)>=0), in_cu = (cl_cu_idx(&C,Lu,ku,hu)>=0);
            double *xs_b=malloc(N*sizeof(double)), *wxs_b=malloc(N*sizeof(double)), *ys_b=malloc(N*sizeof(double)), *wys_b=malloc(N*sizeof(double));
            if (!in_hu && hu==1 && ku>=1) {
                if (do_opt) {
                    int nx=0, ny=0;
                    for (int ip=0; ip<C.npars; ip++) { int p=C.pars[ip]; double wt=emp_edge[p*N+u]+1e-4/(C.gKind[p]==0?(double)C.gD[p]*C.gD[p]:1e4); xs_b[nx]=gInvD[p]; wxs_b[nx]=wt; nx++; }
                    for (int iw=0; iw<C.nchs; iw++) { int w=C.chs[iw]; double wt=emp_edge[u*N+w]+1e-4/(C.gKind[w]==0?(double)C.gD[w]*C.gD[w]:1e4); ys_b[ny]=g[up*N+w]-g[u*N+w]; wys_b[ny]=wt; ny++; }
                    opt_knots(KNOTS, k_pc[u], sgn, xs_b, wxs_b, nx, ys_b, wys_b, ny);
                }
                for (int ip=0; ip<C.npars; ip++) { int p=C.pars[ip]; A_in[p*N+u] += ku * eval_pl(KNOTS, k_pc[u], gInvD[p]); }
                if (do_opt) {
                    int nx=0, ny=0;
                    for (int iw=0; iw<C.nchs; iw++) { int w=C.chs[iw]; double wt=emp_edge[u*N+w]+1e-4/(C.gKind[w]==0?(double)C.gD[w]*C.gD[w]:1e4); xs_b[nx]=gInvD[w]; wxs_b[nx]=wt; nx++; }
                    for (int ip=0; ip<C.npars; ip++) { int p=C.pars[ip]; double wt=emp_edge[p*N+u]+1e-4/(C.gKind[p]==0?(double)C.gD[p]*C.gD[p]:1e4); ys_b[ny]=g[p*N+up]-g[p*N+u]; wys_b[ny]=wt; ny++; }
                    opt_knots(KNOTS, k_cp[u], sgn, xs_b, wxs_b, nx, ys_b, wys_b, ny);
                }
                for (int iw=0; iw<C.nchs; iw++) { int w=C.chs[iw]; A_out[u*N+w] += eval_pl(KNOTS, k_cp[u], gInvD[w]); }
            }
            if (!in_cu && ku>=2) {
                if (do_opt) {
                    int nx=0, ny=0;
                    for (int iw=0; iw<C.nchs; iw++) { int w=C.chs[iw]; double wt=emp_edge[u*N+w]+1e-4/(C.gKind[w]==0?(double)C.gD[w]*C.gD[w]:1e4); xs_b[nx]=gInvD[w]; wxs_b[nx]=wt; nx++; ys_b[ny]=g[up*N+w]-g[u*N+w]; wys_b[ny]=wt; ny++; }
                    opt_knots(KNOTS, k_cc[u], sgn, xs_b, wxs_b, nx, ys_b, wys_b, ny);
                }
                for (int iw=0; iw<C.nchs; iw++) { int w=C.chs[iw]; A_out[u*N+w] += (ku-1) * eval_pl(KNOTS, k_cc[u], gInvD[w]); }
            }
            free(xs_b); free(wxs_b); free(ys_b); free(wys_b);
        }

        solve_g_hv_hc(&C, KNOTS, mode, lam_target, rew_g, rew_h,
                      gInvD, rpar, rch, xpar, xch, xcc, dLf_g, A_in, A_out,
                      k_pc, k_cp, k_cc, p_by_r, p_rcnt, c_by_r, c_rcnt,
                      hV, hC, gn, hVn, hCn);

        rho = fn[v_nb] / f[v_nb];
        double inv_s = 1.0 / fn[v_nb], diff = 0.0;
        for (int i=0; i<C.NV; i++) { double v=fn[i]*inv_s; if (fabs(v-f[i])>diff) diff=fabs(v-f[i]); f[i]=v; }
        for (size_t i=0; i<(size_t)N*N; i++) { double v=gn[i]*inv_s; if (fabs(v-g[i])>diff) diff=fabs(v-g[i]); g[i]=v; }
        for (size_t i=0; i<hv_len; i++) { double v=hVn[i]*inv_s; if (fabs(v-hV[i])>diff) diff=fabs(v-hV[i]); hV[i]=v; }
        for (size_t i=0; i<hc_len; i++) { double v=hCn[i]*inv_s; if (fabs(v-hC[i])>diff) diff=fabs(v-hC[i]); hC[i]=v; }
        if (it % 4 == 0 || it == 15 || diff < 1e-9) {
            printf("  it=%2d rho=%.8f diff=%.2e\n", it, rho, diff); fflush(stdout);
        }
        if (diff < 1e-9 && it >= 12) break;
    }

    /* Step B: Quantize knots to exact int64_t / S and freeze! */
    int64_t (*k_pc_q)[CERT_NKNOTS] = calloc(C.INNER_N, sizeof(*k_pc_q));
    int64_t (*k_cp_q)[CERT_NKNOTS] = calloc(C.INNER_N, sizeof(*k_cp_q));
    int64_t (*k_cc_q)[CERT_NKNOTS] = calloc(C.INNER_N, sizeof(*k_cc_q));
    for (int u = 0; u < C.INNER_N; u++) {
        for (int k = 0; k < CERT_NKNOTS; k++) {
            k_pc_q[u][k] = llround(k_pc[u][k] * H.S); k_pc[u][k] = (double)k_pc_q[u][k] / H.S;
            k_cp_q[u][k] = llround(k_cp[u][k] * H.S); k_cp[u][k] = (double)k_cp_q[u][k] / H.S;
            k_cc_q[u][k] = llround(k_cc[u][k] * H.S); k_cc[u][k] = (double)k_cc_q[u][k] / H.S;
        }
    }
    memset(A_in, 0, (size_t)N*N*sizeof(double)); memset(A_out, 0, (size_t)N*N*sizeof(double));
    for (int u=0; u<C.INNER_N; u++) {
        int Lu=C.gL[u], ku=C.gK[u], hu=C.gH[u];
        int in_hu = (cl_hu_idx(&C,Lu,ku,hu)>=0), in_cu = (cl_cu_idx(&C,Lu,ku,hu)>=0);
        if (!in_hu && hu==1 && ku>=1) {
            for (int ip=0; ip<C.npars; ip++) { int p=C.pars[ip]; A_in[p*N+u] += ku * eval_pl(KNOTS, k_pc[u], gInvD[p]); }
            for (int iw=0; iw<C.nchs; iw++)  { int w=C.chs[iw];  A_out[u*N+w] += eval_pl(KNOTS, k_cp[u], gInvD[w]); }
        }
        if (!in_cu && ku>=2) {
            for (int iw=0; iw<C.nchs; iw++)  { int w=C.chs[iw];  A_out[u*N+w] += (ku-1) * eval_pl(KNOTS, k_cc[u], gInvD[w]); }
        }
    }

    /* Step C: Unnormalized final solve of f (at lam_f), then 2 exact DAG passes of (g, hV, hC) at lam_target! */
    if (mode == 2) { fn[C.VSTAR0] = fn[C.VSTAR1] = 0.0; }
    else {
        double g0_inf = g[C.INF_ID*N + iw_nb];
        double g0_hi  = (1.0 - (double)C.D0/(C.KF+1)) * g0_inf + ((double)C.D0/(C.KF+1)) * g[cl_t1_idx(&C,C.KT+1,0)*N + iw_nb];
        double g1_hi  = (1.0 - (double)C.D0/(C.KF+2)) * g0_inf + ((double)C.D0/(C.KF+2)) * g[cl_t1_idx(&C,C.KT+1,1)*N + iw_nb];
        fn[C.VSTAR0] = fmax(0.0, f[v_nb] + fmax(g0_inf, g0_hi)) / lam_f;
        fn[C.VSTAR1] = fmax(0.0, f[v_nb] + fmax(g0_inf, g1_hi)) / lam_f;
    }
    for (int hb=1; hb>=0; hb--) for (int k=C.KF; k>=0; k--) for (int L=C.LF; L>=0; L--) {
        int u = cl_vidx(&C,L,k,hb), D = L+k+hb; if (D < 2) { fn[u] = 0; continue; }
        double bL = (double)L / D;
        if (L == C.LF) {
            int ub = cl_vidx(&C,L-1,k+1,hb), ubt = cl_vidx(&C,C.LF,k+1,hb);
            double snb = f[v_nb] + eval_g_p(&C, g, L-1, k+1, hb, iw_nb);
            double t_hi = (double)C.D0 / (C.LF + 1 + k + hb);
            double g_t0 = g[C.INF_ID*N + iw_nb];
            double g_t1 = (1.0 - t_hi)*g_t0 + t_hi*g[cl_t1_idx(&C, k+1, hb)*N + iw_nb];
            double snbt = f[v_nb] + ((mode==1) ? fmax(g_t0, g_t1) : fmin(g_t0, g_t1));
            double blo = (double)(C.LF+1) / (C.LF + 1 + k + hb);
            double v1=(bL*(fn[ub]+snb))/(lam_f+bL), v2=(blo*(fn[ubt]+snbt))/(lam_f+blo), v3=(fn[ubt]+snbt)/(lam_f+1.0);
            fn[u] = (mode==1) ? fmax(v1, fmax(v2, v3)) : fmin(v1, fmin(v2, v3));
            continue;
        }
        double f_up = fn[cl_vidx(&C,L+1,k,hb)];
        int ub = (L>=1) ? cl_vidx(&C,L-1,k+1,hb) : u;
        double bv = (L>=1) ? (fn[ub] + f[v_nb] + eval_g_p(&C, g, L-1, k+1, hb, iw_nb)) : 0.0;
        if (D < C.D0) {
            int gu = cl_inner_idx(&C,L,k,hb);
            double up = L + hb*rpar[gu] + k*rch[gu];
            fn[u] = (up*f_up + (L>=1?bL*bv:0.0)) / (lam_f + up + (L>=1?bL:0.0));
        } else {
            double r0=L, r1=L+0.5*(hb+k);
            double v0=(r0*f_up+(L>=1?bL*bv:0))/(lam_f+r0+(L>=1?bL:0));
            double v1=(r1*f_up+(L>=1?bL*bv:0))/(lam_f+r1+(L>=1?bL:0));
            fn[u] = (mode==1) ? fmax(v0, v1) : fmin(v0, v1);
        }
    }
    /* Quantize final f to exact int64_t / S BEFORE solving final (g, hV, hC)! */
    int64_t *f_q = calloc(C.NV, sizeof(int64_t));
    for (int i = 0; i < C.NV; i++) { f_q[i] = llround(fn[i] * H.S); f[i] = (double)f_q[i] / H.S; }
    memset(dLf_g, 0, N*sizeof(double));
    for (int u = 0; u < C.INNER_N; u++) dLf_g[u] = f[cl_vidx(&C,C.gL[u]+1,C.gK[u],C.gH[u])] - f[cl_vidx(&C,C.gL[u],C.gK[u],C.gH[u])];

    /* 2 passes of (g, hV, hC) at EXACT lam_target using EXACT quantized f and EXACT quantized knots! */
    for (int pass = 0; pass < 2; pass++) {
        solve_g_hv_hc(&C, KNOTS, mode, lam_target, rew_g, rew_h,
                      gInvD, rpar, rch, xpar, xch, xcc, dLf_g, A_in, A_out,
                      k_pc, k_cp, k_cc, p_by_r, p_rcnt, c_by_r, c_rcnt,
                      hV, hC, gn, hVn, hCn);
        memcpy(g, gn, (size_t)N*N*sizeof(double));
        memcpy(hV, hVn, hv_len*sizeof(double));
        memcpy(hC, hCn, hc_len*sizeof(double));
    }

    int64_t *g_q  = calloc((size_t)N*N, sizeof(int64_t));
    int64_t *hV_q = calloc(hv_len, sizeof(int64_t));
    int64_t *hC_q = calloc(hc_len, sizeof(int64_t));
    for (size_t i = 0; i < (size_t)N*N; i++) g_q[i] = llround(g[i] * H.S);
    for (size_t i = 0; i < hv_len; i++) hV_q[i] = llround(hV[i] * H.S);
    for (size_t i = 0; i < hc_len; i++) hC_q[i] = llround(hC[i] * H.S);

    FILE *fp = fopen(out_path, "wb");
    if (!fp) { perror("fopen"); exit(1); }
    fwrite(&H, sizeof(H), 1, fp);
    fwrite(rpar_q, sizeof(int64_t), C.INNER_N, fp);
    fwrite(rch_q,  sizeof(int64_t), C.INNER_N, fp);
    fwrite(xpar_q, sizeof(int64_t), C.INNER_N, fp);
    fwrite(xch_q,  sizeof(int64_t), C.INNER_N, fp);
    fwrite(xcc_q,  sizeof(int64_t), C.INNER_N, fp);
    fwrite(k_pc_q, sizeof(int64_t), (size_t)C.INNER_N * CERT_NKNOTS, fp);
    fwrite(k_cp_q, sizeof(int64_t), (size_t)C.INNER_N * CERT_NKNOTS, fp);
    fwrite(k_cc_q, sizeof(int64_t), (size_t)C.INNER_N * CERT_NKNOTS, fp);
    fwrite(f_q,    sizeof(int64_t), C.NV, fp);
    fwrite(g_q,    sizeof(int64_t), (size_t)N * N, fp);
    fwrite(hV_q,   sizeof(int64_t), hv_len, fp);
    fwrite(hC_q,   sizeof(int64_t), hc_len, fp);
    fclose(fp);

    printf("Wrote %s (f=%d, g=%d*%d, hV=%zu, hC=%zu, Phi(T3)=%lld/%lld=%.6f)\n",
           out_path, C.NV, N, N, hv_len, hc_len,
           (long long)f_q[cl_vidx(&C,2,0,0)], (long long)H.S, (double)f_q[cl_vidx(&C,2,0,0)]/H.S);
    fflush(stdout);

    for (int r=0; r<=max_rank; r++) { free(p_by_r[r]); free(c_by_r[r]); }
    free(gInvD); free(emp_rpar); free(emp_rch); free(emp_edge);
    free(rpar_q); free(rch_q); free(xpar_q); free(xch_q); free(xcc_q);
    free(rpar); free(rch); free(xpar); free(xch); free(xcc);
    free(f); free(fn); free(g); free(gn); free(hV); free(hVn); free(hC); free(hCn);
    free(dLf_g); free(A_in); free(A_out); free(k_pc); free(k_cp); free(k_cc);
    free(k_pc_q); free(k_cp_q); free(k_cc_q); free(f_q); free(g_q); free(hV_q); free(hC_q);
    cl_free(&C);
}

int main(int argc, char **argv) {
    if (argc < 5) {
        fprintf(stderr, "Usage: %s <mode:1=UPPER|2=LOWER> <P> <Q> <out.bin>\n", argv[0]);
        return 1;
    }
    build_cert(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), argv[4]);
    return 0;
}
