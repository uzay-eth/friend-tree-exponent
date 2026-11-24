/*
 * Measure exact 1-body, 2-body, and 3-body (2-hop) statistics & correlations
 * in the TRUE undirected Random Friend Tree (isotropic redirection r=1),
 * plus test the EXACT 1-body master equation identity:
 *   mu * q(L,k) = (L-1)q(L-1,k) - (L + L/D)q(L,k) + ((L+1)/D)q(L+1,k-1) + mu*delta_{(1,0)} + R(L-1,k) - R(L,k)
 */
#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <math.h>

static uint64_t xs[4];
static inline uint64_t rotl(uint64_t x, int k) { return (x << k) | (x >> (64 - k)); }
static inline uint64_t next64(void) {
    uint64_t r = rotl(xs[1] * 5, 7) * 9, t = xs[1] << 17;
    xs[2] ^= xs[0]; xs[3] ^= xs[1]; xs[1] ^= xs[2]; xs[0] ^= xs[3];
    xs[2] ^= t; xs[3] = rotl(xs[3], 45); return r;
}
static uint64_t sm(uint64_t *x) {
    uint64_t z = (*x += 0x9e3779b97f4a7c15ULL);
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    return z ^ (z >> 31);
}
static inline uint32_t unif(uint32_t m) {
    uint64_t x = next64() >> 32, prod = x * (uint64_t)m;
    uint32_t low = (uint32_t)prod;
    if (low < m) {
        uint32_t thr = (uint32_t)(-m) % m;
        while (low < thr) { x = next64() >> 32; prod = x * (uint64_t)m; low = (uint32_t)prod; }
    }
    return (uint32_t)(prod >> 32);
}

int main(void) {
    const int NMAX = 2000000;
    const int RUNS = 40;
    uint64_t seed = 2026112401ULL;
    for (int i = 0; i < 4; i++) xs[i] = sm(&seed);

    uint32_t *par = malloc(NMAX * sizeof(uint32_t));
    uint32_t *deg = malloc(NMAX * sizeof(uint32_t));
    uint32_t *kcnt = calloc(NMAX, sizeof(uint32_t));
    uint32_t *kcap = calloc(NMAX, sizeof(uint32_t));
    uint32_t **kids = calloc(NMAX, sizeof(uint32_t *));
    int *Larr = malloc(NMAX * sizeof(int));
    int *karr = malloc(NMAX * sizeof(int));

    /* Track 2-hop correlations conditioned on middle vertex state (L, k) for D <= 6 */
    double tot_C = 0, tot_S = 0;
    double sum_invDp = 0, sum_invDw = 0, sum_invDp_invDw = 0, cnt_edge = 0;
    /* Conditioned on (L,k) = (0,1), (1,1), (0,2) for vertical path p -> u -> w */
    const int NS = 4;
    const int sL[4] = {0, 1, 0, 2};
    const int sK[4] = {1, 1, 2, 1};
    double vp_cnt[4]={0}, vp_ip[4]={0}, vp_iw[4]={0}, vp_ipiw[4]={0}, vp_ip2[4]={0}, vp_iw2[4]={0};
    /* Sibling cherry u -> {w1, w2} conditioned on u in (0,2), (1,2), (2,2) */
    const int cL[3] = {0, 1, 2};
    double ch_cnt[3]={0}, ch_i1[3]={0}, ch_i1i2[3]={0}, ch_i12[3]={0};
    /* 1-body stats for small (L,k) */
    double q_cnt[8][8] = {{0}}, rpar_sum[8][8] = {{0}}, rch_sum[8][8] = {{0}};

    for (int run = 0; run < RUNS; run++) {
        for (int i = 0; i < NMAX; i++) {
            kcnt[i] = 0;
            if (kcap[i] > 64) { free(kids[i]); kids[i] = NULL; kcap[i] = 0; }
        }
        par[0] = UINT32_MAX; deg[0] = 2; par[1] = 0; deg[1] = 1; par[2] = 0; deg[2] = 1;
        if (kcap[0] < 2) { kids[0] = realloc(kids[0], 2 * sizeof(uint32_t)); kcap[0] = 2; }
        kids[0][0] = 1; kids[0][1] = 2; kcnt[0] = 2;

        for (uint32_t n = 3; n < (uint32_t)NMAX; n++) {
            uint32_t V = unif(n), d = deg[V], r = unif(d);
            uint32_t W = (par[V] != UINT32_MAX) ? ((r == 0) ? par[V] : kids[V][r - 1]) : kids[V][r];
            par[n] = W; deg[n] = 1; deg[W]++;
            if (kcnt[W] == kcap[W]) {
                uint32_t nc = kcap[W] ? 2 * kcap[W] : 2;
                kids[W] = realloc(kids[W], nc * sizeof(uint32_t));
                kcap[W] = nc;
            }
            kids[W][kcnt[W]++] = n;
        }

        for (int u = 0; u < NMAX; u++) {
            if (deg[u] < 2) { Larr[u] = -1; karr[u] = -1; continue; }
            int L = 0;
            for (uint32_t i = 0; i < kcnt[u]; i++) if (deg[kids[u][i]] == 1) L++;
            int h = (par[u] != UINT32_MAX) ? 1 : 0;
            int k = (int)deg[u] - L - h;
            Larr[u] = L; karr[u] = k;
            tot_C += 1.0;
            tot_S += (double)L / deg[u];
            if (h == 1 && L < 8 && k < 8) {
                q_cnt[L][k] += 1.0;
                rpar_sum[L][k] += 1.0 / deg[par[u]];
            }
        }
        for (int w = 1; w < NMAX; w++) {
            if (deg[w] < 2) continue;
            int p = (int)par[w];
            double ip = 1.0 / deg[p], iw = 1.0 / deg[w];
            sum_invDp += ip; sum_invDw += iw; sum_invDp_invDw += ip * iw; cnt_edge += 1.0;
            if (par[p] != UINT32_MAX && Larr[p] < 8 && karr[p] < 8) {
                rch_sum[Larr[p]][karr[p]] += iw;
            }
            /* Vertical path pp -> p -> w */
            if (par[p] != UINT32_MAX) {
                int pp = (int)par[p];
                double ipp = 1.0 / deg[pp];
                for (int si = 0; si < NS; si++) {
                    if (Larr[p] == sL[si] && karr[p] == sK[si]) {
                        vp_cnt[si] += 1.0;
                        vp_ip[si] += ipp; vp_iw[si] += iw;
                        vp_ipiw[si] += ipp * iw;
                        vp_ip2[si] += ipp * ipp; vp_iw2[si] += iw * iw;
                    }
                }
            }
        }
        /* Cherries at u with k=2 */
        for (int u = 1; u < NMAX; u++) {
            if (deg[u] < 2 || karr[u] != 2 || Larr[u] > 2) continue;
            int ci = Larr[u];
            double ws[2]; int wn = 0;
            for (uint32_t i = 0; i < kcnt[u]; i++) {
                uint32_t c = kids[u][i];
                if (deg[c] >= 2) ws[wn++] = 1.0 / deg[c];
            }
            if (wn == 2) {
                ch_cnt[ci] += 1.0;
                ch_i1[ci] += 0.5 * (ws[0] + ws[1]);
                ch_i1i2[ci] += ws[0] * ws[1];
                ch_i12[ci] += 0.5 * (ws[0] * ws[0] + ws[1] * ws[1]);
            }
        }
    }

    printf("N=%d, RUNS=%d, total_C=%.0f, empirical S/C = %.7f\n", NMAX, RUNS, tot_C, tot_S / tot_C);
    double mp = sum_invDp / cnt_edge, mw = sum_invDw / cnt_edge;
    printf("Global edge: E[1/Dp]=%.6f, E[1/Dw]=%.6f, Cov(1/Dp,1/Dw)=%.6f\n",
           mp, mw, sum_invDp_invDw / cnt_edge - mp * mw);

    printf("\nSmall (L, k) empirical q, E[1/Dp | L,k], E[1/Dw | L,k]:\n");
    for (int D = 2; D <= 4; D++) {
        for (int k = 0; k <= D - 1; k++) {
            int L = D - 1 - k;
            double q = q_cnt[L][k] / tot_C;
            double rp = rpar_sum[L][k] / q_cnt[L][k];
            double rc = (k > 0) ? (rch_sum[L][k] / (k * q_cnt[L][k])) : 0.0;
            printf("  (L=%d, k=%d, D=%d): q=%.6f  rpar=%.6f  rch=%.6f\n", L, k, D, q, rp, rc);
        }
    }

    printf("\n2-hop VERTICAL PATH correlation Corr(1/D_pp, 1/D_w | s_u=(L,k)):\n");
    for (int si = 0; si < NS; si++) {
        double n = vp_cnt[si];
        double ep = vp_ip[si]/n, ew = vp_iw[si]/n;
        double cov = vp_ipiw[si]/n - ep * ew;
        double vp = vp_ip2[si]/n - ep*ep, vw = vp_iw2[si]/n - ew*ew;
        double corr = cov / sqrt(vp * vw);
        printf("  u=(L=%d,k=%d): n=%.0f  E[1/Dpp]=%.5f  E[1/Dw]=%.5f  Cov=%.6f  Corr=%.5f\n",
               sL[si], sK[si], n, ep, ew, cov, corr);
    }

    printf("\n2-hop SIBLING CHERRY correlation Corr(1/D_w1, 1/D_w2 | s_u=(L,k=2)):\n");
    for (int ci = 0; ci < 3; ci++) {
        double n = ch_cnt[ci];
        double e1 = ch_i1[ci]/n;
        double cov = ch_i1i2[ci]/n - e1*e1;
        double v1 = ch_i12[ci]/n - e1*e1;
        double corr = cov / v1;
        printf("  u=(L=%d,k=2): n=%.0f  E[1/Dw]=%.5f  Cov=%.6f  Corr=%.5f\n",
               cL[ci], n, e1, cov, corr);
    }

    for (int i = 0; i < NMAX; i++) free(kids[i]);
    free(par); free(deg); free(kcnt); free(kcap); free(kids); free(Larr); free(karr);
    return 0;
}
