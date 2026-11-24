#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>

enum { N=10000, R=20000, U=10, J=6 };
static uint64_t state = UINT64_C(240320185);
static uint64_t rng(void) {
    uint64_t z = (state += UINT64_C(0x9e3779b97f4a7c15));
    z = (z ^ (z >> 30))*UINT64_C(0xbf58476d1ce4e5b9);
    z = (z ^ (z >> 27))*UINT64_C(0x94d049bb133111eb);
    return z ^ (z >> 31);
}
static unsigned uniform(unsigned b) {
    uint64_t x, m;
    do {
        x = (uint32_t)rng();
        m = x*b;
    } while ((uint32_t)m < (uint32_t)(-b)%b);
    return (unsigned)(m >> 32);
}

int main(void) {
    static int head[N], degree[N], to[2*N], next[2*N];
    static int hist[U][N], hits[U][J];
    const int cuts[J]={1,10,50,100,500,1000};
    double diff[U]={0}, diff2[U]={0};
    for(int r=0;r<R;r++) {
        head[0]=0; head[1]=1; to[0]=1; to[1]=0;
        next[0]=next[1]=-1; degree[0]=degree[1]=1;
        int e=2;
        for(int n=2;n<N;n++) {
            int v=uniform(n), a=head[v];
            for (unsigned j=uniform(degree[v]);j>0;j--) a=next[a];
            assert(a>=0);
            int w=to[a];
            to[e]=n; next[e]=head[w]; head[w]=e++; degree[w]++;
            to[e]=w; next[e]=-1; head[n]=e++; degree[n]=1;
        }
        assert(e==2*(N-1));
        for(int u=0;u<U;u++) {
            int d=degree[u];
            hist[u][d]++;
            for(int j=0;j<J;j++) hits[u][j] += d>cuts[j];
            if(u+1<U) {
                double x=(d-degree[u+1])/(double)N;
                diff[u]+=x; diff2[u]+=x*x;
            }
        }
    }
    printf("R=%d n=%d SplitMix64 seed=240320185; single thread\n",R,N);
    puts("u mean mean_SE P(>0.0001) P(>0.001) P(>0.005) P(>0.01) P(>0.05) P(>0.1)");
    for(int u=0;u<U;u++) {
        double s=0, q=0;
        for(int d=1;d<N;d++) {
            double x=d/(double)N;
            s+=hist[u][d]*x; q+=hist[u][d]*x*x;
        }
        printf("%2d %.8f %.8f",u+1,s/R,sqrt((q-s*s/R)/(R*(R-1.0))));
        for(int j=0;j<J;j++) printf(" %.5f",hits[u][j]/(double)R);
        putchar('\n');
    }
    puts("u paired mean gap (u minus u+1), paired MC SE");
    for(int u=0;u<U-1;u++) {
        printf("%2d %.8f %.8f\n",u+1,diff[u]/R,sqrt((diff2[u]-diff[u]*diff[u]/R)/(R*(R-1.0))));
    }
    puts("u max empirical survival reversal (next minus current), cut D>");
    for(int u=0;u<U-1;u++) {
        int gap=0,max=0,cut=0;
        for(int d=N-1;d>0;d--) {
            gap+=hist[u+1][d]-hist[u][d];
            if(gap>max) {max=gap; cut=d-1;}
        }
        printf("%2d %.5f %d\n",u+1,max/(double)R,cut);
    }
}
