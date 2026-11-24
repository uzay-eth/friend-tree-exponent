/* Independently reduce packed motif rows to charge extrema, then sweep f.
 * Production only: bounds and sizes are intentionally pinned to the audited files.
 * Values returned are integers in units 1/S (upper) and 1/(4S) (lower).
 */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include <limits.h>
#include "../cert/cert_b3f_format.h"

typedef struct { int64_t l,r; } Ext;
enum { N=210, NG=435, P=397, HU=187, CU=352, P4=35, YT=20, BU=60, BW=35, S3=24 };
static void see(Ext *x,int64_t v) { if(v<x->l)x->l=v; if(v>x->r)x->r=v; }
static Ext init(void) { return (Ext){INT64_MAX,INT64_MIN}; }
static void join(Ext *a, Ext b) { see(a,b.l);see(a,b.r); }
static void sum(Ext *a, Ext b,int k) { a->l+=b.l*k; a->r+=b.r*k; }
static int64_t *load(FILE *f,size_t n) {
    int64_t *a=malloc(n*8);
    if(!a || fread(a,8,n,f)!=n)exit(2);
    return a;
}
static void endpoint(Ext *e,const int64_t *v,size_t n0,int scale) {
    for(size_t i=0;i<n0;i++)see(&e[i%N],scale*v[i]);
}
static void apex(Ext *e,const int64_t *v,size_t n0) {
    for(size_t i=0;i<n0;i++)see(e,v[i]);
}
int main(int argc,char **argv) {
    for(int arg=1;arg<argc;arg++) {
        FILE *fp=fopen(argv[arg],"rb");
        CertB3FHeader h; B3FLayout c;
        if(!fp || fread(&h,sizeof h,1,fp)!=1 || b3f_validate_header(&h))return 1;
        b3f_layout_init(&c,&h);
        if(c.D0!=18||c.KT!=17||c.LF!=14000||c.KF!=180||c.NG!=NG||c.npars!=P||c.nchs!=N||c.NHU!=HU||c.NCU!=CU||c.NP4!=P4||c.NYT!=YT||c.NYBU!=BU||c.NYBW!=BW||c.NS3!=S3)return 3;
        if(fseek(fp,(long)(96+8*(2*c.INNER_N+4*c.NT+c.yp_len+c.yw_len)),SEEK_SET))return 4;
        int64_t *fv=load(fp,c.NV), *v;
        Ext g[N],hv[N]={0},path[N]={0},yb[N]={0},hc[CU],yt[YT],s3[S3],bottom[NG]={{0}};
        size_t np=N*(N+1)/2, ns=(size_t)N*(N+1)*(N+2)/6;
        for(int w=0;w<N;w++)g[w]=init();
        v=load(fp,c.g_len);
        for(int p=0;p<P;p++)for(int w=0;w<N;w++)see(g+w,v[c.pars[p]*NG+c.chs[w]]);
        free(v);
        v=load(fp,c.hv_len); endpoint(hv,v,c.hv_len,1); free(v);
        v=load(fp,c.hc_len); for(int i=0;i<CU;i++){hc[i]=init();apex(hc+i,v+i*np,np);} free(v);
        v=load(fp,c.p4_len); endpoint(path,v,c.p4_len,1); free(v);
        v=load(fp,c.yt_len);
        for(int i=0;i<YT;i++) { yt[i]=init(); for(int p=0;p<P;p++)apex(yt+i,v+(p*YT+i)*np,np); }
        free(v);
        v=load(fp,c.yb_len);
        for(int i=0;i<BU;i++)endpoint(yb,v+(size_t)i*BW*N*N,(size_t)BW*N*N,c.ybu_K[i]-1);
        free(v);
        v=load(fp,c.s3_len); for(int i=0;i<S3;i++){s3[i]=init();apex(s3+i,v+i*ns,ns);} free(v);
        if(fgetc(fp)!=EOF || ferror(fp))return 5;
        fclose(fp);
        for(int i=0;i<N;i++) {
            int j=c.chs[i]; bottom[j]=g[i];sum(bottom+j,hv[i],1);sum(bottom+j,path[i],1);sum(bottom+j,yb[i],1);
        }
        /* Mix extrema separately, as in the theorem: extremizers of the four
         * independent charged terms can be different. */
        Ext tail[19],ta[2][18];
        for(int k=0;k<=18;k++) {
            Ext temp[4]={init(),init(),init(),init()};
            int ids[]={c.INF_ID,b3f_t1(&c,k,1),b3f_t2(&c,k,1),b3f_t3(&c,k,1)};
            for(int j=0;j<4;j++){int i=c.w_to_wi[ids[j]];join(temp,g[i]);join(temp+1,hv[i]);join(temp+2,path[i]);join(temp+3,yb[i]);}
            tail[k]=temp[0]; for(int j=1;j<4;j++)sum(tail+k,temp[j],1);
        }
        for(int h0=0;h0<=1;h0++)for(int k=2;k<=17;k++) {
            ta[h0][k]=(Ext){0,0};
            join(&ta[h0][k],hc[c.cu_of_g[b3f_t1(&c,k,h0)]]);
            join(&ta[h0][k],hc[c.cu_of_g[b3f_t2(&c,k,h0)]]);
            join(&ta[h0][k],hc[c.cu_of_g[b3f_t3(&c,k,h0)]]);
        }
        int64_t lo=INT64_MAX,hi=INT64_MIN;
        for(int h0=0;h0<=1;h0++)for(int k=0;k<=181;k++)for(int L=0;L<=(k==181?0:14000);L++) {
            int D=L+k+h0; if(D<2)continue;
            int i=b3f_vidx(&c,L,k,h0),j=D<18?c.lut[h0][k][L]:-1;
            Ext x={fv[i],fv[i]};
            if(h0)sum(&x,j>=0?bottom[j]:tail[k<18?k:18],1);
            if(k>=2&&k<=17)sum(&x,j>=0?hc[c.cu_of_g[j]]:ta[h0][k],k*(k-1)/2);
            if(h0&&k>=2&&k<=6&&D<9)sum(&x,yt[c.yt_lut[k][L]],k*(k-1)/2);
            if(k>=3&&k<=6&&D<8)sum(&x,s3[c.s3_lut[h0][k][L]],k*(k-1)*(k-2)/6);
            int64_t y=4*x.l+(k-h0)*h.S;
            if(y<lo)lo=y;
            if(x.r>hi)hi=x.r;
        }
        printf("%s: 4Sc0=%" PRId64 " SC0=%" PRId64 " SPhi(T3)=%" PRId64 "\n",argv[arg],lo,hi,fv[b3f_vidx(&c,2,0,0)]);
        free(fv); b3f_layout_free(&c);
    }
    return 0;
}
