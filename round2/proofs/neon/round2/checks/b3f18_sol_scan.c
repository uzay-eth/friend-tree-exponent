/* Standalone sequential scan of all production table bytes, no GMP/verifier. */
#include <stdio.h>
#include <stdint.h>
#include <inttypes.h>
#include "../cert/cert_b3f_format.h"

int main(int argc, char **argv) {
    static int64_t buf[1<<16];
    for (int i=1;i<argc;i++) {
        FILE *fp=fopen(argv[i],"rb");
        CertB3FHeader h;
        B3FLayout c;
        if (!fp || fread(&h,sizeof h,1,fp)!=1 || b3f_validate_header(&h)) return 1;
        b3f_layout_init(&c,&h);
        size_t counts[] = {c.INNER_N,c.INNER_N,c.NT,c.NT,c.NT,c.NT,c.yp_len,c.yw_len,c.NV,c.g_len,c.hv_len,c.hc_len,c.p4_len,c.yt_len,c.yb_len,c.s3_len};
        const char *names[] = {"rpar","rch","trpar","trch","c2","c3","Yp","Yw","f","g","hV","hC","P4","Yt","Yb","S3"};
        printf("%s: inner=%d NG=%d pars=%d chs=%d HU=%d CU=%d P4=%d Yt=%d Ybu=%d Ybw=%d S3=%d\n",argv[i],c.INNER_N,c.NG,c.npars,c.nchs,c.NHU,c.NCU,c.NP4,c.NYT,c.NYBU,c.NYBW,c.NS3);
        for (unsigned t=0;t<sizeof counts/sizeof *counts;t++) {
            int64_t lo=INT64_MAX,hi=INT64_MIN;
            for (size_t left=counts[t];left;) {
                size_t n=left < (1<<16) ? left : (1<<16);
                if (fread(buf,sizeof *buf,n,fp)!=n) return 2;
                for (size_t j=0;j<n;j++) {
                    if (buf[j]<lo) lo=buf[j];
                    if (buf[j]>hi) hi=buf[j];
                }
                left-=n;
            }
            printf("%s %zu [%" PRId64 ", %" PRId64 "]\n",names[t],counts[t],lo,hi);
        }
        if (fgetc(fp)!=EOF || ferror(fp)) return 3;
        fclose(fp);
        b3f_layout_free(&c);
    }
    return 0;
}
