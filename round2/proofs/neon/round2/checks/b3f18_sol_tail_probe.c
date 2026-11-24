/* Evaluate the original add_tail_gmp in isolation on exact synthetic inputs. */
#define main original_verifier_main
#include "../cert/cert_b3f_verify_gmp.c"
#undef main
#include <inttypes.h>

int main(void) {
    GmpWork w;
    gw_init(&w);
    for(;;) {
        int d,kt,kind,k,h,mode,n,ext,transfer;
        int64_t q[17];
        int rc=scanf("%d %d %d %d %d %d %d %d %d",&d,&kt,&kind,&k,&h,&mode,&n,&ext,&transfer);
        if(rc==EOF)break;
        if(rc!=9)return 2;
        for(int j=0;j<17;j++)if(scanf("%" SCNd64,q+j)!=1)return 3;
        zf_set_si(&w.acc,0,1);
        add_tail_gmp(&w.acc,&w,mode,d,kt,1000000000,kind,k,h,q[0],q[1],q[2],q[3],ext,transfer,q[4],q[5],
                     q[6],q[7],q[8],q[9],q[10],q[11],q[12],n,q[13],q[14],q[15],q[16]);
        gmp_printf("%Zd/%Zd\n",w.acc.num,w.acc.den);
    }
    gw_clear(&w);
    return 0;
}
