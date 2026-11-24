#pragma once
#include "../motif/cert4v_format.h"
#include <algorithm>
#include <array>
#include <bit>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

namespace mc {
using i64 = std::int64_t;
using wide = __int128_t;
inline void check(bool ok, const char *msg) {
    if (!ok) throw std::runtime_error(msg);
}
inline i64 narrow(wide x) {
    check(x >= INT64_MIN && x <= INT64_MAX, "signed 64-bit result outside range");
    return (i64)x;
}
inline i64 ceildiv(wide n, i64 d) {
    check(d > 0, "nonpositive divisor");
    return narrow(n/d + (n%d>0));
}
inline i64 floordiv(wide n, i64 d) { return -ceildiv(-n,d); }
inline i64 mag(i64 n) { return n<0?-n:n; }
using Level = Cert4vLevel;
using Layout = Cert4vLayout;

struct Snapshot {
    int fd;
    size_t bytes;
    void *base;
    Cert4vHeader H;
    Layout C;
    Cert4vSizes Z;
    i64 *cent[5], *fy[3], *f, *g, *v, *c, *p4, *f4, *y4, *t4;
    i64 *e23[4], *ep[4], *ef[4], *eyv[2], *eyc[2];

    explicit Snapshot(const char *path,bool writable=false) {
        static_assert(std::endian::native==std::endian::little);
        static_assert(sizeof(Cert4vHeader)%8==0);
        fd=open(path,writable?O_RDWR:O_RDONLY);
        check(fd>=0,"snapshot open failed");
        struct stat st;
        check(fstat(fd,&st)==0 && st.st_size>=(off_t)sizeof(H) && st.st_size<=(off_t(2)<<30),"snapshot length outside supported range");
        bytes=st.st_size;
        check(pread(fd,&H,sizeof(H),0)==(ssize_t)sizeof(H),"short header");
        /* Bounds precede frozen code, whose index/size arithmetic assumes modest geometry. */
        check(H.DF>=4 && H.DF<=2048 && H.D0_2>=3 && H.D0_2<=72 &&
              H.D0_3>=3 && H.D0_3<=28 && H.D0_4>=3 && H.D0_4<=14,"unsupported degree geometry");
        for (int x:{H.KT_2,H.KT_3,H.KT_4,H.DHU,H.KHU,H.DCU,H.KCU,H.DP4_1,H.DP4_2,H.DP4_S,H.KP4_1,H.KP4_2,
                    H.DF4_1,H.DF4_2,H.DF4_S,H.KF4_1,H.KF4_2,H.DY4,H.KY4,H.DT4,H.KT4})
            check(0<=x && x<=72,"geometry field out of range");
        check(H.nknots>=2 && H.nknots<=16 && H.Q_K>0 && H.Q_K<=1000000,"unsupported knot geometry");
        check(H.S==(i64(1)<<40) && H.S_r==(i64(1)<<30) && H.Q>0 && H.Q<=1000000000 && H.P>0 && H.P<=H.Q,"unsupported scales or target");
        for (i64 x:{H.f_p0,H.f_p1,H.f_p2}) check(x>=0 && x<=16*H.S,"analytic coefficient range");
        check(H.alpha_lo_num>=0 && H.alpha_lo_num<=H.S_r,"tree shift range");
        check(!c4_validate_header(&H),"invalid frozen geometry");
        /* Strict nesting ensures the coarse representative at its D0 exists as a fine inner state;
           permits an immediately adjacent cutoff in the fine successor routine. */
        check(H.D0_4<H.D0_3 && H.D0_3<H.D0_2 && H.KT_4+1<H.D0_4 && H.KT_3+1<H.D0_3,"unsupported projection geometry");
        c4_init_layout(&C,&H); c4_compute_sizes(&C,&Z);
        check(Z.total_bytes==bytes,"payload size mismatch");
        base=mmap(nullptr,bytes,PROT_READ|(writable?PROT_WRITE:0),writable?MAP_SHARED:MAP_PRIVATE,fd,0);
        check(base!=MAP_FAILED,"mmap failed");
        auto cur=(i64*)((char*)base+sizeof(H));
        auto section=[&](size_t n){ auto start=cur; cur+=n; return start; };
        for (auto &a:cent) a=section(Z.n_cent/5);
        for (auto &a:fy) a=section(Z.n_fy/3);
        f=section(Z.n_f);g=section(Z.n_g);v=section(Z.n_hV);c=section(Z.n_hC);
        p4=section(Z.n_hP4);f4=section(Z.n_hF4);y4=section(Z.n_hY4);t4=section(Z.n_hT4);
        for(auto &a:e23) a=section(Z.n_env_23/4);
        for(auto &a:ep) a=section(Z.n_env_p4/4);
        for(auto &a:ef) a=section(Z.n_env_f4/4);
        for(auto &a:eyv) a=section(Z.n_env_y4_hv/2);
        for(auto &a:eyc) a=section(Z.n_env_y4_hc/2);
        check((char*)cur==(char*)base+bytes,"reader offsets");
        for (auto a:cent) for (int i=0;i<C.L2.INNER_N;i++) check(a[i]>=0 && a[i]<=H.S_r/2,"centering range");
        for(auto a=fy[0];a<cur;++a) check(*a>=-16*H.S && *a<=16*H.S,"entry outside arithmetic proof range");
        check(G(C.L2.INF_ID,C.L2.INF_ID)==0,"nonzero g(INF,INF)");
        for(int p=0;p<C.L2.NG;p++) if(C.L2.p_to_pi[p]<0)
            for(int w=0;w<C.L2.NG;w++) check(!G(p,w),"nonzero invalid parent row");
    }
    Snapshot(const Snapshot&)=delete;
    ~Snapshot() { munmap(base,bytes); close(fd); c4_free_layout(&C); }
    i64 G(int p,int w) const { return g[(size_t)p*C.L2.NG+w]; }
    i64 V(int a,int u,int b) const { return u<0?0:v[c4_hv_off(&C,a,u,b)]; }
    i64 Ch(int u,int a,int b) const { return u<0?0:c[c4_hc_off(&C,u,a,b)]; }
};

struct Mix {
    int first, second, n, den;  // weights n/den and (den-n)/den
};
inline int corner(const Level &s,int L,int k) {
    return L+k+1<s.D0?c4_lvl_inner(&s,L,k):c4_lvl_tail(&s,k);
}
inline Mix project(const Level &f,const Level &c,int i) {
    if(i==f.INF_ID) return {c.INF_ID,c.INF_ID,1,1};
    int d=f.gD[i],k=f.gK[i];
    if(d<c.D0) { int v=c4_lvl_inner(&c,f.gL[i],k); return {v,v,1,1}; }
    return {c4_lvl_tail(&c,k),c.INF_ID,c.D0,d};
}
inline int representative(const Level &f,const Level &c,int i) {
    if(i==c.INF_ID) return f.INF_ID;
    int d=c.gD[i],k=c.gK[i];
    return c4_lvl_inner(&f,d-1-k,k);
}
} // namespace mc
