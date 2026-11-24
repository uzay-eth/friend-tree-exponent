#include "snapshot.hpp"
#include <fstream>
#include <iostream>
#include <string>
using namespace mc;

struct Bounds {
    i64 lo=0,hi=0;
    void include(i64 v) { lo=std::min(lo,v);hi=std::max(hi,v); }
    void add(Bounds b) { lo+=b.lo;hi+=b.hi; }
    Bounds times(int n) const { return {n*lo,n*hi}; }
};
static Bounds range(const i64 *p,size_t n) {
    check(n>0,"empty active comparison block");
    auto mm=std::minmax_element(p,p+n);return {*mm.first,*mm.second};
}
static Bounds mixture(const std::vector<Bounds> &v,const Level &l,int D,int k) {
    if(D<l.D0) return v[c4_lvl_inner(&l,D-1-k,k)];
    Bounds a=v[c4_lvl_tail(&l,k)], b=v[l.INF_ID];
    return {floordiv(wide(l.D0)*a.lo+wide(D-l.D0)*b.lo,D),
            ceildiv(wide(l.D0)*a.hi+wide(D-l.D0)*b.hi,D)};
}
struct Compare {
    Snapshot &s;
    Layout &c;
    std::vector<Bounds> edge,vert,chain,fork,cherry,y,t;
    i64 lo=INT64_MAX,hi=INT64_MIN;
    std::array<int,2> low_at{},high_at{};
    Compare(Snapshot &snap):s(snap),c(s.C),edge(c.L2.NG,{INT64_MAX,INT64_MIN}),
      vert(c.L3.NG),chain(c.L4.NG),fork(c.L4.NG),cherry(c.NCU_ACT),y(c.NY4_ACT),t(c.NT4_ACT) {
        auto &g=c.L2,&v=c.L3,&q=c.L4;
        for(int i=0;i<g.npars;i++) for(int b=0;b<g.NG;b++) edge[b].include(s.G(g.pars[i],b));
        for(size_t i=0;i<s.Z.n_hV;i++) vert[i%v.NG].include(s.v[i]);
        for(size_t i=0;i<s.Z.n_hP4;i++) chain[i%q.NG].include(s.p4[i]);
        size_t nf=(size_t)q.NG*q.NG;
        for(int a=0;a<c.NF4_ACT;a++) {
            int copies=q.gK[c.f4_u1[a]]-1;
            for(size_t i=0;i<nf;i++) fork[i%q.NG].include(copies*s.f4[(size_t)a*nf+i]);
        }
        for(int a=0;a<c.NCU_ACT;a++) cherry[a]=range(s.c+(size_t)a*v.pair_n,v.pair_n);
        size_t ny=(size_t)q.npars*q.pair_n;
        for(int a=0;a<c.NY4_ACT;a++) y[a]=range(s.y4+(size_t)a*ny,ny);
        for(int a=0;a<c.NT4_ACT;a++) t[a]=range(s.t4+(size_t)a*q.tri_n,q.tri_n);
    }
    Bounds bottom(int D,int k) const {
        Bounds b=mixture(edge,c.L2,D,k);
        b.add(mixture(vert,c.L3,D,k));
        b.add(mixture(chain,c.L4,D,k));
        b.add(mixture(fork,c.L4,D,k));return b;
    }
    void account(int D,int k,Bounds w) {
        // Sum(k-1)=-1. The lower shifted local weight bounds Phi from below.
        i64 shifted=w.lo+(s.H.S/s.H.S_r)*s.H.alpha_lo_num*(k-1);
        if(shifted<lo) {lo=shifted;low_at={D,k};}
        if(w.hi>hi) {hi=w.hi;high_at={D,k};}
    }
    void run(const char *profile) {
        std::ofstream rows,tail;
        if(profile) {
            rows.open(std::string(profile)+".i64",std::ios::binary);
            tail.open(std::string(profile)+"-tail.tsv");
            check(rows.good() && tail.good(),"profile create");
        }
        for(int D=2;D<s.H.DF;D++) for(int k=0;k<D;k++) {
            int L=D-1-k;
            Bounds w=bottom(D,k);i64 f=s.f[c4_vidx(&c,L,k)];
            w.add({f,f});
            int a=c4_cu_idx(&c,L,k), ya=c4_y4_act(&c,L,k), ta=c4_t4_act(&c,L,k);
            if(a>=0) w.add(cherry[c.cu_to_act[a]].times(k*(k-1)/2));
            if(ya>=0) w.add(y[ya].times(k*(k-1)/2));
            if(ta>=0) w.add(t[ta].times(k*(k-1)*(k-2)/6));
            account(D,k,w);
            if(profile) rows.write((const char*)&w.hi,sizeof(w.hi));
        }
        // For D>=DF every bottom interval is affine in 1/D. Bound it by its two endpoints.
        // K covers every child's tail class; the final iteration includes all larger k.
        const int K=c.L2.D0;
        Bounds inf=edge[c.L2.INF_ID];inf.add(vert[c.L3.INF_ID]);
        inf.add(chain[c.L4.INF_ID]);inf.add(fork[c.L4.INF_ID]);
        if(profile) {
            tail<<"clamped_k\tintercept\tslope\n";
            for(int k=0;k<=c.L2.KT+1;k++) {
                i64 slope=0;
                auto take=[&](const std::vector<Bounds> &a,const Level &l) {
                    slope+=l.D0*(a[c4_lvl_tail(&l,k)].hi-a[l.INF_ID].hi);
                };
                take(edge,c.L2);take(vert,c.L3);take(chain,c.L4);take(fork,c.L4);
                tail<<k<<'\t'<<inf.hi<<'\t'<<slope<<'\n';
            }
            rows.close();tail.close();
            check(rows.good() && tail.good(),"writing comparison profile");
        }
        for(int k=0;k<=K+1;k++) {
            Bounds w=bottom(s.H.DF,k);
            w.lo=std::min(w.lo,inf.lo);w.hi=std::max(w.hi,inf.hi);
            if(s.H.mode==1) w.add({s.H.f_p0,s.H.f_p0+s.H.f_p1});
            else {
                // q*(B1 - B2*q - C/D) >= q*(B1-B2-C/DF).
                i64 coeff=s.H.f_p0-s.H.f_p1-ceildiv(s.H.f_p2,s.H.DF);
                check(coeff>=0,"lower tail positivity comparison");
                i64 qmin=(k==K+1)?0:std::max(0,s.H.DF-k-1);
                w.add({floordiv(wide(qmin)*coeff,s.H.DF),s.H.f_p0});
            }
            account(-1,k,w);
        }
        std::cout<<"{\"comparison\":true,\"lower\":"<<lo<<",\"upper\":"<<hi<<",\"unit\":"<<s.H.S
                 <<",\"lo_at\":["<<low_at[0]<<","<<low_at[1]<<"],\"hi_at\":["<<high_at[0]<<","<<high_at[1]
                 <<"],\"finite_rows\":"<<s.Z.n_f<<"}\n";
        check(lo>0 && hi>=lo,"positive comparison not proved");
    }
};
int main(int argc,char **argv) try {
    check(argc==2 || (argc==4 && std::string(argv[2])=="--profile"),"usage: comparison FILE [--profile PREFIX]");
    Snapshot s(argv[1]);Compare(s).run(argc==4?argv[3]:nullptr);return 0;
} catch(const std::exception &e){std::cerr<<e.what()<<"\n";return 2;}
