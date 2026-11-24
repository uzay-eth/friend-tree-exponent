#include "snapshot.hpp"
#include <iostream>
#include <string>
using namespace mc;

struct EnvelopeAudit {
    Snapshot &s;
    bool build;
    uint64_t n1=0,n2=0;
    static int count(const Level &l,bool parent) { return parent?l.npars:l.nchs; }
    static int id(const Level &l,bool parent,int ix) { return parent?l.pars[ix]:l.chs[ix]; }
    static int index(const Level &l,bool parent,int i) { return parent?l.p_to_pi[i]:l.w_to_wi[i]; }

    template<class Difference>
    void one(const Level &fine,const Level &coarse,bool parent,i64 *m,i64 *r,Difference delta) {
        if(build) {
            std::fill(m,m+coarse.NG,0);
            std::fill(r,r+coarse.NG,0);
            for(int j=0;j<count(coarse,parent);j++) {
                int v=id(coarse,parent,j);
                m[v]=delta(representative(fine,coarse,v));
            }
            std::vector<i64> lo(coarse.NG,INT64_MAX),hi(coarse.NG,INT64_MIN);
            for(int j=0;j<count(fine,parent);j++) {
                int v=id(fine,parent,j);
                Mix p=project(fine,coarse,v);
                if(p.first==p.second) continue;
                i64 y=p.den*delta(v)-(p.den-p.n)*m[p.second];
                lo[p.first]=std::min(lo[p.first],y);
                hi[p.first]=std::max(hi[p.first],y);
            }
            for(int v=coarse.TAIL_BASE;v<coarse.NG;v++) if(lo[v]!=INT64_MAX) {
                const int d=coarse.D0;
                m[v]=floordiv(wide(lo[v])+hi[v],2*d);
                r[v]=ceildiv(std::max(mag(hi[v]-d*m[v]),mag(lo[v]-d*m[v])),d);
            }
        }
        for(int j=0;j<coarse.NG;j++) check(r[j]>=0,"negative 1D radius");
        for(int j=0;j<count(fine,parent);j++) {
            int v=id(fine,parent,j);
            Mix p=project(fine,coarse,v);
            i64 err=p.den*delta(v)-p.n*m[p.first]-(p.den-p.n)*m[p.second];
            check(mag(err)<=p.n*r[p.first]+(p.den-p.n)*r[p.second],"1D envelope failure");
            ++n1;
        }
    }

    template<class Difference>
    void two(const Level &fine,const Level &coarse,bool parent,i64 *m,i64 *r,Difference delta) {
        const int n=count(coarse,parent), nb=coarse.NG, nf=count(fine,parent);
        const auto at=[&](int a,int b){ return (size_t)index(coarse,parent,a)*nb+b; };
        std::vector<Mix> A(nf),B(fine.NG);
        for(int i=0;i<nf;i++) A[i]=project(fine,coarse,id(fine,parent,i));
        for(int i=0;i<fine.NG;i++) B[i]=project(fine,coarse,i);
        if(build) {
            std::fill(r,r+(size_t)n*nb,0);
            for(int ia=0;ia<n;ia++) for(int b=0;b<nb;b++) {
                int a=id(coarse,parent,ia);
                m[at(a,b)]=delta(representative(fine,coarse,a),representative(fine,coarse,b));
            }
        }
        // A one-tail radius is established before the residual at a two-tail cell.
        for(int pass=build?1:3;pass<=3;pass++) {
            for(int ia=0;ia<nf;ia++) for(int ib=0;ib<fine.NG;ib++) {
                auto a=A[ia],b=B[ib];
                int tails=(a.first!=a.second)+(b.first!=b.second);
                if(pass<3 && tails!=pass) continue;
                size_t ix[]={at(a.first,b.first),at(a.first,b.second),at(a.second,b.first),at(a.second,b.second)};
                int w[]={a.n*b.n,a.n*(b.den-b.n),(a.den-a.n)*b.n,(a.den-a.n)*(b.den-b.n)};
                wide err=wide(a.den)*b.den*delta(id(fine,parent,ia),ib), rad=0;
                for(int x=0;x<4;x++) { err-=wide(w[x])*m[ix[x]]; rad+=wide(w[x])*r[ix[x]]; }
                wide distance=err<0?-err:err;
                if(pass<3 && distance>rad)
                    r[ix[0]]=narrow(wide(r[ix[0]])+ceildiv(distance-rad,w[0]));
                if(pass==3) { check(distance<=rad,"2D envelope failure"); ++n2; }
            }
        }
        for(size_t j=0;j<(size_t)n*nb;j++) check(r[j]>=0,"negative 2D radius");
    }
    void run() {
        auto &c=s.C;
        auto &l2=c.L2; auto &l3=c.L3; auto &l4=c.L4;
        for(int a=0;a<c.NHU_ACT;a++) {
            int u=c.hu_list[a], L=u%c.H.DHU,k=u/c.H.DHU+1;
            int i=corner(l2,L,k),j=corner(l2,L+1,k);
            size_t off=(size_t)a*l3.NG;
            one(l2,l3,false,s.e23[0]+off,s.e23[1]+off,[&](int b){ return s.G(j,b)-s.G(i,b); });
            one(l2,l3,true,s.e23[2]+off,s.e23[3]+off,[&](int p){ return s.G(p,j)-s.G(p,i); });
        }
        for(bool fork:{false,true}) {
            int total=fork?c.NF4_ACT:c.NP4_ACT;
            auto env=fork?s.ef:s.ep;
            for(int a=0;a<total;a++) {
                int u=(fork?c.f4_u1:c.p4_u1)[a],v=(fork?c.f4_u2:c.p4_u2)[a];
                int L=l4.gL[u],k=l4.gK[u],M=l4.gL[v],j=l4.gK[v];
                int u3=corner(l3,L,k),up=corner(l3,L+1,k),v3=corner(l3,M,j),vp=corner(l3,M+1,j);
                int mid=c4_hu_idx(&c,M,j);
                size_t o=(size_t)a*l4.NG;
                one(l3,l4,false,env[0]+o,env[1]+o,[&](int b){return s.V(up,mid,b)-s.V(u3,mid,b);});
                if(fork) {
                    int center=c4_cu_idx(&c,L,k);
                    one(l3,l4,false,env[2]+o,env[3]+o,[&](int b){return s.Ch(center,b,vp)-s.Ch(center,b,v3);});
                } else {
                    int center=c4_hu_idx(&c,L,k);
                    one(l3,l4,true,env[2]+o,env[3]+o,[&](int p){return s.V(p,center,vp)-s.V(p,center,v3);});
                }
            }
        }
        for(int a=0;a<c.NY4_ACT;a++) {
            int u=c.y4_u[a],L=l4.gL[u],k=l4.gK[u];
            int mid=c4_hu_idx(&c,L,k),up=c4_hu_idx(&c,L+1,k);
            size_t o=(size_t)a*l4.npars*l4.NG;
            two(l3,l4,true,s.eyv[0]+o,s.eyv[1]+o,[&](int p,int w){return s.V(p,up,w)-s.V(p,mid,w);});
            mid=c4_cu_idx(&c,L,k); up=c4_cu_idx(&c,L+1,k);
            o=(size_t)a*l4.NG*l4.NG;
            two(l3,l4,false,s.eyc[0]+o,s.eyc[1]+o,[&](int x,int y){return s.Ch(up,x,y)-s.Ch(mid,x,y);});
        }
        std::cout<<"envelopes 1D="<<n1<<" 2D="<<n2<<" pass\n";
    }
};
int main(int argc,char **argv) try {
    check(argc==3,"usage: envelopes {build|check} FILE");
    bool build=std::string(argv[1])=="build";
    check(build || std::string(argv[1])=="check","unknown command");
    Snapshot s(argv[2],build);
    EnvelopeAudit{s,build}.run();
    if(build) check(msync(s.base,s.bytes,MS_SYNC)==0,"envelope writeback");
    return 0;
} catch(const std::exception &e) { std::cerr<<e.what()<<"\n"; return 2; }
