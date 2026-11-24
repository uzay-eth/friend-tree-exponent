#pragma once
#include "../snapshot.hpp"
#include <gmpxx.h>
#include <algorithm>
#include <array>
#include <cerrno>
#include <cstring>
#include <string>
#include <sys/wait.h>

namespace motif_reserve {
using mc::i64;
using mc::wide;
using mc::check;
using Rational=mpq_class;
struct Core { int D,c; }; // Ordinary degree and TOTAL core degree, including the parent.
enum class Slot { ParentOfParent, OtherChildOfParent, ChildOfChild };
enum class Exterior { Physical, Continuous };
struct Donor {
    int D;
    bool pooled;
    static Donor exact(int D) { return {D,false}; }
    static Donor pool(int minimum_degree) { return {minimum_degree,true}; }
};
struct Jump { i64 n; int d; };
inline Rational fraction(i64 n,i64 d) {
    check(d>0,"invalid reserve denominator");
    Rational q(n,d);q.canonicalize();return q;
}
inline Rational number(Jump j) { return fraction(j.n,j.d); }
inline i64 floor_ratio(wide n,wide d) {
    check(d>0,"nonpositive interval denominator");return mc::narrow(n/d-(n%d<0));
}
inline i64 floor_q(const Rational &q) {
    mpz_class v;mpz_fdiv_q(v.get_mpz_t(),q.get_num_mpz_t(),q.get_den_mpz_t());
    check(v.fits_slong_p(),"reserve floor outside signed 64 bits");return v.get_si();
}
inline void valid(Core x) { check(x.D>=2 && x.c>=1 && x.c<=x.D,"invalid physical D,c"); }

// Hash the opened file (stdin on an inherited descriptor), without quoting a pathname or using a shell.
inline std::string snapshot_hash(const mc::Snapshot &s) {
    int fd[2];check(pipe(fd)==0,"digest pipe");
    check(lseek(s.fd,0,SEEK_SET)==0,"digest rewind");
    pid_t pid=fork();check(pid>=0,"digest fork");
    if(!pid) {
        close(fd[0]);
        if(dup2(s.fd,STDIN_FILENO)<0 || dup2(fd[1],STDOUT_FILENO)<0) _exit(125);
        close(fd[1]);execlp("sha256sum","sha256sum",(char*)nullptr);_exit(126);
    }
    close(fd[1]);char buf[128];std::string out;
    for(;;) {
        ssize_t n=read(fd[0],buf,sizeof(buf));
        if(n<0 && errno==EINTR) continue;
        check(n>=0,"digest read");if(!n)break;out.append(buf,n);
    }
    close(fd[0]);int status=0;pid_t waited;
    do { waited=waitpid(pid,&status,0); } while(waited<0 && errno==EINTR);
    check(waited==pid && WIFEXITED(status) && WEXITSTATUS(status)==0 &&
          out.size()==68 && out.substr(64)=="  -\n","sha256sum failed");
    return out.substr(0,64);
}

class Scores {
    enum class Kind { Missing, Centered, Plain, FY };
    struct Component {
        Kind kind=Kind::Missing;
        int count=0,family=0,state=0,center=0,real=0,continuous=0;
        Jump jump{0,1};
        std::array<uint8_t,CERT4V_MAX_KNOTS> prefix_best{};
    };
    struct Edge {
        int key=-1;
        std::array<Component,3> component{};
        std::array<i64,2> residual_ceiling{};
        bool unused_ready=false;
    };
    struct Singleton { uint64_t key=0; i64 value=0; };
    mc::Snapshot s;
    int sign;
    Rational lambda;
    i64 fcap;
    std::vector<bool> real_knot;
    std::array<int,73> interval_for_degree{};
    std::vector<Edge> edges;
    std::vector<Singleton> singletons;
    static constexpr size_t cache_size=1<<15;
    static size_t bucket(uint64_t k) {return (k*UINT64_C(11400714819323198485))>>(64-15);}
    const mc::Level &l2() const {return s.C.L2;}
    int id(Core x) const {return mc::corner(l2(),x.D-x.c,x.c-1);}
    i64 G(int p,int w) const {return sign*s.G(p,w);}
    i64 knot(int family,int state,int i) const {
        return sign*s.fy[family][(size_t)state*s.H.nknots+i];
    }
    int interval(int D) const {
        return D>l2().D0?0:interval_for_degree[D];
    }
    Rational signed_primal(int family,int state,int D) const {
        int j=interval(D),a=s.H.knot_num[j],b=s.H.knot_num[j+1];
        i64 u=(i64)b*D-s.H.Q_K,v=s.H.Q_K-(i64)a*D;
        return (fraction(u,(i64)D*(b-a))*knot(family,state,j)
               +fraction(v,(i64)D*(b-a))*knot(family,state,j+1));
    }
    Rational w(i64 D,i64 c) const {
        i64 L=D-c,k=c-1;
        if(D<s.H.DF) return sign*s.f[c4_vidx(&s.C,(int)L,(int)k)];
        Rational q=fraction(L,D);
        if(sign==1) return s.H.f_p0+s.H.f_p1*q;
        return -(s.H.f_p0*q-s.H.f_p1*q*q-s.H.f_p2*q/D);
    }
    Rational newborn_edge(i64 D,i64 k) const {
        const auto &l=l2();
        if(D<l.D0) return G(c4_lvl_inner(&l,(int)(D-1-k),(int)k),0);
        int t=c4_lvl_tail(&l,(int)std::min<i64>(k,l.KT+1));
        return G(l.INF_ID,0)+fraction(l.D0,D)*(G(t,0)-G(l.INF_ID,0));
    }
    i64 singleton_cold(Core x) const {
        i64 D=x.D,c=x.c,L=D-c,k=c-1;
        Rational old=w(D,c),delta=w(D+1,c)-old,r;
        if(D<l2().D0) {
            int u=id(x);
            r=L+fraction(s.cent[0][u]+k*s.cent[1][u],s.H.S_r);
        } else r=L+(delta>0?fraction(c,2):Rational(0));
        Rational residual=r*delta-lambda*old;
        if(L) residual+=fraction(L,D)*(w(D,c+1)-old+w(2,1)+newborn_edge(D,k+1));
        i64 score=floor_q(fcap-residual);
        check(score>=0,"singleton exceeded frozen universal cap");return score;
    }
    Jump difference(Core x,int p,int b,bool parent) const {
        int k=x.c-1,L=x.D-x.c;
        auto at=[&](int i){return parent?G(i,b):G(p,i);};
        int from=parent?p:b;
        if(x.D<l2().D0) return {at(mc::corner(l2(),L+1,k))-at(from),1};
        return {at(l2().INF_ID)-at(from),x.D+1}; // D=d2, so pre-gain t=1
    }
    template<class Value> Rational project3(Core x,Value value) const {
        const auto &l=s.C.L3;
        if(x.D<l.D0) return value(mc::corner(l,x.D-x.c,x.c-1));
        i64 v=value(l.INF_ID);
        return v+fraction(l.D0,x.D)*(value(c4_lvl_tail(&l,x.c-1))-v);
    }
    wide knot_score(const Component &z,int j) const {
        return wide(s.H.knot_num[j])*z.jump.n-wide(s.H.Q_K)*z.jump.d*knot(z.family,z.state,j);
    }
    Component classify(Core x,Slot slot,Jump d) const {
        Component z;z.jump=d;z.state=id(x);
        z.count=slot==Slot::ParentOfParent?1:slot==Slot::OtherChildOfParent?x.c-2:x.c-1;
        if(!z.count)return z;
        if(x.D>=l2().D0){z.kind=Kind::Plain;return z;}
        int L=x.D-x.c,k=x.c-1;
        bool sibling=slot==Slot::OtherChildOfParent;
        z.center=slot==Slot::ParentOfParent?2:sibling?4:3;
        bool covered=(sibling?c4_cu_idx(&s.C,L,k):c4_hu_idx(&s.C,L,k))>=0;
        if(covered){z.kind=Kind::Centered;return z;}
        z.kind=Kind::FY;z.family=slot==Slot::ParentOfParent?0:sibling?2:1;
        wide continuous_best=knot_score(z,0),physical_best=continuous_best;
        for(int j=1;j<s.H.nknots;j++) {
            wide q=knot_score(z,j);
            if(q>continuous_best) {continuous_best=q;z.continuous=j;}
            if(real_knot[j] && q>physical_best) {physical_best=q;z.real=j;}
            z.prefix_best[j]=z.continuous;
        }
        return z;
    }
    Rational maximum(const Component &z,Exterior mode) const {
        if(z.kind==Kind::Missing)return 0;
        if(z.kind==Kind::Centered)return number(z.jump)*fraction(s.cent[z.center][z.state],s.H.S_r);
        if(z.kind==Kind::Plain)return z.jump.n>0?number(z.jump)/2:Rational(0);
        int j=mode==Exterior::Physical?z.real:z.continuous;
        return fraction(s.H.knot_num[j],s.H.Q_K)*number(z.jump)-knot(z.family,z.state,j);
    }
    Rational pair_base(Core p,Core b,Jump dp,Jump db) const {
        const int a=id(p),z=id(b),kp=p.c-1,kb=b.c-1,L=p.D-p.c,M=b.D-b.c;
        i64 old=G(a,z);
        Rational r=-lambda*old+(L+fraction(1,b.D))*number(dp)+(M+fraction(1,p.D))*number(db);
        if(p.D<l2().D0) {
            r+=(fraction(1,b.D)-fraction(s.cent[1][a],s.H.S_r))*(w(p.D+1,p.c)-w(p.D,p.c));
            if(c4_hu_idx(&s.C,L,kp)<0) r+=signed_primal(1,a,b.D);
            if(kp>=2 && c4_cu_idx(&s.C,L,kp)<0) r+=(kp-1)*signed_primal(2,a,b.D);
        }
        if(b.D<l2().D0) {
            r+=(fraction(1,p.D)-fraction(s.cent[0][z],s.H.S_r))*(w(b.D+1,b.c)-w(b.D,b.c));
            if(kb && c4_hu_idx(&s.C,M,kb)<0)r+=kb*signed_primal(0,z,p.D);
        }
        if(L) {
            int ci=c4_cu_idx(&s.C,L-1,kp+1);
            Rational born=ci<0?Rational(0):project3(b,[&](int v){return sign*s.Ch(ci,v,0);});
            r+=fraction(L,p.D)*(G(id({p.D,p.c+1}),z)-old+born);
        }
        if(M) {
            int ui=c4_hu_idx(&s.C,M-1,kb+1);
            Rational born=ui<0?Rational(0):project3(p,[&](int v){return sign*s.V(v,ui,0);});
            r+=fraction(M,b.D)*(G(a,id({b.D,b.c+1}))-old+born);
        }
        return r;
    }
    Edge &edge(Core p,Core b) {
        valid(p);valid(b);check(p.c>=2,"an edge parent has a core child");
        check(p.D<=l2().D0 && b.D<=l2().D0,"outside finite edge domain");
        auto ix=[](Core x){return x.D*(x.D-1)/2-1+x.c-1;};
        int key=ix(p)*2627+ix(b);
        Edge &e=edges[bucket(key)];
        if(e.key==key)return e;
        Jump dp=difference(p,id(p),id(b),true),db=difference(b,id(p),id(b),false);
        e.component={classify(p,Slot::ParentOfParent,dp),classify(p,Slot::OtherChildOfParent,dp),classify(b,Slot::ChildOfChild,db)};
        e.unused_ready=false;
        e.key=key;return e;
    }
    void compute_unused(Edge &e,Core p,Core b) const {
        Rational base=pair_base(p,b,e.component[0].jump,e.component[2].jump);
        for(auto mode:{Exterior::Physical,Exterior::Continuous}) {
            Rational r=base;
            for(const auto &z:e.component)if(z.count)r+=z.count*maximum(z,mode);
            e.residual_ceiling[(size_t)mode]=mc::narrow(-wide(floor_q(-r)));
        }
        e.unused_ready=true;
    }
    i64 known_fy_gap(const Component &z,int best,int D,int j) const {
        int a=s.H.knot_num[j],b=s.H.knot_num[j+1];
        i64 u=(i64)b*D-s.H.Q_K,v=s.H.Q_K-(i64)a*D;
        wide primal=wide(u)*knot(z.family,z.state,j)+wide(v)*knot(z.family,z.state,j+1);
        wide at=wide(z.jump.n)*(b-a)-wide(z.jump.d)*primal;
        wide n=knot_score(z,best)*D*(b-a)-wide(s.H.Q_K)*at;
        return floor_ratio(n,wide(s.H.Q_K)*z.jump.d*D*(b-a));
    }
public:
    explicit Scores(const char *file):s(file),sign(s.H.mode==1?1:-1),
        lambda(sign==1?fraction(56581165,100000000):fraction(5656453,10000000)),
        fcap(sign==1?154457:95343),edges(cache_size),singletons(cache_size) {
        const char *expected=sign==1?"cae73eb2d4241ace1e566e257c052e4b908a31b9346b624cf63bb38cec4f8f91":
                                      "5030698d2ef127d1f57684b8f570021fa726e19e91199ecee0319b3b3d994c04";
        check(snapshot_hash(s)==expected,"reserve snapshot differs from frozen certificate");
        for(int j=0;j<s.H.nknots;j++) {
            int a=s.H.knot_num[j];real_knot.push_back(a==0 || s.H.Q_K%a==0);
        }
        check((i64)l2().D0*s.H.knot_num[1]==s.H.Q_K,"unexpected initial reciprocal segment");
        // Between adjacent knots the conjugate is affine. Check that every non-knot
        // realizable reciprocal lies between two realizable knots. Below 1/d2 this holds too.
        for(int D=2;D<=l2().D0;D++) {
            int j=0;
            while(j+2<s.H.nknots && D*s.H.knot_num[j+1]<s.H.Q_K)++j;
            interval_for_degree[D]=j;
            if(D*s.H.knot_num[j]!=s.H.Q_K && D*s.H.knot_num[j+1]!=s.H.Q_K)
                check(real_knot[j] && real_knot[j+1],"discrete FY grid unsupported");
        }
    }
    i64 singleton_floor(Core x) {
        valid(x);uint64_t key=(uint64_t(x.D)<<32)|(unsigned)x.c;
        auto &record=singletons[bucket(key)];
        if(record.key!=key){record.value=singleton_cold(x);record.key=key;}
        return record.value;
    }
    i64 edge_residual_ceiling(Core p,Core b,Exterior mode=Exterior::Physical) {
        check((unsigned)mode<=1,"invalid exterior domain");
        Edge &e=edge(p,b);
        if(!e.unused_ready)compute_unused(e,p,b);
        return e.residual_ceiling[(size_t)mode];
    }
    i64 edge_unused_floor(Core p,Core b,Exterior mode=Exterior::Physical) {
        i64 residual=edge_residual_ceiling(p,b,mode);
        return residual>=0?0:mc::narrow(-wide(residual));
    }
    i64 edge_donor_floor(Core p,Core b,Slot slot,Donor donor,Exterior mode=Exterior::Physical) {
        check(donor.D>=2 && (unsigned)slot<3 && (unsigned)mode<=1,"invalid external donor");
        const auto &z=edge(p,b).component[(size_t)slot];
        if(z.kind==Kind::Missing || z.kind==Kind::Centered)return 0;
        if(z.kind==Kind::Plain) {
            wide n=z.jump.n;
            if(n>=0)return floor_ratio(n*(donor.D-2),wide(2)*z.jump.d*donor.D);
            return donor.pooled?0:floor_ratio(-n,wide(z.jump.d)*donor.D);
        }
        int best=mode==Exterior::Physical?z.real:z.continuous;
        int j=interval(donor.D);
        i64 gap=known_fy_gap(z,best,donor.D,j);
        if(donor.pooled) {
            wide prefix_gap=knot_score(z,best)-knot_score(z,z.prefix_best[j]);
            gap=std::min(gap,floor_ratio(prefix_gap,wide(s.H.Q_K)*z.jump.d));
            // The interval can contain an unrealizable knot. Zero is still a valid physical credit.
            return std::max<i64>(0,gap);
        }
        check(gap>=0,"negative exact FY donor credit");return gap;
    }
    int direction() const {return sign;}
    i64 singleton_cap() const {return fcap;}
};
} // namespace motif_reserve
