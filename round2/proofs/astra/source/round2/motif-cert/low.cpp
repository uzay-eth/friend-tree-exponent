#include "snapshot.hpp"
#include <charconv>
#include <iostream>
#include <string>
using namespace mc;

struct Maximum {
    i64 cap=INT64_MIN;
    uint64_t rows=0;
    int i=-1,j=-1;
    void add(i64 x,int a,int b) {
        ++rows;
        if(x>cap) { cap=x;i=a;j=b; }
    }
    void print(const char *name,int P,int Q) const {
        std::cout<<"{\"family\":\""<<name<<"\",\"rows\":"<<rows<<",\"cap\":"<<cap
                 <<",\"witness\":["<<i<<","<<j<<"],\"P\":"<<P<<",\"Q\":"<<Q<<"}\n";
    }
};
struct Frac { wide n; i64 d; };
struct EdgeState {
    int id,k,L,D,kind,up,bud,invden; // invden=0 at INF
    int hu,cu,hu_bud,cu_bud;
    i64 df;
    Mix proj;
};
class Low {
    Snapshot &s;
    const Layout &c;
    const Level &l;
    const int sign,P,Q,nb2,nb3;
    std::vector<EdgeState> st;
    // For any corner 1/d or zero, the two barycentric weights of its FY knot interval.
    struct Knot { int j,n,d; };
    std::vector<Knot> interp;
    i64 signed_f(int L,int k) const { return sign*s.f[c4_vidx(&c,L,k)]; }
    Frac f(int L,int k) const {
        int D=L+k+1;
        if(D<s.H.DF) return {signed_f(L,k),1};
        if(sign==1) return {wide(s.H.f_p0)*D+wide(s.H.f_p1)*L,D};
        return {-(wide(s.H.f_p0)*L*D-wide(s.H.f_p1)*L*L-wide(s.H.f_p2)*L),(i64)D*D};
    }
    i64 g(int p,int w) const { return sign*s.G(p,w); }
    i64 term(i64 x,i64 n,i64 d) const { return ceildiv(wide(x)*n,d); }
    i64 inverse(i64 delta,int d) const { return d?ceildiv(delta,d):0; }
    i64 center(i64 delta,int t,int i,int copies=1) const {
        return term(delta,s.cent[t][i]*copies,s.H.S_r);
    }
    i64 fy(int family,int u,int d,int copies) const {
        auto knot=interp[d];
        const i64 *a=s.fy[family]+(size_t)u*s.H.nknots+knot.j;
        return ceildiv(wide(sign)*copies*(wide(knot.d-knot.n)*a[0]+wide(knot.n)*a[1]),knot.d);
    }
    i64 conjugate(int family,int u,i64 jump,int copies) const {
        const i64 *a=s.fy[family]+(size_t)u*s.H.nknots;
        // Header bounds Q_K <= 65536, |a| <= 16*S, |jump| <= 32*S and knots <= Q_K/2.
        i64 m=INT64_MIN;
        for(int j=0;j<s.H.nknots;j++)
            m=std::max(m,s.H.knot_num[j]*jump-sign*(i64)s.H.Q_K*a[j]);
        return ceildiv(wide(m)*copies,s.H.Q_K);
    }
    i64 tail_gain(const EdgeState &v,const EdgeState &other,i64 jump) const {
        if(jump>=0) return jump;
        if(v.k>l.KT) return 0;
        int d=other.kind==0?other.D:1;
        int numerator=std::max(0,(l.D0-1-v.k)*d+(other.kind==0));
        return term(jump,numerator,d*(l.D0+1));
    }
    i64 tail_bud(const EdgeState &v,i64 jump) const {
        if(v.k>l.KT) return 0;
        return jump>=0?jump:term(jump,std::max(0,l.D0-1-v.k),l.D0);
    }
    template<class Fn> i64 bud(const EdgeState &v,i64 jump,Mix p,Fn newborn) const {
        wide value=wide(p.den)*jump+wide(p.n)*newborn(p.first)+wide(p.den-p.n)*newborn(p.second);
        return ceildiv(value*v.L,p.den*v.D);
    }
public:
    Low(Snapshot &snap,int num,int den):s(snap),c(s.C),l(c.L2),sign(s.H.mode==1?1:-1),
        P(num),Q(den),nb2(c4_lvl_inner(&l,1,0)),nb3(c4_lvl_inner(&c.L3,1,0)) {
        check(s.H.Q_K<=65536,"knot scale outside fast arithmetic range");
        for(int den=0;den<=l.D0;den++) {
            int inv=den?1:0,d=den?den:1,j=0;
            // The formal denominator 1 is unused (no core has reciprocal 1).
            if(den==1) inv=0;
            while(j+2<s.H.nknots && inv*s.H.Q_K>d*s.H.knot_num[j+1]) ++j;
            interp.push_back({j,inv*s.H.Q_K-d*s.H.knot_num[j],d*(s.H.knot_num[j+1]-s.H.knot_num[j])});
        }
        for(int i=0;i<l.NG;i++) {
            EdgeState v{};
            v.id=i;v.k=l.gK[i];v.L=l.gL[i];v.D=l.gD[i];v.kind=l.gKind[i];v.invden=v.D;
            v.up=l.INF_ID;v.bud=i;v.hu=v.cu=v.hu_bud=v.cu_bud=-1;
            if(v.kind==0) {
                v.up=corner(l,v.L+1,v.k);
                v.hu=c4_hu_idx(&c,v.L,v.k); v.cu=c4_cu_idx(&c,v.L,v.k);
                v.df=signed_f(v.L+1,v.k)-signed_f(v.L,v.k);
                if(v.L>0) {
                    v.bud=corner(l,v.L-1,v.k+1);
                    v.hu_bud=c4_hu_idx(&c,v.L-1,v.k+1);
                    v.cu_bud=c4_cu_idx(&c,v.L-1,v.k+1);
                }
            } else if(v.kind==2 && v.k<=l.KT) v.bud=c4_lvl_tail(&l,v.k+1);
            v.proj=project(l,c.L3,i);
            st.push_back(v);
        }
    }
    i64 vertex(int L,int k) const {
        int D=L+k+1; i64 v=signed_f(L,k);
        Frac next=f(L+1,k); next.n-=wide(v)*next.d;
        i64 total=term(v,-P,Q);
        if(D<l.D0) {
            int u=corner(l,L,k);
            i64 r=(i64)L*s.H.S_r+s.cent[0][u]+k*s.cent[1][u];
            total+=ceildiv(next.n*r,s.H.S_r*next.d);
        } else {
            i64 rate=2*L+(next.n>0?k+1:0);
            total+=ceildiv(next.n*rate,2*next.d);
        }
        if(L) {
            i64 born=signed_f(L-1,k+1)-v+signed_f(1,0);
            int p=corner(l,L-1,k+1);
            if(D<l.D0) total+=term(born+g(p,nb2),L,D);
            else total+=ceildiv(wide(L)*(wide(D)*born+wide(l.D0)*g(p,nb2)+wide(D-l.D0)*g(l.INF_ID,nb2)),(i64)D*D);
        }
        return total;
    }
    i64 edge(int p,int w) const {
        const auto &a=st[p],&b=st[w];
        i64 v=g(p,w),A=g(a.up,w)-v,B=g(p,b.up)-v;
        i64 total=term(v,-P,Q);
        if(a.kind==0) {
            total+=inverse(a.df,b.invden)+center(-a.df,1,p);
            total+=a.L*A;
            if(b.kind==0) total+=inverse(A,b.invden);
            if(a.hu>=0) total+=center(A,2,p);
            else total+=fy(1,p,b.invden,1)+conjugate(0,p,A,1);
            if(a.cu>=0) total+=center(A,4,p,a.k-1);
            else if(a.k>=2) total+=fy(2,p,b.invden,a.k-1)+conjugate(2,p,A,a.k-1);
            if(a.L) total+=bud(a,g(a.bud,w)-v,b.proj,[&](int i){return sign*s.Ch(a.cu_bud,i,nb3);});
        }
        if(b.kind==0) {
            total+=inverse(b.df,a.invden)+center(-b.df,0,w);
            total+=b.L*B;
            if(a.kind==0) total+=inverse(B,a.invden);
            if(b.hu>=0) total+=center(B,3,w,b.k);
            else if(b.k) total+=fy(0,w,a.invden,b.k)+conjugate(1,w,B,b.k);
            if(b.L) total+=bud(b,g(p,b.bud)-v,a.proj,[&](int i){return sign*s.V(i,b.hu_bud,nb3);});
        }
        if(a.kind==2) {
            total+=tail_gain(a,b,A)+tail_bud(a,g(a.bud,w)-v);
            if(b.kind==0) total+=inverse(std::max(B,g(l.INF_ID,b.up)-g(l.INF_ID,w)),l.D0);
        }
        if(b.kind==2) {
            total+=tail_gain(b,a,B)+tail_bud(b,g(p,b.bud)-v);
            if(a.kind==0) total+=inverse(std::max(A,g(a.up,l.INF_ID)-g(p,l.INF_ID)),l.D0);
        }
        return total;
    }
    void scan() const {
        Maximum inner,table,pair;
        for(int D=2;D<s.H.DF;D++) for(int k=0;k<D;k++)
            (D<l.D0?inner:table).add(vertex(D-1-k,k),D,k);
        for(int i=0;i<l.npars;i++) for(int w=0;w<l.NG;w++)
            pair.add(edge(l.pars[i],w),l.pars[i],w);
        inner.print("f_inner",P,Q);table.print("f_table",P,Q);pair.print("g",P,Q);
    }
};
static int integer(const char *p) {
    int i=0; auto end=p+strlen(p);auto ans=std::from_chars(p,end,i);
    check(ans.ec==std::errc{} && ans.ptr==end,"invalid integer");return i;
}
int main(int argc,char **argv) try {
    check(argc==3 || argc==5,"usage: low {scan|probe} FILE [P Q]");
    Snapshot s(argv[2]);
    int P=argc==5?integer(argv[3]):s.H.P,Q=argc==5?integer(argv[4]):s.H.Q;
    check(P>0 && Q>=P && Q<=1000000000,"invalid verification exponent");
    Low low(s,P,Q);
    if(std::string(argv[1])=="scan") low.scan();
    else {
        check(std::string(argv[1])=="probe","unknown command");
        std::string kind;int i,j;
        while(std::cin>>kind>>i>>j) {
            if(kind=="v") {
                check(i>=0 && j>=0 && (i64)i+j+1>=2 && (i64)i+j+1<s.H.DF,"vertex probe range");
                std::cout<<low.vertex(i,j)<<"\n";
            } else {
                check(kind=="g" && i>=0 && i<s.C.L2.NG && s.C.L2.p_to_pi[i]>=0 && j>=0 && j<s.C.L2.NG,"edge probe range");
                std::cout<<low.edge(i,j)<<"\n";
            }
        }
        check(std::cin.eof(),"malformed probe");
    }
    return 0;
} catch(const std::exception &e) { std::cerr<<e.what()<<"\n";return 2; }
