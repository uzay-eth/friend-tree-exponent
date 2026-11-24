// Standalone exhaustive interval checker. No floating-point inputs or operations.
// Fixed-point endpoints are integers in units of 2^-40, with 128-bit intermediates.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <map>
#include <utility>
#include <vector>
using I64=std::int64_t;
using W=__int128_t;
constexpr I64 Q=I64(1)<<40;
struct Overflow {};
I64 checked(W n) { if(!(n>=-(__int128(1)<<62)&&n<(__int128(1)<<62)))throw Overflow{}; return I64(n); }
I64 floor_q(W n,W d) { assert(d>0); return checked(n/d-(n%d<0)); }
I64 ceil_q(W n,W d) { return -floor_q(-n,d); }
struct I {
    I64 l,h;
    I(int n=0):l(I64(n)*Q),h(l) {}
    I(I64 lo,I64 hi):l(lo),h(hi) { assert(lo<=hi); }
    bool pos() const {return l>0;}
    bool neg() const {return h<0;}
};
I rat(int n,int d) {return {floor_q(W(n)*Q,d),ceil_q(W(n)*Q,d)};}
I operator+(I a,I b){return {checked(W(a.l)+b.l),checked(W(a.h)+b.h)};}
I operator-(I a){return {-a.h,-a.l};}
I operator-(I a,I b){return a+-b;}
I operator*(I a,I b){
    auto [l,h]=std::minmax({W(a.l)*b.l,W(a.l)*b.h,W(a.h)*b.l,W(a.h)*b.h});
    return {floor_q(l,Q),ceil_q(h,Q)};
}
I times(int n,I x){assert(n>=0);return {checked(W(n)*x.l),checked(W(n)*x.h)};}
I over(I x,int n){return {floor_q(x.l,n),ceil_q(x.h,n)};}
I operator/(I a,I b){
    assert(b.pos());
    return {floor_q(W(a.l)*Q,a.l<0?b.l:b.h),ceil_q(W(a.h)*Q,a.h>0?b.l:b.h)};
}
I abs(I a){return {a.pos()?a.l:a.neg()?-a.h:0,std::max(-a.l,a.h)};}
I sq(I a){a=abs(a);return {floor_q(W(a.l)*a.l,Q),ceil_q(W(a.h)*a.h,Q)};}
I hull(I a,I b){return {std::min(a.l,b.l),std::max(a.h,b.h)};}
I inflate(I x,I err){assert(err.h>=0);return x+I(-err.h,err.h);}
I center(I a){I64 mid=a.l+(a.h-a.l)/2;return {mid,mid};}
struct Z {I x,y; Z(I r=0,I i=0):x(r),y(i) {}};
Z operator+(Z a,Z b){return {a.x+b.x,a.y+b.y};}
Z operator-(Z a){return {-a.x,-a.y};}
Z operator-(Z a,Z b){return a+-b;}
Z operator*(Z a,Z b){return {a.x*b.x-a.y*b.y,a.x*b.y+a.y*b.x};}
Z operator*(I a,Z b){return {a*b.x,a*b.y};}
Z times(int n,Z b){return {times(n,b.x),times(n,b.y)};}
Z over(Z z,int n){return {over(z.x,n),over(z.y,n)};}
Z inv(Z b){I den=sq(b.x)+sq(b.y);return {b.x/den,-b.y/den};}
Z inflate(Z z,I e){return {inflate(z.x,e),inflate(z.y,e)};}
I norm1(Z z){return abs(z.x)+abs(z.y);}
int bin(int n,int k){if(k>n)return 0;int r=1;for(int j=1;j<=k;++j)r=r*(n-j+1)/j;return r;}
I power(I x,int n){I r=1;while(n--)r=r*x;return r;}

constexpr int D=512, C=36, K=6, P=K+1, A0=20;
constexpr int E=12, A=40;
const I a=rat(9,16); // continuation proved for Re s > -2/3
const I beta=rat(2,3);
const I betaP=rat(3,4);
I rate_scaled(int d,int c){return I(d+1)-times(d,over(a,d-c));}
I forcing(int d,int c){return rat(c,d+1)-rat(1,d);}
// R(s)=b+ integral e^-st E Gb(X_t), so F=(R-s)/s.
I major[D+1][D+1];
I weight(int d,int c){return rat(2*(c+17),2*(d-c)-1);}
I w(int d,int c){return power(rat(c+A0,d),P);}
I ew(int d,int c){return power(rat(c+A,d),E);}
void setup(){
    // Check baseline drift for L=1,...,15.  For larger L use the coarse inequality in the report.
    for(int l=2;l<16;++l){
        I gm=-rat(2*l,2*l+1);
        I correction=rat(l,D)*(rat(2,2*l-3)*(I(1)+rat(1,18))+rat(1,18));
        assert((gm+correction+beta).h<=0);
    }
    assert(( -rat(32,33)+rat(2,29)*(I(1)+rat(1,18))+rat(1,18)+beta).neg());
    // Strong drift before exiting c/d<=2/3 at d>=D.
    I t=power(rat(A0+1,A0),P);
    assert((rat(P,3)-rat(P*(P+1)/2,D)-(t-I(1))-betaP).pos());
    I te=power(rat(A+1,A),E);
    assert((rat(E,3)-rat(E*(E+1)/2,D)-(te-I(1))-betaP).pos());
    // Full DP majorizes |integral e^-st E Gb|, with W/(beta-a) at d=D.
    for(int c=1;c<D;++c) major[D][c]=weight(D,c)/(beta-a);
    for(int d=D-1;d>=2;--d) for(int c=d-1;c>=1;--c){
        I den=rate_scaled(d,c);
        major[d][c]=(abs(forcing(d,c))+major[d][c+1])/den+(I(d)/den)*major[d+1][c];
    }
    std::cout<<"global majorant upper endpoint (units 2^-40): "<<major[2][1].h<<std::endl;
}

// q at the degree boundary for R=sH. Polynomial coefficient vectors in c.
using Poly=std::array<Z,K+2>;
Poly delta(const Poly& p){Poly r{};for(int j=1;j<K+2;++j)for(int i=0;i<j;++i)r[i]=r[i]+times(bin(j,i),p[j]);return r;}
std::array<Poly,K+1> polys(Z s){
    std::array<Poly,K+1> q{};q[0][0]=Z(1);
    for(int n=1;n<=K;++n){
        Poly before{};if(n==1)before[1]=s;
        for(int i=1;i<n;++i){
            int j=n-i+1;
            int u=(j%2?-1:1)*bin(i+j-1,j), v=(j%2?1:-1)*bin(i+j-2,j-1);
            Poly dq=delta(q[i]);
            for(int k=0;k<=i;++k){
                before[k]=before[k]-Z(u)*q[i][k];
                before[k+1]=before[k+1]+Z(v)*q[i][k];
                if(i+1==n)before[k+1]=before[k+1]+dq[k];
            }
        }
        Z div=inv(s+Z(n)), invpow=div;
        for(int j=0;j<=n;++j){
            for(int k=0;k<=n-j;++k)q[n][k]=q[n][k]-invpow*before[k];
            if(j<n){invpow=invpow*div;before=delta(before);}
        }
    }
    return q;
}
I residual(const std::array<Poly,K+1>& q,Z s){
    // Exact cancellations in degrees <=K are proved algebraically in the report.
    std::array<Poly,2*K+2> r{};
    r[K+1][1]=s;
    for(int i=1;i<=K;++i){
        Poly dq=delta(q[i]);
        for(int j=0;j<=K;++j){
            int kj=bin(K,j), gj=bin(K-i,j)-kj;
            for(int k=0;k<=i;++k){
                if(i+j>=P){r[i+j][k]=r[i+j][k]+times(kj,s*q[i][k]-dq[k]);
                           r[i+j][k+1]=r[i+j][k+1]+Z(gj)*q[i][k];}
                if(i+j-1>=P)r[i+j-1][k]=r[i+j-1][k]-Z(gj)*q[i][k];
                if(i+j+1>=P)r[i+j+1][k+1]=r[i+j+1][k+1]+times(kj,dq[k]);
            }
        }
    }
    I M=0;
    for(int j=P;j<=2*K+1;++j)for(int k=0;k<=P;++k){
        I v=norm1(r[j][k]);
        for(int u=P;u<j;++u)v=over(v,D);
        for(int u=k;u<P;++u)v=over(v,A0);
        M=M+v;
    }
    return M;
}

// For |Im(s)| high the baseline global bound handles the unbounded region.
// On a bounded contour use midpoint evaluation and enclose the nearby curve
// by difference-quotient propagation, with rational modulus bounds on steps.
struct Eval {Z value; I error; I slope;};
I offBoundary(const std::array<Poly,K+1>& q){
    // On exit c/d>2/3 (at a c increment) c/d <= 2/3+1/D.
    // Bound W, and the exit Lyapunov function is >= (2/3)^E.
    I value=I(1)+((rat(2,3)+rat(18,D))/(rat(1,3)-rat(2,D)))/(beta-a);
    for(int i=0;i<=K;++i)for(int j=0;j<=i;++j){
        I u=norm1(q[i][j])*power(rat(2,3)+rat(1,D),j);
        for(int k=j;k<i;++k)u=over(u,D);
        value=value+u;
    }
    return value / power(rat(2,3),E);
}

// Evaluate R-s, with the finite truncation error propagated through a
// point/segment-specific *complex* discount modulus bound (|z| >= Re z also).
Eval evaluate(Z s,Z segment){
    assert((segment.x+a).l>=0);
    auto q=polys(s); auto boxq=polys(segment);
    I tail=residual(boxq,segment)/(betaP-a), exit=offBoundary(boxq);
    // Approximate R via its directly forced recursion, (s-G)R=s b.
    std::array<Z,C+1> nxt{};
    std::array<I,C+1> budget{}, slope{}, rnd{};
    // Bound degree-boundary derivative of q with Cauchy's estimate on a
    // box of radius 1/8 (the poles of q are at negative integers only).
    const I eps=rat(1,8);
    auto wideq=polys(inflate(segment,eps));
    for(int c=1;c<C;++c){
        for(int i=K;i>=0;--i){
            for(int k=0;k<=i;++k){
                Z v=power(rat(c,D),k)*q[i][k];
                for(int j=k;j<i;++j)v=over(v,D);
                nxt[c]=nxt[c]+v;
            }
        }
        Z mid(center(nxt[c].x),center(nxt[c].y));
        rnd[c]=norm1(nxt[c]-mid);nxt[c]=mid;
        budget[c]=tail*w(D,c)+exit*ew(D,c);
        for(int i=1;i<=K;++i)for(int k=0;k<=i;++k){
            I u=norm1(wideq[i][k])*power(rat(c,D),k);
            for(int j=k;j<i;++j)u=over(u,D);
            slope[c]=slope[c]+times(8,u);
        }
    }
    for(int d=D-1;d>=2;--d){
        std::array<Z,C+1> row{};
        std::array<I,C+1> bb{},ss{},rr{};
        if(C<d)bb[C]=I(1)+major[d][C];
        for(int c=std::min(d-1,C-1);c>=1;--c){
            int l=d-c;
            Z denom=Z(d+1)+over(times(d,s),l);
            Z calc=(s+times(d,nxt[c])+row[c+1])*over(inv(over(denom,d+1)),d+1);
            row[c]=Z(center(calc.x),center(calc.y));
            Z boxden=Z(d+1)+over(times(d,segment),l);
            // Rational lower bound on |z|: max(x, |y|, (x+|y|)*7/10).
            I ay=abs(boxden.y);
            I m=rat(7,10)*(boxden.x+ay);
            I64 lb=std::max({boxden.x.l,ay.l,m.l});
            I md(lb,lb);
            I gain=I(d)/md;
            bb[c]=gain*budget[c]+bb[c+1]/md;
            rr[c]=gain*rnd[c]+rr[c+1]/md+norm1(calc-row[c]);
            // Difference quotient, subtracting the midpoint recurrence
            // from the recurrence at an arbitrary point on the segment.
            I fac=rat(d,l);
            ss[c]=gain*slope[c]+(I(1)+ss[c+1]+fac*(norm1(row[c])+rr[c]))/md;
        }
        nxt=row;budget=bb;slope=ss;rnd=rr;
    }
    return {nxt[1]-s,budget[1]+rnd[1],slope[1]+I(1)};
}

using Pt=std::pair<I64,I64>;
Z point(Pt p){return {I(p.first,p.first),I(p.second,p.second)};}
Z box(Pt a,Pt b){Z x=point(a),y=point(b);return {hull(x.x,y.x),hull(x.y,y.y)};}
bool misses(Z z){return z.x.pos()||z.x.neg()||z.y.pos()||z.y.neg();}
I special_factor_x=rat(1,4),special_factor_y=rat(3,8);
Z planted(Z s){return s-Z(special_factor_x,special_factor_y);}
void krawczyk(Z z){
    // Complex disk Krawczyk via Banach, with rational upper bounds on moduli.
    // Cauchy on |u|=R bounds all nonlinear Taylor coefficients of G0(z+u).
    const I R=rat(1,128), h=rat(1,2048), r=rat(1,131072);
    Eval e=evaluate(z,inflate(z,R));
    I M=e.error+times(2,R)*e.slope;
    Eval p=evaluate(z+Z(h),z+Z(h)), m=evaluate(z-Z(h),z-Z(h));
    Z derivative=times(1024,p.value-m.value);
    Z Ainv=inv(derivative);
    I hR=h/R, rR=r/R;
    I error=(p.error+m.error)/(times(2,h))+(M/R)*sq(hR)/(I(1)-sq(hR));
    error=error+(M/R)*rR*(I(2)-rR)/sq(I(1)-rR);
    I contraction=norm1(Z(1)-Ainv*derivative)+norm1(Ainv)*error;
    I displacement=norm1(Ainv*e.value)+norm1(Ainv)*evaluate(z,z).error;
    I image=displacement+r*contraction;
    std::cout<<"Krawczyk disk contraction, image radius upper (units 2^-40): "<<contraction.h<<" "<<image.h<<std::endl;
    assert((I(1)-contraction).pos()&&(r-image).pos());
}
struct Contour {
    bool mutant;
    int segments=0;
    std::vector<Pt> polygon;
    std::map<Pt,Pt> cache;
    Pt endpoint(Pt p){
        auto it=cache.find(p); if(it!=cache.end())return it->second;
        Z s=point(p);Eval e=evaluate(s,s); Z v=e.value;
        if(mutant)v=v*planted(s);
        Pt rep={v.x.l+(v.x.h-v.x.l)/2,v.y.l+(v.y.h-v.y.l)/2};
        return cache.emplace(p,rep).first->second;
    }
    void edge(Pt p,Pt q,int depth=0){
        Pt mid={(p.first+q.first)/2,(p.second+q.second)/2};
        Z z=point(mid), seg=box(p,q);
        I radius=norm1(seg-z);
        if(radius.h<=Q/4)try {
            Eval e=evaluate(z,seg);
            Z enclosure=inflate(e.value,e.error+radius*e.slope);
            if(mutant)enclosure=enclosure*planted(seg);
            Pt pp=endpoint(p), qq=endpoint(q);
            Z ep=box(pp,qq);
            enclosure={hull(enclosure.x,ep.x),hull(enclosure.y,ep.y)};
            if(misses(enclosure)){polygon.push_back(pp);++segments;return;}
        }catch(const Overflow&){} // Coarse candidate enclosure; subdivide.
        if(depth==36){Eval e=evaluate(z,seg);std::cerr<<"unresolved segment, midpoint units 2^-40: "<<mid.first<<" "<<mid.second<<"; image "<<e.value.x.l<<" "<<e.value.y.l<<" error "<<e.error.h<<std::endl;std::abort();}
        assert(p!=mid&&q!=mid);
        edge(p,mid,depth+1);edge(mid,q,depth+1);
    }
    int rect(I64 l,I64 r,I64 b,I64 t){
        edge({l,b},{r,b});edge({r,b},{r,t});edge({r,t},{l,t});edge({l,t},{l,b});
        int winding=0;
        for(size_t i=0;i<polygon.size();++i){
            auto [x,y]=polygon[i];auto [u,v]=polygon[(i+1)%polygon.size()];
            W cp=W(x)*v-W(y)*u;
            if(y<=0&&v>0&&cp>0)++winding;
            if(y>0&&v<=0&&cp<0)--winding;
        }
        return winding;
    }
};
int main(int argc,char** argv){
    setup();
    int Y=1;while(!((I(Y)-rat(1,2)-major[2][1]).pos()))++Y;
    std::cout<<"high frequency cutoff "<<Y<<std::endl;
    int mode=argc>1?std::atoi(argv[1]):0;
    const I64 L=-9*Q/16;
    auto run=[](const char* label,I64 l,I64 r,I64 b,I64 t,int expected,bool mutant=false){
        Contour c{mutant}; int w=c.rect(l,r,b,t);
        std::cout<<label<<": winding="<<w<<", segments="<<c.segments<<std::endl;
        assert(w==expected);
    };
    if(mode==0||mode==1)run("right half-plane Re >= -1/2",-Q/2,Q,-Y*Q,Y*Q,1);
    if(mode==0||mode==2)run("right half-plane Re >= -9/16",L,Q,-Y*Q,Y*Q,3);
    if(mode==0||mode==3){
        constexpr I64 X=W(-564548)*Q/1024000, V=W(513564)*Q/1024000;
        Z z=point({X,V}); Eval e=evaluate(z,z);
        std::cout<<"scan point, midpoint residual and error units 2^-40: "<<e.value.x.l<<" "<<e.value.y.l<<" "<<e.error.h<<std::endl;
        run("upper root",X-Q/131072,X+Q/131072,V-Q/131072,V+Q/131072,1);
        krawczyk(z);
    }
    if(mode==0||mode==4)run("planted extra zero",-Q/2,Q,-Y*Q,Y*Q,2,true);
    if(mode==0||mode==5){
        Eval e=evaluate(Z(0),Z(0));Z v=inflate(e.value,e.error);
        assert(v.x.pos());
        std::cout<<"pole residue real endpoints, units 2^-40: "<<v.x.l<<" "<<v.x.h<<std::endl;
    }
    if(mode==0||mode==6){
        I64 num=argc>2?std::strtoll(argv[2],nullptr,10):551315;
        assert(num>=500000&&num<=562500);
        I64 left=floor_q(-W(num)*Q,1000000);
        Contour c{false}; int n=c.rect(left,Q,-Y*Q,Y*Q);
        std::cout<<"right half-plane Re >= -"<<num<<"/1000000 (left endpoint "<<left<<" * 2^-40): winding="<<n<<", segments="<<c.segments<<std::endl;
        assert(n==1);
    }
}
