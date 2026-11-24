// Exact upward integer Bellman solve. See README for state and ordering.
#include <algorithm>
#include <array>
#include <cassert>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <vector>
using U=uint64_t;
using I=__int128_t;
constexpr U Q=10000000000ULL, den=1000000, sn=565658;
using S=std::array<int,6>;
U ceilq(I a,I b) { assert(a>=0&&b>0); return (a+b-1)/b; }
U exp_up(U x) {
    // x is an upper argument in Q units; x<1. Taylor terms through 12,
    // and at most twice the first omitted term.
    assert(x<Q);
    U term=Q,sum=Q;
    for(int j=1;j<=12;j++){term=ceilq(I(term)*x,I(Q)*j);sum+=term;}
    return sum+2*ceilq(I(term)*x,I(Q)*13);
}

std::vector<std::vector<U>> coupon(int B) {
    constexpr U steps=64;
    U disc=ceilq(I(Q)*Q,Q+I(sn)*Q/(den*steps)); // exp(-s/steps) <= 1/(1+s/steps)
    std::vector<std::vector<U>> a(B+2);
    for(int d=2;d<=B+1;d++) {
        U grow=exp_up(ceilq(I(Q)*(d+1),steps*d));
        a[d].assign(d+1,0);
        for(int m=1;m<=d;m++) {
            U age=Q, toll=ceilq(2*I(m)*Q,d);
            I total=0;
            while(toll<2*Q) {
                toll=ceilq(I(toll)*grow,Q);
                total+=I(age)*std::min(toll,2*Q);
                age=ceilq(I(age)*disc,Q);
            }
            a[d][m]=ceilq(total,I(Q)*steps)+ceilq(2*I(age)*den,sn);
        }
    }
    return a;
}
U choose[300][7];
U rank(S a) {
    // Lexicographic among compositions of a fixed total, with e as last index.
    int t=0; for(int x:a)t+=x;
    U r=0;
    for (int j=0;j<5;j++) {
        int k=5-j;
        r+=choose[t+k][k]-choose[t-a[j]+k][k];
        t-=a[j];
    }
    return r;
}
// Storage order l,u,v,c,e,k; displayed order is u,v,c,e,l,k.
enum {l,u,v,c,e,k};
struct Edge { U rate,val; };
struct Mix {U rate=0; I rhs=0; void add(U r,U w) {rate+=r;rhs+=I(r)*w;} };

struct Options {std::array<Mix,3> x; unsigned n;};
template<unsigned K>
U solve(Mix fixed, const std::array<Options,K>& opts,U reward,U dm, bool check_all) {
    fixed.rhs+=I(reward)*Q;
    auto ratio=[](Mix z,U d) {return ceilq(z.rhs*den, I(z.rate)*den+I(d)*sn);};
    U y=0; // Omitting a mandatory LC group is not a feasible policy.
    for(;;) {
        Mix z=fixed;
        for (const auto& group:opts) {
            const Mix* best=&group.x[0];
            for(unsigned j=1;j<group.n;j++) {
                const auto& a=group.x[j];
                if(a.rhs-I(a.rate)*y>best->rhs-I(best->rate)*y)best=&a;
            }
            z.rate+=best->rate;z.rhs+=best->rhs;
        }
        U next=ratio(z,dm);
        if(next<=y) break;
        y=next;
    }
    // Independent rule: re-enumerate endpoint products and check every maximum.
    if(check_all) {
        auto scan=[&](auto&& self,unsigned i,Mix z)->void {
            if(i==opts.size()) {if(ratio(z,dm)>y)std::abort();return;}
            for(unsigned j=0;j<opts[i].n;j++) {auto a=opts[i].x[j];self(self,i+1,{z.rate+a.rate,z.rhs+a.rhs});}
        };
        scan(scan,0,fixed);
    }
    return y;
}

std::vector<std::vector<U>> side(int B) {
    std::vector<std::vector<U>> A(B+2);
    for(int d=B+1;d>=2;--d) {
        A[d].resize(d);
        for(int j=d-1;j>=0;--j) {
            if(d==B+1) {
                A[d][j]=std::min(ceilq(I(Q)*den,sn),ceilq(I(Q)*den*(I(j)*sn+den),I(sn)*(den+I(sn)*d)));
                continue;
            }
            Mix base,more;
            if(d-1-j) base.add(2*(d-1-j),A[d][j+1]);
            base.add(2*d*(d-1-j),A[d+1][j]);
            more.add(d*j,A[d+1][j]);
            A[d][j]=solve<1>(base,{{{{Mix{},more,Mix{}},2}}},2*j,2*d,false);
        }
    }
    return A;
}

int main(int argc,char**argv) {
    if(argc<2) return 2;
    int cap=std::atoi(argv[1]); bool check=argc>=3;
    if(cap<4||cap>150) return 2;
    for(int n=0;n<300;n++) {choose[n][0]=1;for(int j=1;j<=6&&j<=n;j++)choose[n][j]=choose[n-1][j-1]+choose[n-1][j];}
    auto one=side(cap);
    auto cp=coupon(std::max(cap,128));
    U twice=ceilq(I(Q)*2*den,sn);
    std::vector<U> prev,cur;
    U rows=0;
    U maximum=0;int argd=0,argl=0;
    for(int d=cap;d>=3;--d) {
        int t=d-1;
        cur.assign(choose[t+5][5],0);
        U at=0;
        for(int L=0;L<=t;L++)for(int U0=0;U0<=t-L;U0++)for(int V=0;V<=t-L-U0;V++)
        for(int C=0;C<=t-L-U0-V;C++)for(int E=0;E<=t-L-U0-V-C;E++) {
            int K=t-L-U0-V-C-E;
            S a{L,U0,V,C,E,K};
            int b=d-L-K;
            auto term=[&](const S& z) {
                int db=1+z[u]+z[v]+z[c]+z[e],da=db+z[l]+z[k];
                return std::min({twice,one[db][z[v]+z[c]+z[e]]+one[da][z[c]+z[e]+z[k]],
                                 cp[db][std::min(db,z[v]+z[e]+z[l]+z[k])]});
            };
            if(b<2){++at;continue;}
            auto get=[&](S z) {
                if(z==a) return U(0);
                int tot=0;for(int x:z)tot+=x;
                if(tot==t) {auto i=rank(z);assert(i<at);return cur[i];}
                assert(tot==t+1);
                return d==cap?term(z):prev[rank(z)];
            };
            auto ch=a;ch[u]++;U common=get(ch);ch=a;ch[l]++;U priv=get(ch);
            U m=2*d*b;
            Mix f;
            f.add(U0*m,common);f.add((L+V)*m,priv);
            // LC source budget uses common b_i and private (1-b_i).
            Mix lc0,lc1;lc1.add(V*m/2,common);
            // Remove LC baseline from f to express shared alternatives positively.
            f.rate-=V*m;f.rhs-=I(V*m)*priv;
            lc0.add(V*m,priv);lc1.add(V*m/2,priv);
            Mix cc,ep,ec,kp;
            cc.add(C*m/2,common);ep.add(E*m/2,priv);ec.add(E*m/2,common);kp.add(K*m/2,priv);
            if(U0) {ch=a;ch[u]--;ch[c]++;f.add(2*b*U0,get(ch));}
            if(V) {ch=a;ch[v]--;ch[e]++;f.add(2*b*V,get(ch));}
            for(int A : {l,k}) for(int B : std::array<int,5>{u,c,e,v,-1}) {
                int na=a[A],nb=B<0?1:a[B];
                if(!na||!nb)continue;
                ch=a;
                if(A==l){ch[l]--;ch[k]++;}
                if(B==u){ch[u]--;ch[v]++;}
                if(B==c){ch[c]--;ch[e]++;}
                if(ch!=a) f.add(2*na*nb,get(ch));
            }
            U reward=2*(C*(d-b)+E*(d+b)+V*d+K*b);
            cur[at++]=solve<4>(f,{{{{lc0,lc1,Mix{}},2},{{Mix{},cc,Mix{}},2},
                                 {{Mix{},ep,ec},3},{{Mix{},kp,Mix{}},2}}},reward,m,check);
            ++rows;
            if(L==1&&!K&&!V&&!E) {
                U hit=std::min(cur[at-1],cp[b][1]);
                if(hit>maximum){maximum=hit;argd=b;argl=U0;}
                std::cout<<b<<' '<<U0<<' '<<hit<<'\n';
            }
        }
        assert(at==cur.size());
        std::cerr<<"slice DA="<<d<<" states="<<cur.size()<<" cumulative="<<rows<<'\n';
        prev.swap(cur);
    }
    std::cerr<<"coupon d=32,40,48,64,128 "<<cp[32][1]<<' '<<cp[40][1]<<' '<<cp[48][1]<<' '<<cp[64][1]<<' '<<cp[128][1]<<'\n';
    std::cerr<<"max initial finite "<<maximum<<" at "<<argd<<' '<<argl<<"; tail "<<cp[cap][1]<<'\n';
    if(cap==48&&(maximum>=Q||cp[cap][1]>=Q))return 3;
}
