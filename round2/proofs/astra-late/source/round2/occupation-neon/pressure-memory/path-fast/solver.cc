#include <algorithm>
#include <array>
#include <cstdint>
#include <cstdio>
#include <exception>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <vector>

namespace {
using I = std::int64_t;
using W = __int128_t;
using Value = std::array<I,4>;
using Wire = int (*)(const int*, const int*, I*);
constexpr I limit = I(1) << 60;
void check(bool ok, const char* why) { if (!ok) throw std::runtime_error(why); }
I narrow(W x) {
    check(-limit < x && x < limit, "fixed-point range exceeded");
    return I(x);
}
// C++ integer division truncates toward zero; explicitly round outwards.
I ratio(W a, I b, bool upper) {
    check(b > 0, "nonpositive denominator");
    auto q = a/b, r = a%b;
    return narrow(q + (upper ? int(r > 0) : -int(r < 0)));
}
struct Census { int d, k; };
constexpr Census virgin{2,0};
struct State {
    std::array<Census,4> path{};
    int n = 1;
    Census y{-1,-1};   // absent, then finite degrees, then overflow
    int m = 0;         // finite degrees, then overflow
};
struct Reward {
    int status = 0;
    Value minimum{};
    std::array<Value,5> source{}, slope{};
};
#include "charges.hh"
struct Edge { I lo, hi; Value value; };
    struct Solver {
    std::array<int,4> cap, arity;
    std::array<std::array<std::array<int,16>,17>,4> code{};
    std::array<std::size_t,5> stride{1};
    struct Block {
        std::array<std::array<int,17*16>,4> digit{};
        std::array<std::size_t,5> stride{1};
        std::size_t offset{};
    };
    std::array<std::array<Block,2>,5> blocks{};
    int side, mother, ym, mm;
    I unit = 1000000, s0, ds;
    std::array<I,26> r{};
    const Wire load;
    const Charges charges;
    // Slab uses an independent rank for complete path censuses only.
    std::vector<Reward> rewards;
    int part = -2;     // -2 minimum, -1 mother, 0..3 donor attached at that role
    std::vector<Value> value;
    std::vector<unsigned char> seen;
    std::vector<unsigned char> paid;
    std::uint64_t calls{}, paths{}, edges{}, variable{}, ratios{}, budget{}, hash{};
    I min_value{}, max_value{};

    Solver(const int* c, const int* a, int sc, int mc, const Tables* t,bool exclude_old,Wire wire): side(sc), mother(mc),
        ym(2 + (sc ? sc*(sc+1)/2-1 : 0)), mm(std::max(1,mc)), load(wire), charges(t,exclude_old) {
        int top = std::max({2, side+1, mother+1});
        for (int i=0;i<4;++i) {
            cap[i]=c[i]; arity[i]=a[i];
            check(c[i]>=2 && c[i]<=16 && a[i]>=int(i<3) && a[i]<c[i], "invalid path caps");
            top=std::max(top,c[i]);
            int digit=0;
            for (int d=2;d<=c[i];++d)
                for (int k=int(i<3);k<=a[i] && k<d;++k) code[i][d][k]=digit++;
            stride[i+1]=stride[i]*digit;
        }
        check(sc>=0 && sc!=1 && mc>=0 && mc!=1 && top<=24, "invalid memory caps");
        for (int d=2;d<=top;++d) unit=std::lcm(unit,I(d));
        // At most 20 competing arrows, each at most 32*unit. With the
        // narrowed values below this bounds every widened product and sum.
        check(unit <= (I(1)<<48), "rate unit out of range");
        for (int d=2;d<=top;++d) { check(unit%d==0,"inexact rate scale"); r[d]=unit/d; }
        s0=(unit/1000000)*565673;
        ds=(unit/1000000)*102;
        rewards.resize(stride[4]);
    }

    bool live(const State& s) const {
        for (int i=0;i<s.n;++i) {
            auto [d,k]=s.path[i];
            int demand=int(s.n<4 && i==s.n-1)+int(i==part && s.y.d<0);
            if (d>cap[i] || k+demand>arity[i]) return false;
            check(d>=2 && k>=int(i+1<s.n)+int(i==part && s.y.d>=0) && k<d,
                  "invalid census or missing registered child");
        }
        return true;
    }
    std::size_t path_id(const State& s) const {
        std::size_t id=0;
        for (int i=0;i<s.n;++i) id += stride[i]*code[i][s.path[i].d][s.path[i].k];
        return id;
    }
    std::size_t rank(const State& s) const {
        int y = s.y.d<0 ? 0 : s.y.d==0 ? ym-1 : s.y.d*(s.y.d-1)/2+s.y.k;
        int m = s.m ? s.m-2 : mm-1;
        check(y>=0 && y<ym && m>=0 && m<mm,"invalid exterior rank digit");
        const auto& b=blocks[s.n][y!=0];
        std::size_t id=0;
        for (int i=0;i<s.n;++i) {
            auto [d,k]=s.path[i];
            id+=b.stride[i]*b.digit[i][d*16+k];
        }
        return b.offset+(id*(y?ym-1:1)+(y?y-1:0))*mm+m;
    }
    const Reward& reward(const State& s, std::size_t id) {
        auto& b=rewards[id];
        if (!b.status) {
            b.status=charges(s.path,b);
            if (load) {
                int d[4],k[4];
                for (int i=0;i<4;++i) { d[i]=s.path[i].d; k[i]=s.path[i].k; }
                I raw[44]{};
                int status=load(d,k,raw);
                check(status==b.status,"native/Python reward membership mismatch");
                if (status==1) {
                    if (!std::equal(b.minimum.begin(),b.minimum.end(),raw)) {
                        for (auto c:s.path) std::fprintf(stderr," %d/%d",c.d,c.k);
                        for (int j=0;j<4;++j) std::fprintf(stderr," %lld/%lld",(long long)b.minimum[j],(long long)raw[j]);
                        std::fprintf(stderr,"\n");
                        check(false,"native/Python minimum mismatch");
                    }
                    for (int i=0;i<5;++i) {
                        check(std::equal(b.source[i].begin(),b.source[i].end(),raw+4+8*i),"native/Python source mismatch");
                        check(std::equal(b.slope[i].begin(),b.slope[i].end(),raw+8+8*i),"native/Python slope mismatch");
                    }
                }
            }
            ++request;
            supported+=b.status==1;
        }
        return b;
    }

    Value solve(const State& s, std::size_t prev=0, bool entry=false) {
        if (!live(s)) return {};
        auto index=rank(s);
        check(index<seen.size() && (entry || index>prev),"non-increasing mixed-radix DAG edge");
        if (seen[index]) return value[index];
        check(++calls<=budget,"live-state budget exceeded");
        Value source{}, slope{};
        const bool tagged=s.y.d>=0;
        if (s.n==4 && (part<0 || tagged)) {
            auto id=path_id(s);
            const auto& b=reward(s,id);
            if (b.status==1) {
                if (!paid[id]) { paid[id]=1; ++paths; }
                source=part==-2 ? b.minimum : b.source[part+1];
                if (part!=-2) slope=b.slope[part+1];
            }
        }
        std::array<Edge,20> e;
        int ne=0;
        auto add=[&](I lo,I hi,Value z,const Value& w=Value{}) {
            check(lo>=0 && hi>=lo && hi<=32*unit,"arrow rate outside range");
            if (!hi) return;
            check(ne<int(e.size()),"too many arrows");
            for (int j=0;j<4;++j) z[j]=narrow(W(z[j])+w[j]);
            e[ne++]={lo,hi,z};
        };
        for (int i=0;i<s.n;++i) {
            auto [d,k]=s.path[i];
            int leaves=d-k-1;
            State to=s;
            ++to.path[i].d;
            Value grown=solve(to,index);
            I fixed=leaves*unit+(i?r[s.path[i-1].d]:0)+(i+1<s.n?r[s.path[i+1].d]:0);
            I high=fixed;
            if (i==0) {
                I lo=r[s.m], hi=s.m?lo:r[std::max(2,mother+1)];
                if (part==-1 && slope!=Value{}) add(lo,hi,grown,slope);
                else { fixed+=lo; high+=hi; }
            }
            bool known=tagged && i==part;
            high+=(k-int(i+1<s.n)-int(known))*r[2];
            if (known) {
                I lo=r[s.y.d], hi=s.y.d?lo:r[std::max(2,side+1)];
                if (slope!=Value{}) add(lo,hi,grown,slope);
                else { fixed+=lo; high+=hi; }
            }
            add(fixed,high,grown);
            if (leaves) {
                to=s; ++to.path[i].k;
                auto combined=solve(to,index);
                auto accumulate=[&](const State& more) {
                    auto z=solve(more,index);
                    for (int j=0;j<4;++j) combined[j]=narrow(W(combined[j])+z[j]);
                };
                // One physical jump can register both alternatives. Never
                // emit the P4 path and donor simultaneously as the same child.
                if (i==s.n-1 && s.n<4) {
                    State longer=to; longer.path[longer.n++]=virgin; accumulate(longer);
                }
                if (i==part && !tagged) { to.y=side?virgin:Census{0,0}; accumulate(to); }
                add(leaves*r[d],leaves*r[d],combined);
            }
        }
        if (s.y.d>0) {
            auto [d,k]=s.y;
            int leaves=d-k-1;
            State to=s; to.y=d<side?Census{d+1,k}:Census{0,0};
            I fixed=leaves*unit+r[s.path[part].d];
            add(fixed,fixed+k*r[2],solve(to,index));
            if (leaves) { to.y={d,k+1}; add(leaves*r[d],leaves*r[d],solve(to,index)); }
        }
        if (s.m) {
            State to=s; to.m=s.m<mother?s.m+1:0;
            I child=r[s.path[0].d];
            add(child,child+s.m*unit-3*r[2],solve(to,index));
        }

        edges+=ne;
        I baseline=s0;
        for (int i=0;i<ne;++i) { baseline+=e[i].lo; variable+=e[i].lo!=e[i].hi; }
        Value out{};
        for (int j=0;j<4;++j) {
            bool upper=j%2;
            I den=baseline;
            W num=W(source[j])*unit;
            using Choice=std::pair<I,I>;
            std::array<Choice,21> choices;
            int nc=1;
            choices[0]={0,ds};
            for (int i=0;i<ne;++i) {
                const auto& a=e[i];
                num+=W(a.lo)*a.value[j];
                if (a.hi>a.lo) choices[nc++]={a.value[j],a.hi-a.lo};
            }
            std::sort(choices.begin(),choices.begin()+nc,[upper](Choice a,Choice b) {
                return upper?a.first>b.first:a.first<b.first;
            });
            for (int i=0;i<nc;++i) {
                auto [v,inc]=choices[i];
                auto cmp=W(v)*den;
                if (upper?cmp<=num:cmp>=num) break;
                num+=W(v)*inc; den+=inc;
            }
            I v=out[j]=ratio(num,den,upper);
            ++ratios;
            W residual=W(source[j])*unit-W(s0+(((v<0)==upper)?ds:0))*v;
            for (int i=0;i<ne;++i) {
                const auto& a=e[i];
                I delta=a.value[j]-v;
                residual+=W((delta>0)==upper?a.hi:a.lo)*delta;
            }
            check(upper?residual<=0:residual>=0,"signed drift check failed");
            min_value=std::min(min_value,v); max_value=std::max(max_value,v);
        }
        check(out[0]<=out[1] && out[2]<=out[3],"crossed bounds");
        // FNV-1a-style word checksum of ranked reachable states and endpoints.
        auto mix=[this](std::uint64_t x) { hash=(hash^x)*1099511628211ULL; };
        mix(index);
        for (auto x:out) mix(std::uint64_t(x));
        value[index]=out; seen[index]=1;
        return out;
    }

    void run(int component, std::uint64_t bound, I* entries, std::uint64_t* stats, I* range) {
        part=component; budget=bound;
        check(part>=-2 && part<=3,"invalid component");
        std::size_t slots=0;
        for (int n=1;n<=4;++n) for (int tag=0;tag<=int(part>=0);++tag) {
            auto& b=blocks[n][tag];
            b.offset=slots;
            if (tag && part>=n) continue;
            for (int i=0;i<n;++i) {
                int min_k=int(i+1<n)+int(i==part && tag);
                int max_k=arity[i]-int(n<4 && i==n-1)-int(i==part && !tag);
                int count=0;
                for (int d=2;d<=cap[i];++d) for (int k=min_k;k<d && k<=max_k;++k)
                    b.digit[i][d*16+k]=count++;
                b.stride[i+1]=b.stride[i]*count;
            }
            slots+=b.stride[n]*(tag?ym-1:1)*mm;
        }
        check(slots<=160000000,"compact slot budget exceeded");
        value.resize(slots); seen.assign(slots,0); paid.assign(stride[4],0);
        calls=paths=edges=variable=ratios=0;
        min_value=limit; max_value=-limit; hash=14695981039346656037ULL;
        State s;
        s.path[0]=virgin;
        for (int i=0;i<mm;++i) {
            s.m=i?i+1:0;
            auto z=solve(s,0,true);
            std::copy(z.begin(),z.end(),entries+4*i);
        }
        check(calls && ratios==4*calls,"missing range audit");
        auto stats_out=std::array<std::uint64_t,9>{paths,calls,edges,variable,ratios,slots,hash,request,supported};
        std::copy(stats_out.begin(),stats_out.end(),stats);
        range[0]=min_value; range[1]=max_value;
    }
    std::uint64_t request=0,supported=0;
};
} // namespace

extern "C" void* prepare(const int* cap,const int* arity,int side,int mother,const Tables* t,int exclude_old,Wire load) {
    try { return new Solver(cap,arity,side,mother,t,exclude_old,load); }
    catch (const std::exception& e) { std::fprintf(stderr,"path-fast: %s\n",e.what()); return nullptr; }
}
extern "C" int price(void* handle,int part,std::uint64_t budget,I* entries,std::uint64_t* stats,I* range) {
    try { static_cast<Solver*>(handle)->run(part,budget,entries,stats,range); return 0; }
    catch (const std::exception& e) { std::fprintf(stderr,"path-fast: %s\n",e.what()); return 1; }
}
extern "C" void release(void* handle) { delete static_cast<Solver*>(handle); }
