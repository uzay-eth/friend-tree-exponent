#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <iostream>
#include <random>
#include <vector>

struct Vertex { uint64_t leaves; double time; std::vector<size_t> core; };
struct Sim {
    std::mt19937_64 rng;
    std::vector<Vertex> v;
    std::vector<int64_t> bit;
    double next;
    double u() { return (double(rng() >> 12)+0.5)*0x1p-52; }
    double e() { return -std::log(u()); }
    uint64_t idx(uint64_t n) {
        __uint128_t x=__uint128_t(rng())*n;
        while (uint64_t(x)<n && uint64_t(x)<(uint64_t(0)-n)%n) x=__uint128_t(rng())*n;
        return uint64_t(x>>64);
    }
    void advance(size_t i,double t) {
        auto &w=v[i]; double dt=t-w.time; w.time=t;
        auto k=w.leaves;
        if (dt<=0 || !k) return;
        if (k<8) {
            double q=-std::expm1(-dt);
            for (uint64_t j=0;j<k;++j) { double a=u(); if (a<q) w.leaves+=uint64_t(std::log(a)/std::log(q)); }
        } else {
            double intensity=std::gamma_distribution<double>(double(k),std::expm1(dt))(rng);
            w.leaves+=std::poisson_distribution<uint64_t>(intensity)(rng);
        }
    }
    uint64_t observe(double t) {
        while (next<=t) {
            size_t i=idx(v.size()); advance(i,next);
            size_t j=idx(v[i].leaves+v[i].core.size());
            if (j<v[i].core.size()) { auto w=v[i].core[j]; advance(w,next); ++v[w].leaves; }
            else { --v[i].leaves; v[i].core.push_back(v.size()); v.push_back({1,next,{i}}); }
            next+=e()/v.size();
        }
        uint64_t n=v.size();
        for (size_t i=0;i<v.size();++i) { advance(i,t); n+=v[i].leaves; }
        return n;
    }
    void add(size_t i,int delta) { for (++i;i<bit.size();i+=i&-i) bit[i]+=delta; }
    size_t find(int64_t k) {
        size_t i=0,s=1; while (s<bit.size()) s*=2;
        for (;s;s/=2) if (i+s<bit.size() && bit[i+s]<=k) { k-=bit[i+s]; i+=s; }
        return i;
    }
    bool grow(uint64_t n,uint64_t target) {
        double t=0;
        for (auto &w:v) w.time=0;
        next=e()/v.size();
        while (n+400<target) {
            double r=target-n, margin=12*std::sqrt(r)+50;
            t+=std::min(1.0,std::log1p((r-margin)/n));
            n=observe(t);
            if (n>target) return false;
        }
        bit.assign(v.size()+target-n+1,0);
        for(size_t i=0;i<v.size();++i) bit[i+1]=v[i].leaves;
        for(size_t i=1;i<bit.size();++i) if (i+(i&-i)<bit.size()) bit[i+(i&-i)]+=bit[i];
        while(n<target) {
            size_t i=idx(n);
            if (i>=v.size()) { auto w=find(i-v.size()); ++v[w].leaves; add(w,1); }
            else {
                auto j=idx(v[i].leaves+v[i].core.size());
                if (j<v[i].core.size()) { auto w=v[i].core[j]; ++v[w].leaves; add(w,1); }
                else { --v[i].leaves; add(i,-1); add(v.size(),1); v[i].core.push_back(v.size()); v.push_back({1,t,{i}}); }
            }
            ++n;
        }
        return true;
    }
};
int main(int argc,char **argv) {
    if(argc!=5) return 2;
    Sim s; s.rng.seed(std::strtoull(argv[1],nullptr,10));
    auto reps=std::strtoull(argv[2],nullptr,10), lo=std::strtoull(argv[3],nullptr,10), hi=std::strtoull(argv[4],nullptr,10);
    if(lo<3 || hi<lo) return 2;
    std::cout<<"c_lo,c_hi\n";
    uint64_t rejected=0;
    for(uint64_t j=0;j<reps;++j) {
        s.v.clear(); s.v.push_back({2,0,{}});
        if(!s.grow(3,lo)) { ++rejected; continue; }
        auto a=s.v.size();
        if(!s.grow(lo,hi)) { ++rejected; continue; }
        std::cout<<a<<','<<s.v.size()<<'\n';
    }
    std::cerr<<"attempted "<<reps<<", rejected "<<rejected<<'\n';
}
