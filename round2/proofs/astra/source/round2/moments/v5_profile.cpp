// Use the nominal integer row evaluator unchanged. Aggregate only after evaluating each row.
#define main archived_low_main
#include "../motif-cert/low.cpp"
#undef main
#include <climits>
#include <fstream>
#include <map>
#include <set>

using Key = std::array<int, 3>;
struct Values {
    i64 f = INT64_MIN, score = INT64_MIN, small = INT64_MAX, big = INT64_MIN;
    int dmin = INT_MAX, dmax = 0, kmin = INT_MAX, kmax = 0;
};
static constexpr int degree_end[] = {8,16,32,48,71,95,127,191,255,383,511,767,1023,1399};
static constexpr int child_end[] = {0,1,2,3,7,15,31,63,INT_MAX};

static int band(int x, const auto &ends) {
    return int(std::lower_bound(std::begin(ends),std::end(ends),x)-std::begin(ends));
}
static int edge_group(int d,int k,int cut) {
    if(d<=cut) return d*(d-1)/2-1+k;
    return 10000+100*std::min(5,band(d,degree_end))+band(k,child_end);
}

int main(int argc,char **argv) try {
    check(argc==7,"usage: v5_profile SNAPSHOT MBAR EDGE_CUTOFF VERTEX_CUTOFF NUM DEN");
    Snapshot s(argv[1]); int K=std::stoi(argv[3]), M=std::stoi(argv[4]);
    check(K>=3 && K<=7 && K<M && M<s.H.DF-1 && s.H.D0_2==72 && s.H.DF==1400,"unsupported cutoffs");
    int num=std::stoi(argv[5]),den=std::stoi(argv[6]);
    check(0<num && num<=den && den<=1000000000,"invalid nominal exponent");
    Low low(s,num,den);
    const auto &l=s.C.L2;
    check(l.KT<64,"unrepresented exact child-count tail label");
    std::ifstream in(argv[2],std::ios::binary);
    check(in.good(),"mbar open");
    std::map<Key,Values> tails;
    const int sign=s.H.mode==1?1:-1;
    std::cout<<"header "<<s.H.mode<<" "<<num<<" "<<den<<" "<<s.H.S<<" "<<s.H.DF
             <<" "<<s.H.f_p0<<" "<<s.H.f_p1<<" "<<s.H.f_p2<<" "<<l.KT<<"\n";
    for(int d=2;d<s.H.DF;d++) for(int k=0;k<d;k++) {
        int L=d-k-1;
        i64 score=0; in.read(reinterpret_cast<char*>(&score),sizeof(score));
        check(in.good(),"short mbar profile");
        i64 f=low.vertex(L,k), delta=0;
        if(d>=l.D0 && d<s.H.DF-1)
            delta=sign*(s.f[c4_vidx(&s.C,L+1,k)]-s.f[c4_vidx(&s.C,L,k)]);
        if(d<=M) {
            std::cout<<"v "<<d<<" "<<k<<" "<<f<<" "<<score<<" "<<delta<<"\n";
        } else {
            Key key{band(d,degree_end),band(k,child_end),16*L/d};
            auto &v=tails[key];
            v.f=std::max(v.f,f);v.score=std::max(v.score,score);
            v.small=std::min(v.small,delta);v.big=std::max(v.big,delta);
            v.dmin=std::min(v.dmin,d);v.dmax=std::max(v.dmax,d);
            v.kmin=std::min(v.kmin,k);v.kmax=std::max(v.kmax,k);
        }
    }
    check(in.peek()==std::char_traits<char>::eof(),"long mbar profile");
    for(const auto &[key,v]:tails) {
        std::cout<<"b "<<key[0]<<" "<<key[1]<<" "<<key[2]<<" "<<v.dmin<<" "<<v.dmax
                 <<" "<<v.kmin<<" "<<v.kmax<<" "<<v.f<<" "<<v.score<<" "<<v.small<<" "<<v.big<<"\n";
    }
    std::vector<i64> edge_max(l.NG,INT64_MIN);
    for(int a=0;a<l.npars;a++) for(int w=0;w<l.NG;w++)
        edge_max[w]=std::max(edge_max[w],s.G(l.pars[a],w));
    auto groups=[&](int i,bool parent) {
        std::set<int> out;
        if(l.gKind[i]==0) out.insert(edge_group(l.gD[i],l.gK[i],K));
        else if(l.gKind[i]==1 || l.gK[i]==l.KT+1) {
            int kmin=l.gKind[i]==1?int(parent):l.KT+1;
            for(int k=kmin;k<=std::max(64,kmin);k++) out.insert(edge_group(l.D0,k,K));
        } else out.insert(edge_group(l.D0,l.gK[i],K));
        return std::vector<int>(out.begin(),out.end());
    };
    std::map<std::pair<int,int>,std::pair<i64,i64>> edges;
    std::vector<std::vector<int>> pg(l.NG),wg(l.NG);
    for(int i=0;i<l.NG;i++) {pg[i]=groups(i,true);wg[i]=groups(i,false);}
    for(int a=0;a<l.npars;a++) for(int w=0;w<l.NG;w++) {
        int p=l.pars[a];i64 drift=low.edge(p,w), defect=edge_max[w]-s.G(p,w);
        for(int gi:pg[p]) for(int gj:wg[w]) {
            auto [it,inserted]=edges.try_emplace({gi,gj},drift,defect);
            if(!inserted) {
                it->second.first=std::max(it->second.first,drift);
                it->second.second=std::min(it->second.second,defect);
            }
        }
    }
    for(const auto &[pair,v]:edges)
        std::cout<<"e "<<pair.first<<" "<<pair.second<<" "<<v.first<<" "<<v.second<<"\n";
    return 0;
} catch(const std::exception &e) {std::cerr<<e.what()<<"\n";return 2;}
