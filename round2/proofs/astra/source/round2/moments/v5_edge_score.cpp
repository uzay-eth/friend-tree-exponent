// Unsigned comparison conditioned on the parent class of each bottom vertex.
#define main archived_edge_comparison_main
#include "../motif-cert/comparison.cpp"
#undef main
#include <climits>
#include <map>
#include <set>

static int group(int d,int k,int K) {
    if(d<=K) return d*(d-1)/2-1+k;
    constexpr int degrees[]={8,16,32,48,71},children[]={0,1,2,3,7,15,31,63};
    return 10000+100*int(std::lower_bound(std::begin(degrees),std::end(degrees),d)-std::begin(degrees))
           +int(std::lower_bound(std::begin(children),std::end(children),k)-std::begin(children));
}
using Groups=std::map<int,std::set<int>>;
static void record(Groups& groups,const Level& l,int g,int d,int k) {
    groups[g].insert(corner(l,d-k-1,k));
    if(d>l.D0) groups[g].insert(l.INF_ID);
}
using Table=std::vector<std::vector<i64>>;
static Table table(int rows,int cols,i64 init=0) {
    return Table(rows,std::vector<i64>(cols,init));
}
struct Conditional {
    Snapshot &s;
    Compare cmp;
    const Layout &c;
    const Level &v,&q;
    int K;
    Table vert,cherry,chain,fork,chain_high,fork_high,y,trident;
    Conditional(Snapshot &snapshot,int cutoff):s(snapshot),cmp(s),c(s.C),v(c.L3),q(c.L4),K(cutoff),
        vert(table(c.NHU_ACT,v.NG)),cherry(table(c.NCU_ACT,v.NG,INT64_MIN)),
        chain(table(q.NG,q.NG)),fork(chain),chain_high(chain),fork_high(chain),
        y(table(c.NY4_ACT,q.NG,INT64_MIN)),trident(table(c.NT4_ACT,q.NG,INT64_MIN)) {
        for(int h=0;h<c.NHU_ACT;h++) for(int a=0;a<v.npars;a++) for(int w=0;w<v.NG;w++)
            vert[h][w]=std::max(vert[h][w],s.v[(size_t(a)*c.NHU_ACT+h)*v.nchs+v.w_to_wi[w]]);
        for(int h=0;h<c.NCU_ACT;h++) for(int a=0;a<v.NG;a++) for(int b=a;b<v.NG;b++) {
            i64 z=s.c[size_t(h)*v.pair_n+c4_lvl_pair(&v,a,b)];
            cherry[h][a]=std::max(cherry[h][a],z);cherry[h][b]=std::max(cherry[h][b],z);
        }
        for(int h=0;h<c.NP4_ACT;h++) for(int a=0;a<q.npars;a++) for(int w=0;w<q.NG;w++) {
            int m=c.p4_u2[h]; i64 z=s.p4[c4_p4_off(&c,h,q.pars[a],w)];
            chain[m][w]=std::max(chain[m][w],z);
            if(q.gD[c.p4_u1[h]]>K) chain_high[m][w]=std::max(chain_high[m][w],z);
        }
        for(int h=0;h<c.NF4_ACT;h++) for(int a=0;a<q.NG;a++) for(int w=0;w<q.NG;w++) {
            int top=c.f4_u1[h],m=c.f4_u2[h];
            i64 z=(q.gK[top]-1)*s.f4[c4_f4_off(&c,h,a,w)];
            fork[m][w]=std::max(fork[m][w],z);
            if(q.gD[top]>K) fork_high[m][w]=std::max(fork_high[m][w],z);
        }
        for(int h=0;h<c.NY4_ACT;h++) for(int p=0;p<q.npars;p++)
            for(int a=0;a<q.NG;a++) for(int b=a;b<q.NG;b++) {
                i64 z=s.y4[c4_y4_off(&c,h,q.pars[p],a,b)];
                y[h][a]=std::max(y[h][a],z);y[h][b]=std::max(y[h][b],z);
            }
        for(int h=0;h<c.NT4_ACT;h++) for(int a=0;a<q.NG;a++)
            for(int b=a;b<q.NG;b++) for(int e=b;e<q.NG;e++) {
                i64 z=s.t4[c4_t4_off(&c,h,a,b,e)];
                for(int w:{a,b,e}) trident[h][w]=std::max(trident[h][w],z);
            }
    }
    void run() {
        Groups w3,w4,outer4;
        constexpr int last=73;
        check(q.D0<last && v.D0<last && q.KT+1<last && v.KT+1<last
              && 3<=K && K<q.D0 && K<=7,"unsupported score cutoff");
        for(int d=2;d<=last;d++) for(int k=0;k<d;k++) {
            int g=group(d,k,K);
            record(w3,v,g,d,k);record(w4,q,g,d,k);record(outer4,q,d<=K?g:-1,d,k);
        }
        std::map<std::pair<int,int>,i64> edge;
        for(int d=2;d<=last;d++) for(int k=1;k<d;k++) {
            int L=d-k-1,m=group(d,k,K),u=d<q.D0?c4_lvl_inner(&q,L,k):-1;
            if(d<=K) {
                for(const auto &[b,ws]:outer4) {
                    i64 price=INT64_MAX;
                    for(int w:ws)
                        price=std::min(price,cmp.chain[w].hi-chain_high[u][w]+cmp.fork[w].hi-fork_high[u][w]);
                    check(price>=0,"negative high-parent comparison");
                    if(price) std::cout<<"t "<<m<<' '<<b<<' '<<price<<'\n';
                }
                continue;
            }
            int hu=c4_hu_idx(&c,L,k),cu=c4_cu_idx(&c,L,k);
            int ya=c4_y4_act(&c,L,k),ta=c4_t4_act(&c,L,k);
            int h=hu<0?-1:c.hu_to_act[hu],ch=cu<0?-1:c.cu_to_act[cu];
            for(const auto &[b,ws]:w3) {
                i64 lo3=INT64_MAX,lo4=INT64_MAX;
                for(int w:ws) {
                    i64 z=cmp.vert[w].hi-(h<0?0:vert[h][w]);
                    if(ch>=0) z+=(k-1)*(cmp.cherry[ch].hi-cherry[ch][w])/2;
                    lo3=std::min(lo3,z);
                }
                for(int w:w4.at(b)) {
                    i64 z=cmp.chain[w].hi+cmp.fork[w].hi;
                    if(u>=0) z-=chain[u][w]+fork[u][w];
                    if(ya>=0) z+=(k-1)*(cmp.y[ya].hi-y[ya][w])/2;
                    if(ta>=0) z+=(k-1)*(k-2)*(cmp.t[ta].hi-trident[ta][w])/6;
                    lo4=std::min(lo4,z);
                }
                check(lo3>=0 && lo4>=0,"negative conditional edge defect");
                auto [it,fresh]=edge.try_emplace({m,b},lo3+lo4);
                if(!fresh) it->second=std::min(it->second,lo3+lo4);
            }
        }
        for(const auto &[key,z]:edge) if(z)
            std::cout<<"e "<<key.first<<' '<<key.second<<' '<<z<<'\n';
    }
};
int main(int argc,char**argv)try {
    check(argc==3,"usage: v5_edge_score SNAPSHOT JOINT_CUTOFF|--probe");
    Snapshot s(argv[1]);
    if(std::string(argv[2])=="--probe") {
        Compare cmp(s);
        int d,k;
        while(std::cin>>d>>k) {
            check(d>=2 && 0<=k && k<d,"invalid probe state");
            int L=d-k-1,cu=c4_cu_idx(&s.C,L,k),ya=c4_y4_act(&s.C,L,k),ta=c4_t4_act(&s.C,L,k);
            std::cout<<mixture(cmp.vert,s.C.L3,d,k).hi<<' '<<mixture(cmp.chain,s.C.L4,d,k).hi
                     <<' '<<mixture(cmp.fork,s.C.L4,d,k).hi<<' '<<(cu<0?0:cmp.cherry[s.C.cu_to_act[cu]].hi)
                     <<' '<<(ya<0?0:cmp.y[ya].hi)<<' '<<(ta<0?0:cmp.t[ta].hi)<<'\n';
        }
    } else Conditional(s,std::stoi(argv[2])).run();
    return 0;
} catch(const std::exception& e) {std::cerr<<e.what()<<'\n';return 2;}
