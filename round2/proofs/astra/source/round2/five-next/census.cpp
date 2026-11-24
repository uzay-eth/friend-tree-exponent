// Independent discrete physical cohorts. Both original quartet tables use
// their full multilinear endpoint interpolation on each observed tree.
#include "../numerics4/physical.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <vector>

using namespace numeric4;
struct CoreTree {
    std::vector<int> parent;
    std::vector<Count> degree;
    std::vector<std::vector<int>> children;
};
struct Totals {
    uint64_t sets=0,incidences=0;
    Sum all,small,tail,old_tail,seed_donor,small_all,up,down;
};
Pair coupon(Pair delta,Count degree){
    auto j=[&](Real a){return std::max(a,Real(0))/2-a/degree;};
    return {j(-delta.lo),j(delta.hi)};
}
void census(const Tables& tables,const CoreTree& t,uint64_t seed,Count size){
    std::array<Totals,NF> stats;
    auto quartet=[&](Family f,std::array<int,4> ids){
        Nodes a{};
        Count max_degree=0;bool baseline=true;
        for(int j=0;j<4;++j){
            int x=ids[j],c=1+int(t.children[x].size());Count d=t.degree[x];
            a[j].s={Count(d-c),U(c-1)};max_degree=std::max(max_degree,d);
            baseline&=d<=6&&c<=5;
        }
        if(tables.active(f,a)<0)return;
        auto& s=stats[f];++s.sets;
        Pair old=tables.value(f,a);
        for(int j=0;j<4;++j){
            const int x=ids[j];auto after=a;after[j].s=after[j].s.gain();
            const Pair delta=tables.value(f,after)-old;
            auto donor=[&](int y,bool parent){
                if(std::find(ids.begin(),ids.end(),y)!=ids.end())return;
                Pair q=coupon(delta,t.degree[y]);
                if(y==0){s.seed_donor.add(q);return;}
                ++s.incidences;s.all.add(q);(parent?s.up:s.down).add(q);
                if(baseline)s.small.add(q);
                if(baseline&&t.degree[y]<=6&&t.children[y].size()<5)s.small_all.add(q);
                if(max_degree>=Count(tables.c.L4.D0))s.tail.add(q);
                if(t.degree[ids[0]]>=Count(tables.c.L4.D0))s.old_tail.add(q);
            };
            donor(t.parent[x],true);
            for(int y:t.children[x])donor(y,false);
        }
    };
    for(int u=1;u<int(t.parent.size());++u){
        const auto& ch=t.children[u];
        for(int mid:ch)for(int below:t.children[mid]){
            for(int end:t.children[below])quartet(Path,{u,mid,below,end});
            for(int sibling:ch)if(sibling!=mid)quartet(Fork,{u,sibling,mid,below});
        }
        for(size_t i=0;i<ch.size();++i)for(size_t j=i+1;j<ch.size();++j){
            if(t.parent[u]>0)quartet(Stem,{t.parent[u],u,ch[i],ch[j]});
            for(size_t k=j+1;k<ch.size();++k)quartet(Trident,{u,ch[i],ch[j],ch[k]});
        }
    }
    for(int f=Path;f<NF;++f){
        const auto& s=stats[f];
        for(int side=0;side<2;++side){
            std::cout<<(side?"upper":"lower")<<" "<<seed<<" "<<size<<" "<<t.parent.size()<<" "<<names[f]
                     <<" "<<s.sets<<" "<<s.incidences;
            for(const Sum* a:{&s.all,&s.small,&s.small_all,&s.tail,&s.old_tail,&s.up,&s.down,&s.seed_donor})
                std::cout<<" "<<(side?a->value.hi:a->value.lo);
            std::cout<<"\n";
        }
    }
    std::cout<<std::flush;
}
void grow(const Tables& tables,int n,uint64_t seed){
    std::mt19937_64 rng(seed);
    auto pick=[&](uint64_t bound){uint64_t r,threshold=(-bound)%bound;do r=rng();while(r<threshold);return int(r%bound);};
    std::vector<int> parent(n),cid(n,-1);cid[0]=0;
    std::vector<std::vector<int>> incident{{1,2}};
    CoreTree tree{{-1},{2},{{}}};
    for(int v=3;v<n;++v){
        int w=pick(v),i=cid[w];
        int u=i<0?parent[w]:incident[i][pick(incident[i].size())];
        parent[v]=u;i=cid[u];
        if(i<0){
            int a=parent[u],p=cid[a],j=int(tree.parent.size());check(p>=0,"missing parent core");
            cid[u]=j;incident.push_back({a,v});tree.children[p].push_back(j);
            tree.parent.push_back(p);tree.degree.push_back(2);tree.children.emplace_back();
        }else{incident[i].push_back(v);++tree.degree[i];}
        if(v+1==n||v+1==n/8)census(tables,tree,seed,v+1);
    }
}
void read_tree(const Tables& tables,const char* filename){
    std::ifstream input(filename);CoreTree t;int p;Count leaves,all_leaves=0;
    while(input>>p>>leaves){
        int id=int(t.parent.size());check((id==0&&p==-1)||(id>0&&p>=0&&p<id),"invalid core parent");
        t.parent.push_back(p);t.degree.push_back(leaves+(id!=0));t.children.emplace_back();
        if(id){t.children[p].push_back(id);++t.degree[p];}
        all_leaves+=leaves;
    }
    check(input.eof()&&!t.parent.empty(),"invalid shared tree");
    census(tables,t,0,all_leaves+t.parent.size());
}
int main(int argc,char** argv)try{
    check(argc==5,"usage: census LOWER UPPER vertex-count replicates | --tree core-file");Tables t(argv[1],argv[2]);
    std::cout<<std::setprecision(20)<<"side seed n cores family active incidences J small small_all tail oldest_tail parent child founder\n";
    if(std::string(argv[3])=="--tree"){read_tree(t,argv[4]);return 0;}
    int n=std::stoi(argv[3]),count=std::stoi(argv[4]);check(n>=1024&&n<=1<<25&&count>0,"invalid census size");
    for(int k=0;k<count;++k)grow(t,n,UINT64_C(0x8c4872a8f51b3d07)+uint64_t(k)*UINT64_C(0xd1342543de82ef95));
    return 0;
}catch(const std::exception& e){std::cerr<<e.what()<<"\n";return 1;}
