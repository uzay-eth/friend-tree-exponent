// Original signed quartet gradients, including the 1/D endpoint interpolation.
#pragma once
#include "flow.hpp"
#include <map>

namespace next5 {
struct Deletion {
    int outside,site,old_shape;
    std::array<int,5> roles{}; // surviving members, time order
    mh::Family family=mh::P;std::array<int,4> canon{};int canon_site=-1;
};
inline std::vector<Deletion> deletions(const Parent& par,int n){
    std::vector<Deletion> out;
    for(int z=0;z<n;++z){
        int neighbors=(z!=0),touch=z?par[z]:-1;
        for(int j=z+1;j<n;++j)if(par[j]==z){++neighbors;touch=j;}
        if(neighbors!=1)continue;
        Deletion part{};part.outside=z;Parent index{},pr{};int left=0;
        for(int j=0;j<n;++j)if(j!=z){index[j]=left;part.roles[left++]=j;}
        part.site=index[touch];pr[0]=-1;
        for(int j=1;j<left;++j){pr[j]=index[par[part.roles[j]]];part.old_shape=part.old_shape*j+pr[j];}
        if(n==5){
            bool found=false;
            for(int f=mh::P;f<=mh::T&&!found;++f){
                std::array<int,4> perm{0,1,2,3};
                do{bool match=true;
                    for(int j=1;j<4;++j)match&=pr[perm[j]]==perm[mh::shapes[f].parent[j]];
                    if(match){
                        part.family=mh::Family(f);
                        for(int j=0;j<4;++j){part.canon[j]=part.roles[perm[j]];if(part.canon[j]==touch)part.canon_site=j;}
                        found=true;break;
                    }
                }while(std::next_permutation(perm.begin()+1,perm.end()));
            }
            need(found&&part.canon_site>=0,"quartet deletion lost its rooted shape");
        }
        out.push_back(part);
    }
    return out;
}
struct Mixture {int a,b,na,den;};
struct Fraction {W num;I den;};
struct Quartets {
    const mh::Snapshot& tab;const Cert4vLevel& level;
    bool multiple=false;
    std::map<std::array<int,3>,std::vector<Mixture>> pools;
    explicit Quartets(const mh::Snapshot& s):tab(s),level(s.c.L4){}
    Mixture exact(Cell x)const{
        int d=x.d,k=x.c-1;
        if(d<level.D0){int v=c4_lvl_inner(&level,d-x.c,k);return {v,v,1,1};}
        int v=c4_lvl_tail(&level,k);
        if(d==level.D0)return {v,v,1,1};
        return {v,level.INF_ID,level.D0,d};
    }
    const std::vector<Mixture>& pool(Cell x){
        auto [it,new_key]=pools.try_emplace({x.d,x.c,int(x.tail)});
        auto& options=it->second;
        if(!new_key)return options;
        for(int d=x.d;d<level.D0;++d)
            for(int c=x.c;c<=(x.tail?d:std::min(x.c,d));++c)options.push_back(exact({d,c}));
        int first=std::min(x.c-1,level.KT+1),last=x.tail?level.KT+1:first;
        for(int k=first;k<=last;++k){
            int d=std::max({level.D0,x.d,k+1});
            options.push_back(exact({d,k+1}));
        }
        options.push_back({level.INF_ID,level.INF_ID,1,1});return options;
    }
    Fraction weight(mh::Family f,const std::array<Mixture,4>& v)const{
        mh::Labels ids{};std::array<int,3> flex{};int n=0;I den=1;
        for(int j=0;j<4;++j){ids[j]=v[j].a;if(v[j].den>1){need(n<3,"too many quartet tail coordinates");flex[n++]=j;den*=v[j].den;}}
        int active=tab.active(f,ids);
        if(active<0)return {0,1};
        if(!n)return {tab.value(f,active,ids),1};
        W num=0;
        for(int mask=0;mask<(1<<n);++mask){I m=1;
            for(int k=0;k<n;++k){int j=flex[k];bool inf=mask&(1<<k);
                ids[j]=inf?v[j].b:v[j].a;m*=inf?v[j].den-v[j].na:v[j].na;}
            num+=W(m)*tab.value(f,active,ids);
        }
        return {num,den};
    }
    static bool endpoint(mh::Family f,int j){
        if(f==mh::P)return j==0||j==3;
        if(f==mh::F)return j==1||j==3;
        return f==mh::Y?j!=1:j!=0;
    }
    I slope(mh::Family f,int site,const std::array<Cell,4>& v){
        if(v[site].high)return 0;
        // A quartet with a large interior member is inactive, and stays so on gain.
        for(int j=0;j<4;++j)if(!endpoint(f,j)&&(v[j].high||v[j].d>=level.D0))return 0;
        std::array<Mixture,4> a{},b{};
        std::array<int,3> unseen{};int num=0;
        for(int j=0;j<4;++j){
            if(v[j].high){if(num&&!multiple)return 0;unseen[num++]=j;}
            else a[j]=exact(v[j]);
        }
        b=a;b[site]=exact(v[site].gain());
        if(!endpoint(f,site)&&v[site].d+1>=level.D0){int t=c4_lvl_tail(&level,v[site].c-1);b[site]={t,t,1,1};}
        auto change=[&](){auto p=weight(f,a),q=weight(f,b);
            return I((q.num*p.den-p.num*q.den)/(W(q.den)*p.den));};
        if(!num)return change();
        I low=INT64_MAX,high=INT64_MIN;
        std::array<const std::vector<Mixture>*,3> options{};size_t product=1;
        for(int j=0;j<num;++j){options[j]=&pool(v[unseen[j]]);product*=options[j]->size();}
        if(product>4096)return 0;
        auto explore=[&](auto&& self,int k)->void{
            if(k==num){I z=change();low=std::min(low,z);high=std::max(high,z);return;}
            for(auto mix:*options[k]){a[unseen[k]]=b[unseen[k]]=mix;self(self,k+1);
                if(low<=0&&high>=0)break;}
        };
        explore(explore,0);
        return low>0?low:(high<0?high:0);
    }
};
// One gradient table per leaf deletion. Its entries have four coordinates,
// while the donor supplies only its reciprocal interval in the five-set row.
struct Sources {
    struct Part {Deletion def;Parent stride{};std::vector<I> delta;};
    const Mesh& mesh;
    std::vector<Part> parts;
    template<class Delta>
    Sources(const Mesh& m,Delta calc):mesh(m){
        for(const auto& d:deletions(m.par,m.n)){
            if(m.donor_piece>=0&&m.donor_piece!=d.outside)continue;
            Part p;p.def=d;int length=1;
            for(int k=m.n-1;k-->0;){p.stride[k]=length;length*=m.width[d.roles[k]];}
            p.delta.resize(length);Five x{};
            for(int ix=0;ix<length;++ix){
                for(int k=0;k<m.n-1;++k){int j=d.roles[k];x[j]=m.axes[j].values[(ix/p.stride[k])%m.width[j]];}
                p.delta[ix]=calc(d,x);
            }
            parts.push_back(std::move(p));
        }
    }
    I operator()(const Parent& pos,const Five& x)const{
        I sum=0;
        for(const auto& p:parts){int code=0;
            for(int k=0;k<mesh.n-1;++k)code+=p.stride[k]*pos[p.def.roles[k]];
            sum+=donor(p.delta[code],x[p.def.outside]);
        }
        return sum;
    }
};
struct ActiveInterior {
    std::array<int,2> roles{};int width=1,n=0;
    std::vector<char> valid;
    ActiveInterior(const Mesh& m,const Deletion& d,const mh::Snapshot& witness){
        std::array<int,2> canonical{};
        for(int j=0;j<4;++j)if(!Quartets::endpoint(d.family,j)){
            canonical[n]=j;roles[n++]=d.canon[j];
        }
        if(n==2)width=m.width[roles[1]];
        valid.resize(m.width[roles[0]]*width);
        for(int i=0;i<int(valid.size());++i){mh::Labels ids{};
            for(int j=0;j<n;++j){
                auto x=m.axes[roles[j]].values[j?i%width:i/width];
                need(!x.high&&x.d<witness.c.L4.D0,"specialized internal position overflow");
                ids[canonical[j]]=c4_lvl_inner(&witness.c.L4,x.d-x.c,x.c-1);
            }
            valid[i]=witness.active(d.family,ids)>=0;
        }
    }
    bool operator()(const Parent& pos)const{return valid[pos[roles[0]]*width+(n==2?pos[roles[1]]:0)];}
};
} // namespace next5
