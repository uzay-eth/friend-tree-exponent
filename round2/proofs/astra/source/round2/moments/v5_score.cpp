// Comparison defects for known three-site overlaps of the unsigned potential.
#define main archived_comparison_main
#include "../motif-cert/comparison.cpp"
#undef main
#include <map>

using Classes=std::map<int,std::vector<int>>;
static int id(int d,int k){return d*(d-1)/2-1+k;}
static Classes classes(const Level& l,int K,bool parent){
    Classes z;
    for(int i=0;i<l.NG;i++){
        if(parent && l.p_to_pi[i]<0)continue;
        int group=l.gKind[i]==0 && l.gD[i]<=K?id(l.gD[i],l.gK[i]):-1;
        z[group].push_back(i);
    }
    return z;
}
template <class Fn> static i64 least(const std::vector<int>& xs,const std::vector<int>& ys,Fn price){
    i64 val=INT64_MAX;
    for(int x:xs)for(int y:ys)val=std::min(val,price(x,y));
    check(val>=0 && val<INT64_MAX,"invalid potential defect");return val;
}
int main(int argc,char**argv)try{
    check(argc==3,"usage: v5_score SNAPSHOT JOINT_CUTOFF");
    Snapshot s(argv[1]);Compare cmp(s);int K=std::stoi(argv[2]);check(K>=3 && K<=7,"joint cutoff");
    auto &c=s.C;auto &l=c.L4;
    auto p3=classes(c.L3,K,true),w3=classes(c.L3,K,false);
    auto p4=classes(c.L4,K,true),w4=classes(c.L4,K,false);
    // Best chain for each internal pair and bottom with its top still free.
    std::vector<i64> best(size_t(c.NP4_ACT)*l.NG,INT64_MIN);
    for(int z=0;z<c.NP4_ACT;z++)for(int a=0;a<l.npars;a++)for(int w=0;w<l.NG;w++){
        auto& value=best[size_t(z)*l.NG+w];
        value=std::max(value,s.p4[c4_p4_off(&c,z,l.pars[a],w)]);
    }
    for(int d=2;d<=K;d++)for(int k=1;k<d;k++){
        int L=d-1-k,m=id(d,k),hu=c4_hu_idx(&c,L,k);
        for(const auto &[a,aa]:p3)for(const auto &[b,bb]:w3){
            i64 price=least(aa,bb,[&](int x,int w){return cmp.vert[w].hi-s.V(x,hu,w);});
            if(a>=0){
                int ix=aa[0],ak=c.L3.gK[ix],aL=c.L3.gL[ix];
                int pa=c4_p4_act(&c,aL,ak,L,k),fa=c4_f4_act(&c,aL,ak,L,k);
                i64 low=INT64_MAX;
                for(int w:w4[b]){
                    i64 fmax=0;
                    if(fa>=0){
                        fmax=INT64_MIN;
                        for(int j=0;j<l.NG;j++)fmax=std::max(fmax,s.f4[c4_f4_off(&c,fa,j,w)]);
                        fmax*=ak-1;
                    }
                    i64 pc=cmp.chain[w].hi-(pa<0?0:best[size_t(pa)*l.NG+w]);
                    i64 fc=cmp.fork[w].hi-fmax;
                    check(pc>=0 && fc>=0,"bottom comparison defect");
                    low=std::min(low,pc+fc);
                }
                price+=low;
            }
            if(price)std::cout<<"t "<<a<<" "<<m<<" "<<b<<" "<<price<<"\n";
        }
        int pa=c4_p4_act(&c,L,k,0,1);
        if(pa>=0)for(const auto &[a,aa]:p4)for(const auto &[b,bb]:w4){
            i64 price=least(aa,bb,[&](int x,int w){return best[size_t(pa)*l.NG+w]-s.p4[c4_p4_off(&c,pa,x,w)];});
            if(price)std::cout<<"q "<<a<<" "<<m<<" "<<b<<" "<<price<<"\n";
        }
        int cu=c4_cu_idx(&c,L,k);
        if(cu<0)continue;
        int ya=c4_y4_act(&c,L,k),ta=c4_t4_act(&c,L,k);
        for(const auto &[a,aa]:w3)for(const auto &[b,bb]:w3){
            if(a>b)continue;
            i64 price=least(aa,bb,[&](int x,int w){return cmp.cherry[c.cu_to_act[cu]].hi-s.Ch(cu,x,w);});
            price+=least(w4[a],w4[b],[&](int x,int w){
                i64 ygap=0,tgap=0;
                if(ya>=0){
                    i64 hi=INT64_MIN;
                    for(int j=0;j<l.npars;j++)hi=std::max(hi,s.y4[c4_y4_off(&c,ya,l.pars[j],x,w)]);
                    ygap=cmp.y[ya].hi-hi;
                }
                if(ta>=0){
                    i64 hi=INT64_MIN;
                    for(int j=0;j<l.NG;j++)hi=std::max(hi,s.t4[c4_t4_off(&c,ta,x,w,j)]);
                    tgap=(k-2)*(cmp.t[ta].hi-hi)/3;
                }
                return ygap+tgap;
            });
            if(price)std::cout<<"c "<<m<<" "<<a<<" "<<b<<" "<<price<<"\n";
        }
    }
    return 0;
}catch(const std::exception &e){std::cerr<<e.what()<<"\n";return 2;}
