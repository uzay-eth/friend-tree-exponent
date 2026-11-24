#include "quartet.hpp"
#include "projection.hpp"
#include <fstream>
#include <memory>

namespace next5 {
using Profile=std::array<std::vector<Cuts>,6>;
Profile profile(std::string mode,int pass){
    Profile p;
    if(mode=="small"||mode=="small-split"){
        p[1]={{{13,6,3}}};p[2]={{{8,4,3},{6,4,3}}};
        p[3]={{{7,4,3},{5,3,3},{4,3,3}}};p[4]={Cuts(4,{4,3,3})};
        p[5]={Cuts(5,{4,3,3}),{{6,4,3},{4,3,3},{4,3,3},{3,3,3},{3,3,3}}};
    }else if(mode=="pilot"||pass){
        p[1]={{{1600,33,7}}};p[2]={{{140,12,6},{70,10,6}}};
        p[3]={{{37,8,5},{22,6,5},{14,5,4}}};
        p[4]={{{15,6,5},{10,5,4},{8,4,4},{7,4,3}}};
        p[5]={Cuts(5,{6,5,4}),{{14,5,5},{8,4,4},{6,4,3},{5,3,3},{5,3,3}}};
        if(pass){
            p[5]={{{9,4,4},{6,4,3},{5,3,3},{5,3,3},{5,3,3}}};
            if(pass>=2){p[4]={{{10,4,4},{7,4,3},{6,3,3},{5,3,3}}};p[5]={Cuts(5,{5,3,3})};}
        }
    }else if(mode=="wide"||mode=="split"){
        p[1]={{{14000,76,10}}};p[2]={{{600,21,9},{200,14,7}}};
        p[3]={{{120,11,7},{40,7,6},{24,6,5}}};
        p[4]={{{32,7,6},{15,5,5},{10,5,4},{8,4,4}}};
        p[5]={Cuts(5,{7,5,4}),{{24,6,6},{10,5,4},{7,4,4},{6,4,3},{5,3,3}}};
    }else throw std::runtime_error("unknown mesh");
    return p;
}
struct Hierarchy {
    std::array<std::vector<std::vector<Mesh>>,6> levels;
    I newborn=0;bool split;
    static Cuts adapted(const Deletion& def,bool small){
        Cuts window(5);
        int donor=def.outside;
        window[donor]=small?Cut{4,3,3}:(donor==0?Cut{26,6,8}:Cut{7,4,5});
        for(int j=0;j<4;++j){
            int pos=def.canon[j];bool outside=Quartets::endpoint(def.family,j);
            if(!outside){
                int limit=def.family==mh::T?7:9;
                window[pos]=small?Cut{4,3,3,true}:Cut{limit,5,3,true};
            }else if(small)window[pos]={pos==0?6:4,3,4};
            else if(pos==0)window[pos]={38,7,10};
            else if(pos==1)window[pos]={10,5,7};
            else window[pos]={pos==2?8:6,4,6};
        }
        return window;
    }
    template<class Get>
    I combine(int n,int sh,Get get)const{
        I full=0,by_donor=0;
        for(const auto& m:levels[n][sh]){
            I z=get(m);
            if(m.donor_piece<0)full=std::max(full,z);else by_donor+=z;
        }
        return std::max(full,by_donor);
    }
    I value(int n,int sh,const Five& x)const{
        return combine(n,sh,[&](const Mesh& m){return m.value(x);});
    }
    I entry(int n,int sh,const Five& x)const{
        return combine(n,sh,[&](const Mesh& m){return split?m.bounded_entry(x):m.entry(x);});
    }
    I correction(int n,const Deletion& part,const Five& x)const{
        Five old{};
        for(int k=0;k<n;++k)old[k]=x[part.roles[k]];
        if(old[part.site].high)return 0;
        for(const auto& candidate:levels[n][part.old_shape])for(int j=0;j<n;++j){
            auto a=old[j];const auto& axis=candidate.axes[j];
            if(a.high&&(a.d<=axis.cut.d||(a.tail&&a.c<axis.cut.high_c)))return 0;
        }
        I prev=value(n,part.old_shape,old);old[part.site]=old[part.site].gain();
        return prev-value(n,part.old_shape,old);
    }
    using Project=std::function<void(const Mesh&)>;
    Hierarchy(Quartets& q,const Profile& cuts,const Hierarchy* prev,bool retain_five,bool adaptive,bool small,
              const Project& project={},bool only_five=false):split(adaptive){
        const auto tree=shapes();
        for(int n=5;n>=1;--n){
            uint64_t rows=0;
            for(int j=0;j<int(tree[n].size());++j){
                levels[n].emplace_back();auto& group=levels[n].back();
                std::vector<std::pair<int,Cuts>> local;
                for(const auto& windows:cuts[n])local.emplace_back(-1,windows);
                if(split&&!prev&&n==5)for(const auto& d:deletions(tree[n][j],n))
                    local.emplace_back(d.outside,adapted(d,small));
                for(const auto& [donor,windows]:local){
                    auto& m=group.emplace_back(tree[n][j],j,windows);m.donor_piece=donor;
                    Mesh::Entry birth;
                    if(n<5)birth=[&,n,j](int site,const Five& x){return entry(n+1,j*n+site,x);};
                    if((prev&&n>1)||(!prev&&n==5)){
                        Sources source(m,[&](const Deletion& part,const Five& x){
                            if(prev)return prev->correction(n-1,part,x);
                            std::array<Cell,4> a{};for(int k=0;k<4;++k)a[k]=x[part.canon[k]];
                            return q.slope(part.family,part.canon_site,a);
                        });
                        if(donor>=0){
                            need(source.parts.size()==1,"specialized donor decomposition");
                            ActiveInterior live(m,source.parts[0].def,q.tab);
                            m.solve(std::cref(source),birth,retain_five||bool(project),std::cref(live));
                        }else m.solve(std::cref(source),birth,n<5||retain_five||bool(project));
                    }else m.solve({},birth,true);
                    if(n==5&&project){project(m);if(!retain_five)std::vector<I>().swap(m.weights);}
                    rows+=m.count;
                }
            }
            std::cout<<" rank "<<n<<" rows "<<rows<<"\n";
            if(only_five)break;
        }
        if(!only_five){Five nb{Cell{2,1}};newborn=value(1,0,nb);}
    }
};
} // namespace next5

#ifndef FIVE_NEXT_NO_MAIN
int main(int argc,char** argv)try{
    using namespace next5;
    need(argc==4||argc==6,"usage: reserve SNAPSHOT small|pilot|wide PASSES [requests output]");
    std::cout.setf(std::ios::unitbuf);mh::Snapshot snap(argv[1]);
    need(snap.h.S==S,"incorrect fixed-point scale");
    std::string mode=argv[2];int passes=std::stoi(argv[3]);need(passes>=1&&passes<=5,"unsupported number of passes");
    const bool small=mode=="small"||mode=="small-split";
    const bool split=mode=="split"||mode=="small-split";
    const bool requests=argc==6;need(!requests||small,"query mode retains full small score tables");
    struct Query{int n,sh;Five x;I score=0;};std::vector<Query> probe;
    if(requests){std::ifstream input(argv[4]);int n,j;
        while(input>>n>>j){Query r{n,j,{}};for(int k=0;k<n;++k)need(bool(input>>r.x[k].d>>r.x[k].c),"bad query");
            probe.push_back(r);}need(input.eof()&&!probe.empty(),"empty queries");}
    Quartets quartets(snap);quartets.multiple=split;std::unique_ptr<Hierarchy> prior;I total=0;
    const char* side=snap.sign==1?"upper":"lower";
    for(int pass=0;pass<passes;++pass){
        std::cout<<side<<" "<<mode<<" pass "<<pass<<"\n";
        std::vector<Compact> projection;Hierarchy::Project capture;
        if(mode=="split"&&pass==0){
            projection=projection_tables(false);
            capture=[&](const Mesh& m){projection[m.id].receive(m);};
        }
        auto next=std::make_unique<Hierarchy>(quartets,profile(mode,pass),prior.get(),requests,split,small,capture);
        if(capture){
            for(auto& t:projection)t.seal();
            write_projection(std::string("round2/five-next/build/primary-")+side+".bin",snap.sign,next->newborn,projection);
        }
        total+=next->newborn;
        for(auto& r:probe)r.score+=next->value(r.n,r.sh,r.x);
        std::cout<<side<<" quartet pass "<<pass<<" newborn "<<next->newborn<<"/"<<S<<"\n";
        prior=std::move(next);
    }
    std::cout<<side<<" quartet replacement "<<total<<"/"<<S<<"\n";
    if(requests){std::ofstream out(argv[5]);out<<total<<"\n";for(const auto& r:probe)out<<r.score<<"\n";need(out.good(),"query output failed");}
    return 0;
}catch(const std::exception&e){std::cerr<<e.what()<<"\n";return 1;}
#endif
