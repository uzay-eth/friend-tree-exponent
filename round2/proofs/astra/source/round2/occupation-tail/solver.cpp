// Lower discounted occupation using interpolation in each observed reciprocal degree.
#include "../motif-cert/reserve/reserve.hpp"
#include "../motif-high/occupation/credit.hpp"
#include <array>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {
using I = std::int64_t;
using Wide = __int128_t;
constexpr I unit = I(1)<<40;
constexpr int max_arity = 4;
using Indices = std::array<int,max_arity>;
void check(bool v,const char *message) {if (!v) throw std::runtime_error(message);}
I narrow(Wide v) {
    check(v>=INT64_MIN && v<=INT64_MAX,"integer out of range");
    return static_cast<I>(v);
}
Wide floor_ratio(Wide n,Wide d) {
    check(d>0,"nonpositive denominator");
    return n/d-(n%d<0);
}
I down(Wide n,Wide d) { return narrow(floor_ratio(n,d)); }

enum class Kind { inner, tail, infinity };
struct State { int d, c; Kind kind; };
struct Mixture { int first=-1, second=-1, weight=0, denom=1; };

struct Axis {
    int threshold, cmax, first_c, tbase, ibase;
    std::vector<State> labels;
    std::vector<std::vector<int>> inner;
    Axis(int degree,int cores,int children):threshold(degree),cmax(cores),first_c(children+1),
        inner(degree,std::vector<int>(cores+1,-1)) {
        check(2<=threshold && 1<=cmax && cmax<=threshold && threshold<=2048,"axis cutoff");
        for (int d=2;d<threshold;++d) for(int c=first_c;c<=std::min(d,cmax);++c) {
            inner[d][c]=int(labels.size()); labels.push_back({d,c,Kind::inner});
        }
        tbase=int(labels.size());
        for(int c=first_c;c<=cmax;++c) labels.push_back({threshold,c,Kind::tail});
        ibase=int(labels.size());
        for(int c=first_c;c<=cmax;++c) labels.push_back({threshold,c,Kind::infinity});
    }
    Mixture project(State s) const {
        if(s.c>cmax) return {};
        check(s.c>=first_c && s.d>=2,"infeasible successor");
        int ti=tbase+s.c-first_c, inf=ibase+s.c-first_c;
        if(s.kind==Kind::infinity) return {inf,-1,1,1};
        if(s.d<threshold) {
            check(s.kind==Kind::inner && s.c<=s.d,"fine/coarse ordering");
            return {inner[s.d][s.c],-1,1,1};
        }
        if(s.d==threshold) return {ti,-1,1,1};
        return {ti,inf,threshold,s.d};
    }
    int old_id(State s) const {
        if (s.c>cmax) return -1;
        if(s.kind==Kind::infinity) return ibase+s.c-first_c;
        if(s.d==threshold) return tbase+s.c-first_c;
        return inner[s.d][s.c];
    }
};

using Members = std::array<State,max_arity>;
struct Grid {
    int n;
    Indices parent, child_count{}, stride{};
    std::vector<Axis> axis;
    std::vector<I> values;
    Grid(int r,Indices p,const Indices& d,const Indices& c):n(r),parent(p) {
        for(int i=1;i<n;++i) ++child_count[p[i]];
        for(int i=0;i<n;++i) axis.emplace_back(d[i],c[i],child_count[i]);
        size_t count=1;
        for(int i=n-1;i>=0;--i) {
            stride[i]=int(count);
            count*=axis[i].labels.size();
            check(count<100000000,"grid exceeds the single-table memory limit");
        }
        values.resize(count);
    }
    I interpolate(const Members& x) const {
        // Exact tensor interpolation followed by one downward rounding.
        std::array<int,16> offsets{};
        std::array<Wide,16> mass{};
        mass[0]=1;
        int count=1;
        Wide total=1;
        for(int i=0;i<n;++i) {
            Mixture m=axis[i].project(x[i]);
            if(m.first<0)return 0;
            if(m.second<0) {
                for(int j=0;j<count;++j) offsets[j]+=stride[i]*m.first;
            } else {
                total*=m.denom;
                for(int j=0;j<count;++j) {
                    offsets[j+count]=offsets[j]+stride[i]*m.second;
                    mass[j+count]=mass[j]*(m.denom-m.weight);
                    offsets[j]+=stride[i]*m.first;
                    mass[j]*=m.weight;
                }
                count*=2;
            }
        }
        Wide sum=0;
        for(int j=0;j<count;++j) sum+=mass[j]*values[offsets[j]];
        return down(sum,total);
    }
};

struct Range {
    I low, high, denom;
    I destination;
    Wide term(I z,I precision) const {
        I diff=destination-z;
        return floor_ratio(Wide(precision)*diff*(diff>=0?low:high),denom);
    }
};
struct Slot { int donor, parent, child; motif_reserve::Slot position; };
struct Pruned {
    int outside, contact;
    std::array<int,3> retained;
    mh_occ::Shape triple;
    mh_occ::Direction direction;
};

std::vector<Pruned> quartet_donors(const Indices& parent) {
    std::vector<Pruned> cuts;
    for(int outside=0;outside<4;++outside) {
        int contact=parent[outside],degree=(outside!=0);
        for(int j=1;j<4;++j)if(parent[j]==outside){++degree;contact=j;}
        if(degree!=1)continue;
        Indices renumber{},p{-1,-1,-1,-1};
        std::array<int,3> retained{};
        int used=0;
        for(int i=0;i<4;++i)if(i!=outside){renumber[i]=used;retained[used++]=i;}
        for(int i=1;i<3;++i)p[i]=renumber[parent[retained[i]]];
        auto dir=parent[contact]==outside?mh_occ::Direction::Parent:mh_occ::Direction::Child;
        cuts.push_back({outside,renumber[contact],retained,mh_occ::Shape::from_parents(3,p.data()),dir});
    }
    return cuts;
}

class Solver {
    motif_reserve::Scores low_;
    mh_occ::Source high_;
    std::array<Indices,max_arity> degrees_,cores_;
    int mask_, rmax_;
    I num_,den_;
    std::ofstream dump_;
    std::array<uint64_t,5> count_{},active_{};
    static constexpr I precision=4096;
    static bool inner(State x) {return x.kind==Kind::inner;}
    static motif_reserve::Core loc(State x) {return {x.d,x.c};}
    static std::vector<Slot> directed_slots(const Indices& p) {
        using motif_reserve::Slot;
        if(p[2]==1) return {{0,1,2,Slot::ParentOfParent},{2,0,1,Slot::ChildOfChild}};
        return {{1,0,2,Slot::OtherChildOfParent},{2,0,1,Slot::OtherChildOfParent}};
    }
    I credit(const Grid& g,const Members& x,const mh_occ::Shape& sh,
             const std::vector<Slot>& slots,const std::vector<Pruned>& cuts) {
        if(g.n==1) return (mask_&1) && inner(x[0]) ? low_.singleton_floor(loc(x[0])):0;
        if(g.n==2) {
            if(!(mask_&2) || !inner(x[0]) || !inner(x[1]) || x[0].d>72 || x[1].d>72)return 0;
            return low_.edge_unused_floor(loc(x[0]),loc(x[1]));
        }
        I f=0;
        if(g.n==3 && (mask_&4)) for(auto s:slots) {
            State a=x[s.parent], b=x[s.child], z=x[s.donor];
            if(inner(a) && inner(b) && a.d<=72 && b.d<=72)
                f+=low_.edge_donor_floor(loc(a),loc(b),s.position,{z.d,!inner(z)});
        }
        mh_occ::Members m{};
        for(int i=0;i<g.n;++i)m[i]={x[i].d,inner(x[i])?x[i].c:0};
        if(mask_&8) f+=high_.row(sh,m.data());
        if(g.n==4 && (mask_&16)) {
            if(!(mask_&32))f+=high_.triple_donors_on_four(sh,m.data());
            else for(const auto& cut:cuts) {
                mh_occ::Members triple{};
                for(int j=0;j<3;++j)triple[j]=m[cut.retained[j]];
                State donor=x[cut.outside];
                bool far=donor.kind==Kind::infinity;
                // J(delta,p) is affine in p. The near corner uses its value
                // at 1/d and the far corner a lower bound valid at p=0.
                f+=high_.triple_donor(cut.triple,triple.data(),cut.contact,cut.direction,
                                      far?INT_MAX-1:donor.d,far);
            }
        }
        check(f>=0 && f<=8*unit,"source magnitude");
        return f;
    }

    Grid solve(int n,Indices parent) {
        Grid g(n,parent,degrees_[n-1],cores_[n-1]);
        if(g.values.empty())return g;
        mh_occ::Shape sh{};
        if(n>=3)sh=mh_occ::Shape::from_parents(n,parent.data());
        std::vector<Slot> slots;
        if(n==3)slots=directed_slots(parent);
        std::vector<Pruned> cuts;
        if(n==4)cuts=quartet_donors(parent);
        std::vector<Grid> extended;
        if(n<rmax_)for(int i=0;i<n;++i){
            parent.at(n)=i; extended.push_back(solve(n+1,parent));
        }
        for(int row=int(g.values.size())-1;row>=0;--row) {
            Members x{};
            Indices pos{};
            for(int i=0;i<n;++i) {
                pos[i]=(row/g.stride[i])%g.axis[i].labels.size();
                x[i]=g.axis[i].labels[pos[i]];
            }
            I f=credit(g,x,sh,slots,cuts);
            if(f)++active_[n];
            std::array<Range,2*max_arity> ranges{};
            int used=0;
            for(int i=0;i<n;++i) {
                State current=x[i];
                auto old=[&](State to) {
                    int j=g.axis[i].old_id(to);
                    check(j<0 || j>pos[i],"old state does not advance");
                    return j<0?I(0):g.values[row+(j-pos[i])*g.stride[i]];
                };
                if(inner(current)) {
                    I q=2;
                    auto selected=[&](int j){return (j>0 && parent[j]==i)||(i>0 && parent[i]==j);};
                    for(int j=0;j<n;++j)if(selected(j))q*=x[j].d;
                    I lo=(current.d-current.c)*q, hi=lo;
                    int outside=current.c-(i!=0)-g.child_count[i];
                    for(int j=0;j<n;++j) if(selected(j)) {
                        if(inner(x[j]))lo+=q/x[j].d;
                        hi+=q/x[j].d;
                    }
                    hi+=outside*(q/2);
                    ranges[used++]={lo,hi,q,old({current.d+1,current.c,Kind::inner})};
                } else if(current.kind==Kind::tail) {
                    ranges[used++]={current.d-current.c,current.d+1,current.d+1,
                                    old({current.d,current.c,Kind::infinity})};
                }
                if(inner(current) && current.c==current.d)continue;
                State bud{current.d,current.c+1,current.kind};
                I post=old(bud);
                if(n<rmax_) {
                    auto newset=x;newset[i]=bud;newset.at(n)={2,1,Kind::inner};
                    post+=extended[i].interpolate(newset);
                }
                I lo=current.d-current.c, hi=inner(current)?lo:current.d;
                ranges[used++]={lo,hi,current.d,post};
            }
            // Each nonnegative interval acts on its own signed successor jump.
            // The endpoints used for the proposal have no role in acceptance.
            double z=0;
            for(int it=0;it<=used+1;++it) {
                double a=double(f),b=double(num_)/den_;
                for(int i=0;i<used;++i) {
                    auto r=ranges[i];
                    double rate=double(r.destination>=z?r.low:r.high)/r.denom;
                    a+=rate*r.destination; b+=rate;
                }
                double next=a/b;
                if(std::abs(z-next)<.1){z=next;break;}
                z=next;
            }
            check(std::isfinite(z) && z>=0 && z<1e18,"invalid floating proposal");
            I candidate=std::min<I>(I(z),unit);
            auto residual=[&](I v) {
                Wide a=Wide(f)*precision+down(-Wide(precision)*num_*v,den_);
                for(int i=0;i<used;++i)a+=ranges[i].term(v,precision);
                return a;
            };
            Wide gap=residual(candidate);
            if(gap<0)candidate=std::max<I>(0,candidate+down(gap*den_,Wide(precision)*num_));
            check(residual(candidate)>=0,"outward interval Bellman row failed");
            g.values[row]=candidate;
            if(dump_.is_open()) {
                if(n==1)dump_<<"-";
                for(int i=1;i<n;++i)dump_<<parent[i];
                for(int i=0;i<n;++i)
                    dump_<<" "<<x[i].d<<" "<<x[i].c<<" "<<int(x[i].kind);
                dump_<<" "<<candidate<<" "<<f<<"\n";
            }
        }
        count_[n]+=g.values.size();
        return g;
    }
public:
    Solver(const char* file,int side,int mask,int rmax,const std::array<Indices,4>& d,
           const std::array<Indices,4>& c,I num,I den)
       :low_(file),high_(file,side==1?56581165:5656453,side==1?100000000:10000000),
        degrees_(d),cores_(c),mask_(mask),rmax_(rmax),num_(num),den_(den) {
        check((side==1||side==2) && (side==1?1:-1)==high_.sign(),"snapshot side");
        check(1<=rmax && rmax<=4 && mask>0 && mask<64 && num>0 && den>0 &&
              num<=1000000 && den<=1000000,"invalid options");
        for(int n=1;n<rmax;++n)for(int i=0;i<n;++i)
            check(d[n][i]<=d[n-1][i] && c[n][i]<=c[n-1][i],"nonnested thresholds");
        if(const char *path=std::getenv("FT2_TAIL_DUMP")) {
            dump_.open(path);check(dump_.good(),"cannot create fixture");
        }
    }
    void run() {
        Grid root=solve(1,{-1,0,0,0});
        Members newborn{};newborn[0]={2,1,Kind::inner};
        std::cout<<"mask "<<mask_<<" at "<<num_<<"/"<<den_<<" rows";
        for(int i=1;i<=rmax_;++i)std::cout<<" "<<count_[i];
        std::cout<<" positive sources";
        for(int i=1;i<=rmax_;++i)std::cout<<" "<<active_[i];
        std::cout<<" newborn "<<root.interpolate(newborn)<<"/"<<unit<<"\n";
    }
};

std::array<Indices,4> config(const std::string &text) {
    std::array<Indices,4> a{};
    std::istringstream in(text);
    for(int n=1;n<=4;++n)for(int i=0;i<n;++i) {
        check(bool(in>>a[n-1][i]),"expected 10 ordered layer cutoffs");
        if(n<4||i<n-1){char c;check(bool(in>>c) && c==',',"cutoff separator");}
    }
    std::string rest;check(!(in>>rest),"excess cutoffs");
    return a;
}
} // namespace

int main(int argc,char**argv) try {
    check(argc==7 || argc==9,"usage: solver FILE SIDE(1|2) MASK ARITY D1,D21,D22,... C1,C21,C22,... [discount numerator denominator]");
    Solver solver(argv[1],std::stoi(argv[2]),std::stoi(argv[3]),std::stoi(argv[4]),
                  config(argv[5]),config(argv[6]),argc==9?std::stoll(argv[7]):1133,
                  argc==9?std::stoll(argv[8]):2000);
    solver.run();
} catch(const std::exception& e) {std::cerr<<e.what()<<"\n";return 1;}
