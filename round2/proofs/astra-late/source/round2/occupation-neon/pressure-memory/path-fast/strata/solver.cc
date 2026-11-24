// Borrow the frozen D18 table decoder, integer helpers and ABI types.
// No calls to the dense solver's allocator or transition code.
#include "../solver.cc"
#include <utility>

namespace {
using ID = std::uint32_t;
constexpr ID dead = std::numeric_limits<ID>::max();
struct Role {
    std::vector<Census> list;
    std::array<std::array<int,16>,17> digit{};
};
struct Layer {
    std::array<Role,4> role;
    std::array<ID,5> stride{1};
    ID offset=0;
};
struct Arc {
    ID grow=dead, stay=dead, extend=dead, attach=dead;
    int internal=0, unseen=0, leaves=0, d=0;
};
struct Pay {
    Value source{}, slope{};
};
struct Model {
    const Charges charges;
    std::array<int,4> cap, arity;
    int side, mother, mm, part;
    I unit=1000000, s0, ds;
    std::array<I,25> r{};
    std::array<std::array<Layer,2>,5> lay;
    ID size=0, tag_count=0;
    std::vector<ID> compact, face;
    std::vector<Pay> toll;
    std::vector<std::array<Value,4>> joint;
    std::array<std::vector<ID>,2> active;
    bool prune, merge;
    std::array<std::array<std::pair<std::size_t,std::size_t>,2>,5> segment{};
    std::uint64_t candidates=0, payable=0, calls=0, arrows=0, varying=0, checksum=0, budget;
    std::uint64_t residuals=0, slab_values=0;
    I range_lo=limit,range_hi=-limit;
    std::vector<Value> previous, overflow, absent;

    Model(const int* c,const int* a,int sc,int mc,const Tables* t,bool old,int p,int flags,std::uint64_t max)
        :charges(t,!old),side(sc),mother(mc),mm(std::max(1,mc)),part(p),prune(flags&1),merge(flags&2),budget(max) {
        check(p>=-3 && p<=3 && (side==0 || (side>=2 && side<=23)) &&
              (mother==0 || (mother>=2 && mother<=23)),"invalid strata parameters");
        int top=std::max({2,side+1,mother+1});
        for(int i=0;i<4;++i) {
            cap[i]=c[i]; arity[i]=a[i];
            check(c[i]>=2 && c[i]<=16 && a[i]>=int(i<3) && a[i]<c[i],"invalid path caps");
            top=std::max(top,c[i]);
        }
        for(int d=2;d<=top;++d) unit=std::lcm(unit,I(d));
        check(unit<=(I(1)<<48),"strata rate unit out of range");
        for(int d=2;d<=top;++d) r[d]=unit/d;
        s0=(unit/1000000)*565673; ds=(unit/1000000)*102;
        for(int n=1;n<=4;++n) for(int tag=0;tag<=int(p>=0);++tag) {
            auto& b=lay[n][tag]; b.offset=size;
            if(tag && p>=n) { b.stride[n]=0; continue; }
            for(int i=0;i<n;++i) {
                auto& z=b.role[i];
                for(auto& row:z.digit) row.fill(-1);
                int low=int(i+1<n)+int(i==p && tag);
                int high=arity[i]-int(n<4 && i==n-1)-int(i==p && !tag);
                for(int d=2;d<=cap[i];++d) for(int k=low;k<d && k<=high;++k) {
                    z.digit[d][k]=int(z.list.size()); z.list.push_back({d,k});
                }
                auto next=std::uint64_t(b.stride[i])*z.list.size();
                check(next<=40000000,"strata skeleton exceeds budget");
                b.stride[i+1]=ID(next);
            }
            check(std::uint64_t(size)+b.stride[n]<=40000000,"strata skeleton exceeds budget");
            size+=b.stride[n];
        }
        compact.assign(size,dead);
        face.assign(size,dead);
    }

    ID id(const State& s,bool tag) const {
        if(tag && part>=s.n) return dead;
        const auto& b=lay[s.n][tag];
        ID j=b.offset;
        for(int i=0;i<s.n;++i) {
            auto [d,k]=s.path[i];
            if(d>cap[i] || k>arity[i]) return dead;
            int t=b.role[i].digit[d][k];
            if(t<0) return dead;
            j+=b.stride[i]*ID(t);
        }
        return j;
    }

    State decode(ID j,int n,bool tag) const {
        const auto& b=lay[n][tag];
        j-=b.offset;
        State s; s.n=n;
        for(int i=0;i<n;++i) {
            const auto& z=b.role[i].list;
            s.path[i]=z[(j/b.stride[i])%z.size()];
        }
        return s;
    }

    std::array<Arc,4> topology(ID j,const State& s,bool tag) const {
        std::array<Arc,4> out;
        for(int i=0;i<s.n;++i) {
            auto [d,k]=s.path[i];
            auto& a=out[i];
            State z=s;
            ++z.path[i].d;
            a.grow=id(z,tag); a.d=d; a.leaves=d-k-1;
            a.unseen=k-int(i+1<s.n)-int(tag && i==part);
            check(a.unseen>=0,"missing physical child");
            if(a.leaves) {
                z=s; ++z.path[i].k;
                a.stay=id(z,tag);
                if(i==s.n-1 && s.n<4) { State q=z; q.path[q.n++]=virgin; a.extend=id(q,tag); }
                if(i==part && !tag) a.attach=id(z,true);
            }
            for(ID next:{a.grow,a.stay,a.extend,a.attach})
                check(next==dead || next>j,"strata skeleton is not topological");
        }
        return out;
    }

    void support() {
        // Backward topological sweep; live but unpayable rows map to the
        // same exact-zero continuation. Rate exists for every listed arc.
        for(int n=4;n>=1;--n) for(int t=int(part>=0);t>=0;--t) {
            const auto& b=lay[n][t];
            auto start=active[t].size();
            for(ID q=b.stride[n];q-->0;) {
                ID j=b.offset+q;
                auto s=decode(j,n,t);
                bool yes=false;
                if(n==4 && (part<0 || t)) {
                    ++candidates;
                    Reward z;
                    if(charges(s.path,z)==1) {
                        Pay a;
                        if(part<=-2) a.source=z.minimum;
                        else { a.source=z.source[part+1]; a.slope=z.slope[part+1]; }
                        if(part==-3) {
                            // Full affine row. Every pressure multiplies its
                            // own slope and drives the same receiver gain.
                            a.slope=z.slope[0];
                            for(int h=0;h<4;++h) {
                                W source=W(a.source[h])+z.source[0][h];
                                for(int i=0;i<4;++i) {
                                    int copies=s.path[i].k-int(i<3);
                                    source+=W(copies)*z.source[i+1][h];
                                }
                                a.source[h]=narrow(source);
                            }
                        }
                        bool nonzero=a.source!=Value{} || a.slope!=Value{};
                        if(part==-3) for(int i=1;i<5;++i) nonzero|=z.slope[i]!=Value{};
                        if(nonzero) {
                            // Source/slope are the exact outward ticks for both
                            // witnesses: keep the union of their nonzero faces.
                            yes=true; ++payable; face[j]=ID(toll.size()); toll.push_back(a);
                            if(part==-3) joint.push_back({z.slope[1],z.slope[2],z.slope[3],z.slope[4]});
                        }
                    }
                }
                if(prune && !yes) {
                    auto adj=topology(j,s,t);
                    for(int i=0;i<n && !yes;++i) for(ID w:{adj[i].grow,adj[i].stay,adj[i].extend,adj[i].attach})
                        if(w!=dead && compact[w]!=dead) { yes=true; break; }
                }
                if(yes || !prune) { compact[j]=ID(active[t].size()); active[t].push_back(j); }
            }
            segment[n][t]={start,active[t].size()};
        }
    }

    Value find(const std::vector<Value>& v,ID j,int m) const {
        return j==dead || compact[j]==dead ? Value{} : v[std::size_t(compact[j])*mm+m];
    }

    Value bellman(const Value& source,const std::array<Edge,20>& e,int ne) {
        I base=s0;
        arrows+=ne;
        for(int i=0;i<ne;++i) { base+=e[i].lo; varying+=e[i].lo!=e[i].hi; }
        Value out{};
        for(int j=0;j<4;++j) {
            bool upper=j%2;
            I den=base;
            W num=W(source[j])*unit;
            using Choice=std::pair<I,I>;
            std::array<Choice,21> order;
            int count=1; order[0]={0,ds};
            for(int i=0;i<ne;++i) {
                auto a=e[i]; num+=W(a.lo)*a.value[j];
                if(a.hi>a.lo) order[count++]={a.value[j],a.hi-a.lo};
            }
            std::sort(order.begin(),order.begin()+count,[upper](Choice a,Choice b){return upper?a.first>b.first:a.first<b.first;});
            for(int i=0;i<count;++i) {
                auto [v,inc]=order[i];
                W cmp=W(v)*den;
                if(upper?cmp<=num:cmp>=num) break;
                num+=W(v)*inc; den+=inc;
            }
            I v=out[j]=ratio(num,den,upper);
            W drift=W(source[j])*unit-W(s0+(((v<0)==upper)?ds:0))*v;
            for(int i=0;i<ne;++i) {
                auto a=e[i]; I delta=a.value[j]-v;
                drift+=W((delta>0)==upper?a.hi:a.lo)*delta;
            }
            check(upper?drift<=0:drift>=0,"strata drift inequality failed");
            ++residuals;
            range_lo=std::min(range_lo,v); range_hi=std::max(range_hi,v);
        }
        check(out[0]<=out[1] && out[2]<=out[3],"crossed strata enclosure");
        return out;
    }

    void slab(bool tag,int yd,int yk,std::vector<Value>& here,const std::vector<Value>& grown_donor,
              const std::vector<Value>& born_donor,const std::vector<Value>& bud_donor) {
        for(int n=4;n>=1;--n) {
            auto [start,end]=segment[n][tag];
            for(auto irow=start;irow<end;++irow) {
                ID j=active[tag][irow];
                auto s=decode(j,n,tag);
                auto adj=topology(j,s,tag);
                auto pay=face[j]==dead?Pay{}:toll[face[j]];
                for(int m=mm-1;m>=0;--m) {
                    if(++calls>budget) throw std::runtime_error("strata row budget exceeded");
                    int md=m==mm-1?0:m+2;
                    std::array<Edge,20> e;
                    int ne=0;
                    auto add=[&](I lo,I hi,Value z,const Value& w=Value{}) {
                        check(lo>=0 && lo<=hi && hi<=32*unit,"strata arrow range");
                        if(!hi) return;
                        check(ne<int(e.size()),"strata arrow limit");
                        for(int h=0;h<4;++h) z[h]=narrow(W(z[h])+w[h]);
                        e[ne++]={lo,hi,z};
                    };
                    for(int i=0;i<n;++i) {
                        auto a=adj[i];
                        auto g=find(here,a.grow,m);
                        I lo=a.leaves*unit+(i?r[s.path[i-1].d]:0)+(i+1<n?r[s.path[i+1].d]:0), hi=lo;
                        auto rate=[&](I l,I h,const Value& slope) {
                            if(merge && slope==Value{}) {lo+=l; hi+=h;}
                            else add(l,h,g,slope);
                        };
                        if(!i) rate(r[md],md?r[md]:r[std::max(2,mother+1)],(part==-1 || part==-3)?pay.slope:Value{});
                        rate(0,a.unseen*r[2],part==-3 && face[j]!=dead?joint[face[j]][i]:Value{});
                        if(tag && part==i) rate(r[yd],yd?r[yd]:r[std::max(2,side+1)],pay.slope);
                        add(lo,hi,g);
                        if(a.leaves) {
                            auto z=find(here,a.stay,m);
                            auto w=find(here,a.extend,m);
                            auto y=find(born_donor,a.attach,m);
                            for(int h=0;h<4;++h) z[h]=narrow(W(z[h])+w[h]+y[h]);
                            add(a.leaves*r[a.d],a.leaves*r[a.d],z);
                        }
                    }
                    if(yd>0) {
                        int leaves=yd-yk-1;
                        I fixed=leaves*unit+r[s.path[part].d];
                        add(fixed,fixed+yk*r[2],find(grown_donor,j,m));
                        if(leaves) add(leaves*r[yd],leaves*r[yd],find(bud_donor,j,m));
                    }
                    if(md) {
                        I child=r[s.path[0].d];
                        add(child,child+md*unit-3*r[2],find(here,j,m+1));
                    }
                    auto v=bellman(pay.source,e,ne);
                    here[std::size_t(compact[j])*mm+m]=v;
                    auto mix=[&](std::uint64_t x){checksum=(checksum^x)*1099511628211ULL;};
                    mix(j); mix(yd); mix(yk); mix(m); for(auto x:v) mix(std::uint64_t(x));
                }
            }
        }
    }

    void run(I* entries,std::uint64_t* stats,I* span) {
        support();
        auto nt=std::size_t(mm)*active[1].size(), na=std::size_t(mm)*active[0].size();
        // Generation axis: at degree d a gain points to d+1 at same
        // arity; a promotion points within d at arity k+1. Never transpose
        // them: previous[k] must still be intact when computing current[k].
        const std::vector<Value> none;
        checksum=14695981039346656037ULL;
        if(part>=0 && nt) {
            overflow.resize(nt);
            slab(true,0,0,overflow,none,none,none);
            if(side) {
                using Degree=std::vector<std::vector<Value>>;
                Degree layer(side);
                for(int d=side;d>=2;--d) {
                    for(int k=d-1;k>=0;--k) {
                        layer[k].resize(nt);
                        // The d+1 successor has exactly this path and mother
                        // index, so its old slot is read before we overwrite it.
                        slab(true,d,k,layer[k],d==side?overflow:layer[k],none,k+1<d?layer[k+1]:none);
                    }
                    if(d==side) std::vector<Value>().swap(overflow);
                    std::vector<Value>().swap(layer[d-1]);
                    // At d=2 only k=0 is used by registration.
                    if(d==2) previous=std::move(layer[0]);
                }
                slab_values=nt*(1+std::size_t(side)); // exact peak: d=side plus overflow
            } else previous=std::move(overflow);
        }
        absent.resize(na);
        slab(false,0,-1,absent,none,previous,none);
        slab_values=std::max(slab_values,na+previous.size());
        State s; s.path[0]=virgin;
        for(int m=0;m<mm;++m) {
            int lane=m?m-1:mm-1;
            auto z=find(absent,id(s,false),lane);
            std::copy(z.begin(),z.end(),entries+4*m);
        }
        check(residuals==4*calls,"strata residual coverage mismatch");
        std::array<std::uint64_t,11> out{size,candidates,payable,active[0].size(),active[1].size(),
                                         calls,arrows,varying,residuals,slab_values,checksum};
        std::copy(out.begin(),out.end(),stats);
        span[0]=calls?range_lo:0; span[1]=calls?range_hi:0;
    }
};
} // namespace

extern "C" int strata(const int* c,const int* a,int sc,int mc,const Tables* t,int old,int p,int flags,
                      std::uint64_t budget,I* entries,std::uint64_t* stats,I* range) {
    try { Model model(c,a,sc,mc,t,old,p,flags,budget); model.run(entries,stats,range); return 0; }
    catch(const std::exception& e) { std::fprintf(stderr,"strata: %s\n",e.what()); return 1; }
}
