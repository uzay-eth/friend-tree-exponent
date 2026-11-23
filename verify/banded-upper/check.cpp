// Reconstruct a positive local potential and verify its upper drift by integer arithmetic.
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>
using std::vector;
using I = std::int64_t;
using W = __int128_t;
using Bag = vector<int>;
using Edge = std::array<int,4>;

void require(bool value, const char* message) {
    if (!value) throw std::runtime_error(message);
}
I narrow(W v) {
    require(v >= std::numeric_limits<I>::min() && v <= std::numeric_limits<I>::max(), "integer out of range");
    return I(v);
}
I up(W n, I d) {
    require(d > 0, "nonpositive denominator");
    return narrow(n/d + (n%d > 0));
}
struct Recipe {
    I scale=0;
    int threshold=0;
    vector<I> rewards;
    std::map<Edge,I> gains;
    std::map<std::pair<int,int>,I> arrows;
    explicit Recipe(const std::string& path) {
        std::ifstream input(path);
        require(input.good(),"cannot read recipe");
        std::string line,kind;
        while(std::getline(input,line)) {
            std::istringstream row(line);
            if(!(row>>kind) || kind[0]=='#') continue;
            int D,c,d,k,i; I w;
            if(kind=="scale") row>>scale;
            else if(kind=="threshold") row>>threshold;
            else if(kind=="reward") {
                row>>i>>w; require(i==int(rewards.size()),"reward order"); rewards.push_back(w);
            } else if(kind=="edge") {
                row>>D>>c>>d>>k>>w; require(gains.emplace(Edge{D,c,d,k},w).second,"duplicate edge");
            } else if(kind=="arrow") {
                row>>D>>c>>w; require(c>0 && c<=D && w>=0 && arrows.emplace(std::make_pair(D,c),w).second,"bad arrow");
            } else throw std::runtime_error("unknown recipe record");
            require(!row.fail(),"invalid recipe field");
            std::string extra; require(!(row>>extra),"unexpected recipe field");
        }
        require(scale>0 && scale<=100000000000000LL && threshold>=2 && !rewards.empty(),"bad recipe header");
        for(I v:rewards) require(v>0 && v<=16*scale,"reward out of range");
        for(auto it:gains) require(it.second>=-scale && it.second<=scale,"gain transfer out of range");
        for(auto it:arrows) require(it.second<=scale,"arrow out of range");
    }
    I arrow(int D,int c) const {
        auto pos=arrows.find({D,c}); return pos==arrows.end()?0:pos->second;
    }
    I gain(int D,int c,int d,int k) const {
        auto pos=gains.find({D,c,d,k}); return pos==gains.end()?0:pos->second;
    }
};
struct Type {int d,k,e=0,j=0;};
enum Event { GAIN, CONVERT, FAR_GAIN, FAR_CONVERT };
struct Change {int destination,label,count; Event event;};
struct Group {int label,count;};
struct State {
    vector<Group> groups;
    vector<Change> changes;
    // Non-newborn successors have degree D or D+1. Keep just those two slices.
    std::array<I,2> weights{0,0};
    std::array<int,2> degrees{-1,-1};
    int promoted=-1, coarser=-1;
    I outgoing_arrows=0;
};
struct Layer {
    int c,cap,kcap,first_pool,start,stop,global_limit;
    vector<Type> types;
    vector<vector<int>> offsets;
    vector<State> states;
    std::map<Bag,int> indices;
    Layer(int c_,int h,int low,int high,int limit):c(c_),cap(h),kcap(c==1?20:c==2?7:5),
        start(std::max(low,c)),stop(high),global_limit(limit) {
        offsets.resize(cap);
        for(int d=2;d<cap;++d) {
            offsets[d].resize(std::min(d,kcap)+1);
            for(int k=1;k<=std::min(d,kcap);++k) {
                offsets[d][k]=int(types.size());
                if(witness(d,k)) {
                    for(int e=2;e<6;++e) for(int j=1;j<=e;++j) types.push_back({d,k,e,j});
                    types.push_back({d,k,6,0});
                } else types.push_back({d,k});
            }
        }
        first_pool=int(types.size()); types.push_back({cap,0});
        for(int d=4;d<limit;d*=2) if(d>cap) types.push_back({d,0});
        if(limit>cap) types.push_back({limit,0});
        Bag bag; enumerate(bag,0);
    }
    bool witness(int d,int k) const {return c==1 && d<cap && d<=14 && k>=2 && k<=4;}
    int label(int d,int k,int e=0,int j=0) const {
        if(d>=cap) {
            int t=first_pool;
            while(t+1<int(types.size()) && types[t+1].d<=d) ++t;
            return t;
        }
        require(k>=1 && k<=d,"invalid known neighbour");
        int pos=offsets[d][std::min(k,kcap)];
        if(!witness(d,k)) return pos;
        require(e>=2,"missing witness");
        if(e>=6) return pos+14; // 2+3+4+5 exact witness types
        require(j>=1 && j<=e,"invalid witness");
        return pos+(e-2)*(e+1)/2+j-1;
    }
    void enumerate(Bag& bag,int lo) {
        if(int(bag.size())==c) {
            int id=int(states.size()); indices.emplace(bag,id);
            State s;
            for(int t:bag) {
                if(s.groups.empty() || s.groups.back().label!=t) s.groups.push_back({t,1});
                else ++s.groups.back().count;
            }
            states.push_back(std::move(s)); return;
        }
        for(int t=lo;t<int(types.size());++t) {
            if(t>first_pool && int(bag.size())+1<c) break;
            bag.push_back(t); enumerate(bag,t); bag.pop_back();
        }
    }
    int index(Bag bag) const {
        std::sort(bag.begin(),bag.end());
        for(size_t j=0;j+1<bag.size();++j) if(bag[j]>first_pool) bag[j]=first_pool;
        return indices.at(bag);
    }
    I get(int id,int D) const {
        const auto& row=states.at(id);
        require(D>=start && D<stop && row.degrees[D%2]==D,"uncomputed successor");
        return row.weights[D%2];
    }
};
struct Band {
    int low,high;
    vector<Layer> layers;
    Band(int lo,int hi,int H,const vector<int>& caps):low(lo),high(hi) {
        layers.reserve(hi);
        for(int c=0;c<hi;++c) {
            int hc=(c>=1 && c<=int(caps.size()))?caps[c-1]:2;
            require(hc>=2 && hc<=H && (c<=1 || hc<=layers.back().cap),"invalid cutoffs");
            layers.emplace_back(c,hc,lo,hi,H);
        }
    }
};
struct Term {I target,lo,hi,den;};
struct Checker {
    I p,q;
    Recipe recipe;
    vector<std::unique_ptr<Band>> bands;
    vector<int> ends;
    int H=0;
    I A=0,B=0,amax=0,minimum=std::numeric_limits<I>::max(),maximum=0,slack_min=std::numeric_limits<I>::max();
    std::uint64_t rows=0;
    Checker(const std::string& geometry,const std::string& input,I numerator,I denominator):
        p(numerator),q(denominator),recipe(input) {
        require(p>0 && p<q && q<=1000000,"invalid test exponent");
        std::ifstream f(geometry); require(f.good(),"cannot read band geometry");
        vector<vector<int>> caps;
        std::string line;
        while(std::getline(f,line)) {
            std::istringstream row(line); int end;
            if(!(row>>end)) continue;
            require(end>=4 && end<=8000 && (ends.empty() || end>ends.back()),"invalid band endpoints");
            vector<int> h; int v;
            while(row>>v) h.push_back(v);
            require(h.size()>=3 && h[2]>=6,"witness needs exact degrees below six");
            if(!caps.empty()) require(h[0]==caps[0][0],"inconsistent newborn alphabet");
            ends.push_back(end); caps.push_back(h);
        }
        require(!ends.empty(),"empty geometry"); H=ends.back();
        int low=2;
        for(size_t i=0;i<ends.size();++i) {
            bands.push_back(std::make_unique<Band>(low,ends[i],H,caps[i])); low=ends[i];
        }
        auto& first=bands.front()->layers[1];
        require(recipe.rewards.size()==first.states.size(),"wrong number of rewards");
        require(recipe.threshold<=H,"bad arrow threshold");
        for(auto [dc,a]:recipe.arrows) {
            require(dc.first>=2 && dc.first<H,"arrow at tail or non-core"); amax=std::max(amax,a);
        }
        for(auto [key,g]:recipe.gains) {
            auto [D,c,d,k]=key;
            require(2<=D && D<H && 2<=d && d<H && 1<=c && c<=D && 1<=k && k<=d,"invalid edge key");
            auto& u=band(D).layers[c]; auto& w=band(d).layers[k];
            require(d<u.cap && k<u.kcap && D<w.cap && c<w.kcap,"unobservable edge key");
        }
        I s=recipe.rewards.at(first.label(H,0));
        B=std::max<I>(0,narrow(W(s-amax)*q/(q+p)));
        A=std::max(up(W(B+amax)*q,p),up(W(s)*q,p)-B);
        require(A>0 && B>=0 && W(B+amax)*q<=W(p)*A && W(s)*q<=W(p)*(A+B),"tail inequalities failed");
        minimum=A; maximum=A+B;
    }
    Band& band(int degree) {return *bands.at(std::upper_bound(ends.begin(),ends.end(),degree)-ends.begin());}
    static Bag expand(const State& s) {
        Bag b; for(auto g:s.groups) b.insert(b.end(),g.count,g.label); return b;
    }
    I receiver_max(const Layer& l,int i) const {
        Type w=l.types[i]; I best=0;
        for(auto [dc,a]:recipe.arrows) {
            auto [d,k]=dc;
            const bool possible=i>=l.first_pool ? d>=w.d :
                (d==w.d && (w.k==l.kcap ? k>=w.k : k==w.k));
            if(possible) best=std::max(best,a);
        }
        return best;
    }
    void transitions(size_t bi) {
        auto& b=*bands[bi];
        for(auto& l:b.layers) {
            Layer* next=bi+1<bands.size()? &bands[bi+1]->layers[l.c] : nullptr;
            require(!next || next->cap<=l.cap,"increasing cutoff between bands");
            vector<I> charges;
            for(int i=0;i<int(l.types.size());++i) charges.push_back(receiver_max(l,i));
            for(auto& s:l.states) {
                Bag bag=expand(s), fresh{ }, projected;
                auto* after=l.c+1<b.high ? &b.layers[l.c+1] : nullptr;
                if(after) fresh.push_back(after->label(2,1));
                for(int t:bag) {
                    Type w=l.types[t];
                    if(next) projected.push_back(next->label(w.d,w.k,w.e,w.j));
                    if(after) fresh.push_back(after->label(w.d,w.k));
                }
                if(next) s.coarser=next->index(projected);
                if(after) s.promoted=after->index(fresh);
                size_t slot=0;
                for(auto group:s.groups) {
                    int t=group.label,n=group.count;
                    s.outgoing_arrows=narrow(W(s.outgoing_arrows)+W(n)*charges[t]);
                    if(t<l.first_pool) {
                        auto w=l.types[t];
                        auto add=[&](int target,Event kind) {
                            int before=bag[slot]; bag[slot]=target;
                            s.changes.push_back({l.index(bag),t,n,kind}); bag[slot]=before;
                        };
                        add(l.label(w.d+1,w.k,w.e,w.j),GAIN);
                        if(w.k<w.d && w.k<l.kcap) add(l.label(w.d,w.k+1,2,1),CONVERT);
                        if(w.j) {
                            add(l.label(w.d,w.k,w.e+1,w.j),FAR_GAIN);
                            if(w.j<w.e) add(l.label(w.d,w.k,w.e,w.j+1),FAR_CONVERT);
                        }
                    }
                    slot+=n;
                }
            }
        }
    }
    int newborn(int D,const Layer& l,const State& s) const {
        auto &first=bands.front()->layers[1];
        if(first.witness(D,l.c+1)) {
            require(!s.groups.empty(),"missing old neighbour at birth");
            Type w=l.types[s.groups.front().label];
            return first.label(D,l.c+1,w.d,w.k);
        }
        return first.label(D,l.c+1);
    }
    void solve(size_t bi) {
        auto& b=*bands[bi]; vector<Term> terms;
        for(int D=b.high-1;D>=b.low;--D) for(int c=D;c>=0;--c) {
            auto& l=b.layers[c]; int L=D-c;
            I arrow=recipe.arrow(D,c);
            vector<I> gain_in(l.types.size()),gain_out(l.types.size());
            for(int t=0;t<l.first_pool;++t) {
                auto w=l.types[t];
                if(w.k<l.kcap) {
                    gain_in[t]=recipe.gain(D,c,w.d,w.k);
                    gain_out[t]=recipe.gain(w.d,w.k,D,c);
                }
            }
            for(int i=int(l.states.size())-1;i>=0;--i) {
                auto& s=l.states[i];
                I target=D+1<b.high?l.get(i,D+1):
                    bi+1<bands.size()?bands[bi+1]->layers[c].get(s.coarser,D+1):A+up(W(B)*(L+1),H);
                for(auto g:s.groups) target=narrow(W(target)+W(g.count)*gain_out[g.label]);
                terms.clear();
                if(L) terms.push_back({target,L,L,1});
                for(auto g:s.groups) {
                    int d=l.types[g.label].d;
                    I adjusted=target-(d>=recipe.threshold?arrow:0);
                    terms.push_back({adjusted,g.label<l.first_pool?g.count:0,g.count,d});
                }
                if(L) {
                    I at_birth=recipe.rewards.at(newborn(D,l,s));
                    terms.push_back({narrow(W(b.layers[c+1].get(s.promoted,D))+at_birth),L,L,D});
                }
                for(auto mv:s.changes) {
                    auto w=l.types[mv.label];
                    I x=l.get(mv.destination,D);
                    I lo=0,hi=0,den=1;
                    if(mv.event==FAR_GAIN) {
                        den=2*w.d; lo=den*(w.e-w.j)+2; hi=lo+w.d*(w.j-1);
                    } else if(mv.event==FAR_CONVERT) {
                        den=w.e; hi=lo=w.e-w.j;
                    } else if(mv.event==CONVERT) {
                        den=w.d; hi=lo=w.d-w.k;
                    } else {
                        x-=gain_in[mv.label];
                        if(w.e) {
                            den=2LL*D*w.e;
                            I base=den*(w.d-w.k)+2*w.e;
                            lo=base+(w.j?2*D:0); hi=base+D*(w.k-2)*w.e+2*D;
                        } else {
                            den=2*D; lo=2+(w.k<l.kcap?den*(w.d-w.k):0);
                            hi=2+D*(2*w.d-w.k-1);
                        }
                    }
                    terms.push_back({x,lo*mv.count,hi*mv.count,den});
                }
                I fee=D>=recipe.threshold?up(s.outgoing_arrows,D):0;
                // The floating proposal only selects the first integer to try.
                double v=target;
                for(int iter=0;iter<8;++iter) {
                    double num=double(fee),den=double(p)/q;
                    for(auto t:terms) {
                        double rate=double(t.target>v?t.hi:t.lo)/t.den;
                        num+=rate*t.target; den+=rate;
                    }
                    double candidate=std::max(1.,num/den);
                    if(std::abs(v-candidate)<0.25) {v=candidate;break;}
                    v=candidate;
                }
                require(std::isfinite(v) && v>=0 && v<=32.*recipe.scale,"invalid proposal");
                I z=I(std::ceil(v))+up(W(terms.size()+2)*q,p);
                auto slack=[&](I value) -> W {
                    W a=W(p)*value/q-fee;
                    for(auto t:terms) {
                        I diff=t.target-value;
                        a-=up(W(diff)*(diff>0?t.hi:t.lo),t.den);
                    }
                    return a;
                };
                W residual=slack(z);
                if(residual<0) {
                    z=narrow(W(z)+up((-residual+1)*q,p)); residual=slack(z);
                }
                require(residual>=0 && z>0 && z<=32*recipe.scale,"invalid row inequality");
                s.weights[D%2]=z; s.degrees[D%2]=D;
#ifdef FRIEND_TREE_OBSERVE_ROW
                FRIEND_TREE_OBSERVE_ROW(D, c, i, z);
#endif
                minimum=std::min(minimum,z); maximum=std::max(maximum,z);
                slack_min=std::min(slack_min,narrow(residual)); ++rows;
            }
        }
    }
    void run() {
        for(size_t i=0;i<bands.size();++i) {
            transitions(i);
            for(auto& layer:bands[i]->layers) layer.indices.clear();
        }
        for(size_t j=bands.size();j-->0;) {
            solve(j);
            std::cout<<"checked central degrees ["<<bands[j]->low<<","<<bands[j]->high<<")"<<std::endl;
        }
        auto& first=bands[0]->layers[1]; I gap=std::numeric_limits<I>::max();
        for(size_t i=0;i<recipe.rewards.size();++i) gap=std::min(gap,recipe.rewards[i]-first.get(i,2));
        require(gap>=0,"completed newborn weight exceeds its trial reward");
        std::cout<<"verified rows "<<rows<<", signed transfers "<<recipe.gains.size()
                 <<", minimum integer drift slack "<<slack_min<<"\n"
                 <<"scale "<<recipe.scale<<", minimum newborn reserve "<<gap<<", weight range ["<<minimum<<", "<<maximum
                 <<"], centre "<<bands[0]->layers[0].get(0,2)<<"\n"
                 <<"PASS upper "<<p<<"/"<<q<<std::endl;
    }
};
#ifndef FRIEND_TREE_NO_MAIN
int main(int argc,char** argv) {
    if(argc!=5) { std::cerr<<"usage: check geometry recipe numerator denominator\n"; return 2; }
    try { Checker(argv[1],argv[2],std::stoll(argv[3]),std::stoll(argv[4])).run(); }
    catch(const std::exception& e) { std::cerr<<"FAIL "<<e.what()<<"\n"; return 1; }
}
#endif
