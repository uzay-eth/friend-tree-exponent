// C++17. Each state is a central (L,c) and a sorted multiset of neighbour types.
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <iomanip>
#include <iostream>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>
using namespace std;
using Int = int64_t;
const Int scale = 10000000000LL;
const Int weight_limit = 16*scale;
Int quantize(double value) {
    if(!isfinite(value) || value<0 || value>double(weight_limit)/scale)
        throw runtime_error("invalid floating proposal");
    return llround(value*scale);
}
void weight_check(Int z) {
    if(z<0 || z>weight_limit) throw runtime_error("weight outside checked range");
}
using Bag = vector<int>;

constexpr int witness_cap=6;
struct Type { int d,k; int wd=0,wk=0; };
struct Move { int target, count, type; bool grow, witness=false; };
struct Row {
    Bag bag;
    vector<Move> moves;
    int convert = -1;
    double rho = 0;
    double extra = 0;
    vector<double> f;
    vector<Int> z;
};
struct Layer {
    int c, cap, pool, k_limit;
    vector<Type> types;
    vector<int> offset;
    vector<Row> rows;
    map<Bag,int> index;
    Layer(int c_,int cap_,int H):c(c_),cap(cap_),k_limit(c==1?20:c==2?7:5) {
        offset.resize(cap);
        for(int d=2;d<cap;++d) {
            offset[d]=types.size();
            for(int k=1;k<=min(d,k_limit);++k) {
                if(c==1 && d==2 && k==2) {
                    for(int w=2;w<witness_cap;++w)
                        for(int j=1;j<=w;++j) types.push_back({2,2,w,j});
                    types.push_back({2,2,witness_cap,0});
                } else types.push_back({d,k});
            }
        }
        pool=types.size();
        types.push_back({cap,0}); // pooled with only a lower degree bound
        for(int d=4;d<H;d*=2)
            if(d>cap) types.push_back({d,0});
        if(H>cap) types.push_back({H,0});
        Bag bag;
        enumerate(bag,0,H);
    }
    void enumerate(Bag& bag,int start,int H) {
        if((int)bag.size()==c) {
            index[bag]=rows.size();
            Row r;
            r.bag=bag;
            for(int t:bag) {
                if(t>=pool) r.extra+=1.0/types[t].d;
                else r.rho+=1.0/types[t].d;
            }
            r.f.resize(H-c);
            rows.push_back(move(r));
            return;
        }
        for(int t=start;t<(int)types.size();++t) {
            if(t>pool && (int)bag.size()+1<c) break;
            bag.push_back(t);
            enumerate(bag,t,H);
            bag.pop_back();
        }
    }
    int label(int d,int k) const {
        if(d>=cap) {
            int t=pool;
            while(t+1<(int)types.size() && types[t+1].d<=d) ++t;
            return t;
        }
        return offset[d]+min(k,k_limit)-1;
    }
    int bridge(int d,int k) const {
        // A degree-two, core-degree-two neighbour with an observed far endpoint.
        return d>=witness_cap ? offset[3]-1 : 1+(d-2)*(d+1)/2+k-1;
    }
    int lookup(Bag b) const {
        sort(b.begin(),b.end());
        // Keep only the greatest certified degree among pooled neighbours.
        for(int j=0;j+1<(int)b.size();++j) if(b[j]>pool) b[j]=pool;
        return index.at(b);
    }
};

struct Model {
    bool upper;
    int H;
    vector<Layer> layers;
    vector<double> rewards;
    vector<Int> iz;
    Int par_num,par_den;
    double par;
    Model(bool up,int h,const vector<int>& caps,Int a,Int b):
        upper(up), H(h),par_num(a),par_den(b),par(double(a)/b) {
        if(H<4 || H>512 || caps.empty() || H<caps[0] || caps[0]<witness_cap || a<=0 || a>=b || b>10000 || 4*a<b)
            throw runtime_error("invalid parameters");
        for(int c=0;c<H;++c) {
            int cap=c==0 ? 2 : c<=(int)caps.size()?caps[c-1]:2;
            if(cap<2 || (c>1 && cap>layers[c-1].cap)) throw runtime_error("caps must decrease");
            layers.emplace_back(c,cap,H);
        }
        for(auto &g:layers) for(auto &r:g.rows) {
            if(g.c+1<H) {
                auto &next=layers[g.c+1];
                Bag nb{next.label(2,1)};
                for(int t:r.bag) {
                    Type v=g.types[t];
                    nb.push_back(next.label(v.d,v.k));
                }
                r.convert=next.lookup(nb);
            }
            for(int j=0;j<g.c;) {
                int t=r.bag[j], end=j+1;
                while(end<g.c && r.bag[end]==t) ++end;
                if(t<g.pool) {
                    auto [d,k,wd,wk]=g.types[t];
                    for(bool grow:{true,false}) {
                        if(!grow && (k==d || k==g.k_limit)) continue;
                        Bag nb=r.bag;
                        nb[j]=grow?g.label(d+1,k):g.label(d,k+1);
                        r.moves.push_back({g.lookup(nb),end-j,t,grow});
                    }
                    if(wk) for(bool grow:{true,false}) {
                        if(!grow && wd==wk) continue;
                        Bag nb{grow?g.bridge(wd+1,wk):g.bridge(wd,wk+1)};
                        r.moves.push_back({g.lookup(nb),1,t,grow,true});
                    }
                }
                j=end;
            }
        }
        rewards.assign(layers[1].rows.size(),1);
    }
    int newborn(int D,int c,const Row& r) const {
        if(D==2 && c==1) {
            auto w=layers[1].types[r.bag[0]];
            return layers[1].bridge(w.d,w.k);
        }
        return layers[1].label(D,c+1);
    }
    // All D>=H have pooled newborn reward s. Upper tail is A*(1+par*L/D), A=s/(par*(1+par)).
    // Lower tail is B*L/D with B=s/(par+1/H).
    double tail(int L,int D) const {
        double s=rewards.back();
        return upper ? s/(par*(1+par))*(1+par*L/D) : s/(par+1.0/H)*L/D;
    }
    double opt(double lo,double hi,double target,double v) const {
        return ((target>v)==upper)?hi:lo;
    }
    void sweep() {
        for(int D=H-1;D>=2;--D) for(int c=D;c>=0;--c) {
            int L=D-c;
            auto &g=layers[c];
            for(int i=int(g.rows.size())-1;i>=0;--i) {
                auto &r=g.rows[i];
                double y=1.0/D, birth=L*y;
                double grow=D+1==H?tail(L+1,D+1):r.f[L+1];
                double made=L?layers[c+1].rows[r.convert].f[L-1]+rewards[newborn(D,c,r)]:0;
                double v=grow;
                for(int iter=0;;++iter) {
                    double h=L+r.rho+opt(0,r.extra,grow,v);
                    double num=h*grow+birth*made,den=par+h+birth;
                    for(auto m:r.moves) {
                        auto [d,k,wd,wk]=g.types[m.type];
                        double dest=g.rows[m.target].f[L];
                        double lo=k==g.k_limit?0:d-k;
                        double rate=m.grow ? y+opt(lo,d-k+(k-1)*0.5,dest,v) : double(d-k)/d;
                        if(m.witness) {
                            rate=m.grow ? wd-wk+0.5+opt(0,(wk-1)*0.5,dest,v) : double(wd-wk)/wd;
                        } else if(wd) {
                            rate=y+opt(wk?1.0/wd:0,1.0/wd,dest,v);
                        }
                        num+=m.count*rate*dest;
                        den+=m.count*rate;
                    }
                    double nv=num/den;
                    if(upper) nv=max(nv,1e-5);
                    if(abs(nv-v)<2e-14 || iter>=12) { r.f[L]=nv; break; }
                    v=nv;
                }
            }
        }
    }
    vector<double> births() const {
        vector<double> out;
        for(auto &r:layers[1].rows) out.push_back(r.f[1]);
        return out;
    }
    double iterate(int n) {
        double min_ratio=0,max_ratio=0;
        for(int j=0;j<n;++j) {
            sweep();
            auto next=births();
            min_ratio=1e10; max_ratio=0;
            for(int i=0;i<(int)next.size();++i) {
                double ratio=next[i]/rewards[i];
                min_ratio=min(min_ratio,ratio); max_ratio=max(max_ratio,ratio);
            }
            if(j+1==n || max_ratio-min_ratio<1e-10) break;
            double norm=next.back();
            for(int i=0;i<(int)next.size();++i) rewards[i]=next[i]/norm;
        }
        cout<<setprecision(12)<<"newborn ratios ["<<min_ratio<<", "<<max_ratio<<"]\n";
        return upper?max_ratio:min_ratio;
    }
    // Directed integer rounding of a signed rational. All stored potentials
    // share one denominator (scale), so it cancels out of drift inequalities.
    Int divide(__int128 a,Int b,bool ceiling) const {
        if(b<=0) throw runtime_error("nonpositive denominator");
        __int128 q=a/b,r=a%b;
        if(ceiling && r>0) ++q;
        if(!ceiling && r<0) --q;
        if(q>INT64_MAX || q<INT64_MIN) throw runtime_error("out of range");
        return Int(q);
    }
    Int coeff(Int x,Int a,Int b=1) const {
        return divide((__int128)x*a,b,upper);
    }
    Int tail_coeff=0;
    Int bound(int D,int c,int i,Int z,const vector<Int>& sv) const {
        int L=D-c;
        const auto &g=layers[c];
        const auto &r=g.rows[i];
        Int up=D+1<H?r.z[L+1]:(upper?tail_coeff+divide((__int128)tail_coeff*par_num*(L+1),par_den*H,true)
                                                       :coeff(tail_coeff,L+1,H));
        Int diff=up-z;
        __int128 drift=coeff(diff,L);
        for(int t:r.bag) {
            if(t>=g.pool) {
                if((diff>0)==upper) drift+=coeff(diff,1,g.types[t].d);
            } else drift+=coeff(diff,1,g.types[t].d);
        }
        if(L) drift+=coeff(layers[c+1].rows[r.convert].z[L-1]+sv[newborn(D,c,r)]-z,L,D);
        for(auto m:r.moves) {
            auto [d,k,wd,wk]=g.types[m.type];
            Int v=g.rows[m.target].z[L]-z;
            if(m.witness) {
                if(m.grow) {
                    Int twice=2*(wd-wk)+1+(((v>0)==upper)?wk-1:0);
                    drift+=coeff(v,twice,2);
                } else drift+=coeff(v,wd-wk,wd);
            } else if(wd) {
                drift+=coeff(v,1,D);
                if(wk || (v>0)==upper) drift+=coeff(v,1,wd);
            } else if(m.grow) {
                Int twice_D_rate=2;
                if((v>0)==upper) twice_D_rate+=(2*d-k-1)*D;
                else if(k!=g.k_limit) twice_D_rate+=2*(d-k)*D;
                drift+=coeff(v,twice_D_rate*m.count,2*D);
            } else drift+=coeff(v,m.count*(d-k),d);
        }
        Int compare=divide((__int128)par_num*z,par_den,!upper);
        return divide(upper?compare-drift:drift-compare,1,false);
    }
    void certify() {
        for(double v:rewards) iz.push_back(quantize(v));
        if(upper) tail_coeff=divide((__int128)iz.back()*par_den*par_den,par_num*(par_num+par_den),true);
        else tail_coeff=divide((__int128)iz.back()*par_den*H,par_num*H+par_den,false);
        weight_check(tail_coeff);
        for(auto &g:layers) for(auto &r:g.rows) r.z.resize(H-g.c);
        for(int D=H-1;D>=2;--D) for(int c=D;c>=0;--c) {
            auto &g=layers[c];
            for(int i=int(g.rows.size())-1;i>=0;--i) {
                auto &r=g.rows[i];
                Int &z=r.z[D-c];
                z=quantize(r.f[D-c]);
                for(;;) {
                    Int slack=bound(D,c,i,z,iz);
                    if(slack>=0) break;
                    Int step=divide((__int128)(1-slack)*par_den,par_num,true);
                    Int far=upper?min<Int>(weight_limit,z+step):max<Int>(0,z-step);
                    if(bound(D,c,i,far,iz)<0) throw runtime_error("weight limit in repair");
                    Int near=z;
                    while(abs(far-near)>1) {
                        Int mid=near+(far-near)/2;
                        if(bound(D,c,i,mid,iz)>=0) far=mid;
                        else near=mid;
                    }
                    z=far;
                }
            }
        }
        vector<Int> actual;
        Int consistency=INT64_MAX;
        for(int i=0;i<(int)iz.size();++i) {
            Int value=layers[1].rows[i].z[1];
            actual.push_back(value);
            consistency=min(consistency,upper?iz[i]-value:value-iz[i]);
        }
        cout<<"minimum newborn consistency "<<consistency<<" scale units\n";
        if(consistency<0) throw runtime_error("newborn consistency failed");
        Int min_slack=INT64_MAX, min_z=INT64_MAX, max_z=0, low_core_min=INT64_MAX;
        if(layers[2].cap<=2) throw runtime_error("c=2 must track degree two");
        int64_t count=0;
        for(int D=2;D<H;++D) for(int c=0;c<=D;++c) {
            auto &g=layers[c];
            for(int i=0;i<(int)g.rows.size();++i) {
                Int z=g.rows[i].z[D-c];
                if(z<0 || (upper && z==0)) throw runtime_error("positivity failed");
                Int residual=bound(D,c,i,z,actual);
                if(residual<0) throw runtime_error("actual newborn drift failed");
                min_z=min(min_z,z); max_z=max(max_z,z);
                bool tracked_B=c==2 && find(g.rows[i].bag.begin(),g.rows[i].bag.end(),
                                             g.label(2,2))!=g.rows[i].bag.end();
                if(c<=2 && (D>c || tracked_B)) low_core_min=min(low_core_min,z);
                min_slack=min(min_slack,residual);
                ++count;
            }
        }
        auto tail_lhs=(__int128)tail_coeff*(upper?par_num*(par_num+par_den):par_num*H+par_den);
        auto tail_rhs=(__int128)actual.back()*par_den*(upper?par_den:H);
        if(upper?tail_lhs<tail_rhs:tail_lhs>tail_rhs)
            throw runtime_error("hub-tail check failed");
        Int center=layers[0].rows[0].z[2];
        Int norm=upper?min(min_z,tail_coeff):max(max_z,tail_coeff);
        if(center<=0 || norm<=0) throw runtime_error("zero comparison constant");
        cout<<"verified "<<count<<" finite states and analytic D>="<<H
            <<" tail; minimum finite slack "<<min_slack<<" scale units\n";
        if(low_core_min<=0 || tail_coeff<=0) throw runtime_error("low-core minimum");
        Int tail_min_den=H*(upper?par_den:1);
        Int tail_min_num=upper?tail_coeff*(par_den*H+par_num*(H-2)):tail_coeff*(H-2);
        Int low_min_num=min(low_core_min*tail_min_den,tail_min_num);
        Int tail_max=upper?divide((__int128)tail_coeff*(par_den+par_num),par_den,true):tail_coeff;
        cout<<"weights in units of 1/"<<scale<<": finite ["<<min_z<<", "<<max_z
            <<"], tail coefficient "<<tail_coeff<<", global ["
            <<(upper?min(min_z,tail_coeff):0)<<", "<<max(max_z,tail_max)<<"]\n";
        cout<<"c<=2 specified subset: cap(2)="<<layers[2].cap<<", finite minimum "
            <<low_core_min<<"/"<<scale<<", tail minimum "<<tail_min_num<<"/"<<tail_min_den*scale
            <<", combined minimum "<<low_min_num<<"/"<<tail_min_den*scale<<"\n";
        cout<<"PASS "<<(upper?"upper":"lower")<<" "<<par_num<<"/"<<par_den
            <<" product prefactor "<<center<<"/"<<norm<<"\n";
    }

};

int main(int argc,char **argv) {
    if(argc<7) { cerr<<"usage: hub upper|lower numerator denominator H iterations cap1 [cap2 ...]\n"; return 2; }
    if(string(argv[1])!="upper" && string(argv[1])!="lower") return 2;
    bool up=string(argv[1])=="upper";
    vector<int> caps;
    for(int i=6;i<argc;++i) caps.push_back(stoi(argv[i]));
    Model m(up,stoi(argv[4]),caps,stoll(argv[2]),stoll(argv[3]));
    int64_t count=0;
    for(auto &g:m.layers) count+=g.rows.size()*(m.H-max(2,g.c));
    cout<<"states "<<count<<"\n";
    double ratio=m.iterate(stoi(argv[5]));
    if(up?ratio<=1:ratio>=1) m.certify();
    else { cout<<"screen failed\n"; return 1; }
}
