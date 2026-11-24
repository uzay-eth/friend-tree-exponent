// Discounted lifetime credit for real donors adjacent to four-core motifs.
#include "../motif/cert4v_format.h"
#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

using I = std::int64_t;
using Wide = __int128_t;
constexpr int Max = 5;
using Five = std::array<int, Max>;
using Four = std::array<int, 4>;
using Ids = std::array<I, Max>;

void require(bool condition, const char* msg) {
    if (!condition) throw std::runtime_error(msg);
}

I down(Wide n, Wide d) {
    require(d > 0, "nonpositive denominator");
    Wide k = n / d - (n % d < 0);
    require(k >= std::numeric_limits<I>::min() && k <= std::numeric_limits<I>::max(),
            "fixed-point integer overflow");
    return static_cast<I>(k);
}

struct Type {
    int degree;
    int cores;
    int leaves() const { return degree - cores; }
    int children() const { return cores - 1; }
    bool pooled() const { return cores == 0; }
};

enum Family { Path, Fork, Y, Trident };
struct Piece {
    int donor, contact;
    Family family;
    Four roles;
};

// The member order is chronological, and every noninitial member has an
// earlier parent in this subset. A fifth donor must be an undirected leaf.
std::vector<Piece> pieces(const Five& parent) {
    std::vector<Piece> result;
    for (int del = 0; del < 5; ++del) {
        std::vector<int> adj;
        if (del) adj.push_back(parent[del]);
        for (int j = del + 1; j < 5; ++j) if (parent[j] == del) adj.push_back(j);
        if (adj.size() != 1) continue;
        std::array<std::vector<int>, 5> kids;
        for (int j = 1; j < 5; ++j)
            if (j != del && parent[j] != del) kids[parent[j]].push_back(j);
        int top = del == 0 ? adj[0] : 0;
        Four r{};
        Family f;
        if (kids[top].size() == 3) {
            f = Trident;
            r = {top, kids[top][0], kids[top][1], kids[top][2]};
        } else if (kids[top].size() == 2) {
            f = Fork;
            int mid = kids[top][0], side = kids[top][1];
            if (kids[mid].empty()) std::swap(mid, side);
            require(kids[mid].size() == 1, "fork middle");
            r = {top, mid, side, kids[mid][0]};
        } else {
            require(kids[top].size() == 1, "root ordering");
            int mid = kids[top][0];
            if (kids[mid].size() == 2) {
                f = Y;
                r = {top, mid, kids[mid][0], kids[mid][1]};
            } else {
                f = Path;
                require(kids[mid].size() == 1 && kids[kids[mid][0]].size() == 1, "chain");
                r = {top, mid, kids[mid][0], kids[kids[mid][0]][0]};
            }
        }
        auto where = std::find(r.begin(), r.end(), adj[0]);
        require(where != r.end(), "lost contact");
        result.push_back({del, static_cast<int>(where - r.begin()), f, r});
    }
    return result;
}

class Tables {
    const unsigned char* mapping_;
    size_t bytes_;
    Cert4vLayout layout_{};
    std::array<const I*, 4> weights_{};
    const I* vertex_ = nullptr;
public:
    Cert4vHeader h{};
    explicit Tables(const char* path) {
        int fd = open(path, O_RDONLY);
        require(fd >= 0, "cannot open four-core witness");
        struct stat st{};
        require(fstat(fd, &st) == 0 && st.st_size >= static_cast<I>(sizeof h), "invalid witness size");
        bytes_ = st.st_size;
        void* ptr = mmap(nullptr, bytes_, PROT_READ, MAP_PRIVATE, fd, 0);
        close(fd);
        require(ptr != MAP_FAILED, "cannot map witness");
        mapping_ = static_cast<const unsigned char*>(ptr);
        std::memcpy(&h, mapping_, sizeof h);
        require(c4_validate_header(&h) == 0 && h.S == (I(1) << 40), "invalid witness header");
        c4_init_layout(&layout_, &h);
        Cert4vSizes sizes;
        c4_compute_sizes(&layout_, &sizes);
        require(bytes_ == sizes.total_bytes, "witness body has wrong size");
        const I* weights = reinterpret_cast<const I*>(mapping_ + sizeof h);
        vertex_ = weights + sizes.n_cent + sizes.n_fy;
        size_t at = sizes.n_cent + sizes.n_fy + sizes.n_f + sizes.n_g + sizes.n_hV + sizes.n_hC;
        weights_[Path] = weights + at; at += sizes.n_hP4;
        weights_[Fork] = weights + at; at += sizes.n_hF4;
        weights_[Y] = weights + at; at += sizes.n_hY4;
        weights_[Trident] = weights + at;
    }
    Tables(const Tables&) = delete;
    ~Tables() { c4_free_layout(&layout_); munmap(const_cast<unsigned char*>(mapping_), bytes_); }

    I vertex_gain(Type t) const {
        if (t.pooled() || t.degree < h.D0_2 || t.degree + 1 >= h.DF) return 0;
        int k = t.children(), leaves = t.leaves();
        I a = vertex_[c4_vidx(&layout_, leaves, k)];
        I b = vertex_[c4_vidx(&layout_, leaves + 1, k)];
        require(-16*h.S<a && a<16*h.S && -16*h.S<b && b<16*h.S, "vertex magnitude");
        return b-a;
    }

    I value(Family f, const std::array<Type, 4>& t) const {
        const auto& l = layout_;
        const auto& s = l.L4;
        auto id = [&](int a) {
            require(!t[a].pooled() && t[a].degree < s.D0, "input outside intended interior");
            return c4_lvl_inner(&s, t[a].leaves(), t[a].children());
        };
        auto L = [&](int a) { return t[a].leaves(); };
        auto k = [&](int a) { return t[a].children(); };
        int j;
        size_t address;
        switch (f) {
            case Path:
                j = c4_p4_act(&l, L(1), k(1), L(2), k(2));
                if (j < 0) return 0;
                address = c4_p4_off(&l, j, id(0), id(3));
                break;
            case Fork:
                j = c4_f4_act(&l, L(0), k(0), L(1), k(1));
                if (j < 0) return 0;
                address = c4_f4_off(&l, j, id(2), id(3));
                break;
            case Y:
                j = c4_y4_act(&l, L(1), k(1));
                if (j < 0) return 0;
                address = c4_y4_off(&l, j, id(0), id(2), id(3));
                break;
            case Trident:
                j = c4_t4_act(&l, L(0), k(0));
                if (j < 0) return 0;
                address = c4_t4_off(&l, j, id(1), id(2), id(3));
                break;
            default: throw std::runtime_error("invalid four-core shape");
        }
        I v = weights_[f][address];
        require(v > -16*h.S && v < 16*h.S, "witness magnitude outside safe range");
        return v;
    }
};

struct Axis {
    int dlim, clim;
    std::vector<Type> values;
    std::vector<std::vector<int>> lookup;
    int pool_id = -1;
    Axis(int dmax, int cmax, int children) :
        dlim(dmax), clim(cmax), lookup(dlim+1, std::vector<int>(clim+1, -1)) {
        for (int d=2; d<=dlim; ++d) {
            for (int c=1+children; c<=std::min(d,clim); ++c) {
                lookup[d][c] = static_cast<int>(values.size());
                values.push_back({d,c});
            }
        }
        if (children <= 1) {
            pool_id = static_cast<int>(values.size());
            values.push_back({std::min(dlim,clim)+1, 0});
        }
    }
    int index(Type t) const {
        if (t.pooled() || t.degree > dlim || t.cores > clim) return pool_id;
        return lookup[t.degree][t.cores];
    }
};

struct Grid {
    int n;
    Five parent, children{}, stride{};
    std::vector<Axis> axes;
    std::vector<I> val;
    Grid(int size, Five tree, Five dlim, Five clim) : n(size), parent(tree) {
        for (int j=1; j<n; ++j) ++children[parent[j]];
        for (int i=0; i<n; ++i) axes.emplace_back(dlim[i],clim[i],children[i]);
        size_t length = 1;
        for (int i=n-1; i>=0; --i) {
            stride[i] = static_cast<int>(length);
            length *= axes[i].values.size();
            require(length < 120000000, "grid exceeds memory allowance");
        }
        val.assign(length,0);
    }
    I get(const std::array<Type,5>& entry) const {
        int index = 0;
        for (int i=0; i<n; ++i) {
            int a = axes[i].index(entry[i]);
            if (a < 0) return 0;
            index += a*stride[i];
        }
        return val[index];
    }
};

class Reserves {
    const Tables& tab_;
    Five dmax_, cmax_;
    int sign_, sp_, sq_, last_, depth_, target_;
    std::array<std::uint64_t, 6> rows_{};
    // Deltas have identical values whenever the shape, role and its four
    // (degree, core-degree) pairs agree. One compact array handles every donor.
    int packed_width_;
    std::vector<std::array<I,4>> deltas_;
    std::vector<char> ready_;
    std::ofstream dump_;
    int code(Type t) const {
        return (t.degree-1)*t.degree/2 - 1 + t.cores-1;
    }
    I difference(const Piece& p, const std::array<Type,5>& x) {
        std::array<Type,4> types{};
        size_t key = p.family;
        for (int i=0;i<4;++i) {
            types[i] = x[p.roles[i]];
            if (types[i].pooled()) return 0;
            key = key * packed_width_ + code(types[i]);
        }
        if (!ready_[key]) {
            I base = tab_.value(p.family, types);
            for (int i=0;i<4;++i) {
                ++types[i].degree;
                deltas_[key][i] = sign_*(tab_.value(p.family,types)-base);
                --types[i].degree;
            }
            ready_[key]=true;
        }
        return deltas_[key][p.contact];
    }
    I force(const std::vector<Piece>& parts, const std::array<Type,5>& x) {
        I sum=0;
        for (const auto& p:parts) {
            I d=difference(p,x);
            const Type donor=x[p.donor];
            if (d>=0) sum+=down(Wide(d)*(donor.degree-2),2*donor.degree);
            else if (!donor.pooled()) sum+=down(-Wide(d),donor.degree);
        }
        return sum;
    }

    Grid solve(int n,Five tree) {
        Five dlimits{}, climits{};
        dlimits.fill(dmax_[n-1]); climits.fill(cmax_[n-1]);
        if (depth_ == 2 && n == 2) {
            dlimits[1-target_] = dmax_[2]; climits[1-target_] = cmax_[2];
        }
        Grid grid(n,tree,dlimits,climits);
        if (grid.val.empty()) return grid;
        std::vector<Grid> expanded;
        if (n<depth_) for (int i=0;i<n;++i) {
            tree[n]=i;
            expanded.push_back(solve(n+1,tree));
        }
        std::vector<Piece> parts;
        if (n==5) parts=pieces(tree);
        for (int row=static_cast<int>(grid.val.size())-1;row>=0;--row) {
            std::array<Type,5> x{};
            Five pos{};
            I q=2;
            for (int i=0;i<n;++i) {
                pos[i]=(row/grid.stride[i])%grid.axes[i].values.size();
                x[i]=grid.axes[i].values[pos[i]];
                q*=x[i].degree;
            }
            require(Wide(q)*(sp_+Wide(sq_)*(5*dmax_[0]+10)) < (Wide(1)<<62),
                    "rate denominator exceeds supported integer range");
            I source=(n==5)?force(parts,x):0;
            if (depth_ == 2 && n == 2) {
                const Type donor=x[1-target_];
                I change=sign_*tab_.vertex_gain(x[target_]);
                if (change>=0) source=down(Wide(change)*(donor.degree-2),2*donor.degree);
                else if (!donor.pooled()) source=down(-Wide(change),donor.degree);
            }
            I base_den=sp_*q;
            Wide base_num=Wide(sq_)*q*source;
            Ids gain_to{}, lo{}, hi{};
            for (int i=0;i<n;++i) {
                const Type t=x[i];
                if (t.pooled()) continue;
                I known=t.leaves()*q, unobserved=0;
                auto neighbor=[&](int j) {
                    if(x[j].pooled()) unobserved+=q/x[j].degree;
                    else known+=q/x[j].degree;
                };
                if (i>0) neighbor(tree[i]);
                for (int j=i+1;j<n;++j) if (tree[j]==i) neighbor(j);
                int free=t.cores-grid.children[i]-(i!=0);
                require(free>=0,"negative remaining core slots");
                lo[i]=known;
                hi[i]=known+unobserved+free*q/2;
                auto old=[&](Type t2) {
                    int target=grid.axes[i].index(t2);
                    require(target==-1 || target>pos[i],"nonadvancing old state");
                    return target<0?I(0):grid.val[row+(target-pos[i])*grid.stride[i]];
                };
                gain_to[i]=old({t.degree+1,t.cores});
                if (t.leaves()) {
                    Type after{t.degree,t.cores+1};
                    I value=old(after);
                    if (n<depth_) {
                        auto entry=x;
                        entry[i]=after;
                        entry[n]={2,1};
                        value+=expanded[i].get(entry);
                    }
                    I bud=t.leaves()*(q/t.degree);
                    base_den+=sq_*bud;
                    base_num+=Wide(sq_)*bud*value;
                }
            }
            // The proposal selects a nearby integer. The last residual is
            // checked at every unknown input's true worst endpoint.
            double trial=static_cast<double>(base_num)/base_den;
            for (int it=0;it<=n;++it) {
                double num=static_cast<double>(base_num), den=base_den;
                for(int i=0;i<n;++i) {
                    double r=static_cast<double>(sq_)*(gain_to[i]<trial?hi[i]:lo[i]);
                    num+=r*gain_to[i]; den+=r;
                }
                double next=num/den;
                if (std::abs(next-trial)<0.1) {trial=next;break;}
                trial=next;
            }
            require(std::isfinite(trial) && trial>=0 && trial<1.e18,"unsafe trial");
            I z=std::min<I>(tab_.h.S,static_cast<I>(trial));
            auto residual=[&](I candidate) {
                Wide v=base_num-Wide(base_den)*candidate;
                for(int i=0;i<n;++i) {
                    I jump=gain_to[i]-candidate;
                    v+=Wide(sq_)*jump*(jump<0?hi[i]:lo[i]);
                }
                return v;
            };
            Wide r=residual(z);
            if(r<0) z=std::max<I>(0,z+down(r,base_den));
            require(residual(z)>=0,"failed exact occupation inequality");
            grid.val[row]=z;
            if (dump_.is_open()) {
                for (int j=1;j<n;++j) dump_ << tree[j];
                if(n==1) dump_ << '-';
                for(int i=0;i<n;++i) dump_ << " " << x[i].degree << " " << x[i].cores;
                dump_ << " " << z << " " << source << "\n";
            }
        }
        rows_[n]+=grid.val.size();
        return grid;
    }
public:
    Reserves(const Tables& t,Five dlim,Five clim,int sp,int sq,int target):
        tab_(t),dmax_(dlim),cmax_(clim),sign_(t.h.mode==1?1:-1),sp_(sp),sq_(sq),
        last_(dlim[4]),depth_(target<0?5:2),target_(target) {
        for(int i=0;i<5;++i) {
            require(dmax_[i]>=2 && dmax_[i]<=1400 && cmax_[i]>=2 && cmax_[i]<=dmax_[i],
                    "invalid state cutoffs");
            if(i) require(dmax_[i-1]>=dmax_[i] && cmax_[i-1]>=cmax_[i],"cutoffs must be nonincreasing");
        }
        require(sp>0 && sq>0 && sp<=1000000 && sq<=1000000,"invalid discount");
        if (const char* dump_path=std::getenv("FT2_RESERVE_DUMP")) {
            dump_.open(dump_path);
            require(dump_.good(), "cannot create diagnostic dump");
        }
        if (depth_==2) return;
        require(last_+1<t.h.D0_4 && last_<=9,"last layer must be strictly interior");
        packed_width_=last_*(last_+1)/2-1;
        size_t n=4;
        for(int j=0;j<4;++j) n*=packed_width_;
        require(n<=50000000,"delta cache exceeds memory allowance");
        deltas_.resize(n);
        ready_.resize(n);
    }
    void run() {
        Five root{-1,0,0,0,0};
        auto g=solve(1,root);
        std::array<Type,5> initial{};
        initial[0]={2,1};
        I v=g.get(initial);
        std::cout<<"mode "<<tab_.h.mode;
        if (depth_==2) std::cout<<(target_==0?" parent":" child")<<" vertex slack";
        std::cout<<" degree caps";
        for(auto x:dmax_)std::cout<<" "<<x;
        std::cout<<" core caps";
        for(auto x:cmax_)std::cout<<" "<<x;
        std::cout<<" rows";
        for(int i=1;i<=5;++i)std::cout<<" "<<rows_[i];
        std::cout<<" newborn reserve "<<v<<"/"<<tab_.h.S<<" at discount "<<sp_<<"/"<<sq_<<"\n";
    }
};

Five parse(const char* list) {
    Five v{};
    std::string input(list);
    for(int i=0;i<5;++i) {
        size_t pos=input.find(',');
        require((pos==std::string::npos)==(i==4),"expected five comma-separated cutoffs");
        v[i]=std::stoi(input.substr(0,pos));
        if(i<4)input.erase(0,pos+1);
    }
    return v;
}

int main(int argc,char** argv) {
    try {
        require(argc==2 || argc==4 || argc==6 || argc==7,"usage: five_reserve WITNESS [D1,..,D5 C1,..,C5 [s_num s_den [parent|child]]]");
        Tables table(argv[1]);
        Five d{32,14,7,5,4},c{12,7,5,4,4};
        if(argc>=4) {d=parse(argv[2]);c=parse(argv[3]);}
        int sp=argc>=6?std::stoi(argv[4]):1133, sq=argc>=6?std::stoi(argv[5]):2000;
        int target=-1;
        if(argc==7) {
            const std::string side(argv[6]);
            require(side=="parent" || side=="child", "invalid high vertex");
            target=side=="parent"?0:1;
        }
        Reserves reserves(table,d,c,sp,sq,target);
        reserves.run();
    } catch (const std::exception& error) {
        std::cerr<<error.what()<<"\n";return 1;
    }
}
