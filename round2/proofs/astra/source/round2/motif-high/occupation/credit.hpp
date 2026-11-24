// Additional occupation sources on physically specified noninitial cores.
// Header-only C++17. Link with -lgmpxx -lgmp. The snapshot is opened read-only.
#pragma once
#include "../rows.hpp"
#include <memory>
#include <type_traits>

namespace mh_occ {
using mh::Family;
using mh::I;
using mh::Labels;
using mh::Wide;
using Q = mpq_class;

// core_degree counts ALL core neighbours, including the parent. core_degree=0
// marks a pooled member and makes every source on the selected set zero.
struct Core { int degree, core_degree; };
using Members = std::array<Core,4>;
enum class Direction { Parent, Child };

struct Shape {
    Family family;
    Labels role;   // role[j] = input position of frozen table coordinate j
    static Shape from_parents(int n, const int *parent) {
        mh::require(n == 3 || n == 4, "occupation shape must have 3 or 4 members");
        for (int i = 1; i < n; ++i)
            mh::require(parent[i] >= 0 && parent[i] < i,"invalid chronological parent");
        for (int f=mh::V; f<=mh::T; ++f) {
            auto &form=mh::shapes[f];
            if (form.n != n) continue;
            // The chronological apex is 0. Symmetric children have the first matching order.
            Labels m{0,1,2,3};
            do {
                bool same=true;
                for (int j=1;j<n;++j) if (parent[m[j]] != m[form.parent[j]]) same=false;
                if (same) return {Family(f),m};
            } while (std::next_permutation(m.begin()+1,m.begin()+n));
        }
        throw std::logic_error("no rooted shape");
    }
    static Shape canonical(Family f) {
        mh::require(f>=mh::V && f<=mh::T,"invalid higher family");
        return {f,{0,1,2,3}};
    }
};

inline I floor_q(const Q &a) {
    mpz_class n;
    mpz_fdiv_q(n.get_mpz_t(),a.get_num_mpz_t(),a.get_den_mpz_t());
    mh::require(n.fits_slong_p(),"occupation source exceeds int64 range");
    return n.get_si();
}
inline Q fraction(I n,I d) {
    Q v(n,d);
    v.canonicalize();
    return v;
}

inline I donor_floor(I d,int degree,bool pooled) {
    if (d<0) return pooled ? 0 : (-d)/degree;
    return I(Wide(d)*(degree-2)/(2LL*degree));
}
inline I donor_floor(const Q &d,int degree,bool pooled) {
    if (d<0) return pooled ? 0 : floor_q(-d/degree);
    return floor_q(d*fraction(degree-2,2LL*degree));
}

struct Prepared {
    bool active = false;
    int arity = 0;
    I residual_ceiling = 0;  // exact ceil of signed R_phys,max in integer snapshot units
    I intrinsic = 0;         // exact floor of max(0,-R_phys,max)
    Labels free_parent{}, free_children{};
    std::array<I,4> gains{};
    std::unique_ptr<std::array<Q,4>> fractional_gains;
    // One additional J from an uncovered outside donor of a TRIPLE. For a
    // pooled donor, degree is a proven lower bound on its actual degree.
    // At quartets both free arrays stay zero: their outside J is already booked.
    I triple_donor(int contact, Direction dir, int degree, bool pooled=false) const {
        if (!active) return 0;
        mh::require(contact>=0 && contact<arity && degree>=2,"invalid donor descriptor");
        int slots=(dir==Direction::Parent?free_parent:free_children)[contact];
        if (!slots) return 0;
        return fractional_gains ? donor_floor((*fractional_gains)[contact],degree,pooled)
                                : donor_floor(gains[contact],degree,pooled);
    }
};

class Source {
    mh::Snapshot snap;
    int P,Qden;
    struct Geometry {
        int n=0, ncut=0;
        Labels children{},degree{};
        std::array<Labels,4> neighbors{};
        std::array<mh::Cut,3> cuts{};
        std::array<mh::Extension,4> child{};
        mh::Extension parent{};
        explicit Geometry(Family f) {
            const auto &shape=mh::shapes[f];
            n=shape.n;
            for (int i=1;i<n;++i) {
                int p=shape.parent[i];
                ++children[p];
                neighbors[i][degree[i]++]=p; neighbors[p][degree[p]++]=i;
            }
            for (int i=0;i<n;++i) if (degree[i]==1) cuts[ncut++]={i,neighbors[i][0]};
            if (f==mh::V) {
                parent={mh::P,{3,0,1,2}};
                child[0]={mh::F,{0,3,1,2}}; child[1]={mh::Y,{0,1,2,3}}; child[2]={mh::P,{0,1,2,3}};
            } else if (f==mh::C) {
                parent={mh::Y,{3,0,1,2}};
                child[0]={mh::T,{0,1,2,3}}; child[1]={mh::F,{0,2,1,3}}; child[2]={mh::F,{0,1,2,3}};
            }
        }
    };
    std::array<Geometry,6> plans{Geometry(mh::V),Geometry(mh::C),Geometry(mh::P),
                                 Geometry(mh::F),Geometry(mh::Y),Geometry(mh::T)};
    const Cert4vLevel &level(Family f) const {
        return f==mh::Edge ? snap.c.L2 : snap.level(f);
    }
    static int L(Core v) { return v.degree-v.core_degree; }
    static int k(Core v) { return v.core_degree-1; }
    static constexpr Core newborn{2,1};
    static int point(const Cert4vLevel &g,Core v) {
        return v.degree<g.D0 ? c4_lvl_inner(&g,L(v),k(v)) : c4_lvl_tail(&g,k(v));
    }
    static Labels points(const Cert4vLevel &g, const Members &s, int n) {
        Labels x{};
        for (int i=0;i<n;++i) x[i]=point(g,s[i]);
        return x;
    }
    int active(Family f,const Members &s) const {
        if (f==mh::Edge) return 0;
        return snap.active(f,points(level(f),s,mh::shapes[f].n));
    }
    I table(Family f, int a, const Labels &x) const {
        if (f==mh::Edge) return snap.sign*snap.weight[f][size_t(x[0])*snap.c.L2.NG+x[1]];
        return snap.value(f,a,x);
    }
    static Members extend(const mh::Extension &e,const Members &s) {
        Members out{};
        for (int j=0;j<4;++j) out[j]=e.from[j]==3 ? newborn : s[e.from[j]];
        return out;
    }
    struct Prune { Family f; Labels take; };
    static Prune prune(Family f, int cut) {
        switch (f) {
        case mh::V: return cut==0 ? Prune{mh::Edge,{1,2}} : Prune{mh::Edge,{0,1}};
        case mh::C: return {mh::Edge,{0,3-cut}};
        case mh::P: return cut==0 ? Prune{mh::V,{1,2,3}} : Prune{mh::V,{0,1,2}};
        case mh::F: return cut==1 ? Prune{mh::V,{0,2,3}} : Prune{mh::C,{0,1,2}};
        case mh::Y: return cut==0 ? Prune{mh::C,{1,2,3}} : Prune{mh::V,{0,1,5-cut}};
        case mh::T: return {mh::C,{0,cut==1?2:1,cut==3?2:3}};
        default: throw std::logic_error("bad pruning family");
        }
    }
    // With all original members strictly inside their level, every retained
    // self and lower weight uses a point label; a gain at the seam lands on T1.
    I exact_point(Family f, const Members &s) const {
        const auto &g=level(f);
        auto x=points(g,s,mh::shapes[f].n);
        return table(f,f==mh::Edge?0:snap.active(f,x),x);
    }
    template<class Numerator>
    static Numerator numerator(I value) {
        if constexpr (std::is_same_v<Numerator,Wide>) return Wide(value);
        else return mpz_class(value);
    }
    // Tensor interpolation. The fast path is used only for active newborn
    // quartets with <=2 nontrivial factors of denominator <=127.
    template<class Z> std::pair<Z,Z> average(Family f, const Members &s) const {
        const auto &g=level(f);
        int n=mh::shapes[f].n;
        Labels x=points(g,s,n);
        int a=f==mh::Edge ? 0 : snap.active(f,x);
        if (a<0) return {numerator<Z>(0),numerator<Z>(1)};
        int var[3]{}, nv=0;
        for (int j=0;j<n;++j) if (s[j].degree>g.D0) {
            if (nv==3) throw std::logic_error("active motif has more than three tail coordinates");
            var[nv++]=j;
        }
        if (!nv) return {numerator<Z>(table(f,a,x)),numerator<Z>(1)};
        Labels original=x;
        Z sum=numerator<Z>(0), den=numerator<Z>(1);
        for (int j=0;j<nv;++j) den*=s[var[j]].degree;
        for (int mask=0;mask<(1<<nv);++mask) {
            Z mult=numerator<Z>(1);
            for (int j=0;j<nv;++j) {
                int v=var[j]; bool inf=(mask & (1<<j)) !=0;
                x[v]=inf ? g.INF_ID : original[v];
                mult*=inf ? s[v].degree-g.D0 : g.D0;
            }
            sum+=mult*table(f,a,x);
        }
        return {sum,den};
    }
    Q exact_weight(Family f,const Members &s) const {
        auto [n,d]=average<mpz_class>(f,s);
        Q v(n,d); v.canonicalize(); return v;
    }
    template<bool Fast> using Value=std::conditional_t<Fast,I,Q>;
    template<bool Fast> Value<Fast> weight(Family f,const Members &s) const {
        if constexpr (Fast) return exact_point(f,s);
        else return exact_weight(f,s);
    }
    template<bool Fast,class Sum>
    static void add(Sum &sum, const Value<Fast> &v, I n=1,I d=1) {
        if constexpr (Fast) sum.add(Wide(v)*n,d);
        else if (n && v!=0) sum+=v*fraction(n,d);
    }
    template<bool Fast,class Sum>
    void add_newborn(Sum &sum, const mh::Extension &e,const Members &s,Core old) const {
        Members grown=extend(e,s);
        if constexpr (Fast) {
            auto [n,d]=average<Wide>(e.f,grown);
            sum.add(n*L(old),I(d)*old.degree);
        } else add<false>(sum,exact_weight(e.f,grown),L(old),old.degree);
    }

    template<bool Fast,class Sum>
    void drift(const Shape &shape, Members s, Prepared &out, Sum &sum) const {
        Family f=shape.family;
        const auto &p=plans[f-mh::V];
        auto w=weight<Fast>(f,s);
        add<Fast>(sum,w,-P,Qden);
        for (int j=0;j<p.ncut;++j) {
            auto [y,x]=p.cuts[j];
            auto small=prune(f,y);
            Members r{}; int at=-1;
            for (int z=0;z<p.n-1;++z) {
                r[z]=s[small.take[z]];
                if (small.take[z]==x) at=z;
            }
            if (at<0) throw std::logic_error("pruning lost the contact");
            auto old=weight<Fast>(small.f,r);
            ++r[at].degree;
            Value<Fast> diff=weight<Fast>(small.f,r)-old;
            bool parent=mh::shapes[f].parent[x]==y;
            int cent=f==mh::C ? 4 : f==mh::V ? (parent?2:3) : (parent?0:1);
            I z=snap.cent[cent][point(snap.c.L2,s[x])];
            add<Fast>(sum,diff,1,s[y].degree);
            add<Fast>(sum,diff,-z,snap.h.S_r);
        }
        for (int i=0;i<p.n;++i) {
            Core old=s[i];
            ++s[i].degree;
            Value<Fast> d=weight<Fast>(f,s)-w;
            s[i]=old;
            if constexpr (Fast) out.gains[shape.role[i]]=d;
            else (*out.fractional_gains)[shape.role[i]]=d;
            add<Fast>(sum,d,L(old));
            for (int j=0;j<p.degree[i];++j) add<Fast>(sum,d,1,s[p.neighbors[i][j]].degree);
            int open=0;
            if (i==0) {
                if (p.n==3 && active(p.parent.f,extend(p.parent,s))>=0)
                    add<Fast>(sum,d,snap.cent[0][point(snap.c.L2,old)],snap.h.S_r);
                else {
                    ++open;
                    if (p.n==3) out.free_parent[shape.role[i]]=1;
                }
            }
            int children=k(old)-p.children[i];
            if (children) {
                if (p.n==3 && active(p.child[i].f,extend(p.child[i],s))>=0) {
                    // Covered contacts are internal to a supported quartet and
                    // their child counts are below the finite interior cutoff.
                    add<Fast>(sum,d,I(children)*snap.cent[1][point(snap.c.L2,old)],snap.h.S_r);
                } else {
                    open+=children;
                    if (p.n==3) out.free_children[shape.role[i]]=children;
                }
            }
            if (d>0) add<Fast>(sum,d,open,2);
            if (L(old)) {
                ++s[i].core_degree;  // degree is fixed, L goes down and k goes up
                Value<Fast> db=weight<Fast>(f,s)-w;
                add<Fast>(sum,db,L(old),old.degree);
                if (p.n==3) add_newborn<Fast>(sum,p.child[i],s,old);
                s[i]=old;
            }
        }
    }
    static void validate_shape(const Shape &sh) {
        mh::require(sh.family>=mh::V && sh.family<=mh::T,"invalid occupation family");
        int n=mh::shapes[sh.family].n,used=0;
        for (int j=0;j<n;++j) {
            int role=sh.role[j];
            mh::require(role>=0 && role<n && !(used & (1<<role)),"invalid occupation role map");
            used |= 1<<role;
        }
    }
    bool setup(const Shape &sh,const Core *in,Members &s,bool &fast) const {
        validate_shape(sh);
        const auto &p=plans[sh.family-mh::V];
        fast=true;
        for (int i=0;i<p.n;++i) {
            s[i]=in[sh.role[i]];
            auto v=s[i];
            if (!v.core_degree) return false;
            mh::require(v.degree>=2 && v.degree<INT_MAX && v.core_degree>=1+p.children[i] &&
                        v.core_degree<=v.degree,"invalid physical motif state");
            fast=fast && v.degree<level(sh.family).D0;
        }
        return active(sh.family,s)>=0;
    }
public:
    Source(const char *path,int num,int den):snap(path),P(num),Qden(den) {
        mh::require(P>0 && Qden>0,"nonpositive occupation target");
        I limit=I(1)<<52;
        for (size_t i=0;i<snap.z.n_g;++i) {
            I w=snap.weight[mh::Edge][i];
            if (w < -limit || w > limit) throw std::runtime_error("edge weight outside arithmetic range");
        }
    }
    I scale() const { return snap.h.S; }
    int sign() const { return snap.sign; }
    Prepared prepare(const Shape &sh,const Core *in) const {
        Prepared out;
        Members s{}; bool fast;
        if (!setup(sh,in,s,fast)) return out;
        out.active=true; out.arity=mh::shapes[sh.family].n;
        if (fast) {
            mh::Grid g;
            drift<true>(sh,s,out,g);
            I hi=mh::Grid::ceil_grid(g.up), lo=mh::Grid::ceil_grid(g.up-g.noninteger);
            if (hi!=lo) {
                mh::Rational r;
                drift<true>(sh,s,out,r);
                I exact=r.ceil();
                mh::require(lo<=exact && exact<=hi,"occupation filter disagrees with GMP");
                hi=exact;
            }
            out.residual_ceiling=hi;
        } else {
            out.fractional_gains=std::make_unique<std::array<Q,4>>();
            Q r=0; drift<false>(sh,s,out,r);
            I neg=floor_q(-r);
            mh::require(neg>INT64_MIN,"occupation residual outside int64 range");
            out.residual_ceiling=-neg;
        }
        mh::require(out.residual_ceiling>INT64_MIN,"occupation slack outside int64 range");
        out.intrinsic=std::max<I>(0,-out.residual_ceiling);
        return out;
    }
    I row(const Shape &sh,const Core *in) const { return prepare(sh,in).intrinsic; }
    I triple_donor(const Shape &sh,const Core *in,int contact,Direction dir,int degree,bool pooled=false) const {
        if (sh.family>mh::C) return 0;  // quartet outside credits belong to the existing five-core source
        mh::require(contact>=0 && contact<3 && degree>=2,"invalid triple donor");
        Members s{}; bool fast;
        if (!setup(sh,in,s,fast)) return 0;
        int i=0;
        while (sh.role[i]!=contact) ++i;
        const auto &p=plans[sh.family-mh::V];
        if (dir==Direction::Parent ? i!=0 : k(s[i])==p.children[i]) return 0;
        const auto &e=dir==Direction::Parent ? p.parent : p.child[i];
        if (active(e.f,extend(e,s))>=0) return 0;
        if (fast) {
            I old=exact_point(sh.family,s);
            ++s[i].degree;
            return donor_floor(exact_point(sh.family,s)-old,degree,pooled);
        }
        Q old=exact_weight(sh.family,s);
        ++s[i].degree;
        Q d=exact_weight(sh.family,s)-old;
        return donor_floor(d,degree,pooled);
    }
    // Book one copy of every eligible triple/donor incidence on its unique
    // four-set. The four-set itself may be outside the quartet weight support.
    I triple_donors_on_four(const Shape &sh,const Core *in) const {
        validate_shape(sh);
        mh::require(sh.family>=mh::P,"four members required");
        const auto &p=plans[sh.family-mh::V];
        Wide sum=0;
        for (int j=0;j<p.ncut;++j) {
            auto [donor,contact]=p.cuts[j];
            const Core z=in[sh.role[donor]];
            mh::require(z.degree>=2,"donor needs an actual degree or a pool lower bound");
            Prune small=prune(sh.family,donor);
            Members triple{}; int at=-1;
            for (int i=0;i<3;++i) {
                triple[i]=in[sh.role[small.take[i]]];
                if (small.take[i]==contact) at=i;
            }
            if(at<0) throw std::logic_error("four-set pruning lost contact");
            Direction dir=mh::shapes[sh.family].parent[contact]==donor ? Direction::Parent : Direction::Child;
            sum+=triple_donor(Shape::canonical(small.f),triple.data(),at,dir,z.degree,z.core_degree==0);
        }
        mh::require(sum<=INT64_MAX,"four-set credit outside int64 range");
        return I(sum);
    }
    // Ready-to-book additional source for the chronological layer of THIS
    // three- or four-set. Book this once, or book its components separately.
    I additional(const Shape &sh,const Core *in) const {
        I v=row(sh,in);
        if (sh.family<=mh::C) return v;
        Wide total=Wide(v)+triple_donors_on_four(sh,in);
        mh::require(total<=INT64_MAX,"additional source outside int64 range");
        return I(total);
    }
    // Independent audits and individual probes can recover the full rational.
    // For an absent embedding this returns zero.
    Q exact_row(const Shape &sh,const Core *in) const {
        Members s{}; bool fast;
        if (!setup(sh,in,s,fast)) return 0;
        Prepared out; out.fractional_gains=std::make_unique<std::array<Q,4>>();
        Q a=0; drift<false>(sh,s,out,a); return a;
    }
};
} // namespace mh_occ
