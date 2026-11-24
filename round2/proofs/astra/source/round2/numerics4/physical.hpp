// Full unrelaxed physical observable on an actual rooted core tree.
#pragma once
#include "../motif/cert4v_format.h"
#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstring>
#include <fcntl.h>
#include <limits>
#include <stdexcept>
#include <sys/mman.h>
#include <sys/stat.h>
#include <unistd.h>

namespace numeric4 {
using Real = long double;
using U = uint32_t;
using Count = uint64_t;
constexpr U none = UINT32_MAX;
static_assert(std::numeric_limits<Real>::digits >= 64, "requires extended precision");
static_assert(sizeof(Cert4vHeader)==248 && __BYTE_ORDER__==__ORDER_LITTLE_ENDIAN__);
inline void check(bool ok, const char *why) { if (!ok) throw std::runtime_error(why); }

struct Pair { Real lo=0, hi=0; };
inline Pair operator+(Pair a, Pair b) { return {a.lo+b.lo,a.hi+b.hi}; }
inline Pair operator-(Pair a, Pair b) { return {a.lo-b.lo,a.hi-b.hi}; }
inline Pair operator*(Real a, Pair b) { return {a*b.lo,a*b.hi}; }
struct Sum {
    Pair value{}, err{};
    void add(Pair x) {
        Pair y=x-err, t=value+y;
        err=(t-value)-y;
        value=t;
    }
};

enum Family { Vertex, Edge, Vertical, Cherry, Path, Fork, Stem, Trident, NF };
constexpr const char *names[NF]={"f","E","V","C","P","F","Y","T"};
constexpr int npos[NF]={1,2,3,3,4,4,4,4};
constexpr int nend[NF]={0,2,2,2,2,2,3,3};
constexpr int ends[NF][3]={{0},{0,1},{0,2},{1,2},{0,3},{1,3},{0,2,3},{1,2,3}};

struct File {
    Cert4vHeader h{};
    void *mem=MAP_FAILED;
    size_t bytes=0;
    const int64_t *w[NF]{};
    explicit File(const char *path) {
        int fd=open(path,O_RDONLY);
        check(fd>=0,"cannot open four-vertex snapshot");
        struct stat st{};
        bool ok=fstat(fd,&st)==0 && st.st_size>=(off_t)sizeof h;
        if (ok) { bytes=st.st_size; mem=mmap(nullptr,bytes,PROT_READ,MAP_PRIVATE,fd,0); }
        close(fd);
        check(ok && mem!=MAP_FAILED,"cannot mmap four-vertex snapshot");
        memcpy(&h,mem,sizeof h);
        check(c4_validate_header(&h)==0 && h.DF<65536 && h.KHU<=128,"invalid snapshot");
    }
    File(const File &)=delete;
    ~File(){if(mem!=MAP_FAILED) munmap(mem,bytes);}
    void bind(const Cert4vSizes &z) {
        check(bytes==z.total_bytes,"wrong snapshot length");
        auto *p=reinterpret_cast<const int64_t *>(static_cast<const char *>(mem)+sizeof h)
                +z.n_cent+z.n_fy;
        const size_t lengths[NF]={z.n_f,z.n_g,z.n_hV,z.n_hC,z.n_hP4,z.n_hF4,z.n_hY4,z.n_hT4};
        for(int f=0;f<NF;++f) { w[f]=p; p+=lengths[f]; }
    }
};

struct State {
    Count l=0; U k=0; // table state: h=1, including the founder extension
    Count d() const { return l+k+1; }
    State gain() const { return {l+1,k}; }
    State bud() const { return {l-1,k+1}; }
};
struct Node { State s; Real r,b; };
using Nodes=std::array<Node,4>;
struct Mix { int a=-1,b=-1; Real t=0; }; // value at a + t(value at b - value at a)
using Mixes=std::array<Mix,3>;

struct Tables {
    File low,high;
    Cert4vLayout c{};
    explicit Tables(const char *l,const char *u):low(l),high(u) {
        check(low.h.mode==2 && high.h.mode==1 && low.h.S==high.h.S &&
              memcmp(&low.h.DF,&high.h.DF,sizeof(Cert4vHeader)-offsetof(Cert4vHeader,DF))==0,
              "incompatible nominal snapshots");
        c4_init_layout(&c,&low.h);
        Cert4vSizes z{}; c4_compute_sizes(&c,&z);
        low.bind(z); high.bind(z);
    }
    ~Tables(){c4_free_layout(&c);}
    const Cert4vLevel &level(Family f) const { return f==Edge?c.L2:f<=Cherry?c.L3:c.L4; }
    Pair corner(Family f,size_t off) const {
        return {Real(low.w[f][off]),Real(high.w[f][off])};
    }
    int active(Family f,const Nodes &a) const {
        auto inner=[&](int j){return a[j].s.d()<U(c.H.DHU);};
        auto l=[&](int j){return int(a[j].s.l);};
        auto k=[&](int j){return int(a[j].s.k);};
        switch(f) {
        case Edge: return 0;
        case Vertical: { if(!inner(1)) return -1; int u=c4_hu_idx(&c,l(1),k(1)); return u<0?-1:c.hu_to_act[u]; }
        case Cherry: { if(!inner(0)) return -1; int u=c4_cu_idx(&c,l(0),k(0)); return u<0?-1:c.cu_to_act[u]; }
        case Path: return inner(1)&&inner(2)?c4_p4_act(&c,l(1),k(1),l(2),k(2)):-1;
        case Fork: return inner(0)&&inner(2)?c4_f4_act(&c,l(0),k(0),l(2),k(2)):-1;
        case Stem: return inner(1)?c4_y4_act(&c,l(1),k(1)):-1;
        case Trident: return inner(0)?c4_t4_act(&c,l(0),k(0)):-1;
        default: throw std::logic_error("invalid interaction");
        }
    }
    static int tail(const Cert4vLevel &g,U k) { return g.TAIL_BASE+int(std::min(k,U(g.KT+1))); }
    static Mix mass(const Cert4vLevel &g,State s) {
        auto d=s.d();
        if(d<U(g.D0)) return {c4_lvl_inner(&g,int(s.l),int(s.k)),-1,0};
        int v=tail(g,s.k);
        if(d==U(g.D0)) return {v,-1,0};
        return {g.INF_ID,v,g.D0/Real(d)};
    }
    Pair at(Family f,int act,const std::array<int,3> &x) const {
        size_t off;
        switch(f) {
        case Edge: off=size_t(x[0])*c.L2.NG+x[1]; break;
        case Vertical: off=(size_t(c.L3.p_to_pi[x[0]])*c.NHU_ACT+act)*c.L3.NG+x[1]; break;
        case Cherry: off=size_t(act)*c.L3.pair_n+c4_lvl_pair(&c.L3,x[0],x[1]); break;
        case Path: off=c4_p4_off(&c,act,x[0],x[1]); break;
        case Fork: off=c4_f4_off(&c,act,x[0],x[1]); break;
        case Stem: off=c4_y4_off(&c,act,x[0],x[1],x[2]); break;
        case Trident: off=c4_t4_off(&c,act,x[0],x[1],x[2]); break;
        default: throw std::logic_error("invalid corner");
        }
        return corner(f,off);
    }
    template<int J>
    Pair interpol(Family f,int act,const Mixes &m,std::array<int,3> x) const {
        if constexpr(J==3) return at(f,act,x);
        else {
            if(J==nend[f]) return at(f,act,x);
            x[J]=m[J].a; Pair v=interpol<J+1>(f,act,m,x);
            if(m[J].b<0) return v;
            x[J]=m[J].b;
            return v+m[J].t*(interpol<J+1>(f,act,m,x)-v);
        }
    }
    Pair value(Family f,int a,const Mixes &m) const {
        return a<0?Pair{}:interpol<0>(f,a,m,{});
    }
    Mixes masses(Family f,const Nodes &a) const {
        Mixes m{}; const auto &g=level(f);
        for(int e=0;e<nend[f];++e) m[e]=mass(g,a[ends[f][e]].s);
        return m;
    }
    Pair value(Family f,const Nodes &a) const {
        int dom=active(f,a); return dom<0?Pair{}:value(f,dom,masses(f,a));
    }
    Pair singleton(State s) const {
        if(s.d()<U(c.H.DF)) return corner(Vertex,c4_vidx(&c,int(s.l),int(s.k)));
        Real x=1/Real(s.d()),q=s.l*x;
        return {q*(low.h.f_p0-low.h.f_p1*q-low.h.f_p2*x),high.h.f_p0+high.h.f_p1*q};
    }
    Pair singleton_drift(Node v,Pair base) const {
        State s=v.s;
        if(s.d()<U(c.H.DF)) {
            Pair dg=v.r*(singleton(s.gain())-base);
            return v.b>0?dg+v.b*(singleton(s.bud())-base):dg;
        }
        Real d=s.d(),x=1/d,xp=1/(d+1),q=s.l*x, dq=(Real(s.k)+1)*x*xp;
        Real a=low.h.f_p0,b=low.h.f_p1,h=low.h.f_p2, slope=high.h.f_p1;
        Pair plus{dq*(a-b*(2*q+dq)-h*xp)+h*q*x*xp,slope*dq};
        Pair bud{-x*(a-b*(2*q-x)-h*x),-slope*x};
        return v.r*plus+v.b*bud;
    }
};

struct Component { Sum phi,drift; };
struct Observation { Component parts[2][NF]{}; }; // bulk, terms containing actual founder
struct Tree {
    size_t count;
    const Count *l;
    const U *c,*par,*first,*next;
    const double *neighbor;
    U only=none; // optional restriction of generator to changes at one site in the audit
    unsigned events=3; // gain=1, bud=2
    Node node(U i) const {
        Real d=Real(l[i]+c[i]);
        bool live=(only==none || only==i);
        return {{l[i],c[i]-U(i!=0)},live&&(events&1)?Real(l[i])+neighbor[i]:0,live&&(events&2)?l[i]/d:0};
    }
};

class Observer {
    const Tables &t;
    const Tree &g;
    Observation out{};
    const Node baby{{1,0},0,0};
    void existing(Family f,const Nodes &ns,bool root) {
        int a=t.active(f,ns);
        if(a<0) return; // degree and child count only increase for this fixed set
        const auto &lev=t.level(f);
        Mixes m=t.masses(f,ns);
        Pair base=t.value(f,a,m),drift{};
        out.parts[root][f].phi.add(base);
        for(int j=0;j<npos[f];++j) {
            Node v=ns[j]; if(v.r==0 && v.b==0) continue;
            int e=0; while(e<nend[f] && ends[f][e]!=j) ++e;
            if(e<nend[f]) {
                auto eval=[&](Mix changed){ Mixes ms=m; ms[e]=changed; return t.value(f,a,ms); };
                if(v.s.d()>=U(lev.D0)) {
                    Real d=v.s.d(),w=lev.D0/d;
                    int j0=Tables::tail(lev,v.s.k), jb=Tables::tail(lev,v.s.k+1);
                    Pair knot=eval({j0,-1,0});
                    drift=drift+(w*(v.r/(d+1)))*(eval({lev.INF_ID,-1,0})-knot);
                    if(v.b && jb!=j0) drift=drift+(w*v.b)*(eval({jb,-1,0})-knot);
                } else {
                    drift=drift+v.r*(eval(Tables::mass(lev,v.s.gain()))-base);
                    if(v.b) drift=drift+v.b*(eval(Tables::mass(lev,v.s.bud()))-base);
                }
            } else {
                Nodes changed=ns; changed[j].s=v.s.gain();
                drift=drift+v.r*(t.value(f,t.active(f,changed),m)-base);
                if(v.b) {
                    changed[j].s=v.s.bud();
                    drift=drift+v.b*(t.value(f,t.active(f,changed),m)-base);
                }
            }
        }
        out.parts[root][f].drift.add(drift);
    }
    // Removing the newborn from its post-bud embedding gives exactly one old set.
    void born(Family f,Nodes ns,int promoter,bool root) {
        Real b=ns[promoter].b;
        if(!b) return;
        ns[promoter].s=ns[promoter].s.bud();
        out.parts[root][f].drift.add(b*t.value(f,ns));
    }
public:
    Observer(const Tables &tables,const Tree &tree):t(tables),g(tree){}
    Observation run() {
        Pair fb=t.singleton(baby.s);
        for(U w=0;w<g.count;++w) {
            Node v=g.node(w);
            Pair base=t.singleton(v.s);
            out.parts[w==0][Vertex].phi.add(base);
            out.parts[w==0][Vertex].drift.add(t.singleton_drift(v,base));
            out.parts[0][Vertex].drift.add(v.b*fb); // every promoter, including founder
            born(Edge,{v,baby},0,w==0);
            U p=g.par[w];
            if(p!=none) {
                Node a=g.node(p); bool root=(p==0);
                existing(Edge,{a,v},root);
                born(Vertical,{a,v,baby},1,root);
                born(Cherry,{a,v,baby},0,root);
                U gp=g.par[p];
                if(gp!=none && a.s.d()<U(t.c.H.DHU) && a.s.k<=U(t.c.H.KHU)) {
                    Node b=g.node(gp); bool r=(gp==0);
                    existing(Vertical,{b,a,v},r);
                    born(Path,{b,a,v,baby},2,r);
                    born(Stem,{b,a,v,baby},1,r);
                    born(Fork,{b,baby,a,v},0,r);
                    U ggp=g.par[gp];
                    if(ggp!=none) existing(Path,{g.node(ggp),b,a,v},ggp==0);
                    if(b.s.d()<U(t.c.H.DF4_1) && b.s.k<=U(t.c.H.KF4_1))
                        for(U side=g.first[gp];side!=none;side=g.next[side])
                            if(side!=p) existing(Fork,{b,g.node(side),a,v},r);
                }
            }
            if(v.s.d()>=U(t.c.H.DCU) || v.s.k>U(t.c.H.KCU)) continue;
            for(U x=g.first[w];x!=none;x=g.next[x]) {
                Node a=g.node(x);
                for(U y=g.next[x];y!=none;y=g.next[y]) {
                    Node b=g.node(y); bool root=(w==0);
                    existing(Cherry,{v,a,b},root);
                    born(Trident,{v,a,b,baby},0,root);
                    born(Fork,{v,b,a,baby},2,root);
                    born(Fork,{v,a,b,baby},2,root);
                    if(p!=none) existing(Stem,{g.node(p),v,a,b},p==0);
                    if(v.s.d()<U(t.c.H.DT4) && v.s.k<=U(t.c.H.KT4))
                        for(U z=g.next[y];z!=none;z=g.next[z])
                            existing(Trident,{v,a,b,g.node(z)},root);
                }
            }
        }
        return out;
    }
};
} // namespace numeric4
