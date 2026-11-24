#include "quintet_format.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
using R=long double;
using I=int64_t;
using namespace std;
enum {F,G,H,P,S,Q,B,T};
const char *names[]={"f","g","h3","p4","s3","p5","b5","s4"};
struct Mark {I d; int c;};
struct Set {
    array<int,6> v{};
    int n=0;
    void push(int i){assert(n<6);v.at(n++)=i;}
    bool has(int i)const{return find(v.begin(),v.begin()+n,i)!=v.begin()+n;}
};
struct Shape {int k=-1; Set in,out; size_t block=0;};
struct Tree {
    vector<int> par;
    vector<vector<int>> a;
    vector<Mark> m;
    vector<R> gain;
    I pop=0;
    Tree(const char* path){
        ifstream f(path); assert(f); int p; I l;
        while(f>>p>>l){int n=int(m.size());assert(l>=0 && (n ? 0<=p && p<n : p==-1));
            par.push_back(p);a.emplace_back();m.push_back({l,0});pop+=l+1;
            if(n){a[p].push_back(n);a[n].push_back(p);}}
        assert(f.eof() && !m.empty());
        for(size_t i=0;i<m.size();++i){m[i].c=int(a[i].size());m[i].d+=m[i].c;assert(m[i].d>=2);}
        for(size_t i=0;i<m.size();++i){R g=m[i].d-m[i].c;for(int j:a[i])g+=1.L/m[j].d;gain.push_back(g);}
    }
    bool link(int u,int v)const{return (v==(int)par.size() || u==(int)par.size()) ? false : par[u]==v || par[v]==u;}
};
struct Witness {
    int fd; const void* map; size_t len; const int64_t* w;
    QuintetHeader h; Q5Layout z; R sign, lam;
    Witness(const char* f){fd=open(f,O_RDONLY);assert(fd>=0);struct stat st{};assert(fstat(fd,&st)==0);
         len=st.st_size;assert(len==6344902168);map=mmap(nullptr,len,PROT_READ,MAP_PRIVATE,fd,0);assert(map!=MAP_FAILED);
         memcpy(&h,map,512); assert(h.S==(1LL<<40) && h.Q==100000000 && h.DF==4000);
         assert(h.P==(h.mode==2?56571368:56573686));
         constexpr array<I,20> expected{0,40,80,120,168,252,360,504,672,840,1120,1440,1680,2016,2520,2880,3360,3780,4200,5040};
         assert(h.nknots==20&&h.Q_K==10080&&equal(expected.begin(),expected.end(),h.knot_num));
         z.init(h);assert(z.total_bytes==len);
        w=(const int64_t*)map+64; sign=h.mode==2?-1:1;lam=R(h.P)/h.Q;
    }
    ~Witness(){munmap(const_cast<void*>(map),len);close(fd);}
    int exact(Mark a)const{return a.d<z.G2.dcut && a.c>=1 && a.c<=z.G2.Kcut ? z.G2.idx[a.d][a.c] : -1;}
    R center(Mark a)const{int i=exact(a);assert(i>=0);return R(w[z.off_r1+i])/h.S_r;}
    R knot(int i)const{return R(h.knot_num[i])/h.Q_K;}
    I fyword(int i,int j)const{return w[z.off_knots+i*h.nknots+j];}
    R fy(int i,R x)const{
        int j=1;while(j+1<h.nknots && knot(j)<x)++j;
        R t=(x-knot(j-1))/(knot(j)-knot(j-1));return (1-t)*fyword(i,j-1)+t*fyword(i,j);
    }
    static int index(const Q5InnerSpace& s,Mark a){return a.d<s.dcut && a.c<=s.maxC && a.c>=s.minC?s.idx[a.d][a.c]:-1;}
    using Corn=array<pair<int,R>,2>;
    static Corn corner(const Q5EndSpace& s,Mark a){
        if(a.d<s.dcut)return {{{a.c<=s.Kcut?s.idx[a.d][a.c]:s.ovf_idx[a.d],1},{0,0}}};
        R t=R(s.dcut)/a.d;return {{{s.t_id(a.c),t},{s.INF,1-t}}};
    }
    R weight(const Shape& q,const Tree& tr,int who=-1,int how=0)const{
        auto mark=[&](int v){if(v==(int)tr.m.size())return Mark{2,1};Mark a=tr.m[v];
            if(v==who){a.d+=how==1;a.c+=how==2;}return a;};
        if(q.k==F){auto [d,c]=mark(q.in.v[0]);if(d<h.DF)return w[z.off_f+z.f_id(d,c)];
            R x=1.L/d, b=R(d-c)/d;return h.mode==1?h.f_p0+h.f_p1*b:b*(h.f_p0-h.f_p1*b-h.f_p2*x);}
        array<int,3> j{};
        auto id=[&](int i,const Q5InnerSpace& sp){return j[i]=index(sp,mark(q.in.v[i]));};
        switch(q.k){
            case H:if(id(0,z.M3)<0)return 0;break;
            case P:if(id(0,z.M4v)<0||id(1,z.M4v)<0||z.p4_ord_map[j[0]*z.M4v.N+j[1]]<0)return 0;break;
            case S:if(id(0,z.MS3)<0)return 0;break;
            case Q:{int n=z.I5_2.N;if(id(0,z.I5_2)<0||id(1,z.I5_2)<0||id(2,z.I5_2)<0||z.p5_ord_map[(j[0]*n+j[1])*n+j[2]]<0)return 0;break;}
            case B:if(id(0,z.I5_3)<0||id(1,z.I5_2c)<0||(j[0]=z.b5_int_map[j[0]*z.I5_2c.N+j[1]])<0)return 0;break;
            case T:if(id(0,z.MS4)<0)return 0;break;
        }
        const auto& sp=q.k==G?z.G2:q.k==H?z.G3:q.k<=S?z.G4:z.G5;
        array<Corn,4> cc{};
        for(int i=0;i<q.out.n;++i)cc[i]=corner(sp,mark(q.out.v[i]));
        R ans=0;
        for(int b=0;b<(1<<q.out.n);++b){array<int,4>a{};R mass=1;
            for(int i=0;i<q.out.n;++i){auto [e,t]=cc[i][(b>>i)&1];a[i]=e;mass*=t;}if(mass==0)continue;
            size_t ix=0;
            switch(q.k){
                case G:ix=z.off_g+z.g_off(a[0],a[1]);break;
                case H:ix=z.off_h3+z.h3_off(j[0],a[0],a[1]);break;
                case P:ix=z.off_p4+z.p4_off(j[0],j[1],a[0],a[1]);break;
                case S:ix=z.off_s3+z.s3_off(j[0],a[0],a[1],a[2]);break;
                case Q:ix=z.off_p5+z.p5_off(j[0],j[1],j[2],a[0],a[1]);break;
                case B:ix=z.off_b5+z.b5_off(j[0],a[0],a[1],a[2]);break;
                case T:ix=z.off_s4+z.s4_off(j[0],a[0],a[1],a[2],a[3]);break;
            }
            ans+=mass*w[ix];
        }return ans;
    }
    bool supported(const Shape& q,const Tree& tr)const{
        auto x=[&](int i){return tr.m[q.in.v[i]];};
        switch(q.k){
            case H:return index(z.M3,x(0))>=0;
            case P:return index(z.M4v,x(0))>=0&&index(z.M4v,x(1))>=0&&x(0).d+x(1).d<=h.SP4;
            case S:return index(z.MS3,x(0))>=0;
            case Q:return index(z.I5_2,x(0))>=0&&index(z.I5_2,x(1))>=0&&index(z.I5_2,x(2))>=0&&
                max(x(0).d+x(1).d,x(1).d+x(2).d)<=h.SP4&&x(0).d+x(1).d+x(2).d<=h.SP5;
            case B:return index(z.I5_3,x(0))>=0&&index(z.I5_2c,x(1))>=0&&x(0).d+x(1).d<=h.SB5;
            case T:return index(z.MS4,x(0))>=0;
            default:return true;
        }
    }
};

Shape classify(const Set& s,const Tree& tr,int nbhost=-1){
    Shape q;
    if(s.n==1){q.k=F;q.in=s;return q;}
    if(s.n==2){q.k=G;q.out=s;return q;}
    int nb=tr.m.size();
    auto link=[&](int u,int v){return tr.link(u,v)||(nbhost>=0&&((u==nb&&v==nbhost)||(v==nb&&u==nbhost)));};
    array<int,6> ar{};int hub=-1;
    for(int i=0;i<s.n;++i){for(int j=0;j<s.n;++j)ar[i]+=link(s.v[i],s.v[j]);
        assert(ar[i]>0);if(ar[i]==1)q.out.push(s.v[i]);if(ar[i]>=3)hub=s.v[i];}
    if(hub<0){q.k=s.n==3?H:s.n==4?P:Q;int u=q.out.v[0],prev=-1;
        for(int i=1;i<s.n-1;++i){int v=-1;for(int j=0;j<s.n;++j)if(s.v[j]!=prev&&link(u,s.v[j])){v=s.v[j];break;}
            assert(v>=0);q.in.push(v);prev=u;u=v;}return q;}
    q.in.push(hub);
    if(q.out.n==s.n-1){q.k=s.n==4?S:T;return q;}
    q.k=B;int arm=-1;
    for(int i=0;i<s.n;++i)if(ar[i]==2)arm=s.v[i];
    assert(arm>=0);q.in.push(arm);
    for(int i=0;i<2;++i)if(link(arm,q.out.v[i]))swap(q.out.v[i],q.out.v[2]);
    return q;
}

struct Agg {I n=0;R amount=0;};
struct Ledger {
    const Witness& p;const Tree& t;
    array<I,8> count{};
    array<R,8> mass{},raw{},charged{},negative{},positive{},slots{},lowcap{};
    array<array<R,8>,8> transfer{};
    vector<R> dg,db;
    map<array<int,8>,Agg> groups;
    Ledger(const Witness& pp,const Tree& tt):p(pp),t(tt),dg(t.m.size()),db(t.m.size()){}
    static int degree(I d){return d<=5?0:d<=11?1:d<=18?2:d<=37?3:d<=179?4:d<=1998?5:d==1999?6:d<=2011?7:d<=9999?8:9;}
    static int arity(int c){return c<=1?c:c<=4?2:c<=8?3:c<=15?4:c<=26?5:6;}
    void log(array<int,8> key,R amount){auto& a=groups[key];++a.n;a.amount+=amount;}
    void row(const Set& s,const Shape& q){
        const int k=q.k; ++count[k]; R old=p.weight(q,t), r=-p.sign*p.lam*old, moves=0, gap=0;
        mass[k]+=old;
        I maxd=0;int maxc=0;for(int i=0;i<s.n;++i){auto a=t.m[s.v[i]];maxd=max(maxd,a.d);maxc=max(maxc,a.c);}
        for(int i=0;i<s.n;++i){
            int u=s.v[i];auto a=t.m[u]; R dz=p.sign*(p.weight(q,t,u,1)-old), b=R(a.d-a.c)/a.d;
            r+=t.gain[u]*dz;dg[u]+=dz;
            if(b){R jump=p.weight(q,t,u,2)-old;
                if(s.n<5){Set add=s;add.push(t.m.size());jump+=p.weight(classify(add,t,u),t,u,2);}
                if(k==F){Set j;j.push(t.m.size());jump+=p.weight(classify(j,t),t);}
                r+=b*p.sign*jump;db[u]+=p.sign*jump;
            }
            if(k==F){if(p.exact(a)>=0){R sent=(t.gain[u]-(a.d-a.c)-a.c*p.center(a))*dz;
                    moves-=sent;transfer[F][G]+=sent;}continue;}
            int ex=p.exact(a);
            if(k==G&&ex>=0){R inc=p.sign*(1.L/t.m[s.v[1-i]].d-p.center(a))*(p.weight(Shape{F,Set{{u},1},{},0},t,u,1)-p.weight(Shape{F,Set{{u},1},{},0},t));moves+=inc;}
            bool fy=k==G&&ex>=0&&a.c>=2&&Witness::index(p.z.M3,a)<0;
            R hi=0, fyto=0;
            if(fy){hi=-INFINITY;for(int j=0;j<p.h.nknots;++j)hi=max(hi,p.knot(j)*dz-p.sign*p.fyword(ex,j));
                fyto=p.sign*p.fy(ex,1.L/t.m[s.v[1-i]].d);moves+=(a.c-1)*fyto;}
            for(int v:t.a[u])if(!s.has(v)){
                R x=1.L/t.m[v].d;
                if(fy){R from=p.sign*p.fy(ex,x);moves-=from;R j=hi-(x*dz-from);gap+=j;
                    int pay=s.v[1-i];auto d=t.m[pay],e=t.m[v];
                    int dir=t.par[u]==pay?0:t.par[u]==v?1:(pay<v?2:3);
                    bool low=max(a.d,d.d)<=18;
                    bool main=!low&&a.d<=175&&d.d<=5&&e.d<=5&&a.c<=8;
                    bool ext=!low&&a.d<=175&&d.d>5&&d.d<=36&&e.d<=5&&a.c<=10&&d.c<=6;
                    bool older=!low&&a.d<=90&&a.c<=5&&d.d<=5&&e.d>5&&e.d<=24&&e.c<=5;
                    log({1,dir,degree(a.d),degree(d.d),degree(e.d),arity(a.c),low?1:main?2:ext?3:older?4:0,arity(d.c)},j);
                }else{
                    bool covered=false;
                    if(s.n<5){Set add=s;add.push(v);Shape dst=classify(add,t);
                        if(p.supported(dst,t)){R sent=(x-p.center(a))*dz;moves-=sent;transfer[k][dst.k]+=sent;covered=true;}}
                    if(!covered){R j=max(R(0),dz/2)-x*dz;gap+=j;
                        log({2,k,t.par[u]==v?0:1,degree(a.d),degree(t.m[v].d),arity(a.c),degree(maxd),arity(maxc)},j);}
                }
            }
            // Each deleted eligible tip sends an actual lower-level gain term.
            if(k>G){for(int j=0;j<q.out.n;++j){int v=q.out.v[j];if(!t.link(u,v))continue;
                Set src;for(int l=0;l<s.n;++l)if(s.v[l]!=v)src.push(s.v[l]);
                Shape pre=classify(src,t);assert(p.supported(pre,t));
                moves+=p.sign*(1.L/t.m[v].d-p.center(a))*(p.weight(pre,t,u,1)-p.weight(pre,t));}}
        }
        raw[k]+=r;R ch=r+moves;charged[k]+=ch;negative[k]+=max(R(0),-ch);positive[k]+=max(R(0),ch);slots[k]+=gap;
        log({0,k,degree(maxd),arity(maxc),q.in.n?degree(t.m[q.in.v[0]].d):-1,
             q.in.n?arity(t.m[q.in.v[0]].c):-1,s.has(0),0},ch);
        log({3,k,degree(maxd),arity(maxc),0,0,s.has(0),0},max(R(0),-ch));
         int mask=0;
         if(k==H||k==P){
            int apex=*min_element(s.v.begin(),s.v.begin()+s.n);
            bool strict=apex!=0&&(k==H?q.in.v[0]!=apex&&maxc<=2:
                apex!=q.in.v[0]&&apex!=q.in.v[1]&&maxc<=2);
             if(strict){int tip=q.out.v[0]==apex?q.out.v[1]:q.out.v[0];
                 if(t.m[tip].c==1&&((k==H&&t.m[q.in.v[0]].d<=5&&maxd<=65)||
                     (k==P&&t.m[q.in.v[0]].d<=5&&t.m[q.in.v[1]].d<=5&&maxd<=36)))mask=1;}
            // Also exclude the separately tested high-middle donor lane.
            auto m=t.m[q.in.v[0]];
             if(k==H&&m.c>=9&&m.c<=12&&m.d<=24&&t.m[q.out.v[0]].d<=6&&t.m[q.out.v[1]].d<=6)mask=2;
            I d0=t.m[q.out.v[0]].d,d1=t.m[q.out.v[1]].d;
             log({5,k,degree(min(d0,d1)),degree(max(d0,d1)),degree(m.d),arity(m.c),mask,s.has(0)},max(R(0),-ch));
        }
        if(k>=Q){I tipd=0;for(int j=0;j<q.out.n;++j)tipd=max(tipd,t.m[q.out.v[j]].d);
            log({6,k,degree(tipd),arity(maxc),degree(t.m[q.in.v[0]].d),arity(t.m[q.in.v[0]].c),s.has(0),0},max(R(0),-ch));}
         if(k==G){I lo=min(t.m[s.v[0]].d,t.m[s.v[1]].d);R cap=p.h.mode==2?792:463;
             log({7,k,degree(lo),degree(maxd),arity(maxc),s.has(0),0,0},cap-ch-gap);}
         if(k==F){auto x=t.m[s.v[0]];R c=p.h.mode==2?55:2359;
             log({4,k,degree(x.d),arity(x.c),0,0,s.has(0),0},c-ch);}
         if(k==G&&maxd<=18)lowcap[k]+=(p.h.mode==2?792:463)-ch;
         if(k==G){
             int y=max(s.v[0],s.v[1]),h=t.par[y];assert(s.has(h));
             auto a=t.m[h],b=t.m[y];
             log({8,k,degree(a.d),degree(b.d),arity(a.c),arity(b.c),s.has(0),0},max(R(0),-ch));
         }
         if(k==H){
             int m=q.in.v[0], x=q.out.v[0],y=q.out.v[1];
             if(t.m[x].d<t.m[y].d)swap(x,y);
             auto a=t.m[x],b=t.m[y],c=t.m[m];
             int orient=a.d==b.d?3:t.par[m]==x?0:t.par[m]==y?1:2;
             if(!s.has(0)&&!mask&&a.d>=19&&a.d<=1998&&c.c<=4)
                 log({8,k,orient,degree(a.d),degree(b.d),c.d<=5?5:c.d<=24?24:25,arity(c.c),0},max(R(0),-ch));
             if(!s.has(0)&&orient==0&&c.d<=5&&c.c<=4&&b.d<=18)
                 log({9,k,a.d<=37?37:a.d<=65?65:a.d<=100?100:a.d<=179?179:1998,
                       b.d<=5?5:b.d<=11?11:18,a.c<=7?7:8,b.c<=3?3:4,
                       a.c==2&&c.c==2&&b.c==1&&a.d<=65,int(a.d>=19&&a.d<=1998)},max(R(0),-ch));
         }
    }
    void run(){
        // Unique connected subsets: root at the minimum label; expand the
        // frontier in fixed order. A monotone unsupported set cannot extend
        // to a supported higher motif by downward closure.
        auto walk=[&](auto&& self,Set s,vector<int> todo)->void{
            Shape sh=classify(s,t); if(!p.supported(sh,t))return;
            row(s,sh); if(s.n==5)return;
            for(size_t i=0;i<todo.size();++i){int u=todo[i];Set next=s;next.push(u);
                vector<int> children(todo.begin()+i+1,todo.end());
                if(Witness::index(p.z.M3,t.m[u])>=0){
                    for(int v:t.a[u])if(v>s.v[0]&&!s.has(v))children.push_back(v);
                }
                self(self,next,move(children));
            }
        };
        for(int i=0;i<(int)t.m.size();++i){Set s;s.push(i);vector<int> a;
            for(int j:t.a[i])if(j>i)a.push_back(j);
            // At a high-degree root only singletons and individual edges
            // can use more than one of its incident edges as candidates.
            if(Witness::index(p.z.M3,t.m[i])>=0)walk(walk,s,move(a));
            else{row(s,classify(s,t));for(int j:a){Set e=s;e.push(j);vector<int> f;
                if(Witness::index(p.z.M3,t.m[j])>=0)for(int v:t.a[j])if(v>i&&v!=i)f.push_back(v);
                walk(walk,e,move(f));}}
        }
    }
    void print(const char* path){
        R sr=0,sc=0,sd=0,m=0;
        for(int k=0;k<8;++k){sr+=raw[k];sc+=charged[k];m+=mass[k];}
        for(size_t i=0;i<t.m.size();++i){auto [d,c]=t.m[i];sd+=t.gain[i]*dg[i]+R(d-c)/d*db[i];}
        sd-=p.sign*p.lam*m;R scale=p.h.S;
        assert(abs(sr-sc)/scale<2e-9L&&abs(sr-sd)/scale<2e-9L);
        cout<<"tree\t"<<path<<'\t'<<p.h.mode<<'\t'<<t.m.size()<<'\t'<<t.pop<<'\t'<<m/scale<<'\t'<<sc/scale<<'\t'<<sr/scale<<'\t'<<sd/scale<<'\n';
        for(int k=0;k<8;++k){cout<<"family\t"<<names[k]<<'\t'<<count[k]<<'\t'<<mass[k]/scale<<'\t'<<raw[k]/scale<<'\t'<<charged[k]/scale<<'\t'<<negative[k]/scale<<'\t'<<positive[k]/scale<<'\t'<<slots[k]/scale<<'\t'<<lowcap[k]/scale<<'\n';
            for(int j=k+1;j<8;++j)if(transfer[k][j])cout<<"move\t"<<k<<'\t'<<j<<'\t'<<transfer[k][j]/scale<<'\n';}
        for(auto [key,a]:groups){cout<<"bin";for(int i:key)cout<<'\t'<<i;cout<<'\t'<<a.n<<'\t'<<a.amount/scale<<'\n';}
    }
};

int main(int argc,char**argv){assert(argc>=3);Witness w(argv[1]);cout<<setprecision(18);
    for(int i=2;i<argc;++i){Tree t(argv[i]);Ledger l(w,t);l.run();l.print(argv[i]);}}
