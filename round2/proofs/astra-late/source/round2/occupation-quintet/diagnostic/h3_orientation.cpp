#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main census_main
#include "census.cpp"
#undef main
#pragma GCC diagnostic pop

// One row per physical H3, with separate actual arrows (no clipping per arrow).
struct Measure {I n=0;R ch=0,neg=0,gap=0;};
struct Orientation {
    const Witness& p;const Tree& t;
    map<array<int,8>,Measure> rows;
    map<array<int,8>,Measure> ports;
    map<array<int,4>,Measure> exposure;
    R sum=0, negative=0, gaps=0;
    I count=0;
    static int band(I d){return int(lower_bound(limits.begin(),limits.end(),d)-limits.begin());}
    static constexpr array<I,17> limits{3,5,8,11,18,24,37,65,100,179,360,1988,1998,1999,2000,2011,9999};
    Orientation(const Witness& w,const Tree& tr):p(w),t(tr){}
    void add(const Set& s,const Shape& q,int m,int l,int h){
        auto a=t.m[m],lo=t.m[l],hi=t.m[h];
        int o=t.par[m]==h?0:t.par[m]==l?1:2;
        bool founder=s.has(0);
        bool strict=o!=2&&!founder&&a.d<=5&&a.c==2&&hi.d<=65&&
            (o==0?hi.c==2&&lo.c==1:lo.c==2&&hi.c==1);
        bool highmid=a.c>=9&&a.c<=12&&a.d<=24&&hi.d<=6;
        int flags=founder+2*strict+4*highmid;
        bool tie=lo.d==hi.d;
        int stock=10000*hi.c+100*lo.c+int(min<I>(lo.d-lo.c,99));
        array<int,8> key{o,band(hi.d),band(lo.d),stock,int(a.d),a.c,flags,tie};
        bool small=!flags&&!tie&&o==0&&a.c<=4&&a.d<=5&&
            ((hi.d>=19&&hi.d<=65&&lo.d<=18)||(hi.d>=19&&hi.d<=100&&lo.d<=11));
        R old=p.weight(q,t), score=-p.sign*p.lam*old, allgap=0;
        for(int u:{m,l,h}){
            auto [d,c]=t.m[u];R dz=p.sign*(p.weight(q,t,u,1)-old);
            score+=(d-c)*dz;
            if(d>c){Set b=s;b.push(t.m.size());
                score+=p.sign*R(d-c)/d*(p.weight(q,t,u,2)-old+p.weight(classify(b,t,u),t,u,2));}
            int dst=u==m?0:u==l?1:2;
            auto arrow=[&](int v,int type,R value,R gap){
                auto reg=[&](int f){auto& z=ports[{o,band(hi.d),band(lo.d),dst,type,band(t.m[v].d),f,tie}];
                    ++z.n;z.ch+=value;z.gap+=gap;};
                reg(flags);if(small)reg(8);};
            for(int v:t.a[u]){
                R x=1.L/t.m[v].d;
                if(s.has(v)){
                    int src=v==m?0:v==l?1:2;
                    score+=x*dz;arrow(v,src,x*dz,0);
                }else{
                    Set b=s;b.push(v);
                    bool covered=p.supported(classify(b,t),t);
                    R val=(covered?p.center(t.m[u]):x)*dz;
                    R gap=covered?0:max(R(0),dz/2)-x*dz;
                    score+=val;allgap+=gap;
                    arrow(v,(covered?9:3)+(t.par[u]!=v),x*dz,gap);
                    if(covered)arrow(v,5+(t.par[u]!=v),val-x*dz,0);
                }
            }
            if(u==m){
                for(int v:{l,h}){
                    Set e;e.push(m);e.push(v==l?h:l);auto pre=classify(e,t);
                    R moved=p.sign*(1.L/t.m[v].d-p.center(a))*(p.weight(pre,t,m,1)-p.weight(pre,t));
                    score+=moved;arrow(v,v==l?7:8,moved,0);
                }
            }
        }
        R neg=max(R(0),-score);
        auto reg=[&](Measure& z){++z.n;z.ch+=score;z.neg+=neg;z.gap+=allgap;};
        reg(rows[key]);++count;sum+=score;negative+=neg;gaps+=allgap;
        if(small&&hi.d<=65&&lo.d<=18){
            int ancestor=t.par[h];assert(ancestor>=0);
            I parentD=t.m[ancestor].d;R parent=1.L/parentD, kids=0;
            for(int j:t.a[h])if(j!=m&&j!=ancestor)kids+=1.L/t.m[j].d;
            auto& z=exposure[{hi.c,band(parentD),int(kids*4),band(hi.d)}];
            ++z.n;z.ch+=parent;z.neg+=neg;z.gap+=kids;
        }
    }
    template<class M> void print(const char* tag,const M& v){for(auto [k,x]:v){
        cout<<tag;for(int j:k)cout<<'\t'<<j;
        cout<<'\t'<<x.n<<'\t'<<x.ch/p.h.S<<'\t'<<x.neg/p.h.S<<'\t'<<x.gap/p.h.S<<'\n';}}
    void run(const char* path){
        for(int m=0;m<(int)t.m.size();++m)if(Witness::index(p.z.M3,t.m[m])>=0){
            const auto& nb=t.a[m];
            for(size_t i=0;i<nb.size();++i)for(size_t j=i+1;j<nb.size();++j){
                int l=nb[i],h=nb[j];
                if(make_pair(t.m[l].d,l)>make_pair(t.m[h].d,h))swap(l,h);
                Set s;s.push(m);s.push(l);s.push(h);add(s,classify(s,t),m,l,h);
            }}
        cout<<"tree\t"<<path<<'\t'<<p.h.mode<<'\t'<<t.m.size()<<'\t'<<count<<'\t'<<sum/p.h.S<<'\t'<<negative/p.h.S<<'\t'<<gaps/p.h.S<<'\n';
        print("row",rows);print("port",ports);
        for(auto [k,z]:exposure){cout<<"rate";for(int i:k)cout<<'\t'<<i;
            cout<<'\t'<<z.n<<'\t'<<z.ch<<'\t'<<z.gap<<'\t'<<z.neg/p.h.S<<'\n';}
    }
};
int main(int argc,char** argv){assert(argc>=3);Witness w(argv[1]);cout<<setprecision(18);
    for(int i=2;i<argc;++i){Tree t(argv[i]);Orientation(w,t).run(argv[i]);}return 0;}
