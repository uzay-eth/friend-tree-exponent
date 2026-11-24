#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main census_main
#include "../census.cpp"
#undef main
#pragma GCC diagnostic pop

struct Sample { I n=0; R direct=0, box=0, increment=0, one=0, up=0; };

void sample(const Witness& w,const Tree& t,const char* name){
    array<Sample,2> sums{};
    for(int m=1;m<(int)t.m.size();++m){
        int h=t.par[m];
        if(h<=0)continue;
        auto H=t.m[h], M=t.m[m];
        if(H.d<19||H.d>65||H.c<3||H.c>7||M.d>5||M.c<2||M.c>4)continue;
        for(int y:t.a[m]){
            if(t.par[y]!=m||t.m[y].d>18||t.m[y].c>3)continue;
            Set s;s.push(h);s.push(m);s.push(y);
            const auto q=classify(s,t);assert(q.k==::H);
            R old=w.weight(q,t),score=-w.sign*w.lam*old;
            vector<R> z;
            for(int u:{h,m,y}){
                auto [d,c]=t.m[u];
                R jump=w.sign*(w.weight(q,t,u,1)-old);
                score+=(d-c)*jump;
                if(d>c){Set e=s;e.push(t.m.size());
                    score+=w.sign*R(d-c)/d*(w.weight(q,t,u,2)-old+w.weight(classify(e,t,u),t,u,2));}
                for(int v:t.a[u]){
                    R x=1.L/t.m[v].d;
                    if(s.has(v)){score+=x*jump;continue;}
                    Set e=s;e.push(v);
                    bool covered=w.supported(classify(e,t),t);
                    score+=(covered?w.center(t.m[u]):x)*jump;
                    if(u==h&&t.par[v]==h)z.push_back(covered?0:max(R(0),jump/2)-x*jump);
                }
            }
            for(int v:{h,y}){
                Set e;e.push(m);e.push(v==h?y:h);auto q2=classify(e,t);
                score+=w.sign*(1.L/t.m[v].d-w.center(M))*(w.weight(q2,t,m,1)-w.weight(q2,t));
            }
            assert(z.size()==size_t(H.c-2));
            R b=score;for(R x:z)b+=x;
            auto clip=[](R x){return max(R(0),-x);};
            R base=clip(b),delta=clip(score)-base,one=0,up=0;
            for(R x:z){one+=clip(b-x)-base;up+=(clip(b-z.size()*x)-base)/z.size();}
            assert(one<=delta+1e-5L&&delta<=up+1e-5L);
            auto& a=sums[H.d>=38];
            ++a.n;a.direct+=clip(score);a.box+=base;a.increment+=delta;a.one+=one;a.up+=up;
        }
    }
    for(size_t i=0;i<sums.size();++i){auto a=sums[i];
        cout<<name<<'\t'<<i<<'\t'<<t.m.size()<<'\t'<<a.n;
        for(R x:{a.direct,a.box,a.increment,a.one,a.up})cout<<'\t'<<x/w.h.S;
        cout<<'\n';}
}
int main(int argc,char** argv){assert(argc>=3);Witness w(argv[1]);cout<<setprecision(18);
    for(int i=2;i<argc;++i){Tree t(argv[i]);sample(w,t,argv[i]);}}
