#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-type"
#define main census_main
#include "../census.cpp"
#undef main
#pragma GCC diagnostic pop

struct Tot { I n=0; R physical=0, base=0, joint=0, sib=0, up=0, gap=0, compressed=0; };
void sample(const Witness& w,const Tree& t,const char* name){
    map<pair<int,int>,Tot> sums;
    for(int m=1;m<(int)t.m.size();++m){
        int h=t.par[m];if(h<=0)continue;
        auto H=t.m[h], M=t.m[m];
        if(H.d<19||H.d>1998||M.c<2||M.c>4||Witness::index(w.z.M3,M)<0)continue;
        for(int y:t.a[m]){
            auto Y=t.m[y];
            if(t.par[y]!=m||Y.d>1998||Y.d>=H.d)continue;
            if(H.d<=65&&Y.d<=65&&M.d<=5&&H.c==2&&M.c==2&&Y.c==1)continue;
            Set s;s.push(h);s.push(m);s.push(y);
            const auto q=classify(s,t);
            R old=w.weight(q,t),score=-w.sign*w.lam*old;
            R sib=0,up=0,other=0,slope=0;
            int k=0;
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
                    if(covered)continue;
                    R z=max(R(0),jump/2)-x*jump;
                    if(u!=h){other+=z;continue;}
                    slope=jump;
                    if(v==t.par[h])up+=z;
                    else {sib+=z;k+=t.m[v].d==2;}
                }
            }
            for(int v:{h,y}){
                Set e;e.push(m);e.push(v==h?y:h);auto q2=classify(e,t);
                score+=w.sign*(1.L/t.m[v].d-w.center(M))*(w.weight(q2,t,m,1)-w.weight(q2,t));
            }
            auto clip=[](R x){return max(R(0),-x);};
            R B=score+sib+up, direct=clip(score),base=clip(B),between=clip(score+sib);
            auto tally=[&](int cut,R comp){
                int band=H.d<=37?37:H.d<=65?65:H.d<=100?100:H.d<=179?179:H.d<=360?360:1998;
                auto& a=sums[{cut,band}];
                ++a.n;a.physical+=direct;a.base+=base;a.joint+=direct-base;a.sib+=direct-between;
                a.up+=between-base;a.gap+=sib+up;a.compressed+=comp;
            };
            tally(0,0);
            if(M.d<=5&&Y.d<=18&&Y.c<=3&&H.c<=8){
                tally(1,0);
                if(t.m[t.par[h]].d>=3){
                    R mn=slope>0?slope*(R(H.c-1)/2-(R(k)/2+R(H.c-2-k)/3+1.L/3)):-slope*R(k)/2;
                    tally(2,max(R(0),mn-max(R(0),B+other)));
                }
            }
        }
    }
    for(auto [key,a]:sums){
        cout<<name<<'\t'<<key.first<<'\t'<<key.second<<'\t'<<t.m.size()<<'\t'<<a.n;
        for(R x:{a.physical,a.base,a.joint,a.sib,a.up,a.gap,a.compressed})cout<<'\t'<<x/w.h.S;
        cout<<'\n';}
}
int main(int argc,char** argv){assert(argc>=3);Witness w(argv[1]);cout<<setprecision(18);
    for(int i=2;i<argc;++i){Tree t(argv[i]);sample(w,t,argv[i]);}}
