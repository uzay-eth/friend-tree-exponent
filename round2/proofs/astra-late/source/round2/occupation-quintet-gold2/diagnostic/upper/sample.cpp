#define main diagnostic_main
#include "../census.cpp"
#undef main

int main(int argc,char** argv){
    assert(argc==11); Witness w(argv[1]);assert(w.h.mode==1);cout<<setprecision(18);
    for(int i=2;i<argc;++i){
        Tree t(argv[i]);Ledger l(w,t),wide(w,t);
        for(int m=1;m<(int)t.m.size();++m){
            int h=t.par[m];if(h<=0||t.m[h].d<80||t.m[h].d>1998||t.m[m].d!=2)continue;
            for(int y:t.a[m]){
                if(t.par[y]!=m||t.m[y].d>=t.m[h].d)continue;
                Set s;s.push(h);s.push(m);s.push(y);auto q=classify(s,t);assert(q.k==H);
                if(t.m[h].d<=420&&t.m[h].c<=12&&t.m[y].d<=8&&t.m[y].c<=2)l.row(s,q);
                if(t.m[h].d<=480&&t.m[h].c<=12&&t.m[y].d<=12&&t.m[y].c<=3)wide.row(s,q);
            }
        }
        cout<<argv[i]<<'\t'<<t.m.size()<<'\t'<<l.count[H]<<'\t'<<l.negative[H]/w.h.S
            <<'\t'<<wide.count[H]<<'\t'<<wide.negative[H]/w.h.S<<'\n';
    }
}
