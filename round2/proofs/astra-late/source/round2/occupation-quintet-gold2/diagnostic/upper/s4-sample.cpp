#define main diagnostic_main
#include "../census.cpp"
#undef main

int band(I d){return d<=5?5:d<=11?11:d<=18?18:d<=37?37:d<=179?179:d<=1998?1998:9999;}

int main(int argc,char** argv){
    assert(argc==11);Witness w(argv[1]);assert(w.h.mode==1);cout<<setprecision(18);
    for(int i=2;i<argc;++i){
        Tree t(argv[i]);Ledger l(w,t);
        map<array<int,13>,Agg> groups;
        for(int m=0;m<(int)t.m.size();++m){
            auto hub=t.m[m];if(Witness::index(w.z.MS4,hub)<0)continue;
            const auto& a=t.a[m];
            for(size_t j=0;j<a.size();++j)for(size_t k=j+1;k<a.size();++k)
                for(size_t n=k+1;n<a.size();++n)for(size_t o=n+1;o<a.size();++o){
                    array<int,4> tips{a[j],a[k],a[n],a[o]};
                    Set s;s.push(m);for(int v:tips)s.push(v);
                    Shape q=classify(s,t);assert(q.k==T&&q.in.v[0]==m);
                    R prior=l.negative[T];l.row(s,q);R reward=l.negative[T]-prior;
                    auto p=find(tips.begin(),tips.end(),t.par[m]);bool has=p!=tips.end();
                    array<I,4> ds;transform(tips.begin(),tips.end(),ds.begin(),[&](int v){return t.m[v].d;});
                    sort(ds.begin(),ds.end());
                    int old=ds[2]==ds[3]?2:(has&&t.m[*p].d==ds[3]?0:1);
                    int rank=0;
                    if(old==1){
                        auto v=*find_if(tips.begin(),tips.end(),[&](int x){return t.m[x].d==ds[3];});
                        rank=1+count_if(tips.begin(),tips.end(),[&](int x){return x<v&&x!=t.par[m];});
                    }
                    I young=999999999;int num2=0,young_child=0;
                    for(int v:tips)if(v!=t.par[m]){
                        young=min(young,t.m[v].d);num2+=t.m[v].d==2;young_child+=t.m[v].d<=5;}
                    auto& g=groups[{s.has(0),has,old,rank,int(hub.d),hub.c,band(ds[3]),band(ds[2]),
                        band(young),num2,young_child,has?band(t.m[*p].d):0,int(ds[3]>=19)}];
                    ++g.n;g.amount+=reward;
                }
        }
        cout<<"tree\t"<<argv[i]<<'\t'<<t.m.size()<<'\t'<<l.count[T]<<'\t'<<l.negative[T]/w.h.S<<'\n';
        for(auto [key,g]:groups){cout<<"bin";for(int v:key)cout<<'\t'<<v;
            cout<<'\t'<<g.n<<'\t'<<g.amount/w.h.S<<'\n';}
    }
}
