#define main diagnostic_main
#include "../census.cpp"
#undef main

int band(I d){return d<=5?5:d<=11?11:d<=18?18:d<=37?37:d<=179?179:d<=1998?1998:9999;}

int main(int argc,char** argv){
    assert(argc==11);Witness w(argv[1]);assert(w.h.mode==1);cout<<setprecision(18);
    for(int i=2;i<argc;++i){
        Tree t(argv[i]);Ledger l(w,t);
        map<array<int,13>,Agg> groups;
        for(int h=0;h<(int)t.m.size();++h){
            auto Hm=t.m[h];if(Witness::index(w.z.I5_3,Hm)<0)continue;
            const auto& a=t.a[h];
            for(int m:a){
                auto M=t.m[m];if(Witness::index(w.z.I5_2c,M)<0||Hm.d+M.d>w.h.SB5)continue;
                for(int y:t.a[m])if(y!=h){
                    for(size_t j=0;j<a.size();++j)for(size_t k=j+1;k<a.size();++k){
                        int x=a[j],z=a[k];if(x==m||z==m)continue;
                        Set s;for(int v:{h,m,y,x,z})s.push(v);
                        Shape q=classify(s,t);assert(q.k==B&&q.in.v[0]==h&&q.in.v[1]==m);
                        R prior=l.negative[B];l.row(s,q);R reward=l.negative[B]-prior;
                        assert(w.supported(q,t));
                        int root=t.par[h]==m?(t.par[m]==y?3:2):(t.par[h]==x||t.par[h]==z?1:0);
                        I dx=t.m[x].d,dz=t.m[z].d,dy=t.m[y].d;
                        I mx=max({dx,dz,dy}),mn=min({dx,dz,dy});
                        int high=(dx==mx)+(dz==mx)+(dy==mx)>1?3:dy==mx?2:0;
                        // Identify the far/short arm's direct relation to its own parent.
                        if(high==0){int u=dx==mx?x:z;high=t.par[h]==u?1:0;}
                        int young=(dx<=5)+(dz<=5)+(dy<=5);
                        auto& g=groups[{s.has(0),root,high,int(Hm.d),Hm.c,int(M.d),M.c,band(mx),band(mn),
                            band(dx+dz+dy-mx-mn),young,int(mx>=19),int(s.has(t.par[m]))}];
                        ++g.n;g.amount+=reward;
                    }
                }
            }
        }
        cout<<"tree\t"<<argv[i]<<'\t'<<t.m.size()<<'\t'<<l.count[B]<<'\t'<<l.negative[B]/w.h.S<<'\n';
        for(auto [key,g]:groups){cout<<"bin";for(int v:key)cout<<'\t'<<v;
            cout<<'\t'<<g.n<<'\t'<<g.amount/w.h.S<<'\n';}
    }
}
